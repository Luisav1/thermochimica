!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleMQMQAResponseMapping.f90
!> \brief   Build, aggregate, apply, and transactionally solve verified phase-local MQMQA GEM corrections.
!>
!> \details MQ-4B packages the diagnostic MQ-4A derivation as reusable software. MQ-4C adds strict active-phase
!!          routing, all-or-nothing aggregation, and a baseline-preserving fixed-alpha linear transaction used
!!          by a default-off GEMNewton path. For one active, uncharged, interior plain-SUBG or SUBQ phase, a
!!          model-specific builder
!!          upgrades the historical ideal composition response to the complete supported MQMQA response and
!!          returns only their difference:
!!
!!              deltaA = N*S^T*(Rcorr-Rbase)
!!              deltaB = N*S^T*(rMuCorr-rMuBase).
!!
!!          The builder never changes Thermochimica global state or a GEM matrix. The separate applicator changes
!!          only the caller-owned element block and element residual. Alpha weights the completed correction,
!!          never the local Hessian. The separate default-off MQ-4D trust layer selects how much of the
!!          completed correction pair to apply.
!-------------------------------------------------------------------------------------------------------------

module ModuleMQMQAResponseMapping

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleThermo, ONLY: cSolnPhaseType, dChemicalPotential, dMolesPhase, dMolFraction, dStoichSpecies, &
        iAssemblage, iParticlesPerMole, iPhaseElectronID, nElements, nSpeciesPhase
    USE ModuleMQMQAUnconstrained, ONLY: MQMQAModelData, MQMQAInteractionTerm, &
        CompMQMQAHessianUnconstrained
    USE ModuleMQMQAProductionAdapter, ONLY: DecodeProductionSUBGPhase, DecodeProductionSUBQPhase
    USE ModuleConstrainedResponse, ONLY: SolveConstrainedResponse
    USE ModuleGEMNewtonDiagnosticCapture, ONLY: CaptureGEMNewtonCorrectedSystem

    implicit none
    private

    integer, parameter, public :: MQMQA_MAP_SUCCESS = 0
    integer, parameter, public :: MQMQA_MAP_NOT_APPLICABLE = 1
    integer, parameter, public :: MQMQA_MAP_UNSUPPORTED_CHARGED_PHASE = 2
    integer, parameter, public :: MQMQA_MAP_INVALID_INPUT = 3
    integer, parameter, public :: MQMQA_MAP_DECODE_FAILURE = 4
    integer, parameter, public :: MQMQA_MAP_HESSIAN_FAILURE = 5
    integer, parameter, public :: MQMQA_MAP_CORRECTED_ELEMENT_RESPONSE_FAILURE = 6
    integer, parameter, public :: MQMQA_MAP_BASELINE_ELEMENT_RESPONSE_FAILURE = 7
    integer, parameter, public :: MQMQA_MAP_CORRECTED_RESIDUAL_RESPONSE_FAILURE = 8
    integer, parameter, public :: MQMQA_MAP_BASELINE_RESIDUAL_RESPONSE_FAILURE = 9
    integer, parameter, public :: MQMQA_MAP_INVALID_CORRECTION = 10
    integer, parameter, public :: MQMQA_MAP_INVALID_APPLICATION = 11
    integer, parameter, public :: MQMQA_MAP_OUTSIDE_INTERIOR = 12
    real(8), parameter, public :: MQMQA_INTERIOR_MINIMUM = 1D-12
    integer, parameter, public :: MQMQA_AGGREGATE_SUCCESS = 0
    integer, parameter, public :: MQMQA_AGGREGATE_NO_APPLICABLE_PHASE = 1
    integer, parameter, public :: MQMQA_AGGREGATE_INVALID_INPUT = 2
    integer, parameter, public :: MQMQA_AGGREGATE_PHASE_FAILURE = 3
    integer, parameter, public :: MQMQA_TRIAL_ACCEPTED = 0
    integer, parameter, public :: MQMQA_TRIAL_APPLICATION_FALLBACK = 1
    integer, parameter, public :: MQMQA_TRIAL_DGESV_FALLBACK = 2
    integer, parameter, public :: MQMQA_TRIAL_NONFINITE_FALLBACK = 3
    integer, parameter, public :: MQMQA_TRIAL_BASELINE_FAILURE = 4

    public :: BuildMQMQAGEMCorrection, BuildMQMQASUBQGEMCorrection
    public :: BuildMQMQAReducedCorrection, ApplyMQMQAGEMCorrection
    public :: BuildActiveMQMQAGEMCorrection, AggregateMQMQACorrectionPairs
    public :: SolveMQMQACorrectionTrial

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Route every active SUBG/SUBQ phase and build one all-or-nothing reduced-GEM correction.
    !>
    !> \details Unrelated models are ignored before routing. Unsupported charged phases are counted as deliberate
    !!          exclusions. Every other failure from a routed model invalidates the complete aggregate, preventing
    !!          a partially corrected GEM system when one eligible phase cannot provide its response.
    !---------------------------------------------------------------------------------------------------------
    subroutine BuildActiveMQMQAGEMCorrection(nActive,dDeltaA,dDeltaB,lSupportedFound,lEligible, &
        nAccepted,nCharged,iFailurePhase,iFailureStatus,dFailureMinimumFraction,iStatus)

        integer, intent(in) :: nActive
        real(8), intent(out) :: dDeltaA(:,:), dDeltaB(:)
        logical, intent(out) :: lSupportedFound, lEligible
        integer, intent(out) :: nAccepted, nCharged, iFailurePhase, iFailureStatus, iStatus
        real(8), intent(out) :: dFailureMinimumFraction

        integer :: iPhase, iSlot, k, nRouted
        integer, allocatable :: iPhaseID(:), iPhaseStatus(:)
        logical, allocatable :: lApplicable(:)
        real(8), allocatable :: dPhaseA(:,:,:), dPhaseB(:,:)

        dDeltaA = 0D0
        dDeltaB = 0D0
        lSupportedFound = .FALSE.
        lEligible = .FALSE.
        nAccepted = 0
        nCharged = 0
        iFailurePhase = 0
        iFailureStatus = MQMQA_MAP_SUCCESS
        dFailureMinimumFraction = HUGE(1D0)
        iStatus = MQMQA_AGGREGATE_INVALID_INPUT
        if ((nActive < 0) .OR. (SIZE(dDeltaA,1) /= nElements) .OR. &
            (SIZE(dDeltaA,2) /= nElements) .OR. (SIZE(dDeltaB) /= nElements)) return
        if ((.NOT. ALLOCATED(iAssemblage)) .OR. (.NOT. ALLOCATED(cSolnPhaseType))) return

        if (nActive == 0) then
            iStatus = MQMQA_AGGREGATE_NO_APPLICABLE_PHASE
            return
        end if
        allocate(dPhaseA(nElements,nElements,nActive),dPhaseB(nElements,nActive), &
            lApplicable(nActive),iPhaseStatus(nActive),iPhaseID(nActive))
        dPhaseA = 0D0
        dPhaseB = 0D0
        lApplicable = .FALSE.
        iPhaseStatus = MQMQA_MAP_NOT_APPLICABLE
        iPhaseID = 0

        nRouted = 0
        do k = 1, nActive
            iSlot = nElements-k+1
            if ((iSlot < 1) .OR. (iSlot > SIZE(iAssemblage))) cycle
            iPhase = -iAssemblage(iSlot)
            if ((iPhase < 1) .OR. (iPhase > SIZE(cSolnPhaseType))) cycle
            if ((cSolnPhaseType(iPhase) /= 'SUBG') .AND. &
                (cSolnPhaseType(iPhase) /= 'SUBQ')) cycle

            lSupportedFound = .TRUE.
            nRouted = nRouted+1
            iPhaseID(nRouted) = iPhase
            select case(cSolnPhaseType(iPhase))
            case('SUBG')
                call BuildMQMQAGEMCorrection(iPhase,iSlot,dPhaseA(:,:,nRouted),dPhaseB(:,nRouted), &
                    lApplicable(nRouted),iPhaseStatus(nRouted))
            case('SUBQ')
                call BuildMQMQASUBQGEMCorrection(iPhase,iSlot,dPhaseA(:,:,nRouted),dPhaseB(:,nRouted), &
                    lApplicable(nRouted),iPhaseStatus(nRouted))
            end select
        end do

        if (nRouted == 0) then
            iStatus = MQMQA_AGGREGATE_NO_APPLICABLE_PHASE
            return
        end if
        call AggregateMQMQACorrectionPairs(dPhaseA(:,:,1:nRouted),dPhaseB(:,1:nRouted), &
            lApplicable(1:nRouted),iPhaseStatus(1:nRouted),dDeltaA,dDeltaB,nAccepted,nCharged, &
            k,iFailureStatus,iStatus)
        if (iStatus == MQMQA_AGGREGATE_PHASE_FAILURE) then
            iFailurePhase = iPhaseID(k)
            if ((iFailureStatus == MQMQA_MAP_OUTSIDE_INTERIOR) .AND. (iFailurePhase > 0)) then
                dFailureMinimumFraction = MINVAL(dMolFraction( &
                    nSpeciesPhase(iFailurePhase-1)+1:nSpeciesPhase(iFailurePhase)))
            end if
        end if
        lEligible = (iStatus == MQMQA_AGGREGATE_SUCCESS) .AND. (nAccepted > 0)

    end subroutine BuildActiveMQMQAGEMCorrection


    !---------------------------------------------------------------------------------------------------------
    !> \brief Combine prebuilt phase-local correction pairs using MQ-4C transactional failure semantics.
    !>
    !> \details This state-free helper is shared by live routing and controlled aggregation tests. A charged
    !!          exclusion is skipped deliberately; every other unsuccessful routed pair erases all prior sums.
    !---------------------------------------------------------------------------------------------------------
    subroutine AggregateMQMQACorrectionPairs(dPhaseA,dPhaseB,lApplicable,iPhaseStatus,dDeltaA,dDeltaB, &
        nAccepted,nCharged,iFailurePair,iFailureStatus,iStatus)

        real(8), intent(in) :: dPhaseA(:,:,:), dPhaseB(:,:)
        logical, intent(in) :: lApplicable(:)
        integer, intent(in) :: iPhaseStatus(:)
        real(8), intent(out) :: dDeltaA(:,:), dDeltaB(:)
        integer, intent(out) :: nAccepted, nCharged, iFailurePair, iFailureStatus, iStatus

        integer :: k, nPair

        dDeltaA = 0D0
        dDeltaB = 0D0
        nAccepted = 0
        nCharged = 0
        iFailurePair = 0
        iFailureStatus = MQMQA_MAP_SUCCESS
        iStatus = MQMQA_AGGREGATE_INVALID_INPUT
        nPair = SIZE(dPhaseA,3)
        if ((nPair < 1) .OR. (SIZE(dPhaseA,1) /= SIZE(dDeltaA,1)) .OR. &
            (SIZE(dPhaseA,2) /= SIZE(dDeltaA,2)) .OR. (SIZE(dPhaseB,1) /= SIZE(dDeltaB)) .OR. &
            (SIZE(dPhaseB,2) /= nPair) .OR. (SIZE(lApplicable) /= nPair) .OR. &
            (SIZE(iPhaseStatus) /= nPair)) return

        do k = 1, nPair
            if (lApplicable(k) .AND. (iPhaseStatus(k) == MQMQA_MAP_SUCCESS)) then
                if ((.NOT. ALL(IEEE_IS_FINITE(dPhaseA(:,:,k)))) .OR. &
                    (.NOT. ALL(IEEE_IS_FINITE(dPhaseB(:,k))))) then
                    iFailurePair = k
                    iFailureStatus = MQMQA_MAP_INVALID_CORRECTION
                    iStatus = MQMQA_AGGREGATE_PHASE_FAILURE
                    dDeltaA = 0D0
                    dDeltaB = 0D0
                    return
                end if
                dDeltaA = dDeltaA+dPhaseA(:,:,k)
                dDeltaB = dDeltaB+dPhaseB(:,k)
                nAccepted = nAccepted+1
            else if ((.NOT. lApplicable(k)) .AND. &
                (iPhaseStatus(k) == MQMQA_MAP_UNSUPPORTED_CHARGED_PHASE)) then
                nCharged = nCharged+1
            else
                iFailurePair = k
                iFailureStatus = iPhaseStatus(k)
                iStatus = MQMQA_AGGREGATE_PHASE_FAILURE
                dDeltaA = 0D0
                dDeltaB = 0D0
                nAccepted = 0
                return
            end if
        end do

        if (nAccepted > 0) then
            iStatus = MQMQA_AGGREGATE_SUCCESS
        else
            iStatus = MQMQA_AGGREGATE_NO_APPLICABLE_PHASE
        end if

    end subroutine AggregateMQMQACorrectionPairs

    !---------------------------------------------------------------------------------------------------------
    !> \brief Build the unscaled reduced-GEM correction for one live plain-SUBG phase.
    !>
    !> \param[in]  iPhaseIndex Absolute production solution-phase index.
    !> \param[in]  iPhaseSlot  Slot in `iAssemblage` and `dMolesPhase` containing this active phase.
    !> \param[out] dDeltaA     Element-by-element response correction.
    !> \param[out] dDeltaB     Element residual correction.
    !> \param[out] lApplicable True once the phase passes model/domain eligibility, even if mapping later fails.
    !> \param[out] iStatus     `MQMQA_MAP_SUCCESS`, an inapplicability code, or the precise failed stage.
    !>
    !> \details `dChemicalPotential` must already contain the complete production partial molars produced by
    !!          `CompChemicalPotential`; adding `dPartialExcessGibbs` here would double count excess terms.
    !!          Outputs are cleared before every check, so no failure can expose a partial correction.
    !---------------------------------------------------------------------------------------------------------
    subroutine BuildMQMQAGEMCorrection(iPhaseIndex,iPhaseSlot,dDeltaA,dDeltaB,lApplicable,iStatus)

        integer, intent(in) :: iPhaseIndex, iPhaseSlot
        real(8), intent(out) :: dDeltaA(:,:), dDeltaB(:)
        logical, intent(out) :: lApplicable
        integer, intent(out) :: iStatus

        call BuildMQMQAProductionCorrection(iPhaseIndex,iPhaseSlot,'SUBG',dDeltaA,dDeltaB, &
            lApplicable,iStatus)

    end subroutine BuildMQMQAGEMCorrection


    !---------------------------------------------------------------------------------------------------------
    !> \brief Build the unscaled reduced-GEM correction for one live SUBQ phase.
    !>
    !> \param[in]  iPhaseIndex Absolute production solution-phase index.
    !> \param[in]  iPhaseSlot  Slot in `iAssemblage` and `dMolesPhase` containing this active phase.
    !> \param[out] dDeltaA     Element-by-element response correction.
    !> \param[out] dDeltaB     Element residual correction.
    !> \param[out] lApplicable True once the phase passes model/domain eligibility, even if mapping later fails.
    !> \param[out] iStatus     `MQMQA_MAP_SUCCESS`, an inapplicability code, or the precise failed stage.
    !>
    !> \details This entry point has the same output-clearing, eligibility, state-preservation, and response
    !!          contract as the plain-SUBG builder, but requires `cSolnPhaseType == 'SUBQ'` and invokes only
    !!          `DecodeProductionSUBQPhase`. It cannot silently reinterpret a SUBG phase as SUBQ.
    !---------------------------------------------------------------------------------------------------------
    subroutine BuildMQMQASUBQGEMCorrection(iPhaseIndex,iPhaseSlot,dDeltaA,dDeltaB,lApplicable,iStatus)

        integer, intent(in) :: iPhaseIndex, iPhaseSlot
        real(8), intent(out) :: dDeltaA(:,:), dDeltaB(:)
        logical, intent(out) :: lApplicable
        integer, intent(out) :: iStatus

        call BuildMQMQAProductionCorrection(iPhaseIndex,iPhaseSlot,'SUBQ',dDeltaA,dDeltaB, &
            lApplicable,iStatus)

    end subroutine BuildMQMQASUBQGEMCorrection


    !---------------------------------------------------------------------------------------------------------
    !> \brief Perform model-independent live-state checks and dispatch to one strict production decoder.
    !>
    !> \details `dChemicalPotential` must already contain the complete production partial molars produced by
    !!          `CompChemicalPotential`; adding `dPartialExcessGibbs` here would double count excess terms.
    !!          Outputs are cleared before every check, so no failure can expose a partial correction.
    !---------------------------------------------------------------------------------------------------------
    subroutine BuildMQMQAProductionCorrection(iPhaseIndex,iPhaseSlot,cExpectedType,dDeltaA,dDeltaB, &
        lApplicable,iStatus)

        integer, intent(in) :: iPhaseIndex, iPhaseSlot
        character(len=*), intent(in) :: cExpectedType
        real(8), intent(out) :: dDeltaA(:,:), dDeltaB(:)
        logical, intent(out) :: lApplicable
        integer, intent(out) :: iStatus

        integer :: e, iFirst, iInfo, iLast, nQuad
        real(8) :: dN
        real(8), allocatable :: dHessian(:,:), dHx(:,:), dMoles(:), dS(:,:), dX(:)
        type(MQMQAModelData) :: tModel
        type(MQMQAInteractionTerm), allocatable :: tInteraction(:)

        dDeltaA = 0D0
        dDeltaB = 0D0
        lApplicable = .FALSE.
        iStatus = MQMQA_MAP_NOT_APPLICABLE

        if ((.NOT. ALLOCATED(cSolnPhaseType)) .OR. (.NOT. ALLOCATED(dChemicalPotential)) .OR. &
            (.NOT. ALLOCATED(dMolesPhase)) .OR. (.NOT. ALLOCATED(dMolFraction)) .OR. &
            (.NOT. ALLOCATED(dStoichSpecies)) .OR. (.NOT. ALLOCATED(iAssemblage)) .OR. &
            (.NOT. ALLOCATED(iParticlesPerMole)) .OR. (.NOT. ALLOCATED(iPhaseElectronID)) .OR. &
            (.NOT. ALLOCATED(nSpeciesPhase))) then
            iStatus = MQMQA_MAP_INVALID_INPUT
            return
        end if
        if ((SIZE(dDeltaA,1) /= nElements) .OR. (SIZE(dDeltaA,2) /= nElements) .OR. &
            (SIZE(dDeltaB) /= nElements)) then
            iStatus = MQMQA_MAP_INVALID_INPUT
            return
        end if
        if ((iPhaseIndex < 1) .OR. (iPhaseIndex > SIZE(cSolnPhaseType)) .OR. &
            (iPhaseSlot < 1) .OR. (iPhaseSlot > SIZE(iAssemblage))) return
        if (iAssemblage(iPhaseSlot) /= -iPhaseIndex) return
        if ((cExpectedType /= 'SUBG') .AND. (cExpectedType /= 'SUBQ')) then
            iStatus = MQMQA_MAP_INVALID_INPUT
            return
        end if
        if (cSolnPhaseType(iPhaseIndex) /= cExpectedType) return
        if (iPhaseElectronID(iPhaseIndex) /= 0) then
            iStatus = MQMQA_MAP_UNSUPPORTED_CHARGED_PHASE
            return
        end if

        iFirst = nSpeciesPhase(iPhaseIndex-1)+1
        iLast = nSpeciesPhase(iPhaseIndex)
        nQuad = iLast-iFirst+1
        if ((nElements < 1) .OR. (nQuad < 2) .OR. (iLast > SIZE(dMolFraction)) .OR. &
            (iLast > SIZE(dChemicalPotential)) .OR. (iLast > SIZE(iParticlesPerMole)) .OR. &
            (iLast > SIZE(dStoichSpecies,1)) .OR. (nElements > SIZE(dStoichSpecies,2)) .OR. &
            (iPhaseSlot > SIZE(dMolesPhase))) then
            iStatus = MQMQA_MAP_INVALID_INPUT
            return
        end if

        dN = dMolesPhase(iPhaseSlot)
        allocate(dX(nQuad),dMoles(nQuad),dS(nQuad,nElements))
        dX = dMolFraction(iFirst:iLast)
        if ((.NOT. IEEE_IS_FINITE(dN)) .OR. (dN <= 0D0) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dX))) .OR. &
            (ABS(SUM(dX)-1D0) > 1D-10) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dChemicalPotential(iFirst:iLast)))) .OR. &
            ANY(iParticlesPerMole(iFirst:iLast) <= 0)) then
            iStatus = MQMQA_MAP_INVALID_INPUT
            return
        end if
        if (MINVAL(dX) <= MQMQA_INTERIOR_MINIMUM) then
            iStatus = MQMQA_MAP_OUTSIDE_INTERIOR
            return
        end if

        lApplicable = .TRUE.
        do e = 1, nElements
            dS(:,e) = dStoichSpecies(iFirst:iLast,e)/DFLOAT(iParticlesPerMole(iFirst:iLast))
        end do
        if (.NOT. ALL(IEEE_IS_FINITE(dS))) then
            iStatus = MQMQA_MAP_INVALID_INPUT
            return
        end if

        select case (cExpectedType)
        case ('SUBG')
            call DecodeProductionSUBGPhase(iPhaseIndex,tModel,tInteraction,iInfo)
        case ('SUBQ')
            call DecodeProductionSUBQPhase(iPhaseIndex,tModel,tInteraction,iInfo)
        end select
        if (iInfo /= 0) then
            iStatus = MQMQA_MAP_DECODE_FAILURE
            return
        end if

        dMoles = dN*dX
        allocate(dHessian(nQuad,nQuad),dHx(nQuad,nQuad))
        call CompMQMQAHessianUnconstrained(tModel,dMoles,1D0,tInteraction,dHessian,iInfo)
        if (iInfo /= 0) then
            iStatus = MQMQA_MAP_HESSIAN_FAILURE
            return
        end if

        dHx = dN*dHessian
        call BuildMQMQAReducedCorrection(dN,dX,dChemicalPotential(iFirst:iLast),dS,dHx, &
            dDeltaA,dDeltaB,iStatus)

    end subroutine BuildMQMQAProductionCorrection


    !---------------------------------------------------------------------------------------------------------
    !> \brief Condense supplied local response data into the verified reduced-GEM correction.
    !>
    !> \details This state-based kernel contains no Thermochimica decoding or global reads. Besides keeping the
    !!          production adapter separate from the response algebra, it permits direct failure-path tests with
    !!          singular or invalid local curvature matrices. Every failure returns zero corrections.
    !---------------------------------------------------------------------------------------------------------
    subroutine BuildMQMQAReducedCorrection(dN,dX,dMu,dS,dHx,dDeltaA,dDeltaB,iStatus)

        real(8), intent(in) :: dN, dX(:), dMu(:), dS(:,:), dHx(:,:)
        real(8), intent(out) :: dDeltaA(:,:), dDeltaB(:)
        integer, intent(out) :: iStatus

        integer :: e, iInfo, nElement, nQuad
        real(8) :: dConstraintResidual, dKKTResidual, dSymmetryResidual
        real(8), allocatable :: dConstraint(:,:), dForceMu(:,:), dHbase(:,:), dMultiplier(:,:)
        real(8), allocatable :: dMuResponse(:,:), dMuResponseBase(:,:), dResponse(:,:), dResponseBase(:,:)

        dDeltaA = 0D0
        dDeltaB = 0D0
        iStatus = MQMQA_MAP_INVALID_INPUT
        nQuad = SIZE(dX)
        nElement = SIZE(dS,2)
        if ((nQuad < 2) .OR. (nElement < 1) .OR. (SIZE(dMu) /= nQuad) .OR. &
            (SIZE(dS,1) /= nQuad) .OR. (SIZE(dHx,1) /= nQuad) .OR. (SIZE(dHx,2) /= nQuad) .OR. &
            (SIZE(dDeltaA,1) /= nElement) .OR. (SIZE(dDeltaA,2) /= nElement) .OR. &
            (SIZE(dDeltaB) /= nElement)) return
        if ((.NOT. IEEE_IS_FINITE(dN)) .OR. (dN <= 0D0) .OR. (.NOT. ALL(IEEE_IS_FINITE(dX))) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dMu))) .OR. (.NOT. ALL(IEEE_IS_FINITE(dS))) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dHx))) .OR. &
            (ABS(SUM(dX)-1D0) > 1D-10)) return
        if (MINVAL(dX) <= MQMQA_INTERIOR_MINIMUM) then
            iStatus = MQMQA_MAP_OUTSIDE_INTERIOR
            return
        end if

        allocate(dHbase(nQuad,nQuad),dConstraint(1,nQuad),dResponse(nQuad,nElement), &
            dResponseBase(nQuad,nElement),dMultiplier(1,nElement),dForceMu(nQuad,1), &
            dMuResponse(nQuad,1),dMuResponseBase(nQuad,1))
        dHbase = 0D0
        do e = 1, nQuad
            dHbase(e,e) = 1D0/dX(e)
        end do
        dConstraint = 1D0

        call SolveConstrainedResponse(dHx,dConstraint,dS,dResponse,iInfo,dMultiplier)
        if (iInfo /= 0) then
            iStatus = MQMQA_MAP_CORRECTED_ELEMENT_RESPONSE_FAILURE
            return
        end if
        dKKTResidual = KKTResidual(dHx,dConstraint,dS,dResponse,dMultiplier)
        dConstraintResidual = MAXVAL(ABS(MATMUL(dConstraint,dResponse)))
        if ((dKKTResidual > 1D-10) .OR. (dConstraintResidual > 1D-10)) then
            iStatus = MQMQA_MAP_CORRECTED_ELEMENT_RESPONSE_FAILURE
            return
        end if

        call SolveConstrainedResponse(dHbase,dConstraint,dS,dResponseBase,iInfo,dMultiplier)
        if (iInfo /= 0) then
            iStatus = MQMQA_MAP_BASELINE_ELEMENT_RESPONSE_FAILURE
            return
        end if
        dKKTResidual = KKTResidual(dHbase,dConstraint,dS,dResponseBase,dMultiplier)
        dConstraintResidual = MAXVAL(ABS(MATMUL(dConstraint,dResponseBase)))
        if ((dKKTResidual > 1D-10) .OR. (dConstraintResidual > 1D-10)) then
            iStatus = MQMQA_MAP_BASELINE_ELEMENT_RESPONSE_FAILURE
            return
        end if

        dForceMu(:,1) = dMu-1D0
        deallocate(dMultiplier)
        allocate(dMultiplier(1,1))
        call SolveConstrainedResponse(dHx,dConstraint,dForceMu,dMuResponse,iInfo,dMultiplier)
        if (iInfo /= 0) then
            iStatus = MQMQA_MAP_CORRECTED_RESIDUAL_RESPONSE_FAILURE
            return
        end if
        dKKTResidual = KKTResidual(dHx,dConstraint,dForceMu,dMuResponse,dMultiplier)
        dConstraintResidual = MAXVAL(ABS(MATMUL(dConstraint,dMuResponse)))
        if ((dKKTResidual > 1D-10) .OR. (dConstraintResidual > 1D-10)) then
            iStatus = MQMQA_MAP_CORRECTED_RESIDUAL_RESPONSE_FAILURE
            return
        end if

        call SolveConstrainedResponse(dHbase,dConstraint,dForceMu,dMuResponseBase,iInfo,dMultiplier)
        if (iInfo /= 0) then
            iStatus = MQMQA_MAP_BASELINE_RESIDUAL_RESPONSE_FAILURE
            return
        end if
        dKKTResidual = KKTResidual(dHbase,dConstraint,dForceMu,dMuResponseBase,dMultiplier)
        dConstraintResidual = MAXVAL(ABS(MATMUL(dConstraint,dMuResponseBase)))
        if ((dKKTResidual > 1D-10) .OR. (dConstraintResidual > 1D-10)) then
            iStatus = MQMQA_MAP_BASELINE_RESIDUAL_RESPONSE_FAILURE
            return
        end if

        dDeltaA = dN*MATMUL(TRANSPOSE(dS),dResponse-dResponseBase)
        dDeltaB = dN*MATMUL(TRANSPOSE(dS),dMuResponse(:,1)-dMuResponseBase(:,1))
        dSymmetryResidual = SQRT(SUM((dDeltaA-TRANSPOSE(dDeltaA))**2))/ &
            DMAX1(1D0,SQRT(SUM(dDeltaA*dDeltaA)))
        if ((.NOT. ALL(IEEE_IS_FINITE(dDeltaA))) .OR. (.NOT. ALL(IEEE_IS_FINITE(dDeltaB))) .OR. &
            (dSymmetryResidual > 1D-10)) then
            dDeltaA = 0D0
            dDeltaB = 0D0
            iStatus = MQMQA_MAP_INVALID_CORRECTION
            return
        end if

        ! The response derivation is symmetric. Remove only the roundoff-level skew that remains after the raw
        ! matrix has passed the symmetry gate, so callers receive an exactly symmetric GEM correction.
        dDeltaA = 0.5D0*(dDeltaA+TRANSPOSE(dDeltaA))

        iStatus = MQMQA_MAP_SUCCESS

    end subroutine BuildMQMQAReducedCorrection


    !---------------------------------------------------------------------------------------------------------
    !> \brief Add a caller-selected fraction of one verified correction to caller-owned GEM arrays.
    !>
    !> \details This routine has no access to GEMNewton globals. It modifies only `A(1:nElement,1:nElement)` and
    !!          `B(1:nElement)`. The weight is a future globalization choice, not a scale applied to the MQMQA
    !!          Hessian itself.
    !---------------------------------------------------------------------------------------------------------
    subroutine ApplyMQMQAGEMCorrection(dA,dB,nElement,dDeltaA,dDeltaB,dAlpha,iStatus)

        real(8), intent(inout) :: dA(:,:), dB(:)
        integer, intent(in) :: nElement
        real(8), intent(in) :: dDeltaA(:,:), dDeltaB(:), dAlpha
        integer, intent(out) :: iStatus

        real(8) :: dSymmetryResidual

        iStatus = MQMQA_MAP_INVALID_APPLICATION
        if ((nElement < 1) .OR. (SIZE(dA,1) < nElement) .OR. (SIZE(dA,2) < nElement) .OR. &
            (SIZE(dB) < nElement) .OR. (SIZE(dDeltaA,1) /= nElement) .OR. &
            (SIZE(dDeltaA,2) /= nElement) .OR. (SIZE(dDeltaB) /= nElement)) return
        if ((.NOT. IEEE_IS_FINITE(dAlpha)) .OR. (dAlpha < 0D0) .OR. (dAlpha > 1D0) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dDeltaA))) .OR. (.NOT. ALL(IEEE_IS_FINITE(dDeltaB)))) return
        dSymmetryResidual = SQRT(SUM((dDeltaA-TRANSPOSE(dDeltaA))**2))/ &
            DMAX1(1D0,SQRT(SUM(dDeltaA*dDeltaA)))
        if (dSymmetryResidual > 1D-10) return

        dA(1:nElement,1:nElement) = dA(1:nElement,1:nElement)+dAlpha*dDeltaA
        dB(1:nElement) = dB(1:nElement)+dAlpha*dDeltaB
        iStatus = MQMQA_MAP_SUCCESS

    end subroutine ApplyMQMQAGEMCorrection


    !---------------------------------------------------------------------------------------------------------
    !> \brief Solve one fixed-alpha correction trial while preserving a complete historical fallback system.
    !>
    !> \details The supplied baseline arrays are never passed to the corrected `DGESV` call because LAPACK
    !!          overwrites both matrix and right-hand side. A valid correction is first applied to private trial
    !!          copies. If application fails, the corrected matrix is singular, or the resulting update is not
    !!          finite, untouched baseline copies are solved instead. This state-free transaction is shared by
    !!          live GEM integration and controlled failure-path tests; model routing remains outside it.
    !---------------------------------------------------------------------------------------------------------
    subroutine SolveMQMQACorrectionTrial(dABase,dBBase,nElement,dDeltaA,dDeltaB,dAlpha, &
        dASolved,dBSolved,iPiv,iInfo,lApplied,lAccepted,iStatus)

        integer, intent(in) :: nElement
        real(8), intent(in) :: dABase(:,:), dBBase(:), dDeltaA(:,:), dDeltaB(:), dAlpha
        real(8), intent(out) :: dASolved(:,:), dBSolved(:)
        integer, intent(out) :: iPiv(:), iInfo, iStatus
        logical, intent(out) :: lApplied, lAccepted

        integer :: iApplyStatus, iTrialInfo, nVar
        integer, allocatable :: iTrialPiv(:)
        real(8), allocatable :: dATrial(:,:), dBTrial(:)

        nVar = SIZE(dBBase)
        dASolved = 0D0
        dBSolved = 0D0
        iPiv = 0
        iInfo = -1
        iStatus = MQMQA_TRIAL_APPLICATION_FALLBACK
        lApplied = .FALSE.
        lAccepted = .FALSE.
        if ((nVar < 1) .OR. (nElement < 1) .OR. (nElement > nVar) .OR. &
            (SIZE(dABase,1) /= nVar) .OR. (SIZE(dABase,2) /= nVar) .OR. &
            (SIZE(dASolved,1) /= nVar) .OR. (SIZE(dASolved,2) /= nVar) .OR. &
            (SIZE(dBSolved) /= nVar) .OR. (SIZE(iPiv) /= nVar) .OR. &
            (SIZE(dDeltaA,1) /= nElement) .OR. (SIZE(dDeltaA,2) /= nElement) .OR. &
            (SIZE(dDeltaB) /= nElement) .OR. (.NOT. ALL(IEEE_IS_FINITE(dABase))) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dBBase)))) return

        allocate(dATrial(nVar,nVar),dBTrial(nVar),iTrialPiv(nVar))
        dASolved = dABase
        dBSolved = dBBase
        dATrial = dABase
        dBTrial = dBBase
        iTrialPiv = 0
        call ApplyMQMQAGEMCorrection(dATrial,dBTrial,nElement,dDeltaA,dDeltaB,dAlpha,iApplyStatus)
        if (iApplyStatus /= MQMQA_MAP_SUCCESS) then
            iStatus = MQMQA_TRIAL_APPLICATION_FALLBACK
            call dgesv(nVar,1,dASolved,nVar,iPiv,dBSolved,nVar,iInfo)
            if ((iInfo /= 0) .OR. (.NOT. ALL(IEEE_IS_FINITE(dBSolved)))) &
                iStatus = MQMQA_TRIAL_BASELINE_FAILURE
            return
        end if

        lApplied = .TRUE.
        call CaptureGEMNewtonCorrectedSystem(dATrial,dBTrial,nVar)
        call dgesv(nVar,1,dATrial,nVar,iTrialPiv,dBTrial,nVar,iTrialInfo)
        if (iTrialInfo /= 0) then
            iStatus = MQMQA_TRIAL_DGESV_FALLBACK
        else if (.NOT. ALL(IEEE_IS_FINITE(dBTrial))) then
            iStatus = MQMQA_TRIAL_NONFINITE_FALLBACK
        else
            dASolved = dATrial
            dBSolved = dBTrial
            iPiv = iTrialPiv
            iInfo = 0
            iStatus = MQMQA_TRIAL_ACCEPTED
            lAccepted = .TRUE.
            return
        end if

        ! The corrected trial has consumed only private arrays, so fallback still sees the exact historical system.
        dASolved = dABase
        dBSolved = dBBase
        iPiv = 0
        call dgesv(nVar,1,dASolved,nVar,iPiv,dBSolved,nVar,iInfo)
        if ((iInfo /= 0) .OR. (.NOT. ALL(IEEE_IS_FINITE(dBSolved)))) &
            iStatus = MQMQA_TRIAL_BASELINE_FAILURE

    end subroutine SolveMQMQACorrectionTrial


    !> Return a scale-normalized residual for the complete bordered response equations.
    real(8) function KKTResidual(dH,dC,dF,dR,dLambda)

        real(8), intent(in) :: dH(:,:), dC(:,:), dF(:,:), dR(:,:), dLambda(:,:)
        real(8) :: dScale
        real(8), allocatable :: dTop(:,:)

        allocate(dTop(SIZE(dF,1),SIZE(dF,2)))
        dTop = MATMUL(dH,dR)+MATMUL(TRANSPOSE(dC),dLambda)-dF
        dScale = DMAX1(1D0,MAXVAL(ABS(dF)),MAXVAL(ABS(dH))*MAXVAL(ABS(dR)),MAXVAL(ABS(dLambda)))
        KKTResidual = MAXVAL(ABS(dTop))/dScale
        deallocate(dTop)

    end function KKTResidual

end module ModuleMQMQAResponseMapping
