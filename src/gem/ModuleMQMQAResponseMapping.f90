!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleMQMQAResponseMapping.f90
!> \brief   Build and apply the verified phase-local MQMQA correction to the reduced GEM equations.
!>
!> \details MQ-4B packages the diagnostic MQ-4A derivation as reusable software without activating it in
!!          GEMNewton. For one active, uncharged, interior plain-SUBG phase, the builder upgrades the historical
!!          ideal composition response to the complete supported MQMQA response and returns only their difference:
!!
!!              deltaA = N*S^T*(Rcorr-Rbase)
!!              deltaB = N*S^T*(rMuCorr-rMuBase).
!!
!!          The builder never changes Thermochimica global state or a GEM matrix. The separate applicator changes
!!          only the caller-owned element block and element residual. Solver activation, trust selection, and
!!          globalization remain MQ-4C work.
!-------------------------------------------------------------------------------------------------------------

module ModuleMQMQAResponseMapping

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleThermo, ONLY: cSolnPhaseType, dChemicalPotential, dMolesPhase, dMolFraction, dStoichSpecies, &
        iAssemblage, iParticlesPerMole, iPhaseElectronID, nElements, nSpeciesPhase
    USE ModuleMQMQAUnconstrained, ONLY: MQMQAModelData, MQMQAInteractionTerm, &
        CompMQMQAHessianUnconstrained
    USE ModuleMQMQAProductionAdapter, ONLY: DecodeProductionSUBGPhase
    USE ModuleConstrainedResponse, ONLY: SolveConstrainedResponse

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

    public :: BuildMQMQAGEMCorrection, BuildMQMQAReducedCorrection, ApplyMQMQAGEMCorrection

contains

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
        if (cSolnPhaseType(iPhaseIndex) /= 'SUBG') return
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
            (.NOT. ALL(IEEE_IS_FINITE(dX))) .OR. (MINVAL(dX) <= 1D-12) .OR. &
            (ABS(SUM(dX)-1D0) > 1D-10) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dChemicalPotential(iFirst:iLast)))) .OR. &
            ANY(iParticlesPerMole(iFirst:iLast) <= 0)) then
            iStatus = MQMQA_MAP_INVALID_INPUT
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

        call DecodeProductionSUBGPhase(iPhaseIndex,tModel,tInteraction,iInfo)
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

    end subroutine BuildMQMQAGEMCorrection


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
            (.NOT. ALL(IEEE_IS_FINITE(dHx))) .OR. (MINVAL(dX) <= 1D-12) .OR. &
            (ABS(SUM(dX)-1D0) > 1D-10)) return

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
