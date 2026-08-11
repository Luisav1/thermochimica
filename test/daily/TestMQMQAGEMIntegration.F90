!-------------------------------------------------------------------------------------------------------------
!> \file    TestMQMQAGEMIntegration.F90
!> \brief   Verify default-off fixed-alpha MQMQA correction plumbing in GEMNewton.
!>
!> \details MQ-4C starts from the already verified local Hessian, constrained response, and reduced mapping.
!!          This test checks only their connection to the live GEM linear system: persistent controls, strict
!!          SUBG/SUBQ routing, all-or-nothing aggregation, baseline-preserving linear fallback, and simultaneous
!!          application of deltaA and deltaB. Fixed alpha one is a linear-solve experiment, not a claim of
!!          nonlinear globalization or guaranteed equilibrium convergence.
!-------------------------------------------------------------------------------------------------------------
program TestMQMQAGEMIntegration

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE, IEEE_VALUE, IEEE_QUIET_NAN, &
        IEEE_POSITIVE_INF, IEEE_NEGATIVE_INF
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleParseCS
    USE ModuleGEMSolver
    USE ModuleMQMQAUnconstrained, ONLY: MQMQAModelData, MQMQA_MODEL_UNSET
    USE ModuleMQMQAResponseMapping
    USE ModuleGEMNewtonDiagnosticCapture

    implicit none

    interface
        subroutine SetMQMQAHessianControls(lEnable,dAlpha,iInfo)
            logical, intent(in) :: lEnable
            real(8), intent(in) :: dAlpha
            integer, intent(out) :: iInfo
        end subroutine SetMQMQAHessianControls
        subroutine ResetMQMQAHessianControls
        end subroutine ResetMQMQAHessianControls
        subroutine SetRKMPHessianControls(lEnable,dAlphaMax,lDebug,iInfo)
            logical, intent(in) :: lEnable, lDebug
            real(8), intent(in) :: dAlphaMax
            integer, intent(out) :: iInfo
        end subroutine SetRKMPHessianControls
        subroutine ResetRKMPHessianControls
        end subroutine ResetRKMPHessianControls
        subroutine GEMNewton(iInfo)
            integer, intent(out) :: iInfo
        end subroutine GEMNewton
        subroutine CompChemicalPotential(lCompEverything)
            logical, intent(in) :: lCompEverything
        end subroutine CompChemicalPotential
    end interface

    logical :: lAlphaOneReport, lHaveCalculation, lPass, lReport
    character(len=32) :: cArgument

    lPass = .TRUE.
    lReport = .FALSE.
    lAlphaOneReport = .FALSE.
    lHaveCalculation = .FALSE.
    if (COMMAND_ARGUMENT_COUNT() > 0) then
        call GET_COMMAND_ARGUMENT(1,cArgument)
        lReport = TRIM(cArgument) == '--report'
        lAlphaOneReport = TRIM(cArgument) == '--alpha-one-report'
    end if

    call CheckControlContract(lPass,lReport)
    call CheckStateFreeTransactions(lPass,lReport)
    call CheckDefaultAndAlphaZero(lPass,lReport)
    call CheckStrictSUBGRouting(lPass,lReport)
    call CheckLiveSUBQIntegration(lPass,lReport)
    call CheckRemainingDualControlCases(lPass,lReport)
    if (lAlphaOneReport) call ReportFullAlphaOneEvidence(lPass)
    call FinishTest(lPass)

contains

    ! Verify that an explicit formulation and valid controls are required, and invalid setters are atomic.
    subroutine CheckControlContract(lAllPass,lDetailed)

        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lDetailed
        integer :: i, iInfo
        real(8) :: dInvalidAlpha(5)
        type(MQMQAModelData) :: tDefaultModel

        call ResetMQMQAHessianControls
        lAllPass = lAllPass .AND. (.NOT. lMQMQAHessianControlsConfigured) .AND. &
            (.NOT. lMQMQAHessianRequestedEnable) .AND. (dMQMQAHessianRequestedAlpha == 0D0) .AND. &
            (tDefaultModel%iModelType == MQMQA_MODEL_UNSET)

        call SetMQMQAHessianControls(.TRUE.,0.35D0,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0) .AND. lMQMQAHessianControlsConfigured .AND. &
            lMQMQAHessianRequestedEnable .AND. (dMQMQAHessianRequestedAlpha == 0.35D0)
        dInvalidAlpha(1:2) = [-0.1D0,1.1D0]
        dInvalidAlpha(3) = IEEE_VALUE(0D0,IEEE_QUIET_NAN)
        dInvalidAlpha(4) = IEEE_VALUE(0D0,IEEE_POSITIVE_INF)
        dInvalidAlpha(5) = IEEE_VALUE(0D0,IEEE_NEGATIVE_INF)
        do i = 1, SIZE(dInvalidAlpha)
            call SetMQMQAHessianControls(.FALSE.,dInvalidAlpha(i),iInfo)
            lAllPass = lAllPass .AND. (iInfo /= 0) .AND. lMQMQAHessianRequestedEnable .AND. &
                (dMQMQAHessianRequestedAlpha == 0.35D0)
        end do
        call SetMQMQAHessianControls(.FALSE.,0.20D0,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0) .AND. (.NOT. lMQMQAHessianRequestedEnable) .AND. &
            (dMQMQAHessianRequestedAlpha == 0.20D0)
        call ResetMQMQAHessianControls
        lAllPass = lAllPass .AND. (.NOT. lMQMQAHessianControlsConfigured) .AND. &
            (.NOT. lMQMQAHessianRequestedEnable) .AND. (dMQMQAHessianRequestedAlpha == 0D0)

        if (lDetailed) write(*,'(A)') 'controls: valid/persistent/atomic-invalid/reset and unset-model checks passed'

    end subroutine CheckControlContract


    ! Exercise aggregation and corrected-solve fallback without relying on a particular assessed assemblage.
    subroutine CheckStateFreeTransactions(lAllPass,lDetailed)

        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lDetailed
        integer :: iFailurePair, iFailureStatus, iInfo, iStatus, nAccepted, nCharged
        integer :: iPiv(2), iPhaseStatus(2)
        logical :: lAccepted, lApplicable(2), lApplied
        real(8) :: dABase(2,2), dASolved(2,2), dDeltaA(2,2), dDeltaB(2)
        real(8) :: dBBase(2), dBSolved(2), dExpectedA(2,2), dExpectedB(2)
        real(8) :: dBoundaryA(1,1), dBoundaryB(1), dBoundaryHx(2,2), dBoundaryMu(2)
        real(8) :: dBoundaryS(2,1), dBoundaryX(2)
        real(8) :: dPhaseA(2,2,2), dPhaseB(2,2), dAggregateA(2,2), dAggregateB(2)

        dPhaseA = 0D0
        dPhaseB = 0D0
        dPhaseA(1,1,1) = 1D0
        dPhaseA(2,2,1) = 2D0
        dPhaseA(1,2,2) = 0.5D0
        dPhaseA(2,1,2) = 0.5D0
        dPhaseB(:,1) = [3D0,4D0]
        dPhaseB(:,2) = [-1D0,2D0]
        lApplicable = .TRUE.
        iPhaseStatus = MQMQA_MAP_SUCCESS
        call AggregateMQMQACorrectionPairs(dPhaseA,dPhaseB,lApplicable,iPhaseStatus,dAggregateA,dAggregateB, &
            nAccepted,nCharged,iFailurePair,iFailureStatus,iStatus)
        lAllPass = lAllPass .AND. (iStatus == MQMQA_AGGREGATE_SUCCESS) .AND. (nAccepted == 2) .AND. &
            (nCharged == 0) .AND. (MAXVAL(ABS(dAggregateA-SUM(dPhaseA,DIM=3))) == 0D0) .AND. &
            (MAXVAL(ABS(dAggregateB-SUM(dPhaseB,DIM=2))) == 0D0)

        iPhaseStatus(2) = MQMQA_MAP_HESSIAN_FAILURE
        call AggregateMQMQACorrectionPairs(dPhaseA,dPhaseB,lApplicable,iPhaseStatus,dAggregateA,dAggregateB, &
            nAccepted,nCharged,iFailurePair,iFailureStatus,iStatus)
        lAllPass = lAllPass .AND. (iStatus == MQMQA_AGGREGATE_PHASE_FAILURE) .AND. &
            (iFailurePair == 2) .AND. (iFailureStatus == MQMQA_MAP_HESSIAN_FAILURE) .AND. &
            (MAXVAL(ABS(dAggregateA)) == 0D0) .AND. (MAXVAL(ABS(dAggregateB)) == 0D0)

        lApplicable = [.FALSE.,.TRUE.]
        iPhaseStatus = [MQMQA_MAP_UNSUPPORTED_CHARGED_PHASE,MQMQA_MAP_SUCCESS]
        call AggregateMQMQACorrectionPairs(dPhaseA,dPhaseB,lApplicable,iPhaseStatus,dAggregateA,dAggregateB, &
            nAccepted,nCharged,iFailurePair,iFailureStatus,iStatus)
        lAllPass = lAllPass .AND. (iStatus == MQMQA_AGGREGATE_SUCCESS) .AND. &
            (nAccepted == 1) .AND. (nCharged == 1) .AND. &
            (MAXVAL(ABS(dAggregateA-dPhaseA(:,:,2))) == 0D0)

        dABase = 0D0
        dABase(1,1) = 1D0
        dABase(2,2) = 1D0
        dBBase = [2D0,3D0]
        dDeltaA = 0D0
        dDeltaA(1,1) = 0.25D0
        dDeltaA(2,2) = 0.50D0
        dDeltaB = [0.5D0,-0.5D0]
        dExpectedA = dABase+dDeltaA
        dExpectedB = dBBase+dDeltaB
        call SolveMQMQACorrectionTrial(dABase,dBBase,2,dDeltaA,dDeltaB,1D0,dASolved,dBSolved,iPiv, &
            iInfo,lApplied,lAccepted,iStatus)
        lAllPass = lAllPass .AND. lApplied .AND. lAccepted .AND. (iInfo == 0) .AND. &
            (iStatus == MQMQA_TRIAL_ACCEPTED) .AND. &
            (MAXVAL(ABS(MATMUL(dExpectedA,dBSolved)-dExpectedB)) <= 1D-14)

        ! A materially nonsymmetric supplied correction must fail before touching the trial and use baseline.
        dDeltaA = 0D0
        dDeltaA(1,2) = 1D0
        dDeltaB = 0D0
        call SolveMQMQACorrectionTrial(dABase,dBBase,2,dDeltaA,dDeltaB,1D0,dASolved,dBSolved,iPiv, &
            iInfo,lApplied,lAccepted,iStatus)
        lAllPass = lAllPass .AND. (.NOT. lApplied) .AND. (.NOT. lAccepted) .AND. (iInfo == 0) .AND. &
            (iStatus == MQMQA_TRIAL_APPLICATION_FALLBACK) .AND. (MAXVAL(ABS(dBSolved-dBBase)) == 0D0)

        ! Make only the corrected matrix singular. The untouched identity baseline must still solve exactly.
        dDeltaA = -dABase
        dDeltaB = 0D0
        call SolveMQMQACorrectionTrial(dABase,dBBase,2,dDeltaA,dDeltaB,1D0,dASolved,dBSolved,iPiv, &
            iInfo,lApplied,lAccepted,iStatus)
        lAllPass = lAllPass .AND. lApplied .AND. (.NOT. lAccepted) .AND. (iInfo == 0) .AND. &
            (iStatus == MQMQA_TRIAL_DGESV_FALLBACK) .AND. (MAXVAL(ABS(dBSolved-dBBase)) == 0D0)

        ! A finite but extreme corrected right-hand side overflows only the trial update; baseline remains finite.
        dABase = 0D0
        dABase(1,1) = 1D-200
        dABase(2,2) = 1D-200
        dBBase = 1D0
        dDeltaA = 0D0
        dDeltaB = 0.5D0*HUGE(1D0)
        call SolveMQMQACorrectionTrial(dABase,dBBase,2,dDeltaA,dDeltaB,1D0,dASolved,dBSolved,iPiv, &
            iInfo,lApplied,lAccepted,iStatus)
        lAllPass = lAllPass .AND. lApplied .AND. (.NOT. lAccepted) .AND. (iInfo == 0) .AND. &
            (iStatus == MQMQA_TRIAL_NONFINITE_FALLBACK) .AND. ALL(IEEE_IS_FINITE(dBSolved))

        ! Freeze the strict-interior contract immediately below and above the unchanged 1E-12 threshold.
        dBoundaryMu = 0D0
        dBoundaryS(:,1) = [1D0,-1D0]
        dBoundaryX = [0.5D0*MQMQA_INTERIOR_MINIMUM,1D0-0.5D0*MQMQA_INTERIOR_MINIMUM]
        dBoundaryHx = 0D0
        dBoundaryHx(1,1) = 1D0/dBoundaryX(1)
        dBoundaryHx(2,2) = 1D0/dBoundaryX(2)
        call BuildMQMQAReducedCorrection(1D0,dBoundaryX,dBoundaryMu,dBoundaryS,dBoundaryHx, &
            dBoundaryA,dBoundaryB,iStatus)
        lAllPass = lAllPass .AND. (iStatus == MQMQA_MAP_OUTSIDE_INTERIOR) .AND. &
            (MAXVAL(ABS(dBoundaryA)) == 0D0) .AND. (MAXVAL(ABS(dBoundaryB)) == 0D0)
        dBoundaryX = [1.5D0*MQMQA_INTERIOR_MINIMUM,1D0-1.5D0*MQMQA_INTERIOR_MINIMUM]
        dBoundaryHx = 0D0
        dBoundaryHx(1,1) = 1D0/dBoundaryX(1)
        dBoundaryHx(2,2) = 1D0/dBoundaryX(2)
        call BuildMQMQAReducedCorrection(1D0,dBoundaryX,dBoundaryMu,dBoundaryS,dBoundaryHx, &
            dBoundaryA,dBoundaryB,iStatus)
        lAllPass = lAllPass .AND. (iStatus == MQMQA_MAP_SUCCESS)

        if (lDetailed) write(*,'(A)') &
            'transactions: aggregation, three fallbacks, and interior-threshold classification passed'

    end subroutine CheckStateFreeTransactions


    ! Alpha zero is a control-level bypass, so it must reproduce a default-off calculation exactly.
    subroutine CheckDefaultAndAlphaZero(lAllPass,lDetailed)

        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lDetailed
        integer :: iInfo, iInfoBase
        real(8) :: dGibbsBase
        real(8), allocatable :: dFractionBase(:), dMolesBase(:), dPhaseBase(:)

        call ResetMQMQAHessianControls
        call PrepareFeTiVO
        call Thermochimica
        iInfoBase = INFOThermo
        lAllPass = lAllPass .AND. (iInfoBase == 0) .AND. (.NOT. lUseMQMQAExactHessian) .AND. &
            (nMQMQAHessianApplyCount == 0) .AND. (nMQMQAHessianAcceptedSolveCount == 0)
        if (iInfoBase /= 0) return
        dGibbsBase = dGibbsEnergySys
        allocate(dFractionBase(SIZE(dMolFraction)),dMolesBase(SIZE(dMolesSpecies)), &
            dPhaseBase(SIZE(dMolesPhase)))
        dFractionBase = dMolFraction
        dMolesBase = dMolesSpecies
        dPhaseBase = dMolesPhase

        call SetMQMQAHessianControls(.TRUE.,0D0,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0)
        call SetRKMPHessianControls(.TRUE.,1D0,.FALSE.,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0)
        call PrepareFeTiVO
        call Thermochimica
        lAllPass = lAllPass .AND. (INFOThermo == iInfoBase) .AND. lUseMQMQAExactHessian .AND. &
            (dMQMQAHessianAlpha == 0D0) .AND. (dGibbsEnergySys == dGibbsBase) .AND. &
            (SIZE(dMolFraction) == SIZE(dFractionBase)) .AND. (SIZE(dMolesSpecies) == SIZE(dMolesBase)) .AND. &
            (SIZE(dMolesPhase) == SIZE(dPhaseBase))
        if ((SIZE(dMolFraction) == SIZE(dFractionBase)) .AND. (SIZE(dMolesSpecies) == SIZE(dMolesBase)) .AND. &
            (SIZE(dMolesPhase) == SIZE(dPhaseBase))) then
            lAllPass = lAllPass .AND. (MAXVAL(ABS(dMolFraction-dFractionBase)) == 0D0) .AND. &
                (MAXVAL(ABS(dMolesSpecies-dMolesBase)) == 0D0) .AND. &
                (MAXVAL(ABS(dMolesPhase-dPhaseBase)) == 0D0)
        end if
        lAllPass = lAllPass .AND. (nMQMQAHessianApplyCount == 0) .AND. &
            (nMQMQAHessianAcceptedSolveCount == 0) .AND. (.NOT. lMQMQAHessianSupportedPhaseFound)

        if (lDetailed) write(*,'(A)') 'alpha zero: exact default-off final-state and zero-metric identity passed'

    end subroutine CheckDefaultAndAlphaZero


    ! A real CuFeC calculation establishes that the active scanner reaches the strict SUBG builder.
    subroutine CheckStrictSUBGRouting(lAllPass,lDetailed)

        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lDetailed
        integer :: iFailurePhase, iFailureStatus, iStatus, nAccepted, nCharged
        logical :: lEligible, lSupported
        real(8) :: dFailureMinimumFraction
        real(8), allocatable :: dDeltaA(:,:), dDeltaB(:)

        call ResetMQMQAHessianControls
        call PrepareFreshCalculation
        cThermoFileName = DATA_DIRECTORY // 'CuFeC-Kang.dat'
        dTemperature = 1400D0
        dPressure = 1D0
        dElementMass = 0D0
        dElementMass(6) = 1D0
        dElementMass(26) = 1D0
        dElementMass(29) = 1D0
        call ParseCSDataFile(cThermoFileName)
        if (INFOThermo == 0) call Thermochimica
        lAllPass = lAllPass .AND. (INFOThermo == 0)
        if (INFOThermo /= 0) return
        allocate(dDeltaA(nElements,nElements),dDeltaB(nElements))
        call BuildActiveMQMQAGEMCorrection(nSolnPhases,dDeltaA,dDeltaB,lSupported,lEligible,nAccepted,nCharged, &
            iFailurePhase,iFailureStatus,dFailureMinimumFraction,iStatus)
        lAllPass = lAllPass .AND. lSupported .AND. lEligible .AND. &
            (iStatus == MQMQA_AGGREGATE_SUCCESS) .AND. (nAccepted >= 1)
        if (lDetailed) write(*,'(A,I0)') 'strict SUBG routing: accepted active phase corrections = ',nAccepted

    end subroutine CheckStrictSUBGRouting


    ! Install a positive interior FeTiVO SUBQ state and verify the exact pre-solve correction seen by GEMNewton.
    subroutine CheckLiveSUBQIntegration(lAllPass,lDetailed)

        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lDetailed
        integer :: i, iFailurePhase, iFailureStatus, iFirst, iInfo, iLast, iNewtonInfo
        integer :: iBaselineInfo, iCorrectedInfo, iPhase, iSlot, iStatus, nAccepted, nCharged, nQuad, nVar
        integer, allocatable :: iBaselinePivot(:), iCorrectedPivot(:)
        logical :: lEligible, lRevertSave, lSupported, lUpdateDiagnosticsValid
        real(8) :: dBaselineUpdateNorm, dConflictUpdateError, dCorrectedReplayError
        real(8) :: dCorrectedUpdateNorm, dElementError, dFailureMinimumFraction
        real(8) :: dN, dOutsideA, dOutsideB, dStateDifference
        real(8) :: dRawUpdateCosine, dUpdateAngleDegrees, dUpdateCosine, dUpdateDifference, dUpdateNormRatio
        real(8), allocatable :: dBaselineSolveA(:,:), dBaselineUpdate(:)
        real(8), allocatable :: dChemicalSave(:), dDeltaA(:,:), dDeltaB(:), dFractionSave(:)
        real(8), allocatable :: dCorrectedSolveA(:,:), dCorrectedUpdate(:)
        real(8), allocatable :: dEffStoichSave(:,:), dMolesSave(:), dPhaseSave(:)
        real(8), allocatable :: dUpdateRKMP(:), dUpdateSave(:), dX(:)

        call SetMQMQAHessianControls(.TRUE.,0D0,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0)
        call PrepareFeTiVO
        call Thermochimica
        lAllPass = lAllPass .AND. (INFOThermo == 0)
        if (INFOThermo /= 0) return

        iPhase = 0
        iSlot = 0
        do i = 1, nSolnPhases
            if (cSolnPhaseType(-iAssemblage(nElements-i+1)) == 'SUBQ') then
                iPhase = -iAssemblage(nElements-i+1)
                iSlot = nElements-i+1
                exit
            end if
        end do
        lAllPass = lAllPass .AND. (iPhase > 0) .AND. (iSlot > 0)
        if (iPhase <= 0) return
        iFirst = nSpeciesPhase(iPhase-1)+1
        iLast = nSpeciesPhase(iPhase)
        nQuad = iLast-iFirst+1
        dN = dMolesPhase(iSlot)
        allocate(dX(nQuad))
        dX = dMolFraction(iFirst:iLast)
        dX = 0.8D0*dX/SUM(dX)+0.2D0/DFLOAT(nQuad)
        dMolesSpecies(iFirst:iLast) = dN*dX
        call CompChemicalPotential(.FALSE.)
        lAllPass = lAllPass .AND. (INFOThermo == 0) .AND. (MINVAL(dMolFraction(iFirst:iLast)) > 1D-12)
        if (INFOThermo /= 0) return

        allocate(dDeltaA(nElements,nElements),dDeltaB(nElements), &
            dChemicalSave(SIZE(dChemicalPotential)),dFractionSave(SIZE(dMolFraction)), &
            dMolesSave(SIZE(dMolesSpecies)),dPhaseSave(SIZE(dMolesPhase)), &
            dEffStoichSave(SIZE(dEffStoichSolnPhase,1),SIZE(dEffStoichSolnPhase,2)), &
            dUpdateSave(SIZE(dUpdateVar)),dUpdateRKMP(SIZE(dUpdateVar)))
        dChemicalSave = dChemicalPotential
        dFractionSave = dMolFraction
        dMolesSave = dMolesSpecies
        dPhaseSave = dMolesPhase
        dEffStoichSave = dEffStoichSolnPhase
        dUpdateSave = dUpdateVar
        lRevertSave = lRevertSystem
        call BuildActiveMQMQAGEMCorrection(nSolnPhases,dDeltaA,dDeltaB,lSupported,lEligible,nAccepted,nCharged, &
            iFailurePhase,iFailureStatus,dFailureMinimumFraction,iStatus)
        dStateDifference = DMAX1(MAXVAL(ABS(dChemicalPotential-dChemicalSave)), &
            MAXVAL(ABS(dMolFraction-dFractionSave)),MAXVAL(ABS(dMolesSpecies-dMolesSave)), &
            MAXVAL(ABS(dMolesPhase-dPhaseSave)))
        lAllPass = lAllPass .AND. lSupported .AND. lEligible .AND. (iStatus == MQMQA_AGGREGATE_SUCCESS) .AND. &
            (nAccepted >= 1) .AND. (dStateDifference == 0D0)

        call ResetGEMNewtonDiagnosticCapture
        lCaptureGEMNewtonSystem = .TRUE.
        lCaptureGEMNewtonCorrectedSystem = .TRUE.
        lUseMQMQAExactHessian = .TRUE.
        dMQMQAHessianAlpha = 1D0
        lUseRKMPExactHessian = .TRUE.
        lRKMPHessianActive = .FALSE.
        call GEMNewton(iNewtonInfo)
        lCaptureGEMNewtonSystem = .FALSE.
        lCaptureGEMNewtonCorrectedSystem = .FALSE.
        lAllPass = lAllPass .AND. (iNewtonInfo == 0) .AND. lGEMNewtonSystemCaptured .AND. &
            lGEMNewtonCorrectedSystemCaptured .AND. lMQMQAHessianCorrectionApplied .AND. &
            lMQMQAHessianCorrectedSolveAccepted .AND. (nMQMQAHessianApplyCount > 0) .AND. &
            (nMQMQAHessianAcceptedSolveCount > 0) .AND. (nMQMQAHessianRKMPConflictCount == 0) .AND. &
            ALL(IEEE_IS_FINITE(dUpdateVar))
        if ((.NOT. lGEMNewtonSystemCaptured) .OR. (.NOT. lGEMNewtonCorrectedSystemCaptured)) return

        nVar = SIZE(dCapturedGEMNewtonA,1)
        dElementError = DMAX1(MAXVAL(ABS((dCapturedGEMNewtonCorrectedA(1:nElements,1:nElements)- &
            dCapturedGEMNewtonA(1:nElements,1:nElements))-dDeltaA)), &
            MAXVAL(ABS((dCapturedGEMNewtonCorrectedB(1:nElements)- &
            dCapturedGEMNewtonB(1:nElements))-dDeltaB)))
        dOutsideA = 0D0
        dOutsideB = 0D0
        if (nVar > nElements) then
            dOutsideA = DMAX1(MAXVAL(ABS(dCapturedGEMNewtonCorrectedA(nElements+1:nVar,:)- &
                dCapturedGEMNewtonA(nElements+1:nVar,:))), &
                MAXVAL(ABS(dCapturedGEMNewtonCorrectedA(:,nElements+1:nVar)- &
                dCapturedGEMNewtonA(:,nElements+1:nVar))))
            dOutsideB = MAXVAL(ABS(dCapturedGEMNewtonCorrectedB(nElements+1:nVar)- &
                dCapturedGEMNewtonB(nElements+1:nVar)))
        end if
        lAllPass = lAllPass .AND. (dElementError <= 1D-13) .AND. &
            (dOutsideA == 0D0) .AND. (dOutsideB == 0D0)

        ! Replay both captured systems independently. This checks the behavioral consequence of the correction:
        ! it must not merely appear in A/B, but must produce a numerically resolved change in the solved Newton
        ! direction while reproducing the update accepted by the live corrected path.
        allocate(dBaselineSolveA(nVar,nVar),dBaselineUpdate(nVar),iBaselinePivot(nVar), &
            dCorrectedSolveA(nVar,nVar),dCorrectedUpdate(nVar),iCorrectedPivot(nVar))
        dBaselineSolveA = dCapturedGEMNewtonA
        dBaselineUpdate = dCapturedGEMNewtonB
        dCorrectedSolveA = dCapturedGEMNewtonCorrectedA
        dCorrectedUpdate = dCapturedGEMNewtonCorrectedB
        call DGESV(nVar,1,dBaselineSolveA,nVar,iBaselinePivot,dBaselineUpdate,nVar,iBaselineInfo)
        call DGESV(nVar,1,dCorrectedSolveA,nVar,iCorrectedPivot,dCorrectedUpdate,nVar,iCorrectedInfo)
        lUpdateDiagnosticsValid = .FALSE.
        dBaselineUpdateNorm = HUGE(1D0)
        dCorrectedUpdateNorm = HUGE(1D0)
        dUpdateDifference = HUGE(1D0)
        dUpdateNormRatio = HUGE(1D0)
        dUpdateCosine = HUGE(1D0)
        dUpdateAngleDegrees = HUGE(1D0)
        dCorrectedReplayError = HUGE(1D0)
        if ((iBaselineInfo == 0) .AND. (iCorrectedInfo == 0) .AND. &
            ALL(IEEE_IS_FINITE(dBaselineUpdate)) .AND. ALL(IEEE_IS_FINITE(dCorrectedUpdate))) then
            dBaselineUpdateNorm = SQRT(SUM(dBaselineUpdate**2))
            dCorrectedUpdateNorm = SQRT(SUM(dCorrectedUpdate**2))
            if (IEEE_IS_FINITE(dBaselineUpdateNorm) .AND. IEEE_IS_FINITE(dCorrectedUpdateNorm) .AND. &
                (dBaselineUpdateNorm > 0D0) .AND. (dCorrectedUpdateNorm > 0D0)) then
                dRawUpdateCosine = SUM(dCorrectedUpdate*dBaselineUpdate)/ &
                    (dCorrectedUpdateNorm*dBaselineUpdateNorm)
                if (IEEE_IS_FINITE(dRawUpdateCosine)) then
                    dUpdateDifference = SQRT(SUM((dCorrectedUpdate-dBaselineUpdate)**2))/ &
                        DMAX1(dCorrectedUpdateNorm,dBaselineUpdateNorm)
                    dUpdateNormRatio = dCorrectedUpdateNorm/dBaselineUpdateNorm
                    dUpdateCosine = DMAX1(-1D0,DMIN1(1D0,dRawUpdateCosine))
                    dUpdateAngleDegrees = DACOS(dUpdateCosine)*180D0/DACOS(-1D0)
                    dCorrectedReplayError = SQRT(SUM((dCorrectedUpdate-dUpdateVar(1:nVar))**2))/ &
                        DMAX1(1D0,dCorrectedUpdateNorm,SQRT(SUM(dUpdateVar(1:nVar)**2)))
                    lUpdateDiagnosticsValid = IEEE_IS_FINITE(dUpdateDifference) .AND. &
                        IEEE_IS_FINITE(dUpdateNormRatio) .AND. IEEE_IS_FINITE(dUpdateAngleDegrees) .AND. &
                        IEEE_IS_FINITE(dCorrectedReplayError)
                end if
            end if
        end if
        lAllPass = lAllPass .AND. lUpdateDiagnosticsValid .AND. (dUpdateDifference > 1D-8) .AND. &
            (dCorrectedReplayError <= 1D-13)

        ! First solve the identical restored state with RKMP ownership alone. This supplies the update oracle for
        ! the controlled both-eligible conflict without fabricating a multiphase assessed database case.
        dChemicalPotential = dChemicalSave
        dMolFraction = dFractionSave
        dMolesSpecies = dMolesSave
        dMolesPhase = dPhaseSave
        dEffStoichSolnPhase = dEffStoichSave
        dUpdateVar = dUpdateSave
        lRevertSystem = lRevertSave
        lUseMQMQAExactHessian = .FALSE.
        lUseRKMPExactHessian = .TRUE.
        lRKMPHessianActive = .TRUE.
        lRKMPHessianNonlinearReady = .FALSE.
        call GEMNewton(iNewtonInfo)
        lAllPass = lAllPass .AND. (iNewtonInfo == 0) .AND. ALL(IEEE_IS_FINITE(dUpdateVar))
        dUpdateRKMP = dUpdateVar

        ! Restore again, request both paths, and require the conflict update to equal the RKMP-owned update.
        dChemicalPotential = dChemicalSave
        dMolFraction = dFractionSave
        dMolesSpecies = dMolesSave
        dMolesPhase = dPhaseSave
        dEffStoichSolnPhase = dEffStoichSave
        dUpdateVar = dUpdateSave
        lRevertSystem = lRevertSave
        nMQMQAHessianRKMPConflictCount = 0
        nMQMQAHessianApplyCount = 0
        lMQMQAHessianFallbackUsed = .FALSE.
        lUseMQMQAExactHessian = .TRUE.
        dMQMQAHessianAlpha = 1D0
        lRKMPHessianActive = .TRUE.
        lUseRKMPExactHessian = .TRUE.
        lRKMPHessianNonlinearReady = .FALSE.
        call ResetGEMNewtonDiagnosticCapture
        lCaptureGEMNewtonCorrectedSystem = .TRUE.
        call GEMNewton(iNewtonInfo)
        lCaptureGEMNewtonCorrectedSystem = .FALSE.
        dConflictUpdateError = SQRT(SUM((dUpdateVar-dUpdateRKMP)**2))/ &
            DMAX1(1D0,SQRT(SUM(dUpdateRKMP**2)))
        lAllPass = lAllPass .AND. (iNewtonInfo == 0) .AND. (nMQMQAHessianRKMPConflictCount == 1) .AND. &
            (nMQMQAHessianApplyCount == 0) .AND. lMQMQAHessianFallbackUsed .AND. &
            (.NOT. lGEMNewtonCorrectedSystemCaptured) .AND. (dConflictUpdateError <= 1D-14)
        lRKMPHessianActive = .FALSE.
        lUseRKMPExactHessian = .FALSE.

        if (lDetailed) then
            write(*,'(A,ES12.4)') 'live SUBQ captured application error = ',dElementError
            write(*,'(A,2ES12.4)') 'maximum correction ratios rhoA/rhoB = ', &
                dMQMQAHessianMaxRatioA,dMQMQAHessianMaxRatioB
            write(*,'(A,ES12.4)') 'pre-line-search GEM solution-vector separation = ',dUpdateDifference
            write(*,'(A,2ES12.4)') 'historical/corrected GEM solution-vector norms = ', &
                dBaselineUpdateNorm,dCorrectedUpdateNorm
            write(*,'(A,ES12.4)') 'corrected-to-historical solution norm ratio = ',dUpdateNormRatio
            write(*,'(A,ES12.4)') 'GEM solution-vector cosine similarity = ',dUpdateCosine
            write(*,'(A,ES12.4)') 'unscaled-coordinate solution-vector angle (degrees) = ',dUpdateAngleDegrees
            write(*,'(A,ES12.4)') 'corrected solution-vector replay error = ',dCorrectedReplayError
            write(*,'(A)') 'controlled RKMP/MQMQA ownership conflict preserved the untouched RKMP path'
            write(*,'(A,ES12.4)') 'conflict versus equivalent RKMP-owned update error = ',dConflictUpdateError
        end if

    end subroutine CheckLiveSUBQIntegration


    ! Complete the public-control routing matrix with real RKMP-only and neither-eligible calculations.
    subroutine CheckRemainingDualControlCases(lAllPass,lDetailed)

        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lDetailed
        integer :: iInfo
        logical :: lNeitherPass, lRKMPOnlyPass

        call SetRKMPHessianControls(.TRUE.,1D0,.FALSE.,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0)
        call SetMQMQAHessianControls(.TRUE.,1D0,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0)

        call PrepareFreshCalculation
        cThermoFileName = DATA_DIRECTORY // 'WAuArO-1.dat'
        dPressure = 1D0
        dTemperature = 1455D0
        dElementMass = 0D0
        dElementMass(74) = 1.95D0
        dElementMass(79) = 1D0
        dElementMass(18) = 2D0
        dElementMass(8) = 10D0
        call ParseCSDataFile(cThermoFileName)
        if (INFOThermo == 0) call Thermochimica
        lRKMPOnlyPass = (INFOThermo == 0) .AND. lRKMPHessianWasActive .AND. &
            (nRKMPHessianApplyCount > 0) .AND. (.NOT. lMQMQAHessianSupportedPhaseFound) .AND. &
            (nMQMQAHessianApplyCount == 0) .AND. (nMQMQAHessianRKMPConflictCount == 0)
        lAllPass = lAllPass .AND. lRKMPOnlyPass

        call PrepareFreshCalculation
        cThermoFileName = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'
        dPressure = 1D0
        dTemperature = 2250D0
        dElementMass = 0D0
        dElementMass(42) = 0.8D0
        dElementMass(44) = 0.2D0
        call ParseCSDataFile(cThermoFileName)
        if (INFOThermo == 0) call Thermochimica
        lNeitherPass = (INFOThermo == 0) .AND. (.NOT. lRKMPHessianWasActive) .AND. &
            (.NOT. lMQMQAHessianSupportedPhaseFound) .AND. (nMQMQAHessianApplyCount == 0) .AND. &
            (nMQMQAHessianRKMPConflictCount == 0)
        lAllPass = lAllPass .AND. lNeitherPass

        if (lDetailed) write(*,'(A,L1,A,L1)') &
            'dual-control real cases: only-RKMP=',lRKMPOnlyPass,', neither-eligible=',lNeitherPass

    end subroutine CheckRemainingDualControlCases


    ! Full fixed-alpha-one behavior is recorded for MQ-4D; convergence is deliberately not an MQ-4C gate.
    subroutine ReportFullAlphaOneEvidence(lAllPass)

        logical, intent(inout) :: lAllPass
        integer :: i, iInfo
        real(8) :: dFractionDifference, dGibbsBase, dGibbsDifference
        real(8) :: dMolesDifference, dPhaseDifference
        real(8), allocatable :: dFractionBase(:), dMolesBase(:), dPhaseBase(:)

        call ResetMQMQAHessianControls
        call PrepareFeTiVO
        call Thermochimica
        lAllPass = lAllPass .AND. (INFOThermo == 0)
        if (INFOThermo /= 0) return
        dGibbsBase = dGibbsEnergySys
        allocate(dFractionBase(SIZE(dMolFraction)),dMolesBase(SIZE(dMolesSpecies)), &
            dPhaseBase(SIZE(dMolesPhase)))
        dFractionBase = dMolFraction
        dMolesBase = dMolesSpecies
        dPhaseBase = dMolesPhase

        call SetMQMQAHessianControls(.TRUE.,1D0,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0)
        call PrepareFeTiVO
        call Thermochimica
        dGibbsDifference = DABS(dGibbsEnergySys-dGibbsBase)/DMAX1(1D0,DABS(dGibbsBase))
        dFractionDifference = MAXVAL(ABS(dMolFraction-dFractionBase))/ &
            DMAX1(1D0,MAXVAL(ABS(dFractionBase)))
        dMolesDifference = MAXVAL(ABS(dMolesSpecies-dMolesBase))/ &
            DMAX1(1D0,MAXVAL(ABS(dMolesBase)))
        dPhaseDifference = MAXVAL(ABS(dMolesPhase-dPhaseBase))/ &
            DMAX1(1D0,MAXVAL(ABS(dPhaseBase)))
        write(*,'(A)') 'alpha-one nonlinear evidence (not an MQ-4C acceptance gate)'
        write(*,'(A,I0)') '  INFOThermo = ',INFOThermo
        write(*,'(A,L1,A,I0,A,I0)') '  converged=',lConverged,', iterations=',iterGlobal,', reversions=',iterRevert
        write(*,'(A,I0,A,I0)') '  applications=',nMQMQAHessianApplyCount, &
            ', accepted solves=',nMQMQAHessianAcceptedSolveCount
        write(*,'(A,6(I0,1X))') '  fallbacks interior/aggregate/application/DGESV/nonfinite/conflict=', &
            nMQMQAHessianInteriorFallbackCount,nMQMQAHessianAggregateFailureCount,nMQMQAHessianApplicationFailureCount, &
            nMQMQAHessianDGESVFallbackCount,nMQMQAHessianNonfiniteFallbackCount,nMQMQAHessianRKMPConflictCount
        if (nMQMQAHessianInteriorFallbackCount > 0) write(*,'(A,ES12.4)') &
            '  minimum rejected fraction=',dMQMQAHessianMinimumRejectedFraction
        write(*,'(A,I0,A,I0)') '  last failed phase/status=',iMQMQAHessianLastFailurePhase,'/', &
            iMQMQAHessianLastFailureStatus
        write(*,'(A,2ES12.4)') '  max rhoA/rhoB=',dMQMQAHessianMaxRatioA,dMQMQAHessianMaxRatioB
        write(*,'(A,4ES12.4)') '  scaled final differences G/x/n/N=',dGibbsDifference, &
            dFractionDifference,dMolesDifference,dPhaseDifference
        write(*,'(A,I0,A,I0)') '  active solution phases=',nSolnPhases,', pure phases=',nConPhases
        do i = 1, nSolnPhases
            write(*,'(A,I0,A,A,A,ES12.4)') '    solution ',i,': ', &
                TRIM(cSolnPhaseName(-iAssemblage(nElements-i+1))),'  moles=',dMolesPhase(nElements-i+1)
        end do
        do i = 1, nConPhases
            write(*,'(A,I0,A,A,A,ES12.4)') '    pure ',i,': ',TRIM(cSpeciesName(iAssemblage(i))), &
                '  moles=',dMolesPhase(i)
        end do

    end subroutine ReportFullAlphaOneEvidence


    subroutine PrepareFeTiVO

        call PrepareFreshCalculation
        cThermoFileName = DATA_DIRECTORY // 'FeTiVO.dat'
        dTemperature = 2000D0
        dPressure = 1D0
        dElementMass = 0D0
        dElementMass(8) = 2D0
        dElementMass(22) = 0.5D0
        dElementMass(23) = 0.5D0
        dElementMass(26) = 0.5D0
        call ParseCSDataFile(cThermoFileName)

    end subroutine PrepareFeTiVO


    subroutine PrepareFreshCalculation

        if (lHaveCalculation) call ResetThermoAll
        lHaveCalculation = .TRUE.
        cInputUnitTemperature = 'K'
        cInputUnitPressure = 'atm'
        cInputUnitMass = 'moles'
        dElementMass = 0D0

    end subroutine PrepareFreshCalculation


    subroutine FinishTest(lSucceeded)

        logical, intent(in) :: lSucceeded

        call ResetGEMNewtonDiagnosticCapture
        call ResetMQMQAHessianControls
        call ResetRKMPHessianControls
        if (lHaveCalculation) call ResetThermoAll
        if (lSucceeded) then
            write(*,'(A)') 'TestMQMQAGEMIntegration: PASS'
            call EXIT(0)
        else
            write(*,'(A)') 'TestMQMQAGEMIntegration: FAIL <---'
            call EXIT(1)
        end if

    end subroutine FinishTest

end program TestMQMQAGEMIntegration
