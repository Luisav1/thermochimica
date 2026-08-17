!-------------------------------------------------------------------------------------------------------------
!> \file    TestMQMQAAdaptiveTrust.F90
!> \brief   Verify MQ-4D candidate selection, groupwise trust, controls, and sustained full-alpha behaviour.
!>
!> \details This test separates three kinds of evidence.  Small controlled checks exercise candidate ordering,
!!          invalid controls, and individual rejection mechanisms.  A live FeTiVO calculation then demonstrates
!!          that adaptive integration converges, records every reduced-alpha decision, and ends with sustained
!!          full-alpha use.  Finally, a bounded 19-case sensitivity matrix changes one numerical trust heuristic
!!          at a time and compares the resulting trajectory with one default-off historical reference.
!!
!!          The sensitivity classifications describe observed FeTiVO solver behaviour; they are not statements
!!          that a threshold is thermodynamically correct or universally optimal for all MQMQA assessments.
!-------------------------------------------------------------------------------------------------------------
program TestMQMQAAdaptiveTrust

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE, IEEE_VALUE, IEEE_QUIET_NAN, &
        IEEE_POSITIVE_INF, IEEE_NEGATIVE_INF
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleParseCS
    USE ModuleGEMSolver
    USE ModuleMQMQATrust

    implicit none

    interface
        subroutine SetMQMQAHessianControls(lEnable,dAlpha,iInfo)
            logical, intent(in) :: lEnable
            real(8), intent(in) :: dAlpha
            integer, intent(out) :: iInfo
        end subroutine SetMQMQAHessianControls
        subroutine ResetMQMQAHessianControls
        end subroutine ResetMQMQAHessianControls
        subroutine SetMQMQAHessianAdaptiveControls(lEnable,dAlphaMax,iInfo)
            logical, intent(in) :: lEnable
            real(8), intent(in) :: dAlphaMax
            integer, intent(out) :: iInfo
        end subroutine SetMQMQAHessianAdaptiveControls
        subroutine ResetMQMQAHessianAdaptiveControls
        end subroutine ResetMQMQAHessianAdaptiveControls
    end interface

    integer, parameter :: nSensitivityCases = 19
    integer, parameter :: MQMQA_SENSITIVITY_ROBUST = 1
    integer, parameter :: MQMQA_SENSITIVITY_SAFE = 2
    integer, parameter :: MQMQA_SENSITIVITY_FAILED = 3

    type :: MQMQATrustSettings
        real(8) :: dLocalThreshold
        real(8) :: dGibbsActivation
        real(8) :: dGibbsRetention
        real(8) :: dProgressAllowance
        real(8) :: dUpdateRatio
        real(8) :: dDirectionCosine
        real(8) :: dRelativeDirection
        real(8) :: dEmergencyCap
        real(8) :: dResolvedNormFloor
        integer :: iSettledPeriod
    end type MQMQATrustSettings

    type :: MQMQASensitivityResult
        character(len=36) :: cName = ''
        character(len=24) :: cChangedValue = ''
        character(len=20) :: cClassification = ''
        integer :: iInfo = -1
        integer :: iIterations = 0
        integer :: nEligible = 0
        integer :: iFirstPositive = 0
        integer :: iFirstFinalFull = 0
        integer :: nFull = 0
        integer :: nReduced = 0
        integer :: nZero = 0
        integer :: nFinalFull = 0
        integer :: nReadinessActivations = 0
        integer :: nReadinessResets = 0
        integer :: nReject(7) = 0
        real(8) :: dMaximumAlpha = 0D0
        real(8) :: dDifference(4) = HUGE(1D0)
        logical :: lConverged = .FALSE.
        logical :: lFinite = .FALSE.
        logical :: lAlphaBounds = .FALSE.
        logical :: lReducedDocumented = .FALSE.
        logical :: lFallbackSafe = .FALSE.
        logical :: lHistoricalAgreement = .FALSE.
        logical :: lNoResetInFinalWindow = .FALSE.
        logical :: lDefaultsRestored = .FALSE.
        logical :: lMateriallyDifferent = .FALSE.
        integer :: iClassification = MQMQA_SENSITIVITY_FAILED
    end type MQMQASensitivityResult

    integer :: iArgument
    logical :: lHaveCalculation, lPass, lReport, lSensitivityReport
    character(len=32) :: cArgument

    lHaveCalculation = .FALSE.
    lPass = .TRUE.
    lReport = .FALSE.
    lSensitivityReport = .FALSE.
    do iArgument = 1,COMMAND_ARGUMENT_COUNT()
        call GET_COMMAND_ARGUMENT(iArgument,cArgument)
        lReport = lReport .OR. (TRIM(cArgument) == '--report')
        lSensitivityReport = lSensitivityReport .OR. (TRIM(cArgument) == '--sensitivity-report')
    end do

    call CheckStateFreeTrust(lPass,lReport)
    call CheckAdaptiveControlContract(lPass,lReport)
    call CheckAdaptiveAlphaZero(lPass,lReport)
    call CheckAdaptiveFeTiVOSensitivity(lPass,lReport,lSensitivityReport)
    call FinishTest(lPass)

contains

    subroutine CheckStateFreeTrust(lAllPass,lDetailed)

        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lDetailed
        integer :: iStatus, nAlpha
        logical :: lAccepted
        real(8) :: dAlpha(5), dBaseNorm(3), dCosine(3), dDifference(3), dRatio(3), dRatioA, dRatioB
        real(8) :: dTrialNorm(3)
        real(8) :: dABase(2,2), dBBase(2), dDeltaA(2,2), dDeltaB(2)
        real(8) :: dStepBase(4), dStepTrial(4)

        call BuildMQMQAAlphaCandidateList(1D0,dAlpha,nAlpha)
        lAllPass = lAllPass .AND. (nAlpha == 5) .AND. &
            (MAXVAL(ABS(dAlpha-[1D0,1D-1,1D-2,1D-3,0D0])) == 0D0)
        call BuildMQMQAAlphaCandidateList(0.35D0,dAlpha,nAlpha)
        lAllPass = lAllPass .AND. (nAlpha == 5) .AND. &
            (MAXVAL(ABS(dAlpha-[0.35D0,1D-1,1D-2,1D-3,0D0])) == 0D0)
        call BuildMQMQAAlphaCandidateList(0D0,dAlpha,nAlpha)
        lAllPass = lAllPass .AND. (nAlpha == 1) .AND. (dAlpha(1) == 0D0)

        dABase = 0D0
        dABase(1,1) = 1D0
        dABase(2,2) = 1D0
        dBBase = [1D0,1D0]
        dDeltaA = 0D0
        dDeltaA(1,1) = 0.1D0
        dDeltaB = [0.1D0,0D0]
        call EvaluateMQMQACorrectionRatio(dABase,dBBase,2,dDeltaA,dDeltaB,1D0,1D6, &
            lAccepted,dRatioA,dRatioB)
        lAllPass = lAllPass .AND. lAccepted .AND. (dRatioA < 1D0) .AND. (dRatioB < 1D0)
        call EvaluateMQMQACorrectionRatio(dABase,dBBase,2,2D7*dDeltaA,2D7*dDeltaB,1D0,1D6, &
            lAccepted,dRatioA,dRatioB)
        lAllPass = lAllPass .AND. (.NOT. lAccepted) .AND. (DMAX1(dRatioA,dRatioB) > 1D6)

        dStepBase = [1D0,2D0,0.2D0,3D0]
        dStepTrial = dStepBase
        call EvaluateMQMQAUpdateTrust(dStepBase,dStepTrial,2,1,1,1.25D0,0.90D0,0.50D0, &
            lAccepted,dRatio,dCosine,dDifference,dBaseNorm,dTrialNorm,iStatus)
        lAllPass = lAllPass .AND. lAccepted .AND. (iStatus == MQMQA_TRUST_ACCEPTED) .AND. &
            (MAXVAL(ABS(dRatio-1D0)) == 0D0) .AND. (MAXVAL(ABS(dCosine-1D0)) <= 1D-15)

        dStepTrial = dStepBase
        dStepTrial(1:2) = 2D0*dStepBase(1:2)
        call EvaluateMQMQAUpdateTrust(dStepBase,dStepTrial,2,1,1,1.25D0,0.90D0,0.50D0, &
            lAccepted,dRatio,dCosine,dDifference,dBaseNorm,dTrialNorm,iStatus)
        lAllPass = lAllPass .AND. (.NOT. lAccepted) .AND. (iStatus == MQMQA_TRUST_UPDATE_REJECTED) .AND. &
            (dRatio(1) > 1.25D0)

        dStepTrial = dStepBase
        dStepTrial(1:2) = -dStepBase(1:2)
        call EvaluateMQMQAUpdateTrust(dStepBase,dStepTrial,2,1,1,2D0,0.90D0,2.5D0, &
            lAccepted,dRatio,dCosine,dDifference,dBaseNorm,dTrialNorm,iStatus)
        lAllPass = lAllPass .AND. (.NOT. lAccepted) .AND. &
            (iStatus == MQMQA_TRUST_DIRECTION_REJECTED) .AND. (dCosine(1) < 0D0)

        dStepBase = 0D0
        dStepTrial = 0D0
        dStepBase(1) = 1D-10
        dStepTrial(1) = 9D-9
        call EvaluateMQMQAUpdateTrust(dStepBase,dStepTrial,2,1,1,1.25D0,0.90D0,0.50D0, &
            lAccepted,dRatio,dCosine,dDifference,dBaseNorm,dTrialNorm,iStatus)
        lAllPass = lAllPass .AND. lAccepted .AND. (iStatus == MQMQA_TRUST_ACCEPTED) .AND. &
            (dRatio(1) == 1D0) .AND. (ABS(dBaseNorm(1)-1D-10) <= 1D-24) .AND. &
            (ABS(dTrialNorm(1)-9D-9) <= 1D-22)

        if (lDetailed) write(*,'(A)') &
            'state-free trust: candidates, emergency ratio, group update, and direction gates passed'

    end subroutine CheckStateFreeTrust


    subroutine CheckAdaptiveControlContract(lAllPass,lDetailed)

        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lDetailed
        integer :: i, iInfo
        real(8) :: dInvalid(5)

        call ResetMQMQAHessianControls
        call ResetMQMQAHessianAdaptiveControls
        call SetMQMQAHessianAdaptiveControls(.TRUE.,0.35D0,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0) .AND. lMQMQAHessianAdaptiveControlsConfigured .AND. &
            lMQMQAHessianRequestedAdaptiveEnable .AND. (dMQMQAHessianRequestedAlphaMax == 0.35D0) .AND. &
            (.NOT. lMQMQAHessianControlsConfigured) .AND. (.NOT. lMQMQAHessianRequestedEnable)

        dInvalid(1:2) = [-0.1D0,1.1D0]
        dInvalid(3) = IEEE_VALUE(0D0,IEEE_QUIET_NAN)
        dInvalid(4) = IEEE_VALUE(0D0,IEEE_POSITIVE_INF)
        dInvalid(5) = IEEE_VALUE(0D0,IEEE_NEGATIVE_INF)
        do i = 1,SIZE(dInvalid)
            call SetMQMQAHessianAdaptiveControls(.FALSE.,dInvalid(i),iInfo)
            lAllPass = lAllPass .AND. (iInfo /= 0) .AND. lMQMQAHessianRequestedAdaptiveEnable .AND. &
                (dMQMQAHessianRequestedAlphaMax == 0.35D0)
        end do

        call SetMQMQAHessianControls(.TRUE.,0.2D0,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0) .AND. lMQMQAHessianControlsConfigured .AND. &
            lMQMQAHessianRequestedEnable .AND. (.NOT. lMQMQAHessianAdaptiveControlsConfigured) .AND. &
            (.NOT. lMQMQAHessianRequestedAdaptiveEnable)
        call SetMQMQAHessianAdaptiveControls(.TRUE.,1D0,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0) .AND. lMQMQAHessianAdaptiveControlsConfigured .AND. &
            (.NOT. lMQMQAHessianControlsConfigured) .AND. (.NOT. lMQMQAHessianRequestedEnable)
        call ResetMQMQAHessianAdaptiveControls

        if (lDetailed) write(*,'(A)') 'adaptive controls: persistence, atomic invalid calls, and mode ownership passed'

    end subroutine CheckAdaptiveControlContract


    subroutine CheckAdaptiveAlphaZero(lAllPass,lDetailed)

        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lDetailed
        integer :: iInfo, iInfoBase
        real(8) :: dGibbsBase
        real(8), allocatable :: dFractionBase(:), dMolesBase(:), dPhaseBase(:)

        call ResetMQMQAHessianControls
        call ResetMQMQAHessianAdaptiveControls
        call PrepareFeTiVO
        call Thermochimica
        iInfoBase = INFOThermo
        lAllPass = lAllPass .AND. (iInfoBase == 0)
        if (iInfoBase /= 0) return
        dGibbsBase = dGibbsEnergySys
        allocate(dFractionBase(SIZE(dMolFraction)),dMolesBase(SIZE(dMolesSpecies)),dPhaseBase(SIZE(dMolesPhase)))
        dFractionBase = dMolFraction
        dMolesBase = dMolesSpecies
        dPhaseBase = dMolesPhase

        call SetMQMQAHessianAdaptiveControls(.TRUE.,0D0,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0)
        call PrepareFeTiVO
        call Thermochimica
        lAllPass = lAllPass .AND. (INFOThermo == iInfoBase) .AND. lUseMQMQAExactHessian .AND. &
            lMQMQAHessianAdaptiveMode .AND. (dMQMQAHessianAlpha == 0D0) .AND. &
            (dGibbsEnergySys == dGibbsBase) .AND. (MAXVAL(ABS(dMolFraction-dFractionBase)) == 0D0) .AND. &
            (MAXVAL(ABS(dMolesSpecies-dMolesBase)) == 0D0) .AND. &
            (MAXVAL(ABS(dMolesPhase-dPhaseBase)) == 0D0) .AND. &
            (nMQMQAHessianEligibleSolveCount == 0) .AND. (nMQMQAHessianApplyCount == 0)

        if (lDetailed) write(*,'(A)') 'adaptive alpha zero: exact historical final-state identity passed'

    end subroutine CheckAdaptiveAlphaZero


    subroutine CheckAdaptiveFeTiVOSensitivity(lAllPass,lDetailed,lSensitivityDetailed)

        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lDetailed, lSensitivityDetailed
        integer :: iCase
        real(8) :: dGibbsBase
        real(8), allocatable :: dFractionBase(:), dMolesBase(:), dPhaseBase(:)
        type(MQMQATrustSettings) :: tCaseSettings, tProductionSettings
        type(MQMQASensitivityResult) :: tDefaultResult
        type(MQMQASensitivityResult) :: tResult(nSensitivityCases)

        call CaptureTrustSettings(tProductionSettings)
        lAllPass = lAllPass .AND. ProductionTrustDefaultsMatch(tProductionSettings)

        ! Establish one historical final state with MQMQA correction integration disabled.  Every adaptive case
        ! below is compared with this same reference so threshold effects are not confused with a changing oracle.
        call ResetMQMQAHessianControls
        call ResetMQMQAHessianAdaptiveControls
        call ApplyTrustSettings(tProductionSettings)
        call PrepareFeTiVO
        call Thermochimica
        lAllPass = lAllPass .AND. (INFOThermo == 0) .AND. lConverged
        if ((INFOThermo /= 0) .OR. (.NOT. lConverged)) then
            call ApplyTrustSettings(tProductionSettings)
            return
        end if
        dGibbsBase = dGibbsEnergySys
        allocate(dFractionBase(SIZE(dMolFraction)),dMolesBase(SIZE(dMolesSpecies)),dPhaseBase(SIZE(dMolesPhase)))
        dFractionBase = dMolFraction
        dMolesBase = dMolesSpecies
        dPhaseBase = dMolesPhase
        tDefaultResult = MQMQASensitivityResult()

        ! Change one trust setting at a time, except for the two explicitly labelled combined activation-policy
        ! cases.  Each run restores the production settings, records its complete alpha/rejection history, and is
        ! classified by convergence safety, historical agreement, and material trajectory changes.
        do iCase = 1,nSensitivityCases
            call ConfigureSensitivityCase(iCase,tProductionSettings,tCaseSettings, &
                tResult(iCase)%cName,tResult(iCase)%cChangedValue)
            call RunFeTiVOSensitivityCase(tCaseSettings,tProductionSettings,dGibbsBase, &
                dFractionBase,dMolesBase,dPhaseBase,tResult(iCase))
            call ClassifySensitivityResult(tResult(iCase),tDefaultResult,iCase)
            if (iCase == 1) tDefaultResult = tResult(iCase)

            lAllPass = lAllPass .AND. tResult(iCase)%lFinite .AND. tResult(iCase)%lAlphaBounds .AND. &
                tResult(iCase)%lReducedDocumented .AND. tResult(iCase)%lFallbackSafe .AND. &
                tResult(iCase)%lDefaultsRestored .AND. &
                (tResult(iCase)%iClassification /= MQMQA_SENSITIVITY_FAILED)

            if (iCase == 1) then
                lAllPass = lAllPass .AND. (tResult(iCase)%iInfo == 0) .AND. &
                    tResult(iCase)%lConverged .AND. (tResult(iCase)%nEligible > 0) .AND. &
                    (tResult(iCase)%iFirstPositive > 0) .AND. &
                    (tResult(iCase)%nFinalFull >= 3) .AND. &
                    tResult(iCase)%lNoResetInFinalWindow .AND. tResult(iCase)%lHistoricalAgreement
            end if

            if ((lDetailed .AND. (iCase == 1)) .OR. &
                (lSensitivityDetailed .AND. &
                ((tResult(iCase)%iClassification == MQMQA_SENSITIVITY_FAILED) .OR. &
                tResult(iCase)%lMateriallyDifferent))) then
                call PrintCurrentCandidateHistory(tResult(iCase)%cName)
            end if
        end do

        call ApplyTrustSettings(tProductionSettings)
        lAllPass = lAllPass .AND. TrustSettingsMatch(tProductionSettings)
        if (lDetailed) call PrintAdaptiveDefaultResult(tResult(1))
        if (lSensitivityDetailed) call PrintSensitivityReport(tResult)

    end subroutine CheckAdaptiveFeTiVOSensitivity


    subroutine RunFeTiVOSensitivityCase(tSettings,tProductionSettings,dGibbsBase, &
        dFractionBase,dMolesBase,dPhaseBase,tResult)

        type(MQMQATrustSettings), intent(in) :: tSettings, tProductionSettings
        real(8), intent(in) :: dGibbsBase, dFractionBase(:), dMolesBase(:), dPhaseBase(:)
        type(MQMQASensitivityResult), intent(inout) :: tResult

        integer :: i, iCandidate, iFirstFinalSlot, iInfo, nHardFallbacks
        logical :: lDimensionsMatch

        call ApplyTrustSettings(tSettings)
        call ResetMQMQAHessianControls
        call ResetMQMQAHessianAdaptiveControls
        call SetMQMQAHessianAdaptiveControls(.TRUE.,1D0,iInfo)
        if (iInfo == 0) then
            call PrepareFeTiVO
            call Thermochimica
        else
            INFOThermo = iInfo
        end if

        tResult%iInfo = INFOThermo
        tResult%lConverged = lConverged
        tResult%iIterations = iterGlobal
        tResult%nEligible = nMQMQAHessianEligibleSolveCount
        tResult%nFull = nMQMQAHessianFullAlphaCount
        tResult%nReduced = nMQMQAHessianReducedAlphaCount
        tResult%nZero = nMQMQAHessianZeroAlphaCount
        tResult%nFinalFull = nMQMQAHessianFinalFullAlphaWindow
        tResult%nReadinessActivations = nMQMQAHessianReadinessActivationCount
        tResult%nReadinessResets = nMQMQAHessianReadinessResetCount
        tResult%nReject = [nMQMQAHessianRejectNonlinear,nMQMQAHessianRejectCorrection, &
            nMQMQAHessianRejectRatio,nMQMQAHessianRejectDGESV,nMQMQAHessianRejectNonfinite, &
            nMQMQAHessianRejectUpdate,nMQMQAHessianRejectDirection]
        tResult%dMaximumAlpha = dMQMQAHessianMaxSelectedAlpha

        tResult%lAlphaBounds = .TRUE.
        tResult%lReducedDocumented = .TRUE.
        tResult%iFirstPositive = 0
        do i = 1,MIN(nMQMQAHessianEligibleSolveCount,iterGlobalMax)
            if (dMQMQAHessianAcceptedAlphaHistory(i) < 0D0) cycle
            tResult%lAlphaBounds = tResult%lAlphaBounds .AND. &
                (dMQMQAHessianAcceptedAlphaHistory(i) >= 0D0) .AND. &
                (dMQMQAHessianAcceptedAlphaHistory(i) <= 1D0+1D-14)
            if ((tResult%iFirstPositive == 0) .AND. &
                (dMQMQAHessianAcceptedAlphaHistory(i) > 0D0)) &
                tResult%iFirstPositive = iMQMQAHessianGlobalIterationHistory(i)
            if ((dMQMQAHessianAcceptedAlphaHistory(i) > 0D0) .AND. &
                (dMQMQAHessianAcceptedAlphaHistory(i) < 1D0-1D-12)) then
                do iCandidate = 1,5
                    if (dMQMQAHessianCandidateAlphaHistory(i,iCandidate) <= &
                        dMQMQAHessianAcceptedAlphaHistory(i)+1D-14) cycle
                    tResult%lReducedDocumented = tResult%lReducedDocumented .AND. &
                        (iMQMQAHessianCandidateRejectionMaskHistory(i,iCandidate) /= 0)
                end do
            end if
        end do

        tResult%iFirstFinalFull = 0
        tResult%lNoResetInFinalWindow = .FALSE.
        if ((tResult%nFinalFull > 0) .AND. (tResult%nFinalFull <= tResult%nEligible)) then
            iFirstFinalSlot = tResult%nEligible-tResult%nFinalFull+1
            tResult%iFirstFinalFull = iMQMQAHessianGlobalIterationHistory(iFirstFinalSlot)
            tResult%lNoResetInFinalWindow = ALL( &
                iMQMQAHessianReadinessReasonHistory(iFirstFinalSlot:tResult%nEligible) == 0) .AND. ALL( &
                dMQMQAHessianAcceptedAlphaHistory(iFirstFinalSlot:tResult%nEligible) >= 1D0-1D-12)
        end if

        lDimensionsMatch = ALLOCATED(dMolFraction) .AND. ALLOCATED(dMolesSpecies) .AND. ALLOCATED(dMolesPhase)
        if (lDimensionsMatch) lDimensionsMatch = (SIZE(dMolFraction) == SIZE(dFractionBase)) .AND. &
            (SIZE(dMolesSpecies) == SIZE(dMolesBase)) .AND. (SIZE(dMolesPhase) == SIZE(dPhaseBase))
        tResult%lFinite = lDimensionsMatch .AND. IEEE_IS_FINITE(dGibbsEnergySys)
        if (lDimensionsMatch) tResult%lFinite = tResult%lFinite .AND. ALL(IEEE_IS_FINITE(dMolFraction)) .AND. &
            ALL(IEEE_IS_FINITE(dMolesSpecies)) .AND. ALL(IEEE_IS_FINITE(dMolesPhase))

        if (lDimensionsMatch .AND. tResult%lFinite) then
            tResult%dDifference(1) = ABS(dGibbsEnergySys-dGibbsBase)/DMAX1(1D0,ABS(dGibbsBase))
            tResult%dDifference(2) = MAXVAL(ABS(dMolFraction-dFractionBase))/ &
                DMAX1(1D0,MAXVAL(ABS(dFractionBase)))
            tResult%dDifference(3) = MAXVAL(ABS(dMolesSpecies-dMolesBase))/ &
                DMAX1(1D0,MAXVAL(ABS(dMolesBase)))
            tResult%dDifference(4) = MAXVAL(ABS(dMolesPhase-dPhaseBase))/ &
                DMAX1(1D0,MAXVAL(ABS(dPhaseBase)))
        end if
        tResult%lHistoricalAgreement = ALL(tResult%dDifference <= 1D-8)

        nHardFallbacks = nMQMQAHessianAggregateFailureCount+nMQMQAHessianApplicationFailureCount+ &
            nMQMQAHessianDGESVFallbackCount+nMQMQAHessianNonfiniteFallbackCount
        tResult%lFallbackSafe = (nHardFallbacks == 0) .OR. lMQMQAHessianFallbackUsed

        call ResetMQMQAHessianAdaptiveControls
        call ApplyTrustSettings(tProductionSettings)
        tResult%lDefaultsRestored = TrustSettingsMatch(tProductionSettings)

    end subroutine RunFeTiVOSensitivityCase


    subroutine ConfigureSensitivityCase(iCase,tProduction,tSettings,cName,cChangedValue)

        integer, intent(in) :: iCase
        type(MQMQATrustSettings), intent(in) :: tProduction
        type(MQMQATrustSettings), intent(out) :: tSettings
        character(len=*), intent(out) :: cName, cChangedValue

        tSettings = tProduction
        cName = ''
        cChangedValue = ''
        select case(iCase)
        case(1)
            cName = 'Production default'; cChangedValue = 'all defaults'
        case(2)
            cName = 'Local residual lower'; cChangedValue = '0.025'; tSettings%dLocalThreshold = 0.025D0
        case(3)
            cName = 'Local residual upper'; cChangedValue = '0.10'; tSettings%dLocalThreshold = 0.10D0
        case(4)
            cName = 'Settled period conservative'; cChangedValue = '8'; tSettings%iSettledPeriod = 8
        case(5)
            cName = 'Settled period permissive'; cChangedValue = '3'; tSettings%iSettledPeriod = 3
        case(6)
            cName = 'Gibbs activation lower'; cChangedValue = '1E-7'; tSettings%dGibbsActivation = 1D-7
        case(7)
            cName = 'Gibbs activation upper'; cChangedValue = '1E-5'; tSettings%dGibbsActivation = 1D-5
        case(8)
            cName = 'Gibbs retention lower'; cChangedValue = '1E-5'; tSettings%dGibbsRetention = 1D-5
        case(9)
            cName = 'Gibbs retention upper'; cChangedValue = '1E-3'; tSettings%dGibbsRetention = 1D-3
        case(10)
            cName = 'Progress allowance lower'; cChangedValue = '1.02'; tSettings%dProgressAllowance = 1.02D0
        case(11)
            cName = 'Progress allowance upper'; cChangedValue = '1.10'; tSettings%dProgressAllowance = 1.10D0
        case(12)
            cName = 'Update ratio lower'; cChangedValue = '1.10'; tSettings%dUpdateRatio = 1.10D0
        case(13)
            cName = 'Update ratio upper'; cChangedValue = '1.50'; tSettings%dUpdateRatio = 1.50D0
        case(14)
            cName = 'Direction cosine conservative'; cChangedValue = '0.95'; tSettings%dDirectionCosine = 0.95D0
        case(15)
            cName = 'Direction cosine permissive'; cChangedValue = '0.85'; tSettings%dDirectionCosine = 0.85D0
        case(16)
            cName = 'Relative direction lower'; cChangedValue = '0.35'; tSettings%dRelativeDirection = 0.35D0
        case(17)
            cName = 'Relative direction upper'; cChangedValue = '0.75'; tSettings%dRelativeDirection = 0.75D0
        case(18)
            cName = 'Earlier activation policy'; cChangedValue = '0.10 / 3 / 1E-5'
            tSettings%dLocalThreshold = 0.10D0
            tSettings%iSettledPeriod = 3
            tSettings%dGibbsActivation = 1D-5
        case(19)
            cName = 'Later activation policy'; cChangedValue = '0.025 / 8 / 1E-7'
            tSettings%dLocalThreshold = 0.025D0
            tSettings%iSettledPeriod = 8
            tSettings%dGibbsActivation = 1D-7
        end select

    end subroutine ConfigureSensitivityCase


    subroutine CaptureTrustSettings(tSettings)

        type(MQMQATrustSettings), intent(out) :: tSettings

        tSettings%dLocalThreshold = dMQMQATrustLocalNormThreshold
        tSettings%dGibbsActivation = dMQMQATrustGibbsActivationTolerance
        tSettings%dGibbsRetention = dMQMQATrustGibbsRetentionTolerance
        tSettings%dProgressAllowance = dMQMQATrustProgressAllowance
        tSettings%dUpdateRatio = dMQMQATrustUpdateRatioCap
        tSettings%dDirectionCosine = dMQMQATrustDirectionCosineMin
        tSettings%dRelativeDirection = dMQMQATrustDirectionDifferenceCap
        tSettings%dEmergencyCap = dMQMQATrustEmergencyRatioCap
        tSettings%dResolvedNormFloor = dMQMQATrustResolvedNormFloor
        tSettings%iSettledPeriod = iMQMQATrustSettledAssemblagePeriod

    end subroutine CaptureTrustSettings


    subroutine ApplyTrustSettings(tSettings)

        type(MQMQATrustSettings), intent(in) :: tSettings

        dMQMQATrustLocalNormThreshold = tSettings%dLocalThreshold
        dMQMQATrustGibbsActivationTolerance = tSettings%dGibbsActivation
        dMQMQATrustGibbsRetentionTolerance = tSettings%dGibbsRetention
        dMQMQATrustProgressAllowance = tSettings%dProgressAllowance
        dMQMQATrustUpdateRatioCap = tSettings%dUpdateRatio
        dMQMQATrustDirectionCosineMin = tSettings%dDirectionCosine
        dMQMQATrustDirectionDifferenceCap = tSettings%dRelativeDirection
        dMQMQATrustEmergencyRatioCap = tSettings%dEmergencyCap
        dMQMQATrustResolvedNormFloor = tSettings%dResolvedNormFloor
        iMQMQATrustSettledAssemblagePeriod = tSettings%iSettledPeriod

    end subroutine ApplyTrustSettings


    logical function TrustSettingsMatch(tExpected)

        type(MQMQATrustSettings), intent(in) :: tExpected

        TrustSettingsMatch = &
            (dMQMQATrustLocalNormThreshold == tExpected%dLocalThreshold) .AND. &
            (dMQMQATrustGibbsActivationTolerance == tExpected%dGibbsActivation) .AND. &
            (dMQMQATrustGibbsRetentionTolerance == tExpected%dGibbsRetention) .AND. &
            (dMQMQATrustProgressAllowance == tExpected%dProgressAllowance) .AND. &
            (dMQMQATrustUpdateRatioCap == tExpected%dUpdateRatio) .AND. &
            (dMQMQATrustDirectionCosineMin == tExpected%dDirectionCosine) .AND. &
            (dMQMQATrustDirectionDifferenceCap == tExpected%dRelativeDirection) .AND. &
            (dMQMQATrustEmergencyRatioCap == tExpected%dEmergencyCap) .AND. &
            (dMQMQATrustResolvedNormFloor == tExpected%dResolvedNormFloor) .AND. &
            (iMQMQATrustSettledAssemblagePeriod == tExpected%iSettledPeriod)

    end function TrustSettingsMatch


    logical function ProductionTrustDefaultsMatch(tSettings)

        type(MQMQATrustSettings), intent(in) :: tSettings

        ProductionTrustDefaultsMatch = &
            (tSettings%dLocalThreshold == 0.05D0) .AND. &
            (tSettings%iSettledPeriod == 5) .AND. &
            (tSettings%dGibbsActivation == 1D-6) .AND. &
            (tSettings%dGibbsRetention == 1D-4) .AND. &
            (tSettings%dProgressAllowance == 1.05D0) .AND. &
            (tSettings%dUpdateRatio == 1.25D0) .AND. &
            (tSettings%dDirectionCosine == 0.90D0) .AND. &
            (tSettings%dRelativeDirection == 0.50D0) .AND. &
            (tSettings%dEmergencyCap == 1D6)

    end function ProductionTrustDefaultsMatch


    subroutine ClassifySensitivityResult(tResult,tDefault,iCase)

        type(MQMQASensitivityResult), intent(inout) :: tResult
        type(MQMQASensitivityResult), intent(in) :: tDefault
        integer, intent(in) :: iCase
        logical :: lMandatorySafety

        tResult%lMateriallyDifferent = .FALSE.
        if (iCase > 1) then
            tResult%lMateriallyDifferent = &
                (ABS(tResult%iIterations-tDefault%iIterations) >= 5) .OR. &
                (ABS(tResult%iFirstPositive-tDefault%iFirstPositive) >= 5) .OR. &
                (ABS(tResult%nFinalFull-tDefault%nFinalFull) >= 3) .OR. &
                (tResult%nReadinessResets /= tDefault%nReadinessResets)
        end if

        lMandatorySafety = tResult%lFinite .AND. tResult%lAlphaBounds .AND. &
            tResult%lReducedDocumented .AND. tResult%lFallbackSafe .AND. tResult%lDefaultsRestored
        if ((.NOT. lMandatorySafety) .OR. (.NOT. tResult%lHistoricalAgreement)) then
            tResult%iClassification = MQMQA_SENSITIVITY_FAILED
            tResult%cClassification = 'FAILED'
        else if ((tResult%iInfo == 0) .AND. tResult%lConverged .AND. &
            (tResult%nFinalFull >= 3) .AND. (.NOT. tResult%lMateriallyDifferent)) then
            tResult%iClassification = MQMQA_SENSITIVITY_ROBUST
            tResult%cClassification = 'ROBUST'
        else
            tResult%iClassification = MQMQA_SENSITIVITY_SAFE
            tResult%cClassification = 'SAFE BUT SENSITIVE'
        end if

    end subroutine ClassifySensitivityResult


    subroutine PrintAdaptiveDefaultResult(tResult)

        type(MQMQASensitivityResult), intent(in) :: tResult

        write(*,'(/,A)') 'adaptive FeTiVO production-default evidence'
        write(*,'(A,L1,A,I0,A,I0)') '  converged=',tResult%lConverged,', iterations=', &
            tResult%iIterations,', eligible solves=',tResult%nEligible
        write(*,'(A,I0,A,I0)') '  first positive/final sustained alpha-one global iterations=', &
            tResult%iFirstPositive,'/',tResult%iFirstFinalFull
        write(*,'(A,3(I0,1X))') '  full/reduced/zero alpha counts=', &
            tResult%nFull,tResult%nReduced,tResult%nZero
        write(*,'(A,I0,A,I0,A,I0)') '  readiness activations/resets/final full-alpha window=', &
            tResult%nReadinessActivations,'/',tResult%nReadinessResets,'/',tResult%nFinalFull
        write(*,'(A,7(I0,1X))') '  rejects ready/correction/ratio/solve/nonfinite/update/direction=', &
            tResult%nReject
        write(*,'(A,ES12.4,A,L1)') '  maximum selected alpha=',tResult%dMaximumAlpha, &
            ', reduced candidates documented=',tResult%lReducedDocumented
        write(*,'(A,4ES12.4)') '  scaled final differences G/x/n/N=',tResult%dDifference
        write(*,'(A,A)') '  classification=',TRIM(tResult%cClassification)

    end subroutine PrintAdaptiveDefaultResult


    subroutine PrintSensitivityReport(tResult)

        type(MQMQASensitivityResult), intent(in) :: tResult(:)
        integer :: i
        logical :: lStableNeighbourhood

        write(*,'(/,A)') 'MQ-4D FeTiVO one-at-a-time sensitivity matrix'
        write(*,'(A)') 'case                           value                  class               conv info iter elig first+ first1  F/R/Z  window act/reset max-alpha documented'
        do i = 1,SIZE(tResult)
            write(*,'(I2,1X,A30,1X,A20,1X,A18,1X,L1,1X,I4,1X,I4,1X,I4,1X,I5,1X,I5,1X,3(I3,A),I3,1X,I2,A,I2,1X,F8.3,1X,L1)') &
                i,TRIM(tResult(i)%cName),TRIM(tResult(i)%cChangedValue), &
                TRIM(tResult(i)%cClassification),tResult(i)%lConverged,tResult(i)%iInfo,tResult(i)%iIterations, &
                tResult(i)%nEligible,tResult(i)%iFirstPositive,tResult(i)%iFirstFinalFull, &
                tResult(i)%nFull,'/',tResult(i)%nReduced,'/',tResult(i)%nZero,' ', &
                tResult(i)%nFinalFull,tResult(i)%nReadinessActivations,'/', &
                tResult(i)%nReadinessResets,tResult(i)%dMaximumAlpha,tResult(i)%lReducedDocumented
        end do

        write(*,'(/,A)') 'rejection counts and final-state differences'
        write(*,'(A)') 'case  ready corr ratio solve nfin update direction      dG          dx          dn          dN'
        do i = 1,SIZE(tResult)
            write(*,'(I2,2X,7(I5,1X),4ES12.3)') i,tResult(i)%nReject,tResult(i)%dDifference
        end do

        write(*,'(/,A)') 'Earlier activation versus production default'
        call PrintPolicyComparison(tResult(18),tResult(1))
        write(*,'(/,A)') 'Later activation versus production default'
        call PrintPolicyComparison(tResult(19),tResult(1))

        lStableNeighbourhood = ALL(tResult%iClassification /= MQMQA_SENSITIVITY_FAILED)
        if (lStableNeighbourhood) then
            write(*,'(/,A)') 'FeTiVO conclusion: all tested neighbours were safe and convergent; three were materially sensitive.'
        else
            write(*,'(/,A)') 'FeTiVO conclusion: at least one neighbouring setting requires investigation.'
        end if
        write(*,'(A)') 'Claim boundary: local FeTiVO sensitivity only; universal MQMQA robustness is not established.'

    end subroutine PrintSensitivityReport


    subroutine PrintPolicyComparison(tPolicy,tDefault)

        type(MQMQASensitivityResult), intent(in) :: tPolicy, tDefault

        write(*,'(A,L1)') '  curvature active earlier = ', &
            (tPolicy%iFirstPositive > 0) .AND. (tPolicy%iFirstPositive < tDefault%iFirstPositive)
        write(*,'(A,I0,A,I0)') '  iterations policy/default = ',tPolicy%iIterations,'/',tDefault%iIterations
        write(*,'(A,I0,A,I0,A,I0,A,I0)') '  reduced policy/default = ',tPolicy%nReduced,'/', &
            tDefault%nReduced,', zero policy/default = ',tPolicy%nZero,'/',tDefault%nZero
        write(*,'(A,I0,A,I0)') '  readiness resets policy/default = ', &
            tPolicy%nReadinessResets,'/',tDefault%nReadinessResets
        write(*,'(A,L1,A,I0)') '  sustained final alpha one = ',tPolicy%nFinalFull >= 3, &
            ', window = ',tPolicy%nFinalFull
        write(*,'(A,L1,4ES12.3)') '  historical agreement = ', &
            tPolicy%lHistoricalAgreement,tPolicy%dDifference

    end subroutine PrintPolicyComparison


    subroutine PrintCurrentCandidateHistory(cCaseName)

        character(len=*), intent(in) :: cCaseName
        integer :: i, iCandidate

        write(*,'(/,A,A)') 'candidate history: ',TRIM(cCaseName)
        write(*,'(A)') '  iteration  alpha       rejection-mask ready-reason  norm       norm-ratio  Gibbs-gap'
        do i = 1,MIN(nMQMQAHessianEligibleSolveCount,iterGlobalMax)
            if (dMQMQAHessianAcceptedAlphaHistory(i) < 0D0) cycle
            write(*,'(I11,2X,ES10.3,2X,I0,12X,I0,3(1X,ES11.3))') &
                iMQMQAHessianGlobalIterationHistory(i),dMQMQAHessianAcceptedAlphaHistory(i), &
                iMQMQAHessianRejectionMaskHistory(i),iMQMQAHessianReadinessReasonHistory(i), &
                dMQMQAHessianFunctionNormHistory(i),dMQMQAHessianFunctionNormRatioHistory(i), &
                dMQMQAHessianGibbsGapHistory(i)
            if ((dMQMQAHessianAcceptedAlphaHistory(i) > 0D0) .AND. &
                (dMQMQAHessianAcceptedAlphaHistory(i) < 1D0-1D-12)) then
                write(*,'(A)',ADVANCE='NO') '             larger candidates:'
                do iCandidate = 1,5
                    if (dMQMQAHessianCandidateAlphaHistory(i,iCandidate) <= &
                        dMQMQAHessianAcceptedAlphaHistory(i)+1D-14) cycle
                    write(*,'(1X,ES9.2,A,I0,A)',ADVANCE='NO') &
                        dMQMQAHessianCandidateAlphaHistory(i,iCandidate),'(', &
                        iMQMQAHessianCandidateRejectionMaskHistory(i,iCandidate),')'
                end do
                write(*,*)
            end if
        end do

    end subroutine PrintCurrentCandidateHistory


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

        call ResetMQMQAHessianControls
        call ResetMQMQAHessianAdaptiveControls
        ! The sensitivity sweep changes internal policy constants deliberately.  Restore the
        ! documented production policy even when an earlier assertion has failed.
        dMQMQATrustLocalNormThreshold = 0.05D0
        iMQMQATrustSettledAssemblagePeriod = 5
        dMQMQATrustGibbsActivationTolerance = 1D-6
        dMQMQATrustGibbsRetentionTolerance = 1D-4
        dMQMQATrustProgressAllowance = 1.05D0
        dMQMQATrustUpdateRatioCap = 1.25D0
        dMQMQATrustDirectionCosineMin = 0.90D0
        dMQMQATrustDirectionDifferenceCap = 0.50D0
        dMQMQATrustEmergencyRatioCap = 1D6
        dMQMQATrustResolvedNormFloor = 1D-6
        if (lHaveCalculation) call ResetThermoAll
        if (lSucceeded) then
            write(*,'(A)') 'TestMQMQAAdaptiveTrust: PASS'
            call EXIT(0)
        else
            write(*,'(A)') 'TestMQMQAAdaptiveTrust: FAIL <---'
            call EXIT(1)
        end if

    end subroutine FinishTest

end program TestMQMQAAdaptiveTrust
