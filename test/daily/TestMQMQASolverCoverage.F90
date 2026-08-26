!-------------------------------------------------------------------------------------------------------------
!> \file    TestMQMQASolverCoverage.F90
!> \brief   Public-database solver coverage for the MQ-4E-A checkpoint.
!>
!> \details This test compares the historical GEM solve, fixed full MQMQA curvature, and adaptive MQMQA
!!          curvature for one public plain-SUBG state and a twelve-state public SUBQ trajectory.  It records
!!          whether each corrected run actually applied curvature, reached only a documented boundary fallback,
!!          or reached only a documented adaptive-trust fallback.  A fallback-only result is never reported as
!!          curvature evidence.  Pass --report to print the complete coverage matrix.
!-------------------------------------------------------------------------------------------------------------
program TestMQMQASolverCoverage

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleParseCS
    USE ModuleGEMSolver
    USE ModuleGEMNewtonDiagnosticCapture, ONLY: ResetGEMNewtonDiagnosticCapture, &
        lMQMQADiagnosticPhasePathStudy, &
        nMQMQADiagnosticPhasePathCandidates, nMQMQADiagnosticLeadingIdentityChanges, &
        nMQMQADiagnosticEligibilityCrossings, nMQMQADiagnosticOrderingReversals, &
        nMQMQADiagnosticRemovalCrossings, &
        dMQMQADiagnosticMaxScaledForceShift, dMQMQADiagnosticMaxActiveAmountDisplacement

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

    integer, parameter :: nCases = 13, nCorrectedModes = 2
    integer, parameter :: MODE_FIXED = 1, MODE_ADAPTIVE = 2
    integer, parameter :: COVERAGE_FULL_CURVATURE = 1
    integer, parameter :: COVERAGE_REDUCED_CURVATURE = 2
    integer, parameter :: COVERAGE_BOUNDARY_FALLBACK = 3
    integer, parameter :: COVERAGE_TRUST_FALLBACK = 4
    integer, parameter :: COVERAGE_REFERENCE_DIFFERENCE = 5
    integer, parameter :: COVERAGE_FAILURE = 6
    real(8), parameter :: GIBBS_AGREEMENT_TOLERANCE = 1D-8
    real(8), parameter :: STATE_AGREEMENT_TOLERANCE = 1D-6

    type :: ReferenceState
        integer :: iInfo = -1
        integer :: iIterations = 0
        logical :: lConverged = .FALSE.
        real(8) :: dGibbs = 0D0
        real(8), allocatable :: dFraction(:), dMoles(:), dPhase(:)
    end type ReferenceState

    type :: CoverageResult
        character(len=20) :: cCase = ''
        character(len=10) :: cModel = ''
        character(len=12) :: cMode = ''
        character(len=28) :: cClassification = ''
        integer :: iClassification = COVERAGE_FAILURE
        integer :: iInfo = -1
        integer :: iIterations = 0
        integer :: iHistoricalIterations = 0
        integer :: nEligible = 0
        integer :: nApplied = 0
        integer :: nAccepted = 0
        integer :: nFull = 0
        integer :: nReduced = 0
        integer :: nZero = 0
        integer :: nBoundary = 0
        integer :: nAggregateFailure = 0
        integer :: nApplicationFailure = 0
        integer :: nSolveFallback = 0
        integer :: nNonfiniteFallback = 0
        integer :: nFinalFull = 0
        integer :: nReadinessActivations = 0
        integer :: nReadinessResets = 0
        integer :: nReject(7) = 0
        integer :: nPath(5) = 0
        real(8) :: dMinimumBoundaryFraction = 1D0
        real(8) :: dIterationRatio = HUGE(1D0)
        real(8) :: dMaximumAlpha = 0D0
        real(8) :: dMaximumRatioA = 0D0
        real(8) :: dMaximumRatioB = 0D0
        real(8) :: dDifference(4) = HUGE(1D0)
        real(8) :: dPathMaximum(2) = 0D0
        logical :: lConverged = .FALSE.
        logical :: lFinite = .FALSE.
        logical :: lHistoricalAgreement = .FALSE.
        logical :: lReducedDocumented = .FALSE.
    end type CoverageResult

    type(CoverageResult) :: tResult(nCases,nCorrectedModes)
    type(ReferenceState) :: tReference
    integer :: iArgument, iCase, iMode
    logical :: lHaveCalculation, lPass, lReport
    character(len=32) :: cArgument

    lHaveCalculation = .FALSE.
    lPass = .TRUE.
    lReport = .FALSE.
    do iArgument = 1,COMMAND_ARGUMENT_COUNT()
        call GET_COMMAND_ARGUMENT(iArgument,cArgument)
        lReport = lReport .OR. (TRIM(cArgument) == '--report')
    end do

    do iCase = 1,nCases
        call RunHistoricalReference(iCase,tReference,lPass)
        do iMode = 1,nCorrectedModes
            call RunCorrectedCase(iCase,iMode,tReference,tResult(iCase,iMode))
            lPass = lPass .AND. (tResult(iCase,iMode)%iClassification /= COVERAGE_FAILURE)
            if (iMode == MODE_ADAPTIVE) lPass = lPass .AND. tResult(iCase,iMode)%lHistoricalAgreement
        end do
        call DestroyReference(tReference)
    end do

    ! The portable coverage set must include actual full-curvature solver evidence for both supported model
    ! routes.  Other individual states may legitimately be boundary- or trust-fallback evidence.
    lPass = lPass .AND. ANY(tResult(:,MODE_FIXED)%nAccepted > 0 .AND. &
        tResult(:,MODE_FIXED)%cModel == 'SUBG')
    lPass = lPass .AND. ANY(tResult(:,MODE_FIXED)%nAccepted > 0 .AND. &
        tResult(:,MODE_FIXED)%cModel == 'SUBQ')
    lPass = lPass .AND. ANY(tResult(:,MODE_ADAPTIVE)%iClassification == COVERAGE_FULL_CURVATURE .AND. &
        tResult(:,MODE_ADAPTIVE)%cModel == 'SUBG')
    lPass = lPass .AND. ANY(tResult(:,MODE_ADAPTIVE)%iClassification == COVERAGE_FULL_CURVATURE .AND. &
        tResult(:,MODE_ADAPTIVE)%cModel == 'SUBQ')

    if (lReport) call PrintCoverageReport(tResult)
    call FinishTest(lPass)

contains

    subroutine RunHistoricalReference(iCase,tState,lAllPass)

        integer, intent(in) :: iCase
        type(ReferenceState), intent(inout) :: tState
        logical, intent(inout) :: lAllPass

        call ResetMQMQAHessianControls
        call ResetMQMQAHessianAdaptiveControls
        call PrepareCase(iCase)
        if (INFOThermo == 0) call Thermochimica

        tState%iInfo = INFOThermo
        tState%iIterations = iterGlobal
        tState%lConverged = lConverged
        lAllPass = lAllPass .AND. (INFOThermo == 0) .AND. lConverged
        if ((INFOThermo /= 0) .OR. (.NOT. lConverged)) return

        tState%dGibbs = dGibbsEnergySys
        allocate(tState%dFraction(SIZE(dMolFraction)),tState%dMoles(SIZE(dMolesSpecies)), &
            tState%dPhase(SIZE(dMolesPhase)))
        tState%dFraction = dMolFraction
        tState%dMoles = dMolesSpecies
        tState%dPhase = dMolesPhase

        ! Default-off coverage is exact only when no MQMQA diagnostic path was entered.
        lAllPass = lAllPass .AND. (.NOT. lUseMQMQAExactHessian) .AND. &
            (nMQMQAHessianApplyCount == 0) .AND. (nMQMQAHessianEligibleSolveCount == 0) .AND. &
            (nMQMQAHessianInteriorFallbackCount == 0)

    end subroutine RunHistoricalReference


    subroutine RunCorrectedCase(iCase,iMode,tState,tCoverage)

        integer, intent(in) :: iCase, iMode
        type(ReferenceState), intent(in) :: tState
        type(CoverageResult), intent(out) :: tCoverage
        integer :: i, iCandidate, iInfo
        logical :: lDimensionsMatch

        tCoverage = CoverageResult()
        call ResetGEMNewtonDiagnosticCapture
        call DescribeCase(iCase,tCoverage%cCase,tCoverage%cModel)
        if (iMode == MODE_FIXED) then
            tCoverage%cMode = 'fixed-one'
            call ResetMQMQAHessianAdaptiveControls
            call SetMQMQAHessianControls(.TRUE.,1D0,iInfo)
        else
            tCoverage%cMode = 'adaptive'
            call ResetMQMQAHessianControls
            call SetMQMQAHessianAdaptiveControls(.TRUE.,1D0,iInfo)
            lMQMQADiagnosticPhasePathStudy = .TRUE.
        end if

        if (iInfo == 0) then
            call PrepareCase(iCase)
            if (INFOThermo == 0) call Thermochimica
        else
            INFOThermo = iInfo
        end if

        tCoverage%iInfo = INFOThermo
        tCoverage%lConverged = lConverged
        tCoverage%iIterations = iterGlobal
        tCoverage%iHistoricalIterations = tState%iIterations
        if (tState%iIterations > 0) tCoverage%dIterationRatio = &
            DFLOAT(tCoverage%iIterations)/DFLOAT(tState%iIterations)
        tCoverage%nEligible = nMQMQAHessianEligibleSolveCount
        tCoverage%nApplied = nMQMQAHessianApplyCount
        tCoverage%nAccepted = nMQMQAHessianAcceptedSolveCount
        tCoverage%nFull = nMQMQAHessianFullAlphaCount
        tCoverage%nReduced = nMQMQAHessianReducedAlphaCount
        tCoverage%nZero = nMQMQAHessianZeroAlphaCount
        tCoverage%nBoundary = nMQMQAHessianInteriorFallbackCount
        tCoverage%nAggregateFailure = nMQMQAHessianAggregateFailureCount
        tCoverage%nApplicationFailure = nMQMQAHessianApplicationFailureCount
        tCoverage%nSolveFallback = nMQMQAHessianDGESVFallbackCount
        tCoverage%nNonfiniteFallback = nMQMQAHessianNonfiniteFallbackCount
        tCoverage%nFinalFull = nMQMQAHessianFinalFullAlphaWindow
        tCoverage%nReadinessActivations = nMQMQAHessianReadinessActivationCount
        tCoverage%nReadinessResets = nMQMQAHessianReadinessResetCount
        tCoverage%nReject = [nMQMQAHessianRejectNonlinear,nMQMQAHessianRejectCorrection, &
            nMQMQAHessianRejectRatio,nMQMQAHessianRejectDGESV,nMQMQAHessianRejectNonfinite, &
            nMQMQAHessianRejectUpdate,nMQMQAHessianRejectDirection]
        tCoverage%nPath = [nMQMQADiagnosticPhasePathCandidates,nMQMQADiagnosticLeadingIdentityChanges, &
            nMQMQADiagnosticEligibilityCrossings,nMQMQADiagnosticOrderingReversals, &
            nMQMQADiagnosticRemovalCrossings]
        tCoverage%dPathMaximum = [dMQMQADiagnosticMaxScaledForceShift, &
            dMQMQADiagnosticMaxActiveAmountDisplacement]
        tCoverage%dMinimumBoundaryFraction = dMQMQAHessianMinimumRejectedFraction
        tCoverage%dMaximumAlpha = dMQMQAHessianMaxSelectedAlpha
        tCoverage%dMaximumRatioA = dMQMQAHessianMaxRatioA
        tCoverage%dMaximumRatioB = dMQMQAHessianMaxRatioB
        tCoverage%lReducedDocumented = .TRUE.
        if (iMode == MODE_ADAPTIVE) then
            do i = 1,MIN(nMQMQAHessianEligibleSolveCount,iterGlobalMax)
                if ((dMQMQAHessianAcceptedAlphaHistory(i) <= 0D0) .OR. &
                    (dMQMQAHessianAcceptedAlphaHistory(i) >= 1D0-1D-12)) cycle
                do iCandidate = 1,5
                    if (dMQMQAHessianCandidateAlphaHistory(i,iCandidate) <= &
                        dMQMQAHessianAcceptedAlphaHistory(i)+1D-14) cycle
                    tCoverage%lReducedDocumented = tCoverage%lReducedDocumented .AND. &
                        (iMQMQAHessianCandidateRejectionMaskHistory(i,iCandidate) /= 0)
                end do
            end do
        end if

        lDimensionsMatch = ALLOCATED(tState%dFraction) .AND. ALLOCATED(tState%dMoles) .AND. &
            ALLOCATED(tState%dPhase) .AND. ALLOCATED(dMolFraction) .AND. ALLOCATED(dMolesSpecies) .AND. &
            ALLOCATED(dMolesPhase)
        if (lDimensionsMatch) lDimensionsMatch = (SIZE(dMolFraction) == SIZE(tState%dFraction)) .AND. &
            (SIZE(dMolesSpecies) == SIZE(tState%dMoles)) .AND. (SIZE(dMolesPhase) == SIZE(tState%dPhase))
        tCoverage%lFinite = lDimensionsMatch .AND. IEEE_IS_FINITE(dGibbsEnergySys)
        if (lDimensionsMatch) tCoverage%lFinite = tCoverage%lFinite .AND. &
            ALL(IEEE_IS_FINITE(dMolFraction)) .AND. ALL(IEEE_IS_FINITE(dMolesSpecies)) .AND. &
            ALL(IEEE_IS_FINITE(dMolesPhase))
        if (tCoverage%lFinite) then
            tCoverage%dDifference(1) = ABS(dGibbsEnergySys-tState%dGibbs)/DMAX1(1D0,ABS(tState%dGibbs))
            tCoverage%dDifference(2) = MAXVAL(ABS(dMolFraction-tState%dFraction))/ &
                DMAX1(1D0,MAXVAL(ABS(tState%dFraction)))
            tCoverage%dDifference(3) = MAXVAL(ABS(dMolesSpecies-tState%dMoles))/ &
                DMAX1(1D0,MAXVAL(ABS(tState%dMoles)))
            tCoverage%dDifference(4) = MAXVAL(ABS(dMolesPhase-tState%dPhase))/ &
                DMAX1(1D0,MAXVAL(ABS(tState%dPhase)))
        end if
        tCoverage%lHistoricalAgreement = (tCoverage%dDifference(1) <= GIBBS_AGREEMENT_TOLERANCE) .AND. &
            ALL(tCoverage%dDifference(2:4) <= STATE_AGREEMENT_TOLERANCE)

        call ClassifyCoverage(iMode,tCoverage)
        call ResetMQMQAHessianControls
        call ResetMQMQAHessianAdaptiveControls

    end subroutine RunCorrectedCase


    subroutine ClassifyCoverage(iMode,tCoverage)

        integer, intent(in) :: iMode
        type(CoverageResult), intent(inout) :: tCoverage
        integer :: nHardFailures, nPositive
        logical :: lSafeResult

        nHardFailures = tCoverage%nAggregateFailure+tCoverage%nApplicationFailure+ &
            tCoverage%nSolveFallback+tCoverage%nNonfiniteFallback
        if (iMode == MODE_FIXED) then
            nPositive = tCoverage%nAccepted
        else
            nPositive = tCoverage%nFull+tCoverage%nReduced
        end if
        lSafeResult = (tCoverage%iInfo == 0) .AND. tCoverage%lConverged .AND. tCoverage%lFinite .AND. &
            tCoverage%lReducedDocumented .AND. (nHardFailures == 0)

        if (.NOT. lSafeResult) then
            tCoverage%iClassification = COVERAGE_FAILURE
            tCoverage%cClassification = 'FAILURE'
        else if (.NOT. tCoverage%lHistoricalAgreement) then
            tCoverage%iClassification = COVERAGE_REFERENCE_DIFFERENCE
            tCoverage%cClassification = 'REFERENCE_STATE_DIFFERENCE'
        else if ((iMode == MODE_FIXED) .AND. (nPositive > 0)) then
            tCoverage%iClassification = COVERAGE_FULL_CURVATURE
            tCoverage%cClassification = 'FULL_CURVATURE_EVIDENCE'
        else if ((iMode == MODE_ADAPTIVE) .AND. (tCoverage%nFull > 0)) then
            tCoverage%iClassification = COVERAGE_FULL_CURVATURE
            tCoverage%cClassification = 'FULL_CURVATURE_EVIDENCE'
        else if ((iMode == MODE_ADAPTIVE) .AND. (tCoverage%nReduced > 0)) then
            tCoverage%iClassification = COVERAGE_REDUCED_CURVATURE
            tCoverage%cClassification = 'REDUCED_CURVATURE_EVIDENCE'
        else if (tCoverage%nBoundary > 0) then
            tCoverage%iClassification = COVERAGE_BOUNDARY_FALLBACK
            tCoverage%cClassification = 'SAFE_BOUNDARY_FALLBACK'
        else if ((iMode == MODE_ADAPTIVE) .AND. (tCoverage%nEligible > 0) .AND. (tCoverage%nZero > 0)) then
            tCoverage%iClassification = COVERAGE_TRUST_FALLBACK
            tCoverage%cClassification = 'SAFE_TRUST_FALLBACK'
        else
            tCoverage%iClassification = COVERAGE_FAILURE
            tCoverage%cClassification = 'FAILURE'
        end if

    end subroutine ClassifyCoverage


    subroutine PrepareCase(iCase)

        integer, intent(in) :: iCase
        real(8) :: dPath

        call PrepareFreshCalculation
        dPressure = 1D0
        dElementMass = 0D0
        if (iCase == 1) then
            cThermoFileName = DATA_DIRECTORY // 'CuFeC-Kang.dat'
            dTemperature = 1400D0
            dElementMass(6) = 1D0
            dElementMass(26) = 1D0
            dElementMass(29) = 1D0
        else
            dPath = DFLOAT(iCase-2)/11D0
            cThermoFileName = DATA_DIRECTORY // 'FeTiVO.dat'
            dTemperature = 1900D0+200D0*dPath
            dElementMass(8) = 2D0
            dElementMass(22) = 0.45D0+0.10D0*dPath
            dElementMass(23) = 0.5D0
            dElementMass(26) = 0.55D0-0.10D0*dPath
        end if
        call ParseCSDataFile(cThermoFileName)

    end subroutine PrepareCase


    subroutine DescribeCase(iCase,cCase,cModel)

        integer, intent(in) :: iCase
        character(len=*), intent(out) :: cCase, cModel

        if (iCase == 1) then
            cCase = 'CuFeC 1400 K'
            cModel = 'SUBG'
        else
            write(cCase,'(A,I2.2)') 'FeTiVO step ',iCase-1
            cModel = 'SUBQ'
        end if

    end subroutine DescribeCase


    subroutine PrepareFreshCalculation

        if (lHaveCalculation) call ResetThermoAll
        lHaveCalculation = .TRUE.
        cInputUnitTemperature = 'K'
        cInputUnitPressure = 'atm'
        cInputUnitMass = 'moles'
        dElementMass = 0D0

    end subroutine PrepareFreshCalculation


    subroutine DestroyReference(tState)

        type(ReferenceState), intent(inout) :: tState

        if (ALLOCATED(tState%dFraction)) deallocate(tState%dFraction)
        if (ALLOCATED(tState%dMoles)) deallocate(tState%dMoles)
        if (ALLOCATED(tState%dPhase)) deallocate(tState%dPhase)
        tState = ReferenceState()

    end subroutine DestroyReference


    subroutine PrintCoverageReport(tCoverage)

        type(CoverageResult), intent(in) :: tCoverage(:,:)
        integer :: i, j

        write(*,'(/,A)') 'MQ-4E-A public-database solver coverage'
        write(*,'(A)') 'case                 model mode         class                       iter elig apply accept  F/R/Z  boundary final1'
        do i = 1,SIZE(tCoverage,1)
            do j = 1,SIZE(tCoverage,2)
                write(*,'(A20,1X,A5,1X,A12,1X,A27,1X,I4,1X,I4,1X,I5,1X,I6,1X,3(I3,A),I5,1X,I5)') &
                    TRIM(tCoverage(i,j)%cCase),TRIM(tCoverage(i,j)%cModel),TRIM(tCoverage(i,j)%cMode), &
                    TRIM(tCoverage(i,j)%cClassification),tCoverage(i,j)%iIterations, &
                    tCoverage(i,j)%nEligible,tCoverage(i,j)%nApplied,tCoverage(i,j)%nAccepted, &
                    tCoverage(i,j)%nFull,'/',tCoverage(i,j)%nReduced,'/',tCoverage(i,j)%nZero,' ', &
                    tCoverage(i,j)%nBoundary,tCoverage(i,j)%nFinalFull
            end do
        end do
        write(*,'(/,A)') 'scaled final differences from each historical reference'
        write(*,'(A)') 'case                 mode               dG          dx          dn          dN      min boundary x'
        do i = 1,SIZE(tCoverage,1)
            do j = 1,SIZE(tCoverage,2)
                write(*,'(A20,1X,A12,5ES12.3)') TRIM(tCoverage(i,j)%cCase),TRIM(tCoverage(i,j)%cMode), &
                    tCoverage(i,j)%dDifference,tCoverage(i,j)%dMinimumBoundaryFraction
            end do
        end do
        write(*,'(/,A)') 'iteration cost relative to the historical solve'
        write(*,'(A)') 'case                 historical fixed ratio adaptive ratio'
        do i = 1,SIZE(tCoverage,1)
            write(*,'(A20,1X,I10,2F12.3)') TRIM(tCoverage(i,MODE_FIXED)%cCase), &
                tCoverage(i,MODE_FIXED)%iHistoricalIterations,tCoverage(i,MODE_FIXED)%dIterationRatio, &
                tCoverage(i,MODE_ADAPTIVE)%dIterationRatio
        end do
        write(*,'(/,A)') 'adaptive trust and correction diagnostics'
        write(*,'(A)') 'case                 act/reset max-alpha  max-rhoA    max-rhoB    documented reject ready/corr/ratio/solve/nfin/update/dir'
        do i = 1,SIZE(tCoverage,1)
            j = MODE_ADAPTIVE
            write(*,'(A20,1X,I3,A,I3,3ES12.3,1X,L1,1X,7(I4,1X))') TRIM(tCoverage(i,j)%cCase), &
                tCoverage(i,j)%nReadinessActivations,'/',tCoverage(i,j)%nReadinessResets, &
                tCoverage(i,j)%dMaximumAlpha,tCoverage(i,j)%dMaximumRatioA,tCoverage(i,j)%dMaximumRatioB, &
                tCoverage(i,j)%lReducedDocumented,tCoverage(i,j)%nReject
        end do
        write(*,'(/,A)') 'diagnostic same-state phase-path metrics (adaptive candidates only)'
        write(*,'(A)') 'case                 candidates identity eligibility ordering removal  max-force-shift max-amount-shift'
        do i = 1,SIZE(tCoverage,1)
            j = MODE_ADAPTIVE
            write(*,'(A20,1X,5I9,2ES17.5)') TRIM(tCoverage(i,j)%cCase), &
                tCoverage(i,j)%nPath,tCoverage(i,j)%dPathMaximum
        end do
        write(*,'(/,A)') 'Classification boundary: FULL_CURVATURE_EVIDENCE requires at least one full-alpha solve.'
        write(*,'(A)') 'REDUCED_CURVATURE_EVIDENCE means only damped positive corrections were used.'
        write(*,'(A)') 'SAFE_* rows are transactional fallback evidence, not curvature-effect evidence.'
        write(*,'(A,ES9.2,A,ES9.2)') 'Historical agreement tolerances: Gibbs=',GIBBS_AGREEMENT_TOLERANCE, &
            ', state arrays=',STATE_AGREEMENT_TOLERANCE

    end subroutine PrintCoverageReport


    subroutine FinishTest(lSucceeded)

        logical, intent(in) :: lSucceeded

        call ResetMQMQAHessianControls
        call ResetMQMQAHessianAdaptiveControls
        if (lHaveCalculation) call ResetThermoAll
        if (lSucceeded) then
            write(*,'(A)') 'TestMQMQASolverCoverage: PASS'
            call EXIT(0)
        else
            write(*,'(A)') 'TestMQMQASolverCoverage: FAIL <---'
            call EXIT(1)
        end if

    end subroutine FinishTest

end program TestMQMQASolverCoverage
