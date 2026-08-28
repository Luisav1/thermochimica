!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleGEMNewtonDiagnosticCapture.f90
!> \brief   Opt-in snapshot of the unreduced GEMNewton linear system for verification tests.
!>
!> \details MQMQA response mapping must be compared with the matrix and right-hand side that GEMNewton actually
!!          assembles, not only with a second copy of the source formulas in a test. This module provides a
!!          deliberately small diagnostic boundary: when capture is requested, GEMNewton copies its completed
!!          baseline system immediately before any experimental curvature correction or linear solve.
!!
!!          Capture is disabled by default, does not print, and never changes the supplied matrix or vector.
!-------------------------------------------------------------------------------------------------------------
module ModuleGEMNewtonDiagnosticCapture

    implicit none
    private

    logical, public :: lCaptureGEMNewtonSystem = .FALSE.
    logical, public :: lGEMNewtonSystemCaptured = .FALSE.
    logical, public :: lCaptureGEMNewtonCorrectedSystem = .FALSE.
    logical, public :: lGEMNewtonCorrectedSystemCaptured = .FALSE.
    ! When enabled, retain the baseline and corrected matrices from the first
    ! GEM call that actually reaches an MQMQA corrected trial.  This keeps the
    ! two snapshots paired at one identical pre-step state.
    logical, public :: lCaptureFirstGEMNewtonCorrectionPair = .FALSE.
    ! Test-only same-matrix comparison of the production LU solution with an
    ! SVD rank-revealing minimum-norm solution.  Results never influence the
    ! live Newton update or phase-assemblage decision.
    logical, public :: lMQMQADiagnosticMinimumNormStudy = .FALSE.
    logical, public :: lMQMQADiagnosticMinimumNormCaptured = .FALSE.
    integer, public :: iMQMQADiagnosticMinimumNormRank(2) = 0
    integer, public :: iMQMQADiagnosticMinimumNormInfo(2,2) = 0
    integer, public :: iMQMQADiagnosticMinimumNormDecision(2,2,2) = 0
    ! Solver residuals are reported both relative to the right-hand side and
    ! as a normwise backward error.  Rows are LU/minimum-norm; columns are
    ! baseline/corrected.
    real(8), public :: dMQMQADiagnosticMinimumNormRHSResidual(2,2) = 0D0
    real(8), public :: dMQMQADiagnosticMinimumNormResidual(2,2) = 0D0
    real(8), public :: dMQMQADiagnosticMinimumNormSolutionNorm(2,2) = 0D0
    real(8), public :: dMQMQADiagnosticMinimumNormRelativeDifference(2) = 0D0
    real(8), public :: dMQMQADiagnosticMinimumNormNullResidual(2) = 0D0
    real(8), public :: dMQMQADiagnosticMinimumNormGammaDifference(2) = 0D0
    real(8), public :: dMQMQADiagnosticMinimumNormDecisionForce(3,2,2) = 0D0
    real(8), public :: dMQMQADiagnosticMinimumNormCrossSystemDifference(2) = 0D0
    real(8), public :: dMQMQADiagnosticMinimumNormCrossGammaDifference(2) = 0D0
    ! Rows are largest singular value, smallest retained value, largest
    ! discarded value, and the numerical-rank tolerance; columns are the
    ! baseline and corrected systems.
    real(8), public :: dMQMQADiagnosticMinimumNormSingularSummary(4,2) = 0D0
    ! Per discarded singular mode, rows are: singular value, relative A*v
    ! residual, relative left-hand-side compatibility projection, required
    ! singular coefficient, and the element/solution/pure variable-block
    ! norms, followed by the largest leading phase-force change under a
    ! unit maximum element-potential perturbation.  The final index selects
    ! the paired baseline/corrected system.
    real(8), allocatable, public :: dMQMQADiagnosticNullModeSummary(:,:,:)
    integer, allocatable, public :: iMQMQADiagnosticNullModeDecisionChanges(:,:)
    ! The left singular vectors identify incompatible equation combinations.
    ! Rows 1:3 are their element/solution/pure equation-block norms; rows 4:6
    ! are the signed contributions of those blocks to u^T B, normalized by
    ! ||B||; row 7 is the signed total projection.  The stored vectors make
    ! the small diagnostic systems auditable against named GEM equations.
    real(8), allocatable, public :: dMQMQADiagnosticLeftNullSummary(:,:,:)
    real(8), allocatable, public :: dMQMQADiagnosticLeftNullVector(:,:,:)
    real(8), allocatable, public :: dMQMQADiagnosticRightNullVector(:,:,:)
    ! Row metadata are [equation class, production index], where classes are
    ! 1=element balance, 2=solution phase, and 3=pure phase.
    integer, allocatable, public :: iMQMQADiagnosticEquationIdentity(:,:)
    ! Baseline-to-corrected comparison: maximum sine of the principal angles
    ! for the left and right numerical null spaces, the direct corrected RHS
    ! forcing projected onto the baseline left-null space, and the corrected
    ! RHS projected onto the baseline left-null space.  The latter two values
    ! are normalized by the baseline RHS norm.
    real(8), public :: dMQMQADiagnosticNullSpaceComparison(4) = 0D0
    ! Basis-independent active-phase structure.  The rank belongs to the
    ! element-by-active-phase stoichiometry block.  Summary rows are largest,
    ! smallest retained, and largest discarded singular values, followed by
    ! the normalized element-inventory and phase-energy closure residuals, the
    ! unnormalized phase-energy closure norm, and the number of dependent
    ! active-phase combinations.
    integer, public :: iMQMQADiagnosticPhaseStructureRank(2) = 0
    ! Per paired baseline/corrected system: number of elements, charged
    ! constraints, active solution phases, active pure phases, total active
    ! phases, nominal phase-rule limit, global iteration, and whether the
    ! count-based CorrectPhaseRule condition would fire.  This distinguishes
    ! the nominal component count from the measured active-stoichiometry rank.
    integer, public :: iMQMQADiagnosticPhaseRuleSummary(8,2) = 0
    real(8), public :: dMQMQADiagnosticPhaseStructureSummary(7,2) = 0D0
    real(8), allocatable, public :: dMQMQADiagnosticPhaseDependencyVector(:,:,:)
    real(8), allocatable, public :: dMQMQADiagnosticPhaseDependencyProjection(:,:)
    ! Rows are the raw stoichiometry-column norm and phase-energy RHS value
    ! for each active solution/pure phase equation.
    real(8), allocatable, public :: dMQMQADiagnosticPhaseEquationData(:,:,:)
    integer, parameter, public :: nMQMQADiagnosticMaxPhaseChangeChecks = 32
    integer, parameter, public :: nMQMQADiagnosticMaxActivePhases = 16
    integer, public :: nMQMQADiagnosticPhaseChangeChecks = 0
    ! Rows are global iteration, elements, charged constraints, active
    ! solution phases, active pure phases, GEMNewton INFO, pass/fail, and
    ! numerical rank of the active element-by-phase stoichiometry block.
    integer, public :: iMQMQADiagnosticPhaseChangeCheck(8,nMQMQADiagnosticMaxPhaseChangeChecks) = 0
    integer, public :: iMQMQADiagnosticPhaseChangeAssemblage( &
        nMQMQADiagnosticMaxActivePhases,nMQMQADiagnosticMaxPhaseChangeChecks) = 0
    ! Rows are maximum absolute update, active acceptance threshold, minimum
    ! active phase amount, and the ordinary amount-removal tolerance.
    real(8), public :: dMQMQADiagnosticPhaseChangeCheck(4,nMQMQADiagnosticMaxPhaseChangeChecks) = 0D0
    ! Candidate 1 is a physically grouped gauge projection of the production
    ! LU result; candidate 2 is an algebraically equivalent row/column-
    ! equilibrated exact solve.  Neither candidate is returned to GEMNewton.
    integer, public :: iMQMQADiagnosticConstrainedInfo(2,2) = 0
    integer, public :: iMQMQADiagnosticConstrainedDecision(2,2,2) = 0
    real(8), public :: dMQMQADiagnosticConstrainedRHSResidual(2,2) = 0D0
    real(8), public :: dMQMQADiagnosticConstrainedResidual(2,2) = 0D0
    real(8), public :: dMQMQADiagnosticConstrainedSolutionNorm(2,2) = 0D0
    real(8), public :: dMQMQADiagnosticConstrainedGroupNorm(3,2,2) = 0D0
    real(8), public :: dMQMQADiagnosticConstrainedDecisionForce(3,2,2) = 0D0
    ! Test-only causal probe: enter the adaptive MQMQA path but always return
    ! its independently solved historical system without testing candidates.
    logical, public :: lMQMQADiagnosticForceZeroAlpha = .FALSE.
    ! Test-only causal probe: evaluate adaptive candidates normally, but return
    ! the untouched baseline even when a positive candidate would be accepted.
    logical, public :: lMQMQADiagnosticRejectAcceptedCandidate = .FALSE.
    ! Opt-in, test-only study of whether a candidate Newton solve immediately
    ! crosses a production phase-addition or phase-removal boundary.
    logical, public :: lMQMQADiagnosticPhasePathStudy = .FALSE.
    integer, public :: nMQMQADiagnosticPhasePathCandidates = 0
    integer, public :: nMQMQADiagnosticLeadingIdentityChanges = 0
    integer, public :: nMQMQADiagnosticEligibilityCrossings = 0
    integer, public :: nMQMQADiagnosticOrderingReversals = 0
    integer, public :: nMQMQADiagnosticRemovalCrossings = 0
    real(8), public :: dMQMQADiagnosticMaxScaledForceShift = 0D0
    real(8), public :: dMQMQADiagnosticMaxActiveAmountDisplacement = 0D0
    integer, public :: nCapturedGEMNewtonVariables = 0
    integer, parameter, public :: nMQMQADiagnosticMaxAttempts = 4
    integer, parameter, public :: nMQMQADiagnosticMaxRemovalEvents = 8
    integer, parameter, public :: nMQMQADiagnosticMaxCandidateSpecies = 64
    logical, public :: lCaptureMQMQATrajectory = .FALSE.
    integer, public :: iMQMQADiagnosticCurrentAttempt = 0
    integer, public :: iMQMQADiagnosticSkipAttempt = -1
    integer, public :: iMQMQADiagnosticSkipIteration = -1
    integer, public :: nCapturedMQMQAAttempts = 0
    integer, public :: iCapturedMQMQAPoints(nMQMQADiagnosticMaxAttempts) = 0
    integer, allocatable, public :: iCapturedMQMQAAssemblage(:,:,:)
    integer, allocatable, public :: iCapturedMQMQAAcceptedCount(:,:)
    integer, allocatable, public :: iCapturedMQMQAEligibleCount(:,:)
    integer, allocatable, public :: iCapturedMQMQABoundaryCount(:,:)
    integer, allocatable, public :: iCapturedMQMQAIterLast(:,:)
    ! Exact phase-addition ranking used by CheckPhaseAssemblage.  Integer
    ! fields are [decision present, pure species index, solution phase index];
    ! real fields are [pure force, solution force, pure-minus-solution margin].
    integer, allocatable, public :: iCapturedMQMQAPhaseDecision(:,:,:)
    real(8), allocatable, public :: dCapturedMQMQAPhaseDecision(:,:,:)
    ! Same-pre-step adaptive comparison.  Candidate index 1 is the untouched
    ! alpha-zero solve and index 2 is the selected positive-alpha solve.
    ! Gamma index 1 is the current state, 2 the alpha-zero target, and 3 the
    ! selected corrected target.  Decision fields are [pure,solution] indices
    ! and [pure force,solution force,pure-minus-solution] values.
    logical, allocatable, public :: lCapturedMQMQACandidateComparison(:,:)
    integer, allocatable, public :: iCapturedMQMQACandidateDecision(:,:,:,:)
    real(8), allocatable, public :: dCapturedMQMQACandidateGamma(:,:,:,:)
    real(8), allocatable, public :: dCapturedMQMQACandidateDecision(:,:,:,:)
    integer, allocatable, public :: iCapturedMQMQACandidateSolutionSpecies(:,:,:,:)
    real(8), allocatable, public :: dCapturedMQMQACandidateSolutionFraction(:,:,:,:)
    real(8), allocatable, public :: dCapturedMQMQACandidateAlpha(:,:)
    ! Removal events preserve the triggering state and result.  Integer fields
    ! are [phase id, phase class, event type, outcome]; real fields are
    ! [current amount, previous amount, change, tolerance, condition limit].
    integer, allocatable, public :: nCapturedMQMQARemovalEvents(:,:)
    integer, allocatable, public :: iCapturedMQMQARemovalEvent(:,:,:,:)
    real(8), allocatable, public :: dCapturedMQMQARemovalEvent(:,:,:,:)
    real(8), allocatable, public :: dCapturedMQMQAAlpha(:,:)
    real(8), allocatable, public :: dCapturedMQMQAFunctionNorm(:,:)
    real(8), allocatable, public :: dCapturedMQMQAMinGibbs(:,:)
    real(8), allocatable, public :: dCapturedMQMQAMinimumBoundaryFraction(:,:)
    real(8), allocatable, public :: dCapturedMQMQAPhaseMoles(:,:,:)
    real(8), allocatable, public :: dCapturedMQMQAElementPotential(:,:,:)
    real(8), allocatable, public :: dCapturedGEMNewtonA(:,:), dCapturedGEMNewtonB(:)
    real(8), allocatable, public :: dCapturedGEMNewtonCorrectedA(:,:), dCapturedGEMNewtonCorrectedB(:)

    public :: CaptureGEMNewtonSystem, CaptureGEMNewtonCorrectedSystem, ResetGEMNewtonDiagnosticCapture
    public :: BeginMQMQATrajectoryAttempt, CaptureMQMQATrajectoryPoint
    public :: CaptureMQMQAPhaseDecision
    public :: CaptureMQMQACandidateComparison, CaptureMQMQAPhaseRemovalEvent
    public :: RecordMQMQAPhasePathCandidate
    public :: CaptureMQMQAPhaseChangeCheck

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Record one initialization/global candidate check without changing its outcome.
    !---------------------------------------------------------------------------------------------------------
    subroutine CaptureMQMQAPhaseChangeCheck(iIteration,nElementCount,nChargeCount,nSolutionCount, &
        nPureCount,iInfo,lPass,iActiveRank,dMaximumUpdate,dThreshold,dMinimumAmount,dRemovalTolerance, &
        iAssemblageIn)

        integer, intent(in) :: iIteration, nElementCount, nChargeCount, nSolutionCount, nPureCount, iInfo
        integer, intent(in) :: iActiveRank
        logical, intent(in) :: lPass
        real(8), intent(in) :: dMaximumUpdate, dThreshold, dMinimumAmount, dRemovalTolerance
        integer, intent(in) :: iAssemblageIn(:)
        integer :: i, iCheck, iWrite

        if (.NOT. lMQMQADiagnosticMinimumNormStudy) return
        if (nMQMQADiagnosticPhaseChangeChecks >= nMQMQADiagnosticMaxPhaseChangeChecks) return
        iCheck = nMQMQADiagnosticPhaseChangeChecks+1
        nMQMQADiagnosticPhaseChangeChecks = iCheck
        iMQMQADiagnosticPhaseChangeCheck(:,iCheck) = [iIteration,nElementCount,nChargeCount, &
            nSolutionCount,nPureCount,iInfo,MERGE(1,0,lPass),iActiveRank]
        dMQMQADiagnosticPhaseChangeCheck(:,iCheck) = &
            [dMaximumUpdate,dThreshold,dMinimumAmount,dRemovalTolerance]
        iWrite = 0
        do i = 1,nPureCount
            if (iWrite >= nMQMQADiagnosticMaxActivePhases) exit
            iWrite = iWrite+1
            iMQMQADiagnosticPhaseChangeAssemblage(iWrite,iCheck) = iAssemblageIn(i)
        end do
        do i = nElementCount,nElementCount-nSolutionCount+1,-1
            if (iWrite >= nMQMQADiagnosticMaxActivePhases) exit
            iWrite = iWrite+1
            iMQMQADiagnosticPhaseChangeAssemblage(iWrite,iCheck) = iAssemblageIn(i)
        end do

    end subroutine CaptureMQMQAPhaseChangeCheck

    !---------------------------------------------------------------------------------------------------------
    !> \brief Accumulate model-independent same-state phase-path diagnostics for one alpha candidate.
    !---------------------------------------------------------------------------------------------------------
    subroutine RecordMQMQAPhasePathCandidate(lIdentityChange,lEligibilityCrossing,lOrderingReversal, &
        lRemovalCrossing,dScaledForceShift,dAmountDisplacement)

        logical, intent(in) :: lIdentityChange, lEligibilityCrossing, lOrderingReversal
        logical, intent(in) :: lRemovalCrossing
        real(8), intent(in) :: dScaledForceShift, dAmountDisplacement

        if (.NOT. lMQMQADiagnosticPhasePathStudy) return
        nMQMQADiagnosticPhasePathCandidates = nMQMQADiagnosticPhasePathCandidates+1
        if (lIdentityChange) nMQMQADiagnosticLeadingIdentityChanges = &
            nMQMQADiagnosticLeadingIdentityChanges+1
        if (lEligibilityCrossing) nMQMQADiagnosticEligibilityCrossings = &
            nMQMQADiagnosticEligibilityCrossings+1
        if (lOrderingReversal) nMQMQADiagnosticOrderingReversals = &
            nMQMQADiagnosticOrderingReversals+1
        if (lRemovalCrossing) nMQMQADiagnosticRemovalCrossings = &
            nMQMQADiagnosticRemovalCrossings+1
        dMQMQADiagnosticMaxScaledForceShift = DMAX1(dMQMQADiagnosticMaxScaledForceShift,dScaledForceShift)
        dMQMQADiagnosticMaxActiveAmountDisplacement = &
            DMAX1(dMQMQADiagnosticMaxActiveAmountDisplacement,dAmountDisplacement)

    end subroutine RecordMQMQAPhasePathCandidate

    !---------------------------------------------------------------------------------------------------------
    !> \brief Copy one completed baseline Newton system when the opt-in request is active.
    !---------------------------------------------------------------------------------------------------------
    subroutine CaptureGEMNewtonSystem(dA,dB,nVar)

        integer, intent(in) :: nVar
        real(8), intent(in) :: dA(:,:), dB(:)

        if (.NOT. lCaptureGEMNewtonSystem) return
        if (lCaptureFirstGEMNewtonCorrectionPair .AND. lGEMNewtonCorrectedSystemCaptured) return
        if ((nVar <= 0) .OR. (SIZE(dA,1) < nVar) .OR. (SIZE(dA,2) < nVar) .OR. &
            (SIZE(dB) < nVar)) return

        if (allocated(dCapturedGEMNewtonA)) deallocate(dCapturedGEMNewtonA)
        if (allocated(dCapturedGEMNewtonB)) deallocate(dCapturedGEMNewtonB)
        allocate(dCapturedGEMNewtonA(nVar,nVar),dCapturedGEMNewtonB(nVar))
        dCapturedGEMNewtonA = dA(1:nVar,1:nVar)
        dCapturedGEMNewtonB = dB(1:nVar)
        nCapturedGEMNewtonVariables = nVar
        lGEMNewtonSystemCaptured = .TRUE.

    end subroutine CaptureGEMNewtonSystem


    !---------------------------------------------------------------------------------------------------------
    !> \brief Copy the MQMQA-corrected trial system before its destructive linear solve.
    !---------------------------------------------------------------------------------------------------------
    subroutine CaptureGEMNewtonCorrectedSystem(dA,dB,nVar)

        integer, intent(in) :: nVar
        real(8), intent(in) :: dA(:,:), dB(:)

        if (.NOT. lCaptureGEMNewtonCorrectedSystem) return
        if (lCaptureFirstGEMNewtonCorrectionPair .AND. lGEMNewtonCorrectedSystemCaptured) return
        if ((nVar <= 0) .OR. (SIZE(dA,1) < nVar) .OR. (SIZE(dA,2) < nVar) .OR. &
            (SIZE(dB) < nVar)) return

        if (allocated(dCapturedGEMNewtonCorrectedA)) deallocate(dCapturedGEMNewtonCorrectedA)
        if (allocated(dCapturedGEMNewtonCorrectedB)) deallocate(dCapturedGEMNewtonCorrectedB)
        allocate(dCapturedGEMNewtonCorrectedA(nVar,nVar),dCapturedGEMNewtonCorrectedB(nVar))
        dCapturedGEMNewtonCorrectedA = dA(1:nVar,1:nVar)
        dCapturedGEMNewtonCorrectedB = dB(1:nVar)
        lGEMNewtonCorrectedSystemCaptured = .TRUE.

    end subroutine CaptureGEMNewtonCorrectedSystem


    !---------------------------------------------------------------------------------------------------------
    !> \brief Start one independently initialized GEM attempt in the opt-in MQMQA trajectory trace.
    !---------------------------------------------------------------------------------------------------------
    subroutine BeginMQMQATrajectoryAttempt(nElements,nIterations)

        integer, intent(in) :: nElements, nIterations

        if (.NOT. lCaptureMQMQATrajectory) return
        if ((nElements <= 0) .OR. (nIterations <= 0)) return
        if (.NOT. allocated(iCapturedMQMQAAssemblage)) then
            allocate(iCapturedMQMQAAssemblage(nElements,nIterations,nMQMQADiagnosticMaxAttempts), &
                iCapturedMQMQAAcceptedCount(nIterations,nMQMQADiagnosticMaxAttempts), &
                iCapturedMQMQAEligibleCount(nIterations,nMQMQADiagnosticMaxAttempts), &
                iCapturedMQMQABoundaryCount(nIterations,nMQMQADiagnosticMaxAttempts), &
                iCapturedMQMQAIterLast(nIterations,nMQMQADiagnosticMaxAttempts), &
                iCapturedMQMQAPhaseDecision(3,nIterations,nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQAPhaseDecision(3,nIterations,nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQAPhaseMoles(nElements,nIterations,nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQAElementPotential(nElements,nIterations,nMQMQADiagnosticMaxAttempts), &
                lCapturedMQMQACandidateComparison(nIterations,nMQMQADiagnosticMaxAttempts), &
                iCapturedMQMQACandidateDecision(2,2,nIterations,nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQACandidateGamma(nElements,3,nIterations,nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQACandidateDecision(3,2,nIterations,nMQMQADiagnosticMaxAttempts), &
                iCapturedMQMQACandidateSolutionSpecies(nMQMQADiagnosticMaxCandidateSpecies,2, &
                    nIterations,nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQACandidateSolutionFraction(nMQMQADiagnosticMaxCandidateSpecies,2, &
                    nIterations,nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQACandidateAlpha(nIterations,nMQMQADiagnosticMaxAttempts), &
                nCapturedMQMQARemovalEvents(nIterations,nMQMQADiagnosticMaxAttempts), &
                iCapturedMQMQARemovalEvent(4,nMQMQADiagnosticMaxRemovalEvents,nIterations, &
                    nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQARemovalEvent(5,nMQMQADiagnosticMaxRemovalEvents,nIterations, &
                    nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQAAlpha(nIterations,nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQAFunctionNorm(nIterations,nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQAMinGibbs(nIterations,nMQMQADiagnosticMaxAttempts), &
                dCapturedMQMQAMinimumBoundaryFraction(nIterations,nMQMQADiagnosticMaxAttempts))
            iCapturedMQMQAAssemblage = 0
            iCapturedMQMQAAcceptedCount = 0
            iCapturedMQMQAEligibleCount = 0
            iCapturedMQMQABoundaryCount = 0
            iCapturedMQMQAIterLast = 0
            iCapturedMQMQAPhaseDecision = 0
            dCapturedMQMQAPhaseDecision = 0D0
            dCapturedMQMQAPhaseMoles = 0D0
            dCapturedMQMQAElementPotential = 0D0
            lCapturedMQMQACandidateComparison = .FALSE.
            iCapturedMQMQACandidateDecision = 0
            dCapturedMQMQACandidateGamma = 0D0
            dCapturedMQMQACandidateDecision = 0D0
            iCapturedMQMQACandidateSolutionSpecies = 0
            dCapturedMQMQACandidateSolutionFraction = 0D0
            dCapturedMQMQACandidateAlpha = 0D0
            nCapturedMQMQARemovalEvents = 0
            iCapturedMQMQARemovalEvent = 0
            dCapturedMQMQARemovalEvent = 0D0
            dCapturedMQMQAAlpha = 0D0
            dCapturedMQMQAFunctionNorm = 0D0
            dCapturedMQMQAMinGibbs = 0D0
            dCapturedMQMQAMinimumBoundaryFraction = 1D0
        end if
        if (nCapturedMQMQAAttempts >= nMQMQADiagnosticMaxAttempts) return
        nCapturedMQMQAAttempts = nCapturedMQMQAAttempts+1
        iMQMQADiagnosticCurrentAttempt = nCapturedMQMQAAttempts

    end subroutine BeginMQMQATrajectoryAttempt


    !---------------------------------------------------------------------------------------------------------
    !> \brief Record the post-line-search, post-assemblage state for one global iteration.
    !---------------------------------------------------------------------------------------------------------
    subroutine CaptureMQMQATrajectoryPoint(iIteration,iLast,iAssemblage,dPhaseMoles,dElementPotentialState, &
        nAccepted,nEligible,nBoundary,dAlpha,dFunctionNorm,dMinGibbs,dMinimumBoundaryFraction)

        integer, intent(in) :: iIteration, iLast, iAssemblage(:), nAccepted, nEligible, nBoundary
        real(8), intent(in) :: dPhaseMoles(:), dElementPotentialState(:)
        real(8), intent(in) :: dAlpha, dFunctionNorm, dMinGibbs, dMinimumBoundaryFraction
        integer :: iAttempt

        if (.NOT. lCaptureMQMQATrajectory) return
        if (.NOT. allocated(iCapturedMQMQAAssemblage)) return
        iAttempt = iMQMQADiagnosticCurrentAttempt
        if ((iAttempt <= 0) .OR. (iAttempt > nMQMQADiagnosticMaxAttempts)) return
        if ((iIteration <= 0) .OR. (iIteration > SIZE(iCapturedMQMQAAcceptedCount,1))) return
        if (SIZE(iAssemblage) /= SIZE(iCapturedMQMQAAssemblage,1)) return
        if (SIZE(dPhaseMoles) < SIZE(iAssemblage)) return
        if (SIZE(dElementPotentialState) /= SIZE(iAssemblage)) return
        iCapturedMQMQAAssemblage(:,iIteration,iAttempt) = iAssemblage
        dCapturedMQMQAPhaseMoles(:,iIteration,iAttempt) = dPhaseMoles(1:SIZE(iAssemblage))
        dCapturedMQMQAElementPotential(:,iIteration,iAttempt) = dElementPotentialState
        iCapturedMQMQAAcceptedCount(iIteration,iAttempt) = nAccepted
        iCapturedMQMQAEligibleCount(iIteration,iAttempt) = nEligible
        iCapturedMQMQABoundaryCount(iIteration,iAttempt) = nBoundary
        iCapturedMQMQAIterLast(iIteration,iAttempt) = iLast
        dCapturedMQMQAAlpha(iIteration,iAttempt) = dAlpha
        dCapturedMQMQAFunctionNorm(iIteration,iAttempt) = dFunctionNorm
        dCapturedMQMQAMinGibbs(iIteration,iAttempt) = dMinGibbs
        dCapturedMQMQAMinimumBoundaryFraction(iIteration,iAttempt) = dMinimumBoundaryFraction
        iCapturedMQMQAPoints(iAttempt) = MAX(iCapturedMQMQAPoints(iAttempt),iIteration)

    end subroutine CaptureMQMQATrajectoryPoint


    !---------------------------------------------------------------------------------------------------------
    !> \brief Record the exact pure/solution driving-force ranking used for one phase-addition decision.
    !>
    !> \details The signed margin is pure minus solution.  A negative value therefore means the pure phase
    !!          is considered first by CheckPhaseAssemblage.  This routine is diagnostic-only and inactive
    !!          unless the existing MQMQA trajectory capture has been requested.
    !---------------------------------------------------------------------------------------------------------
    subroutine CaptureMQMQAPhaseDecision(iIteration,iPureSpecies,dPureForce,iSolutionPhase,dSolutionForce)

        integer, intent(in) :: iIteration, iPureSpecies, iSolutionPhase
        real(8), intent(in) :: dPureForce, dSolutionForce
        integer :: iAttempt

        if (.NOT. lCaptureMQMQATrajectory) return
        if (.NOT. allocated(iCapturedMQMQAPhaseDecision)) return
        iAttempt = iMQMQADiagnosticCurrentAttempt
        if ((iAttempt <= 0) .OR. (iAttempt > nMQMQADiagnosticMaxAttempts)) return
        if ((iIteration <= 0) .OR. (iIteration > SIZE(iCapturedMQMQAPhaseDecision,2))) return
        iCapturedMQMQAPhaseDecision(:,iIteration,iAttempt) = &
            [1,iPureSpecies,iSolutionPhase]
        dCapturedMQMQAPhaseDecision(:,iIteration,iAttempt) = &
            [dPureForce,dSolutionForce,dPureForce-dSolutionForce]

    end subroutine CaptureMQMQAPhaseDecision


    !---------------------------------------------------------------------------------------------------------
    !> \brief Record alpha-zero and selected corrected phase rankings from one identical pre-step state.
    !---------------------------------------------------------------------------------------------------------
    subroutine CaptureMQMQACandidateComparison(iIteration,dAlpha,dGammaCurrent,dGammaBase,dGammaCorrected, &
        iPureBase,dPureBase,iSolutionBase,dSolutionBase,iPureCorrected,dPureCorrected, &
        iSolutionCorrected,dSolutionCorrected,iSpeciesBase,dFractionBase,iSpeciesCorrected,dFractionCorrected)

        integer, intent(in) :: iIteration, iPureBase, iSolutionBase, iPureCorrected, iSolutionCorrected
        real(8), intent(in) :: dAlpha, dGammaCurrent(:), dGammaBase(:), dGammaCorrected(:)
        real(8), intent(in) :: dPureBase, dSolutionBase, dPureCorrected, dSolutionCorrected
        integer, intent(in) :: iSpeciesBase(:), iSpeciesCorrected(:)
        real(8), intent(in) :: dFractionBase(:), dFractionCorrected(:)
        integer :: iAttempt, nCopy

        if (.NOT. lCaptureMQMQATrajectory) return
        if (.NOT. allocated(dCapturedMQMQACandidateGamma)) return
        iAttempt = iMQMQADiagnosticCurrentAttempt
        if ((iAttempt <= 0) .OR. (iAttempt > nMQMQADiagnosticMaxAttempts)) return
        if ((iIteration <= 0) .OR. (iIteration > SIZE(dCapturedMQMQACandidateGamma,3))) return
        if (SIZE(dGammaCurrent) /= SIZE(dCapturedMQMQACandidateGamma,1)) return
        if ((SIZE(dGammaBase) /= SIZE(dGammaCurrent)) .OR. &
            (SIZE(dGammaCorrected) /= SIZE(dGammaCurrent))) return

        lCapturedMQMQACandidateComparison(iIteration,iAttempt) = .TRUE.
        dCapturedMQMQACandidateAlpha(iIteration,iAttempt) = dAlpha
        dCapturedMQMQACandidateGamma(:,1,iIteration,iAttempt) = dGammaCurrent
        dCapturedMQMQACandidateGamma(:,2,iIteration,iAttempt) = dGammaBase
        dCapturedMQMQACandidateGamma(:,3,iIteration,iAttempt) = dGammaCorrected
        iCapturedMQMQACandidateDecision(:,1,iIteration,iAttempt) = [iPureBase,iSolutionBase]
        iCapturedMQMQACandidateDecision(:,2,iIteration,iAttempt) = [iPureCorrected,iSolutionCorrected]
        dCapturedMQMQACandidateDecision(:,1,iIteration,iAttempt) = &
            [dPureBase,dSolutionBase,dPureBase-dSolutionBase]
        dCapturedMQMQACandidateDecision(:,2,iIteration,iAttempt) = &
            [dPureCorrected,dSolutionCorrected,dPureCorrected-dSolutionCorrected]
        nCopy = MIN(nMQMQADiagnosticMaxCandidateSpecies,SIZE(iSpeciesBase),SIZE(dFractionBase))
        iCapturedMQMQACandidateSolutionSpecies(1:nCopy,1,iIteration,iAttempt) = iSpeciesBase(1:nCopy)
        dCapturedMQMQACandidateSolutionFraction(1:nCopy,1,iIteration,iAttempt) = dFractionBase(1:nCopy)
        nCopy = MIN(nMQMQADiagnosticMaxCandidateSpecies,SIZE(iSpeciesCorrected),SIZE(dFractionCorrected))
        iCapturedMQMQACandidateSolutionSpecies(1:nCopy,2,iIteration,iAttempt) = iSpeciesCorrected(1:nCopy)
        dCapturedMQMQACandidateSolutionFraction(1:nCopy,2,iIteration,iAttempt) = dFractionCorrected(1:nCopy)

    end subroutine CaptureMQMQACandidateComparison


    !---------------------------------------------------------------------------------------------------------
    !> \brief Record one solution- or pure-phase removal attempt, result, swap, or reversion.
    !---------------------------------------------------------------------------------------------------------
    subroutine CaptureMQMQAPhaseRemovalEvent(iIteration,iPhase,iPhaseClass,iEventType,iOutcome, &
        dCurrent,dPrevious,dChange,dThreshold,dConditionLimit)

        integer, intent(in) :: iIteration, iPhase, iPhaseClass, iEventType, iOutcome
        real(8), intent(in) :: dCurrent, dPrevious, dChange, dThreshold, dConditionLimit
        integer :: iAttempt, iEvent

        if (.NOT. lCaptureMQMQATrajectory) return
        if (.NOT. allocated(nCapturedMQMQARemovalEvents)) return
        iAttempt = iMQMQADiagnosticCurrentAttempt
        if ((iAttempt <= 0) .OR. (iAttempt > nMQMQADiagnosticMaxAttempts)) return
        if ((iIteration <= 0) .OR. (iIteration > SIZE(nCapturedMQMQARemovalEvents,1))) return
        iEvent = nCapturedMQMQARemovalEvents(iIteration,iAttempt)+1
        if (iEvent > nMQMQADiagnosticMaxRemovalEvents) return
        nCapturedMQMQARemovalEvents(iIteration,iAttempt) = iEvent
        iCapturedMQMQARemovalEvent(:,iEvent,iIteration,iAttempt) = &
            [iPhase,iPhaseClass,iEventType,iOutcome]
        dCapturedMQMQARemovalEvent(:,iEvent,iIteration,iAttempt) = &
            [dCurrent,dPrevious,dChange,dThreshold,dConditionLimit]

    end subroutine CaptureMQMQAPhaseRemovalEvent


    !---------------------------------------------------------------------------------------------------------
    !> \brief Clear captured storage and restore the default-inactive diagnostic state.
    !---------------------------------------------------------------------------------------------------------
    subroutine ResetGEMNewtonDiagnosticCapture

        if (allocated(dCapturedGEMNewtonA)) deallocate(dCapturedGEMNewtonA)
        if (allocated(dCapturedGEMNewtonB)) deallocate(dCapturedGEMNewtonB)
        if (allocated(dCapturedGEMNewtonCorrectedA)) deallocate(dCapturedGEMNewtonCorrectedA)
        if (allocated(dCapturedGEMNewtonCorrectedB)) deallocate(dCapturedGEMNewtonCorrectedB)
        if (allocated(dMQMQADiagnosticNullModeSummary)) deallocate(dMQMQADiagnosticNullModeSummary)
        if (allocated(iMQMQADiagnosticNullModeDecisionChanges)) &
            deallocate(iMQMQADiagnosticNullModeDecisionChanges)
        if (allocated(dMQMQADiagnosticLeftNullSummary)) deallocate(dMQMQADiagnosticLeftNullSummary)
        if (allocated(dMQMQADiagnosticLeftNullVector)) deallocate(dMQMQADiagnosticLeftNullVector)
        if (allocated(dMQMQADiagnosticRightNullVector)) deallocate(dMQMQADiagnosticRightNullVector)
        if (allocated(iMQMQADiagnosticEquationIdentity)) deallocate(iMQMQADiagnosticEquationIdentity)
        if (allocated(dMQMQADiagnosticPhaseDependencyVector)) &
            deallocate(dMQMQADiagnosticPhaseDependencyVector)
        if (allocated(dMQMQADiagnosticPhaseDependencyProjection)) &
            deallocate(dMQMQADiagnosticPhaseDependencyProjection)
        if (allocated(dMQMQADiagnosticPhaseEquationData)) deallocate(dMQMQADiagnosticPhaseEquationData)
        if (allocated(iCapturedMQMQAAssemblage)) deallocate(iCapturedMQMQAAssemblage)
        if (allocated(iCapturedMQMQAAcceptedCount)) deallocate(iCapturedMQMQAAcceptedCount)
        if (allocated(iCapturedMQMQAEligibleCount)) deallocate(iCapturedMQMQAEligibleCount)
        if (allocated(iCapturedMQMQABoundaryCount)) deallocate(iCapturedMQMQABoundaryCount)
        if (allocated(iCapturedMQMQAIterLast)) deallocate(iCapturedMQMQAIterLast)
        if (allocated(iCapturedMQMQAPhaseDecision)) deallocate(iCapturedMQMQAPhaseDecision)
        if (allocated(dCapturedMQMQAPhaseDecision)) deallocate(dCapturedMQMQAPhaseDecision)
        if (allocated(dCapturedMQMQAPhaseMoles)) deallocate(dCapturedMQMQAPhaseMoles)
        if (allocated(dCapturedMQMQAElementPotential)) deallocate(dCapturedMQMQAElementPotential)
        if (allocated(lCapturedMQMQACandidateComparison)) deallocate(lCapturedMQMQACandidateComparison)
        if (allocated(iCapturedMQMQACandidateDecision)) deallocate(iCapturedMQMQACandidateDecision)
        if (allocated(dCapturedMQMQACandidateGamma)) deallocate(dCapturedMQMQACandidateGamma)
        if (allocated(dCapturedMQMQACandidateDecision)) deallocate(dCapturedMQMQACandidateDecision)
        if (allocated(iCapturedMQMQACandidateSolutionSpecies)) &
            deallocate(iCapturedMQMQACandidateSolutionSpecies)
        if (allocated(dCapturedMQMQACandidateSolutionFraction)) &
            deallocate(dCapturedMQMQACandidateSolutionFraction)
        if (allocated(dCapturedMQMQACandidateAlpha)) deallocate(dCapturedMQMQACandidateAlpha)
        if (allocated(nCapturedMQMQARemovalEvents)) deallocate(nCapturedMQMQARemovalEvents)
        if (allocated(iCapturedMQMQARemovalEvent)) deallocate(iCapturedMQMQARemovalEvent)
        if (allocated(dCapturedMQMQARemovalEvent)) deallocate(dCapturedMQMQARemovalEvent)
        if (allocated(dCapturedMQMQAAlpha)) deallocate(dCapturedMQMQAAlpha)
        if (allocated(dCapturedMQMQAFunctionNorm)) deallocate(dCapturedMQMQAFunctionNorm)
        if (allocated(dCapturedMQMQAMinGibbs)) deallocate(dCapturedMQMQAMinGibbs)
        if (allocated(dCapturedMQMQAMinimumBoundaryFraction)) &
            deallocate(dCapturedMQMQAMinimumBoundaryFraction)
        lCaptureGEMNewtonSystem = .FALSE.
        lCaptureGEMNewtonCorrectedSystem = .FALSE.
        lCaptureFirstGEMNewtonCorrectionPair = .FALSE.
        lMQMQADiagnosticMinimumNormStudy = .FALSE.
        lMQMQADiagnosticMinimumNormCaptured = .FALSE.
        iMQMQADiagnosticMinimumNormRank = 0
        iMQMQADiagnosticMinimumNormInfo = 0
        iMQMQADiagnosticMinimumNormDecision = 0
        dMQMQADiagnosticMinimumNormRHSResidual = 0D0
        dMQMQADiagnosticMinimumNormResidual = 0D0
        dMQMQADiagnosticMinimumNormSolutionNorm = 0D0
        dMQMQADiagnosticMinimumNormRelativeDifference = 0D0
        dMQMQADiagnosticMinimumNormNullResidual = 0D0
        dMQMQADiagnosticMinimumNormGammaDifference = 0D0
        dMQMQADiagnosticMinimumNormDecisionForce = 0D0
        dMQMQADiagnosticMinimumNormCrossSystemDifference = 0D0
        dMQMQADiagnosticMinimumNormCrossGammaDifference = 0D0
        dMQMQADiagnosticNullSpaceComparison = 0D0
        iMQMQADiagnosticPhaseStructureRank = 0
        iMQMQADiagnosticPhaseRuleSummary = 0
        dMQMQADiagnosticPhaseStructureSummary = 0D0
        nMQMQADiagnosticPhaseChangeChecks = 0
        iMQMQADiagnosticPhaseChangeCheck = 0
        iMQMQADiagnosticPhaseChangeAssemblage = 0
        dMQMQADiagnosticPhaseChangeCheck = 0D0
        dMQMQADiagnosticMinimumNormSingularSummary = 0D0
        iMQMQADiagnosticConstrainedInfo = 0
        iMQMQADiagnosticConstrainedDecision = 0
        dMQMQADiagnosticConstrainedRHSResidual = 0D0
        dMQMQADiagnosticConstrainedResidual = 0D0
        dMQMQADiagnosticConstrainedSolutionNorm = 0D0
        dMQMQADiagnosticConstrainedGroupNorm = 0D0
        dMQMQADiagnosticConstrainedDecisionForce = 0D0
        lGEMNewtonSystemCaptured = .FALSE.
        lGEMNewtonCorrectedSystemCaptured = .FALSE.
        lMQMQADiagnosticForceZeroAlpha = .FALSE.
        lMQMQADiagnosticRejectAcceptedCandidate = .FALSE.
        lMQMQADiagnosticPhasePathStudy = .FALSE.
        nMQMQADiagnosticPhasePathCandidates = 0
        nMQMQADiagnosticLeadingIdentityChanges = 0
        nMQMQADiagnosticEligibilityCrossings = 0
        nMQMQADiagnosticOrderingReversals = 0
        nMQMQADiagnosticRemovalCrossings = 0
        dMQMQADiagnosticMaxScaledForceShift = 0D0
        dMQMQADiagnosticMaxActiveAmountDisplacement = 0D0
        lCaptureMQMQATrajectory = .FALSE.
        iMQMQADiagnosticCurrentAttempt = 0
        iMQMQADiagnosticSkipAttempt = -1
        iMQMQADiagnosticSkipIteration = -1
        nCapturedMQMQAAttempts = 0
        iCapturedMQMQAPoints = 0
        nCapturedGEMNewtonVariables = 0

    end subroutine ResetGEMNewtonDiagnosticCapture

end module ModuleGEMNewtonDiagnosticCapture
