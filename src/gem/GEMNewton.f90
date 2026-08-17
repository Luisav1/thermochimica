
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GEMNewton.f90
    !> \brief   Compute the direction vector for the GEMSolver using Newton's method.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      GEMSolver.f90
    !> \sa      GEMLineSearch.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code (new GEM solver)
    !   05/25/2012      M.H.A. Piro         Check for a NAN immediately after call to DGESV.
    !   01/31/2013      M.H.A. Piro         Check if a charged phase is contained in the database, but is
    !                                        not represented by the current phase assemblage.
    !   03/04/2013      M.H.A. Piro         Fix bug in correction process when dealing with ionic phases
    !                                        the loop should count back from the number of constraints,
    !                                        not the number of charged phases).
    !   09/06/2021      M. Poschmann        Correct moles of species to be proportional to mole fraction times
    !                                        moles of respective phase before direction vector is computed.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \brief The purpose of this subroutine is to compute the direction vector for the Gibbs energy
    !! minimization (GEM) solver using Newton's method.  The Hessian matrix and its corresponding constraint
    !! vector are first constructed and then the direction vector representing the system parameters is solved
    !! with the DGESV driver routine from LAPACK.  The updated element potentials, adjustments to the number of
    !! moles of solution phases and the number of moles of pure condensed phases are applied in the
    !! GEMLineSearch.f90 subroutine.
    !!
    !! Thermochimica is capable of handling ionic phases, which have an additional charge neutrality
    !! constraint imposed for each ionic phase.  Thus, an electron is added as a system component for every
    !! charged phase in the system.  It may be possible that an ionic phase is not predicted to be stable at
    !! a particular iteration and, thus, there aren't any stable species in the system representing that electron.
    !! To prevent a numerical singularity in the Hessian matrix, a check is performed after the Hessian matrix
    !! has been constructed ensuring that the Hessian does not contain a zero row.  In the event that the
    !! Hessian matrix contains all zeroes in the jth row (and necessarily, the jth column), a unit value is
    !! assigned to A(j,j).  Since the total balance of an electron is necessarily zero (i.e., ensuring charge
    !! neutrality) and there aren't any species for this solution phase, the corresponding value on the b vector
    !! will also be zero.  This procedure effectively ignores the jth row while preventing a numerical
    !! singularity.
    !
    !
    ! References:
    ! ===========
    !
    !> \details For further information regarding this methodology, refer to the following material:
    !! <ul>
    !! <li>  W.B. White, S.M. Johnson, G.B. Dantzig, "Chemical Equilibrium in Complex Mixtures," Journal of
    !!        Chemical Physics, V. 28, N. 5, 1958.
    !!
    !! <li>  G. Eriksson, "Thermodynamic Studies of High Temperature Equilibria," Acta Chemica Scandinavica,
    !!        25, 1971.
    !!
    !! <li>  G. Eriksson, E. Rosen, "General Equations for the Calculation of Equilibria in Multiphase Systems,"
    !!        Chemica Scripta, 4, 1973.
    !! </ul>
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[out]  INFO        An integer scalar used by LAPACK indicating a successful exit or an error.
    !
    ! nVar                      An integer scalar representing the total number of unknowns/linear equations.
    ! nElements                 An integer scalar representing the total number of elements in the system.
    ! nSpeciesPhase             An integer vector representing the number of species in a particular solution
    !                            phase (accumulative indexing)
    ! dStoichSpecies            A double real matrix representing stoichiometry coefficients.
    ! dMolesSpecies             A double real vector representing the number of moles of each species.
    ! dMolesPhase               A double real vector representing the number of moles of each phase.
    ! dMolesElement             A double real vector representing the number of moles of each element.
    ! JacobianLong              A double real matrix representing part of the Jacobian matrix that involves the
    !                            stoichiometry coefficients of solution species.
    ! JacobianShort             A double real vector that incorporates the JacobianLong matrix along with the
    !                            updated number of moles of each solution species.
    ! A                         Hessian matrix
    ! B                         Constraint vector (before call to LAPACK); unknown vector (after call to LAPACK)
    ! dEffStoichSolnPhase       A double real matrix representing the effective stoichiometry of a solution phase.
    ! dUpdateVar                A double real vector represending the updated system variables.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine GEMNewton(INFO)

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo, dTemperature
    USE ModuleGEMSolver
    USE ModuleGEMNewtonDiagnosticCapture, ONLY: CaptureGEMNewtonSystem
    USE ModuleMQMQAResponseMapping, ONLY: BuildActiveMQMQAGEMCorrection, SolveMQMQACorrectionTrial, &
        MQMQA_AGGREGATE_SUCCESS, MQMQA_AGGREGATE_NO_APPLICABLE_PHASE, MQMQA_AGGREGATE_PHASE_FAILURE, &
        MQMQA_MAP_OUTSIDE_INTERIOR, &
        MQMQA_TRIAL_ACCEPTED, MQMQA_TRIAL_APPLICATION_FALLBACK, MQMQA_TRIAL_DGESV_FALLBACK, &
        MQMQA_TRIAL_NONFINITE_FALLBACK, MQMQA_TRIAL_BASELINE_FAILURE
    USE ModuleMQMQATrust, ONLY: BuildMQMQAAlphaCandidateList, EvaluateMQMQACorrectionRatio, &
        EvaluateMQMQAUpdateTrust, &
        MQMQA_TRUST_ACCEPTED, MQMQA_TRUST_UPDATE_REJECTED, MQMQA_TRUST_DIRECTION_REJECTED

    implicit none

    integer                              :: i, j, k, l, m, INFO, nVar, iTry, nMaxTry
    integer, dimension(nElements)        :: iErrCol
    integer, dimension(:),   allocatable :: IPIV
    real(8)                              :: dTemp
    real(8), dimension(:),   allocatable :: B
    real(8), dimension(:,:), allocatable :: A

    ! Count phases:
    j = nConPhases
    nConPhases  = 0
    CountCon: do i = 1, j
        if (iAssemblage(i) > 0) then
            nConPhases  = nConPhases  + 1
        else
            exit CountCon
        end if
    end do CountCon

    j = nSolnPhases
    nSolnPhases = 0
    CountSoln: do i = nElements, nElements + 1 - j, -1
        if (iAssemblage(i) < 0) then
            nSolnPhases = nSolnPhases + 1
        else
            exit CountSoln
        end if
    end do CountSoln

    if ((nConPhases + nSolnPhases) <= 0) return

    ! Determine the number of unknowns/linear equations:
    nVar = nElements + nConPhases + nSolnPhases

    iErrCol = 0
    nMaxTry = nElements - (nConPhases + nSolnPhases)
    if (nMaxTry < 0) nMaxTry = 0
    TryLoop: do iTry = 0, nMaxTry
        ! on retry we are going to use dummy phases
        if (iTry > 0) nVar = nElements * 2

        ! Allocate memory:
        allocate(A(nVar, nVar))
        allocate(B(nVar))
        allocate(IPIV(nVar))

        ! Initialize variables:
        IPIV                = 0
        INFO                = 0
        A                   = 0D0
        B                   = 0D0
        dUpdateVar          = 0D0
        dEffStoichSolnPhase = 0D0

        do k = 1, nSolnPhases
            ! Absolute solution phase index:
            m = -iAssemblage(nElements - k + 1)
            ! Loop through species in phase:
            do l = nSpeciesPhase(m-1) + 1, nSpeciesPhase(m)
                dMolesSpecies(l) = dMolesPhase(nElements - k + 1) * dMolFraction(l)
                dMolesSpecies(l) = DMAX1(dMolesSpecies(l), dTolerance(8))
            end do
        end do

        ! Construct the Hessian matrix (elements):
        do j = 1, nElements
            do i = j, nElements
                do k = 1, nSolnPhases
                    ! Absolute solution phase index:
                    m = -iAssemblage(nElements - k + 1)
                    ! Loop through species in phase:
                    do l = nSpeciesPhase(m-1) + 1, nSpeciesPhase(m)
                        dTemp  = dStoichSpecies(l,i) * dStoichSpecies(l,j) * dMolesSpecies(l)
                        A(i,j) = A(i,j) + dTemp / (DFLOAT(iParticlesPerMole(l))**2)
                    end do
                end do
                ! Apply symmetry:
                A(j,i) = A(i,j)
            end do
        end do

        ! Compute the constraint vector (elements):
        do j = 1, nElements
            B(j) = dMolesElement(j)
            do l = 1, nSolnPhases
                k = -iAssemblage(nElements - l + 1)
                do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    dTemp = dStoichSpecies(i,j) * dMolesSpecies(i) * (dChemicalPotential(i) - 1D0)
                    B(j)  = B(j) + dTemp / DFLOAT(iParticlesPerMole(i))
                end do
            end do
        end do

        ! Construct the Hessian matrix and constraint vector (contribution from solution phases):
        do j = nElements + 1, nElements + nSolnPhases
            l = 2 * nElements - j + 1       ! Relative solution phase index (in iAssemblage vector).
            k = -iAssemblage(l)             ! Absolute solution phase index.

            ! Compute the stoichiometry of this phase:
            call CompStoichSolnPhase(k)

            do i = 1,nElements
                A(i,j) = dEffStoichSolnPhase(k,i) * dMolesPhase(l)
                A(j,i) = A(i,j)
            end do
            B(j) = dGibbsSolnPhase(k)
        end do

        ! Construct the Hessian matrix and constraint vector (contribution from pure condensed phases):
        do j = nElements + nSolnPhases + 1, nElements + nConPhases + nSolnPhases
            k = j - nElements - nSolnPhases
            do i = 1, nElements
                A(i,j) = dStoichSpecies(iAssemblage(k),i)
                A(j,i) = A(i,j)
            end do
            B(j) = dStdGibbsEnergy(iAssemblage(k))
        end do

        do k = 1, iTry
            i = iErrCol(k)
            j = nElements + nSolnPhases + nConPhases + k
            A(i,j) = 1D0
            A(j,i) = A(i,j)
            B(j) = 0D0
        end do

        ! Check if the Hessian is properly structured if the system contains any charged phases:
        if (nCountSublattice > 0) then
            ! Loop through elements
            LOOP_SUB: do j = nElements, nElements - nChargedConstraints + 1, -1
                dTemp = 0D0
                ! Loop through coefficients along column:
                do i = 1, nElements
                    dTemp = dTemp + DABS(A(i,j))
                    if (dTemp > 0D0) cycle LOOP_SUB
                end do
                ! The phase corresponding to this electron is not stable.
                A(j,j) = 1D0
            end do LOOP_SUB
        end if

        ! Verification tests may request the exact baseline system at this point. The capture routine is a
        ! default-inactive copy operation and runs before any experimental model-specific correction.
        call CaptureGEMNewtonSystem(A,B,nVar)

        ! Optionally verify mapped RKMP second-order terms without changing the Newton matrix:
        if (lDebugRKMPHessianFD) then
            call RKMPMappedHessianDiagnostic
            call RKMPResponseDiagnostic
            call RKMPGEMIdealReconstructionDiagnostic(A, B, nVar)
        end if

        ! Call the linear equation solver:
        if ((nConPhases > 1) .OR. (nSolnPhases > 0)) then
        
            ! The system is not purely elemental, so use the RKMP Hessian if it is active and enabled.
            if (lUseRKMPExactHessian .AND. lRKMPHessianActive .AND. &
                lUseMQMQAExactHessian .AND. (dMQMQAHessianAlpha > 0D0)) then
                if (lMQMQAHessianAdaptiveMode) then
                    call SolveMQMQAAlphaTrust(A,B,nVar,IPIV,INFO,.TRUE.)
                else
                    call SolveMQMQAFixedAlpha(A,B,nVar,IPIV,INFO,.TRUE.)
                end if
            else if (lUseRKMPExactHessian .AND. lRKMPHessianActive) then
                call SolveRKMPAlphaTrust(A, B, nVar, IPIV, INFO)
            else if (lUseMQMQAExactHessian .AND. (dMQMQAHessianAlpha > 0D0)) then
                if (lMQMQAHessianAdaptiveMode) then
                    call SolveMQMQAAlphaTrust(A,B,nVar,IPIV,INFO,.FALSE.)
                else
                    call SolveMQMQAFixedAlpha(A,B,nVar,IPIV,INFO,.FALSE.)
                end if
            else
                call dgesv( nVar, 1, A, nVar, IPIV, B, nVar, INFO )
            end if
        else
            do i = 1, nElements
                B(i) = dElementPotential(i)
            end do
            B(nElements + 1) = dMolesPhase(1)
        end if

        do k = 1, iTry
            j = nElements + nSolnPhases + nConPhases + k
            B(j) = 0D0
        end do

        ! Check for a NAN:
        LOOP_CheckNan: do i = 1, nVar
            if (B(i) /= B(i)) then
                INFO = 1
                exit LOOP_CheckNan
            end if
        end do LOOP_CheckNan

        if (iTry < nMaxTry) then
            if ((INFO <= 0) .OR. (INFO > nElements)) then
                exit TryLoop
            else
                iErrCol(iTry+1) = INFO
                INFO = 0
                deallocate(A, B, IPIV)
            end if
        end if
    end do TryLoop

    ! Store the updated variables if LAPACK is successful:
    if (INFO == 0) then
        do j = 1, nVar
            dUpdateVar(j) = B(j)
        end do

        ! Reset:
        lRevertSystem = .FALSE.
    else
        ! The system failed.  Revert to a previous assemblage.
        lRevertSystem = .TRUE.
        dUpdateVar    = 0D0
    end if

    ! Deallocate memory of local variables:
    i = 0
    deallocate(A, B, IPIV, STAT = i)
    if (i /= 0) INFOThermo = 24

    return

contains

    !> \brief Select the largest trustworthy MQMQA response correction while retaining an exact historical solve.
    !!
    !> \details The complete phase-local `(deltaA,deltaB)` aggregate is built once, while an untouched copy of
    !! the historical GEM system is solved as both the trust reference and exact fallback.  Alpha zero selects
    !! that historical Newton solve; it does not turn Thermochimica into a first-order method.  Readiness first
    !! decides whether the current assemblage and nonlinear trajectory are settled enough to try the correction.
    !! If ready, candidates are tested from the requested maximum downward, and the first candidate passing the
    !! correction, linear-solve, and grouped-update checks is selected.
    !!
    !! Element potentials, solution logarithmic increments, and pure-phase amounts are assessed as separate
    !! displacement groups because the raw GEM solution vector mixes different meanings and units.  The existing
    !! GEM line search remains the nonlinear step globalization.  Its observed residual and Gibbs progress feed
    !! the next call's readiness decision, so an apparently safe linear direction can still cause readiness to be
    !! revoked if the nonlinear calculation does not progress.
    subroutine SolveMQMQAAlphaTrust(AIn,BIn,nLocalVar,IPIVIn,INFOOut,lRKMPOwnsSolve)

        integer, intent(in) :: nLocalVar
        integer, intent(inout) :: IPIVIn(:)
        integer, intent(out) :: INFOOut
        real(8), intent(inout) :: AIn(:,:), BIn(:)
        logical, intent(in) :: lRKMPOwnsSolve

        integer :: iAggregateStatus, iAlpha, iFailurePhase, iFailureStatus, iHistorySlot, iInfoBase, iReadinessReason
        integer :: iLocal, iTrialStatus, iTrustStatus, nAccepted, nAlpha, nCharged, nSolutionStep
        integer, allocatable :: IPIVBase(:), IPIVTrial(:)
        logical :: lAccepted, lApplied, lEligible, lOldReady, lSupported, lTrustAccepted
        real(8) :: dAlphaCandidate, dAppliedNormA, dAppliedNormB, dBaseNormA, dBaseNormB
        real(8) :: dCurrentGibbs, dFailureMinimumFraction, dGibbsGap, dGibbsScale, dNormRatio, dRatioA, dRatioB
        real(8) :: dSelectedAlpha
        integer :: iCandidateRejectionMask(5)
        real(8) :: dAlphaList(5), dCandidateAlpha(5), dDirectionCosine(3), dDirectionDifference(3), dUpdateRatio(3)
        real(8) :: dBaseNormGroup(3), dTrialNormGroup(3)
        real(8), allocatable :: ABase(:,:), ABaseSolved(:,:), ATrial(:,:), BBase(:), BBaseSolved(:), BTrial(:)
        real(8), allocatable :: dDeltaA(:,:), dDeltaB(:), dStepBase(:), dStepTrial(:)

        INFOOut = 0

        ! MQ-4D step 1: build one complete correction pair for all currently eligible MQMQA phases.  The
        ! builder is transactional: a failure in any applicable phase rejects the aggregate rather than leaving
        ! a partially corrected GEM system.
        allocate(dDeltaA(nElements,nElements),dDeltaB(nElements))
        call BuildActiveMQMQAGEMCorrection(nSolnPhases,dDeltaA,dDeltaB,lSupported,lEligible, &
            nAccepted,nCharged,iFailurePhase,iFailureStatus,dFailureMinimumFraction,iAggregateStatus)
        lMQMQAHessianSupportedPhaseFound = lMQMQAHessianSupportedPhaseFound .OR. lSupported
        nMQMQAHessianChargedSkipCount = nMQMQAHessianChargedSkipCount+nCharged

        if ((iAggregateStatus /= MQMQA_AGGREGATE_SUCCESS) .AND. &
            (iAggregateStatus /= MQMQA_AGGREGATE_NO_APPLICABLE_PHASE)) then
            if (iFailureStatus == MQMQA_MAP_OUTSIDE_INTERIOR) then
                nMQMQAHessianInteriorFallbackCount = nMQMQAHessianInteriorFallbackCount+1
                dMQMQAHessianMinimumRejectedFraction = DMIN1( &
                    dMQMQAHessianMinimumRejectedFraction,dFailureMinimumFraction)
            else
                nMQMQAHessianAggregateFailureCount = nMQMQAHessianAggregateFailureCount+1
            end if
            lMQMQAHessianFallbackUsed = .TRUE.
            iMQMQAHessianLastFailurePhase = iFailurePhase
            if (iAggregateStatus == MQMQA_AGGREGATE_PHASE_FAILURE) then
                iMQMQAHessianLastFailureStatus = iFailureStatus
            else
                iMQMQAHessianLastFailureStatus = MQMQA_INTEGRATION_AGGREGATE_FAILURE
            end if
            if (lRKMPOwnsSolve) then
                call SolveRKMPAlphaTrust(AIn,BIn,nLocalVar,IPIVIn,INFOOut)
            else
                call dgesv(nLocalVar,1,AIn,nLocalVar,IPIVIn,BIn,nLocalVar,INFOOut)
            end if
            return
        end if
        if ((iAggregateStatus == MQMQA_AGGREGATE_NO_APPLICABLE_PHASE) .OR. (.NOT. lEligible)) then
            if (lRKMPOwnsSolve) then
                call SolveRKMPAlphaTrust(AIn,BIn,nLocalVar,IPIVIn,INFOOut)
            else
                call dgesv(nLocalVar,1,AIn,nLocalVar,IPIVIn,BIn,nLocalVar,INFOOut)
            end if
            return
        end if

        ! The RKMP and MQMQA experimental paths do not yet share a combined trust calculation.  When both are
        ! requested, preserve the established RKMP-owned solve and record the explicit ownership conflict.
        lMQMQAHessianEligibleCorrectionBuilt = .TRUE.
        lMQMQAHessianAggregateBuilt = .TRUE.
        nMQMQAHessianPhaseCorrectionCount = nMQMQAHessianPhaseCorrectionCount+nAccepted
        dMQMQAHessianMaxDeltaA = DMAX1(dMQMQAHessianMaxDeltaA,MAXVAL(DABS(dDeltaA)))
        dMQMQAHessianMaxDeltaB = DMAX1(dMQMQAHessianMaxDeltaB,MAXVAL(DABS(dDeltaB)))

        if (lRKMPOwnsSolve) then
            nMQMQAHessianRKMPConflictCount = nMQMQAHessianRKMPConflictCount+1
            lMQMQAHessianFallbackUsed = .TRUE.
            iMQMQAHessianLastFailurePhase = 0
            iMQMQAHessianLastFailureStatus = MQMQA_INTEGRATION_RKMP_CONFLICT
            call SolveRKMPAlphaTrust(AIn,BIn,nLocalVar,IPIVIn,INFOOut)
            return
        end if

        ! MQ-4D step 2: solve an untouched copy of the historical GEM system before trying any correction.  This
        ! solution is both the reference direction for the trust metrics and the exact alpha-zero fallback.
        allocate(ABase(nLocalVar,nLocalVar),ABaseSolved(nLocalVar,nLocalVar),ATrial(nLocalVar,nLocalVar), &
            BBase(nLocalVar),BBaseSolved(nLocalVar),BTrial(nLocalVar),IPIVBase(nLocalVar),IPIVTrial(nLocalVar), &
            dStepBase(nElements+nSpecies+nConPhases),dStepTrial(nElements+nSpecies+nConPhases))
        ABase = AIn
        BBase = BIn
        ABaseSolved = ABase
        BBaseSolved = BBase
        IPIVBase = 0
        call dgesv(nLocalVar,1,ABaseSolved,nLocalVar,IPIVBase,BBaseSolved,nLocalVar,iInfoBase)
        call CheckSolvedUpdate(BBaseSolved,nLocalVar,iInfoBase)
        if (iInfoBase /= 0) then
            INFOOut = iInfoBase
            return
        end if
        call BuildMQMQAGroupedDisplacement(BBaseSolved,nLocalVar,dStepBase,nSolutionStep)

        nMQMQAHessianEligibleSolveCount = nMQMQAHessianEligibleSolveCount+1
        iHistorySlot = nMQMQAHessianEligibleSolveCount
        iMQMQAHessianLastRejectionMask = 0
        dCandidateAlpha = -1D0
        iCandidateRejectionMask = 0
        dSelectedAlpha = 0D0
        dUpdateRatio = 1D0
        dDirectionCosine = 1D0
        dDirectionDifference = 0D0
        dBaseNormGroup = 0D0
        dTrialNormGroup = 0D0

        ! MQ-4D step 3: decide whether recent nonlinear behavior is settled enough to attempt a nonzero alpha.
        ! Readiness is trajectory state, not a property of the Hessian.  It may be revoked on a later call when
        ! the actual residual or Gibbs progress deteriorates after line search.
        dCurrentGibbs = 0D0
        do iLocal = 1,nElements
            dCurrentGibbs = dCurrentGibbs+dElementPotential(iLocal)*dMolesElement(iLocal)
        end do
        dCurrentGibbs = dCurrentGibbs*dTemperature*dIdealConstant
        dGibbsScale = DMAX1(DABS(dMinGibbs),1D0)
        dNormRatio = dGEMFunctionNorm/DMAX1(dGEMFunctionNormLast,1D-30)
        dGibbsGap = DABS(dCurrentGibbs-dMinGibbs)/dGibbsScale
        iReadinessReason = 0
        if (dMinGibbs >= 0.5D0*1D200) iReadinessReason = IOR(iReadinessReason,1)
        if (dGEMFunctionNorm >= dMQMQATrustLocalNormThreshold) iReadinessReason = IOR(iReadinessReason,2)
        if (iterGlobal-iterLast < iMQMQATrustSettledAssemblagePeriod) &
            iReadinessReason = IOR(iReadinessReason,4)
        if ((dGEMFunctionNorm > dMQMQATrustResolvedNormFloor) .AND. &
            (dNormRatio > dMQMQATrustProgressAllowance)) iReadinessReason = IOR(iReadinessReason,8)
        lOldReady = lMQMQAHessianNonlinearReady
        if (.NOT. lMQMQAHessianNonlinearReady) then
            if (dGibbsGap > dMQMQATrustGibbsActivationTolerance) &
                iReadinessReason = IOR(iReadinessReason,16)
            lMQMQAHessianNonlinearReady = (dMinGibbs < 0.5D0*1D200) .AND. &
                (dGEMFunctionNorm < dMQMQATrustLocalNormThreshold) .AND. &
                (iterGlobal-iterLast >= iMQMQATrustSettledAssemblagePeriod) .AND. &
                ((dGEMFunctionNorm <= dMQMQATrustResolvedNormFloor) .OR. &
                (dGEMFunctionNorm <= dMQMQATrustProgressAllowance*DMAX1(dGEMFunctionNormLast,1D-30))) .AND. &
                (DABS(dCurrentGibbs-dMinGibbs)/dGibbsScale <= dMQMQATrustGibbsActivationTolerance)
        else
            if (dGEMFunctionNorm >= dMQMQATrustProgressAllowance*dMQMQATrustLocalNormThreshold) &
                iReadinessReason = IOR(iReadinessReason,2)
            if (dGibbsGap > dMQMQATrustGibbsRetentionTolerance) &
                iReadinessReason = IOR(iReadinessReason,16)
            lMQMQAHessianNonlinearReady = &
                (dGEMFunctionNorm < dMQMQATrustProgressAllowance*dMQMQATrustLocalNormThreshold) .AND. &
                (iterGlobal-iterLast >= iMQMQATrustSettledAssemblagePeriod) .AND. &
                ((dGEMFunctionNorm <= dMQMQATrustResolvedNormFloor) .OR. &
                (dGEMFunctionNorm <= dMQMQATrustProgressAllowance*DMAX1(dGEMFunctionNormLast,1D-30))) .AND. &
                (DABS(dCurrentGibbs-dMinGibbs)/dGibbsScale <= dMQMQATrustGibbsRetentionTolerance)
        end if
        if ((.NOT. lOldReady) .AND. lMQMQAHessianNonlinearReady) &
            nMQMQAHessianReadinessActivationCount = nMQMQAHessianReadinessActivationCount+1
        if (lOldReady .AND. (.NOT. lMQMQAHessianNonlinearReady)) &
            nMQMQAHessianReadinessResetCount = nMQMQAHessianReadinessResetCount+1

        lAccepted = .FALSE.
        if (.NOT. lMQMQAHessianNonlinearReady) then
            nMQMQAHessianRejectNonlinear = nMQMQAHessianRejectNonlinear+1
            iMQMQAHessianLastRejectionMask = MQMQA_REJECT_NOT_READY
        else
            ! MQ-4D steps 4-5: test the largest permitted correction first.  A candidate must pass the
            ! correction-size check, corrected linear solve, finiteness checks, and grouped comparison with the
            ! historical update.  Rejection evidence is retained separately for every larger candidate.
            call BuildMQMQAAlphaCandidateList(dMQMQAHessianAlpha,dAlphaList,nAlpha)
            LOOP_MQMQA_ALPHA: do iAlpha = 1,nAlpha
                dAlphaCandidate = dAlphaList(iAlpha)
                dCandidateAlpha(iAlpha) = dAlphaCandidate
                if (dAlphaCandidate <= 0D0) exit LOOP_MQMQA_ALPHA

                call EvaluateMQMQACorrectionRatio(ABase,BBase,nElements,dDeltaA,dDeltaB,dAlphaCandidate, &
                    dMQMQATrustEmergencyRatioCap,lTrustAccepted,dRatioA,dRatioB)
                if (.NOT. lTrustAccepted) then
                    nMQMQAHessianRejectRatio = nMQMQAHessianRejectRatio+1
                    iMQMQAHessianLastRejectionMask = IOR(iMQMQAHessianLastRejectionMask,MQMQA_REJECT_RATIO)
                    iCandidateRejectionMask(iAlpha) = MQMQA_REJECT_RATIO
                    cycle LOOP_MQMQA_ALPHA
                end if

                call SolveMQMQACorrectionTrial(ABase,BBase,nElements,dDeltaA,dDeltaB,dAlphaCandidate, &
                    ATrial,BTrial,IPIVTrial,INFOOut,lApplied,lTrustAccepted,iTrialStatus)
                if (iTrialStatus == MQMQA_TRIAL_APPLICATION_FALLBACK) then
                    nMQMQAHessianRejectCorrection = nMQMQAHessianRejectCorrection+1
                    iMQMQAHessianLastRejectionMask = IOR(iMQMQAHessianLastRejectionMask,MQMQA_REJECT_CORRECTION)
                    iCandidateRejectionMask(iAlpha) = MQMQA_REJECT_CORRECTION
                    cycle LOOP_MQMQA_ALPHA
                else if (iTrialStatus == MQMQA_TRIAL_DGESV_FALLBACK) then
                    nMQMQAHessianRejectDGESV = nMQMQAHessianRejectDGESV+1
                    iMQMQAHessianLastRejectionMask = IOR(iMQMQAHessianLastRejectionMask,MQMQA_REJECT_DGESV)
                    iCandidateRejectionMask(iAlpha) = MQMQA_REJECT_DGESV
                    cycle LOOP_MQMQA_ALPHA
                else if (iTrialStatus == MQMQA_TRIAL_NONFINITE_FALLBACK) then
                    nMQMQAHessianRejectNonfinite = nMQMQAHessianRejectNonfinite+1
                    iMQMQAHessianLastRejectionMask = IOR(iMQMQAHessianLastRejectionMask,MQMQA_REJECT_NONFINITE)
                    iCandidateRejectionMask(iAlpha) = MQMQA_REJECT_NONFINITE
                    cycle LOOP_MQMQA_ALPHA
                else if ((iTrialStatus /= MQMQA_TRIAL_ACCEPTED) .OR. (.NOT. lTrustAccepted)) then
                    INFOOut = 1
                    return
                end if

                call BuildMQMQAGroupedDisplacement(BTrial,nLocalVar,dStepTrial,nSolutionStep)
                call EvaluateMQMQAUpdateTrust(dStepBase,dStepTrial,nElements,nSolutionStep,nConPhases, &
                    dMQMQATrustUpdateRatioCap,dMQMQATrustDirectionCosineMin, &
                    dMQMQATrustDirectionDifferenceCap,lTrustAccepted,dUpdateRatio,dDirectionCosine, &
                    dDirectionDifference,dBaseNormGroup,dTrialNormGroup,iTrustStatus)
                if (.NOT. lTrustAccepted) then
                    dMQMQAHessianLastRejectedAlpha = dAlphaCandidate
                    dMQMQAHessianRejectedGroupUpdateRatio = dUpdateRatio
                    dMQMQAHessianRejectedGroupDirectionCosine = dDirectionCosine
                    dMQMQAHessianRejectedGroupDirectionDifference = dDirectionDifference
                    dMQMQAHessianRejectedGroupBaseNorm = dBaseNormGroup
                    dMQMQAHessianRejectedGroupTrialNorm = dTrialNormGroup
                    if (iTrustStatus == MQMQA_TRUST_UPDATE_REJECTED) then
                        nMQMQAHessianRejectUpdate = nMQMQAHessianRejectUpdate+1
                        iMQMQAHessianLastRejectionMask = IOR(iMQMQAHessianLastRejectionMask,MQMQA_REJECT_UPDATE)
                        iCandidateRejectionMask(iAlpha) = MQMQA_REJECT_UPDATE
                    else if (iTrustStatus == MQMQA_TRUST_DIRECTION_REJECTED) then
                        nMQMQAHessianRejectDirection = nMQMQAHessianRejectDirection+1
                        iMQMQAHessianLastRejectionMask = IOR(iMQMQAHessianLastRejectionMask,MQMQA_REJECT_DIRECTION)
                        iCandidateRejectionMask(iAlpha) = MQMQA_REJECT_DIRECTION
                    else
                        nMQMQAHessianRejectNonfinite = nMQMQAHessianRejectNonfinite+1
                        iMQMQAHessianLastRejectionMask = IOR(iMQMQAHessianLastRejectionMask,MQMQA_REJECT_NONFINITE)
                        iCandidateRejectionMask(iAlpha) = MQMQA_REJECT_NONFINITE
                    end if
                    cycle LOOP_MQMQA_ALPHA
                end if

                dSelectedAlpha = dAlphaCandidate
                lAccepted = .TRUE.
                exit LOOP_MQMQA_ALPHA
            end do LOOP_MQMQA_ALPHA
        end if

        ! MQ-4D step 6: record what was selected and why larger candidates were rejected.  Histories use eligible
        ! solve slots and separately store the global Newton iteration because an initialization solve can occur
        ! at iteration zero.
        dMQMQAHessianSelectedAlpha = dSelectedAlpha
        dMQMQAHessianMaxSelectedAlpha = DMAX1(dMQMQAHessianMaxSelectedAlpha,dSelectedAlpha)
        dMQMQAHessianGroupUpdateRatio = dUpdateRatio
        dMQMQAHessianGroupDirectionCosine = dDirectionCosine
        dMQMQAHessianGroupDirectionDifference = dDirectionDifference
        dMQMQAHessianGroupBaseNorm = dBaseNormGroup
        dMQMQAHessianGroupTrialNorm = dTrialNormGroup
        if (iHistorySlot <= iterGlobalMax) then
            dMQMQAHessianAcceptedAlphaHistory(iHistorySlot) = dSelectedAlpha
            iMQMQAHessianGlobalIterationHistory(iHistorySlot) = iterGlobal
            iMQMQAHessianRejectionMaskHistory(iHistorySlot) = iMQMQAHessianLastRejectionMask
            dMQMQAHessianCandidateAlphaHistory(iHistorySlot,:) = dCandidateAlpha
            iMQMQAHessianCandidateRejectionMaskHistory(iHistorySlot,:) = iCandidateRejectionMask
            dMQMQAHessianFunctionNormHistory(iHistorySlot) = dGEMFunctionNorm
            dMQMQAHessianFunctionNormRatioHistory(iHistorySlot) = dNormRatio
            dMQMQAHessianGibbsGapHistory(iHistorySlot) = dGibbsGap
            if (lMQMQAHessianNonlinearReady) then
                iMQMQAHessianReadinessReasonHistory(iHistorySlot) = 0
            else
                iMQMQAHessianReadinessReasonHistory(iHistorySlot) = iReadinessReason
            end if
        end if

        ! MQ-4D step 7: return either the accepted corrected solve or the already solved historical arrays.
        ! GEMLineSearch subsequently globalizes this direction and its observed progress informs the next call's
        ! readiness decision.
        if (lAccepted) then
            AIn = ATrial
            BIn = BTrial
            IPIVIn = IPIVTrial
            INFOOut = 0
            lMQMQAHessianCorrectionApplied = .TRUE.
            lMQMQAHessianCorrectedSolveAccepted = .TRUE.
            nMQMQAHessianApplyCount = nMQMQAHessianApplyCount+1
            nMQMQAHessianAcceptedSolveCount = nMQMQAHessianAcceptedSolveCount+1
            dAppliedNormA = dSelectedAlpha*SQRT(SUM(dDeltaA**2))
            dAppliedNormB = dSelectedAlpha*SQRT(SUM(dDeltaB**2))
            dBaseNormA = SQRT(SUM(ABase(1:nElements,1:nElements)**2))
            dBaseNormB = SQRT(SUM(BBase(1:nElements)**2))
            dMQMQAHessianMaxAppliedA = DMAX1(dMQMQAHessianMaxAppliedA, &
                MAXVAL(DABS(dSelectedAlpha*dDeltaA)))
            dMQMQAHessianMaxAppliedB = DMAX1(dMQMQAHessianMaxAppliedB, &
                MAXVAL(DABS(dSelectedAlpha*dDeltaB)))
            dMQMQAHessianMaxRatioA = DMAX1(dMQMQAHessianMaxRatioA,dAppliedNormA/DMAX1(dBaseNormA,1D-30))
            dMQMQAHessianMaxRatioB = DMAX1(dMQMQAHessianMaxRatioB,dAppliedNormB/DMAX1(dBaseNormB,1D-30))
            if (dSelectedAlpha >= 1D0-1D-12) then
                nMQMQAHessianFullAlphaCount = nMQMQAHessianFullAlphaCount+1
                nMQMQAHessianFinalFullAlphaWindow = nMQMQAHessianFinalFullAlphaWindow+1
            else
                nMQMQAHessianReducedAlphaCount = nMQMQAHessianReducedAlphaCount+1
                nMQMQAHessianFinalFullAlphaWindow = 0
            end if
        else
            AIn = ABaseSolved
            BIn = BBaseSolved
            IPIVIn = IPIVBase
            INFOOut = 0
            nMQMQAHessianZeroAlphaCount = nMQMQAHessianZeroAlphaCount+1
            nMQMQAHessianFinalFullAlphaWindow = 0
        end if

    end subroutine SolveMQMQAAlphaTrust


    !> \brief Convert a solved GEM vector to the three displacement groups used by MQMQA trust.
    !!
    !> \details The solution-phase unknown stored by GEM is not itself a composition increment.  For each
    !! active solution constituent, the actual Newton-coordinate increment used by `GEMLineSearch` is the
    !! phase multiplier plus the stoichiometric element-potential target, less the current dimensionless
    !! chemical potential.  This quantity is the first-order logarithmic constituent-mole increment.  Building
    !! it constituent by constituent prevents a small phase multiplier from hiding a large composition update.
    subroutine BuildMQMQAGroupedDisplacement(BSolved,nLocalVar,dStep,nSolutionStep)

        integer, intent(in) :: nLocalVar
        integer, intent(out) :: nSolutionStep
        real(8), intent(in) :: BSolved(:)
        real(8), intent(out) :: dStep(:)

        integer :: iElement, iPhase, iSpecies, iSystemPhase, iWrite
        real(8) :: dIncrement

        dStep = 0D0
        dStep(1:nElements) = BSolved(1:nElements)-dElementPotential(1:nElements)
        nSolutionStep = 0
        do iPhase = 1,nSolnPhases
            iSystemPhase = -iAssemblage(nElements-iPhase+1)
            do iSpecies = nSpeciesPhase(iSystemPhase-1)+1,nSpeciesPhase(iSystemPhase)
                dIncrement = BSolved(nElements+iPhase)-dChemicalPotential(iSpecies)
                do iElement = 1,nElements
                    dIncrement = dIncrement+BSolved(iElement)*dStoichSpecies(iSpecies,iElement)/ &
                        DFLOAT(iParticlesPerMole(iSpecies))
                end do
                nSolutionStep = nSolutionStep+1
                dStep(nElements+nSolutionStep) = dIncrement
            end do
        end do

        do iPhase = 1,nConPhases
            iWrite = nElements+nSolutionStep+iPhase
            if ((iWrite > SIZE(dStep)) .OR. (nElements+nSolnPhases+iPhase > nLocalVar)) exit
            dStep(iWrite) = BSolved(nElements+nSolnPhases+iPhase)-dMolesPhase(iPhase)
        end do

    end subroutine BuildMQMQAGroupedDisplacement

    !> \brief Apply one fixed-alpha MQMQA aggregate transactionally and solve, with baseline fallback.
    !!
    !> \details `AIn/BIn` remain the only surviving historical baseline until a corrected trial has passed both
    !! application and `DGESV`. On a simultaneous eligible RKMP/MQMQA solve, the untouched baseline is passed
    !! directly to the existing RKMP trust routine and no MQMQA-modified trial can enter that path.
    subroutine SolveMQMQAFixedAlpha(AIn,BIn,nLocalVar,IPIVIn,INFOOut,lRKMPOwnsSolve)

        integer, intent(in) :: nLocalVar
        integer, intent(inout) :: IPIVIn(:)
        integer, intent(out) :: INFOOut
        real(8), intent(inout) :: AIn(:,:), BIn(:)
        logical, intent(in) :: lRKMPOwnsSolve

        integer :: iAggregateStatus, iFailurePhase, iFailureStatus, iTrialStatus
        integer :: nAccepted, nCharged
        integer, allocatable :: IPIVTrial(:)
        logical :: lAccepted, lApplied, lEligible, lSupported
        real(8) :: dAppliedNormA, dAppliedNormB, dBaseNormA, dBaseNormB, dFailureMinimumFraction
        real(8), allocatable :: ATrial(:,:), BTrial(:), dDeltaA(:,:), dDeltaB(:)

        INFOOut = 0
        allocate(dDeltaA(nElements,nElements),dDeltaB(nElements))
        call BuildActiveMQMQAGEMCorrection(nSolnPhases,dDeltaA,dDeltaB,lSupported,lEligible, &
            nAccepted,nCharged,iFailurePhase,iFailureStatus,dFailureMinimumFraction,iAggregateStatus)
        lMQMQAHessianSupportedPhaseFound = lMQMQAHessianSupportedPhaseFound .OR. lSupported
        nMQMQAHessianChargedSkipCount = nMQMQAHessianChargedSkipCount+nCharged

        if ((iAggregateStatus /= MQMQA_AGGREGATE_SUCCESS) .AND. &
            (iAggregateStatus /= MQMQA_AGGREGATE_NO_APPLICABLE_PHASE)) then
            if (iFailureStatus == MQMQA_MAP_OUTSIDE_INTERIOR) then
                nMQMQAHessianInteriorFallbackCount = nMQMQAHessianInteriorFallbackCount+1
                dMQMQAHessianMinimumRejectedFraction = DMIN1( &
                    dMQMQAHessianMinimumRejectedFraction,dFailureMinimumFraction)
            else
                nMQMQAHessianAggregateFailureCount = nMQMQAHessianAggregateFailureCount+1
            end if
            lMQMQAHessianFallbackUsed = .TRUE.
            iMQMQAHessianLastFailurePhase = iFailurePhase
            if (iAggregateStatus == MQMQA_AGGREGATE_PHASE_FAILURE) then
                iMQMQAHessianLastFailureStatus = iFailureStatus
            else
                iMQMQAHessianLastFailureStatus = MQMQA_INTEGRATION_AGGREGATE_FAILURE
            end if
            if (lRKMPOwnsSolve) then
                call SolveRKMPAlphaTrust(AIn,BIn,nLocalVar,IPIVIn,INFOOut)
            else
                call dgesv(nLocalVar,1,AIn,nLocalVar,IPIVIn,BIn,nLocalVar,INFOOut)
            end if
            return
        end if
        if ((iAggregateStatus == MQMQA_AGGREGATE_NO_APPLICABLE_PHASE) .OR. (.NOT. lEligible)) then
            if (lRKMPOwnsSolve) then
                call SolveRKMPAlphaTrust(AIn,BIn,nLocalVar,IPIVIn,INFOOut)
            else
                call dgesv(nLocalVar,1,AIn,nLocalVar,IPIVIn,BIn,nLocalVar,INFOOut)
            end if
            return
        end if

        lMQMQAHessianEligibleCorrectionBuilt = .TRUE.
        lMQMQAHessianAggregateBuilt = .TRUE.
        nMQMQAHessianPhaseCorrectionCount = nMQMQAHessianPhaseCorrectionCount+nAccepted
        dMQMQAHessianMaxDeltaA = DMAX1(dMQMQAHessianMaxDeltaA,MAXVAL(DABS(dDeltaA)))
        dMQMQAHessianMaxDeltaB = DMAX1(dMQMQAHessianMaxDeltaB,MAXVAL(DABS(dDeltaB)))

        if (lRKMPOwnsSolve) then
            nMQMQAHessianRKMPConflictCount = nMQMQAHessianRKMPConflictCount+1
            lMQMQAHessianFallbackUsed = .TRUE.
            iMQMQAHessianLastFailurePhase = 0
            iMQMQAHessianLastFailureStatus = MQMQA_INTEGRATION_RKMP_CONFLICT
            call SolveRKMPAlphaTrust(AIn,BIn,nLocalVar,IPIVIn,INFOOut)
            return
        end if

        allocate(ATrial(nLocalVar,nLocalVar),BTrial(nLocalVar),IPIVTrial(nLocalVar))
        call SolveMQMQACorrectionTrial(AIn,BIn,nElements,dDeltaA,dDeltaB,dMQMQAHessianAlpha, &
            ATrial,BTrial,IPIVTrial,INFOOut,lApplied,lAccepted,iTrialStatus)
        if (lApplied) then
            lMQMQAHessianCorrectionApplied = .TRUE.
            nMQMQAHessianApplyCount = nMQMQAHessianApplyCount+1
            dMQMQAHessianMaxAppliedA = DMAX1(dMQMQAHessianMaxAppliedA, &
                MAXVAL(DABS(dMQMQAHessianAlpha*dDeltaA)))
            dMQMQAHessianMaxAppliedB = DMAX1(dMQMQAHessianMaxAppliedB, &
                MAXVAL(DABS(dMQMQAHessianAlpha*dDeltaB)))
            dAppliedNormA = DSQRT(SUM((dMQMQAHessianAlpha*dDeltaA)**2))
            dAppliedNormB = DSQRT(SUM((dMQMQAHessianAlpha*dDeltaB)**2))
            dBaseNormA = DSQRT(SUM(AIn(1:nElements,1:nElements)**2))
            dBaseNormB = DSQRT(SUM(BIn(1:nElements)**2))
            dMQMQAHessianMaxRatioA = DMAX1(dMQMQAHessianMaxRatioA,dAppliedNormA/DMAX1(dBaseNormA,1D-30))
            dMQMQAHessianMaxRatioB = DMAX1(dMQMQAHessianMaxRatioB,dAppliedNormB/DMAX1(dBaseNormB,1D-30))
        end if

        select case(iTrialStatus)
        case(MQMQA_TRIAL_APPLICATION_FALLBACK)
            nMQMQAHessianApplicationFailureCount = nMQMQAHessianApplicationFailureCount+1
            iMQMQAHessianLastFailureStatus = MQMQA_INTEGRATION_APPLICATION_FAILURE
        case(MQMQA_TRIAL_DGESV_FALLBACK)
            nMQMQAHessianDGESVFallbackCount = nMQMQAHessianDGESVFallbackCount+1
            iMQMQAHessianLastFailureStatus = MQMQA_INTEGRATION_DGESV_FAILURE
        case(MQMQA_TRIAL_NONFINITE_FALLBACK)
            nMQMQAHessianNonfiniteFallbackCount = nMQMQAHessianNonfiniteFallbackCount+1
            iMQMQAHessianLastFailureStatus = MQMQA_INTEGRATION_NONFINITE_UPDATE
        case(MQMQA_TRIAL_BASELINE_FAILURE)
            iMQMQAHessianLastFailureStatus = MQMQA_INTEGRATION_DGESV_FAILURE
        end select

        AIn = ATrial
        BIn = BTrial
        IPIVIn = IPIVTrial
        if (iTrialStatus == MQMQA_TRIAL_ACCEPTED) then
            lMQMQAHessianCorrectedSolveAccepted = .TRUE.
            nMQMQAHessianAcceptedSolveCount = nMQMQAHessianAcceptedSolveCount+1
        else
            lMQMQAHessianFallbackUsed = .TRUE.
            iMQMQAHessianLastFailurePhase = 0
        end if

    end subroutine SolveMQMQAFixedAlpha

    !> \brief Select the largest locally trustworthy RKMP response correction.
    !!
    !> \details The alpha-zero solve remains the reference direction.  RKMP curvature is withheld while the
    !! solver is finding a basin: a feasible Gibbs minimum must already have been recorded, the assemblage must
    !! be settled, and the previous nonlinear step must have maintained residual progress.  Once locally ready,
    !! candidates are tried in descending alpha order.  Floating-point validity, an emergency correction-ratio
    !! guard, DGESV success, update size, and direction agreement with the alpha-zero solve are checked.  The
    !! accepted correction is replayed once with metrics enabled.  Alpha is a correction-blend value selected by
    !! the globalization/trust logic and applied to the fully constructed condensed GEM matrix and residual
    !! correction; it does not scale the RKMP derivatives or the local excess Hessian itself.
    subroutine SolveRKMPAlphaTrust(AIn, BIn, nLocalVar, IPIVIn, INFOOut)

        integer, intent(in)                    :: nLocalVar
        integer, intent(out)                   :: INFOOut
        integer, dimension(:)                  :: IPIVIn
        real(8), dimension(:,:)                :: AIn
        real(8), dimension(:)                  :: BIn

        integer                                :: iAlpha, iLocal, INFOBase, INFOTrial, nAlpha
        integer, dimension(:), allocatable     :: IPIVTrial
        real(8)                                :: dAlphaCandidate, dBestAlpha, dNormBase, dNormTrial
        real(8)                                :: dTrialRatio, dBestUpdateRatio
        real(8)                                :: dCurrentGibbs, dGibbsScale, dDirectionCosine, dDirectionDifference
        real(8)                                :: dNormBase2, dNormTrial2
        real(8), dimension(5)                  :: dAlphaList
        real(8), dimension(:), allocatable     :: BBase, BTrial, BZero, dStepTrial, dStepZero
        real(8), dimension(:,:), allocatable   :: ABase, ATrial
        logical                                :: lCorrectionOK, lAccepted

        call BuildAlphaCandidateList(dRKMPHessianBlendAlpha, dAlphaList, nAlpha)
        dBestAlpha = 0D0
        dBestUpdateRatio = 0D0
        dDirectionCosine = 1D0
        dDirectionDifference = 0D0
        lAccepted = .FALSE.
        INFOOut = 0

        allocate(ABase(nLocalVar,nLocalVar), ATrial(nLocalVar,nLocalVar), &
                 BBase(nLocalVar), BTrial(nLocalVar), BZero(nLocalVar), &
                 dStepTrial(nLocalVar), dStepZero(nLocalVar), IPIVTrial(nLocalVar))

        ABase = AIn
        BBase = BIn

        ATrial = ABase
        BTrial = BBase
        IPIVTrial = 0
        call dgesv(nLocalVar, 1, ATrial, nLocalVar, IPIVTrial, BTrial, nLocalVar, INFOBase)
        call CheckSolvedUpdate(BTrial, nLocalVar, INFOBase)
        if (INFOBase /= 0) then
            INFOOut = INFOBase
            deallocate(ABase, ATrial, BBase, BTrial, BZero, dStepTrial, dStepZero, IPIVTrial)
            return
        end if

        BZero = BTrial
        call BuildSolvedDisplacement(BZero, nLocalVar, dStepZero)
        dNormBase = DMAX1(MAXVAL(DABS(dStepZero)), 1D-30)
        dNormBase2 = DMAX1(SQRT(SUM(dStepZero**2)), 1D-30)

        dCurrentGibbs = 0D0
        do iLocal = 1, nElements
            dCurrentGibbs = dCurrentGibbs + dElementPotential(iLocal) * dMolesElement(iLocal)
        end do
        dCurrentGibbs = dCurrentGibbs * dTemperature * dIdealConstant
        dGibbsScale = DMAX1(DABS(dMinGibbs), 1D0)

        ! The ideal solve owns basin finding.  RKMP curvature becomes eligible only near a feasible Gibbs state
        ! after a settled, non-diverging nonlinear step.
        if (.NOT. lRKMPHessianNonlinearReady) then
            lRKMPHessianNonlinearReady = (dMinGibbs < 0.5D0 * 1D200) .AND. &
                (dGEMFunctionNorm < dRKMPTrustLocalNormThreshold) .AND. (iterGlobal - iterLast >= 5) .AND. &
                (dGEMFunctionNorm <= dRKMPTrustProgressAllowance * DMAX1(dGEMFunctionNormLast,1D-30)) .AND. &
                (DABS(dCurrentGibbs - dMinGibbs) / dGibbsScale <= dRKMPTrustGibbsActivationTolerance)
        else
            ! Retain local trust through small nonlinear oscillations, but return basin control to the ideal
            ! solve if residual or Gibbs behavior leaves the neighborhood where trust was established.
            lRKMPHessianNonlinearReady = &
                (dGEMFunctionNorm < dRKMPTrustProgressAllowance * dRKMPTrustLocalNormThreshold) .AND. &
                (iterGlobal - iterLast >= 5) .AND. &
                (DABS(dCurrentGibbs - dMinGibbs) / dGibbsScale <= dRKMPTrustGibbsRetentionTolerance)
        end if

        if ((dRKMPHessianBlendAlpha <= 0D0) .OR. (.NOT. lRKMPHessianNonlinearReady)) then
            dRKMPHessianSelectedAlpha = 0D0
            dRKMPHessianUpdateNormRatio = 1D0
            dRKMPHessianDirectionCosine = 1D0
            dRKMPHessianDirectionDifference = 0D0
            if ((dRKMPHessianBlendAlpha > 0D0) .AND. (.NOT. lRKMPHessianNonlinearReady)) then
                nRKMPHessianRejectNonlinear = nRKMPHessianRejectNonlinear + 1
            end if
        else
            LOOP_ALPHA_TRUST: do iAlpha = 1, nAlpha
                dAlphaCandidate = dAlphaList(iAlpha)

                ATrial = ABase
                BTrial = BBase
                lCorrectionOK = .TRUE.
                dTrialRatio = 0D0

                ! Apply the candidate RKMP correction to the GEM matrix and residual, but do not solve for the updated direction.
                call MapRKMPHessianToGEMVariables(ATrial, BTrial, nLocalVar, dAlphaCandidate, &
                                                   .FALSE., lCorrectionOK, dTrialRatio)

                if (.NOT. lCorrectionOK) then
                    nRKMPHessianRejectBadDelta = nRKMPHessianRejectBadDelta + 1
                    nRKMPHessianRejectLocalResponse = nRKMPHessianRejectLocalResponse + 1
                    cycle LOOP_ALPHA_TRUST
                end if

                ! This cap catches pathological scaling only; ordinary trust is based on nonlinear state and
                ! solved-direction behavior rather than the entrywise matrix-correction ratio.
                if (dTrialRatio > dRKMPTrustEmergencyRatioCap) then
                    nRKMPHessianRejectRatio = nRKMPHessianRejectRatio + 1
                    cycle LOOP_ALPHA_TRUST
                end if

                IPIVTrial = 0
                call dgesv(nLocalVar, 1, ATrial, nLocalVar, IPIVTrial, BTrial, nLocalVar, INFOTrial)
                call CheckSolvedUpdate(BTrial, nLocalVar, INFOTrial)
                if (INFOTrial /= 0) then
                    nRKMPHessianRejectDGESV = nRKMPHessianRejectDGESV + 1
                    cycle LOOP_ALPHA_TRUST
                end if

                call BuildSolvedDisplacement(BTrial, nLocalVar, dStepTrial)
                dNormTrial = MAXVAL(DABS(dStepTrial))
                dBestUpdateRatio = dNormTrial / dNormBase
                if (dBestUpdateRatio > dRKMPTrustUpdateRatioCap) then
                    nRKMPHessianRejectUpdate = nRKMPHessianRejectUpdate + 1
                    cycle LOOP_ALPHA_TRUST
                end if

                dNormTrial2 = DMAX1(SQRT(SUM(dStepTrial**2)), 1D-30)
                dDirectionCosine = DOT_PRODUCT(dStepZero,dStepTrial) / (dNormBase2*dNormTrial2)
                dDirectionDifference = SQRT(SUM((dStepTrial-dStepZero)**2)) / dNormBase2
                if ((dAlphaCandidate > 0D0) .AND. &
                    ((dDirectionCosine < dRKMPTrustDirectionCosineMin) .OR. &
                     (dDirectionDifference > dRKMPTrustDirectionDifferenceCap))) then
                    nRKMPHessianRejectDirection = nRKMPHessianRejectDirection + 1
                    cycle LOOP_ALPHA_TRUST
                end if

                dBestAlpha = dAlphaCandidate
                lAccepted = .TRUE.
                exit LOOP_ALPHA_TRUST
            end do LOOP_ALPHA_TRUST

            if (.NOT. lAccepted) then
                dBestAlpha = 0D0
                dBestUpdateRatio = 1D0
                dDirectionCosine = 1D0
                dDirectionDifference = 0D0
            end if

            dRKMPHessianSelectedAlpha = dBestAlpha
            dRKMPHessianUpdateNormRatio = dBestUpdateRatio
            dRKMPHessianDirectionCosine = dDirectionCosine
            dRKMPHessianDirectionDifference = dDirectionDifference
        end if
        dRKMPHessianMaxSelectedAlpha = DMAX1(dRKMPHessianMaxSelectedAlpha, dBestAlpha)
        if ((iterGlobal >= 1) .AND. (iterGlobal <= iterGlobalMax)) then
            dRKMPHessianAcceptedAlphaHistory(iterGlobal) = dBestAlpha
        end if

        AIn = ABase
        BIn = BBase
        lCorrectionOK = .TRUE.
        dTrialRatio = 0D0

        ! Finally apply the selected RKMP correction to the GEM matrix and residual, and solve for the updated direction.
        call MapRKMPHessianToGEMVariables(AIn, BIn, nLocalVar, dBestAlpha, .TRUE., lCorrectionOK, dTrialRatio)
        if (.NOT. lCorrectionOK) then
            nRKMPHessianRejectLocalResponse = nRKMPHessianRejectLocalResponse + 1
            AIn = ABase
            BIn = BZero
            IPIVIn = 0
            INFOOut = 0
            deallocate(ABase, ATrial, BBase, BTrial, BZero, dStepTrial, dStepZero, IPIVTrial)
            return
        end if
        IPIVIn = 0
        call dgesv(nLocalVar, 1, AIn, nLocalVar, IPIVIn, BIn, nLocalVar, INFOOut)
        call CheckSolvedUpdate(BIn, nLocalVar, INFOOut)
        if ((INFOOut == 0) .AND. (dBestAlpha >= 1D0)) then
            nRKMPHessianFullAlphaCount = nRKMPHessianFullAlphaCount + 1
        end if

        deallocate(ABase, ATrial, BBase, BTrial, BZero, dStepTrial, dStepZero, IPIVTrial)

    end subroutine SolveRKMPAlphaTrust


    !> \brief Build descending trust candidates bounded by the requested maximum alpha.
    subroutine BuildAlphaCandidateList(dAlphaMaxInput, dCandidates, nCandidates)

        real(8), intent(in)                  :: dAlphaMaxInput
        real(8), dimension(:), intent(out)   :: dCandidates
        integer, intent(out)                 :: nCandidates

        integer                              :: iCandidate
        real(8)                              :: dAlphaMax
        real(8), dimension(4)                :: dStandard

        dStandard = [1D0, 1D-1, 1D-2, 1D-3]
        dCandidates = 0D0
        nCandidates = 0
        dAlphaMax = DMAX1(0D0, DMIN1(1D0, dAlphaMaxInput))

        if (dAlphaMax > 0D0) then
            nCandidates = 1
            dCandidates(nCandidates) = dAlphaMax
            do iCandidate = 1, SIZE(dStandard)
                if (dStandard(iCandidate) >= dAlphaMax * (1D0 - 1D-12)) cycle
                nCandidates = nCandidates + 1
                dCandidates(nCandidates) = dStandard(iCandidate)
            end do
        end if
        nCandidates = nCandidates + 1
        dCandidates(nCandidates) = 0D0

    end subroutine BuildAlphaCandidateList


    !> \brief Convert a solved GEM variable vector into the displacement applied by the line search.
    !!
    !> \details Element potentials and pure condensed phase amounts are absolute targets, whereas solution-phase
    !! entries are incremental logarithmic multipliers.  Alpha trust must compare these displacements rather than
    !! the absolute DGESV output, whose common target values can hide a materially changed Newton direction.
    subroutine BuildSolvedDisplacement(BSolved, nLocalVar, dStep)

        integer, intent(in)                    :: nLocalVar
        real(8), dimension(:), intent(in)      :: BSolved
        real(8), dimension(:), intent(out)     :: dStep

        integer                                :: iLocal, iOffset

        dStep = BSolved
        dStep(1:nElements) = BSolved(1:nElements) - dElementPotential(1:nElements)

        iOffset = nElements + nSolnPhases
        do iLocal = 1, nConPhases
            if (iOffset + iLocal > nLocalVar) exit
            dStep(iOffset + iLocal) = BSolved(iOffset + iLocal) - dMolesPhase(iLocal)
        end do

    end subroutine BuildSolvedDisplacement


    subroutine CheckSolvedUpdate(BLocal, nLocalVar, INFOLocal)

        integer, intent(in)                  :: nLocalVar
        integer, intent(inout)               :: INFOLocal
        real(8), dimension(:), intent(in)    :: BLocal

        integer                              :: iLocal

        if (INFOLocal /= 0) return

        do iLocal = 1, nLocalVar
            if ((BLocal(iLocal) /= BLocal(iLocal)) .OR. &
                (DABS(BLocal(iLocal)) > 0.5D0 * HUGE(1D0))) then
                INFOLocal = 1
                return
            end if
        end do

    end subroutine CheckSolvedUpdate

end subroutine GEMNewton
