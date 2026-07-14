
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

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo, dTemperature
    USE ModuleGEMSolver

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

        ! Optionally validate mapped RKMP second-order terms without changing the Newton matrix:
        if (lDebugRKMPHessianFD) then
            call RKMPMappedHessianDiagnostic
            call RKMPResponseDiagnostic
            call RKMPGEMIdealReconstructionDiagnostic(A, B, nVar)
        end if

        ! Call the linear equation solver:
        if ((nConPhases > 1) .OR. (nSolnPhases > 0)) then
            if (lUseRKMPExactHessian) then
                call SolveRKMPAlphaTrust(A, B, nVar, IPIV, INFO)
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

    !> \brief Select the largest locally trustworthy RKMP response correction.
    !!
    !> \details The alpha-zero solve remains the reference direction.  RKMP curvature is withheld while the
    !! solver is finding a basin: a feasible Gibbs minimum must already have been recorded, the assemblage must
    !! be settled, and the previous nonlinear step must have maintained residual progress.  Once locally ready,
    !! candidates are tried in descending alpha order.  Floating-point validity, an emergency correction-ratio
    !! guard, DGESV success, update size, and direction agreement with the alpha-zero solve are checked.  The
    !! accepted correction is replayed once with metrics enabled.
    subroutine SolveRKMPAlphaTrust(AIn, BIn, nLocalVar, IPIVIn, INFOOut)

        integer, intent(in)                    :: nLocalVar
        integer, intent(out)                   :: INFOOut
        integer, dimension(:)                  :: IPIVIn
        real(8), dimension(:,:)                :: AIn
        real(8), dimension(:)                  :: BIn

        integer                                :: iAlpha, iLocal, INFOBase, INFOTrial
        integer, dimension(:), allocatable     :: IPIVTrial
        real(8)                                :: dAlphaCandidate, dBestAlpha, dNormBase, dNormTrial
        real(8)                                :: dTrialRatio, dBestUpdateRatio
        real(8)                                :: dCurrentGibbs, dGibbsScale, dDirectionCosine, dDirectionDifference
        real(8)                                :: dNormBase2, dNormTrial2
        real(8), parameter                     :: dEmergencyRatioCap = 1D6
        real(8), parameter                     :: dUpdateRatioCap = 1.25D0
        real(8), parameter                     :: dDirectionCosineMin = 0.90D0
        real(8), parameter                     :: dDirectionDifferenceCap = 0.50D0
        real(8), parameter                     :: dLocalNormThreshold = 5D-2
        real(8), parameter                     :: dProgressAllowance = 1.05D0
        real(8), parameter                     :: dGibbsActivationTolerance = 1D-6
        real(8), parameter                     :: dGibbsRetentionTolerance = 1D-4
        real(8), dimension(5)                  :: dAlphaList
        real(8), dimension(:), allocatable     :: BBase, BTrial, BZero
        real(8), dimension(:,:), allocatable   :: ABase, ATrial
        logical                                :: lCorrectionOK, lAccepted

        dAlphaList = [1D0, 1D-1, 1D-2, 1D-3, 0D0]
        dBestAlpha = 0D0
        dBestUpdateRatio = 0D0
        dDirectionCosine = 1D0
        dDirectionDifference = 0D0
        lAccepted = .FALSE.
        INFOOut = 0

        allocate(ABase(nLocalVar,nLocalVar), ATrial(nLocalVar,nLocalVar), &
                 BBase(nLocalVar), BTrial(nLocalVar), BZero(nLocalVar), IPIVTrial(nLocalVar))

        ABase = AIn
        BBase = BIn

        ATrial = ABase
        BTrial = BBase
        IPIVTrial = 0
        call dgesv(nLocalVar, 1, ATrial, nLocalVar, IPIVTrial, BTrial, nLocalVar, INFOBase)
        call CheckSolvedUpdate(BTrial, nLocalVar, INFOBase)
        if (INFOBase /= 0) then
            INFOOut = INFOBase
            deallocate(ABase, ATrial, BBase, BTrial, BZero, IPIVTrial)
            return
        end if

        BZero = BTrial
        dNormBase = DMAX1(MAXVAL(DABS(BTrial)), 1D-30)
        dNormBase2 = DMAX1(SQRT(SUM(BZero**2)), 1D-30)

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
                (dGEMFunctionNorm < dLocalNormThreshold) .AND. (iterGlobal - iterLast >= 5) .AND. &
                (dGEMFunctionNorm <= dProgressAllowance * DMAX1(dGEMFunctionNormLast,1D-30)) .AND. &
                (DABS(dCurrentGibbs - dMinGibbs) / dGibbsScale <= dGibbsActivationTolerance)
        else
            ! Retain local trust through small nonlinear oscillations, but return basin control to the ideal
            ! solve if residual or Gibbs behavior leaves the neighborhood where trust was established.
            lRKMPHessianNonlinearReady = (dGEMFunctionNorm < dProgressAllowance * dLocalNormThreshold) .AND. &
                (iterGlobal - iterLast >= 5) .AND. &
                (DABS(dCurrentGibbs - dMinGibbs) / dGibbsScale <= dGibbsRetentionTolerance)
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
            LOOP_ALPHA_TRUST: do iAlpha = 1, 5
                dAlphaCandidate = dAlphaList(iAlpha)

                ATrial = ABase
                BTrial = BBase
                lCorrectionOK = .TRUE.
                dTrialRatio = 0D0
                call MapRKMPHessianToGEMVariables(ATrial, BTrial, nLocalVar, dAlphaCandidate, &
                                                   .FALSE., lCorrectionOK, dTrialRatio)

                if (.NOT. lCorrectionOK) then
                    nRKMPHessianRejectBadDelta = nRKMPHessianRejectBadDelta + 1
                    cycle LOOP_ALPHA_TRUST
                end if

                ! This cap catches pathological scaling only; ordinary trust is based on nonlinear state and
                ! solved-direction behavior rather than the entrywise matrix-correction ratio.
                if (dTrialRatio > dEmergencyRatioCap) then
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

                dNormTrial = MAXVAL(DABS(BTrial))
                dBestUpdateRatio = dNormTrial / dNormBase
                if (dBestUpdateRatio > dUpdateRatioCap) then
                    nRKMPHessianRejectUpdate = nRKMPHessianRejectUpdate + 1
                    cycle LOOP_ALPHA_TRUST
                end if

                dNormTrial2 = DMAX1(SQRT(SUM(BTrial**2)), 1D-30)
                dDirectionCosine = DOT_PRODUCT(BZero,BTrial) / (dNormBase2*dNormTrial2)
                dDirectionDifference = SQRT(SUM((BTrial-BZero)**2)) / dNormBase2
                if ((dAlphaCandidate > 0D0) .AND. &
                    ((dDirectionCosine < dDirectionCosineMin) .OR. &
                     (dDirectionDifference > dDirectionDifferenceCap))) then
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

        AIn = ABase
        BIn = BBase
        lCorrectionOK = .TRUE.
        dTrialRatio = 0D0
        call MapRKMPHessianToGEMVariables(AIn, BIn, nLocalVar, dBestAlpha, .TRUE., lCorrectionOK, dTrialRatio)
        IPIVIn = 0
        call dgesv(nLocalVar, 1, AIn, nLocalVar, IPIVIn, BIn, nLocalVar, INFOOut)
        call CheckSolvedUpdate(BIn, nLocalVar, INFOOut)
        if ((INFOOut == 0) .AND. (dBestAlpha >= 1D0)) then
            nRKMPHessianFullAlphaCount = nRKMPHessianFullAlphaCount + 1
        end if

        deallocate(ABase, ATrial, BBase, BTrial, BZero, IPIVTrial)

    end subroutine SolveRKMPAlphaTrust


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
