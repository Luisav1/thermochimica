!-------------------------------------------------------------------------------------------------------------
!> \file    RKMPResponseDiagnostic.f90
!> \brief   Diagnostic-only local constrained RKMP composition response.
!>
!> \details This diagnostic advances beyond projected curvature by asking the
!!          response question needed by GEM: how does the equilibrium composition
!!          of one fixed-amount RKMP phase change when element potentials are
!!          perturbed?
!!
!!          It combines ideal and excess mole-fraction curvature and solves a
!!          bordered system whose infinitesimal composition changes sum to zero,
!!          preserving mole-fraction normalization. It then compares the analytic
!!          chemical-potential response with finite differences of the established
!!          production RKMP partial-molar routine. It also reports the difference
!!          from the ideal-only condensed response.
!!
!!          The diagnostic does not modify A, B, the Newton update, or the phase
!!          assemblage. Its verified mathematics is implemented for production
!!          use by MapRKMPHessianToGEMVariables.
!-------------------------------------------------------------------------------------------------------------

subroutine RKMPResponseDiagnostic

    USE ModuleThermo
    USE ModuleGEMSolver, ONLY: dPartialExcessGibbs, lDebugRKMPHessianFD

    implicit none

    interface
        subroutine CompExcessGibbsEnergyRKMP_unconstrained(iSolnIndex,dHess)
            integer, intent(in)                  :: iSolnIndex
            real(8), intent(out), dimension(:,:) :: dHess
        end subroutine CompExcessGibbsEnergyRKMP_unconstrained

        subroutine CompExcessGibbsEnergyRKMP(iSolnIndex)
            integer :: iSolnIndex
        end subroutine CompExcessGibbsEnergyRKMP
    end interface

    integer                              :: i, j, k, m, p, iStep, nLocalSpecies, iFirstSpecies, iLastSpecies
    integer                              :: INFO
    real(8)                              :: dTotalMoles, dEps, dFDScale, dNorm, dMaxCurvEx, dMaxCurvIdeal
    real(8)                              :: dFDResponseRelative, dMaxResponse, dMaxFDResponseErr
    real(8)                              :: dMaxCandidateDelta
    real(8), dimension(:), allocatable   :: dX, dXBase, dGammaDir, dDX, dMuPlus, dMuMinus, dMuAnalytic
    real(8), dimension(:), allocatable   :: dMolFractionSave, dPartialExcessSave
    real(8), dimension(:,:), allocatable :: dHloc, dHx, dC, dResponse, dCandidate, dIdealCandidate

    if (.NOT. lDebugRKMPHessianFD) return

    !=========================================================================================================
    ! SECTION 1: BUILD AND SOLVE THE CORRECTED LOCAL RESPONSE
    !=========================================================================================================
    do k = 1, nSolnPhases
        m = -iAssemblage(nElements - k + 1)
        if (m <= 0) cycle
        if (cSolnPhaseType(m) /= 'RKMP') cycle

        iFirstSpecies = nSpeciesPhase(m-1) + 1
        iLastSpecies  = nSpeciesPhase(m)
        nLocalSpecies = iLastSpecies - iFirstSpecies + 1
        if (nLocalSpecies <= 1) cycle

        allocate(dX(nLocalSpecies), dXBase(nLocalSpecies), dGammaDir(nElements), dDX(nLocalSpecies), &
                 dMuPlus(nLocalSpecies), dMuMinus(nLocalSpecies), dMuAnalytic(nLocalSpecies), &
                 dMolFractionSave(nLocalSpecies), dPartialExcessSave(nLocalSpecies), &
                 dHloc(nLocalSpecies,nLocalSpecies), dHx(nLocalSpecies,nLocalSpecies), &
                 dC(nLocalSpecies,nElements), dResponse(nLocalSpecies,nElements), &
                 dCandidate(nElements,nElements), dIdealCandidate(nElements,nElements))

        dTotalMoles = SUM(dMolesSpecies(iFirstSpecies:iLastSpecies))
        if (dTotalMoles <= 1D-30) then
            deallocate(dX, dXBase, dGammaDir, dDX, dMuPlus, dMuMinus, dMuAnalytic, &
                       dMolFractionSave, dPartialExcessSave, dHloc, dHx, dC, &
                       dResponse, dCandidate, dIdealCandidate)
            cycle
        end if

        dXBase = dMolesSpecies(iFirstSpecies:iLastSpecies) / dTotalMoles
        dX = DMAX1(dXBase, 1D-30)
        dX = dX / SUM(dX)

        do j = 1, nElements
            do i = 1, nLocalSpecies
                p = iFirstSpecies + i - 1
                dC(i,j) = dStoichSpecies(p,j) / DFLOAT(iParticlesPerMole(p))
            end do
        end do

        call CompExcessGibbsEnergyRKMP_unconstrained(m, dHloc)

        ! Convert mole-number curvature to fixed-phase-amount mole-fraction
        ! curvature. If total phase moles N remain fixed, a mole-fraction change
        ! dx corresponds to a mole change N*dx. Multiplying Hloc by N therefore
        ! predicts the resulting excess chemical-potential change.
        dHx = dTotalMoles * dHloc
        dMaxCurvEx = MAXVAL(DABS(dHx))
        dMaxCurvIdeal = MAXVAL(1D0 / dX)

        do i = 1, nLocalSpecies
            dHx(i,i) = dHx(i,i) + 1D0 / dX(i)
        end do

        call SolveRKMPResponse(nLocalSpecies, nElements, dHx, dC, dResponse, INFO)
        if (INFO /= 0) then
            write(*,'(A,1X,A,1X,I0,1X,A,1X,I0)') &
                'RKMP_RESPONSE_DEBUG', TRIM(cSolnPhaseName(m)), m, 'solve_info=', INFO
            deallocate(dX, dXBase, dGammaDir, dDX, dMuPlus, dMuMinus, dMuAnalytic, &
                       dMolFractionSave, dPartialExcessSave, dHloc, dHx, dC, &
                       dResponse, dCandidate, dIdealCandidate)
            cycle
        end if

        dCandidate = MATMUL(TRANSPOSE(dC), dTotalMoles * dResponse)
        call SolveIdealResponse(nLocalSpecies, nElements, dX, dC, dIdealCandidate, INFO)
        dMaxCandidateDelta = MAXVAL(DABS(dCandidate - dIdealCandidate))
        dMaxResponse = MAXVAL(DABS(dResponse))

        !=====================================================================================================
        ! SECTION 2: PRODUCTION PARTIAL-MOLAR FINITE-DIFFERENCE CHECK
        !
        ! Verify the predicted excess chemical-potential change along a
        ! normalized composition direction generated by the constrained response.
        !=====================================================================================================
        ! Perturb production mole fractions along that direction and compare
        ! finite-difference partial molars with the analytic prediction.
        dGammaDir = 0D0
        do j = 1, nElements
            dGammaDir(j) = 1D0 / DFLOAT(j)
        end do
        dDX = MATMUL(dResponse, dGammaDir)
        dNorm = MAXVAL(DABS(dDX))
        if (dNorm > 0D0) dDX = dDX / dNorm

        ! Sweep the established production partial-molar RKMP routine over notebook-style perturbations.
        ! Reduce the common scale only when a dilute component would otherwise become non-positive.
        dFDScale = 1D0
        do i = 1, nLocalSpecies
            if (DABS(dDX(i)) > 0D0) then
                dFDScale = DMIN1(dFDScale, 25D0 * dX(i) / DABS(dDX(i)))
            end if
        end do

        dMolFractionSave = dMolFraction(iFirstSpecies:iLastSpecies)
        dPartialExcessSave = dPartialExcessGibbs(iFirstSpecies:iLastSpecies)
        dMuAnalytic = MATMUL(dTotalMoles * dHloc, dDX)

        do iStep = 2, 6
            dEps = dFDScale * 10D0**(-iStep)

            dMolFraction(iFirstSpecies:iLastSpecies) = dX + dEps * dDX
            dPartialExcessGibbs(iFirstSpecies:iLastSpecies) = 0D0
            call CompExcessGibbsEnergyRKMP(m)
            dMuPlus = dPartialExcessGibbs(iFirstSpecies:iLastSpecies)

            dMolFraction(iFirstSpecies:iLastSpecies) = dX - dEps * dDX
            dPartialExcessGibbs(iFirstSpecies:iLastSpecies) = 0D0
            call CompExcessGibbsEnergyRKMP(m)
            dMuMinus = dPartialExcessGibbs(iFirstSpecies:iLastSpecies)

            dMaxFDResponseErr = MAXVAL(DABS(((dMuPlus-dMuMinus)/(2D0*dEps))-dMuAnalytic))
            dFDResponseRelative = dMaxFDResponseErr / DMAX1(MAXVAL(DABS(dMuAnalytic)),1D-30)

            write(*,'(A,1X,A,1X,I0,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6)') &
                'RKMP_PRODUCTION_MU_FD', TRIM(cSolnPhaseName(m)), m, 'eps=', dEps, &
                'analyticMax=', MAXVAL(DABS(dMuAnalytic)), 'absErr=', dMaxFDResponseErr, &
                'relerr=', dFDResponseRelative
        end do

        dMolFraction(iFirstSpecies:iLastSpecies) = dMolFractionSave
        dPartialExcessGibbs(iFirstSpecies:iLastSpecies) = dPartialExcessSave

        write(*,'(A,1X,A,1X,I0,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6)') &
            'RKMP_RESPONSE_DEBUG', TRIM(cSolnPhaseName(m)), m, &
            'maxCurvEx=', dMaxCurvEx, 'maxCurvIdeal=', dMaxCurvIdeal, &
            'maxResponse=', dMaxResponse, 'fdMuErr=', dMaxFDResponseErr, &
            'maxCandidate=', MAXVAL(DABS(dCandidate)), &
            'maxCandidateDelta=', dMaxCandidateDelta

        deallocate(dX, dXBase, dGammaDir, dDX, dMuPlus, dMuMinus, dMuAnalytic, &
                   dMolFractionSave, dPartialExcessSave, dHloc, dHx, dC, &
                   dResponse, dCandidate, dIdealCandidate)
    end do

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Solve local response while preserving mole-fraction normalization.
    !>
    !> \details The final row and column introduce one normalization multiplier.
    !!          Requiring all infinitesimal mole-fraction changes to sum to zero
    !!          keeps their total equal to one. This is a coordinate constraint
    !!          within an active phase, not a constraint on the global assemblage.
    !---------------------------------------------------------------------------------------------------------
    subroutine SolveRKMPResponse(nSpecies, nElem, dHxLocal, dCLocal, dResp, INFO)

        integer, intent(in)                    :: nSpecies, nElem
        integer, intent(out)                   :: INFO
        real(8), intent(in), dimension(:,:)    :: dHxLocal, dCLocal
        real(8), intent(out), dimension(:,:)   :: dResp

        integer                                :: i, j, nEqn
        integer, dimension(:), allocatable     :: IPIV
        real(8), dimension(:,:), allocatable   :: dKKT, dRHS

        nEqn = nSpecies + 1
        allocate(dKKT(nEqn,nEqn), dRHS(nEqn,nElem), IPIV(nEqn))
        dKKT = 0D0
        dRHS = 0D0

        dKKT(1:nSpecies,1:nSpecies) = dHxLocal(1:nSpecies,1:nSpecies)
        do i = 1, nSpecies
            dKKT(i,nEqn) = -1D0
            dKKT(nEqn,i) = 1D0
        end do
        dRHS(1:nSpecies,1:nElem) = dCLocal(1:nSpecies,1:nElem)

        call DGESV(nEqn, nElem, dKKT, nEqn, IPIV, dRHS, nEqn, INFO)
        if (INFO == 0) then
            do j = 1, nElem
                dResp(1:nSpecies,j) = dRHS(1:nSpecies,j)
            end do
        else
            dResp = 0D0
        end if

        deallocate(dKKT, dRHS, IPIV)

    end subroutine SolveRKMPResponse

    !---------------------------------------------------------------------------------------------------------
    !> \brief Condense the ideal-only response for comparison with RKMP-corrected response.
    !---------------------------------------------------------------------------------------------------------
    subroutine SolveIdealResponse(nSpecies, nElem, dXLocal, dCLocal, dIdeal, INFO)

        integer, intent(in)                    :: nSpecies, nElem
        integer, intent(out)                   :: INFO
        real(8), intent(in), dimension(:)      :: dXLocal
        real(8), intent(in), dimension(:,:)    :: dCLocal
        real(8), intent(out), dimension(:,:)   :: dIdeal

        integer                                :: i
        real(8), dimension(:,:), allocatable   :: dHideal, dRespIdeal

        allocate(dHideal(nSpecies,nSpecies), dRespIdeal(nSpecies,nElem))
        dHideal = 0D0
        do i = 1, nSpecies
            dHideal(i,i) = 1D0 / dXLocal(i)
        end do

        call SolveRKMPResponse(nSpecies, nElem, dHideal, dCLocal, dRespIdeal, INFO)
        if (INFO == 0) then
            dIdeal = MATMUL(TRANSPOSE(dCLocal), dTotalMoles * dRespIdeal)
        else
            dIdeal = 0D0
        end if

        deallocate(dHideal, dRespIdeal)

    end subroutine SolveIdealResponse

end subroutine RKMPResponseDiagnostic
