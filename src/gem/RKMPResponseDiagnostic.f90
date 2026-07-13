!-------------------------------------------------------------------------------------------------------------
!> \file    RKMPResponseDiagnostic.f90
!> \brief   Diagnostic-only local constrained RKMP composition response.
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

    integer                              :: i, j, k, m, p, nLocalSpecies, iFirstSpecies, iLastSpecies
    integer                              :: INFO
    real(8)                              :: dTotalMoles, dEps, dNorm, dMaxCurvEx, dMaxCurvIdeal
    real(8)                              :: dMaxResponse, dMaxFDResponseErr, dMaxCandidateDelta
    real(8), dimension(:), allocatable   :: dX, dXBase, dGammaDir, dDX, dMuPlus, dMuMinus, dMuAnalytic
    real(8), dimension(:), allocatable   :: dMolFractionSave, dPartialExcessSave
    real(8), dimension(:,:), allocatable :: dHloc, dHx, dC, dResponse, dCandidate, dIdealCandidate

    if (.NOT. lDebugRKMPHessianFD) return

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

        ! Convert mole-number curvature to fixed-phase-amount mole-fraction curvature:
        ! dmu_ex = Hn dn = (N * Hn) dx when sum(dx)=0 and dN=0.
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

        ! Finite-difference check: validate dmu_ex/dx from the unconstrained RKMP Hessian conversion.
        dGammaDir = 0D0
        do j = 1, nElements
            dGammaDir(j) = 1D0 / DFLOAT(j)
        end do
        dDX = MATMUL(dResponse, dGammaDir)
        dNorm = MAXVAL(DABS(dDX))
        if (dNorm > 0D0) dDX = dDX / dNorm

        dEps = 1D-6
        do i = 1, nLocalSpecies
            if (dDX(i) < 0D0) dEps = DMIN1(dEps, 0.25D0 * dX(i) / DABS(dDX(i)))
        end do
        dEps = DMAX1(dEps, 1D-10)

        dMolFractionSave = dMolFraction(iFirstSpecies:iLastSpecies)
        dPartialExcessSave = dPartialExcessGibbs(iFirstSpecies:iLastSpecies)

        dMolFraction(iFirstSpecies:iLastSpecies) = dX + dEps * dDX
        dPartialExcessGibbs(iFirstSpecies:iLastSpecies) = 0D0
        call CompExcessGibbsEnergyRKMP(m)
        dMuPlus = dPartialExcessGibbs(iFirstSpecies:iLastSpecies)

        dMolFraction(iFirstSpecies:iLastSpecies) = dX - dEps * dDX
        dPartialExcessGibbs(iFirstSpecies:iLastSpecies) = 0D0
        call CompExcessGibbsEnergyRKMP(m)
        dMuMinus = dPartialExcessGibbs(iFirstSpecies:iLastSpecies)

        dMolFraction(iFirstSpecies:iLastSpecies) = dMolFractionSave
        dPartialExcessGibbs(iFirstSpecies:iLastSpecies) = dPartialExcessSave

        dMuAnalytic = MATMUL(dTotalMoles * dHloc, dDX)
        dMaxFDResponseErr = MAXVAL(DABS(((dMuPlus - dMuMinus) / (2D0 * dEps)) - dMuAnalytic))

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
