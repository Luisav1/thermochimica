!> \file    RKMPGEMIdealReconstructionDiagnostic.f90
!> \brief   Diagnostic-only reconstruction of GEMNewton's ideal RKMP phase contribution.

subroutine RKMPGEMIdealReconstructionDiagnostic(A, B, nVar)

    USE ModuleThermo
    USE ModuleGEMSolver, ONLY: lDebugRKMPHessianFD

    implicit none

    interface
        subroutine CompExcessGibbsEnergyRKMP_unconstrained(iSolnIndex,dHess)
            integer, intent(in)                  :: iSolnIndex
            real(8), intent(out), dimension(:,:) :: dHess
        end subroutine CompExcessGibbsEnergyRKMP_unconstrained
    end interface

    integer, intent(in)                    :: nVar
    real(8), intent(in), dimension(nVar,nVar) :: A
    real(8), intent(in), dimension(nVar)    :: B

    integer                                :: i, j, k, l, p, iFirst, iLast, nPhaseSpecies
    integer                                :: iSolnPhases, iAssemblageSlot, iPhaseVar, INFO
    real(8)                                :: dPhaseMoles, dCbar, dBpRecon
    real(8)                                :: dMaxAeeDiff, dMaxAepDiff, dMaxBeDiff, dBpDiff
    real(8)                                :: dMaxActualAepDiff, dActualBpDiff
    real(8)                                :: dMaxAee, dMaxAep, dMaxBe, dMaxTotalAee, dMaxIdealResponseAeeDiff
    real(8)                                :: dMaxDeltaAee, dMaxDeltaBe, dMaxNewAee, dDeltaRatio
    real(8), allocatable, dimension(:)     :: dX, dLocalMoles, dMu
    real(8), allocatable, dimension(:)     :: dDirectAep, dReconAep, dDirectBe, dReconBe, dNewBe, dDeltaBe
    real(8), allocatable, dimension(:,:)   :: dC, dDirectAee, dReconAee
    real(8), allocatable, dimension(:,:)   :: dHloc, dHx, dResponse, dIdealResponse
    real(8), allocatable, dimension(:,:)   :: dMuRHS, dMuResponse, dIdealMuResponse
    real(8), allocatable, dimension(:,:)   :: dCenteredAee, dIdealCenteredAee, dOuterAee, dNewAee, dDeltaAee

    if (.NOT. lDebugRKMPHessianFD) return
    if (nSolnPhases <= 0) return

    LOOP_SOLN: do iSolnPhases = 1, nSolnPhases

        iAssemblageSlot = nElements - iSolnPhases + 1
        k = -iAssemblage(iAssemblageSlot)
        if (k <= 0) cycle LOOP_SOLN
        if (cSolnPhaseType(k) /= 'RKMP') cycle LOOP_SOLN

        iPhaseVar = nElements + iSolnPhases
        if (iPhaseVar > nVar) cycle LOOP_SOLN

        iFirst = nSpeciesPhase(k-1) + 1
        iLast  = nSpeciesPhase(k)
        nPhaseSpecies = iLast - iFirst + 1
        if (nPhaseSpecies <= 0) cycle LOOP_SOLN

        allocate(dX(nPhaseSpecies), dLocalMoles(nPhaseSpecies), dMu(nPhaseSpecies), &
                 dC(nPhaseSpecies,nElements), dDirectAee(nElements,nElements), &
                 dReconAee(nElements,nElements), dDirectAep(nElements), dReconAep(nElements), &
                 dDirectBe(nElements), dReconBe(nElements), dNewBe(nElements), dDeltaBe(nElements), &
                 dHloc(nPhaseSpecies,nPhaseSpecies), dHx(nPhaseSpecies,nPhaseSpecies), &
                 dResponse(nPhaseSpecies,nElements), dIdealResponse(nPhaseSpecies,nElements), &
                 dMuRHS(nPhaseSpecies,1), dMuResponse(nPhaseSpecies,1), &
                 dIdealMuResponse(nPhaseSpecies,1), &
                 dCenteredAee(nElements,nElements), dIdealCenteredAee(nElements,nElements), &
                 dOuterAee(nElements,nElements), dNewAee(nElements,nElements), dDeltaAee(nElements,nElements))

        dLocalMoles = dMolesSpecies(iFirst:iLast)
        dPhaseMoles = SUM(dLocalMoles)
        if (dPhaseMoles <= 1D-30) then
            deallocate(dX, dLocalMoles, dMu, dC, dDirectAee, dReconAee, dDirectAep, &
                       dReconAep, dDirectBe, dReconBe, dNewBe, dDeltaBe, dHloc, dHx, &
                       dResponse, dIdealResponse, dCenteredAee, dIdealCenteredAee, &
                       dMuRHS, dMuResponse, dIdealMuResponse, dOuterAee, dNewAee, dDeltaAee)
            cycle LOOP_SOLN
        end if

        dX  = dLocalMoles / dPhaseMoles
        dMu = dChemicalPotential(iFirst:iLast)
        dMuRHS(:,1) = dMu

        do j = 1, nElements
            do i = 1, nPhaseSpecies
                p = iFirst + i - 1
                dC(i,j) = dStoichSpecies(p,j) / DFLOAT(iParticlesPerMole(p))
            end do
        end do

        dDirectAee = 0D0
        dReconAee  = 0D0
        dDirectAep = 0D0
        dReconAep  = 0D0
        dDirectBe  = 0D0
        dReconBe   = 0D0
        dBpRecon   = 0D0
        dOuterAee  = 0D0

        ! Direct GEMNewton phase contribution, written in the same species-mole form as GEMNewton.
        do l = 1, nPhaseSpecies
            do j = 1, nElements
                dDirectAep(j) = dDirectAep(j) + dLocalMoles(l) * dC(l,j)
                dDirectBe(j)  = dDirectBe(j)  + dLocalMoles(l) * dC(l,j) * (dMu(l) - 1D0)
                do i = 1, nElements
                    dDirectAee(i,j) = dDirectAee(i,j) + dLocalMoles(l) * dC(l,i) * dC(l,j)
                end do
            end do
            dBpRecon = dBpRecon + dLocalMoles(l) * dMu(l)
        end do

        ! Same contribution reconstructed from local ideal-condensed quantities: N, x, C, and mu.
        do j = 1, nElements
            dCbar = 0D0
            do l = 1, nPhaseSpecies
                dCbar = dCbar + dX(l) * dC(l,j)
            end do
            dReconAep(j) = dPhaseMoles * dCbar
            do l = 1, nPhaseSpecies
                dReconBe(j) = dReconBe(j) + dPhaseMoles * dX(l) * dC(l,j) * (dMu(l) - 1D0)
            end do

            do i = 1, nElements
                do l = 1, nPhaseSpecies
                    dReconAee(i,j) = dReconAee(i,j) + dPhaseMoles * dX(l) * dC(l,i) * dC(l,j)
                end do
            end do
        end do

        do j = 1, nElements
            do i = 1, nElements
                dOuterAee(i,j) = dReconAep(i) * dReconAep(j) / dPhaseMoles
            end do
        end do

        call CompExcessGibbsEnergyRKMP_unconstrained(k, dHloc)
        dHx = dPhaseMoles * dHloc
        do i = 1, nPhaseSpecies
            dHx(i,i) = dHx(i,i) + 1D0 / dX(i)
        end do

        call SolveLocalResponse(nPhaseSpecies, nElements, dHx, dC, dResponse, INFO)
        if (INFO /= 0) then
            write(*,'(A,1X,A,1X,I0,1X,A,1X,I0)') &
                'RKMP_GEM_RESPONSE_CONDENSE', TRIM(cSolnPhaseName(k)), k, 'response_solve_info=', INFO
            deallocate(dX, dLocalMoles, dMu, dC, dDirectAee, dReconAee, dDirectAep, &
                       dReconAep, dDirectBe, dReconBe, dNewBe, dDeltaBe, dHloc, dHx, &
                       dResponse, dIdealResponse, dCenteredAee, dIdealCenteredAee, &
                       dMuRHS, dMuResponse, dIdealMuResponse, dOuterAee, dNewAee, dDeltaAee)
            cycle LOOP_SOLN
        end if

        call SolveLocalResponse(nPhaseSpecies, 1, dHx, dMuRHS, dMuResponse, INFO)
        if (INFO /= 0) then
            deallocate(dX, dLocalMoles, dMu, dC, dDirectAee, dReconAee, dDirectAep, &
                       dReconAep, dDirectBe, dReconBe, dNewBe, dDeltaBe, dHloc, dHx, &
                       dResponse, dIdealResponse, dCenteredAee, dIdealCenteredAee, &
                       dMuRHS, dMuResponse, dIdealMuResponse, dOuterAee, dNewAee, dDeltaAee)
            cycle LOOP_SOLN
        end if

        call SolveIdealResponse(nPhaseSpecies, nElements, dX, dC, dIdealResponse, INFO)
        if (INFO /= 0) then
            write(*,'(A,1X,A,1X,I0,1X,A,1X,I0)') &
                'RKMP_GEM_RESPONSE_CONDENSE', TRIM(cSolnPhaseName(k)), k, 'ideal_solve_info=', INFO
            deallocate(dX, dLocalMoles, dMu, dC, dDirectAee, dReconAee, dDirectAep, &
                       dReconAep, dDirectBe, dReconBe, dNewBe, dDeltaBe, dHloc, dHx, &
                       dResponse, dIdealResponse, dCenteredAee, dIdealCenteredAee, &
                       dMuRHS, dMuResponse, dIdealMuResponse, dOuterAee, dNewAee, dDeltaAee)
            cycle LOOP_SOLN
        end if

        call SolveIdealResponse(nPhaseSpecies, 1, dX, dMuRHS, dIdealMuResponse, INFO)
        if (INFO /= 0) then
            deallocate(dX, dLocalMoles, dMu, dC, dDirectAee, dReconAee, dDirectAep, &
                       dReconAep, dDirectBe, dReconBe, dNewBe, dDeltaBe, dHloc, dHx, &
                       dResponse, dIdealResponse, dCenteredAee, dIdealCenteredAee, &
                       dMuRHS, dMuResponse, dIdealMuResponse, dOuterAee, dNewAee, dDeltaAee)
            cycle LOOP_SOLN
        end if

        dCenteredAee      = MATMUL(TRANSPOSE(dC), dPhaseMoles * dResponse)
        dIdealCenteredAee = MATMUL(TRANSPOSE(dC), dPhaseMoles * dIdealResponse)
        dNewAee           = dOuterAee + dCenteredAee
        dDeltaAee         = dNewAee - dReconAee
        ! Re-condense the current local stationarity forcing through each response.  Using deltaA*gamma here
        ! would agree only after the species chemical-potential residual is already zero.
        dDeltaBe          = MATMUL(TRANSPOSE(dC), dPhaseMoles * &
            (dMuResponse(:,1) - dIdealMuResponse(:,1)))
        dNewBe            = dReconBe + dDeltaBe

        dMaxAeeDiff = MAXVAL(DABS(dDirectAee - dReconAee))
        dMaxAepDiff = MAXVAL(DABS(dDirectAep - dReconAep))
        dMaxBeDiff  = MAXVAL(DABS(dDirectBe  - dReconBe))
        dBpDiff     = DABS(dGibbsSolnPhase(k) - dBpRecon)

        dMaxActualAepDiff = MAXVAL(DABS(A(1:nElements,iPhaseVar) - dReconAep))
        dActualBpDiff     = DABS(B(iPhaseVar) - dBpRecon)

        dMaxAee      = MAXVAL(DABS(dReconAee))
        dMaxAep      = MAXVAL(DABS(dReconAep))
        dMaxBe       = MAXVAL(DABS(dReconBe))
        dMaxTotalAee = MAXVAL(DABS(A(1:nElements,1:nElements)))
        dMaxIdealResponseAeeDiff = MAXVAL(DABS((dOuterAee + dIdealCenteredAee) - dReconAee))
        dMaxDeltaAee = MAXVAL(DABS(dDeltaAee))
        dMaxDeltaBe  = MAXVAL(DABS(dDeltaBe))
        dMaxNewAee   = MAXVAL(DABS(dNewAee))
        dDeltaRatio  = dMaxDeltaAee / DMAX1(dMaxAee, 1D-30)

        write(*,'(A,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6)') &
            'RKMP_GEM_IDEAL_RECON', TRIM(cSolnPhaseName(k)), k, 'phaseVar=', iPhaseVar, &
            'N=', dPhaseMoles, 'maxAee=', dMaxAee, 'maxAep=', dMaxAep, 'maxBe=', dMaxBe, &
            'maxTotalAee=', dMaxTotalAee, 'reconAeeDiff=', dMaxAeeDiff, &
            'reconAepDiff=', dMaxAepDiff, 'reconBeDiff=', dMaxBeDiff, &
            'reconBpDiff=', dBpDiff, 'actualAepDiff=', dMaxActualAepDiff, &
            'actualBpDiff=', dActualBpDiff

        write(*,'(A,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6)') &
            'RKMP_GEM_RESPONSE_CONDENSE', TRIM(cSolnPhaseName(k)), k, 'phaseVar=', iPhaseVar, &
            'idealResponseAeeDiff=', dMaxIdealResponseAeeDiff, &
            'maxNewAee=', dMaxNewAee, 'maxDeltaAee=', dMaxDeltaAee, &
            'deltaAeeRatio=', dDeltaRatio, 'maxDeltaBe=', dMaxDeltaBe, &
            'maxNewBe=', MAXVAL(DABS(dNewBe))

        deallocate(dX, dLocalMoles, dMu, dC, dDirectAee, dReconAee, dDirectAep, &
                   dReconAep, dDirectBe, dReconBe, dNewBe, dDeltaBe, dHloc, dHx, &
                   dResponse, dIdealResponse, dCenteredAee, dIdealCenteredAee, &
                   dMuRHS, dMuResponse, dIdealMuResponse, dOuterAee, dNewAee, dDeltaAee)

    end do LOOP_SOLN

contains

    subroutine SolveLocalResponse(nSpecies, nElem, dHxLocal, dCLocal, dResp, INFO)

        integer, intent(in)                    :: nSpecies, nElem
        integer, intent(out)                   :: INFO
        real(8), intent(in), dimension(:,:)    :: dHxLocal, dCLocal
        real(8), intent(out), dimension(:,:)   :: dResp

        integer                                :: ii, jj, nEqn
        integer, dimension(:), allocatable     :: IPIV
        real(8), dimension(:,:), allocatable   :: dKKT, dRHS

        nEqn = nSpecies + 1
        allocate(dKKT(nEqn,nEqn), dRHS(nEqn,nElem), IPIV(nEqn))
        dKKT = 0D0
        dRHS = 0D0

        dKKT(1:nSpecies,1:nSpecies) = dHxLocal(1:nSpecies,1:nSpecies)
        do ii = 1, nSpecies
            dKKT(ii,nEqn) = -1D0
            dKKT(nEqn,ii) = 1D0
        end do
        dRHS(1:nSpecies,1:nElem) = dCLocal(1:nSpecies,1:nElem)

        call DGESV(nEqn, nElem, dKKT, nEqn, IPIV, dRHS, nEqn, INFO)
        if (INFO == 0) then
            do jj = 1, nElem
                dResp(1:nSpecies,jj) = dRHS(1:nSpecies,jj)
            end do
        else
            dResp = 0D0
        end if

        deallocate(dKKT, dRHS, IPIV)

    end subroutine SolveLocalResponse

    subroutine SolveIdealResponse(nSpecies, nElem, dXLocal, dCLocal, dResp, INFO)

        integer, intent(in)                    :: nSpecies, nElem
        integer, intent(out)                   :: INFO
        real(8), intent(in), dimension(:)      :: dXLocal
        real(8), intent(in), dimension(:,:)    :: dCLocal
        real(8), intent(out), dimension(:,:)   :: dResp

        integer                                :: ii
        real(8), dimension(:,:), allocatable   :: dHideal

        allocate(dHideal(nSpecies,nSpecies))
        dHideal = 0D0
        do ii = 1, nSpecies
            dHideal(ii,ii) = 1D0 / dXLocal(ii)
        end do

        call SolveLocalResponse(nSpecies, nElem, dHideal, dCLocal, dResp, INFO)

        deallocate(dHideal)

    end subroutine SolveIdealResponse

end subroutine RKMPGEMIdealReconstructionDiagnostic
