!> \file    RKMPGEMIdealReconstructionDiagnostic.f90
!> \brief   Diagnostic-only reconstruction of GEMNewton's ideal RKMP phase contribution.

subroutine RKMPGEMIdealReconstructionDiagnostic(A, B, nVar)

    USE ModuleThermo
    USE ModuleGEMSolver, ONLY: lDebugRKMPHessianFD

    implicit none

    integer, intent(in)                    :: nVar
    real(8), intent(in), dimension(nVar,nVar) :: A
    real(8), intent(in), dimension(nVar)    :: B

    integer                                :: i, j, k, l, p, iFirst, iLast, nPhaseSpecies
    integer                                :: iSolnPhases, iAssemblageSlot, iPhaseVar
    real(8)                                :: dPhaseMoles, dCbar, dBpRecon
    real(8)                                :: dMaxAeeDiff, dMaxAepDiff, dMaxBeDiff, dBpDiff
    real(8)                                :: dMaxActualAepDiff, dActualBpDiff
    real(8)                                :: dMaxAee, dMaxAep, dMaxBe, dMaxTotalAee
    real(8), allocatable, dimension(:)     :: dX, dLocalMoles, dMu
    real(8), allocatable, dimension(:)     :: dDirectAep, dReconAep, dDirectBe, dReconBe
    real(8), allocatable, dimension(:,:)   :: dC, dDirectAee, dReconAee

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
                 dDirectBe(nElements), dReconBe(nElements))

        dLocalMoles = dMolesSpecies(iFirst:iLast)
        dPhaseMoles = SUM(dLocalMoles)
        if (dPhaseMoles <= 1D-30) then
            deallocate(dX, dLocalMoles, dMu, dC, dDirectAee, dReconAee, dDirectAep, &
                       dReconAep, dDirectBe, dReconBe)
            cycle LOOP_SOLN
        end if

        dX  = dLocalMoles / dPhaseMoles
        dMu = dChemicalPotential(iFirst:iLast)

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

        write(*,'(A,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6)') &
            'RKMP_GEM_IDEAL_RECON', TRIM(cSolnPhaseName(k)), k, 'phaseVar=', iPhaseVar, &
            'N=', dPhaseMoles, 'maxAee=', dMaxAee, 'maxAep=', dMaxAep, 'maxBe=', dMaxBe, &
            'maxTotalAee=', dMaxTotalAee, 'reconAeeDiff=', dMaxAeeDiff, &
            'reconAepDiff=', dMaxAepDiff, 'reconBeDiff=', dMaxBeDiff, &
            'reconBpDiff=', dBpDiff, 'actualAepDiff=', dMaxActualAepDiff, &
            'actualBpDiff=', dActualBpDiff

        deallocate(dX, dLocalMoles, dMu, dC, dDirectAee, dReconAee, dDirectAep, &
                   dReconAep, dDirectBe, dReconBe)

    end do LOOP_SOLN

end subroutine RKMPGEMIdealReconstructionDiagnostic
