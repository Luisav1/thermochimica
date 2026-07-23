!-------------------------------------------------------------------------------------------------------------
!> \file    RKMPMappedHessianDiagnostic.f90
!> \brief   Diagnostic-only projection of RKMP species curvature to GEM element directions.
!>
!> \details This historical Stage 1A/1B diagnostic constructs a matrix J whose
!!          column for element e describes one local species-mole direction
!!          induced by that element. The phase-average stoichiometry is removed,
!!          so each direction redistributes species without changing total phase
!!          moles.
!!
!!          The projected matrix Kmap=transpose(J)*Hloc*J measures RKMP curvature
!!          between pairs of those element-driven directions. Hloc is the local
!!          species-mole Hessian; multiplying on the right maps an element
!!          direction into species space, and multiplying by transpose(J) maps
!!          the resulting response back to element-direction coordinates.
!!          This established that the local Hessian projection was mathematically
!!          and numerically consistent.
!!
!!          Kmap is diagnostic curvature, not the production GEM correction.
!!          Directly adding it to GEMNewton was rejected because GEM requires the
!!          constrained phase response and matching right-hand-side condensation.
!-------------------------------------------------------------------------------------------------------------

subroutine RKMPMappedHessianDiagnostic

    USE ModuleThermo
    USE ModuleGEMSolver, ONLY: lDebugRKMPHessianFD

    implicit none

    integer                                :: k, m, i, j, p, q, nLocalSpecies, iFirstSpecies
    real(8)                                :: dTotalMoles, dCbar
    real(8), dimension(:), allocatable     :: dLocalMoles
    real(8), dimension(:,:), allocatable   :: dJ, dHloc, dKmap

    interface
        subroutine CompExcessGibbsEnergyRKMP_unconstrained(iSolnIndex,dHess)
            integer, intent(in)                  :: iSolnIndex
            real(8), intent(out), dimension(:,:) :: dHess
        end subroutine CompExcessGibbsEnergyRKMP_unconstrained
    end interface

    do k = 1, nSolnPhases
        m = -iAssemblage(nElements - k + 1)
        if (cSolnPhaseType(m) /= 'RKMP') cycle

        nLocalSpecies = nSpeciesPhase(m) - nSpeciesPhase(m-1)
        if (nLocalSpecies <= 0) cycle
        iFirstSpecies = nSpeciesPhase(m-1) + 1

        allocate(dLocalMoles(nLocalSpecies), dJ(nLocalSpecies,nElements), &
                 dHloc(nLocalSpecies,nLocalSpecies), dKmap(nElements,nElements))

        dLocalMoles = dMolesSpecies(iFirstSpecies:nSpeciesPhase(m))
        dTotalMoles = SUM(dLocalMoles)
        if (dTotalMoles <= 0D0) then
            deallocate(dLocalMoles, dJ, dHloc, dKmap)
            cycle
        end if

        ! Stage 1A: subtract the phase-average stoichiometry so every column of J
        ! changes composition while preserving the total phase amount.
        do j = 1, nElements
            dCbar = 0D0
            do i = 1, nLocalSpecies
                p = iFirstSpecies + i - 1
                dCbar = dCbar + (dLocalMoles(i) / dTotalMoles) * &
                    dStoichSpecies(p,j) / DFLOAT(iParticlesPerMole(p))
            end do
            do i = 1, nLocalSpecies
                p = iFirstSpecies + i - 1
                dJ(i,j) = dLocalMoles(i) * &
                    (dStoichSpecies(p,j) / DFLOAT(iParticlesPerMole(p)) - dCbar)
            end do
        end do

        call CompExcessGibbsEnergyRKMP_unconstrained(m, dHloc)

        ! Stage 1B: project local species curvature between every pair of those
        ! element-driven composition directions.
        dKmap = 0D0
        do j = 1, nElements
            do i = 1, nElements
                do q = 1, nLocalSpecies
                    do p = 1, nLocalSpecies
                        dKmap(i,j) = dKmap(i,j) + dJ(p,i) * dHloc(p,q) * dJ(q,j)
                    end do
                end do
            end do
        end do

        if (lDebugRKMPHessianFD) call DebugRKMPMappedHessianFiniteDifference(m, dJ, dKmap, nLocalSpecies)

        deallocate(dLocalMoles, dJ, dHloc, dKmap)
    end do

    return

end subroutine RKMPMappedHessianDiagnostic
