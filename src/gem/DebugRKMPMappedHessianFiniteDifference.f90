!-------------------------------------------------------------------------------------------------------------
!> \file    DebugRKMPMappedHessianFiniteDifference.f90
!> \brief   Validate mapped RKMP composition curvature by finite differences.
!-------------------------------------------------------------------------------------------------------------

subroutine DebugRKMPMappedHessianFiniteDifference(iSolnIndex, dJ, dKmap, nLocalSpecies)

    USE ModuleThermo

    implicit none

    integer, intent(in)                   :: iSolnIndex, nLocalSpecies
    real(8), dimension(nLocalSpecies,nElements), intent(in) :: dJ
    real(8), dimension(nElements,nElements), intent(in) :: dKmap

    integer                               :: i, j, iFirstSpecies
    real(8)                               :: dAnalytic, dFD, dE0, dEp, dEm, dEps, dSumV, dNormV
    real(8), dimension(:), allocatable    :: dU, dV, dMoles0, dWork

    if (nElements <= 0) return

    allocate(dU(nElements), dV(nLocalSpecies), dMoles0(nLocalSpecies), dWork(nLocalSpecies))

    iFirstSpecies = nSpeciesPhase(iSolnIndex-1) + 1
    dMoles0       = dMolesSpecies(iFirstSpecies:nSpeciesPhase(iSolnIndex))

    dU = 0D0
    do i = 1, nElements
        dU(i) = 1D0 / DFLOAT(i)
    end do

    dV = 0D0
    do j = 1, nElements
        do i = 1, nLocalSpecies
            dV(i) = dV(i) + dJ(i,j) * dU(j)
        end do
    end do

    dAnalytic = 0D0
    do j = 1, nElements
        do i = 1, nElements
            dAnalytic = dAnalytic + dU(i) * dKmap(i,j) * dU(j)
        end do
    end do

    dNormV = MAXVAL(ABS(dV))
    if (dNormV <= 0D0) then
        deallocate(dU, dV, dMoles0, dWork)
        return
    end if

    dEps = 1D-4 * MAX(1D0, MAXVAL(dMoles0) / dNormV)
    do i = 1, nLocalSpecies
        if (dV(i) < 0D0) dEps = MIN(dEps, 0.25D0 * dMoles0(i) / ABS(dV(i)))
    end do
    if (dEps <= 0D0) then
        deallocate(dU, dV, dMoles0, dWork)
        return
    end if

    call CompRKMPBinaryExcessGibbsFromMoles(iSolnIndex, nLocalSpecies, dMoles0, dE0)
    dWork = dMoles0 + dEps * dV
    call CompRKMPBinaryExcessGibbsFromMoles(iSolnIndex, nLocalSpecies, dWork, dEp)
    dWork = dMoles0 - dEps * dV
    call CompRKMPBinaryExcessGibbsFromMoles(iSolnIndex, nLocalSpecies, dWork, dEm)

    dFD   = (dEp - 2D0*dE0 + dEm) / (dEps*dEps)
    dSumV = SUM(dV)

    write(*,'(A,1X,A,1X,I0,1X,A,1X,ES15.6E3,1X,A,1X,ES15.6E3,1X,A,1X,ES15.6E3,1X,A,1X,ES15.6E3)') &
        'RKMP_MAPPED_FD_DEBUG', TRIM(cSolnPhaseName(iSolnIndex)), iSolnIndex, &
        'analytic=', dAnalytic, 'fd=', dFD, 'sumv=', dSumV, 'eps=', dEps

    deallocate(dU, dV, dMoles0, dWork)

    return

end subroutine
