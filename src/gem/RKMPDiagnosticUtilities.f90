!-------------------------------------------------------------------------------------------------------------
!> \file    RKMPDiagnosticUtilities.f90
!> \brief   Diagnostic-only RKMP excess-energy curvature utilities.
!>
!> These routines intentionally do not modify the GEM Newton matrix.  They provide
!> a local binary-RKMP excess Gibbs energy evaluator, a finite-difference local
!> species-space Hessian, and a mapped finite-difference check for projected
!> composition directions.
!-------------------------------------------------------------------------------------------------------------

subroutine CompRKMPBinaryExcessGibbsFromMoles(iSolnIndex, nLocalSpecies, dLocalMoles, dGex)

    USE ModuleThermo

    implicit none

    integer, intent(in)                  :: iSolnIndex
    integer, intent(in)                  :: nLocalSpecies
    real(8), dimension(nLocalSpecies), intent(in) :: dLocalMoles
    real(8), intent(out)                 :: dGex

    integer                              :: iParam, i1, i2, iExponent
    real(8)                              :: dTotalMoles, x1, x2, dx

    dGex         = 0D0
        if ((cSolnPhaseType(iSolnIndex) /= 'RKMP') .AND. (cSolnPhaseType(iSolnIndex) /= 'RKMPM')) return
    if (nParamPhase(iSolnIndex) - nParamPhase(iSolnIndex-1) == 0) return

    dTotalMoles = SUM(dLocalMoles(1:nLocalSpecies))
    if (dTotalMoles <= 0D0) return

    do iParam = nParamPhase(iSolnIndex-1)+1, nParamPhase(iSolnIndex)
        ! Stage-1 diagnostics are restricted to binary RKMP excess terms.  Higher-order
        ! Muggiano and RKMPM magnetic curvature terms are intentionally skipped here.
        if (iRegularParam(iParam,1) /= 2) cycle

        i1 = iRegularParam(iParam,2)
        i2 = iRegularParam(iParam,3)
        if ((i1 < 1) .OR. (i1 > nLocalSpecies)) cycle
        if ((i2 < 1) .OR. (i2 > nLocalSpecies)) cycle

        iExponent = iRegularParam(iParam,4)
        x1        = dLocalMoles(i1) / dTotalMoles
        x2        = dLocalMoles(i2) / dTotalMoles
        dx        = x1 - x2

        dGex = dGex + dTotalMoles * dExcessGibbsParam(iParam) * x1 * x2 * dx**iExponent
    end do

    return

end subroutine CompRKMPBinaryExcessGibbsFromMoles