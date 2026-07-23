!-------------------------------------------------------------------------------------------------------------
!> \file    RKMPDiagnosticUtilities.f90
!> \brief   Independent scalar RKMP excess-energy evaluator for diagnostics.
!>
!> \details This file evaluates the supported extensive binary RKMP excess energy
!!          directly from a supplied local mole vector. It is deliberately
!!          independent of the analytic Hessian expression, allowing finite
!!          differences to test CompExcessGibbsEnergyRKMP_unconstrained.
!!
!!          For each binary parameter, total phase moles N convert supplied
!!          species moles into mole fractions x1 and x2. The scalar contribution
!!          is N times the database coefficient, the binary mixing factor x1*x2,
!!          and the requested power of the composition contrast x1-x2.
!!
!!          It reads production parameter definitions but does not change
!!          Thermochimica state, GEMNewton matrices, or the phase assemblage.
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

    ! Reconstruct the scalar excess energy whose second derivatives with respect
    ! to the supplied species moles define the local Hessian Hloc.
    ! Keeping this expression outside the Hessian routine reduces common-mode
    ! verification risk.
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
