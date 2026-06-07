
    !---------------------------------------------------------------------------------------------------------
    ! Temporary RKMP finite-difference diagnostic.
    !
    ! Select the first supported binary RKMP interaction in this phase and define a composition-changing
    ! direction v with +1 for species ia and -1 for species ib.  Compare the analytic directional curvature
    !
    !     v^T Hloc v
    !
    ! against a central finite difference of the supported binary RKMP excess-energy expression:
    !
    !     (Gex(n + eps*v) - 2*Gex(n) + Gex(n - eps*v)) / eps^2
    !
    ! This helper is diagnostic-only.  It does not modify thermodynamic state or solver matrices.
    !---------------------------------------------------------------------------------------------------------

subroutine DebugRKMPHessianFiniteDifference(iSolnIndex,dHess)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer, intent(in)                  :: iSolnIndex
    real(8), intent(in), dimension(:,:)  :: dHess

    integer :: i, j, iParam, iStep, ia, ib, iaParam, ibParam, iExponentParam
    integer :: iFirstSpecies, iLastSpecies, nPhaseSpecies
    real(8) :: dAnalytic, dAnalyticParam, dEnergyMinus, dEnergyZero, dEnergyPlus
    real(8) :: dEps, dFD, dFDParam, dN, dRelativeError, dScale
    real(8), allocatable, dimension(:) :: dDirection, dMolesLocal, dMolesMinus, dMolesPlus

    iFirstSpecies = nSpeciesPhase(iSolnIndex-1) + 1
    iLastSpecies  = nSpeciesPhase(iSolnIndex)
    nPhaseSpecies = iLastSpecies - iFirstSpecies + 1
    if (nPhaseSpecies <= 0) return

    ia = 0
    ib = 0
    LOOP_FIND_PARAM: do iParam = nParamPhase(iSolnIndex-1)+1, nParamPhase(iSolnIndex)
        if (iRegularParam(iParam,1) /= 2) cycle LOOP_FIND_PARAM
        if (iRegularParam(iParam,4) < 0) cycle LOOP_FIND_PARAM
        ia = iRegularParam(iParam,2)
        ib = iRegularParam(iParam,3)
        if ((ia < 1) .OR. (ia > nPhaseSpecies) .OR. (ib < 1) .OR. (ib > nPhaseSpecies)) then
            ia = 0
            ib = 0
            cycle LOOP_FIND_PARAM
        end if
        exit LOOP_FIND_PARAM
    end do LOOP_FIND_PARAM
    if ((ia == 0) .OR. (ib == 0)) return

    allocate(dDirection(nPhaseSpecies), dMolesLocal(nPhaseSpecies))
    allocate(dMolesMinus(nPhaseSpecies), dMolesPlus(nPhaseSpecies))

    dDirection = 0D0
    dDirection(ia) = 1D0
    dDirection(ib) = -1D0

    dN = 0D0
    do i = 1, nPhaseSpecies
        dMolesLocal(i) = dMolesSpecies(iFirstSpecies + i - 1)
        dN = dN + dMolesLocal(i)
    end do

    dScale = DMIN1(dN,0.25D0*DMIN1(dMolesLocal(ia),dMolesLocal(ib)))
    if (dScale <= 0D0) then
        deallocate(dDirection,dMolesLocal,dMolesMinus,dMolesPlus)
        return
    end if

    dAnalytic = 0D0
    do i = 1, nPhaseSpecies
        do j = 1, nPhaseSpecies
            dAnalytic = dAnalytic + dDirection(i) * dHess(i,j) * dDirection(j)
        end do
    end do

    dEnergyZero = CompSupportedBinaryRKMPExcessEnergy(dMolesLocal)

    do iStep = 1, 5
        dEps = dScale * 10D0**(-iStep)
        dMolesMinus = dMolesLocal - dEps*dDirection
        dMolesPlus  = dMolesLocal + dEps*dDirection
        dEnergyMinus = CompSupportedBinaryRKMPExcessEnergy(dMolesMinus)
        dEnergyPlus  = CompSupportedBinaryRKMPExcessEnergy(dMolesPlus)
        dFD = (dEnergyPlus - 2D0*dEnergyZero + dEnergyMinus) / (dEps*dEps)
        dRelativeError = DABS(dFD-dAnalytic) / DMAX1(DABS(dAnalytic),1D-30)

        write(*,'(A,1X,F12.4,1X,A,1X,I5,1X,A,1X,I5,1X,A,1X,I5,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6)') &
            'RKMP_FD_DEBUG T=', dTemperature, 'k=', iSolnIndex, 'ia=', ia, 'ib=', ib, 'eps=', dEps, &
            'analytic=', dAnalytic, 'fd=', dFD, 'relerr=', dRelativeError
    end do

    ! Per-parameter diagnostic at the middle epsilon scale used above.  This identifies which binary
    ! interaction term contributes most to any aggregate analytic/finite-difference mismatch.
    dEps = dScale * 1D-3
    dMolesMinus = dMolesLocal - dEps*dDirection
    dMolesPlus  = dMolesLocal + dEps*dDirection

    LOOP_PARAM_DIAG: do iParam = nParamPhase(iSolnIndex-1)+1, nParamPhase(iSolnIndex)
        if (iRegularParam(iParam,1) /= 2) cycle LOOP_PARAM_DIAG

        iaParam = iRegularParam(iParam,2)
        ibParam = iRegularParam(iParam,3)
        if ((iaParam < 1) .OR. (iaParam > nPhaseSpecies) .OR. &
            (ibParam < 1) .OR. (ibParam > nPhaseSpecies)) cycle LOOP_PARAM_DIAG

        iExponentParam = iRegularParam(iParam,4)
        if (iExponentParam < 0) cycle LOOP_PARAM_DIAG

        dAnalyticParam = CompBinaryRKMPAnalyticDirectional(iParam,dMolesLocal,dDirection)
        dFDParam = (CompBinaryRKMPExcessEnergy(iParam,dMolesPlus) - &
            2D0*CompBinaryRKMPExcessEnergy(iParam,dMolesLocal) + &
            CompBinaryRKMPExcessEnergy(iParam,dMolesMinus)) / (dEps*dEps)
        dRelativeError = DABS(dFDParam-dAnalyticParam) / DMAX1(DABS(dAnalyticParam),1D-30)

        write(*,'(A,1X,F12.4,1X,A,1X,I5,1X,A,1X,I8,1X,A,1X,I5,1X,A,1X,I5,1X,A,1X,I5,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6)') &
            'RKMP_FD_PARAM_DEBUG T=', dTemperature, 'k=', iSolnIndex, 'param=', iParam, &
            'ia=', iaParam, 'ib=', ibParam, 'exp=', iExponentParam, 'L=', dExcessGibbsParam(iParam), &
            'eps=', dEps, 'analytic=', dAnalyticParam, 'fd=', dFDParam, 'relerr=', dRelativeError
    end do LOOP_PARAM_DIAG

    deallocate(dDirection,dMolesLocal,dMolesMinus,dMolesPlus)

contains

    real(8) function CompSupportedBinaryRKMPExcessEnergy(dMoles)

        implicit none

        real(8), intent(in), dimension(:) :: dMoles

        integer :: iParamLocal, iaLocal, ibLocal, iExponentLocal
        real(8) :: dNLocal, dDeltaLocal

        CompSupportedBinaryRKMPExcessEnergy = 0D0
        dNLocal = SUM(dMoles)
        if (dNLocal <= 1D-30) return

        do iParamLocal = nParamPhase(iSolnIndex-1)+1, nParamPhase(iSolnIndex)
            if (iRegularParam(iParamLocal,1) /= 2) cycle

            iaLocal = iRegularParam(iParamLocal,2)
            ibLocal = iRegularParam(iParamLocal,3)
            if ((iaLocal < 1) .OR. (iaLocal > SIZE(dMoles)) .OR. &
                (ibLocal < 1) .OR. (ibLocal > SIZE(dMoles))) cycle

            iExponentLocal = iRegularParam(iParamLocal,4)
            if (iExponentLocal < 0) cycle

            CompSupportedBinaryRKMPExcessEnergy = CompSupportedBinaryRKMPExcessEnergy + &
                CompBinaryRKMPExcessEnergy(iParamLocal,dMoles)
        end do

    end function CompSupportedBinaryRKMPExcessEnergy


    real(8) function CompBinaryRKMPExcessEnergy(iParamLocal,dMoles)

        implicit none

        integer, intent(in) :: iParamLocal
        real(8), intent(in), dimension(:) :: dMoles

        integer :: iaLocal, ibLocal, iExponentLocal
        real(8) :: dNLocal, dDeltaLocal

        CompBinaryRKMPExcessEnergy = 0D0
        dNLocal = SUM(dMoles)
        if (dNLocal <= 1D-30) return
        if (iRegularParam(iParamLocal,1) /= 2) return

        iaLocal = iRegularParam(iParamLocal,2)
        ibLocal = iRegularParam(iParamLocal,3)
        if ((iaLocal < 1) .OR. (iaLocal > SIZE(dMoles)) .OR. &
            (ibLocal < 1) .OR. (ibLocal > SIZE(dMoles))) return

        iExponentLocal = iRegularParam(iParamLocal,4)
        if (iExponentLocal < 0) return

        dDeltaLocal = (dMoles(iaLocal) - dMoles(ibLocal)) / dNLocal
        CompBinaryRKMPExcessEnergy = dExcessGibbsParam(iParamLocal) * &
            dMoles(iaLocal) * dMoles(ibLocal) / dNLocal * dDeltaLocal**iExponentLocal

    end function CompBinaryRKMPExcessEnergy


    real(8) function CompBinaryRKMPAnalyticDirectional(iParamLocal,dMoles,dDirectionLocal)

        implicit none

        integer, intent(in) :: iParamLocal
        real(8), intent(in), dimension(:) :: dMoles, dDirectionLocal

        integer :: iLocal, jLocal, iaLocal, ibLocal, iExponentLocal
        real(8) :: dNLocal, dNiLocal, dNjLocal, dBLocal, dCLocal, dALocal, dDeltaLocal
        real(8) :: dL0Local, dL1Local, dL2Local, dHijLocal
        real(8) :: dDelta_ai, dDelta_aj, dDelta_bi, dDelta_bj
        real(8) :: dBa, dBb, dBab, dCa, dCb, dAa, dAb, dAab, dDa, dDb, dDab

        CompBinaryRKMPAnalyticDirectional = 0D0
        dNLocal = SUM(dMoles)
        if (dNLocal <= 1D-30) return
        if (iRegularParam(iParamLocal,1) /= 2) return

        iaLocal = iRegularParam(iParamLocal,2)
        ibLocal = iRegularParam(iParamLocal,3)
        if ((iaLocal < 1) .OR. (iaLocal > SIZE(dMoles)) .OR. &
            (ibLocal < 1) .OR. (ibLocal > SIZE(dMoles))) return

        iExponentLocal = iRegularParam(iParamLocal,4)
        if (iExponentLocal < 0) return

        dNiLocal = dMoles(iaLocal)
        dNjLocal = dMoles(ibLocal)
        dBLocal = dNiLocal * dNjLocal
        dCLocal = dNiLocal - dNjLocal
        dALocal = dBLocal / dNLocal
        dDeltaLocal = dCLocal / dNLocal

        dL0Local = dDeltaLocal**iExponentLocal
        if (iExponentLocal == 0) then
            dL1Local = 0D0
            dL2Local = 0D0
        elseif (iExponentLocal == 1) then
            dL1Local = 1D0
            dL2Local = 0D0
        else
            dL1Local = DFLOAT(iExponentLocal) * dDeltaLocal**(iExponentLocal-1)
            if (iExponentLocal == 2) then
                dL2Local = 2D0
            else
                dL2Local = DFLOAT(iExponentLocal*(iExponentLocal-1)) * &
                    dDeltaLocal**(iExponentLocal-2)
            end if
        end if

        do iLocal = 1, SIZE(dMoles)
            do jLocal = 1, SIZE(dMoles)
                dDelta_ai = 0D0
                dDelta_aj = 0D0
                dDelta_bi = 0D0
                dDelta_bj = 0D0
                if (iLocal == iaLocal) dDelta_ai = 1D0
                if (iLocal == ibLocal) dDelta_aj = 1D0
                if (jLocal == iaLocal) dDelta_bi = 1D0
                if (jLocal == ibLocal) dDelta_bj = 1D0

                dBa  = dDelta_ai*dNjLocal + dDelta_aj*dNiLocal
                dBb  = dDelta_bi*dNjLocal + dDelta_bj*dNiLocal
                dBab = dDelta_ai*dDelta_bj + dDelta_aj*dDelta_bi

                dCa  = dDelta_ai - dDelta_aj
                dCb  = dDelta_bi - dDelta_bj

                dAa  = dBa/dNLocal - dBLocal/(dNLocal*dNLocal)
                dAb  = dBb/dNLocal - dBLocal/(dNLocal*dNLocal)
                dAab = dBab/dNLocal - dBa/(dNLocal*dNLocal) - &
                    dBb/(dNLocal*dNLocal) + 2D0*dBLocal/(dNLocal*dNLocal*dNLocal)

                dDa  = dCa/dNLocal - dCLocal/(dNLocal*dNLocal)
                dDb  = dCb/dNLocal - dCLocal/(dNLocal*dNLocal)
                dDab = -dCa/(dNLocal*dNLocal) - dCb/(dNLocal*dNLocal) + &
                    2D0*dCLocal/(dNLocal*dNLocal*dNLocal)

                dHijLocal = dExcessGibbsParam(iParamLocal) * ( dAab*dL0Local + &
                    (dAa*dDb + dAb*dDa + dALocal*dDab)*dL1Local + &
                    dALocal*dDa*dDb*dL2Local )

                CompBinaryRKMPAnalyticDirectional = CompBinaryRKMPAnalyticDirectional + &
                    dDirectionLocal(iLocal) * dHijLocal * dDirectionLocal(jLocal)
            end do
        end do

    end function CompBinaryRKMPAnalyticDirectional

end subroutine DebugRKMPHessianFiniteDifference
