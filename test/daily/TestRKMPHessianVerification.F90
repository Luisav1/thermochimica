!-------------------------------------------------------------------------------------------------------------
!> \file    TestRKMPHessianVerification.F90
!> \brief   Thermochimica-native finite-difference verification of the plain-RKMP excess Hessian.
!>
!> \details Reproduce the TestThermo30 equilibrium, then compare the production analytic RKMP Hessian with
!!          scalar-energy and established partial-molar finite differences at the converged phase state.
!!          Pass --report to print the precision and convergence tables used for numerical evidence.
!-------------------------------------------------------------------------------------------------------------

program TestRKMPHessianVerification

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    interface
        subroutine CompExcessGibbsEnergyRKMP_unconstrained(iSolnIndex,dHess)
            integer, intent(in)                  :: iSolnIndex
            real(8), intent(out), dimension(:,:) :: dHess
        end subroutine CompExcessGibbsEnergyRKMP_unconstrained

        subroutine CompExcessGibbsEnergyRKMP(iSolnIndex)
            integer, intent(in) :: iSolnIndex
        end subroutine CompExcessGibbsEnergyRKMP

        subroutine CompRKMPBinaryExcessGibbsFromMoles(iSolnIndex,nLocalSpecies,dLocalMoles,dGex)
            integer, intent(in) :: iSolnIndex, nLocalSpecies
            real(8), dimension(nLocalSpecies), intent(in) :: dLocalMoles
            real(8), intent(out) :: dGex
        end subroutine CompRKMPBinaryExcessGibbsFromMoles
    end interface

    integer, parameter :: nFDSteps = 5
    integer :: i, iDirection, iEnergyA, iEnergyB, iFirstSpecies, iLastSpecies, iParam
    integer :: iPhaseIndex, iReference
    integer :: iSolnSlot, iStep, nDirections, nLocalSpecies
    integer :: nStorageBits, nDecimalPrecision, nBinaryDigits
    real(8) :: dAnalyticEx, dAnalyticIdeal, dEnergy0Ex, dEnergy0Ideal
    real(8) :: dEnergyM1Ex, dEnergyM2Ex, dEnergyP1Ex, dEnergyP2Ex
    real(8) :: dEnergyM1Ideal, dEnergyM2Ideal, dEnergyP1Ideal, dEnergyP2Ideal
    real(8) :: dEpsilonMachine, dFD3, dFD5, dH, dScale, dTotalMoles
    real(8) :: dRadialResidual, dSymmetryResidual, dScaleHessian
    real(8) :: dOrder3, dOrder5, dWorstMuBest
    logical :: lPass, lReport
    character(len=32) :: cArgument
    real(8), allocatable, dimension(:) :: dDirection, dMoles0, dMolesM1, dMolesM2
    real(8), allocatable, dimension(:) :: dMolesP1, dMolesP2, dMuAnalytic, dMuMinus, dMuPlus
    real(8), allocatable, dimension(:) :: dMolFractionSave, dPartialExcessSave, dX
    real(8), allocatable, dimension(:) :: dEnergyStep, dErr3Ex, dErr3Ideal, dErr5Ex, dErr5Ideal
    real(8), allocatable, dimension(:,:) :: dErrMu, dHessian, dMuStep

    lPass = .TRUE.
    lReport = .FALSE.
    if (COMMAND_ARGUMENT_COUNT() > 0) then
        call GET_COMMAND_ARGUMENT(1,cArgument)
        lReport = TRIM(cArgument) == '--report'
    end if

    nStorageBits = STORAGE_SIZE(1D0)
    nDecimalPrecision = PRECISION(1D0)
    nBinaryDigits = DIGITS(1D0)
    dEpsilonMachine = EPSILON(1D0)
    lPass = lPass .AND. (nStorageBits == 64)
    lPass = lPass .AND. (nDecimalPrecision >= 15)
    lPass = lPass .AND. (nBinaryDigits >= 53)
    lPass = lPass .AND. (dEpsilonMachine <= 3D-16)

    cInputUnitTemperature = 'K'
    cInputUnitPressure = 'atm'
    cInputUnitMass = 'moles'
    cThermoFileName = DATA_DIRECTORY // 'WAuArO-1.dat'
    dPressure = 1D0
    dTemperature = 1455D0
    dElementMass(74) = 1.95D0
    dElementMass(79) = 1D0
    dElementMass(18) = 2D0
    dElementMass(8) = 10D0

    call ParseCSDataFile(cThermoFileName)
    call Thermochimica
    lPass = lPass .AND. (INFOThermo == 0)
    lPass = lPass .AND. (DABS((dGibbsEnergySys+4.620D5)/4.620D5) < 1D-3)

    iPhaseIndex = 0
    if (INFOThermo == 0) then
        do iSolnSlot = 1, nSolnPhases
            i = -iAssemblage(nElements-iSolnSlot+1)
            if (i <= 0) cycle
            if (cSolnPhaseType(i) == 'RKMP') then
                iPhaseIndex = i
                exit
            end if
        end do
    end if
    lPass = lPass .AND. (iPhaseIndex > 0)

    if (iPhaseIndex > 0) then
        iFirstSpecies = nSpeciesPhase(iPhaseIndex-1) + 1
        iLastSpecies = nSpeciesPhase(iPhaseIndex)
        nLocalSpecies = iLastSpecies - iFirstSpecies + 1
        nDirections = nLocalSpecies - 1
        iReference = nLocalSpecies
        iEnergyA = 0
        iEnergyB = 0
        do iParam = nParamPhase(iPhaseIndex-1)+1, nParamPhase(iPhaseIndex)
            if (iRegularParam(iParam,1) /= 2) cycle
            if (iRegularParam(iParam,4) < 0) cycle
            iEnergyA = iRegularParam(iParam,2)
            iEnergyB = iRegularParam(iParam,3)
            exit
        end do
        lPass = lPass .AND. (iEnergyA > 0) .AND. (iEnergyB > 0)

        allocate(dDirection(nLocalSpecies), dMoles0(nLocalSpecies), dMolesM1(nLocalSpecies), &
            dMolesM2(nLocalSpecies), dMolesP1(nLocalSpecies), dMolesP2(nLocalSpecies), &
            dMuAnalytic(nLocalSpecies), dMuMinus(nLocalSpecies), dMuPlus(nLocalSpecies), &
            dMolFractionSave(nLocalSpecies), dPartialExcessSave(nLocalSpecies), dX(nLocalSpecies), &
            dEnergyStep(nFDSteps), dErr3Ex(nFDSteps), dErr3Ideal(nFDSteps), &
            dErr5Ex(nFDSteps), dErr5Ideal(nFDSteps), &
            dErrMu(nDirections,nFDSteps), dHessian(nLocalSpecies,nLocalSpecies), &
            dMuStep(nDirections,nFDSteps))

        dMoles0 = dMolesSpecies(iFirstSpecies:iLastSpecies)
        dTotalMoles = SUM(dMoles0)
        dX = dMoles0 / dTotalMoles
        call CompExcessGibbsEnergyRKMP_unconstrained(iPhaseIndex,dHessian)

        lPass = lPass .AND. ALL(IEEE_IS_FINITE(dHessian))
        dScaleHessian = DMAX1(1D0,MAXVAL(DABS(dHessian)))
        dSymmetryResidual = MAXVAL(DABS(dHessian-TRANSPOSE(dHessian))) / dScaleHessian
        dRadialResidual = MAXVAL(DABS(MATMUL(dHessian,dMoles0))) / &
            DMAX1(1D0,MAXVAL(DABS(dHessian))*SUM(DABS(dMoles0)))
        lPass = lPass .AND. (dSymmetryResidual <= 1D-12)
        lPass = lPass .AND. (dRadialResidual <= 1D-10)

        dMolFractionSave = dMolFraction(iFirstSpecies:iLastSpecies)
        dPartialExcessSave = dPartialExcessGibbs(iFirstSpecies:iLastSpecies)

        ! Use the actual binary interaction pair for scalar-energy differences.  Both species have substantial
        ! phase amounts in TestThermo30, leaving a visible truncation region before roundoff dominates.
        dDirection = 0D0
        dDirection(iEnergyA) = 1D0
        dDirection(iEnergyB) = -1D0
        dScale = DMIN1(dTotalMoles,DMIN1(dMoles0(iEnergyA),dMoles0(iEnergyB)))
        dAnalyticEx = DOT_PRODUCT(dDirection,MATMUL(dHessian,dDirection))
        dAnalyticIdeal = 1D0/dMoles0(iEnergyA) + 1D0/dMoles0(iEnergyB)
        call CompRKMPBinaryExcessGibbsFromMoles(iPhaseIndex,nLocalSpecies,dMoles0,dEnergy0Ex)
        dEnergy0Ideal = CompIdealMixingEnergy(dMoles0)

        do iStep = 1, nFDSteps
            dH = dScale * 10D0**(-iStep)
            dEnergyStep(iStep) = dH
            dMolesM1 = dMoles0 - dH*dDirection
            dMolesP1 = dMoles0 + dH*dDirection
            dMolesM2 = dMoles0 - 2D0*dH*dDirection
            dMolesP2 = dMoles0 + 2D0*dH*dDirection

            call CompRKMPBinaryExcessGibbsFromMoles(iPhaseIndex,nLocalSpecies,dMolesM1,dEnergyM1Ex)
            call CompRKMPBinaryExcessGibbsFromMoles(iPhaseIndex,nLocalSpecies,dMolesP1,dEnergyP1Ex)
            call CompRKMPBinaryExcessGibbsFromMoles(iPhaseIndex,nLocalSpecies,dMolesM2,dEnergyM2Ex)
            call CompRKMPBinaryExcessGibbsFromMoles(iPhaseIndex,nLocalSpecies,dMolesP2,dEnergyP2Ex)
            dEnergyM1Ideal = CompIdealMixingEnergy(dMolesM1)
            dEnergyP1Ideal = CompIdealMixingEnergy(dMolesP1)
            dEnergyM2Ideal = CompIdealMixingEnergy(dMolesM2)
            dEnergyP2Ideal = CompIdealMixingEnergy(dMolesP2)

            dFD3 = (dEnergyP1Ex-2D0*dEnergy0Ex+dEnergyM1Ex)/(dH*dH)
            dFD5 = (-dEnergyP2Ex+16D0*dEnergyP1Ex-30D0*dEnergy0Ex+ &
                16D0*dEnergyM1Ex-dEnergyM2Ex)/(12D0*dH*dH)
            dErr3Ex(iStep) = CompScaledError(dFD3,dAnalyticEx)
            dErr5Ex(iStep) = CompScaledError(dFD5,dAnalyticEx)

            dFD3 = (dEnergyP1Ideal-2D0*dEnergy0Ideal+dEnergyM1Ideal)/(dH*dH)
            dFD5 = (-dEnergyP2Ideal+16D0*dEnergyP1Ideal-30D0*dEnergy0Ideal+ &
                16D0*dEnergyM1Ideal-dEnergyM2Ideal)/(12D0*dH*dH)
            dErr3Ideal(iStep) = CompScaledError(dFD3,dAnalyticIdeal)
            dErr5Ideal(iStep) = CompScaledError(dFD5,dAnalyticIdeal)
        end do

        lPass = lPass .AND. (MINVAL(dErr3Ex) <= 1D-7)
        lPass = lPass .AND. (MINVAL(dErr5Ex) <= 1D-9)
        dOrder3 = CompObservedOrder(dErr3Ideal(1),dErr3Ideal(2),dEnergyStep(1),dEnergyStep(2))
        dOrder5 = CompObservedOrder(dErr5Ideal(1),dErr5Ideal(2),dEnergyStep(1),dEnergyStep(2))
        lPass = lPass .AND. (dOrder3 >= 1.5D0) .AND. (dOrder3 <= 2.5D0)
        lPass = lPass .AND. (dOrder5 >= 3D0) .AND. (dOrder5 <= 5D0)

        ! The production partial-molar comparison is numerically stable enough to cover the full tangent basis,
        ! including directions involving trace species whose scalar-energy differences are cancellation-limited.
        do iDirection = 1, nDirections
            dDirection = 0D0
            dDirection(iDirection) = 1D0
            dDirection(iReference) = -1D0
            dScale = DMIN1(dTotalMoles,0.25D0*DMIN1(dMoles0(iDirection),dMoles0(iReference)))
            lPass = lPass .AND. (dScale > 0D0)

            dMuAnalytic = MATMUL(dTotalMoles*dHessian,dDirection)

            do iStep = 1, nFDSteps
                dH = dScale * 10D0**(-iStep)
                dMuStep(iDirection,iStep) = dH

                dMolFraction(iFirstSpecies:iLastSpecies) = dX + (dH/dTotalMoles)*dDirection
                dPartialExcessGibbs(iFirstSpecies:iLastSpecies) = 0D0
                call CompExcessGibbsEnergyRKMP(iPhaseIndex)
                dMuPlus = dPartialExcessGibbs(iFirstSpecies:iLastSpecies)

                dMolFraction(iFirstSpecies:iLastSpecies) = dX - (dH/dTotalMoles)*dDirection
                dPartialExcessGibbs(iFirstSpecies:iLastSpecies) = 0D0
                call CompExcessGibbsEnergyRKMP(iPhaseIndex)
                dMuMinus = dPartialExcessGibbs(iFirstSpecies:iLastSpecies)

                dErrMu(iDirection,iStep) = MAXVAL(DABS((dMuPlus-dMuMinus)/(2D0*dH/dTotalMoles)- &
                    dMuAnalytic)) / DMAX1(1D0,MAXVAL(DABS(dMuAnalytic)))
            end do

            lPass = lPass .AND. (MINVAL(dErrMu(iDirection,:)) <= 1D-9)
        end do

        dMolFraction(iFirstSpecies:iLastSpecies) = dMolFractionSave
        dPartialExcessGibbs(iFirstSpecies:iLastSpecies) = dPartialExcessSave

        dWorstMuBest = 0D0
        do iDirection = 1, nDirections
            dWorstMuBest = DMAX1(dWorstMuBest,MINVAL(dErrMu(iDirection,:)))
        end do

        if (lReport) call PrintReport

        deallocate(dDirection,dMoles0,dMolesM1,dMolesM2,dMolesP1,dMolesP2,dMuAnalytic, &
            dMuMinus,dMuPlus,dMolFractionSave,dPartialExcessSave,dX,dEnergyStep,dErr3Ex,dErr3Ideal, &
            dErr5Ex,dErr5Ideal,dErrMu,dHessian,dMuStep)
    end if

    if (lPass) then
        print *, 'TestRKMPHessianVerification: PASS'
        call ResetThermo
        call EXIT(0)
    else
        print *, 'TestRKMPHessianVerification: FAIL <---'
        call ResetThermo
        call EXIT(1)
    end if

contains

    real(8) function CompIdealMixingEnergy(dMoles)
        real(8), intent(in), dimension(:) :: dMoles
        integer :: iLocal
        real(8) :: dNLocal

        CompIdealMixingEnergy = 0D0
        dNLocal = SUM(dMoles)
        do iLocal = 1, SIZE(dMoles)
            if (dMoles(iLocal) > 0D0) then
                CompIdealMixingEnergy = CompIdealMixingEnergy + &
                    dMoles(iLocal)*DLOG(dMoles(iLocal)/dNLocal)
            end if
        end do
    end function CompIdealMixingEnergy

    real(8) function CompScaledError(dApproximate,dExact)
        real(8), intent(in) :: dApproximate, dExact
        CompScaledError = DABS(dApproximate-dExact) / DMAX1(DABS(dExact),1D-12)
    end function CompScaledError

    real(8) function CompObservedOrder(dErrorCoarse,dErrorFine,dHCoarse,dHFine)
        real(8), intent(in) :: dErrorCoarse, dErrorFine, dHCoarse, dHFine
        CompObservedOrder = DLOG(dErrorCoarse/dErrorFine) / DLOG(dHCoarse/dHFine)
    end function CompObservedOrder

    subroutine PrintReport
        integer :: iDirLocal, iStepLocal
        real(8) :: dOrder3Local, dOrder5Local

        write(*,'(A)') 'RKMP Hessian verification at converged TestThermo30 state'
        write(*,'(A,I0)') 'storage bits      = ', nStorageBits
        write(*,'(A,I0)') 'decimal precision = ', nDecimalPrecision
        write(*,'(A,I0)') 'binary digits     = ', nBinaryDigits
        write(*,'(A,ES14.6)') 'machine epsilon   = ', dEpsilonMachine
        write(*,'(A,ES14.6)') 'symmetry residual = ', dSymmetryResidual
        write(*,'(A,ES14.6)') 'radial residual   = ', dRadialResidual

        write(*,'(/,A,I0,A,I0)') 'energy direction: species ', iEnergyA, ' minus species ', iEnergyB
        write(*,'(A)') 'h                  RKMP-3pt       RKMP-5pt       ideal-3pt      ideal-5pt'
        do iStepLocal = 1, nFDSteps
            write(*,'(ES14.6,4(2X,ES14.6))') dEnergyStep(iStepLocal), dErr3Ex(iStepLocal), &
                dErr5Ex(iStepLocal), dErr3Ideal(iStepLocal), dErr5Ideal(iStepLocal)
        end do
        dOrder3Local = CompObservedOrder(dErr3Ideal(1),dErr3Ideal(2),dEnergyStep(1),dEnergyStep(2))
        dOrder5Local = CompObservedOrder(dErr5Ideal(1),dErr5Ideal(2),dEnergyStep(1),dEnergyStep(2))
        write(*,'(A,F8.4,A,F8.4)') 'ideal observed orders: 3pt=', dOrder3Local, '  5pt=', dOrder5Local

        write(*,'(/,A)') 'production partial-molar RKMP comparison'
        do iDirLocal = 1, nDirections
            write(*,'(/,A,I0,A,I0)') 'direction: species ', iDirLocal, ' minus species ', iReference
            write(*,'(A)') 'h                  production-mu scaled error'
            do iStepLocal = 1, nFDSteps
                write(*,'(ES14.6,2X,ES14.6)') dMuStep(iDirLocal,iStepLocal), dErrMu(iDirLocal,iStepLocal)
            end do
        end do
        write(*,'(/,A,ES14.6)') 'worst best production-mu error = ', dWorstMuBest
    end subroutine PrintReport

end program TestRKMPHessianVerification
