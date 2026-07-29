!-------------------------------------------------------------------------------------------------------------
!> \file    TestRKMPHessianVerification.F90
!> \brief   Thermochimica-native finite-difference verification of the plain-RKMP excess Hessian.
!>
!> \details Reproduce the TestThermo30 equilibrium, then compare the production analytic RKMP Hessian with
!!          scalar-energy and established partial-molar finite differences at the converged phase state.
!!          Pass --report to print the precision and convergence tables used for numerical evidence.
!!
!!          Verification map:
!!          1. Run the ordinary TestThermo30 thermodynamic calculation with the
!!             curvature-enabled solver behavior left at its default-off setting.
!!          2. Identify the converged active plain-RKMP phase and evaluate Hloc.
!!          3. Check double precision, finite values, symmetry, and the
!!             extensivity identity: uniformly scaling every species amount
!!             changes phase amount but not composition, so Hloc applied to the
!!             current mole vector should be zero.
!!          4. Transfer moles between two interacting species along direction v.
!!             Compare the energy curvature predicted by Hloc along that
!!             direction with three- and five-point scalar-energy differences.
!!          5. Apply the same stencils to ideal mixing as a truncation/roundoff
!!             control with known convergence order.
!!          6. For each independent composition direction, use Hloc to predict
!!             how all excess partial molars change, then compare with finite
!!             differences of the established production RKMP routine.
!!
!!          This test verifies the local production-linked Hessian. Solver
!!          response mapping and nonlinear hardening are covered separately by
!!          diagnostics and the complete regression suite.
!-------------------------------------------------------------------------------------------------------------

program TestRKMPHessianVerification

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleFiniteDifferenceVerification

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

    integer, parameter :: nEnergyFDSteps = 12
    integer, parameter :: nMuFDSteps = 9
    integer :: i, iDirection, iEnergyA, iEnergyB, iFirstSpecies, iLastSpecies, iParam
    integer :: iPhaseIndex, iReference, iControlledParam, iExponentSave
    integer :: iSolnSlot, iStep, nDirections, nLocalSpecies
    integer :: nStorageBits, nDecimalPrecision, nBinaryDigits
    real(8) :: dAnalyticEx, dAnalyticIdeal, dEnergy0Ex, dEnergy0Ideal
    real(8) :: dEnergyM1Ex, dEnergyM2Ex, dEnergyP1Ex, dEnergyP2Ex
    real(8) :: dEnergyM1Ideal, dEnergyM2Ideal, dEnergyP1Ideal, dEnergyP2Ideal
    real(8) :: dEpsilonMachine, dFDForward, dFDBackward, dFD3, dFD5, dH, dScale, dTotalMoles
    real(8) :: dRadialResidual, dSymmetryResidual, dScaleHessian
    real(8) :: dWorstMuBest
    logical :: lPass, lReport
    character(len=32) :: cArgument
    real(8), allocatable, dimension(:) :: dDirection, dMoles0, dMolesM1, dMolesM2
    real(8), allocatable, dimension(:) :: dMolesP1, dMolesP2, dMuAnalytic, dMuMinus, dMuPlus
    real(8), allocatable, dimension(:) :: dMolFractionSave, dPartialExcessSave, dX
    real(8), allocatable, dimension(:) :: dEnergyStep, dErrForwardEx, dErrBackwardEx
    real(8), allocatable, dimension(:) :: dErr3Ex, dErr3Ideal, dErr5Ex, dErr5Ideal
    real(8), allocatable, dimension(:) :: dAbsForwardEx, dAbsBackwardEx
    real(8), allocatable, dimension(:) :: dAbs3Ex, dAbs3Ideal, dAbs5Ex, dAbs5Ideal
    real(8), allocatable, dimension(:) :: dOrderForwardEx, dOrderBackwardEx
    real(8), allocatable, dimension(:) :: dOrder3Ex, dOrder3Ideal, dOrder5Ex, dOrder5Ideal
    logical, allocatable, dimension(:) :: lOrderForwardEx, lOrderBackwardEx
    logical, allocatable, dimension(:) :: lOrder3Ex, lOrder3Ideal, lOrder5Ex, lOrder5Ideal
    real(8), allocatable, dimension(:,:) :: dErrMu, dHessian, dMuStep, dNormAbsMu
    real(8), allocatable, dimension(:,:) :: dMaxAbsMu, dMaxScaledMu, dOrderMu
    real(8), allocatable, dimension(:,:,:) :: dMuFD
    logical, allocatable, dimension(:,:) :: lOrderMu
    integer, allocatable, dimension(:,:) :: iWorstMu
    type(FDSweepAssessment) :: tRKMPForward, tRKMPBackward, tRKMP3, tRKMP5, tIdeal3, tIdeal5

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

    !=========================================================================================================
    ! SECTION 1: CONVERGED THERMOCHIMICA RKMP STATE
    !=========================================================================================================
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
        iControlledParam = 0
        do iParam = nParamPhase(iPhaseIndex-1)+1, nParamPhase(iPhaseIndex)
            if (iRegularParam(iParam,1) /= 2) cycle
            if (iRegularParam(iParam,4) < 0) cycle
            iEnergyA = iRegularParam(iParam,2)
            iEnergyB = iRegularParam(iParam,3)
            iControlledParam = iParam
            exit
        end do
        lPass = lPass .AND. (iEnergyA > 0) .AND. (iEnergyB > 0)

        allocate(dDirection(nLocalSpecies), dMoles0(nLocalSpecies), dMolesM1(nLocalSpecies), &
            dMolesM2(nLocalSpecies), dMolesP1(nLocalSpecies), dMolesP2(nLocalSpecies), &
            dMuAnalytic(nLocalSpecies), dMuMinus(nLocalSpecies), dMuPlus(nLocalSpecies), &
            dMolFractionSave(nLocalSpecies), dPartialExcessSave(nLocalSpecies), dX(nLocalSpecies), &
            dEnergyStep(nEnergyFDSteps), dErrForwardEx(nEnergyFDSteps), &
            dErrBackwardEx(nEnergyFDSteps), dErr3Ex(nEnergyFDSteps), dErr3Ideal(nEnergyFDSteps), &
            dErr5Ex(nEnergyFDSteps), dErr5Ideal(nEnergyFDSteps), dAbs3Ex(nEnergyFDSteps), &
            dAbs3Ideal(nEnergyFDSteps), dAbs5Ex(nEnergyFDSteps), dAbs5Ideal(nEnergyFDSteps), &
            dAbsForwardEx(nEnergyFDSteps), dAbsBackwardEx(nEnergyFDSteps), &
            dOrderForwardEx(nEnergyFDSteps), dOrderBackwardEx(nEnergyFDSteps), &
            dOrder3Ex(nEnergyFDSteps), dOrder3Ideal(nEnergyFDSteps), &
            dOrder5Ex(nEnergyFDSteps), dOrder5Ideal(nEnergyFDSteps), &
            lOrder3Ex(nEnergyFDSteps), lOrder3Ideal(nEnergyFDSteps), &
            lOrder5Ex(nEnergyFDSteps), lOrder5Ideal(nEnergyFDSteps), &
            lOrderForwardEx(nEnergyFDSteps), lOrderBackwardEx(nEnergyFDSteps), &
            dErrMu(nDirections,nMuFDSteps), dHessian(nLocalSpecies,nLocalSpecies), &
            dMuStep(nDirections,nMuFDSteps), dNormAbsMu(nDirections,nMuFDSteps), &
            dMaxAbsMu(nDirections,nMuFDSteps), dMaxScaledMu(nDirections,nMuFDSteps), &
            dOrderMu(nDirections,nMuFDSteps), lOrderMu(nDirections,nMuFDSteps), &
            iWorstMu(nDirections,nMuFDSteps), dMuFD(nDirections,nMuFDSteps,nLocalSpecies))

        dMoles0 = dMolesSpecies(iFirstSpecies:iLastSpecies)
        dTotalMoles = SUM(dMoles0)
        dX = dMoles0 / dTotalMoles
        !=====================================================================================================
        ! SECTION 2: STRUCTURAL HESSIAN IDENTITIES
        !
        ! Symmetry checks that changing species i then j gives the same mixed
        ! derivative as changing j then i. Extensivity means uniformly scaling
        ! every species amount changes phase amount but not composition; Hloc
        ! applied to the current mole vector should therefore be zero.
        !=====================================================================================================
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

        !=====================================================================================================
        ! SECTION 3: SCALAR-ENERGY DIRECTIONAL CURVATURE
        !
        ! Use the actual binary interaction pair for scalar-energy differences. Both species have substantial
        ! phase amounts in TestThermo30, leaving a visible truncation region before roundoff dominates.
        !=====================================================================================================
        ! TestThermo30 contains only a constant binary interaction, whose
        ! directional scalar energy is too low-order to display stencil
        ! truncation. Temporarily exercise the same parsed parameter as an
        ! exponent-four controlled fixture. Both production-linked analytic and
        ! independent scalar paths consume this state.
        iExponentSave = iRegularParam(iControlledParam,4)
        iRegularParam(iControlledParam,4) = 4
        call CompExcessGibbsEnergyRKMP_unconstrained(iPhaseIndex,dHessian)

        dDirection = 0D0
        dDirection(iEnergyA) = 1D0
        dDirection(iEnergyB) = -1D0
        dScale = DMIN1(dTotalMoles,DMIN1(dMoles0(iEnergyA),dMoles0(iEnergyB)))
        dAnalyticEx = DOT_PRODUCT(dDirection,MATMUL(dHessian,dDirection))
        dAnalyticIdeal = 1D0/dMoles0(iEnergyA) + 1D0/dMoles0(iEnergyB)
        call CompRKMPBinaryExcessGibbsFromMoles(iPhaseIndex,nLocalSpecies,dMoles0,dEnergy0Ex)
        dEnergy0Ideal = CompIdealMixingEnergy(dMoles0)

        do iStep = 1, nEnergyFDSteps
            dH = 0.2D0*dScale * 3D0**(-(iStep-1))
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

            ! All four stencils estimate the same curvature v^T H v from
            ! independent scalar-energy evaluations. Forward and backward
            ! differences retain leading +h*f''' and -h*f''' terms, respectively,
            ! and should lose absolute error as h. Their absolute-error curves
            ! may therefore overlap even though the signed errors oppose.
            ! Centering cancels that term, giving h**2 behavior; the wider
            ! five-point stencil cancels more terms and gives h**4.
            ! At very small h, subtracting nearly equal energies amplifies
            ! floating-point roundoff, so an eventual error increase is expected.
            dFDForward = (dEnergyP2Ex-2D0*dEnergyP1Ex+dEnergy0Ex)/(dH*dH)
            dFDBackward = (dEnergy0Ex-2D0*dEnergyM1Ex+dEnergyM2Ex)/(dH*dH)
            dFD3 = (dEnergyP1Ex-2D0*dEnergy0Ex+dEnergyM1Ex)/(dH*dH)
            dFD5 = (-dEnergyP2Ex+16D0*dEnergyP1Ex-30D0*dEnergy0Ex+ &
                16D0*dEnergyM1Ex-dEnergyM2Ex)/(12D0*dH*dH)
            dAbsForwardEx(iStep) = DABS(dFDForward-dAnalyticEx)
            dAbsBackwardEx(iStep) = DABS(dFDBackward-dAnalyticEx)
            dAbs3Ex(iStep) = DABS(dFD3-dAnalyticEx)
            dAbs5Ex(iStep) = DABS(dFD5-dAnalyticEx)
            dErrForwardEx(iStep) = CompScaledError(dFDForward,dAnalyticEx)
            dErrBackwardEx(iStep) = CompScaledError(dFDBackward,dAnalyticEx)
            dErr3Ex(iStep) = CompScaledError(dFD3,dAnalyticEx)
            dErr5Ex(iStep) = CompScaledError(dFD5,dAnalyticEx)

            dFD3 = (dEnergyP1Ideal-2D0*dEnergy0Ideal+dEnergyM1Ideal)/(dH*dH)
            dFD5 = (-dEnergyP2Ideal+16D0*dEnergyP1Ideal-30D0*dEnergy0Ideal+ &
                16D0*dEnergyM1Ideal-dEnergyM2Ideal)/(12D0*dH*dH)
            dAbs3Ideal(iStep) = DABS(dFD3-dAnalyticIdeal)
            dAbs5Ideal(iStep) = DABS(dFD5-dAnalyticIdeal)
            dErr3Ideal(iStep) = CompScaledError(dFD3,dAnalyticIdeal)
            dErr5Ideal(iStep) = CompScaledError(dFD5,dAnalyticIdeal)
        end do

        call AssessFDSweep(dEnergyStep,dErrForwardEx,FD_ORDER_FIRST_MIN,FD_ORDER_FIRST_MAX,1D-6, &
            tRKMPForward,dOrderForwardEx,lOrderForwardEx)
        call AssessFDSweep(dEnergyStep,dErrBackwardEx,FD_ORDER_FIRST_MIN,FD_ORDER_FIRST_MAX,1D-6, &
            tRKMPBackward,dOrderBackwardEx,lOrderBackwardEx)
        call AssessFDSweep(dEnergyStep,dErr3Ex,FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX,1D-7, &
            tRKMP3,dOrder3Ex,lOrder3Ex)
        call AssessFDSweep(dEnergyStep,dErr5Ex,FD_ORDER_FOURTH_MIN,FD_ORDER_FOURTH_MAX,1D-9, &
            tRKMP5,dOrder5Ex,lOrder5Ex)
        call AssessFDSweep(dEnergyStep,dErr3Ideal,FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX,1D-6, &
            tIdeal3,dOrder3Ideal,lOrder3Ideal)
        call AssessFDSweep(dEnergyStep,dErr5Ideal,FD_ORDER_FOURTH_MIN,FD_ORDER_FOURTH_MAX,1D-8, &
            tIdeal5,dOrder5Ideal,lOrder5Ideal)
        lPass = lPass .AND. tRKMPForward%lPassed .AND. tRKMPBackward%lPassed
        lPass = lPass .AND. tRKMP3%lPassed .AND. tRKMP5%lPassed
        lPass = lPass .AND. tIdeal3%lPassed .AND. tIdeal5%lPassed

        ! The order fixture ends here. Restore the parsed TestThermo30 parameter
        ! and rebuild Hloc so the following production comparison is genuinely
        ! native to the database rather than another controlled-exponent test.
        iRegularParam(iControlledParam,4) = iExponentSave
        call CompExcessGibbsEnergyRKMP_unconstrained(iPhaseIndex,dHessian)

        !=====================================================================================================
        ! SECTION 4: NATIVE PRODUCTION PARTIAL-MOLAR DERIVATIVE ORACLE
        !
        ! This section uses the unmodified TestThermo30 database parameter. Its
        ! low-order RKMP polynomial can make a centered derivative exact up to
        ! roundoff from the coarsest step onward, so this is an accuracy and
        ! implementation-consistency check rather than an observed-order test.
        ! The controlled scalar fixture above supplies the truncation-order
        ! evidence.
        !=====================================================================================================
        dOrderMu = 0D0
        lOrderMu = .FALSE.
        do iDirection = 1, nDirections
            dDirection = 0D0
            dDirection(iDirection) = 1D0
            dDirection(iReference) = -1D0
            dScale = DMIN1(dTotalMoles,0.25D0*DMIN1(dMoles0(iDirection),dMoles0(iReference)))
            lPass = lPass .AND. (dScale > 0D0)

            dMuAnalytic = MATMUL(dTotalMoles*dHessian,dDirection)

            do iStep = 1, nMuFDSteps
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

                dMuFD(iDirection,iStep,:) = (dMuPlus-dMuMinus)/(2D0*dH/dTotalMoles)
                call ComputeVectorErrorMetrics(dMuFD(iDirection,iStep,:),dMuAnalytic, &
                    dNormAbsMu(iDirection,iStep),dErrMu(iDirection,iStep), &
                    dMaxAbsMu(iDirection,iStep),dMaxScaledMu(iDirection,iStep), &
                    iWorstMu(iDirection,iStep))
            end do

            if (VectorTwoNormFD(dMuAnalytic) <= 1D-12) then
                ! A direction that leaves this binary interaction unchanged has
                ! an identically zero derivative. Accuracy is meaningful here,
                ! but an observed order formed from zero/roundoff errors is not.
                lPass = lPass .AND. ALL(dNormAbsMu(iDirection,:) <= 1D-10)
            else
                lPass = lPass .AND. (MINVAL(dErrMu(iDirection,:)) <= 1D-9)
            end if
        end do

        dMolFraction(iFirstSpecies:iLastSpecies) = dMolFractionSave
        dPartialExcessGibbs(iFirstSpecies:iLastSpecies) = dPartialExcessSave

        dWorstMuBest = 0D0
        do iDirection = 1, nDirections
            dWorstMuBest = DMAX1(dWorstMuBest,MINVAL(dErrMu(iDirection,:)))
        end do

        if (lReport) call PrintReport

        deallocate(dDirection,dMoles0,dMolesM1,dMolesM2,dMolesP1,dMolesP2,dMuAnalytic, &
            dMuMinus,dMuPlus,dMolFractionSave,dPartialExcessSave,dX,dEnergyStep,dErrForwardEx, &
            dErrBackwardEx,dErr3Ex,dErr3Ideal,dErr5Ex,dErr5Ideal,dAbsForwardEx,dAbsBackwardEx, &
            dAbs3Ex,dAbs3Ideal,dAbs5Ex,dAbs5Ideal,dOrderForwardEx,dOrderBackwardEx,dOrder3Ex, &
            dOrder3Ideal,dOrder5Ex,dOrder5Ideal,lOrderForwardEx,lOrderBackwardEx,lOrder3Ex, &
            lOrder3Ideal,lOrder5Ex,lOrder5Ideal,dErrMu,dHessian, &
            dMuStep,dNormAbsMu,dMaxAbsMu,dMaxScaledMu,dOrderMu,lOrderMu,iWorstMu,dMuFD)
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

    !=========================================================================================================
    ! SECTION 5: NUMERICAL CONTROLS AND REPORTING
    !=========================================================================================================

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
        CompScaledError = DABS(dApproximate-dExact) / &
            DMAX1(1D0,DABS(dApproximate),DABS(dExact))
    end function CompScaledError

    subroutine PrintReport
        integer :: iDirLocal, iStepLocal

        write(*,'(A)') 'RKMP Hessian verification'
        write(*,'(A)') 'native scope: converged TestThermo30 plain-RKMP topology, parameters, and partial molars'
        write(*,'(A)') 'order fixture: parsed TestThermo30 binary parameter temporarily evaluated at exponent four'
        write(*,'(A)') 'excluded: RKMPM magnetism, ternary/Muggiano curvature, GEM response mapping, solver behavior'
        write(*,'(A,I0)') 'storage bits      = ', nStorageBits
        write(*,'(A,I0)') 'decimal precision = ', nDecimalPrecision
        write(*,'(A,I0)') 'binary digits     = ', nBinaryDigits
        write(*,'(A,ES14.6)') 'machine epsilon   = ', dEpsilonMachine
        write(*,'(A,ES14.6)') 'symmetry residual = ', dSymmetryResidual
        write(*,'(A,ES14.6)') 'radial residual   = ', dRadialResidual

        write(*,'(/,A,I0,A,I0)') 'controlled energy direction: species ', iEnergyA, ' minus species ', iEnergyB
        write(*,'(A)') 'h            forward abs  forward scaled order    backward abs backward scaled order'
        do iStepLocal = 1, nEnergyFDSteps
            write(*,'(ES12.4,2(2X,ES12.4,2X,ES12.4,1X,A8))') dEnergyStep(iStepLocal), &
                dAbsForwardEx(iStepLocal),dErrForwardEx(iStepLocal), &
                TRIM(OrderLabel(dOrderForwardEx(iStepLocal),lOrderForwardEx(iStepLocal))), &
                dAbsBackwardEx(iStepLocal),dErrBackwardEx(iStepLocal), &
                TRIM(OrderLabel(dOrderBackwardEx(iStepLocal),lOrderBackwardEx(iStepLocal)))
        end do
        write(*,'(A,L1,A,L1)') 'controlled RKMP one-sided order/accuracy: forward=', &
            tRKMPForward%lPassed,' backward=',tRKMPBackward%lPassed

        write(*,'(/,A)') 'controlled centered-difference comparison'
        write(*,'(A)') 'h            RKMP3 abs    RKMP3 scaled order    RKMP5 abs    RKMP5 scaled order'
        do iStepLocal = 1, nEnergyFDSteps
            write(*,'(ES12.4,2(2X,ES12.4,2X,ES12.4,1X,A8))') dEnergyStep(iStepLocal), &
                dAbs3Ex(iStepLocal),dErr3Ex(iStepLocal),TRIM(OrderLabel(dOrder3Ex(iStepLocal), &
                lOrder3Ex(iStepLocal))),dAbs5Ex(iStepLocal),dErr5Ex(iStepLocal), &
                TRIM(OrderLabel(dOrder5Ex(iStepLocal),lOrder5Ex(iStepLocal)))
        end do
        write(*,'(A,L1,A,L1)') 'controlled RKMP centered order/accuracy: 3pt=',tRKMP3%lPassed, &
            ' 5pt=',tRKMP5%lPassed

        write(*,'(/,A)') 'ideal-mixing numerical control'
        write(*,'(A)') 'h            ideal3 abs   ideal3 scaled order    ideal5 abs   ideal5 scaled order'
        do iStepLocal = 1, nEnergyFDSteps
            write(*,'(ES12.4,2(2X,ES12.4,2X,ES12.4,1X,A8))') dEnergyStep(iStepLocal), &
                dAbs3Ideal(iStepLocal),dErr3Ideal(iStepLocal),TRIM(OrderLabel(dOrder3Ideal(iStepLocal), &
                lOrder3Ideal(iStepLocal))),dAbs5Ideal(iStepLocal),dErr5Ideal(iStepLocal), &
                TRIM(OrderLabel(dOrder5Ideal(iStepLocal),lOrder5Ideal(iStepLocal)))
        end do
        write(*,'(A,L1,A,L1)') 'ideal order/accuracy: 3pt=',tIdeal3%lPassed,' 5pt=',tIdeal5%lPassed

        write(*,'(/,A)') 'native TestThermo30 production partial-molar RKMP comparison'
        write(*,'(A)') 'observed order: N/A when the native low-order polynomial begins in the roundoff regime'
        do iDirLocal = 1, nDirections
            write(*,'(/,A,I0,A,I0)') 'direction: species ', iDirLocal, ' minus species ', iReference
            write(*,'(A)') 'h            norm abs      norm scaled   max abs       max scaled    worst  order'
            do iStepLocal = 1, nMuFDSteps
                write(*,'(ES12.4,4(2X,ES12.4),2X,I5,2X,A8)') dMuStep(iDirLocal,iStepLocal), &
                    dNormAbsMu(iDirLocal,iStepLocal),dErrMu(iDirLocal,iStepLocal), &
                    dMaxAbsMu(iDirLocal,iStepLocal),dMaxScaledMu(iDirLocal,iStepLocal), &
                    iWorstMu(iDirLocal,iStepLocal),TRIM(OrderLabel(dOrderMu(iDirLocal,iStepLocal), &
                    lOrderMu(iDirLocal,iStepLocal)))
            end do
        end do
        write(*,'(/,A,ES14.6)') 'worst best production-mu error = ', dWorstMuBest
    end subroutine PrintReport

    function OrderLabel(dOrder,lAvailable) result(cLabel)
        real(8), intent(in) :: dOrder
        logical, intent(in) :: lAvailable
        character(len=8) :: cLabel

        if (lAvailable) then
            write(cLabel,'(F8.4)') dOrder
        else
            cLabel = 'N/A'
        end if
    end function OrderLabel

end program TestRKMPHessianVerification
