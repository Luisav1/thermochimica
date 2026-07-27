!-------------------------------------------------------------------------------------------------------------
!> \file    TestRKMPExactHessianRegression.F90
!> \brief   Curvature-enabled regression and trust-sensitivity coverage for the RKMP prototype.
!>
!> \details This test enables the experimental path through its public control
!!          routine.  It checks the targeted RKMP-sensitive regression systems,
!!          including TestThermo90's transient RKMP activity, verifies the
!!          configured alpha cap, and repeats TestThermo30 with modestly tighter
!!          and looser numerical trust gates.  The complete exact-on suite,
!!          rather than this targeted executable, supplies regression evidence
!!          for calculations without active RKMP phases. Equilibrium outputs are
!!          never relaxed when the globalization settings change.
!-------------------------------------------------------------------------------------------------------------

program TestRKMPExactHessianRegression

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    interface
        subroutine SetRKMPHessianControls(lEnable, dAlphaMax, lDebug, iInfo)
            logical, intent(in)  :: lEnable, lDebug
            real(8), intent(in)  :: dAlphaMax
            integer, intent(out) :: iInfo
        end subroutine SetRKMPHessianControls

        subroutine SetRKMPHessianTrustThresholds(dUpdateRatioCap, dDirectionCosineMin, &
                                                 dDirectionDifferenceCap, dProgressAllowance, iInfo)
            real(8), intent(in)  :: dUpdateRatioCap, dDirectionCosineMin
            real(8), intent(in)  :: dDirectionDifferenceCap, dProgressAllowance
            integer, intent(out) :: iInfo
        end subroutine SetRKMPHessianTrustThresholds
    end interface

    integer :: iInfo
    logical :: lPass, lReport, lHaveCalculation
    character(len=32) :: cArgument

    lPass = .TRUE.
    lReport = .FALSE.
    lHaveCalculation = .FALSE.
    if (COMMAND_ARGUMENT_COUNT() > 0) then
        call GET_COMMAND_ARGUMENT(1, cArgument)
        lReport = TRIM(cArgument) == '--report'
    end if
    call SetRKMPHessianReportOutput(lReport)

    call SetRKMPHessianControls(.TRUE., 1D0, .FALSE., iInfo)
    lPass = lPass .AND. (iInfo == 0)
    call RunTest30('alpha_max=1', 1D0, lReport, lPass)
    call RunTest33(lReport, lPass)
    call RunTest90(lReport, lPass)

    call SetRKMPHessianControls(.TRUE., 1D-1, .FALSE., iInfo)
    lPass = lPass .AND. (iInfo == 0)
    call RunTest30('alpha_max=0.1', 1D-1, lReport, lPass)

    call SetRKMPHessianControls(.TRUE., 1D0, .FALSE., iInfo)
    call SetRKMPHessianTrustThresholds(1.15D0, 0.93D0, 0.40D0, 1.03D0, iInfo)
    lPass = lPass .AND. (iInfo == 0)
    call RunTest30('tighter trust', 1D0, lReport, lPass)

    call SetRKMPHessianTrustThresholds(1.35D0, 0.87D0, 0.65D0, 1.08D0, iInfo)
    lPass = lPass .AND. (iInfo == 0)
    call RunTest30('looser trust', 1D0, lReport, lPass)

    call ResetRKMPHessianTrustThresholds
    call ResetRKMPHessianControls
    call SetRKMPHessianReportOutput(.FALSE.)

    if (lPass) then
        write(*,'(A)') 'TestRKMPExactHessianRegression: PASS'
        call EXIT(0)
    else
        write(*,'(A)') 'TestRKMPExactHessianRegression: FAIL <---'
        call EXIT(1)
    end if

contains

    subroutine RunTest30(cLabel, dAlphaCap, lDetailed, lAllPass)

        character(len=*), intent(in) :: cLabel
        real(8), intent(in) :: dAlphaCap
        logical, intent(in) :: lDetailed
        logical, intent(inout) :: lAllPass
        logical :: lCasePass
        real(8) :: dMassResidual

        call PrepareFreshCalculation
        cThermoFileName = DATA_DIRECTORY // 'WAuArO-1.dat'
        dPressure = 1D0
        dTemperature = 1455D0
        dElementMass(74) = 1.95D0
        dElementMass(79) = 1D0
        dElementMass(18) = 2D0
        dElementMass(8) = 10D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        dMassResidual = ComputeMaximumElementResidual()
        lCasePass = (INFOThermo == 0) .AND. &
            (DABS(dGibbsEnergySys - (-4.620D5))/4.620D5 < 1D-3) .AND. &
            lRKMPHessianWasActive .AND. (nRKMPHessianApplyCount > 0) .AND. &
            (nRKMPHessianRejectLocalResponse == 0) .AND. &
            (iRKMPHessianLastFailureReason == RKMP_MAP_SUCCESS) .AND. &
            (dRKMPHessianMaxSelectedAlpha <= dAlphaCap + 1D-14) .AND. &
            (dRKMPHessianMaxSelectedAlpha > 0D0) .AND. (dMassResidual < 1D-8)
        if (dAlphaCap >= 1D0) lCasePass = lCasePass .AND. (nRKMPHessianFullAlphaCount > 0)

        if (lDetailed) call ReportCase('TestThermo30 '//TRIM(cLabel), dMassResidual)
        lAllPass = lAllPass .AND. lCasePass

    end subroutine RunTest30


    subroutine RunTest33(lDetailed, lAllPass)

        logical, intent(in) :: lDetailed
        logical, intent(inout) :: lAllPass
        logical :: lCasePass
        real(8) :: dMassResidual

        call PrepareFreshCalculation
        cThermoFileName = DATA_DIRECTORY // 'WAuArNeO-2.dat'
        dPressure = 2D0
        dTemperature = 900D0
        dElementMass(74) = 20D0
        dElementMass(79) = 2D0
        dElementMass(18) = 7D0
        dElementMass(8) = 5D0
        dElementMass(10) = 1D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        dMassResidual = ComputeMaximumElementResidual()
        lCasePass = (INFOThermo == 0) .AND. &
            (DABS(dMolFraction(1)-0.75306881663786374D0)/0.75306881663786374D0 < 1D-3) .AND. &
            (DABS(dMolFraction(9)-3.0917444033201544D-2)/3.0917444033201544D-2 < 1D-3) .AND. &
            (DABS(dGibbsEnergySys-3.06480D6)/3.06480D6 < 1D-3) .AND. &
            lRKMPHessianWasActive .AND. (nRKMPHessianApplyCount > 0) .AND. &
            (nRKMPHessianRejectLocalResponse == 0) .AND. &
            (iRKMPHessianLastFailureReason == RKMP_MAP_SUCCESS) .AND. &
            (dRKMPHessianMaxSelectedAlpha <= 1D0 + 1D-14) .AND. &
            (nRKMPHessianFullAlphaCount > 0) .AND. (dMassResidual < 1D-8)

        if (lDetailed) call ReportCase('TestThermo33 alpha_max=1', dMassResidual)
        lAllPass = lAllPass .AND. lCasePass

    end subroutine RunTest33


    subroutine RunTest90(lDetailed, lAllPass)

        logical, intent(in) :: lDetailed
        logical, intent(inout) :: lAllPass
        logical :: lCasePass
        real(8) :: dMassResidual

        call PrepareFreshCalculation
        cThermoFileName = DATA_DIRECTORY // 'CsI-Pham.dat'
        cInputUnitTemperature = 'C'
        dTemperature = 400D0
        dPressure = 1D-5
        dElementMass(53) = 1D0
        dElementMass(55) = 1D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        dMassResidual = ComputeMaximumElementResidual()
        lCasePass = (INFOThermo == 0) .AND. &
            (DABS(dGibbsEnergySys-(-4.41869D5))/4.41869D5 < 1D-3) .AND. &
            lRKMPHessianWasActive .AND. (nRKMPHessianApplyCount > 0) .AND. &
            (nRKMPHessianFullAlphaCount > 0) .AND. (nSolnPhases == 0) .AND. (nConPhases == 1) .AND. &
            (nRKMPHessianRejectLocalResponse == 0) .AND. (dMassResidual < 1D-8)

        if (lDetailed) call ReportCase('TestThermo90 transient RKMP, final pure assemblage', dMassResidual)
        lAllPass = lAllPass .AND. lCasePass

    end subroutine RunTest90


    subroutine PrepareFreshCalculation

        if (lHaveCalculation) call ResetThermoAll
        lHaveCalculation = .TRUE.
        cInputUnitTemperature = 'K'
        cInputUnitPressure = 'atm'
        cInputUnitMass = 'moles'
        dElementMass = 0D0

    end subroutine PrepareFreshCalculation


    real(8) function ComputeMaximumElementResidual()

        integer :: iElement, iPhase, iSpecies, iSolution
        real(8) :: dCalculated

        ComputeMaximumElementResidual = 0D0
        do iElement = 1, nElements
            dCalculated = 0D0
            do iSolution = 1, nSolnPhases
                iPhase = -iAssemblage(nElements-iSolution+1)
                do iSpecies = nSpeciesPhase(iPhase-1)+1, nSpeciesPhase(iPhase)
                    dCalculated = dCalculated + dMolesSpecies(iSpecies) * &
                        dStoichSpecies(iSpecies,iElement) / DFLOAT(iParticlesPerMole(iSpecies))
                end do
            end do
            do iPhase = 1, nConPhases
                dCalculated = dCalculated + dMolesPhase(iPhase) * &
                    dStoichSpecies(iAssemblage(iPhase),iElement)
            end do
            ComputeMaximumElementResidual = DMAX1(ComputeMaximumElementResidual, &
                DABS(dCalculated-dMolesElement(iElement))/DMAX1(1D0,DABS(dMolesElement(iElement))))
        end do

    end function ComputeMaximumElementResidual


    subroutine ReportCase(cLabel, dMassResidual)

        character(len=*), intent(in) :: cLabel
        real(8), intent(in) :: dMassResidual
        integer :: iIteration, nZero, nPositive

        nZero = 0
        nPositive = 0
        do iIteration = 1, MIN(iterGlobal,iterGlobalMax)
            if (dRKMPHessianAcceptedAlphaHistory(iIteration) == 0D0) nZero = nZero + 1
            if (dRKMPHessianAcceptedAlphaHistory(iIteration) > 0D0) nPositive = nPositive + 1
        end do

        write(*,'(A,1X,A)') 'case:', TRIM(cLabel)
        write(*,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') '  iterations=', iterGlobal, &
            ' final_norm=', dGEMFunctionNorm, ' mass_residual=', dMassResidual, &
            ' max_selected_alpha=', dRKMPHessianMaxSelectedAlpha
        write(*,'(A,I0,A,I0,A,I0,A,I0)') '  positive_alpha_steps=', nPositive, &
            ' zero_alpha_steps=', nZero, ' full_alpha_steps=', nRKMPHessianFullAlphaCount, &
            ' local_failures=', nRKMPHessianRejectLocalResponse

    end subroutine ReportCase

end program TestRKMPExactHessianRegression
