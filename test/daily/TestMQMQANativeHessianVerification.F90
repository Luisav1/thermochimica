!-------------------------------------------------------------------------------------------------------------
!> \file    TestMQMQANativeHessianVerification.F90
!> \brief   Thermochimica-native verification of the disconnected plain-SUBG MQMQA Hessian.
!>
!> \details Reproduce the TestThermo56 equilibrium with ordinary production
!!          Thermochimica, decode the active nonmagnetic Liquid SUBG state, and
!!          compare ModuleMQMQAUnconstrained with established production
!!          partial-molar thermodynamics.
!!
!!          Verification map:
!!          1. Parse CuFeC-Kang.dat and converge the normal 1400 K calculation.
!!          2. Translate the active Liquid's filtered runtime topology,
!!             coordination numbers, zeta values, reference energies, and
!!             active G parameters through ModuleMQMQAProductionAdapter.
!!          3. Compare reference/configurational and excess scalar blocks using
!!             Euler sums of production partial molars.
!!          4. Compare the analytic generic gradient directly with production
!!             partial molars at the same imposed local composition.
!!          5. Transfer quadruplet moles along independent total-preserving
!!             directions. Finite differences of production partial molars
!!             must agree with the analytic Hessian-vector product.
!!          6. Check raw symmetry and the degree-one extensivity identity H*n=0.
!!
!!          This test does not modify GEMNewton, phase selection, SUBQ behavior,
!!          or production thermodynamic formulas. Pass --report to print the
!!          decoded-state and finite-difference evidence.
!-------------------------------------------------------------------------------------------------------------

program TestMQMQANativeHessianVerification

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleMQMQAUnconstrained
    USE ModuleMQMQAProductionAdapter

    implicit none

    interface
        subroutine CompExcessGibbsEnergySUBG(iSolnIndex)
            integer, intent(in) :: iSolnIndex
        end subroutine CompExcessGibbsEnergySUBG
    end interface

    integer, parameter :: nSteps = 9
    integer :: i, iDirection, iFirst, iInfo, iLast, iPhaseIndex, iReference
    integer :: iSlot, iStep, nDirections, nQuad
    logical :: lPass, lReport
    character(len=32) :: cArgument
    real(8) :: dG, dGExcess, dGIdeal, dGReference
    real(8) :: dGradientError, dHomogeneity, dMinMoles
    real(8) :: dProductionExcess, dProductionReferenceIdeal, dSymmetry
    real(8) :: dWorstBestMu
    real(8), allocatable :: dDirection(:), dErrMu(:,:), dGradient(:), dHessian(:,:)
    real(8), allocatable :: dMoles(:), dMuProduction(:), dSteps(:,:)
    type(MQMQAModelData) :: tModel
    type(MQMQAInteractionTerm), allocatable :: tInteraction(:)

    lPass = .TRUE.
    lReport = .FALSE.
    if (COMMAND_ARGUMENT_COUNT() > 0) then
        call GET_COMMAND_ARGUMENT(1,cArgument)
        lReport = TRIM(cArgument) == '--report'
    end if

    lPass = lPass .AND. (STORAGE_SIZE(1D0) == 64)
    lPass = lPass .AND. (PRECISION(1D0) >= 15)
    lPass = lPass .AND. (DIGITS(1D0) >= 53)

    !=========================================================================================================
    ! SECTION 1: CONVERGED PRODUCTION PLAIN-SUBG STATE
    !=========================================================================================================
    cInputUnitTemperature = 'K'
    cInputUnitPressure = 'atm'
    cInputUnitMass = 'moles'
    cThermoFileName = DATA_DIRECTORY // 'CuFeC-Kang.dat'
    dTemperature = 1400D0
    dPressure = 1D0
    dElementMass = 0D0
    dElementMass(6) = 1D0
    dElementMass(26) = 1D0
    dElementMass(29) = 1D0

    call ParseCSDataFile(cThermoFileName)
    if (INFOThermo == 0) call Thermochimica
    lPass = lPass .AND. (INFOThermo == 0)

    iPhaseIndex = 0
    iSlot = 0
    if (INFOThermo == 0) then
        do i = 1, nElements
            if (iAssemblage(i) >= 0) cycle
            if (cSolnPhaseName(-iAssemblage(i)) == 'Liquid' .AND. &
                cSolnPhaseType(-iAssemblage(i)) == 'SUBG') then
                iPhaseIndex = -iAssemblage(i)
                iSlot = i
                exit
            end if
        end do
    end if
    lPass = lPass .AND. (iPhaseIndex > 0) .AND. (iSlot > 0)

    if (iPhaseIndex > 0) then
        iFirst = nSpeciesPhase(iPhaseIndex-1)+1
        iLast = nSpeciesPhase(iPhaseIndex)
        nQuad = iLast-iFirst+1
        allocate(dMoles(nQuad),dGradient(nQuad),dHessian(nQuad,nQuad),dMuProduction(nQuad))
        dMoles = dMolesSpecies(iFirst:iLast)
        dMinMoles = MINVAL(dMoles)
        lPass = lPass .AND. ALL(dMoles > 0D0)

        !=====================================================================================================
        ! SECTION 2: PRODUCTION-TO-GENERIC DECODING AND LOCAL ENERGY
        !=====================================================================================================
        call DecodeProductionSUBGPhase(iPhaseIndex,tModel,tInteraction,iInfo)
        lPass = lPass .AND. (iInfo == 0)
        if (iInfo == 0) then
            call CompMQMQAGibbsEnergyUnconstrained(tModel,dMoles,1D0,tInteraction, &
                dG,dGReference,dGIdeal,dGExcess,iInfo)
            lPass = lPass .AND. (iInfo == 0)
            call CompMQMQAHessianUnconstrained(tModel,dMoles,1D0,tInteraction,dHessian,iInfo, &
                dGibbs=dG,dGradient=dGradient)
            lPass = lPass .AND. (iInfo == 0)
        end if

        if (iInfo == 0) then
            call EvaluateProductionVector(iPhaseIndex,dMoles,dMuProduction,dProductionReferenceIdeal, &
                dProductionExcess,iInfo)
            lPass = lPass .AND. (iInfo == 0)

            lPass = lPass .AND. (ScaledError((dGReference+dGIdeal)/SUM(dMoles), &
                dProductionReferenceIdeal) <= 1D-10)
            lPass = lPass .AND. (ScaledError(dGExcess/SUM(dMoles),dProductionExcess) <= 1D-10)

            dGradientError = VectorTwoNorm(dGradient-dMuProduction) / &
                DMAX1(1D0,VectorTwoNorm(dGradient),VectorTwoNorm(dMuProduction))
            dSymmetry = FrobeniusNorm(dHessian-TRANSPOSE(dHessian)) / &
                DMAX1(1D0,FrobeniusNorm(dHessian))
            dHomogeneity = VectorTwoNorm(MATMUL(dHessian,dMoles)) / &
                DMAX1(1D0,FrobeniusNorm(dHessian)*VectorTwoNorm(dMoles))
            lPass = lPass .AND. ALL(IEEE_IS_FINITE(dHessian)) .AND. ALL(IEEE_IS_FINITE(dGradient))
            lPass = lPass .AND. (dGradientError <= 1D-10)
            lPass = lPass .AND. (dSymmetry <= 1D-12)
            lPass = lPass .AND. (dHomogeneity <= 1D-10)

            !=================================================================================================
            ! SECTION 3: PRODUCTION PARTIAL-MOLAR FINITE DIFFERENCES
            !
            ! The largest quadruplet is the common reference direction. Moving
            ! equal moles from it into another quadruplet preserves total phase
            ! amount and isolates composition response.
            !=================================================================================================
            iReference = MAXLOC(dMoles,1)
            nDirections = nQuad-1
            allocate(dDirection(nQuad),dErrMu(nDirections,nSteps),dSteps(nDirections,nSteps))
            iDirection = 0
            do i = 1, nQuad
                if (i == iReference) cycle
                iDirection = iDirection+1
                dDirection = 0D0
                dDirection(i) = 1D0
                dDirection(iReference) = -1D0
                call VerifyProductionDirection(iPhaseIndex,dMoles,dDirection,dHessian, &
                    dSteps(iDirection,:),dErrMu(iDirection,:),lPass)
            end do

            dWorstBestMu = 0D0
            do iDirection = 1, nDirections
                dWorstBestMu = DMAX1(dWorstBestMu,MINVAL(dErrMu(iDirection,:)))
            end do
            lPass = lPass .AND. (dWorstBestMu <= 1D-8)

            if (lReport) then
                write(*,'(A)') 'MQ-2B Thermochimica-native plain-SUBG Hessian verification'
                write(*,'(A,A)') 'phase = ',TRIM(cSolnPhaseName(iPhaseIndex))
                write(*,'(A,I0)') 'quadruplet count = ',nQuad
                write(*,'(A,I0)') 'active decoded G-family terms = ',SIZE(tInteraction)
                write(*,'(A,ES14.6)') 'phase amount = ',dMolesPhase(iSlot)
                write(*,'(A,ES14.6)') 'minimum quadruplet moles = ',dMinMoles
                write(*,'(A,ES14.6)') 'reference+configurational scalar error = ', &
                    ScaledError((dGReference+dGIdeal)/SUM(dMoles),dProductionReferenceIdeal)
                write(*,'(A,ES14.6)') 'G-family excess scalar error = ', &
                    ScaledError(dGExcess/SUM(dMoles),dProductionExcess)
                write(*,'(A,ES14.6)') 'direct gradient error = ',dGradientError
                write(*,'(A,ES14.6)') 'symmetry residual = ',dSymmetry
                write(*,'(A,ES14.6)') 'homogeneity residual = ',dHomogeneity
                write(*,'(A)') 'direction h                 production-mu error'
                do iDirection = 1, nDirections
                    do iStep = 1, nSteps
                        write(*,'(I5,2ES22.12)') iDirection,dSteps(iDirection,iStep), &
                            dErrMu(iDirection,iStep)
                    end do
                end do
                write(*,'(A,ES14.6)') 'worst best production-mu error = ',dWorstBestMu
            end if

            deallocate(dDirection,dErrMu,dSteps)
        end if
        deallocate(dMoles,dGradient,dHessian,dMuProduction)
    end if

    call ResetThermoAll
    if (lPass) then
        print *, 'TestMQMQANativeHessianVerification: PASS'
        call EXIT(0)
    else
        print *, 'TestMQMQANativeHessianVerification: FAIL <---'
        call EXIT(1)
    end if

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Compare H*v with finite differences of established production partial molars.
    !---------------------------------------------------------------------------------------------------------
    subroutine VerifyProductionDirection(iPhaseLocal,dMolesLocal,dDirectionLocal,dHessianLocal, &
        dStepValues,dErrors,lAllPass)

        integer, intent(in) :: iPhaseLocal
        real(8), intent(in) :: dMolesLocal(:), dDirectionLocal(:), dHessianLocal(:,:)
        real(8), intent(out) :: dStepValues(:), dErrors(:)
        logical, intent(inout) :: lAllPass

        integer :: iInfoLocal, iStepLocal
        real(8) :: dDummyExcess, dDummyReferenceIdeal, dH, dScale
        real(8), allocatable :: dMinus(:), dMuMinus(:), dMuPlus(:), dPlus(:), dPrediction(:), dRatio(:)

        allocate(dMinus(SIZE(dMolesLocal)),dMuMinus(SIZE(dMolesLocal)), &
            dMuPlus(SIZE(dMolesLocal)),dPlus(SIZE(dMolesLocal)), &
            dPrediction(SIZE(dMolesLocal)),dRatio(SIZE(dMolesLocal)))
        where (DABS(dDirectionLocal) > 0D0)
            dRatio = dMolesLocal/DABS(dDirectionLocal)
        elsewhere
            dRatio = HUGE(1D0)
        end where
        dScale = MINVAL(dRatio)
        dPrediction = MATMUL(dHessianLocal,dDirectionLocal)

        ! Begin safely inside the positive-mole domain, then continue far
        ! enough to observe both truncation-error reduction and the eventual
        ! double-precision roundoff floor.  Some converged quadruplet amounts
        ! are very small, so six reductions do not reach that floor.
        do iStepLocal = 1, SIZE(dStepValues)
            dH = 0.1D0*dScale*3D0**(-(iStepLocal-1))
            dStepValues(iStepLocal) = dH
            dMinus = dMolesLocal-dH*dDirectionLocal
            dPlus = dMolesLocal+dH*dDirectionLocal
            call EvaluateProductionVector(iPhaseLocal,dMinus,dMuMinus,dDummyReferenceIdeal, &
                dDummyExcess,iInfoLocal)
            lAllPass = lAllPass .AND. (iInfoLocal == 0)
            call EvaluateProductionVector(iPhaseLocal,dPlus,dMuPlus,dDummyReferenceIdeal, &
                dDummyExcess,iInfoLocal)
            lAllPass = lAllPass .AND. (iInfoLocal == 0)
            dErrors(iStepLocal) = VectorTwoNorm((dMuPlus-dMuMinus)/(2D0*dH)-dPrediction) / &
                DMAX1(1D0,VectorTwoNorm(dPrediction))
        end do

        deallocate(dMinus,dMuMinus,dMuPlus,dPlus,dPrediction,dRatio)

    end subroutine VerifyProductionDirection


    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate production SUBG partial molars at a temporary local mole state.
    !>
    !> \details The routine saves and restores every global vector modified by
    !!          CompExcessGibbsEnergySUBG. This makes production thermodynamics an
    !!          independent oracle without altering the converged calculation.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateProductionVector(iPhaseLocal,dMolesLocal,dMu,dReferenceIdeal,dExcess,iInfoLocal)

        integer, intent(in) :: iPhaseLocal
        real(8), intent(in) :: dMolesLocal(:)
        real(8), intent(out) :: dMu(:), dReferenceIdeal, dExcess
        integer, intent(out) :: iInfoLocal

        integer :: iFirstLocal, iInfoSave, iLastLocal
        real(8), allocatable :: dChemicalSave(:), dFractionSave(:), dPartialSave(:), dX(:)

        iInfoLocal = 0
        iFirstLocal = nSpeciesPhase(iPhaseLocal-1)+1
        iLastLocal = nSpeciesPhase(iPhaseLocal)
        if ((SIZE(dMolesLocal) /= iLastLocal-iFirstLocal+1) .OR. &
            (SIZE(dMu) /= SIZE(dMolesLocal)) .OR. (SUM(dMolesLocal) <= 0D0) .OR. &
            ANY(dMolesLocal <= 0D0)) then
            iInfoLocal = 1
            return
        end if

        allocate(dChemicalSave(SIZE(dChemicalPotential)),dFractionSave(SIZE(dMolFraction)), &
            dPartialSave(SIZE(dPartialExcessGibbs)),dX(SIZE(dMolesLocal)))
        dChemicalSave = dChemicalPotential
        dFractionSave = dMolFraction
        dPartialSave = dPartialExcessGibbs
        iInfoSave = INFOThermo

        dMolFraction(iFirstLocal:iLastLocal) = dMolesLocal/SUM(dMolesLocal)
        call CompExcessGibbsEnergySUBG(iPhaseLocal)
        dX = dMolFraction(iFirstLocal:iLastLocal)
        dReferenceIdeal = DOT_PRODUCT(dX,dChemicalPotential(iFirstLocal:iLastLocal))
        dExcess = DOT_PRODUCT(dX,dPartialExcessGibbs(iFirstLocal:iLastLocal))
        dMu = dChemicalPotential(iFirstLocal:iLastLocal)+dPartialExcessGibbs(iFirstLocal:iLastLocal)
        if (INFOThermo /= iInfoSave) iInfoLocal = 2

        dChemicalPotential = dChemicalSave
        dMolFraction = dFractionSave
        dPartialExcessGibbs = dPartialSave
        INFOThermo = iInfoSave
        deallocate(dChemicalSave,dFractionSave,dPartialSave,dX)

    end subroutine EvaluateProductionVector


    real(8) function ScaledError(dA,dB)

        real(8), intent(in) :: dA, dB
        ScaledError = DABS(dA-dB)/DMAX1(1D0,DABS(dA),DABS(dB))

    end function ScaledError


    real(8) function VectorTwoNorm(dVector)

        real(8), intent(in) :: dVector(:)
        VectorTwoNorm = SQRT(SUM(dVector*dVector))

    end function VectorTwoNorm


    real(8) function FrobeniusNorm(dMatrix)

        real(8), intent(in) :: dMatrix(:,:)
        FrobeniusNorm = SQRT(SUM(dMatrix*dMatrix))

    end function FrobeniusNorm

end program TestMQMQANativeHessianVerification
