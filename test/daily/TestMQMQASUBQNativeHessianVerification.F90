!-------------------------------------------------------------------------------------------------------------
!> \file    TestMQMQASUBQNativeHessianVerification.F90
!> \brief   Native first-derivative and Hessian verification for the assessed FeTiVO SUBQ G/Q case.
!>
!> \details MQ-2B verifies the standalone SUBQ derivatives against established
!!          Thermochimica partial molars without invoking equilibrium during the
!!          perturbation sweep. The production routine stores the derivative of
!!          reference plus configurational energy in dChemicalPotential and the
!!          derivative of G/Q excess energy in dPartialExcessGibbs. Their sum is
!!          therefore the complete unconstrained quadruplet-mole gradient used
!!          by CompExcessGibbsEnergy; no gauge transformation is required.
!!
!!          The real FeTiVO phase supplies assessed topology, reference energies,
!!          and G/Q parameters. Derivative checks use a fixed positive interior
!!          composition blended from the converged SlagBsoln state. This avoids
!!          letting trace quadruplets determine the perturbation scale while
!!          preserving the real production model and parameter decoding.
!!
!!          Scope remains narrow: FeTiVO has uniform zeta=2.4 and no B or R
!!          records. This test does not verify every SUBQ feature, constrained
!!          response, GEM mapping, phase switching, or solver integration.
!-------------------------------------------------------------------------------------------------------------

program TestMQMQASUBQNativeHessianVerification

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleMQMQAUnconstrained
    USE ModuleMQMQAProductionAdapter
    USE ModuleFiniteDifferenceVerification

    implicit none

    interface
        subroutine CompExcessGibbsEnergySUBG(iSolnIndex)
            integer, intent(in) :: iSolnIndex
        end subroutine CompExcessGibbsEnergySUBG
    end interface

    integer, parameter :: nSteps = 9
    integer :: i, iDirection, iFirst, iGradientWorst, iInfo, iLast, iPhaseIndex, iReference
    integer :: iSlot, iStep, nDirections, nG, nQ, nQuad
    logical :: lPass, lReport
    character(len=32) :: cArgument
    real(8) :: dG, dGExcess, dGIdeal, dGReference, dGradientAbsolute, dGradientError
    real(8) :: dGradientMaxAbsolute, dGradientMaxScaled
    real(8) :: dHomogeneity, dInteriorFraction, dMinConverged, dMinInterior
    real(8) :: dProductionExcess, dProductionReferenceIdeal, dSymmetry
    real(8) :: dWorstBestComponent, dWorstBestMu
    real(8), allocatable :: dConvergedMoles(:), dDirection(:), dErrMu(:,:), dGradient(:)
    real(8), allocatable :: dHessian(:,:), dMaxAbsMu(:,:), dMaxScaledMu(:,:), dMoles(:)
    real(8), allocatable :: dMuProduction(:), dNormAbsMu(:,:), dOrders(:,:), dSteps(:,:)
    integer, allocatable :: iWorstComponent(:,:)
    logical, allocatable :: lOrderAvailable(:,:)
    type(FDSweepAssessment) :: tSweep
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
    ! SECTION 1: REAL ASSESSED SUBQ MODEL AND CONVERGED STATE
    !=========================================================================================================
    cInputUnitTemperature = 'K'
    cInputUnitPressure = 'atm'
    cInputUnitMass = 'moles'
    cThermoFileName = DATA_DIRECTORY // 'FeTiVO.dat'
    dTemperature = 2000D0
    dPressure = 1D0
    dElementMass = 0D0
    dElementMass(8) = 2D0
    dElementMass(22) = 0.5D0
    dElementMass(23) = 0.5D0
    dElementMass(26) = 0.5D0

    call ParseCSDataFile(cThermoFileName)
    if (INFOThermo == 0) call Thermochimica
    lPass = lPass .AND. (INFOThermo == 0)

    iPhaseIndex = 0
    iSlot = 0
    if (INFOThermo == 0) then
        do i = 1, nElements
            if (iAssemblage(i) >= 0) cycle
            if ((cSolnPhaseName(-iAssemblage(i)) == 'SlagBsoln') .AND. &
                (cSolnPhaseType(-iAssemblage(i)) == 'SUBQ')) then
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
        allocate(dConvergedMoles(nQuad),dMoles(nQuad),dGradient(nQuad), &
            dHessian(nQuad,nQuad),dMuProduction(nQuad))
        dConvergedMoles = dMolesSpecies(iFirst:iLast)
        dMinConverged = MINVAL(dConvergedMoles)

        call DecodeProductionSUBQPhase(iPhaseIndex,tModel,tInteraction,iInfo)
        lPass = lPass .AND. (iInfo == 0)
        if (iInfo == 0) then
            nG = COUNT(tInteraction%iFamily == MQMQA_TERM_G)
            nQ = COUNT(tInteraction%iFamily == MQMQA_TERM_Q)
            lPass = lPass .AND. (nQuad == 15) .AND. (nG == 6) .AND. (nQ == 8)
            lPass = lPass .AND. (MAXVAL(DABS(tModel%dZeta-2.4D0)) <= 1D-12)
        end if

        ! Blend twenty percent uniform composition into the converged state.
        ! The resulting state remains tied to the assessed phase but gives every
        ! independent mole-transfer direction a useful positivity-safe interval.
        dInteriorFraction = 0.20D0
        dMoles = (1D0-dInteriorFraction)*dConvergedMoles + &
            dInteriorFraction*SUM(dConvergedMoles)/DBLE(nQuad)
        dMinInterior = MINVAL(dMoles)
        lPass = lPass .AND. ALL(dMoles > 0D0)
        lPass = lPass .AND. (ScaledError(SUM(dMoles),SUM(dConvergedMoles)) <= 1D-14)

        !=====================================================================================================
        ! SECTION 2: COORDINATE CONTRACT AND DIRECT FIRST-DERIVATIVE PARITY
        !
        ! Both implementations differentiate the same extensive dimensionless
        ! Gibbs energy with respect to each independent quadruplet mole amount.
        ! Production's two derivative arrays are added exactly once here, just
        ! as CompExcessGibbsEnergy does for a normal SUBQ calculation.
        !=====================================================================================================
        if (iInfo == 0) then
            call CompMQMQAGibbsEnergyUnconstrained(tModel,dMoles,1D0,tInteraction, &
                dG,dGReference,dGIdeal,dGExcess,iInfo)
            lPass = lPass .AND. (iInfo == 0)
            call CompMQMQAHessianUnconstrained(tModel,dMoles,1D0,tInteraction,dHessian,iInfo, &
                dGibbs=dG,dGradient=dGradient)
            lPass = lPass .AND. (iInfo == 0)
        end if

        if (iInfo == 0) then
            call EvaluateProductionVector(iPhaseIndex,dMoles,dMuProduction, &
                dProductionReferenceIdeal,dProductionExcess,iInfo)
            lPass = lPass .AND. (iInfo == 0)

            lPass = lPass .AND. (ScaledError((dGReference+dGIdeal)/SUM(dMoles), &
                dProductionReferenceIdeal) <= 1D-10)
            lPass = lPass .AND. (ScaledError(dGExcess/SUM(dMoles),dProductionExcess) <= 1D-10)
            call ComputeVectorErrorMetrics(dGradient,dMuProduction,dGradientAbsolute,dGradientError, &
                dGradientMaxAbsolute,dGradientMaxScaled,iGradientWorst)
            dSymmetry = FrobeniusNorm(dHessian-TRANSPOSE(dHessian)) / &
                DMAX1(1D0,FrobeniusNorm(dHessian))
            dHomogeneity = VectorTwoNorm(MATMUL(dHessian,dMoles)) / &
                DMAX1(1D0,FrobeniusNorm(dHessian)*VectorTwoNorm(dMoles))
            lPass = lPass .AND. ALL(IEEE_IS_FINITE(dGradient)) .AND. ALL(IEEE_IS_FINITE(dHessian))
            lPass = lPass .AND. (dGradientError <= 1D-10)
            lPass = lPass .AND. (dGradientMaxScaled <= 1D-10)
            lPass = lPass .AND. (dSymmetry <= 1D-12)
            lPass = lPass .AND. (dHomogeneity <= 1D-10)

            !=================================================================================================
            ! SECTION 3: NATIVE PRODUCTION-PARTIAL-MOLAR FINITE DIFFERENCES
            !
            ! Each direction transfers moles between two quadruplets. Total phase
            ! amount, temperature, pressure, topology, and parameters stay fixed.
            ! Production is reevaluated locally; Thermochimica is not re-solved.
            !=================================================================================================
            iReference = MAXLOC(dMoles,1)
            nDirections = nQuad-1
            lPass = lPass .AND. (nDirections == 14)
            allocate(dDirection(nQuad),dErrMu(nDirections,nSteps),dSteps(nDirections,nSteps), &
                dNormAbsMu(nDirections,nSteps),dMaxAbsMu(nDirections,nSteps), &
                dMaxScaledMu(nDirections,nSteps),dOrders(nDirections,nSteps), &
                lOrderAvailable(nDirections,nSteps),iWorstComponent(nDirections,nSteps))

            iDirection = 0
            dWorstBestComponent = 0D0
            do i = 1, nQuad
                if (i == iReference) cycle
                iDirection = iDirection+1
                dDirection = 0D0
                dDirection(i) = 1D0
                dDirection(iReference) = -1D0
                call VerifyProductionDirection(iPhaseIndex,dMoles,dDirection,dHessian, &
                    dSteps(iDirection,:),dNormAbsMu(iDirection,:),dErrMu(iDirection,:), &
                    dMaxAbsMu(iDirection,:),dMaxScaledMu(iDirection,:), &
                    iWorstComponent(iDirection,:),lPass)
                call AssessFDSweep(dSteps(iDirection,:),dErrMu(iDirection,:), &
                    FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX,1D-8,tSweep, &
                    dOrders(iDirection,:),lOrderAvailable(iDirection,:))
                lPass = lPass .AND. tSweep%lPassed
                if (tSweep%iBest > 0) then
                    dWorstBestComponent = DMAX1(dWorstBestComponent, &
                        dMaxScaledMu(iDirection,tSweep%iBest))
                end if
            end do

            dWorstBestMu = 0D0
            do iDirection = 1, nDirections
                dWorstBestMu = DMAX1(dWorstBestMu,MINVAL(dErrMu(iDirection,:)))
            end do
            lPass = lPass .AND. (dWorstBestMu <= 1D-8)
            lPass = lPass .AND. (dWorstBestComponent <= 1D-8)

            if (lReport) then
                write(*,'(A)') 'MQ-2B Thermochimica-native SUBQ Hessian verification'
                write(*,'(A)') 'native scope: assessed FeTiVO nonmagnetic SUBQ reference/configurational/G/Q'
                write(*,'(A)') 'native exclusions: nonuniform zeta, B, R, magnetism, response mapping, GEM integration'
                write(*,'(A,A)') 'phase = ',TRIM(cSolnPhaseName(iPhaseIndex))
                write(*,'(A,I0,A,I0,A,I0)') 'quadruplets = ',nQuad,', G terms = ',nG,', Q terms = ',nQ
                write(*,'(A,ES14.6)') 'phase amount = ',dMolesPhase(iSlot)
                write(*,'(A,ES14.6)') 'minimum converged quadruplet moles = ',dMinConverged
                write(*,'(A,F7.3)') 'uniform interior blend fraction = ',dInteriorFraction
                write(*,'(A,ES14.6)') 'minimum verification-state quadruplet moles = ',dMinInterior
                write(*,'(A,ES14.6)') 'reference+configurational scalar error = ', &
                    ScaledError((dGReference+dGIdeal)/SUM(dMoles),dProductionReferenceIdeal)
                write(*,'(A,ES14.6)') 'G/Q excess scalar error = ', &
                    ScaledError(dGExcess/SUM(dMoles),dProductionExcess)
                write(*,'(A,ES14.6)') 'direct unconstrained-gradient error = ',dGradientError
                write(*,'(A,ES14.6,A,I0)') 'direct gradient max component error = ', &
                    dGradientMaxScaled,', worst quadruplet = ',iGradientWorst
                write(*,'(A,ES14.6)') 'symmetry residual = ',dSymmetry
                write(*,'(A,ES14.6)') 'homogeneity residual = ',dHomogeneity
                write(*,'(A)') 'central derivative of production partial molars; expected order = 2'
                write(*,'(A)') 'dir  h                 norm abs            norm scaled         max abs             max scaled          worst  order'
                do iDirection = 1, nDirections
                    do iStep = 1, nSteps
                        write(*,'(I3,5ES20.10,I7,2X,A)') iDirection,dSteps(iDirection,iStep), &
                            dNormAbsMu(iDirection,iStep),dErrMu(iDirection,iStep), &
                            dMaxAbsMu(iDirection,iStep),dMaxScaledMu(iDirection,iStep), &
                            iWorstComponent(iDirection,iStep),TRIM(OrderLabel( &
                            dOrders(iDirection,iStep),lOrderAvailable(iDirection,iStep)))
                    end do
                    call AssessFDSweep(dSteps(iDirection,:),dErrMu(iDirection,:), &
                        FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX,1D-8,tSweep)
                    write(*,'(A,I0,A,ES14.6,A,F6.3,A,F6.3,A,L1)') 'direction ',iDirection, &
                        ' best scaled error = ',tSweep%dBestError, &
                        ' order range = ',tSweep%dObservedOrderMin,' to ',tSweep%dObservedOrderMax, &
                        ' small-h upturn = ',tSweep%lRoundoffUpturn
                end do
                write(*,'(A,ES14.6)') 'worst best production-mu error = ',dWorstBestMu
                write(*,'(A,ES14.6)') &
                    'worst componentwise scaled error at normwise-best steps = ',dWorstBestComponent
            end if

            deallocate(dDirection,dErrMu,dSteps,dNormAbsMu,dMaxAbsMu,dMaxScaledMu, &
                dOrders,lOrderAvailable,iWorstComponent)
        end if
        deallocate(dConvergedMoles,dMoles,dGradient,dHessian,dMuProduction)
    end if

    call ResetThermoAll
    if (lPass) then
        print *, 'TestMQMQASUBQNativeHessianVerification: PASS'
        call EXIT(0)
    else
        print *, 'TestMQMQASUBQNativeHessianVerification: FAIL <---'
        call EXIT(1)
    end if

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Compare analytic H*v with centered differences of production partial molars.
    !---------------------------------------------------------------------------------------------------------
    subroutine VerifyProductionDirection(iPhaseLocal,dMolesLocal,dDirectionLocal,dHessianLocal, &
        dStepValues,dNormAbsolute,dErrors,dMaxAbsolute,dMaxScaled,iWorst,lAllPass)

        integer, intent(in) :: iPhaseLocal
        real(8), intent(in) :: dMolesLocal(:), dDirectionLocal(:), dHessianLocal(:,:)
        real(8), intent(out) :: dStepValues(:), dNormAbsolute(:), dErrors(:)
        real(8), intent(out) :: dMaxAbsolute(:), dMaxScaled(:)
        integer, intent(out) :: iWorst(:)
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

        do iStepLocal = 1, SIZE(dStepValues)
            dH = 0.05D0*dScale*3D0**(-(iStepLocal-1))
            dStepValues(iStepLocal) = dH
            dMinus = dMolesLocal-dH*dDirectionLocal
            dPlus = dMolesLocal+dH*dDirectionLocal
            call EvaluateProductionVector(iPhaseLocal,dMinus,dMuMinus,dDummyReferenceIdeal, &
                dDummyExcess,iInfoLocal)
            lAllPass = lAllPass .AND. (iInfoLocal == 0)
            call EvaluateProductionVector(iPhaseLocal,dPlus,dMuPlus,dDummyReferenceIdeal, &
                dDummyExcess,iInfoLocal)
            lAllPass = lAllPass .AND. (iInfoLocal == 0)
            call ComputeVectorErrorMetrics((dMuPlus-dMuMinus)/(2D0*dH),dPrediction, &
                dNormAbsolute(iStepLocal),dErrors(iStepLocal),dMaxAbsolute(iStepLocal), &
                dMaxScaled(iStepLocal),iWorst(iStepLocal))
        end do

        deallocate(dMinus,dMuMinus,dMuPlus,dPlus,dPrediction,dRatio)

    end subroutine VerifyProductionDirection


    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate complete production SUBQ partial molars at one imposed local mole state.
    !>
    !> \details CompExcessGibbsEnergySUBG computes local thermodynamics directly
    !!          from normalized quadruplet moles. It does not run GEM or change
    !!          the phase assemblage. All modified global arrays are restored.
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


    character(len=16) function OrderLabel(dOrder,lAvailable)
        real(8), intent(in) :: dOrder
        logical, intent(in) :: lAvailable
        if (lAvailable) then
            write(OrderLabel,'(F10.4)') dOrder
        else
            OrderLabel = 'N/A'
        end if
    end function OrderLabel


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

end program TestMQMQASUBQNativeHessianVerification
