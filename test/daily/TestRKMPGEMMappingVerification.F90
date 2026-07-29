!-------------------------------------------------------------------------------------------------------------
!> \file    TestRKMPGEMMappingVerification.F90
!> \brief   Independent local-response and reduced GEM-mapping verification for plain RKMP.
!>
!> \details At the converged TestThermo30 state, this test first reconstructs
!!          the ideal and ideal-plus-excess normalized responses as a structural
!!          unit check of the production mapper.  It then supplies independent
!!          evidence by perturbing local thermodynamic driving forces and
!!          re-solving the nonlinear normalized composition problem using the
!!          established production RKMP partial-molar routine and a numerical
!!          Jacobian.  The resulting dx/dGamma, reduced matrix response, and
!!          right-hand-side response are compared with the analytic mapper.
!!          Singular, non-finite, alpha-zero, and no-active-phase controls verify
!!          safe failure and no-op behavior.
!-------------------------------------------------------------------------------------------------------------

program TestRKMPGEMMappingVerification

    USE, INTRINSIC :: IEEE_ARITHMETIC
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleRKMPResponseMapping
    USE ModuleFiniteDifferenceVerification, ONLY: AssessFDSweep, FDSweepAssessment, &
        FD_ORDER_SECOND_MIN, FD_ORDER_SECOND_MAX

    implicit none

    interface
        subroutine CompExcessGibbsEnergyRKMP_unconstrained(iSolnIndex,dHess)
            integer, intent(in) :: iSolnIndex
            real(8), intent(out), dimension(:,:) :: dHess
        end subroutine CompExcessGibbsEnergyRKMP_unconstrained

        subroutine CompExcessGibbsEnergyRKMP(iSolnIndex)
            integer :: iSolnIndex
        end subroutine CompExcessGibbsEnergyRKMP

        subroutine MapRKMPHessianToGEMVariables(A,B,nVar,dAlphaInput,lUpdateMetrics, &
                                                lCorrectionOK,dTrialMaxRatio)
            integer :: nVar
            real(8), dimension(nVar,nVar) :: A
            real(8), dimension(nVar) :: B
            real(8) :: dAlphaInput, dTrialMaxRatio
            logical :: lUpdateMetrics, lCorrectionOK
        end subroutine MapRKMPHessianToGEMVariables
    end interface

    integer, parameter :: nResponseSteps = 14
    integer :: i, j, k, p, iInfo, iFirst, iLast, iStep, iBestStep
    integer :: nSpeciesLocal, nRKMPActive, nSolnSave
    real(8) :: dN, dScaleA, dScaleB, dIdealError, dMapAError, dMapBError
    real(8) :: dConstraintResidual, dStationarityResidual, dFDResponseError, dTrialRatio
    real(8) :: dFullResponseFDError, dReducedAFDError, dReducedBFDError
    real(8) :: dDeltaBMagnitude, dForceScale, dResponseStep
    logical :: lPass, lCorrectionOK, lPlusOK, lMinusOK
    real(8), allocatable :: dX(:), dMu(:), dMuSave(:), dDirection(:)
    real(8), allocatable :: dPartialSave(:), dExpectedB(:), dFDB(:), dAep(:), dPhiBase(:)
    real(8), allocatable :: dGammaDirection(:), dForceGamma(:), dForceMu(:)
    real(8), allocatable :: dXPlus(:), dXMinus(:), dXIdealPlus(:), dXIdealMinus(:)
    real(8), allocatable :: dAnalyticDX(:), dFDDX(:), dResponseSteps(:), dResponseErrors(:)
    real(8), allocatable :: dResponseOrders(:)
    logical, allocatable :: lResponseOrderAvailable(:)
    type(FDSweepAssessment) :: tResponseSweep
    real(8), allocatable :: dC(:,:), dH(:,:), dHx(:,:), dHideal(:,:), dResponse(:,:), dIdealResponse(:,:)
    real(8), allocatable :: dMuRHS(:,:), dMuResponse(:,:), dIdealMuResponse(:,:)
    real(8), allocatable :: dDirectIdeal(:,:), dReconIdeal(:,:), dExpectedA(:,:), dA(:,:), dAOriginal(:,:)
    real(8), allocatable :: dFDResponse(:,:), dFDIdealResponse(:,:), dFDDeltaA(:,:)
    real(8), allocatable :: dB(:), dBOriginal(:)
    real(8) :: dSingularH(3,3), dSingularForce(3,1), dSingularResponse(3,1)

    lPass = .TRUE.
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

    nRKMPActive = 0
    k = 0
    do i = 1, nSolnPhases
        p = -iAssemblage(nElements-i+1)
        if (p > 0) then
            if (cSolnPhaseType(p) == 'RKMP') then
                nRKMPActive = nRKMPActive + 1
                k = p
            end if
        end if
    end do
    lPass = lPass .AND. (nRKMPActive == 1)

    iFirst = nSpeciesPhase(k-1) + 1
    iLast = nSpeciesPhase(k)
    nSpeciesLocal = iLast - iFirst + 1
    allocate(dX(nSpeciesLocal), dMu(nSpeciesLocal), dMuSave(nSpeciesLocal), &
             dDirection(nSpeciesLocal), dPartialSave(nSpeciesLocal), dExpectedB(nElements), &
             dFDB(nElements), &
             dAep(nElements), dPhiBase(nSpeciesLocal), dGammaDirection(nElements), &
             dForceGamma(nSpeciesLocal), dForceMu(nSpeciesLocal), dXPlus(nSpeciesLocal), &
             dXMinus(nSpeciesLocal), dXIdealPlus(nSpeciesLocal), dXIdealMinus(nSpeciesLocal), &
             dAnalyticDX(nSpeciesLocal), dFDDX(nSpeciesLocal), dResponseSteps(nResponseSteps), &
             dResponseErrors(nResponseSteps), dResponseOrders(nResponseSteps), &
             lResponseOrderAvailable(nResponseSteps), &
             dC(nSpeciesLocal,nElements), dH(nSpeciesLocal,nSpeciesLocal), &
             dHx(nSpeciesLocal,nSpeciesLocal), dHideal(nSpeciesLocal,nSpeciesLocal), &
             dResponse(nSpeciesLocal,nElements), dIdealResponse(nSpeciesLocal,nElements), &
             dMuRHS(nSpeciesLocal,1), dMuResponse(nSpeciesLocal,1), &
             dIdealMuResponse(nSpeciesLocal,1), dDirectIdeal(nElements,nElements), &
             dReconIdeal(nElements,nElements), dExpectedA(nElements,nElements), &
             dFDResponse(nSpeciesLocal,nElements), dFDIdealResponse(nSpeciesLocal,nElements), &
             dFDDeltaA(nElements,nElements), &
             dA(nElements,nElements), dAOriginal(nElements,nElements), &
             dB(nElements), dBOriginal(nElements))

    dN = SUM(dMolesSpecies(iFirst:iLast))
    dX = dMolesSpecies(iFirst:iLast) / dN
    dMu = dChemicalPotential(iFirst:iLast)
    dMuRHS(:,1) = dMu
    do j = 1, nElements
        do i = 1, nSpeciesLocal
            p = iFirst + i - 1
            dC(i,j) = dStoichSpecies(p,j) / DFLOAT(iParticlesPerMole(p))
        end do
    end do

    call CompExcessGibbsEnergyRKMP_unconstrained(k, dH)
    dHx = dN * dH
    dHideal = 0D0
    do i = 1, nSpeciesLocal
        dHideal(i,i) = 1D0 / dX(i)
        dHx(i,i) = dHx(i,i) + dHideal(i,i)
    end do

    call SolveRKMPConstrainedResponse(nSpeciesLocal, nElements, dHx, dC, dResponse, iInfo)
    lPass = lPass .AND. (iInfo == 0)
    call SolveRKMPConstrainedResponse(nSpeciesLocal, nElements, dHideal, dC, dIdealResponse, iInfo)
    lPass = lPass .AND. (iInfo == 0)
    call SolveRKMPConstrainedResponse(nSpeciesLocal, 1, dHx, dMuRHS, dMuResponse, iInfo)
    lPass = lPass .AND. (iInfo == 0)
    call SolveRKMPConstrainedResponse(nSpeciesLocal, 1, dHideal, dMuRHS, dIdealMuResponse, iInfo)
    lPass = lPass .AND. (iInfo == 0)

    dConstraintResidual = DMAX1(MAXVAL(DABS(SUM(dResponse,DIM=1))), &
                                DABS(SUM(dMuResponse(:,1))))
    dStationarityResidual = 0D0
    do j = 1, nElements
        dDirection = MATMUL(dHx,dResponse(:,j)) - dC(:,j)
        dStationarityResidual = DMAX1(dStationarityResidual, &
            MAXVAL(DABS(dDirection-dDirection(1))))
    end do

    dDirectIdeal = MATMUL(TRANSPOSE(dC), dN * SPREAD(dX,2,nElements) * dC)
    dAep = MATMUL(TRANSPOSE(dC), dN*dX)
    dReconIdeal = MATMUL(RESHAPE(dAep,[nElements,1]),RESHAPE(dAep,[1,nElements]))/dN + &
        MATMUL(TRANSPOSE(dC),dN*dIdealResponse)
    dScaleA = DMAX1(1D0,MAXVAL(DABS(dDirectIdeal)))
    dIdealError = MAXVAL(DABS(dDirectIdeal-dReconIdeal))/dScaleA

    dExpectedA = MATMUL(TRANSPOSE(dC),dN*(dResponse-dIdealResponse))
    dExpectedB = MATMUL(TRANSPOSE(dC),dN*(dMuResponse(:,1)-dIdealMuResponse(:,1)))

    dA = 0D0
    dB = 0D0
    dAOriginal = dA
    dBOriginal = dB
    call MapRKMPHessianToGEMVariables(dA,dB,nElements,1D0,.FALSE.,lCorrectionOK,dTrialRatio)
    dExpectedA = 0.5D0*(dExpectedA+TRANSPOSE(dExpectedA))
    dScaleA = DMAX1(1D0,MAXVAL(DABS(dExpectedA)))
    dScaleB = DMAX1(1D0,MAXVAL(DABS(dExpectedB)))
    dMapAError = MAXVAL(DABS((dA-dAOriginal)-dExpectedA))/dScaleA
    dMapBError = MAXVAL(DABS((dB-dBOriginal)-dExpectedB))/dScaleB
    lPass = lPass .AND. lCorrectionOK .AND. (dMapAError < 1D-11) .AND. (dMapBError < 1D-11)

    ! Independent nonlinear-response oracle. The base potential is evaluated by
    ! the established production RKMP routine. Perturbed compositions are then
    ! found with a separate finite-difference Newton solve in nSpecies-1
    ! normalized coordinates; this path does not call the analytic response
    ! solver or use Hx while finding the perturbed states.
    dMuSave = dMolFraction(iFirst:iLast)
    dPartialSave = dPartialExcessGibbs(iFirst:iLast)
    call EvaluateLocalPotential(dX,dPhiBase,iInfo)
    lPass = lPass .AND. (iInfo == 0)

    ! A mixed element-potential direction checks the complete composition
    ! response rather than another Hx*v product.
    do j = 1, nElements
        dGammaDirection(j) = 1D0 / DFLOAT(j)
    end do
    dForceGamma = MATMUL(dC,dGammaDirection)
    dAnalyticDX = MATMUL(dResponse,dGammaDirection)
    do iStep = 1, nResponseSteps
        dResponseSteps(iStep) = 4D-2 * 0.5D0**(iStep-1)
        call SolveForcedLocalEquilibrium(dX,dPhiBase,dForceGamma,dResponseSteps(iStep),dXPlus,lPlusOK)
        call SolveForcedLocalEquilibrium(dX,dPhiBase,dForceGamma,-dResponseSteps(iStep),dXMinus,lMinusOK)
        lPass = lPass .AND. lPlusOK .AND. lMinusOK
        if (lPlusOK .AND. lMinusOK) then
            dFDDX = (dXPlus-dXMinus)/(2D0*dResponseSteps(iStep))
            dResponseErrors(iStep) = VectorScaledError(dFDDX,dAnalyticDX)
        else
            dResponseErrors(iStep) = HUGE(1D0)
        end if
    end do
    call AssessFDSweep(dResponseSteps,dResponseErrors,FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX, &
        1D-7,tResponseSweep,dResponseOrders,lResponseOrderAvailable)
    iBestStep = tResponseSweep%iBest
    dResponseStep = dResponseSteps(iBestStep)
    dFDResponseError = dResponseErrors(iBestStep)
    ! Require the independent nonlinear response to approach the analytic
    ! response with the expected second-order central-difference rate. This
    ! prevents one accidentally favourable perturbation from passing the mapper.
    lPass = lPass .AND. tResponseSweep%lPassed
    lPass = lPass .AND. ALL(dResponseErrors(2:nResponseSteps) < &
        dResponseErrors(1:nResponseSteps-1))

    ! Rebuild every column of dx/dGamma independently. Condensing those
    ! finite-difference responses through elemental stoichiometry gives a
    ! reduced delta-A check that does not duplicate the mapper formula's local
    ! response solve.
    do j = 1, nElements
        dForceGamma = dC(:,j)
        call SolveForcedLocalEquilibrium(dX,dPhiBase,dForceGamma,dResponseStep,dXPlus,lPlusOK)
        call SolveForcedLocalEquilibrium(dX,dPhiBase,dForceGamma,-dResponseStep,dXMinus,lMinusOK)
        lPass = lPass .AND. lPlusOK .AND. lMinusOK
        if (lPlusOK .AND. lMinusOK) then
            dFDResponse(:,j) = (dXPlus-dXMinus)/(2D0*dResponseStep)
        else
            dFDResponse(:,j) = 0D0
        end if
        call EvaluateIdealForcedState(dX,dForceGamma,dResponseStep,dXIdealPlus)
        call EvaluateIdealForcedState(dX,dForceGamma,-dResponseStep,dXIdealMinus)
        dFDIdealResponse(:,j) = (dXIdealPlus-dXIdealMinus)/(2D0*dResponseStep)
    end do
    dFullResponseFDError = MatrixScaledError(dFDResponse,dResponse)
    dFDDeltaA = MATMUL(TRANSPOSE(dC),dN*(dFDResponse-dFDIdealResponse))
    dReducedAFDError = MatrixScaledError(dFDDeltaA,dExpectedA)

    ! The mapper's B correction is the reduced response to its current
    ! species-level residual forcing. Remove the irrelevant common component,
    ! normalize it for the nonlinear perturbation, and restore the scale after
    ! differencing. This makes the finite-difference delta-B both independent
    ! and numerically visible.
    dForceMu = dMu - dMu(nSpeciesLocal)
    dForceScale = MAXVAL(DABS(dForceMu))
    lPass = lPass .AND. (dForceScale > 0D0)
    if (dForceScale > 0D0) dForceMu = dForceMu / dForceScale
    call SolveForcedLocalEquilibrium(dX,dPhiBase,dForceMu,dResponseStep,dXPlus,lPlusOK)
    call SolveForcedLocalEquilibrium(dX,dPhiBase,dForceMu,-dResponseStep,dXMinus,lMinusOK)
    call EvaluateIdealForcedState(dX,dForceMu,dResponseStep,dXIdealPlus)
    call EvaluateIdealForcedState(dX,dForceMu,-dResponseStep,dXIdealMinus)
    lPass = lPass .AND. lPlusOK .AND. lMinusOK
    if (lPlusOK .AND. lMinusOK) then
        dFDDX = dForceScale*((dXPlus-dXMinus)-(dXIdealPlus-dXIdealMinus))/(2D0*dResponseStep)
        dFDB = MATMUL(TRANSPOSE(dC),dN*dFDDX)
        dReducedBFDError = VectorScaledError(dFDB,dExpectedB)
    else
        dReducedBFDError = HUGE(1D0)
    end if
    dDeltaBMagnitude = MAXVAL(DABS(dExpectedB))
    dMolFraction(iFirst:iLast) = dMuSave
    dPartialExcessGibbs(iFirst:iLast) = dPartialSave

    ! Alpha zero and a calculation with no active solution phase must leave both reduced objects untouched.
    dA = 3D0
    dB = -2D0
    dAOriginal = dA
    dBOriginal = dB
    call MapRKMPHessianToGEMVariables(dA,dB,nElements,0D0,.FALSE.,lCorrectionOK,dTrialRatio)
    lPass = lPass .AND. lCorrectionOK .AND. ALL(dA == dAOriginal) .AND. ALL(dB == dBOriginal)
    nSolnSave = nSolnPhases
    nSolnPhases = 0
    call MapRKMPHessianToGEMVariables(dA,dB,nElements,1D0,.FALSE.,lCorrectionOK,dTrialRatio)
    lPass = lPass .AND. lCorrectionOK .AND. ALL(dA == dAOriginal) .AND. ALL(dB == dBOriginal)
    nSolnPhases = nSolnSave

    ! A zero-curvature three-species system leaves two unconstrained difference modes and is singular.
    dSingularH = 0D0
    dSingularForce = 1D0
    call SolveRKMPConstrainedResponse(3,1,dSingularH,dSingularForce,dSingularResponse,iInfo)
    lPass = lPass .AND. (iInfo /= 0)

    ! Non-finite local forcing must invalidate the complete mapper result and identify its phase.
    dMuSave = dChemicalPotential(iFirst:iLast)
    dChemicalPotential(iFirst) = IEEE_VALUE(0D0,IEEE_QUIET_NAN)
    dA = dAOriginal
    dB = dBOriginal
    iRKMPHessianLastFailurePhase = 0
    iRKMPHessianLastFailureReason = RKMP_MAP_SUCCESS
    call MapRKMPHessianToGEMVariables(dA,dB,nElements,1D0,.FALSE.,lCorrectionOK,dTrialRatio)
    lPass = lPass .AND. (.NOT. lCorrectionOK) .AND. (iRKMPHessianLastFailurePhase == k) .AND. &
        (iRKMPHessianLastFailureReason == RKMP_MAP_INVALID_CORRECTION) .AND. &
        ALL(dA == dAOriginal) .AND. ALL(dB == dBOriginal)
    dChemicalPotential(iFirst:iLast) = dMuSave

    lPass = lPass .AND. (dConstraintResidual < 1D-12) .AND. &
        (dStationarityResidual < 1D-10) .AND. (dIdealError < 1D-12) .AND. &
        (dFDResponseError < 1D-7) .AND. (dFullResponseFDError < 1D-7) .AND. &
        (dReducedAFDError < 1D-7) .AND. (dReducedBFDError < 1D-7) .AND. &
        (dDeltaBMagnitude > 1D-8)

    if (COMMAND_ARGUMENT_COUNT() > 0) then
        write(*,'(A,ES14.6)') 'ideal reconstruction scaled error = ', dIdealError
        write(*,'(A,ES14.6)') 'composition constraint residual = ', dConstraintResidual
        write(*,'(A,ES14.6)') 'local stationarity residual = ', dStationarityResidual
        write(*,'(A,ES14.6)') 'independent mixed dx/dGamma FD error = ', dFDResponseError
        write(*,'(A,ES14.6)') 'independent full dx/dGamma FD error = ', dFullResponseFDError
        write(*,'(A,ES14.6)') 'independent reduced delta A FD error = ', dReducedAFDError
        write(*,'(A,ES14.6)') 'independent reduced delta B FD error = ', dReducedBFDError
        write(*,'(A,ES14.6)') 'max abs mapped delta B = ', dDeltaBMagnitude
        write(*,'(A,ES14.6)') 'structural mapped delta A error = ', dMapAError
        write(*,'(A,ES14.6)') 'structural mapped delta B error = ', dMapBError
        write(*,'(A,L1,A,F7.4,A,F7.4)') 'mixed-response order-aware pass = ', &
            tResponseSweep%lPassed, ' accepted p range = ',tResponseSweep%dObservedOrderMin, &
            ' to ',tResponseSweep%dObservedOrderMax
        write(*,'(A)') 'mixed-response step        scaled error    observed order'
        do iStep = 1, nResponseSteps
            write(*,'(ES20.10,2X,ES14.6,2X,A14)') dResponseSteps(iStep),dResponseErrors(iStep), &
                TRIM(ResponseOrderLabel(dResponseOrders(iStep),lResponseOrderAvailable(iStep)))
        end do
    end if

    call ResetThermo
    if (lPass) then
        write(*,'(A)') 'TestRKMPGEMMappingVerification: PASS'
        call EXIT(0)
    else
        write(*,'(A)') 'TestRKMPGEMMappingVerification: FAIL <---'
        call EXIT(1)
    end if

contains

    !> Evaluate the normalized local chemical-potential function used by the
    !! independent nonlinear oracle. The logarithm supplies ideal mixing and
    !! the established production routine supplies the RKMP excess part.
    subroutine EvaluateLocalPotential(dXTrial,dPotential,iEvalInfo)

        real(8), intent(in), dimension(:) :: dXTrial
        real(8), intent(out), dimension(:) :: dPotential
        integer, intent(out) :: iEvalInfo
        real(8), allocatable :: dXStored(:), dPartialStored(:)

        iEvalInfo = 0
        if (ANY(dXTrial <= 0D0) .OR. .NOT. ALL(IEEE_IS_FINITE(dXTrial))) then
            iEvalInfo = 1
            dPotential = 0D0
            return
        end if

        allocate(dXStored(nSpeciesLocal),dPartialStored(nSpeciesLocal))
        dXStored = dMolFraction(iFirst:iLast)
        dPartialStored = dPartialExcessGibbs(iFirst:iLast)
        dMolFraction(iFirst:iLast) = dXTrial
        dPartialExcessGibbs(iFirst:iLast) = 0D0
        call CompExcessGibbsEnergyRKMP(k)
        dPotential = LOG(dXTrial) + dPartialExcessGibbs(iFirst:iLast)
        if (.NOT. ALL(IEEE_IS_FINITE(dPotential))) iEvalInfo = 2
        dMolFraction(iFirst:iLast) = dXStored
        dPartialExcessGibbs(iFirst:iLast) = dPartialStored
        deallocate(dXStored,dPartialStored)

    end subroutine EvaluateLocalPotential

    !> Re-equilibrate composition after a species-level thermodynamic forcing.
    !! The final species is the dependent fraction, so every trial composition
    !! remains normalized. A numerical Jacobian of production partial molars
    !! keeps this oracle independent of the analytic Hessian and response solve.
    subroutine SolveForcedLocalEquilibrium(dXBase,dPotentialBase,dForcing,dAmplitude,dXSolved,lSolved)

        real(8), intent(in), dimension(:) :: dXBase, dPotentialBase, dForcing
        real(8), intent(in) :: dAmplitude
        real(8), intent(out), dimension(:) :: dXSolved
        logical, intent(out) :: lSolved

        integer, parameter :: nMaxIteration = 40
        integer :: iIteration, iColumn, iLapack, nIndependent
        integer, allocatable :: iPivot(:)
        real(8) :: dDifferenceStep, dStepFraction
        real(8), allocatable :: dResidual(:), dJacobian(:,:), dNewtonStep(:,:)
        real(8), allocatable :: dPotential(:), dPotentialPlus(:), dPotentialMinus(:)
        real(8), allocatable :: dXPlusLocal(:), dXMinusLocal(:), dXCandidate(:)

        nIndependent = SIZE(dXBase)-1
        allocate(dResidual(nIndependent),dJacobian(nIndependent,nIndependent), &
            dNewtonStep(nIndependent,1),iPivot(nIndependent),dPotential(SIZE(dXBase)), &
            dPotentialPlus(SIZE(dXBase)),dPotentialMinus(SIZE(dXBase)), &
            dXPlusLocal(SIZE(dXBase)),dXMinusLocal(SIZE(dXBase)),dXCandidate(SIZE(dXBase)))

        dXSolved = dXBase
        lSolved = .FALSE.
        do iIteration = 1, nMaxIteration
            dXSolved(SIZE(dXSolved)) = 1D0-SUM(dXSolved(1:nIndependent))
            if (MINVAL(dXSolved) <= 1D-14) exit
            call EvaluateLocalPotential(dXSolved,dPotential,iLapack)
            if (iLapack /= 0) exit
            dResidual = (dPotential(1:nIndependent)-dPotential(SIZE(dPotential))) - &
                (dPotentialBase(1:nIndependent)-dPotentialBase(SIZE(dPotentialBase))) - &
                dAmplitude*(dForcing(1:nIndependent)-dForcing(SIZE(dForcing)))
            if (MAXVAL(DABS(dResidual)) <= 2D-12) then
                lSolved = .TRUE.
                exit
            end if

            do iColumn = 1, nIndependent
                dDifferenceStep = DMIN1(1D-6,0.2D0*DMIN1(dXSolved(iColumn), &
                    dXSolved(SIZE(dXSolved))))
                if (dDifferenceStep <= 1D-14) exit
                dXPlusLocal = dXSolved
                dXMinusLocal = dXSolved
                dXPlusLocal(iColumn) = dXPlusLocal(iColumn) + dDifferenceStep
                dXPlusLocal(SIZE(dXSolved)) = dXPlusLocal(SIZE(dXSolved)) - dDifferenceStep
                dXMinusLocal(iColumn) = dXMinusLocal(iColumn) - dDifferenceStep
                dXMinusLocal(SIZE(dXSolved)) = dXMinusLocal(SIZE(dXSolved)) + dDifferenceStep
                call EvaluateLocalPotential(dXPlusLocal,dPotentialPlus,iLapack)
                if (iLapack /= 0) exit
                call EvaluateLocalPotential(dXMinusLocal,dPotentialMinus,iLapack)
                if (iLapack /= 0) exit
                dJacobian(:,iColumn) = ((dPotentialPlus(1:nIndependent)- &
                    dPotentialPlus(SIZE(dPotentialPlus))) - &
                    (dPotentialMinus(1:nIndependent)-dPotentialMinus(SIZE(dPotentialMinus)))) / &
                    (2D0*dDifferenceStep)
            end do
            if (iColumn <= nIndependent .OR. iLapack /= 0) exit

            dNewtonStep(:,1) = -dResidual
            call DGESV(nIndependent,1,dJacobian,nIndependent,iPivot,dNewtonStep,nIndependent,iLapack)
            if (iLapack /= 0 .OR. .NOT. ALL(IEEE_IS_FINITE(dNewtonStep))) exit

            dStepFraction = 1D0
            do
                dXCandidate = dXSolved
                dXCandidate(1:nIndependent) = dXSolved(1:nIndependent) + &
                    dStepFraction*dNewtonStep(:,1)
                dXCandidate(SIZE(dXCandidate)) = 1D0-SUM(dXCandidate(1:nIndependent))
                if (MINVAL(dXCandidate) > 1D-14) exit
                dStepFraction = 0.5D0*dStepFraction
                if (dStepFraction < 1D-10) exit
            end do
            if (dStepFraction < 1D-10) exit
            dXSolved = dXCandidate
        end do

        deallocate(dResidual,dJacobian,dNewtonStep,iPivot,dPotential,dPotentialPlus, &
            dPotentialMinus,dXPlusLocal,dXMinusLocal,dXCandidate)

    end subroutine SolveForcedLocalEquilibrium

    !> Closed-form ideal-mixing response used as an independent control for the
    !! ideal contribution already present in GEMNewton.
    subroutine EvaluateIdealForcedState(dXBase,dForcing,dAmplitude,dXForced)

        real(8), intent(in), dimension(:) :: dXBase, dForcing
        real(8), intent(in) :: dAmplitude
        real(8), intent(out), dimension(:) :: dXForced

        dXForced = dXBase*EXP(dAmplitude*dForcing)
        dXForced = dXForced/SUM(dXForced)

    end subroutine EvaluateIdealForcedState

    real(8) function VectorScaledError(dApproximate,dReference)

        real(8), intent(in), dimension(:) :: dApproximate, dReference

        VectorScaledError = SQRT(SUM((dApproximate-dReference)**2)) / &
            DMAX1(1D0,SQRT(SUM(dApproximate**2)),SQRT(SUM(dReference**2)))

    end function VectorScaledError

    real(8) function MatrixScaledError(dApproximate,dReference)

        real(8), intent(in), dimension(:,:) :: dApproximate, dReference

        MatrixScaledError = SQRT(SUM((dApproximate-dReference)**2)) / &
            DMAX1(1D0,SQRT(SUM(dApproximate**2)),SQRT(SUM(dReference**2)))

    end function MatrixScaledError

    function ResponseOrderLabel(dOrder,lAvailable) result(cLabel)

        real(8), intent(in) :: dOrder
        logical, intent(in) :: lAvailable
        character(len=14) :: cLabel

        if (lAvailable) then
            write(cLabel,'(F14.6)') dOrder
        else
            cLabel = 'N/A'
        end if

    end function ResponseOrderLabel

end program TestRKMPGEMMappingVerification
