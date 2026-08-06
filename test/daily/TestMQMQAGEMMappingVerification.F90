!-------------------------------------------------------------------------------------------------------------
!> \file    TestMQMQAGEMMappingVerification.F90
!> \brief   Diagnostic-only verification of the plain-SUBG reduced GEM response correction.
!>
!> \details MQ-4A established the mathematical correction; MQ-4B calls the reusable production builder and
!!          compares its output against the retained independent construction and nonlinear finite-difference
!!          oracle. The test also applies the correction only to caller-owned copies of GEM arrays.
!!
!!          The test deliberately moves to a nearby positive off-equilibrium composition. This makes the residual
!!          correction non-vacuous and checks the `mu-1` convention used by GEMNewton. No live solver matrix,
!!          control, phase assemblage, or production thermodynamic formula is changed.
!-------------------------------------------------------------------------------------------------------------
program TestMQMQAGEMMappingVerification

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleMQMQAUnconstrained
    USE ModuleMQMQAProductionAdapter
    USE ModuleMQMQAResponseMapping
    USE ModuleConstrainedResponse
    USE ModuleFiniteDifferenceVerification
    USE ModuleGEMNewtonDiagnosticCapture

    implicit none

    interface
        subroutine CompExcessGibbsEnergySUBG(iPhase)
            integer :: iPhase
        end subroutine CompExcessGibbsEnergySUBG
        subroutine CompStoichSolnPhase(iPhase)
            integer :: iPhase
        end subroutine CompStoichSolnPhase
        subroutine GEMNewton(iInfo)
            integer :: iInfo
        end subroutine GEMNewton
        subroutine CompChemicalPotential(lCompEverything)
            logical :: lCompEverything
        end subroutine CompChemicalPotential
    end interface

    integer, parameter :: nSteps = 8
    character(len=32) :: cArgument
    integer :: e, f, i, iApplyStatus, iBestA, iBestAffine, iBestB, iBuilderStatus
    integer :: iColumn, iDonor, iElectronSave, iFailureStatus, iInfo, iNewtonInfo
    integer :: iFirst, iLast, iPhaseIndex, iReceiver, iSlot, iStep, nQuad
    logical :: lBuilderApplicable, lFailureApplicable, lMinus, lPass, lPlus, lReport
    real(8) :: dAComponent, dAError, dBComponent, dBError
    real(8) :: dBaselineAError, dBaselineBError, dBaselineColumnError, dBaselinePhaseError
    real(8) :: dAggregationError, dApplyError, dBuilderAError, dBuilderBError, dBuilderStateDifference
    real(8) :: dZeroApplyError
    real(8) :: dCaptureControlA(1,1), dCaptureControlB(1)
    real(8) :: dConstantForceError, dConstraintResidual, dDeltaBMagnitude, dFloorDifference
    real(8) :: dAffineComponent, dAffineError, dH, dHColumn
    real(8) :: dHMaximumA, dHMaximumAffine, dHMaximumB, dN, dOffEquilibriumResidual, dSymmetryResidual
    real(8) :: dWorstColumnError, dWorstColumnUncertainty
    real(8), allocatable :: dABaseline(:,:), dABaseCopy(:,:), dACandidate(:,:), dADirect(:,:)
    real(8), allocatable :: dAApplyExpected(:,:), dABuilder(:,:), dATrial(:,:)
    real(8), allocatable :: dAExtracted(:,:), dAOther(:,:), dAep(:), dAErrors(:), dAOrders(:)
    real(8), allocatable :: dAComponents(:), dAUncertainties(:), dAffineComponents(:), dAffineErrors(:)
    real(8), allocatable :: dAffineOrders(:), dAffineSteps(:), dAffineUncertainties(:)
    real(8), allocatable :: dBBaseCopy(:), dBBaseline(:), dBBuilder(:), dBCandidate(:), dBTrial(:)
    real(8), allocatable :: dBComponents(:), dBDirect(:), dBErrors(:), dBOrders(:), dBUncertainties(:)
    real(8), allocatable :: dBestColumnErrors(:), dBestColumnSteps(:), dBestColumnUncertainties(:)
    real(8), allocatable :: dColumnComponents(:,:), dColumnErrors(:,:), dColumnUncertainties(:,:)
    real(8), allocatable :: dChemicalBuilderSave(:), dChemicalSave(:), dEffStoichSave(:,:)
    real(8), allocatable :: dElementBuilderSave(:), dFractionBuilderSave(:), dFractionSave(:)
    real(8), allocatable :: dGibbsBuilderSave(:), dGibbsSave(:)
    real(8), allocatable :: dConstraint(:,:), dDeltaA(:,:), dDeltaB(:), dFDDeltaA(:,:)
    real(8), allocatable :: dAffineReference(:), dFDDeltaB(:), dForceAffine(:), dForcing(:,:)
    real(8), allocatable :: dForceMixed(:), dForceMu(:,:)
    real(8), allocatable :: dForceMuShifted(:,:), dHbase(:,:), dHessian(:,:), dHx(:,:)
    real(8), allocatable :: dKtangent(:,:), dMoles(:), dMolesBuilderSave(:), dMolesSave(:), dMu(:), dMuLowLevel(:)
    real(8), allocatable :: dMuResponse(:,:), dMuResponseBase(:,:), dMuShiftResponse(:,:)
    real(8), allocatable :: dOracleUncertainty(:), dPartialBuilderSave(:), dPartialSave(:)
    real(8), allocatable :: dPhaseMolesBuilderSave(:), dTrialGamma(:), dUpdateSave(:)
    real(8), allocatable :: dResponse(:,:), dResponseBase(:,:), dResponseMinus(:), dResponsePlus(:)
    real(8), allocatable :: dS(:,:), dStepsA(:), dStepsB(:), dX(:), dXMinus(:), dXOff(:), dXPlus(:)
    real(8), allocatable :: dZ(:,:)
    integer, allocatable :: iBestColumn(:)
    logical, allocatable :: lAOrder(:), lAffineOrder(:), lBOrder(:)
    type(FDSweepAssessment) :: tSweepA, tSweepAffine, tSweepB
    type(MQMQAModelData) :: tModel
    type(MQMQAInteractionTerm), allocatable :: tInteraction(:)

    lPass = .TRUE.
    lReport = .FALSE.
    if (COMMAND_ARGUMENT_COUNT() > 0) then
        call GET_COMMAND_ARGUMENT(1,cArgument)
        lReport = TRIM(cArgument) == '--report'
    end if

    !=========================================================================================================
    ! SECTION 1: NATIVE STATE AND LIVE GEMNEWTON BASELINE CAPTURE
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
    if (INFOThermo /= 0) call FinishTest(.FALSE.)

    iPhaseIndex = 0
    iSlot = 0
    do i = 1, nElements
        if (iAssemblage(i) >= 0) cycle
        if ((cSolnPhaseName(-iAssemblage(i)) == 'Liquid') .AND. &
            (cSolnPhaseType(-iAssemblage(i)) == 'SUBG')) then
            iPhaseIndex = -iAssemblage(i)
            iSlot = i
            exit
        end if
    end do
    lPass = lPass .AND. (iPhaseIndex > 0) .AND. (iSlot > 0)
    if (iPhaseIndex <= 0) call FinishTest(lPass)
    ! MQ-4A currently derives only the normalization-constrained response. A charged phase would require the
    ! additional charge constraint used by SubMinNewton and must not enter this prototype silently.
    lPass = lPass .AND. (iPhaseElectronID(iPhaseIndex) == 0)
    if (iPhaseElectronID(iPhaseIndex) /= 0) call FinishTest(.FALSE.)
    iFirst = nSpeciesPhase(iPhaseIndex-1)+1
    iLast = nSpeciesPhase(iPhaseIndex)
    nQuad = iLast-iFirst+1
    dN = dMolesPhase(iSlot)

    ! Exercise the default-disabled diagnostic contract before requesting a real capture.
    call ResetGEMNewtonDiagnosticCapture
    dCaptureControlA = 2D0
    dCaptureControlB = 3D0
    call CaptureGEMNewtonSystem(dCaptureControlA,dCaptureControlB,1)
    lPass = lPass .AND. (.NOT. lCaptureGEMNewtonSystem) .AND. (.NOT. lGEMNewtonSystemCaptured) .AND. &
        (.NOT. ALLOCATED(dCapturedGEMNewtonA)) .AND. (.NOT. ALLOCATED(dCapturedGEMNewtonB))

    allocate(dX(nQuad),dXOff(nQuad),dMoles(nQuad),dMu(nQuad), &
        dMuLowLevel(nQuad),dS(nQuad,nElements),dADirect(nElements,nElements), &
        dAExtracted(nElements,nElements),dAOther(nElements,nElements),dAep(nElements), &
        dBDirect(nElements),dABaseline(nElements,nElements),dBBaseline(nElements), &
        dABuilder(nElements,nElements),dBBuilder(nElements), &
        dChemicalSave(SIZE(dChemicalPotential)),dEffStoichSave(SIZE(dEffStoichSolnPhase,1), &
        SIZE(dEffStoichSolnPhase,2)),dFractionSave(SIZE(dMolFraction)),dGibbsSave(SIZE(dGibbsSolnPhase)), &
        dMolesSave(SIZE(dMolesSpecies)),dPartialSave(SIZE(dPartialExcessGibbs)),dUpdateSave(SIZE(dUpdateVar)), &
        dChemicalBuilderSave(SIZE(dChemicalPotential)),dElementBuilderSave(SIZE(dElementPotential)), &
        dFractionBuilderSave(SIZE(dMolFraction)),dMolesBuilderSave(SIZE(dMolesSpecies)), &
        dPartialBuilderSave(SIZE(dPartialExcessGibbs)),dPhaseMolesBuilderSave(SIZE(dMolesPhase)), &
        dGibbsBuilderSave(SIZE(dGibbsSolnPhase)))
    dX = dMolFraction(iFirst:iLast)
    dX = dX/SUM(dX)
    do e = 1, nElements
        dS(:,e) = dStoichSpecies(iFirst:iLast,e)/DFLOAT(iParticlesPerMole(iFirst:iLast))
    end do

    ! Transfer a small amount between two well-populated quadruplets. Total phase amount and topology remain
    ! fixed, but the production stationarity residual is deliberately no longer at its converged value.
    iDonor = MAXLOC(dX,DIM=1)
    dXOff = dX
    dXOff(iDonor) = -1D0
    iReceiver = MAXLOC(dXOff,DIM=1)
    dXOff = dX
    dH = DMIN1(1D-3,0.02D0*dX(iDonor))
    dXOff(iDonor) = dXOff(iDonor)-dH
    dXOff(iReceiver) = dXOff(iReceiver)+dH
    lPass = lPass .AND. (MINVAL(dXOff) > 1D-12) .AND. (ABS(SUM(dXOff)-1D0) <= 1D-14)

    ! Install the same off-equilibrium state used by the response oracle, recompute production partial molars,
    ! and capture GEMNewton before it solves. All touched global arrays are restored after the phase-local
    ! comparison so this diagnostic cannot alter the converged Thermochimica state seen by later checks.
    dChemicalSave = dChemicalPotential
    dEffStoichSave = dEffStoichSolnPhase
    dFractionSave = dMolFraction
    dGibbsSave = dGibbsSolnPhase
    dMolesSave = dMolesSpecies
    dPartialSave = dPartialExcessGibbs
    dUpdateSave = dUpdateVar
    dMolesSpecies(iFirst:iLast) = dN*dXOff
    call CompChemicalPotential(.FALSE.)
    if (INFOThermo /= 0) call FinishTest(.FALSE.)
    call ResetGEMNewtonDiagnosticCapture
    lCaptureGEMNewtonSystem = .TRUE.
    call GEMNewton(iNewtonInfo)
    lCaptureGEMNewtonSystem = .FALSE.
    lPass = lPass .AND. (iNewtonInfo == 0) .AND. lGEMNewtonSystemCaptured
    if ((iNewtonInfo /= 0) .OR. (.NOT. lGEMNewtonSystemCaptured)) call FinishTest(.FALSE.)
    dFloorDifference = MAXVAL(DABS(dMolesSpecies(iFirst:iLast)-dN*dXOff))/ &
        DMAX1(1D0,MAXVAL(DABS(dN*dXOff)))

    ! The production builder consumes this live off-equilibrium state but must leave every thermodynamic array
    ! untouched. Its outputs are compared below with the independently reconstructed MQ-4A reference and oracle.
    dChemicalBuilderSave = dChemicalPotential
    dElementBuilderSave = dElementPotential
    dFractionBuilderSave = dMolFraction
    dGibbsBuilderSave = dGibbsSolnPhase
    dMolesBuilderSave = dMolesSpecies
    dPartialBuilderSave = dPartialExcessGibbs
    dPhaseMolesBuilderSave = dMolesPhase
    call BuildMQMQAGEMCorrection(iPhaseIndex,iSlot,dABuilder,dBBuilder,lBuilderApplicable,iBuilderStatus)
    dBuilderStateDifference = DMAX1(VectorError(dChemicalPotential,dChemicalBuilderSave), &
        VectorError(dElementPotential,dElementBuilderSave),VectorError(dMolFraction,dFractionBuilderSave), &
        VectorError(dGibbsSolnPhase,dGibbsBuilderSave),VectorError(dMolesSpecies,dMolesBuilderSave), &
        VectorError(dPartialExcessGibbs,dPartialBuilderSave),VectorError(dMolesPhase,dPhaseMolesBuilderSave))
    lPass = lPass .AND. lBuilderApplicable .AND. (iBuilderStatus == MQMQA_MAP_SUCCESS) .AND. &
        (dBuilderStateDifference == 0D0) .AND. &
        (MAXVAL(ABS(dABuilder-TRANSPOSE(dABuilder))) == 0D0)

    ! Remove every other active solution-phase contribution from the captured element equations. The remainder
    ! is a true live snapshot of the selected Liquid contribution, including GEMNewton's species-mole floor.
    dAExtracted = dCapturedGEMNewtonA(1:nElements,1:nElements)
    dBDirect = dCapturedGEMNewtonB(1:nElements)-dMolesElement(1:nElements)
    do i = 1, nSolnPhases
        iColumn = nElements+i
        f = 2*nElements-iColumn+1
        e = -iAssemblage(f)
        if (e == iPhaseIndex) cycle
        call DirectPhaseElementTerms(e,dAOther,dBBaseline)
        dAExtracted = dAExtracted-dAOther
        dBDirect = dBDirect-dBBaseline
    end do
    call DirectPhaseElementTerms(iPhaseIndex,dADirect,dBBaseline)
    dBaselineAError = MatrixError(dAExtracted,dADirect)
    dBaselineBError = VectorError(dBDirect,dBBaseline)

    iColumn = nElements+(nElements-iSlot+1)
    call CompStoichSolnPhase(iPhaseIndex)
    dAep = dN*dEffStoichSolnPhase(iPhaseIndex,1:nElements)
    dBaselineColumnError = DMAX1(VectorError(dCapturedGEMNewtonA(1:nElements,iColumn),dAep), &
        VectorError(dCapturedGEMNewtonA(iColumn,1:nElements),dAep))
    dBaselinePhaseError = ABS(dCapturedGEMNewtonB(iColumn)-dGibbsSolnPhase(iPhaseIndex))/ &
        DMAX1(1D0,ABS(dCapturedGEMNewtonB(iColumn)),ABS(dGibbsSolnPhase(iPhaseIndex)))
    lPass = lPass .AND. (dFloorDifference <= 1D-10) .AND. (dBaselineAError <= 1D-12) .AND. &
        (dBaselineBError <= 1D-12) .AND. (dBaselineColumnError <= 1D-12) .AND. &
        (dBaselinePhaseError <= 1D-12)

    dChemicalPotential = dChemicalSave
    dEffStoichSolnPhase = dEffStoichSave
    dMolFraction = dFractionSave
    dGibbsSolnPhase = dGibbsSave
    dMolesSpecies = dMolesSave
    dPartialExcessGibbs = dPartialSave
    dUpdateVar = dUpdateSave

    !=========================================================================================================
    ! SECTION 2: OFF-EQUILIBRIUM CANDIDATE DELTAS
    !=========================================================================================================
    call DecodeProductionSUBGPhase(iPhaseIndex,tModel,tInteraction,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    if (iInfo /= 0) call FinishTest(.FALSE.)
    dMoles = dN*dXOff
    allocate(dHessian(nQuad,nQuad),dHx(nQuad,nQuad),dHbase(nQuad,nQuad), &
        dConstraint(1,nQuad),dForcing(nQuad,nElements),dResponse(nQuad,nElements), &
        dResponseBase(nQuad,nElements),dForceMu(nQuad,1),dForceMuShifted(nQuad,1), &
        dMuResponse(nQuad,1),dMuResponseBase(nQuad,1),dMuShiftResponse(nQuad,1), &
        dDeltaA(nElements,nElements),dDeltaB(nElements),dACandidate(nElements,nElements), &
        dBCandidate(nElements),dKtangent(nQuad-1,nQuad-1),dZ(nQuad,nQuad-1))
    call CompMQMQAHessianUnconstrained(tModel,dMoles,1D0,tInteraction,dHessian,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    if (iInfo /= 0) call FinishTest(.FALSE.)
    dHx = dN*dHessian
    dHbase = 0D0
    do i = 1, nQuad
        dHbase(i,i) = 1D0/dXOff(i)
    end do
    dConstraint = 1D0
    dForcing = dS
    ! Helmert contrasts provide an orthonormal basis for all normalized composition changes.
    dZ = 0D0
    do i = 1, nQuad-1
        dZ(1:i,i) = 1D0/DSQRT(DFLOAT(i*(i+1)))
        dZ(i+1,i) = -DFLOAT(i)/DSQRT(DFLOAT(i*(i+1)))
    end do
    call EvaluateProductionMu(iPhaseIndex,iFirst,iLast,dXOff,dMuLowLevel,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    if (iInfo /= 0) call FinishTest(.FALSE.)
    dMu = dMuLowLevel
    dOffEquilibriumResidual = TangentNorm(dMu-MATMUL(dS,dElementPotential),dZ)

    call SolveConstrainedResponse(dHx,dConstraint,dForcing,dResponse,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    call SolveConstrainedResponse(dHbase,dConstraint,dForcing,dResponseBase,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    dForceMu(:,1) = dMu-1D0
    dForceMuShifted(:,1) = dMu
    call SolveConstrainedResponse(dHx,dConstraint,dForceMu,dMuResponse,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    call SolveConstrainedResponse(dHbase,dConstraint,dForceMu,dMuResponseBase,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    call SolveConstrainedResponse(dHx,dConstraint,dForceMuShifted,dMuShiftResponse,iInfo)
    lPass = lPass .AND. (iInfo == 0)

    dDeltaA = dN*MATMUL(TRANSPOSE(dS),dResponse-dResponseBase)
    dDeltaB = dN*MATMUL(TRANSPOSE(dS),dMuResponse(:,1)-dMuResponseBase(:,1))
    dBuilderAError = MatrixError(dABuilder,dDeltaA)
    dBuilderBError = VectorError(dBBuilder,dDeltaB)
    lPass = lPass .AND. (dBuilderAError <= 1D-12) .AND. (dBuilderBError <= 1D-12)
    dSymmetryResidual = SQRT(SUM((dDeltaA-TRANSPOSE(dDeltaA))**2))/ &
        DMAX1(1D0,SQRT(SUM(dDeltaA*dDeltaA)))
    dConstraintResidual = DMAX1(MAXVAL(DABS(MATMUL(dConstraint,dResponse))), &
        DABS(SUM(dMuResponse(:,1))))
    dConstantForceError = VectorError(dMuResponse(:,1),dMuShiftResponse(:,1))
    dDeltaBMagnitude = MAXVAL(DABS(dDeltaB))

    ! If corrected curvature is replaced by the historical ideal curvature, both corrections must vanish.
    dACandidate = dN*MATMUL(TRANSPOSE(dS),dResponseBase-dResponseBase)
    dBCandidate = dN*MATMUL(TRANSPOSE(dS),dMuResponseBase(:,1)-dMuResponseBase(:,1))
    lPass = lPass .AND. (MAXVAL(DABS(dACandidate)) == 0D0) .AND. &
        (MAXVAL(DABS(dBCandidate)) == 0D0) .AND. (dSymmetryResidual <= 1D-10) .AND. &
        (dConstraintResidual <= 1D-10) .AND. (dConstantForceError <= 1D-10) .AND. &
        (dDeltaBMagnitude > 1D-10) .AND. (dOffEquilibriumResidual > 1D-8) .AND. &
        ALL(IEEE_IS_FINITE(dDeltaA)) .AND. ALL(IEEE_IS_FINITE(dDeltaB))

    ! Inapplicable, invalid, unsupported, and singular-response cases must leave safe zero outputs with distinct
    ! statuses. The zero-amount kernel call exercises INVALID_INPUT independently of ordinary inapplicability.
    dACandidate = 1D0
    dBCandidate = 1D0
    call BuildMQMQAGEMCorrection(0,iSlot,dACandidate,dBCandidate,lFailureApplicable,iFailureStatus)
    lPass = lPass .AND. (.NOT. lFailureApplicable) .AND. (iFailureStatus == MQMQA_MAP_NOT_APPLICABLE) .AND. &
        (MAXVAL(ABS(dACandidate)) == 0D0) .AND. (MAXVAL(ABS(dBCandidate)) == 0D0)

    dACandidate = 1D0
    dBCandidate = 1D0
    call BuildMQMQAReducedCorrection(0D0,dXOff,dMu,dS,dHx,dACandidate,dBCandidate,iFailureStatus)
    lPass = lPass .AND. (iFailureStatus == MQMQA_MAP_INVALID_INPUT) .AND. &
        (MAXVAL(ABS(dACandidate)) == 0D0) .AND. (MAXVAL(ABS(dBCandidate)) == 0D0)

    iElectronSave = iPhaseElectronID(iPhaseIndex)
    iPhaseElectronID(iPhaseIndex) = 1
    dACandidate = 1D0
    dBCandidate = 1D0
    call BuildMQMQAGEMCorrection(iPhaseIndex,iSlot,dACandidate,dBCandidate,lFailureApplicable,iFailureStatus)
    iPhaseElectronID(iPhaseIndex) = iElectronSave
    lPass = lPass .AND. (.NOT. lFailureApplicable) .AND. &
        (iFailureStatus == MQMQA_MAP_UNSUPPORTED_CHARGED_PHASE) .AND. &
        (MAXVAL(ABS(dACandidate)) == 0D0) .AND. (MAXVAL(ABS(dBCandidate)) == 0D0)

    dACandidate = 1D0
    dBCandidate = 1D0
    call BuildMQMQAReducedCorrection(dN,dXOff,dMu,dS,0D0*dHx,dACandidate,dBCandidate,iFailureStatus)
    lPass = lPass .AND. (iFailureStatus == MQMQA_MAP_CORRECTED_ELEMENT_RESPONSE_FAILURE) .AND. &
        (MAXVAL(ABS(dACandidate)) == 0D0) .AND. (MAXVAL(ABS(dBCandidate)) == 0D0)

    ! Apply only to test-owned copies of the captured GEM system. Full application, alpha-zero, and sequential
    ! correction-pair additivity must preserve every row, column, and residual outside the element equations.
    allocate(dABaseCopy(nCapturedGEMNewtonVariables,nCapturedGEMNewtonVariables), &
        dAApplyExpected(nCapturedGEMNewtonVariables,nCapturedGEMNewtonVariables), &
        dATrial(nCapturedGEMNewtonVariables,nCapturedGEMNewtonVariables), &
        dBBaseCopy(nCapturedGEMNewtonVariables),dBTrial(nCapturedGEMNewtonVariables))
    dABaseCopy = dCapturedGEMNewtonA
    dBBaseCopy = dCapturedGEMNewtonB
    dAApplyExpected = dABaseCopy
    dAApplyExpected(1:nElements,1:nElements) = dAApplyExpected(1:nElements,1:nElements)+dABuilder

    dATrial = dABaseCopy
    dBTrial = dBBaseCopy
    call ApplyMQMQAGEMCorrection(dATrial,dBTrial,nElements,dABuilder,dBBuilder,1D0,iApplyStatus)
    dApplyError = DMAX1(MatrixError(dATrial,dAApplyExpected), &
        VectorError(dBTrial(1:nElements),dBBaseCopy(1:nElements)+dBBuilder), &
        VectorError(dBTrial(nElements+1:),dBBaseCopy(nElements+1:)))
    lPass = lPass .AND. (iApplyStatus == MQMQA_MAP_SUCCESS) .AND. (dApplyError <= 1D-14)

    dATrial = dABaseCopy
    dBTrial = dBBaseCopy
    call ApplyMQMQAGEMCorrection(dATrial,dBTrial,nElements,dABuilder,dBBuilder,0D0,iApplyStatus)
    dZeroApplyError = DMAX1(MatrixError(dATrial,dABaseCopy),VectorError(dBTrial,dBBaseCopy))
    lPass = lPass .AND. (iApplyStatus == MQMQA_MAP_SUCCESS) .AND. (dZeroApplyError == 0D0)

    dATrial = dABaseCopy
    dBTrial = dBBaseCopy
    dACandidate = -0.25D0*dABuilder
    dBCandidate = 0.40D0*dBBuilder
    call ApplyMQMQAGEMCorrection(dATrial,dBTrial,nElements,dABuilder,dBBuilder,1D0,iApplyStatus)
    lPass = lPass .AND. (iApplyStatus == MQMQA_MAP_SUCCESS)
    call ApplyMQMQAGEMCorrection(dATrial,dBTrial,nElements,dACandidate,dBCandidate,1D0,iApplyStatus)
    dAApplyExpected = dABaseCopy
    dAApplyExpected(1:nElements,1:nElements) = dAApplyExpected(1:nElements,1:nElements)+ &
        dABuilder+dACandidate
    dAggregationError = DMAX1(MatrixError(dATrial,dAApplyExpected), &
        VectorError(dBTrial(1:nElements),dBBaseCopy(1:nElements)+dBBuilder+dBCandidate), &
        VectorError(dBTrial(nElements+1:),dBBaseCopy(nElements+1:)))
    lPass = lPass .AND. (iApplyStatus == MQMQA_MAP_SUCCESS) .AND. (dAggregationError <= 1D-14)

    ! The public applicator accepts only a symmetric element-block correction. A materially nonsymmetric caller
    ! input must be rejected before either caller-owned array is changed.
    dATrial = dABaseCopy
    dBTrial = dBBaseCopy
    dACandidate = dABuilder
    dACandidate(1,2) = dACandidate(1,2)+1D0
    call ApplyMQMQAGEMCorrection(dATrial,dBTrial,nElements,dACandidate,dBBuilder,1D0,iApplyStatus)
    lPass = lPass .AND. (iApplyStatus == MQMQA_MAP_INVALID_APPLICATION) .AND. &
        (MatrixError(dATrial,dABaseCopy) == 0D0) .AND. (VectorError(dBTrial,dBBaseCopy) == 0D0)

    dATrial = dABaseCopy
    dBTrial = dBBaseCopy
    call ApplyMQMQAGEMCorrection(dATrial,dBTrial,nElements,dABuilder,dBBuilder,1.5D0,iApplyStatus)
    lPass = lPass .AND. (iApplyStatus == MQMQA_MAP_INVALID_APPLICATION) .AND. &
        (MatrixError(dATrial,dABaseCopy) == 0D0) .AND. (VectorError(dBTrial,dBBaseCopy) == 0D0)

    !=========================================================================================================
    ! SECTION 3: INDEPENDENT NONLINEAR FINITE-DIFFERENCE ORACLES
    !=========================================================================================================
    ! Helmert contrasts span every total-preserving composition change and are used only by the numerical
    ! production-partial-molar oracle. They do not call the analytic response solver while finding states.
    dKtangent = MATMUL(TRANSPOSE(dZ),MATMUL(dHx,dZ))
    allocate(dForceMixed(nQuad),dForceAffine(nQuad),dAffineReference(nElements),dTrialGamma(nElements), &
        dResponsePlus(nQuad),dResponseMinus(nQuad),dOracleUncertainty(nElements), &
        dXPlus(nQuad),dXMinus(nQuad),dStepsA(nSteps),dStepsB(nSteps),dAErrors(nSteps), &
        dBErrors(nSteps),dAComponents(nSteps),dBComponents(nSteps),dAUncertainties(nSteps), &
        dBUncertainties(nSteps),dAOrders(nSteps),dBOrders(nSteps),lAOrder(nSteps),lBOrder(nSteps), &
        dAffineSteps(nSteps),dAffineErrors(nSteps),dAffineComponents(nSteps),dAffineUncertainties(nSteps), &
        dAffineOrders(nSteps),lAffineOrder(nSteps),dFDDeltaA(nElements,nElements),dFDDeltaB(nElements), &
        dColumnErrors(nSteps,nElements),dColumnComponents(nSteps,nElements), &
        dColumnUncertainties(nSteps,nElements),dBestColumnErrors(nElements), &
        dBestColumnUncertainties(nElements),dBestColumnSteps(nElements),iBestColumn(nElements))
    dTrialGamma = [(1D0/DFLOAT(i),i=1,nElements)]
    dForceMixed = MATMUL(dS,dTrialGamma)
    dForceAffine = dForceMixed-dForceMu(:,1)
    dAffineReference = MATMUL(dDeltaA,dTrialGamma)-dDeltaB
    dHMaximumA = PositivityStep(dXOff,MATMUL(dResponse-dResponseBase, &
        dTrialGamma))
    dHMaximumB = PositivityStep(dXOff,dMuResponse(:,1)-dMuResponseBase(:,1))
    dHMaximumAffine = PositivityStep(dXOff,MATMUL(dResponse-dResponseBase,dTrialGamma)- &
        (dMuResponse(:,1)-dMuResponseBase(:,1)))

    do iStep = 1, nSteps
        dStepsA(iStep) = dHMaximumA*3D0**(-(iStep-1))
        call EvaluateDeltaResponse(dXOff,dForceMixed,dStepsA(iStep),dXPlus,dXMinus, &
            dResponsePlus,dResponseMinus,lPlus,lMinus,dOracleUncertainty)
        lPass = lPass .AND. lPlus .AND. lMinus
        dFDDeltaB = dN*MATMUL(TRANSPOSE(dS), &
            ((dXPlus-dXMinus)-(dResponsePlus-dResponseMinus))/(2D0*dStepsA(iStep)))
        call VectorMetrics(dFDDeltaB,MATMUL(dDeltaA,dTrialGamma), &
            dAError,dAComponent)
        dAErrors(iStep) = dAError
        dAComponents(iStep) = dAComponent
        dAUncertainties(iStep) = ScaledUncertainty(dOracleUncertainty,dFDDeltaB, &
            MATMUL(dDeltaA,dTrialGamma))

        dStepsB(iStep) = dHMaximumB*3D0**(-(iStep-1))
        call EvaluateDeltaResponse(dXOff,dForceMu(:,1),dStepsB(iStep),dXPlus,dXMinus, &
            dResponsePlus,dResponseMinus,lPlus,lMinus,dOracleUncertainty)
        lPass = lPass .AND. lPlus .AND. lMinus
        dFDDeltaB = dN*MATMUL(TRANSPOSE(dS), &
            ((dXPlus-dXMinus)-(dResponsePlus-dResponseMinus))/(2D0*dStepsB(iStep)))
        call VectorMetrics(dFDDeltaB,dDeltaB,dBError,dBComponent)
        dBErrors(iStep) = dBError
        dBComponents(iStep) = dBComponent
        dBUncertainties(iStep) = ScaledUncertainty(dOracleUncertainty,dFDDeltaB,dDeltaB)

        ! A combined forcing tests the affine reduced equation in one experiment. Its analytic response is
        ! deltaA*dGamma-deltaB, so a sign error in either correction cannot be hidden by separate checks.
        dAffineSteps(iStep) = dHMaximumAffine*3D0**(-(iStep-1))
        call EvaluateDeltaResponse(dXOff,dForceAffine,dAffineSteps(iStep),dXPlus,dXMinus, &
            dResponsePlus,dResponseMinus,lPlus,lMinus,dOracleUncertainty)
        lPass = lPass .AND. lPlus .AND. lMinus
        dFDDeltaB = dN*MATMUL(TRANSPOSE(dS), &
            ((dXPlus-dXMinus)-(dResponsePlus-dResponseMinus))/(2D0*dAffineSteps(iStep)))
        call VectorMetrics(dFDDeltaB,dAffineReference,dAffineError,dAffineComponent)
        dAffineErrors(iStep) = dAffineError
        dAffineComponents(iStep) = dAffineComponent
        dAffineUncertainties(iStep) = ScaledUncertainty(dOracleUncertainty,dFDDeltaB,dAffineReference)
    end do
    call AssessFDSweep(dStepsA,dAErrors,FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX,1D-6, &
        tSweepA,dAOrders,lAOrder)
    call AssessFDSweep(dStepsB,dBErrors,FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX,1D-6, &
        tSweepB,dBOrders,lBOrder)
    call AssessFDSweep(dAffineSteps,dAffineErrors,FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX,1D-6, &
        tSweepAffine,dAffineOrders,lAffineOrder)
    call BestResolved(dAErrors,dAUncertainties,iBestA,dAError)
    call BestResolved(dBErrors,dBUncertainties,iBestB,dBError)
    call BestResolved(dAffineErrors,dAffineUncertainties,iBestAffine,dAffineError)
    lPass = lPass .AND. tSweepA%lPassed .AND. tSweepB%lPassed .AND. tSweepAffine%lPassed .AND. &
        (iBestA > 0) .AND. (iBestB > 0) .AND. (iBestAffine > 0) .AND. &
        (dAError <= 1D-6) .AND. (dBError <= 1D-6) .AND. (dAffineError <= 1D-6)
    lPass = lPass .AND. HasResolvedSecondOrder(dAErrors,dAUncertainties,dAOrders,lAOrder) .AND. &
        HasResolvedSecondOrder(dBErrors,dBUncertainties,dBOrders,lBOrder) .AND. &
        HasResolvedSecondOrder(dAffineErrors,dAffineUncertainties,dAffineOrders,lAffineOrder)
    if (iBestA > 0) lPass = lPass .AND. (dAComponents(iBestA) <= 1D-6)
    if (iBestB > 0) lPass = lPass .AND. (dBComponents(iBestB) <= 1D-6)
    if (iBestAffine > 0) lPass = lPass .AND. (dAffineComponents(iBestAffine) <= 1D-6)

    ! Reconstruct every element column with its own oracle-resolved step. Storing and checking each column's
    ! mapped uncertainty prevents an unresolved direction from disappearing inside a small matrix norm.
    do e = 1, nElements
        dHColumn = PositivityStep(dXOff,dResponse(:,e)-dResponseBase(:,e))
        do iStep = 1, nSteps
            dH = dHColumn*3D0**(-(iStep-1))
            call EvaluateDeltaResponse(dXOff,dS(:,e),dH,dXPlus,dXMinus, &
                dResponsePlus,dResponseMinus,lPlus,lMinus,dOracleUncertainty)
            lPass = lPass .AND. lPlus .AND. lMinus
            dFDDeltaB = dN*MATMUL(TRANSPOSE(dS), &
                ((dXPlus-dXMinus)-(dResponsePlus-dResponseMinus))/(2D0*dH))
            call VectorMetrics(dFDDeltaB,dDeltaA(:,e),dAError,dAComponent)
            dColumnErrors(iStep,e) = dAError
            dColumnComponents(iStep,e) = dAComponent
            dColumnUncertainties(iStep,e) = ScaledUncertainty(dOracleUncertainty,dFDDeltaB,dDeltaA(:,e))
        end do
        call BestResolved(dColumnErrors(:,e),dColumnUncertainties(:,e),iBestColumn(e),dBestColumnErrors(e))
        if (iBestColumn(e) > 0) then
            dBestColumnUncertainties(e) = dColumnUncertainties(iBestColumn(e),e)
            dBestColumnSteps(e) = dHColumn*3D0**(-(iBestColumn(e)-1))
            dH = dBestColumnSteps(e)
            call EvaluateDeltaResponse(dXOff,dS(:,e),dH,dXPlus,dXMinus, &
                dResponsePlus,dResponseMinus,lPlus,lMinus,dOracleUncertainty)
            dFDDeltaA(:,e) = dN*MATMUL(TRANSPOSE(dS), &
                ((dXPlus-dXMinus)-(dResponsePlus-dResponseMinus))/(2D0*dH))
            lPass = lPass .AND. lPlus .AND. lMinus .AND. (dBestColumnErrors(e) <= 1D-6) .AND. &
                (dColumnComponents(iBestColumn(e),e) <= 1D-6) .AND. &
                (dBestColumnUncertainties(e) <= 0.5D0*dBestColumnErrors(e))
        else
            dBestColumnErrors(e) = HUGE(1D0)
            dBestColumnUncertainties(e) = HUGE(1D0)
            dBestColumnSteps(e) = 0D0
            lPass = .FALSE.
        end if
    end do
    dAError = MatrixError(dFDDeltaA,dDeltaA)
    dAComponent = MAXVAL(ABS(dFDDeltaA-dDeltaA)/MAX(1D0,ABS(dFDDeltaA),ABS(dDeltaA)))
    dWorstColumnError = MAXVAL(dBestColumnErrors)
    dWorstColumnUncertainty = MAXVAL(dBestColumnUncertainties)
    lPass = lPass .AND. (dAError <= 1D-6) .AND. (dAComponent <= 1D-6)

    if (lReport) then
        write(*,'(A)') 'MQ-4A/4B plain-SUBG reduced GEM mapping verification'
        write(*,'(A)') 'scope: reusable builder and copied-array applicator; no live GEM mutation or globalization'
        write(*,'(A,A)') 'phase = ',TRIM(cSolnPhaseName(iPhaseIndex))
        write(*,'(A,ES14.6)') 'authoritative/floored mole discrepancy = ',dFloorDifference
        write(*,'(A,ES14.6)') 'off-equilibrium live phase-local element-block capture error = ',dBaselineAError
        write(*,'(A,ES14.6)') 'off-equilibrium live phase-local element-residual capture error = ',dBaselineBError
        write(*,'(A,ES14.6)') 'off-equilibrium live element-to-phase column capture error = ',dBaselineColumnError
        write(*,'(A,ES14.6)') 'off-equilibrium live solution-phase residual capture error = ',dBaselinePhaseError
        write(*,'(A,ES14.6)') 'off-equilibrium tangent stationarity residual = ',dOffEquilibriumResidual
        write(*,'(A,ES14.6)') 'mu versus mu-1 constrained-response difference = ',dConstantForceError
        write(*,'(A,ES14.6)') 'candidate deltaA symmetry residual = ',dSymmetryResidual
        write(*,'(A,ES14.6)') 'candidate response constraint residual = ',dConstraintResidual
        write(*,'(A,ES14.6)') 'maximum absolute candidate deltaB = ',dDeltaBMagnitude
        write(*,'(A,I0,A,L1)') 'builder status = ',iBuilderStatus,', applicable = ',lBuilderApplicable
        write(*,'(A,ES14.6)') 'builder versus independent deltaA error = ',dBuilderAError
        write(*,'(A,ES14.6)') 'builder versus independent deltaB error = ',dBuilderBError
        write(*,'(A,ES14.6)') 'builder production-state mutation metric = ',dBuilderStateDifference
        write(*,'(A,ES14.6)') 'copied-array full-application error = ',dApplyError
        write(*,'(A,ES14.6)') 'copied-array alpha-zero error = ',dZeroApplyError
        write(*,'(A,ES14.6)') 'copied-array sequential correction-pair additivity error = ',dAggregationError
        write(*,'(A,ES14.6)') 'complete reduced deltaA FD error = ',dAError
        write(*,'(A,ES14.6)') 'complete reduced deltaA maximum scaled component error = ',dAComponent
        write(*,'(A,ES14.6)') 'worst oracle-resolved deltaA column error = ',dWorstColumnError
        write(*,'(A,ES14.6)') 'worst selected deltaA column uncertainty = ',dWorstColumnUncertainty
        do e = 1, nElements
            write(*,'(A,I0,3(A,ES14.6))') 'deltaA column ',e,': step = ',dBestColumnSteps(e), &
                ', error = ',dBestColumnErrors(e),', uncertainty = ',dBestColumnUncertainties(e)
        end do
        call PrintSweep('deltaA mixed element forcing',dStepsA,dAErrors,dAUncertainties,dAOrders,lAOrder,iBestA)
        call PrintSweep('deltaB off-equilibrium residual forcing',dStepsB,dBErrors,dBUncertainties, &
            dBOrders,lBOrder,iBestB)
        call PrintSweep('combined affine forcing: deltaA*dGamma-deltaB',dAffineSteps,dAffineErrors, &
            dAffineUncertainties,dAffineOrders,lAffineOrder,iBestAffine)
        if (iBestA > 0) write(*,'(A,ES14.6)') 'best resolved deltaA maximum scaled component error = ', &
            dAComponents(iBestA)
        if (iBestB > 0) write(*,'(A,ES14.6)') 'best resolved deltaB maximum scaled component error = ', &
            dBComponents(iBestB)
        if (iBestAffine > 0) write(*,'(A,ES14.6)') 'best resolved affine maximum scaled component error = ', &
            dAffineComponents(iBestAffine)
    end if

    call FinishTest(lPass)

contains

    !> Reproduce one phase's exact GEMNewton element-block and element-residual loop from live global arrays.
    subroutine DirectPhaseElementTerms(iPhaseLocal,dAPhase,dBPhase)
        integer, intent(in) :: iPhaseLocal
        real(8), intent(out) :: dAPhase(:,:), dBPhase(:)
        integer :: ee, ff, q, iFirstLocal, iLastLocal
        real(8) :: dSq, dSf
        iFirstLocal = nSpeciesPhase(iPhaseLocal-1)+1
        iLastLocal = nSpeciesPhase(iPhaseLocal)
        dAPhase = 0D0
        dBPhase = 0D0
        do q = iFirstLocal, iLastLocal
            do ee = 1, nElements
                dSq = dStoichSpecies(q,ee)/DFLOAT(iParticlesPerMole(q))
                dBPhase(ee) = dBPhase(ee)+dMolesSpecies(q)*dSq*(dChemicalPotential(q)-1D0)
                do ff = 1, nElements
                    dSf = dStoichSpecies(q,ff)/DFLOAT(iParticlesPerMole(q))
                    dAPhase(ee,ff) = dAPhase(ee,ff)+dMolesSpecies(q)*dSq*dSf
                end do
            end do
        end do
    end subroutine DirectPhaseElementTerms

    !> Obtain the complete production SUBG partial molar from the low-level routine without double counting.
    subroutine EvaluateProductionMu(iPhaseLocal,iFirstLocal,iLastLocal,dXLocal,dMuLocal,iInfoLocal)
        integer, intent(in) :: iPhaseLocal, iFirstLocal, iLastLocal
        real(8), intent(in) :: dXLocal(:)
        real(8), intent(out) :: dMuLocal(:)
        integer, intent(out) :: iInfoLocal
        iInfoLocal = 0
        dMolFraction(iFirstLocal:iLastLocal) = dXLocal
        call CompExcessGibbsEnergySUBG(iPhaseLocal)
        dMuLocal = dChemicalPotential(iFirstLocal:iLastLocal)+dPartialExcessGibbs(iFirstLocal:iLastLocal)
        if ((INFOThermo /= 0) .OR. (.NOT. ALL(IEEE_IS_FINITE(dMuLocal)))) iInfoLocal = 1
    end subroutine EvaluateProductionMu

    !> Compare corrected and ideal nonlinear composition responses for one species-level forcing.
    subroutine EvaluateDeltaResponse(dXBase,dForce,dAmplitude,dCorrPlus,dCorrMinus,dIdealPlus,dIdealMinus, &
        lPlusOK,lMinusOK,dMappedUncertainty)
        real(8), intent(in) :: dXBase(:), dForce(:), dAmplitude
        real(8), intent(out) :: dCorrPlus(:), dCorrMinus(:), dIdealPlus(:), dIdealMinus(:)
        real(8), intent(out) :: dMappedUncertainty(:)
        logical, intent(out) :: lPlusOK, lMinusOK
        real(8), allocatable :: dResidualPlus(:), dResidualMinus(:), dCorrection(:,:), dMappedCorrection(:,:)
        integer :: iInfoLocal
        allocate(dResidualPlus(SIZE(dXBase)),dResidualMinus(SIZE(dXBase)), &
            dCorrection(SIZE(dZ,2),2),dMappedCorrection(nElements,2))
        call SolveProductionForcedState(dXBase,dForce,dAmplitude,dCorrPlus,dResidualPlus,lPlusOK)
        call SolveProductionForcedState(dXBase,dForce,-dAmplitude,dCorrMinus,dResidualMinus,lMinusOK)
        if (lReport .AND. ((.NOT. lPlusOK) .OR. (.NOT. lMinusOK))) then
            write(*,'(A,ES14.6,2(A,L1),2(A,ES14.6))') 'oracle failure at h = ',dAmplitude, &
                ', plus = ',lPlusOK,', minus = ',lMinusOK,', plus residual = ', &
                TangentNorm(dResidualPlus,dZ),', minus residual = ',TangentNorm(dResidualMinus,dZ)
        end if
        dIdealPlus = dXBase*EXP(dAmplitude*dForce)
        dIdealPlus = dIdealPlus/SUM(dIdealPlus)
        dIdealMinus = dXBase*EXP(-dAmplitude*dForce)
        dIdealMinus = dIdealMinus/SUM(dIdealMinus)
        dCorrection(:,1) = MATMUL(TRANSPOSE(dZ),dResidualPlus)
        dCorrection(:,2) = MATMUL(TRANSPOSE(dZ),dResidualMinus)
        call SolveDense(dKtangent,dCorrection,iInfoLocal)
        if (iInfoLocal /= 0) then
            dMappedUncertainty = HUGE(1D0)
        else
            ! Convert the residual-implied composition corrections through the same N*S^T map as the measured
            ! GEM residual. Resolution is therefore judged in element space, not against an unlike x-space norm.
            dMappedCorrection(:,1) = dN*MATMUL(TRANSPOSE(dS),MATMUL(dZ,dCorrection(:,1)))
            dMappedCorrection(:,2) = dN*MATMUL(TRANSPOSE(dS),MATMUL(dZ,dCorrection(:,2)))
            dMappedUncertainty = (DABS(dMappedCorrection(:,1))+DABS(dMappedCorrection(:,2)))/ &
                (2D0*DABS(dAmplitude))
        end if
        deallocate(dResidualPlus,dResidualMinus,dCorrection,dMappedCorrection)
    end subroutine EvaluateDeltaResponse

    !> Re-equilibrate production partial-molar differences without using the analytic MQMQA Hessian.
    subroutine SolveProductionForcedState(dXBase,dForce,dAmplitude,dXSolved,dResidual,lSolved)
        real(8), intent(in) :: dXBase(:), dForce(:), dAmplitude
        real(8), intent(out) :: dXSolved(:), dResidual(:)
        logical, intent(out) :: lSolved
        integer :: iDirection, iInfoLocal, iIteration, iTrial
        real(8) :: dAlpha, dDifference, dNorm, dTrialNorm
        real(8), allocatable :: dDelta(:), dJacobian(:,:), dMuBase(:), dMuMinus(:), dMuPlus(:)
        real(8), allocatable :: dMuTrial(:), dReduced(:,:), dResidualTrial(:), dXTrial(:)
        allocate(dDelta(nQuad),dJacobian(nQuad-1,nQuad-1),dMuBase(nQuad),dMuMinus(nQuad), &
            dMuPlus(nQuad),dMuTrial(nQuad),dReduced(nQuad-1,1),dResidualTrial(nQuad),dXTrial(nQuad))
        call EvaluateProductionMu(iPhaseIndex,iFirst,iLast,dXBase,dMuBase,iInfoLocal)
        dXSolved = dXBase
        lSolved = .FALSE.
        do iIteration = 1, 40
            call EvaluateProductionMu(iPhaseIndex,iFirst,iLast,dXSolved,dMuTrial,iInfoLocal)
            if (iInfoLocal /= 0) exit
            dResidual = dMuTrial-dMuBase-dAmplitude*dForce
            dResidual = dResidual-SUM(dResidual)/DFLOAT(nQuad)
            dReduced(:,1) = MATMUL(TRANSPOSE(dZ),dResidual)
            dNorm = SQRT(SUM(dReduced(:,1)**2))
            if (dNorm <= 1D-12) then
                lSolved = .TRUE.
                exit
            end if
            do iDirection = 1, nQuad-1
                where (DABS(dZ(:,iDirection)) > 0D0)
                    dDelta = dXSolved/DABS(dZ(:,iDirection))
                elsewhere
                    dDelta = HUGE(1D0)
                end where
                dDifference = DMIN1(1D-5,0.1D0*MINVAL(dDelta))
                call EvaluateProductionMu(iPhaseIndex,iFirst,iLast,dXSolved-dDifference*dZ(:,iDirection), &
                    dMuMinus,iInfoLocal)
                if (iInfoLocal /= 0) exit
                call EvaluateProductionMu(iPhaseIndex,iFirst,iLast,dXSolved+dDifference*dZ(:,iDirection), &
                    dMuPlus,iInfoLocal)
                if (iInfoLocal /= 0) exit
                dJacobian(:,iDirection) = MATMUL(TRANSPOSE(dZ),dMuPlus-dMuMinus)/(2D0*dDifference)
            end do
            if (iInfoLocal /= 0) exit
            dReduced = -dReduced
            call SolveDense(dJacobian,dReduced,iInfoLocal)
            if (iInfoLocal /= 0) exit
            dDelta = MATMUL(dZ,dReduced(:,1))
            dAlpha = 1D0
            do i = 1, nQuad
                if (dDelta(i) < 0D0) dAlpha = DMIN1(dAlpha,-0.9D0*dXSolved(i)/dDelta(i))
            end do
            dTrialNorm = HUGE(1D0)
            do iTrial = 1, 20
                dXTrial = dXSolved+dAlpha*dDelta
                call EvaluateProductionMu(iPhaseIndex,iFirst,iLast,dXTrial,dMuTrial,iInfoLocal)
                if (iInfoLocal == 0) then
                    dResidualTrial = dMuTrial-dMuBase-dAmplitude*dForce
                    dResidualTrial = dResidualTrial-SUM(dResidualTrial)/DFLOAT(nQuad)
                    dTrialNorm = SQRT(SUM(MATMUL(TRANSPOSE(dZ),dResidualTrial)**2))
                    if (dTrialNorm < dNorm) exit
                end if
                dAlpha = 0.5D0*dAlpha
            end do
            if ((iInfoLocal /= 0) .OR. (dTrialNorm >= dNorm)) exit
            dXSolved = dXTrial
        end do
        call EvaluateProductionMu(iPhaseIndex,iFirst,iLast,dXSolved,dMuTrial,iInfoLocal)
        dResidual = dMuTrial-dMuBase-dAmplitude*dForce
        dResidual = dResidual-SUM(dResidual)/DFLOAT(nQuad)
        lSolved = (iInfoLocal == 0) .AND. (SQRT(SUM(MATMUL(TRANSPOSE(dZ),dResidual)**2)) <= 1D-11)
        deallocate(dDelta,dJacobian,dMuBase,dMuMinus,dMuPlus,dMuTrial,dReduced,dResidualTrial,dXTrial)
    end subroutine SolveProductionForcedState

    subroutine SolveDense(dMatrix,dRHS,iInfoLocal)
        real(8), intent(in) :: dMatrix(:,:)
        real(8), intent(inout) :: dRHS(:,:)
        integer, intent(out) :: iInfoLocal
        integer, allocatable :: iPivot(:)
        real(8), allocatable :: dWork(:,:)
        allocate(dWork(SIZE(dMatrix,1),SIZE(dMatrix,2)),iPivot(SIZE(dMatrix,1)))
        dWork = dMatrix
        call DGESV(SIZE(dMatrix,1),SIZE(dRHS,2),dWork,SIZE(dMatrix,1),iPivot,dRHS,SIZE(dMatrix,1),iInfoLocal)
        deallocate(dWork,iPivot)
    end subroutine SolveDense

    real(8) function PositivityStep(dXBase,dDirection)
        real(8), intent(in) :: dXBase(:), dDirection(:)
        integer :: q
        PositivityStep = 5D-1
        do q = 1, SIZE(dXBase)
            if (ABS(dDirection(q)) > 0D0) &
                PositivityStep = DMIN1(PositivityStep,0.1D0*dXBase(q)/ABS(dDirection(q)))
        end do
        PositivityStep = DMAX1(PositivityStep,1D-8)
    end function PositivityStep

    real(8) function TangentNorm(dVector,dBasis)
        real(8), intent(in) :: dVector(:), dBasis(:,:)
        TangentNorm = SQRT(SUM(MATMUL(TRANSPOSE(dBasis),dVector)**2))
    end function TangentNorm

    subroutine VectorMetrics(dApprox,dReference,dNormError,dComponentError)
        real(8), intent(in) :: dApprox(:), dReference(:)
        real(8), intent(out) :: dNormError, dComponentError
        dNormError = VectorError(dApprox,dReference)
        dComponentError = MAXVAL(ABS(dApprox-dReference)/ &
            DMAX1(1D0,ABS(dApprox),ABS(dReference)))
    end subroutine VectorMetrics

    !> Scale an absolute element-space oracle uncertainty exactly like the corresponding vector error.
    real(8) function ScaledUncertainty(dAbsolute,dApprox,dReference)
        real(8), intent(in) :: dAbsolute(:), dApprox(:), dReference(:)
        ScaledUncertainty = SQRT(SUM(dAbsolute*dAbsolute))/ &
            DMAX1(1D0,SQRT(SUM(dApprox*dApprox)),SQRT(SUM(dReference*dReference)))
    end function ScaledUncertainty

    real(8) function VectorError(dApprox,dReference)
        real(8), intent(in) :: dApprox(:), dReference(:)
        VectorError = SQRT(SUM((dApprox-dReference)**2))/ &
            DMAX1(1D0,SQRT(SUM(dApprox*dApprox)),SQRT(SUM(dReference*dReference)))
    end function VectorError

    real(8) function MatrixError(dApprox,dReference)
        real(8), intent(in) :: dApprox(:,:), dReference(:,:)
        MatrixError = SQRT(SUM((dApprox-dReference)**2))/ &
            DMAX1(1D0,SQRT(SUM(dApprox*dApprox)),SQRT(SUM(dReference*dReference)))
    end function MatrixError

    subroutine BestResolved(dErrorsLocal,dUncertaintyLocal,iBest,dBest)
        real(8), intent(in) :: dErrorsLocal(:), dUncertaintyLocal(:)
        integer, intent(out) :: iBest
        real(8), intent(out) :: dBest
        integer :: k
        iBest = 0
        dBest = HUGE(1D0)
        do k = 1, SIZE(dErrorsLocal)
            if (dUncertaintyLocal(k) > 0.5D0*dErrorsLocal(k)) cycle
            if (dErrorsLocal(k) < dBest) then
                iBest = k
                dBest = dErrorsLocal(k)
            end if
        end do
    end subroutine BestResolved

    logical function HasResolvedSecondOrder(dError,dUncertainty,dOrder,lOrder)
        real(8), intent(in) :: dError(:), dUncertainty(:), dOrder(:)
        logical, intent(in) :: lOrder(:)
        integer :: k
        HasResolvedSecondOrder = .FALSE.
        do k = 1, SIZE(dError)-2
            if (.NOT. ALL(dUncertainty(k:k+2) <= 0.5D0*dError(k:k+2))) cycle
            if (.NOT. lOrder(k) .OR. .NOT. lOrder(k+1)) cycle
            if (.NOT. ((dError(k) > dError(k+1)) .AND. (dError(k+1) > dError(k+2)))) cycle
            if ((dOrder(k) < FD_ORDER_SECOND_MIN) .OR. (dOrder(k) > FD_ORDER_SECOND_MAX)) cycle
            if ((dOrder(k+1) < FD_ORDER_SECOND_MIN) .OR. (dOrder(k+1) > FD_ORDER_SECOND_MAX)) cycle
            HasResolvedSecondOrder = .TRUE.
            return
        end do
    end function HasResolvedSecondOrder

    subroutine PrintSweep(cName,dStep,dError,dUncertainty,dOrder,lOrder,iBest)
        character(len=*), intent(in) :: cName
        real(8), intent(in) :: dStep(:), dError(:), dUncertainty(:), dOrder(:)
        logical, intent(in) :: lOrder(:)
        integer, intent(in) :: iBest
        integer :: k
        write(*,'(/,A)') TRIM(cName)
        write(*,'(A)') 'step          scaled error   oracle uncertainty  resolved  observed order'
        do k = 1, SIZE(dStep)
            if (lOrder(k)) then
                write(*,'(3ES15.6,3X,L1,3X,F10.4)') dStep(k),dError(k),dUncertainty(k), &
                    dUncertainty(k) <= 0.5D0*dError(k),dOrder(k)
            else
                write(*,'(3ES15.6,3X,L1,3X,A)') dStep(k),dError(k),dUncertainty(k), &
                    dUncertainty(k) <= 0.5D0*dError(k),'N/A'
            end if
        end do
        if (iBest > 0) then
            write(*,'(A,I0,A,ES14.6)') 'best oracle-resolved step = ',iBest,', error = ',dError(iBest)
        else
            write(*,'(A)') 'best oracle-resolved step = none'
        end if
    end subroutine PrintSweep

    subroutine FinishTest(lAllPass)
        logical, intent(in) :: lAllPass
        call ResetGEMNewtonDiagnosticCapture
        call ResetThermoAll
        if (lAllPass) then
            print *, 'TestMQMQAGEMMappingVerification: PASS'
            call EXIT(0)
        else
            print *, 'TestMQMQAGEMMappingVerification: FAIL <---'
            call EXIT(1)
        end if
    end subroutine FinishTest

end program TestMQMQAGEMMappingVerification
