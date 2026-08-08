!-------------------------------------------------------------------------------------------------------------
!> \file    TestMQMQASUBQResponseVerification.F90
!> \brief   Diagnostic-only native verification of the SUBQ constrained composition response.
!>
!> \details This MQ-3B test uses the assessed FeTiVO SlagBsoln SUBQ phase. A
!!          strictly positive blend of its converged composition supplies only
!!          the starting guess; a test-only solver first converges the
!!          production partial-molar stationarity equations at fixed element
!!          potentials. The analytic Hessian and constrained response are then
!!          evaluated at that stationary root. Independently reconverged
!!          compositions at perturbed element potentials provide the nonlinear
!!          finite-difference oracle without using the analytic Hessian.
!!
!!          The test also verifies the existing ideal/simple GEM baseline using
!!          three response/block formulations plus a separate residual
!!          reconstruction. It does not assemble a mapper correction or modify
!!          GEM.
!-------------------------------------------------------------------------------------------------------------

program TestMQMQASUBQResponseVerification

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleMQMQAUnconstrained
    USE ModuleMQMQAProductionAdapter
    USE ModuleConstrainedResponse
    USE ModuleFiniteDifferenceVerification

    implicit none

    interface
        subroutine CompExcessGibbsEnergySUBG(iSolnPhaseIndex)
            integer :: iSolnPhaseIndex
        end subroutine CompExcessGibbsEnergySUBG

        subroutine CompStoichSolnPhase(iSolnPhaseIndex)
            integer :: iSolnPhaseIndex
        end subroutine CompStoichSolnPhase
    end interface

    integer, parameter :: nSteps = 8
    character(len=32) :: cArgument
    integer :: i, iDirection, iFirstLocal, iInfo, iInfoSave, iLastLocal, iPhaseIndex, iRankForcing
    integer :: iSlot, iStep, nG, nQ, nQuad, nTangent
    integer, allocatable :: iBasisColumn(:), iBestResolved(:), iWorstComponent(:,:)
    logical :: lPass, lReport, lZeroOracle
    real(8) :: dAepError, dBaselineClosedError, dBaselineGEMError, dInteriorFraction
    real(8) :: dFloorDifference, dKKTCondition, dKKTConstraintResidual, dKKTTopResidual
    real(8) :: dMaxEigenvalue, dMinEigenvalue, dNullDifference, dNullForcingError, dNullKKTResidual
    real(8) :: dNullReducedResidual, dNullTolerance, dProjectedKKTResidual
    real(8) :: dResidualError, dResponseBest, dResponseRawBest, dStabilityTolerance
    real(8) :: dStationarityBest, dWorstComponentBest, dZeroOracleDrift, dZeroStationarity
    real(8) :: dRootStationarity, dSeedMinimum, dSeedToRootDifference, dStationaryMinimum
    real(8) :: dTangentCondition
    real(8), allocatable :: dAClosed(:,:), dAClosedPhase(:,:), dAFromGEM(:,:), dAOuter(:,:)
    real(8), allocatable :: dAepDirect(:), dAepExpected(:), dConstraint(:,:)
    real(8), allocatable :: dChemicalSave(:), dEffSave(:,:), dEigenvalue(:), dErrors(:,:)
    real(8), allocatable :: dForcing(:,:), dForcingBasis(:,:), dFractionSave(:)
    real(8), allocatable :: dGradient(:), dHbase(:,:), dHessian(:,:), dHx(:,:)
    real(8), allocatable :: dKKT(:,:), dKKTSingular(:), dKtangent(:,:), dMoles(:)
    real(8), allocatable :: dMolesFloor(:), dMolesSave(:), dMuProduction(:), dNullForcing(:,:)
    real(8), allocatable :: dNullForcingMultiplier(:,:), dNullForcingResponse(:,:), dNullResponse(:,:), dOrders(:,:)
    real(8), allocatable :: dKKTMultiplier(:,:), dMaxAbsolute(:,:), dMaxScaled(:,:)
    real(8), allocatable :: dNormAbsolute(:,:), dP(:,:), dProjectedForcing(:,:)
    real(8), allocatable :: dPartialSave(:), dRawMinimum(:), dResolvedBest(:), dResponse(:,:), dResponseBase(:,:)
    real(8), allocatable :: dResponseClosed(:,:), dSteps(:,:), dStationarity(:,:), dX(:), dXConverged(:)
    real(8), allocatable :: dXSeed(:), dZ(:,:), dZeroResidual(:), dZeroX(:)
    logical, allocatable :: lOracleResolved(:,:), lOrderAvailable(:,:)
    type(FDSweepAssessment) :: tSweep
    type(MQMQAModelData) :: tModel
    type(MQMQAInteractionTerm), allocatable :: tInteraction(:)

    lPass = .TRUE.
    lReport = .FALSE.
    if (COMMAND_ARGUMENT_COUNT() > 0) then
        call GET_COMMAND_ARGUMENT(1,cArgument)
        lReport = TRIM(cArgument) == '--report'
    end if

    !=========================================================================================================
    ! SECTION 1: ASSESSED SUBQ PHASE AND POSITIVE STATIONARITY SEED
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
    if (iPhaseIndex <= 0) call FinishTest(lPass)
    lPass = lPass .AND. (iPhaseElectronID(iPhaseIndex) == 0)

    iFirstLocal = nSpeciesPhase(iPhaseIndex-1)+1
    iLastLocal = nSpeciesPhase(iPhaseIndex)
    nQuad = iLastLocal-iFirstLocal+1
    nTangent = nQuad-1
    allocate(dX(nQuad),dXConverged(nQuad),dXSeed(nQuad),dMoles(nQuad),dMolesFloor(nQuad), &
        dGradient(nQuad),dMuProduction(nQuad), &
        dHessian(nQuad,nQuad),dHx(nQuad,nQuad),dHbase(nQuad,nQuad), &
        dConstraint(1,nQuad),dForcing(nQuad,nElements),dZ(nQuad,nTangent), &
        dKtangent(nTangent,nTangent),dEigenvalue(nTangent),iBasisColumn(nElements))

    lPass = lPass .AND. (nQuad == 15)
    dXConverged = dMolFraction(iFirstLocal:iLastLocal)
    dXConverged = dXConverged/SUM(dXConverged)
    dInteriorFraction = 0.20D0
    dXSeed = (1D0-dInteriorFraction)*dXConverged+dInteriorFraction/DFLOAT(nQuad)
    dSeedMinimum = MINVAL(dXSeed)
    lPass = lPass .AND. (dSeedMinimum > 1D-12) .AND. (DABS(SUM(dXSeed)-1D0) <= 1D-12)

    call DecodeProductionSUBQPhase(iPhaseIndex,tModel,tInteraction,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    if (iInfo /= 0) call FinishTest(lPass)
    nG = COUNT(tInteraction%iFamily == MQMQA_TERM_G)
    nQ = COUNT(tInteraction%iFamily == MQMQA_TERM_Q)
    lPass = lPass .AND. (nG == 6) .AND. (nQ == 8)
    lPass = lPass .AND. (MAXVAL(DABS(tModel%dZeta-2.4D0)) <= 1D-12)

    dConstraint = 1D0
    do i = 1, nElements
        dForcing(:,i) = dStoichSpecies(iFirstLocal:iLastLocal,i) / &
            DFLOAT(iParticlesPerMole(iFirstLocal:iLastLocal))
    end do

    ! Helmert contrasts provide an exactly structured orthonormal basis for all
    ! normalized composition changes. Column i moves the first i quadruplets
    ! together against quadruplet i+1 while preserving their total.
    dZ = 0D0
    do i = 1, nTangent
        dZ(1:i,i) = 1D0/DSQRT(DFLOAT(i*(i+1)))
        dZ(i+1,i) = -DFLOAT(i)/DSQRT(DFLOAT(i*(i+1)))
    end do

    ! The blended composition is intentionally only a safe starting guess. A
    ! response derivative is meaningful at the stationary root associated with
    ! the retained production element potentials, not at an arbitrary interior
    ! point chosen for the earlier Hessian test.
    allocate(dZeroX(nQuad),dZeroResidual(nQuad))
    call EvaluateProductionOracleState(iPhaseIndex,iFirstLocal,iLastLocal,dXSeed, &
        dElementPotential*0D0,0D0,dZeroX,dZeroResidual,lZeroOracle)
    lPass = lPass .AND. lZeroOracle
    dSeedToRootDifference = SQRT(SUM((dZeroX-dXSeed)**2))
    dRootStationarity = SQRT(SUM(MATMUL(TRANSPOSE(dZ),dZeroResidual)**2))
    dStationaryMinimum = MINVAL(dZeroX)
    lPass = lPass .AND. (dRootStationarity <= 1D-11) .AND. (dStationaryMinimum > 1D-12)

    dX = dZeroX
    dMoles = dMolesPhase(iSlot)*dX
    call CompMQMQAHessianUnconstrained(tModel,dMoles,1D0,tInteraction,dHessian,iInfo, &
        dGradient=dGradient)
    lPass = lPass .AND. (iInfo == 0)
    dHx = dMolesPhase(iSlot)*dHessian
    dKtangent = MATMUL(TRANSPOSE(dZ),MATMUL(dHx,dZ))
    call SymmetricEigenvalues(dKtangent,dEigenvalue,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    dMinEigenvalue = MINVAL(dEigenvalue)
    dMaxEigenvalue = MAXVAL(dEigenvalue)
    dTangentCondition = dMaxEigenvalue/dMinEigenvalue
    dStabilityTolerance = 1D-12*DMAX1(1D0,dMaxEigenvalue)
    lPass = lPass .AND. (dMinEigenvalue > dStabilityTolerance)

    !=========================================================================================================
    ! SECTION 2: INDEPENDENT ELEMENT-POTENTIAL FORCING AND RESPONSE SOLVES
    !=========================================================================================================
    allocate(dProjectedForcing(nTangent,nElements))
    dProjectedForcing = MATMUL(TRANSPOSE(dZ),dForcing)
    call SelectIndependentColumnsPivotedQR(dProjectedForcing,1D-11,iBasisColumn,iRankForcing)
    lPass = lPass .AND. (iRankForcing > 0)
    allocate(dP(nElements,iRankForcing),dForcingBasis(nQuad,iRankForcing), &
        dResponse(nQuad,iRankForcing),dResponseBase(nQuad,iRankForcing), &
        dResponseClosed(nQuad,iRankForcing),dNullResponse(nQuad,iRankForcing), &
        dKKTMultiplier(1,iRankForcing))
    dP = 0D0
    do i = 1, iRankForcing
        dP(iBasisColumn(i),i) = 1D0
    end do
    dForcingBasis = MATMUL(dForcing,dP)

    call SolveConstrainedResponse(dHx,dConstraint,dForcingBasis,dResponse,iInfo,dKKTMultiplier)
    lPass = lPass .AND. (iInfo == 0)
    dKKTConstraintResidual = MAXVAL(DABS(MATMUL(dConstraint,dResponse)))
    dKKTTopResidual = MAXVAL(DABS(MATMUL(dHx,dResponse)+ &
        MATMUL(TRANSPOSE(dConstraint),dKKTMultiplier)-dForcingBasis))
    dProjectedKKTResidual = ProjectedKKTResidual(dHx,dConstraint,dForcingBasis,dResponse)
    lPass = lPass .AND. (dKKTConstraintResidual <= 1D-10) .AND. &
        (dKKTTopResidual <= 1D-9) .AND. (dProjectedKKTResidual <= 1D-9)

    call SolveNullSpace(dKtangent,dZ,dForcingBasis,dNullResponse,iInfo)
    dNullDifference = MAXVAL(DABS(dResponse-dNullResponse)) / &
        DMAX1(1D0,MAXVAL(DABS(dResponse)),MAXVAL(DABS(dNullResponse)))
    dNullKKTResidual = ProjectedKKTResidual(dHx,dConstraint,dForcingBasis,dNullResponse)
    dNullKKTResidual = dNullKKTResidual/DMAX1(1D0,MAXVAL(DABS(dForcingBasis)), &
        MAXVAL(DABS(MATMUL(dHx,dNullResponse))))
    dNullReducedResidual = MAXVAL(DABS(MATMUL(TRANSPOSE(dZ),dForcingBasis)- &
        MATMUL(dKtangent,MATMUL(TRANSPOSE(dZ),dNullResponse))))
    dNullReducedResidual = dNullReducedResidual/DMAX1(1D0, &
        MAXVAL(DABS(MATMUL(TRANSPOSE(dZ),dForcingBasis))), &
        MAXVAL(DABS(dKtangent))*MAXVAL(DABS(MATMUL(TRANSPOSE(dZ),dNullResponse))))
    ! The tangent system is strongly ill-conditioned. Its condition-number
    ! estimate supplies the standard first-order bound on disagreement between
    ! independently solved forward solutions. Both equation residuals are
    ! scale-normalized because their large intermediate terms cancel on the
    ! tangent space; this distinguishes backward accuracy from forward error.
    dNullTolerance = DMAX1(1D-9,EPSILON(1D0)*dTangentCondition)
    lPass = lPass .AND. (iInfo == 0) .AND. (dNullDifference <= dNullTolerance) .AND. &
        (dNullReducedResidual <= 1D-10) .AND. (dNullKKTResidual <= dNullTolerance)

    ! A forcing entirely in range(C^T) changes only the normalization
    ! multiplier. It must not produce a physical composition response.
    allocate(dNullForcing(nQuad,1),dNullForcingResponse(nQuad,1),dNullForcingMultiplier(1,1))
    dNullForcing = 1D0
    call SolveConstrainedResponse(dHx,dConstraint,dNullForcing,dNullForcingResponse,iInfo, &
        dNullForcingMultiplier)
    dNullForcingError = MAXVAL(DABS(dNullForcingResponse))
    lPass = lPass .AND. (iInfo == 0) .AND. (dNullForcingError <= 1D-10)

    allocate(dKKT(nQuad+1,nQuad+1),dKKTSingular(nQuad+1))
    dKKT = 0D0
    dKKT(1:nQuad,1:nQuad) = dHx
    dKKT(1:nQuad,nQuad+1) = 1D0
    dKKT(nQuad+1,1:nQuad) = 1D0
    call SingularValues(dKKT,dKKTSingular,iInfo)
    lPass = lPass .AND. (iInfo == 0) .AND. (MINVAL(dKKTSingular) > 0D0)
    dKKTCondition = MAXVAL(dKKTSingular)/MINVAL(dKKTSingular)

    !=========================================================================================================
    ! SECTION 3: PHASE-LOCAL IDEAL GEM BASELINE CROSS-CHECKS
    !=========================================================================================================
    dHbase = 0D0
    do i = 1, nQuad
        dHbase(i,i) = 1D0/dX(i)
    end do
    call SolveConstrainedResponse(dHbase,dConstraint,dForcingBasis,dResponseBase,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    dResponseClosed = MATMUL(DiagonalMatrix(dX)-OuterProduct(dX,dX),dForcingBasis)
    dBaselineClosedError = MAXVAL(DABS(dResponseBase-dResponseClosed)) / &
        DMAX1(1D0,MAXVAL(DABS(dResponseBase)),MAXVAL(DABS(dResponseClosed)))

    allocate(dAFromGEM(nElements,nElements),dAOuter(nElements,nElements), &
        dAClosed(nElements,nElements),dAClosedPhase(nElements,nElements), &
        dAepDirect(nElements),dAepExpected(nElements))

    ! Reconstruct the production formulas at the same stationary state used by
    ! the response check. Save and restore every global array installed for this
    ! source-faithful phase-local comparison.
    allocate(dChemicalSave(SIZE(dChemicalPotential)),dPartialSave(SIZE(dPartialExcessGibbs)), &
        dFractionSave(SIZE(dMolFraction)),dMolesSave(SIZE(dMolesSpecies)), &
        dEffSave(SIZE(dEffStoichSolnPhase,1),SIZE(dEffStoichSolnPhase,2)))
    dChemicalSave = dChemicalPotential
    dPartialSave = dPartialExcessGibbs
    dFractionSave = dMolFraction
    dMolesSave = dMolesSpecies
    dEffSave = dEffStoichSolnPhase
    iInfoSave = INFOThermo
    dMolFraction(iFirstLocal:iLastLocal) = dX
    dMolesSpecies(iFirstLocal:iLastLocal) = MAX(dMolesPhase(iSlot)*dX,dTolerance(8))
    dMolesFloor = dMolesSpecies(iFirstLocal:iLastLocal)
    dFloorDifference = MAXVAL(DABS(dMolesFloor-dMoles))/ &
        DMAX1(1D0,MAXVAL(DABS(dMoles)))
    call EvaluateProductionMu(iPhaseIndex,iFirstLocal,iLastLocal,dX,dMuProduction,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    call ReconstructPhaseLocalGEMBaseline(iSlot,dX,dMuProduction,dForcing, &
        dAFromGEM,dAOuter,dAepDirect,dResidualError)
    call CompStoichSolnPhase(iPhaseIndex)
    dAClosedPhase = dAFromGEM-dAOuter
    dAClosed = dMolesPhase(iSlot)*MATMUL(TRANSPOSE(dForcing), &
        MATMUL(DiagonalMatrix(dX)-OuterProduct(dX,dX),dForcing))
    dBaselineGEMError = MAXVAL(DABS(dAClosed-dAClosedPhase)) / &
        DMAX1(1D0,MAXVAL(DABS(dAClosed)),MAXVAL(DABS(dAClosedPhase)))
    dAepExpected = dMolesPhase(iSlot)*MATMUL(TRANSPOSE(dForcing),dX)
    dAepError = DMAX1(MAXVAL(DABS(dAepExpected-dAepDirect)), &
        MAXVAL(DABS(dAepExpected-dMolesPhase(iSlot)*dEffStoichSolnPhase(iPhaseIndex,:)))) / &
        DMAX1(1D0,MAXVAL(DABS(dAepExpected)),MAXVAL(DABS(dAepDirect)))
    lPass = lPass .AND. (dBaselineClosedError <= 1D-12) .AND. &
        (dBaselineGEMError <= 1D-10) .AND. (dAepError <= 1D-10) .AND. &
        (dResidualError <= 1D-10)
    dChemicalPotential = dChemicalSave
    dPartialExcessGibbs = dPartialSave
    dMolFraction = dFractionSave
    dMolesSpecies = dMolesSave
    dEffStoichSolnPhase = dEffSave
    INFOThermo = iInfoSave

    ! Re-solving the unperturbed stationary equations from x0 must return the
    ! same root before any finite-difference response is accepted.
    call EvaluateProductionOracleState(iPhaseIndex,iFirstLocal,iLastLocal,dX, &
        dElementPotential*0D0,0D0,dZeroX,dZeroResidual,lZeroOracle)
    lPass = lPass .AND. lZeroOracle
    dZeroOracleDrift = SQRT(SUM((dZeroX-dX)**2))
    dZeroStationarity = SQRT(SUM(MATMUL(TRANSPOSE(dZ),dZeroResidual)**2))
    lPass = lPass .AND. (dZeroOracleDrift <= 1D-8) .AND. (dZeroStationarity <= 1D-11)

    !=========================================================================================================
    ! SECTION 4: PRODUCTION PARTIAL-MOLAR STATIONARITY ORACLE
    !=========================================================================================================
    allocate(dSteps(iRankForcing,nSteps),dNormAbsolute(iRankForcing,nSteps), &
        dErrors(iRankForcing,nSteps),dMaxAbsolute(iRankForcing,nSteps), &
        dMaxScaled(iRankForcing,nSteps),iWorstComponent(iRankForcing,nSteps), &
        dStationarity(iRankForcing,nSteps),dOrders(iRankForcing,nSteps), &
        lOracleResolved(iRankForcing,nSteps),lOrderAvailable(iRankForcing,nSteps), &
        dRawMinimum(iRankForcing),dResolvedBest(iRankForcing),iBestResolved(iRankForcing))
    dResponseBest = 0D0
    dResponseRawBest = 0D0
    dStationarityBest = 0D0
    dWorstComponentBest = 0D0
    do iDirection = 1, iRankForcing
        call VerifyNativeResponseDirection(iPhaseIndex,iFirstLocal,iLastLocal,dX,dP(:,iDirection), &
            dResponse(:,iDirection),dZ,dSteps(iDirection,:),dErrors(iDirection,:), &
            dNormAbsolute(iDirection,:),dMaxAbsolute(iDirection,:),dMaxScaled(iDirection,:), &
            iWorstComponent(iDirection,:),dStationarity(iDirection,:),lPass)
        call AssessFDSweep(dSteps(iDirection,:),dErrors(iDirection,:), &
            FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX,1D-6,tSweep, &
            dOrders(iDirection,:),lOrderAvailable(iDirection,:))
        lOracleResolved(iDirection,:) = dStationarity(iDirection,:) <= &
            0.5D0*dErrors(iDirection,:)
        dRawMinimum(iDirection) = MINVAL(dErrors(iDirection,:))
        call FindBestResolvedPoint(dErrors(iDirection,:),lOracleResolved(iDirection,:), &
            iBestResolved(iDirection),dResolvedBest(iDirection))
        lPass = lPass .AND. (iBestResolved(iDirection) > 0)
        lPass = lPass .AND. HasResolvedSecondOrderRegion(dErrors(iDirection,:), &
            dOrders(iDirection,:),lOrderAvailable(iDirection,:),lOracleResolved(iDirection,:))
        if (iBestResolved(iDirection) > 0) then
            iStep = iBestResolved(iDirection)
            lPass = lPass .AND. (dResolvedBest(iDirection) <= 1D-6) .AND. &
                (dMaxScaled(iDirection,iStep) <= 1D-6)
            dResponseBest = DMAX1(dResponseBest,dResolvedBest(iDirection))
            dStationarityBest = DMAX1(dStationarityBest,dStationarity(iDirection,iStep))
            dWorstComponentBest = DMAX1(dWorstComponentBest,dMaxScaled(iDirection,iStep))
        end if
        dResponseRawBest = DMAX1(dResponseRawBest,dRawMinimum(iDirection))
    end do

    if (lReport) then
        write(*,'(A)') 'MQ-3B native SUBQ constrained-response verification'
        write(*,'(A)') 'scope: assessed FeTiVO SUBQ G/Q local response and ideal GEM baseline only'
        write(*,'(A)') 'excludes reduced deltaA/deltaB construction and GEM activation'
        write(*,'(A,A)') 'phase = ',TRIM(cSolnPhaseName(iPhaseIndex))
        write(*,'(A,I0)') 'quadruplet count = ',nQuad
        write(*,'(A,I0)') 'decoded G parameter count = ',nG
        write(*,'(A,I0)') 'decoded Q parameter count = ',nQ
        write(*,'(A,ES14.6)') 'minimum pair zeta = ',MINVAL(tModel%dZeta)
        write(*,'(A,ES14.6)') 'maximum pair zeta = ',MAXVAL(tModel%dZeta)
        write(*,'(A,ES14.6)') 'positive seed minimum fraction = ',dSeedMinimum
        write(*,'(A,ES14.6)') 'seed-to-stationary-root composition change = ',dSeedToRootDifference
        write(*,'(A,ES14.6)') 'stationary-root minimum fraction = ',dStationaryMinimum
        write(*,'(A,ES14.6)') 'stationary-root tangent residual = ',dRootStationarity
        write(*,'(A,I0)') 'constraint rank = ',1
        write(*,'(A,I0)') 'supported forcing rank = ',iRankForcing
        write(*,'(A,I0)') 'null/dependent original forcing columns = ',nElements-iRankForcing
        write(*,'(A,*(I0,1X))') 'forcing basis element columns = ',iBasisColumn(1:iRankForcing)
        write(*,'(A)',ADVANCE='NO') 'forcing basis elements = '
        do i = 1, iRankForcing
            write(*,'(I0,A,A,A)',ADVANCE='NO') iBasisColumn(i),' (', &
                TRIM(cElementName(iBasisColumn(i))),')  '
        end do
        write(*,*)
        write(*,'(A,ES14.6)') 'authoritative/floored mole discrepancy = ',dFloorDifference
        write(*,'(A,ES14.6)') 'minimum tangent eigenvalue = ',dMinEigenvalue
        write(*,'(A,ES14.6)') 'maximum tangent eigenvalue = ',dMaxEigenvalue
        write(*,'(A,ES14.6)') 'tangent condition estimate = ',dTangentCondition
        write(*,'(A,ES14.6)') 'positive-definite tolerance = ',dStabilityTolerance
        write(*,'(A,ES14.6)') 'KKT minimum singular value = ',MINVAL(dKKTSingular)
        write(*,'(A,ES14.6)') 'KKT maximum singular value = ',MAXVAL(dKKTSingular)
        write(*,'(A,ES14.6)') 'KKT condition estimate = ',dKKTCondition
        write(*,'(A,ES14.6)') 'KKT top-equation residual = ',dKKTTopResidual
        write(*,'(A,ES14.6)') 'KKT constraint residual = ',dKKTConstraintResidual
        write(*,'(A,ES14.6)') 'projected KKT residual cross-check = ',dProjectedKKTResidual
        write(*,'(A,ES14.6)') 'null-space response difference = ',dNullDifference
        write(*,'(A,ES14.6)') 'condition-aware null-space difference tolerance = ',dNullTolerance
        write(*,'(A,ES14.6)') 'null-space scaled reduced-equation residual = ',dNullReducedResidual
        write(*,'(A,ES14.6)') 'null-space scaled projected equation residual = ',dNullKKTResidual
        write(*,'(A,ES14.6)') 'normalization-only forcing response = ',dNullForcingError
        write(*,'(A,ES14.6)') 'ideal solve/closed response difference = ',dBaselineClosedError
        write(*,'(A,ES14.6)') 'source-faithful phase-local centered-block reconstruction difference = ', &
            dBaselineGEMError
        write(*,'(A,ES14.6)') 'source-faithful phase-local element-phase reconstruction difference = ', &
            dAepError
        write(*,'(A,ES14.6)') 'source-faithful phase-local residual reconstruction difference = ', &
            dResidualError
        write(*,'(A,ES14.6)') 'zero-perturbation production-oracle composition drift = ',dZeroOracleDrift
        write(*,'(A,ES14.6)') 'zero-perturbation tangent stationarity residual = ',dZeroStationarity
        write(*,'(A,*(I0,1X))') 'particles per mole = ',iParticlesPerMole(iFirstLocal:iLastLocal)
        write(*,'(A)') 'dir  h            norm abs       norm scaled    max abs        max scaled' // &
            '     worst  uncertainty    resolved  order'
        do iDirection = 1, iRankForcing
            do iStep = 1, nSteps
                write(*,'(I3,5ES15.6,I7,ES15.6,3X,L1,3X,A)') iDirection,dSteps(iDirection,iStep), &
                    dNormAbsolute(iDirection,iStep),dErrors(iDirection,iStep), &
                    dMaxAbsolute(iDirection,iStep),dMaxScaled(iDirection,iStep), &
                    iWorstComponent(iDirection,iStep),dStationarity(iDirection,iStep), &
                    lOracleResolved(iDirection,iStep), &
                    TRIM(OrderLabel(dOrders(iDirection,iStep),lOrderAvailable(iDirection,iStep)))
            end do
            write(*,'(A,I0,A,ES14.6)') 'direction ',iDirection,' raw minimum error = ', &
                dRawMinimum(iDirection)
            if (iBestResolved(iDirection) > 0) then
                write(*,'(A,I0,A,I0,A,ES14.6,A,I0)') 'direction ',iDirection, &
                    ' best oracle-resolved step = ',iBestResolved(iDirection),', error = ', &
                    dResolvedBest(iDirection),', worst quadruplet = ', &
                    iWorstComponent(iDirection,iBestResolved(iDirection))
            else
                write(*,'(A,I0,A)') 'direction ',iDirection,' has no oracle-resolved point'
            end if
        end do
        write(*,'(A,ES14.6)') 'worst raw minimum native response error = ',dResponseRawBest
        write(*,'(A,ES14.6)') 'worst best oracle-resolved response error = ',dResponseBest
        write(*,'(A,ES14.6)') 'worst best resolved component scaled error = ',dWorstComponentBest
        write(*,'(A,ES14.6)') 'worst best resolved stationarity uncertainty = ',dStationarityBest
    end if

    call ResetThermoAll
    call FinishTest(lPass)

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Compare one analytic forcing response with reconverged production compositions.
    !---------------------------------------------------------------------------------------------------------
    subroutine VerifyNativeResponseDirection(iPhase,iFirst,iLast,dXBase,dGammaDirection,dPrediction, &
        dZLocal,dStepValues,dErrorValues,dNormAbsoluteValues, &
        dMaxAbsoluteValues,dMaxScaledValues,iWorstValues,dStationarityValues,lAllPass)

        integer, intent(in) :: iPhase, iFirst, iLast
        real(8), intent(in) :: dXBase(:), dGammaDirection(:), dPrediction(:)
        real(8), intent(in) :: dZLocal(:,:)
        real(8), intent(out) :: dStepValues(:), dErrorValues(:), dNormAbsoluteValues(:)
        real(8), intent(out) :: dMaxAbsoluteValues(:), dMaxScaledValues(:), dStationarityValues(:)
        integer, intent(out) :: iWorstValues(:)
        logical, intent(inout) :: lAllPass

        integer :: iStepLocal
        logical :: lMinus, lPlus
        real(8) :: dH, dHMaximum
        real(8), allocatable :: dFD(:), dJacobianMinus(:,:), dJacobianPlus(:,:)
        real(8), allocatable :: dResidualMinus(:), dResidualPlus(:), dXMinus(:), dXPlus(:)

        allocate(dFD(SIZE(dXBase)),dResidualMinus(SIZE(dXBase)),dResidualPlus(SIZE(dXBase)), &
            dXMinus(SIZE(dXBase)),dXPlus(SIZE(dXBase)), &
            dJacobianMinus(SIZE(dZLocal,2),SIZE(dZLocal,2)), &
            dJacobianPlus(SIZE(dZLocal,2),SIZE(dZLocal,2)))
        ! Begin with a nonlinear perturbation large enough to expose the
        ! expected truncation region, then refine toward finite precision.
        dHMaximum = 5D-1
        do iStepLocal = 1, SIZE(dXBase)
            if (DABS(dPrediction(iStepLocal)) > 0D0) then
                dHMaximum = DMIN1(dHMaximum,0.1D0*dXBase(iStepLocal)/DABS(dPrediction(iStepLocal)))
            end if
        end do
        dHMaximum = DMAX1(dHMaximum,1D-10)

        do iStepLocal = 1, SIZE(dStepValues)
            dH = dHMaximum*3D0**(-(iStepLocal-1))
            dStepValues(iStepLocal) = dH
            call EvaluateProductionOracleState(iPhase,iFirst,iLast,dXBase,dGammaDirection,-dH, &
                dXMinus,dResidualMinus,lMinus,dJacobianMinus)
            call EvaluateProductionOracleState(iPhase,iFirst,iLast,dXBase,dGammaDirection,dH, &
                dXPlus,dResidualPlus,lPlus,dJacobianPlus)
            lAllPass = lAllPass .AND. lMinus .AND. lPlus
            dFD = (dXPlus-dXMinus)/(2D0*dH)
            call ComputeVectorErrorMetrics(dFD,dPrediction,dNormAbsoluteValues(iStepLocal), &
                dErrorValues(iStepLocal),dMaxAbsoluteValues(iStepLocal), &
                dMaxScaledValues(iStepLocal),iWorstValues(iStepLocal))
            call EstimateStationarityUncertainty(dResidualMinus,dResidualPlus,dH,dZLocal, &
                dJacobianMinus,dJacobianPlus,dStationarityValues(iStepLocal),lAllPass)
        end do

        deallocate(dFD,dJacobianMinus,dJacobianPlus,dResidualMinus,dResidualPlus,dXMinus,dXPlus)

    end subroutine VerifyNativeResponseDirection


    !---------------------------------------------------------------------------------------------------------
    !> \brief Converge the production partial-molar equations at one perturbed element-potential state.
    !>
    !> \details Every modified global array is restored before returning. The
    !!          returned residual is the nonconstant part of mu-S*Gamma; a
    !!          constant component is the normalization multiplier and is not
    !!          a stationarity defect.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateProductionOracleState(iPhase,iFirst,iLast,dXBase,dGammaDirection,dH, &
        dXResult,dResidual,lSuccess,dTangentJacobian)

        integer, intent(in) :: iPhase, iFirst, iLast
        real(8), intent(in) :: dXBase(:), dGammaDirection(:), dH
        real(8), intent(out) :: dXResult(:), dResidual(:)
        logical, intent(out) :: lSuccess
        real(8), intent(out), optional :: dTangentJacobian(:,:)

        integer :: iInfoSave
        logical :: lOracleConverged
        real(8), allocatable :: dChemicalSave(:), dDrivingSave(:), dElementSave(:)
        real(8), allocatable :: dFractionSave(:), dPartialSave(:)

        allocate(dChemicalSave(SIZE(dChemicalPotential)),dDrivingSave(SIZE(dDrivingForceSoln)), &
            dElementSave(SIZE(dElementPotential)),dFractionSave(SIZE(dMolFraction)), &
            dPartialSave(SIZE(dPartialExcessGibbs)))
        dChemicalSave = dChemicalPotential
        dDrivingSave = dDrivingForceSoln
        dElementSave = dElementPotential
        dFractionSave = dMolFraction
        dPartialSave = dPartialExcessGibbs
        iInfoSave = INFOThermo

        dElementPotential = dElementSave+dH*dGammaDirection
        dMolFraction(iFirst:iLast) = dXBase
        call SolveProductionStationarity(iPhase,iFirst,iLast,dXBase,dElementPotential, &
            dXResult,dResidual,lOracleConverged,dTangentJacobian)
        lSuccess = (INFOThermo == iInfoSave) .AND. lOracleConverged .AND. &
            ALL(dXResult > 1D-12) .AND. ALL(IEEE_IS_FINITE(dXResult)) .AND. &
            ALL(IEEE_IS_FINITE(dResidual))

        dChemicalPotential = dChemicalSave
        dDrivingForceSoln = dDrivingSave
        dElementPotential = dElementSave
        dMolFraction = dFractionSave
        dPartialExcessGibbs = dPartialSave
        INFOThermo = iInfoSave
        deallocate(dChemicalSave,dDrivingSave,dElementSave,dFractionSave,dPartialSave)

    end subroutine EvaluateProductionOracleState


    !---------------------------------------------------------------------------------------------------------
    !> \brief Solve the established production partial-molar stationarity equations.
    !>
    !> \details The analytic MQMQA Hessian is deliberately not used. Each
    !!          Jacobian column is instead formed by central differences of
    !!          CompExcessGibbsEnergySUBG partial molars. A positivity-preserving
    !!          line search then drives the tangent stationarity residual to a
    !!          tolerance suitable for differentiating the converged response.
    !---------------------------------------------------------------------------------------------------------
    subroutine SolveProductionStationarity(iPhase,iFirst,iLast,dXStart,dGammaLocal, &
        dXResult,dResidualResult,lConverged,dTangentJacobian)

        integer, intent(in) :: iPhase, iFirst, iLast
        real(8), intent(in) :: dXStart(:), dGammaLocal(:)
        real(8), intent(out) :: dXResult(:), dResidualResult(:)
        logical, intent(out) :: lConverged
        real(8), intent(out), optional :: dTangentJacobian(:,:)

        integer :: iDirectionLocal, iInfoLocal, iIteration, iTrial
        real(8) :: dAlpha, dResidualNorm, dTrialNorm
        real(8), allocatable :: dDelta(:), dJacobian(:,:)
        real(8), allocatable :: dMuTrial(:), dReducedResidual(:,:), dResidualTrial(:)
        real(8), allocatable :: dXTrial(:)

        allocate(dDelta(SIZE(dXStart)),dJacobian(SIZE(dZ,2),SIZE(dZ,2)), &
            dMuTrial(SIZE(dXStart)), &
            dReducedResidual(SIZE(dZ,2),1),dResidualTrial(SIZE(dXStart)), &
            dXTrial(SIZE(dXStart)))
        dXResult = dXStart/SUM(dXStart)
        lConverged = .FALSE.
        do iIteration = 1, 40
            call EvaluateProductionStationarity(iPhase,iFirst,iLast,dXResult,dGammaLocal, &
                dMuTrial,dResidualResult,iInfoLocal)
            if (iInfoLocal /= 0) exit
            dReducedResidual(:,1) = MATMUL(TRANSPOSE(dZ),dResidualResult)
            dResidualNorm = SQRT(SUM(dReducedResidual(:,1)**2))
            if (dResidualNorm <= 1D-12) then
                lConverged = .TRUE.
                exit
            end if

            call BuildProductionTangentJacobian(iPhase,iFirst,iLast,dXResult,dJacobian,iInfoLocal)
            if (iInfoLocal /= 0) exit

            dReducedResidual = -dReducedResidual
            call SolveDense(dJacobian,dReducedResidual,iInfoLocal)
            if (iInfoLocal /= 0) exit
            dDelta = MATMUL(dZ,dReducedResidual(:,1))
            dAlpha = 1D0
            do iDirectionLocal = 1, SIZE(dXResult)
                if (dDelta(iDirectionLocal) < 0D0) then
                    dAlpha = DMIN1(dAlpha,-0.9D0*dXResult(iDirectionLocal)/dDelta(iDirectionLocal))
                end if
            end do

            do iTrial = 1, 20
                dXTrial = dXResult+dAlpha*dDelta
                call EvaluateProductionStationarity(iPhase,iFirst,iLast,dXTrial,dGammaLocal, &
                    dMuTrial,dResidualTrial,iInfoLocal)
                if (iInfoLocal == 0) then
                    dTrialNorm = SQRT(SUM(MATMUL(TRANSPOSE(dZ),dResidualTrial)**2))
                    if (dTrialNorm < dResidualNorm) exit
                end if
                dAlpha = 0.5D0*dAlpha
            end do
            if ((iInfoLocal /= 0) .OR. (dTrialNorm >= dResidualNorm)) exit
            dXResult = dXTrial
        end do

        call EvaluateProductionStationarity(iPhase,iFirst,iLast,dXResult,dGammaLocal, &
            dMuTrial,dResidualResult,iInfoLocal)
        if (iInfoLocal == 0) then
            dResidualNorm = SQRT(SUM(MATMUL(TRANSPOSE(dZ),dResidualResult)**2))
            lConverged = dResidualNorm <= 1D-11
        end if
        if (PRESENT(dTangentJacobian)) then
            if ((SIZE(dTangentJacobian,1) /= SIZE(dZ,2)) .OR. &
                (SIZE(dTangentJacobian,2) /= SIZE(dZ,2))) then
                lConverged = .FALSE.
            else if (iInfoLocal == 0) then
                call BuildProductionTangentJacobian(iPhase,iFirst,iLast,dXResult, &
                    dTangentJacobian,iInfoLocal)
                lConverged = lConverged .AND. (iInfoLocal == 0)
            end if
        end if
        deallocate(dDelta,dJacobian,dMuTrial,dReducedResidual,dResidualTrial,dXTrial)

    end subroutine SolveProductionStationarity


    !---------------------------------------------------------------------------------------------------------
    !> \brief Build the oracle's tangent Jacobian from production partial molars.
    !>
    !> \details This matrix is independent of the analytic SUBQ Hessian. Each
    !!          column differentiates production partial molars along one
    !!          normalized-composition tangent direction at the supplied state.
    !---------------------------------------------------------------------------------------------------------
    subroutine BuildProductionTangentJacobian(iPhase,iFirst,iLast,dXLocal,dJacobianLocal,iInfoLocal)

        integer, intent(in) :: iPhase, iFirst, iLast
        real(8), intent(in) :: dXLocal(:)
        real(8), intent(out) :: dJacobianLocal(:,:)
        integer, intent(out) :: iInfoLocal

        integer :: iDirectionLocal
        real(8) :: dHJacobian
        real(8), allocatable :: dDistance(:), dMuMinus(:), dMuPlus(:), dXMinus(:), dXPlus(:)

        iInfoLocal = 0
        dJacobianLocal = 0D0
        allocate(dDistance(SIZE(dXLocal)),dMuMinus(SIZE(dXLocal)),dMuPlus(SIZE(dXLocal)), &
            dXMinus(SIZE(dXLocal)),dXPlus(SIZE(dXLocal)))
        do iDirectionLocal = 1, SIZE(dZ,2)
            where (DABS(dZ(:,iDirectionLocal)) > 0D0)
                dDistance = dXLocal/DABS(dZ(:,iDirectionLocal))
            elsewhere
                dDistance = HUGE(1D0)
            end where
            dHJacobian = DMIN1(1D-5,0.1D0*MINVAL(dDistance))
            dXMinus = dXLocal-dHJacobian*dZ(:,iDirectionLocal)
            dXPlus = dXLocal+dHJacobian*dZ(:,iDirectionLocal)
            call EvaluateProductionMu(iPhase,iFirst,iLast,dXMinus,dMuMinus,iInfoLocal)
            if (iInfoLocal /= 0) exit
            call EvaluateProductionMu(iPhase,iFirst,iLast,dXPlus,dMuPlus,iInfoLocal)
            if (iInfoLocal /= 0) exit
            dJacobianLocal(:,iDirectionLocal) = MATMUL(TRANSPOSE(dZ), &
                (dMuPlus-dMuMinus))/(2D0*dHJacobian)
        end do
        deallocate(dDistance,dMuMinus,dMuPlus,dXMinus,dXPlus)

    end subroutine BuildProductionTangentJacobian


    subroutine EvaluateProductionStationarity(iPhase,iFirst,iLast,dXLocal,dGammaLocal, &
        dMuLocal,dResidualLocal,iInfoLocal)

        integer, intent(in) :: iPhase, iFirst, iLast
        real(8), intent(in) :: dXLocal(:), dGammaLocal(:)
        real(8), intent(out) :: dMuLocal(:), dResidualLocal(:)
        integer, intent(out) :: iInfoLocal
        real(8) :: dMean

        call EvaluateProductionMu(iPhase,iFirst,iLast,dXLocal,dMuLocal,iInfoLocal)
        if (iInfoLocal /= 0) return
        dResidualLocal = dMuLocal-MATMUL(dForcing,dGammaLocal)
        dMean = SUM(dResidualLocal)/DFLOAT(SIZE(dResidualLocal))
        dResidualLocal = dResidualLocal-dMean

    end subroutine EvaluateProductionStationarity


    subroutine EvaluateProductionMu(iPhase,iFirst,iLast,dXLocal,dMuLocal,iInfoLocal)

        integer, intent(in) :: iPhase, iFirst, iLast
        real(8), intent(in) :: dXLocal(:)
        real(8), intent(out) :: dMuLocal(:)
        integer, intent(out) :: iInfoLocal

        iInfoLocal = 0
        if ((SIZE(dXLocal) /= iLast-iFirst+1) .OR. (SIZE(dMuLocal) /= SIZE(dXLocal)) .OR. &
            (ANY(dXLocal <= 0D0)) .OR. (DABS(SUM(dXLocal)-1D0) > 1D-10)) then
            iInfoLocal = 1
            return
        end if
        dMolFraction(iFirst:iLast) = dXLocal
        call CompExcessGibbsEnergySUBG(iPhase)
        dMuLocal = dChemicalPotential(iFirst:iLast)+dPartialExcessGibbs(iFirst:iLast)
        if ((INFOThermo /= 0) .OR. (.NOT. ALL(IEEE_IS_FINITE(dMuLocal)))) iInfoLocal = 2

    end subroutine EvaluateProductionMu


    subroutine EstimateStationarityUncertainty(dResidualMinus,dResidualPlus,dH,dZLocal, &
        dJacobianMinus,dJacobianPlus,dUncertainty,lAllPass)

        real(8), intent(in) :: dResidualMinus(:), dResidualPlus(:), dH
        real(8), intent(in) :: dZLocal(:,:), dJacobianMinus(:,:), dJacobianPlus(:,:)
        real(8), intent(out) :: dUncertainty
        logical, intent(inout) :: lAllPass

        integer :: iInfoLocal
        real(8), allocatable :: dCorrectionMinus(:), dCorrectionPlus(:)
        real(8), allocatable :: dRHSMinus(:,:), dRHSPlus(:,:)

        ! Translate each root's remaining stationarity defect into a composition
        ! correction with that root's production-FD Jacobian. The analytic SUBQ
        ! Hessian under test does not participate in this resolution decision.
        allocate(dRHSMinus(SIZE(dZLocal,2),1),dRHSPlus(SIZE(dZLocal,2),1), &
            dCorrectionMinus(SIZE(dZLocal,1)),dCorrectionPlus(SIZE(dZLocal,1)))
        dRHSMinus(:,1) = MATMUL(TRANSPOSE(dZLocal),dResidualMinus)
        dRHSPlus(:,1) = MATMUL(TRANSPOSE(dZLocal),dResidualPlus)
        call SolveDense(dJacobianMinus,dRHSMinus,iInfoLocal)
        lAllPass = lAllPass .AND. (iInfoLocal == 0)
        if (iInfoLocal == 0) call SolveDense(dJacobianPlus,dRHSPlus,iInfoLocal)
        lAllPass = lAllPass .AND. (iInfoLocal == 0)
        if (iInfoLocal == 0) then
            dCorrectionMinus = MATMUL(dZLocal,dRHSMinus(:,1))
            dCorrectionPlus = MATMUL(dZLocal,dRHSPlus(:,1))
            dUncertainty = (SQRT(SUM(dCorrectionMinus**2))+ &
                SQRT(SUM(dCorrectionPlus**2)))/(2D0*dH)
        else
            dUncertainty = HUGE(1D0)
        end if
        deallocate(dRHSMinus,dRHSPlus,dCorrectionMinus,dCorrectionPlus)

    end subroutine EstimateStationarityUncertainty


    subroutine ReconstructPhaseLocalGEMBaseline(iAssemblageSlot,dXLocal,dMuLocal,dS, &
        dA,dOuter,dAep,dResidualDifference)

        integer, intent(in) :: iAssemblageSlot
        real(8), intent(in) :: dXLocal(:), dMuLocal(:), dS(:,:)
        real(8), intent(out) :: dA(:,:), dOuter(:,:), dAep(:), dResidualDifference

        integer :: e, f, q
        real(8) :: dPhaseFromSpecies
        real(8), allocatable :: dResidualDirect(:), dResidualExpected(:)

        allocate(dResidualDirect(nElements),dResidualExpected(nElements))
        dA = 0D0
        dAep = 0D0
        dResidualDirect = 0D0
        do q = 1, SIZE(dXLocal)
            do e = 1, nElements
                dAep(e) = dAep(e)+dMolesSpecies(iFirstLocal+q-1)*dS(q,e)
                dResidualDirect(e) = dResidualDirect(e)+ &
                    dMolesSpecies(iFirstLocal+q-1)*dS(q,e)*(dMuLocal(q)-1D0)
                do f = 1, nElements
                    dA(e,f) = dA(e,f)+dMolesSpecies(iFirstLocal+q-1)*dS(q,e)*dS(q,f)
                end do
            end do
        end do
        dPhaseFromSpecies = SUM(dMolesSpecies(iFirstLocal:iLastLocal))
        dOuter = OuterProduct(dAep,dAep)/dPhaseFromSpecies
        dResidualExpected = dMolesPhase(iAssemblageSlot)*MATMUL(TRANSPOSE(dS), &
            dXLocal*(dMuLocal-1D0))
        dResidualDifference = MAXVAL(DABS(dResidualDirect-dResidualExpected))/ &
            DMAX1(1D0,MAXVAL(DABS(dResidualDirect)),MAXVAL(DABS(dResidualExpected)))

        deallocate(dResidualDirect,dResidualExpected)

    end subroutine ReconstructPhaseLocalGEMBaseline


    real(8) function ProjectedKKTResidual(dH,dC,dF,dR)

        real(8), intent(in) :: dH(:,:), dC(:,:), dF(:,:), dR(:,:)
        integer :: j
        real(8), allocatable :: dProjected(:), dResidualLocal(:)

        allocate(dProjected(SIZE(dH,1)),dResidualLocal(SIZE(dH,1)))
        ProjectedKKTResidual = 0D0
        if (SIZE(dC,1) /= 1) then
            ProjectedKKTResidual = HUGE(1D0)
            deallocate(dProjected,dResidualLocal)
            return
        end if
        do j = 1, SIZE(dF,2)
            dResidualLocal = MATMUL(dH,dR(:,j))-dF(:,j)
            dProjected = dResidualLocal-SUM(dResidualLocal)/DFLOAT(SIZE(dResidualLocal))
            ProjectedKKTResidual = DMAX1(ProjectedKKTResidual,MAXVAL(DABS(dProjected)))
        end do
        deallocate(dProjected,dResidualLocal)

    end function ProjectedKKTResidual


    subroutine FindBestResolvedPoint(dError,lResolved,iBest,dBest)

        real(8), intent(in) :: dError(:)
        logical, intent(in) :: lResolved(:)
        integer, intent(out) :: iBest
        real(8), intent(out) :: dBest
        integer :: iLocal

        iBest = 0
        dBest = HUGE(1D0)
        do iLocal = 1, SIZE(dError)
            if (.NOT. lResolved(iLocal)) cycle
            if (dError(iLocal) >= dBest) cycle
            iBest = iLocal
            dBest = dError(iLocal)
        end do

    end subroutine FindBestResolvedPoint


    logical function HasResolvedSecondOrderRegion(dError,dOrder,lOrder,lResolved)

        real(8), intent(in) :: dError(:), dOrder(:)
        logical, intent(in) :: lOrder(:), lResolved(:)
        integer :: iLocal

        HasResolvedSecondOrderRegion = .FALSE.
        do iLocal = 1, SIZE(dError)-2
            if (.NOT. ALL(lResolved(iLocal:iLocal+2))) cycle
            if (.NOT. lOrder(iLocal) .OR. .NOT. lOrder(iLocal+1)) cycle
            if (.NOT. ((dError(iLocal) > dError(iLocal+1)) .AND. &
                (dError(iLocal+1) > dError(iLocal+2)))) cycle
            if ((dOrder(iLocal) < FD_ORDER_SECOND_MIN) .OR. &
                (dOrder(iLocal) > FD_ORDER_SECOND_MAX)) cycle
            if ((dOrder(iLocal+1) < FD_ORDER_SECOND_MIN) .OR. &
                (dOrder(iLocal+1) > FD_ORDER_SECOND_MAX)) cycle
            HasResolvedSecondOrderRegion = .TRUE.
            return
        end do

    end function HasResolvedSecondOrderRegion


    subroutine SolveNullSpace(dReduced,dZLocal,dF,dR,iInfoLocal)

        real(8), intent(in) :: dReduced(:,:), dZLocal(:,:), dF(:,:)
        real(8), intent(out) :: dR(:,:)
        integer, intent(out) :: iInfoLocal
        integer :: iRefine
        real(8), allocatable :: dCorrectionLocal(:,:), dRHSLocal(:,:), dRHSOriginal(:,:)

        allocate(dCorrectionLocal(SIZE(dZLocal,2),SIZE(dF,2)), &
            dRHSLocal(SIZE(dZLocal,2),SIZE(dF,2)), &
            dRHSOriginal(SIZE(dZLocal,2),SIZE(dF,2)))
        dRHSOriginal = MATMUL(TRANSPOSE(dZLocal),dF)
        dRHSLocal = dRHSOriginal
        call SolveDense(dReduced,dRHSLocal,iInfoLocal)
        ! Iterative refinement reduces forward-solution disagreement caused by
        ! the strongly ill-conditioned tangent matrix while preserving the
        ! independently formulated null-space solve.
        do iRefine = 1, 2
            if (iInfoLocal /= 0) exit
            dCorrectionLocal = dRHSOriginal-MATMUL(dReduced,dRHSLocal)
            call SolveDense(dReduced,dCorrectionLocal,iInfoLocal)
            dRHSLocal = dRHSLocal+dCorrectionLocal
        end do
        dR = MATMUL(dZLocal,dRHSLocal)
        deallocate(dCorrectionLocal,dRHSLocal,dRHSOriginal)

    end subroutine SolveNullSpace


    subroutine SolveDense(dMatrix,dRHSLocal,iInfoLocal)

        real(8), intent(in) :: dMatrix(:,:)
        real(8), intent(inout) :: dRHSLocal(:,:)
        integer, intent(out) :: iInfoLocal
        integer :: n
        integer, allocatable :: iPivot(:)
        real(8), allocatable :: dWork(:,:)

        n = SIZE(dMatrix,1)
        allocate(dWork(n,n),iPivot(n))
        dWork = dMatrix
        call DGESV(n,SIZE(dRHSLocal,2),dWork,n,iPivot,dRHSLocal,n,iInfoLocal)
        deallocate(dWork,iPivot)

    end subroutine SolveDense


    subroutine SymmetricEigenvalues(dMatrix,dValue,iInfoLocal)

        real(8), intent(in) :: dMatrix(:,:)
        real(8), intent(out) :: dValue(:)
        integer, intent(out) :: iInfoLocal
        integer :: lWork, n
        real(8) :: dQuery(1)
        real(8), allocatable :: dWork(:), dMatrixWork(:,:)

        n = SIZE(dMatrix,1)
        allocate(dMatrixWork(n,n))
        dMatrixWork = dMatrix
        lWork = -1
        call DSYEV('N','U',n,dMatrixWork,n,dValue,dQuery,lWork,iInfoLocal)
        lWork = MAX(1,INT(dQuery(1)))
        allocate(dWork(lWork))
        dMatrixWork = dMatrix
        call DSYEV('N','U',n,dMatrixWork,n,dValue,dWork,lWork,iInfoLocal)
        deallocate(dWork,dMatrixWork)

    end subroutine SymmetricEigenvalues


    subroutine SingularValues(dMatrix,dValue,iInfoLocal)

        real(8), intent(in) :: dMatrix(:,:)
        real(8), intent(out) :: dValue(:)
        integer, intent(out) :: iInfoLocal

        ! The bordered KKT matrix is symmetric, so its singular values are the
        ! absolute values of its real eigenvalues. This avoids DGESVD's IEEE
        ! capability probe, which intentionally divides by zero and conflicts
        ! with Thermochimica's floating-point trap build.
        call SymmetricEigenvalues(dMatrix,dValue,iInfoLocal)
        dValue = DABS(dValue)

    end subroutine SingularValues


    subroutine SelectIndependentColumnsPivotedQR(dMatrix,dTolerance,iColumn,nRank)

        real(8), intent(in) :: dMatrix(:,:), dTolerance
        integer, intent(out) :: iColumn(:), nRank
        integer :: iColumnLocal, iPivot, j
        real(8) :: dNorm, dScale
        real(8), allocatable :: dCandidate(:), dQ(:,:), dResidualColumn(:,:), dResidualNorm(:)
        integer, allocatable :: iPermutation(:)

        ! Column pivoting chooses the largest remaining projected forcing at
        ! each QR step, so rank and basis selection do not depend on the input
        ! element ordering when columns have very different scales.
        allocate(dCandidate(SIZE(dMatrix,1)),dQ(SIZE(dMatrix,1),SIZE(dMatrix,2)), &
            dResidualColumn(SIZE(dMatrix,1),SIZE(dMatrix,2)),dResidualNorm(SIZE(dMatrix,2)), &
            iPermutation(SIZE(dMatrix,2)))
        dQ = 0D0
        dResidualColumn = dMatrix
        do iColumnLocal = 1, SIZE(dMatrix,2)
            iPermutation(iColumnLocal) = iColumnLocal
        end do
        iColumn = 0
        nRank = 0
        dScale = DMAX1(1D0,MAXVAL(SQRT(SUM(dMatrix*dMatrix,DIM=1))))
        do iColumnLocal = 1, SIZE(dMatrix,2)
            dResidualNorm(iColumnLocal:) = SQRT(SUM( &
                dResidualColumn(:,iColumnLocal:)*dResidualColumn(:,iColumnLocal:),DIM=1))
            iPivot = iColumnLocal-1+MAXLOC(dResidualNorm(iColumnLocal:),DIM=1)
            if (iPivot /= iColumnLocal) then
                dCandidate = dResidualColumn(:,iColumnLocal)
                dResidualColumn(:,iColumnLocal) = dResidualColumn(:,iPivot)
                dResidualColumn(:,iPivot) = dCandidate
                j = iPermutation(iColumnLocal)
                iPermutation(iColumnLocal) = iPermutation(iPivot)
                iPermutation(iPivot) = j
            end if
            dCandidate = dResidualColumn(:,iColumnLocal)
            dNorm = SQRT(SUM(dCandidate*dCandidate))
            if (dNorm <= dTolerance*dScale) exit
            nRank = nRank+1
            iColumn(nRank) = iPermutation(iColumnLocal)
            dQ(:,nRank) = dCandidate/dNorm
            do j = iColumnLocal+1, SIZE(dMatrix,2)
                dResidualColumn(:,j) = dResidualColumn(:,j)- &
                    DOT_PRODUCT(dQ(:,nRank),dResidualColumn(:,j))*dQ(:,nRank)
            end do
        end do
        deallocate(dCandidate,dQ,dResidualColumn,dResidualNorm,iPermutation)

    end subroutine SelectIndependentColumnsPivotedQR


    function DiagonalMatrix(dVector) result(dMatrix)

        real(8), intent(in) :: dVector(:)
        real(8) :: dMatrix(SIZE(dVector),SIZE(dVector))
        integer :: iLocal

        dMatrix = 0D0
        do iLocal = 1, SIZE(dVector)
            dMatrix(iLocal,iLocal) = dVector(iLocal)
        end do

    end function DiagonalMatrix


    function OuterProduct(dLeft,dRight) result(dMatrix)

        real(8), intent(in) :: dLeft(:), dRight(:)
        real(8) :: dMatrix(SIZE(dLeft),SIZE(dRight))
        integer :: iLocal

        do iLocal = 1, SIZE(dLeft)
            dMatrix(iLocal,:) = dLeft(iLocal)*dRight
        end do

    end function OuterProduct


    character(len=16) function OrderLabel(dOrder,lAvailable)

        real(8), intent(in) :: dOrder
        logical, intent(in) :: lAvailable

        if (lAvailable) then
            write(OrderLabel,'(F10.4)') dOrder
        else
            OrderLabel = 'N/A'
        end if

    end function OrderLabel


    subroutine FinishTest(lAllPass)

        logical, intent(in) :: lAllPass

        if (lAllPass) then
            print *, 'TestMQMQASUBQResponseVerification: PASS'
            call EXIT(0)
        else
            print *, 'TestMQMQASUBQResponseVerification: FAIL <---'
            call EXIT(1)
        end if

    end subroutine FinishTest

end program TestMQMQASUBQResponseVerification
