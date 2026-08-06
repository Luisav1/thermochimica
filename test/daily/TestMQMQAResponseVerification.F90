!-------------------------------------------------------------------------------------------------------------
!> \file    TestMQMQAResponseVerification.F90
!> \brief   Diagnostic-only native verification of the plain-SUBG constrained composition response.
!>
!> \details This MQ-3B test starts from the converged TestThermo56 Liquid,
!!          decodes its verified local MQMQA Hessian, imposes the phase
!!          normalization constraint, and predicts how its quadruplet
!!          composition changes when independent element potentials are
!!          perturbed. A test-only nonlinear solver of the established
!!          production partial-molar stationarity equations supplies the
!!          independent oracle without using the new analytic Hessian.
!!
!!          The test also verifies the existing ideal/simple GEM baseline in
!!          three forms: a bordered local solve, its closed normalized formula,
!!          and a phase-local reconstruction from the arrays consumed by
!!          GEMNewton. It does not assemble a mapper correction or modify GEM.
!-------------------------------------------------------------------------------------------------------------

program TestMQMQAResponseVerification

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
    integer :: i, iDirection, iFirstLocal, iInfo, iLastLocal, iPhaseIndex, iRankForcing
    integer :: iSlot, iStep, nQuad, nTangent
    integer, allocatable :: iBasisColumn(:), iBestResolved(:), iWorstComponent(:,:)
    logical :: lPass, lReport, lZeroOracle
    real(8) :: dAepError, dBaselineClosedError, dBaselineGEMError
    real(8) :: dFloorDifference, dKKTCondition, dKKTConstraintResidual, dKKTTopResidual
    real(8) :: dMaxEigenvalue, dMinEigenvalue, dNullDifference, dProjectedKKTResidual
    real(8) :: dResidualError, dResponseBest, dResponseRawBest, dStabilityTolerance
    real(8) :: dStationarityBest, dWorstComponentBest, dZeroOracleDrift, dZeroStationarity
    real(8) :: dTangentCondition
    real(8), allocatable :: dAClosed(:,:), dAClosedPhase(:,:), dAFromGEM(:,:), dAOuter(:,:)
    real(8), allocatable :: dAepDirect(:), dAepExpected(:), dConstraint(:,:)
    real(8), allocatable :: dEigenvalue(:), dErrors(:,:), dForcing(:,:), dForcingBasis(:,:)
    real(8), allocatable :: dGradient(:), dHbase(:,:), dHessian(:,:), dHx(:,:)
    real(8), allocatable :: dKKT(:,:), dKKTSingular(:), dKtangent(:,:), dMoles(:)
    real(8), allocatable :: dMolesFloor(:), dMuProduction(:), dNullResponse(:,:), dOrders(:,:)
    real(8), allocatable :: dKKTMultiplier(:,:), dMaxAbsolute(:,:), dMaxScaled(:,:)
    real(8), allocatable :: dNormAbsolute(:,:), dP(:,:), dProjectedForcing(:,:)
    real(8), allocatable :: dRawMinimum(:), dResolvedBest(:), dResponse(:,:), dResponseBase(:,:)
    real(8), allocatable :: dResponseClosed(:,:), dSteps(:,:), dStationarity(:,:), dX(:), dZ(:,:)
    real(8), allocatable :: dZeroResidual(:), dZeroX(:)
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
    ! SECTION 1: CONVERGED, INTERIOR, UNCHARGED PRODUCTION STATE
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
            if ((cSolnPhaseName(-iAssemblage(i)) == 'Liquid') .AND. &
                (cSolnPhaseType(-iAssemblage(i)) == 'SUBG')) then
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
    allocate(dX(nQuad),dMoles(nQuad),dMolesFloor(nQuad),dGradient(nQuad),dMuProduction(nQuad), &
        dHessian(nQuad,nQuad),dHx(nQuad,nQuad),dHbase(nQuad,nQuad), &
        dConstraint(1,nQuad),dForcing(nQuad,nElements),dZ(nQuad,nTangent), &
        dKtangent(nTangent,nTangent),dEigenvalue(nTangent),iBasisColumn(nElements))

    dX = dMolFraction(iFirstLocal:iLastLocal)
    dX = dX/SUM(dX)
    dMoles = dMolesPhase(iSlot)*dX
    dMolesFloor = dMolesSpecies(iFirstLocal:iLastLocal)
    dFloorDifference = MAXVAL(DABS(dMolesFloor-dMoles))/ &
        DMAX1(1D0,MAXVAL(DABS(dMoles)))
    lPass = lPass .AND. (MINVAL(dX) > 1D-12)

    call DecodeProductionSUBGPhase(iPhaseIndex,tModel,tInteraction,iInfo)
    lPass = lPass .AND. (iInfo == 0)
    if (iInfo == 0) then
        call CompMQMQAHessianUnconstrained(tModel,dMoles,1D0,tInteraction,dHessian,iInfo, &
            dGradient=dGradient)
        lPass = lPass .AND. (iInfo == 0)
    end if
    dHx = dMolesPhase(iSlot)*dHessian
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
    call SelectIndependentColumns(dProjectedForcing,1D-11,iBasisColumn,iRankForcing)
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
    ! The reduced tangent matrix is strongly ill-conditioned; agreement at 1D-9 is
    ! consistent with two independently solved systems whose direct residuals
    ! remain near machine precision.
    lPass = lPass .AND. (iInfo == 0) .AND. (dNullDifference <= 1D-9)

    allocate(dKKT(nQuad+1,nQuad+1),dKKTSingular(nQuad+1))
    dKKT = 0D0
    dKKT(1:nQuad,1:nQuad) = dHx
    dKKT(1:nQuad,nQuad+1) = 1D0
    dKKT(nQuad+1,1:nQuad) = 1D0
    call SingularValues(dKKT,dKKTSingular,iInfo)
    lPass = lPass .AND. (iInfo == 0) .AND. (MINVAL(dKKTSingular) > 0D0)
    dKKTCondition = MAXVAL(dKKTSingular)/MINVAL(dKKTSingular)

    !=========================================================================================================
    ! SECTION 3: THREE-WAY, PHASE-LOCAL IDEAL GEM BASELINE
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

    allocate(dZeroX(nQuad),dZeroResidual(nQuad))
    call EvaluateProductionOracleState(iPhaseIndex,iFirstLocal,iLastLocal,dX, &
        dElementPotential*0D0,0D0,dZeroX,dZeroResidual,lZeroOracle)
    lPass = lPass .AND. lZeroOracle
    dZeroOracleDrift = SQRT(SUM((dZeroX-dX)**2))
    dZeroStationarity = SQRT(SUM(dZeroResidual*dZeroResidual))
    lPass = lPass .AND. (dZeroOracleDrift <= 1D-5) .AND. (dZeroStationarity <= 1D-11)

    ! The stored production composition satisfies Thermochimica's historical
    ! update-size tolerance but is not an exact root of the stationarity
    ! equations. Keep that stored state for the GEM baseline comparison above;
    ! center the differential response check on the nearby, explicitly reported
    ! production-equation root so analytic and numerical derivatives refer to
    ! the same local state.
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
    lPass = lPass .AND. (iInfo == 0) .AND. (dNullDifference <= 1D-9)
    dKKT = 0D0
    dKKT(1:nQuad,1:nQuad) = dHx
    dKKT(1:nQuad,nQuad+1) = 1D0
    dKKT(nQuad+1,1:nQuad) = 1D0
    call SingularValues(dKKT,dKKTSingular,iInfo)
    lPass = lPass .AND. (iInfo == 0) .AND. (MINVAL(dKKTSingular) > 0D0)
    dKKTCondition = MAXVAL(dKKTSingular)/MINVAL(dKKTSingular)

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
            dResponse(:,iDirection),dZ,dKtangent,dSteps(iDirection,:),dErrors(iDirection,:), &
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
        write(*,'(A)') 'MQ-3B native plain-SUBG constrained-response verification'
        write(*,'(A)') 'scope: local response and ideal GEM baseline only; no deltaA/deltaB or GEM activation'
        write(*,'(A,A)') 'phase = ',TRIM(cSolnPhaseName(iPhaseIndex))
        write(*,'(A,I0)') 'quadruplet count = ',nQuad
        write(*,'(A,I0)') 'constraint rank = ',1
        write(*,'(A,I0)') 'supported forcing rank = ',iRankForcing
        write(*,'(A,I0)') 'null/dependent original forcing columns = ',nElements-iRankForcing
        write(*,'(A,*(I0,1X))') 'forcing basis element columns = ',iBasisColumn(1:iRankForcing)
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
        dZLocal,dKtangentLocal,dStepValues,dErrorValues,dNormAbsoluteValues, &
        dMaxAbsoluteValues,dMaxScaledValues,iWorstValues,dStationarityValues,lAllPass)

        integer, intent(in) :: iPhase, iFirst, iLast
        real(8), intent(in) :: dXBase(:), dGammaDirection(:), dPrediction(:)
        real(8), intent(in) :: dZLocal(:,:), dKtangentLocal(:,:)
        real(8), intent(out) :: dStepValues(:), dErrorValues(:), dNormAbsoluteValues(:)
        real(8), intent(out) :: dMaxAbsoluteValues(:), dMaxScaledValues(:), dStationarityValues(:)
        integer, intent(out) :: iWorstValues(:)
        logical, intent(inout) :: lAllPass

        integer :: iStepLocal
        logical :: lMinus, lPlus
        real(8) :: dH, dHMaximum
        real(8), allocatable :: dFD(:), dResidualMinus(:), dResidualPlus(:), dXMinus(:), dXPlus(:)

        allocate(dFD(SIZE(dXBase)),dResidualMinus(SIZE(dXBase)),dResidualPlus(SIZE(dXBase)), &
            dXMinus(SIZE(dXBase)),dXPlus(SIZE(dXBase)))
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
                dXMinus,dResidualMinus,lMinus)
            call EvaluateProductionOracleState(iPhase,iFirst,iLast,dXBase,dGammaDirection,dH, &
                dXPlus,dResidualPlus,lPlus)
            lAllPass = lAllPass .AND. lMinus .AND. lPlus
            dFD = (dXPlus-dXMinus)/(2D0*dH)
            call ComputeVectorErrorMetrics(dFD,dPrediction,dNormAbsoluteValues(iStepLocal), &
                dErrorValues(iStepLocal),dMaxAbsoluteValues(iStepLocal), &
                dMaxScaledValues(iStepLocal),iWorstValues(iStepLocal))
            call EstimateStationarityUncertainty(dResidualMinus,dResidualPlus,dH, &
                dZLocal,dKtangentLocal,dStationarityValues(iStepLocal),lAllPass)
        end do

        deallocate(dFD,dResidualMinus,dResidualPlus,dXMinus,dXPlus)

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
        dXResult,dResidual,lSuccess)

        integer, intent(in) :: iPhase, iFirst, iLast
        real(8), intent(in) :: dXBase(:), dGammaDirection(:), dH
        real(8), intent(out) :: dXResult(:), dResidual(:)
        logical, intent(out) :: lSuccess

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
            dXResult,dResidual,lOracleConverged)
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
        dXResult,dResidualResult,lConverged)

        integer, intent(in) :: iPhase, iFirst, iLast
        real(8), intent(in) :: dXStart(:), dGammaLocal(:)
        real(8), intent(out) :: dXResult(:), dResidualResult(:)
        logical, intent(out) :: lConverged

        integer :: iDirectionLocal, iInfoLocal, iIteration, iTrial
        real(8) :: dAlpha, dHJacobian, dResidualNorm, dTrialNorm
        real(8), allocatable :: dDelta(:), dJacobian(:,:), dMuMinus(:), dMuPlus(:)
        real(8), allocatable :: dMuTrial(:), dReducedResidual(:,:), dResidualTrial(:)
        real(8), allocatable :: dXMinus(:), dXPlus(:), dXTrial(:)

        allocate(dDelta(SIZE(dXStart)),dJacobian(SIZE(dZ,2),SIZE(dZ,2)), &
            dMuMinus(SIZE(dXStart)),dMuPlus(SIZE(dXStart)),dMuTrial(SIZE(dXStart)), &
            dReducedResidual(SIZE(dZ,2),1),dResidualTrial(SIZE(dXStart)), &
            dXMinus(SIZE(dXStart)),dXPlus(SIZE(dXStart)),dXTrial(SIZE(dXStart)))
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

            do iDirectionLocal = 1, SIZE(dZ,2)
                where (DABS(dZ(:,iDirectionLocal)) > 0D0)
                    dDelta = dXResult/DABS(dZ(:,iDirectionLocal))
                elsewhere
                    dDelta = HUGE(1D0)
                end where
                dHJacobian = DMIN1(1D-5,0.1D0*MINVAL(dDelta))
                dXMinus = dXResult-dHJacobian*dZ(:,iDirectionLocal)
                dXPlus = dXResult+dHJacobian*dZ(:,iDirectionLocal)
                call EvaluateProductionMu(iPhase,iFirst,iLast,dXMinus,dMuMinus,iInfoLocal)
                if (iInfoLocal /= 0) exit
                call EvaluateProductionMu(iPhase,iFirst,iLast,dXPlus,dMuPlus,iInfoLocal)
                if (iInfoLocal /= 0) exit
                dJacobian(:,iDirectionLocal) = MATMUL(TRANSPOSE(dZ), &
                    (dMuPlus-dMuMinus))/(2D0*dHJacobian)
            end do
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
        deallocate(dDelta,dJacobian,dMuMinus,dMuPlus,dMuTrial,dReducedResidual, &
            dResidualTrial,dXMinus,dXPlus,dXTrial)

    end subroutine SolveProductionStationarity


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
        dKtangentLocal,dUncertainty,lAllPass)

        real(8), intent(in) :: dResidualMinus(:), dResidualPlus(:), dH
        real(8), intent(in) :: dZLocal(:,:), dKtangentLocal(:,:)
        real(8), intent(out) :: dUncertainty
        logical, intent(inout) :: lAllPass

        integer :: iInfoLocal
        real(8), allocatable :: dCorrection(:,:), dRHSLocal(:,:)

        allocate(dRHSLocal(SIZE(dZLocal,2),2),dCorrection(SIZE(dZLocal,1),2))
        dRHSLocal(:,1) = MATMUL(TRANSPOSE(dZLocal),dResidualMinus)
        dRHSLocal(:,2) = MATMUL(TRANSPOSE(dZLocal),dResidualPlus)
        call SolveDense(dKtangentLocal,dRHSLocal,iInfoLocal)
        lAllPass = lAllPass .AND. (iInfoLocal == 0)
        dCorrection = MATMUL(dZLocal,dRHSLocal)
        dUncertainty = (SQRT(SUM(dCorrection(:,1)**2))+ &
            SQRT(SUM(dCorrection(:,2)**2)))/(2D0*dH)
        deallocate(dRHSLocal,dCorrection)

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


    subroutine SelectIndependentColumns(dMatrix,dTolerance,iColumn,nRank)

        real(8), intent(in) :: dMatrix(:,:), dTolerance
        integer, intent(out) :: iColumn(:), nRank
        integer :: iColumnLocal, j
        real(8) :: dNorm, dScale
        real(8), allocatable :: dCandidate(:), dQ(:,:)

        allocate(dCandidate(SIZE(dMatrix,1)),dQ(SIZE(dMatrix,1),SIZE(dMatrix,2)))
        dQ = 0D0
        iColumn = 0
        nRank = 0
        dScale = DMAX1(1D0,MAXVAL(SQRT(SUM(dMatrix*dMatrix,DIM=1))))
        do iColumnLocal = 1, SIZE(dMatrix,2)
            dCandidate = dMatrix(:,iColumnLocal)
            do j = 1, nRank
                dCandidate = dCandidate-DOT_PRODUCT(dQ(:,j),dCandidate)*dQ(:,j)
            end do
            dNorm = SQRT(SUM(dCandidate*dCandidate))
            if (dNorm <= dTolerance*dScale) cycle
            nRank = nRank+1
            iColumn(nRank) = iColumnLocal
            dQ(:,nRank) = dCandidate/dNorm
        end do
        deallocate(dCandidate,dQ)

    end subroutine SelectIndependentColumns


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
            print *, 'TestMQMQAResponseVerification: PASS'
            call EXIT(0)
        else
            print *, 'TestMQMQAResponseVerification: FAIL <---'
            call EXIT(1)
        end if

    end subroutine FinishTest

end program TestMQMQAResponseVerification
