!-------------------------------------------------------------------------------------------------------------
!> \file    TestMQMQASUBQNativeEnergyVerification.F90
!> \brief   Thermochimica-native scalar-energy verification of the standalone SUBQ model.
!>
!> \details Reproduce the Fe-Ti-V-O calculation used by TestThermo57, decode the
!!          active SlagBsoln phase through ModuleMQMQAProductionAdapter, and
!!          compare the standalone SUBQ scalar energy with established production
!!          Thermochimica at the same converged quadruplet composition.
!!
!!          MQ-2A verifies the production-data translation and scalar equations:
!!          - the decoder selects the SUBQ formulation explicitly;
!!          - real FeTiVO topology, coordination, zeta, reference energies, and
!!            active G/Q parameters reach the standalone model;
!!          - reference/configurational and excess energy blocks agree separately.
!!
!!          FeTiVO uses uniform zeta=2.4. This is therefore a native assessed
!!          SUBQ G/Q case, not verification of every SUBQ feature and not a
!!          native test of the nonuniform pair-specific-zeta distinction. A
!!          controlled modified-runtime section supplies positive nonuniform
!!          zeta values, verifies the corrected weighted S3 scalar block and
!!          complete production partial molars, and restores the parsed row.
!!
!!          The production comparison uses Euler sums of production partial
!!          molars. CompExcessGibbsEnergySUBG returns reference/configurational
!!          terms in dChemicalPotential and excess terms in
!!          dPartialExcessGibbs, so the two blocks remain independently visible.
!!          Hessian-vector finite differences are intentionally deferred to MQ-2B.
!!          Pass --report to print the decoded native-state evidence.
!-------------------------------------------------------------------------------------------------------------

program TestMQMQASUBQNativeEnergyVerification

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleMQMQAUnconstrained
    USE ModuleMQMQAProductionAdapter
    USE ModuleFiniteDifferenceVerification, ONLY: ComputeVectorErrorMetrics

    implicit none

    interface
        subroutine CompExcessGibbsEnergySUBG(iSolnIndex)
            integer, intent(in) :: iSolnIndex
        end subroutine CompExcessGibbsEnergySUBG
    end interface

    integer :: i, iFirst, iInfo, iLast, iPhaseIndex, iSlot, nG, nQ, nQuad
    logical :: lPass, lReport
    character(len=32) :: cArgument
    real(8) :: dG, dGExcess, dGIdeal, dGReference, dMinMoles
    real(8) :: dProductionExcess, dProductionReferenceIdeal
    real(8) :: dErrorExcess, dErrorReferenceIdeal, dErrorTotal
    real(8) :: dUniformLegacyS3, dUniformPairFractionDifference
    real(8) :: dUniformS3Difference, dUniformWeightedS3
    real(8), allocatable :: dMoles(:)
    type(MQMQAModelData) :: tModel, tRejectedModel
    type(MQMQAInteractionTerm), allocatable :: tInteraction(:), tRejectedInteraction(:)

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
    ! SECTION 1: CONVERGED PRODUCTION SUBQ STATE
    !
    ! TestThermo57 supplies all Fe, Ti, and V constituents needed for an
    ! interior SlagBsoln composition and already provides regression evidence
    ! that this phase is stable at the selected temperature and pressure.
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
            if (cSolnPhaseName(-iAssemblage(i)) == 'SlagBsoln' .AND. &
                cSolnPhaseType(-iAssemblage(i)) == 'SUBQ') then
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
        allocate(dMoles(nQuad))
        dMoles = dMolesSpecies(iFirst:iLast)
        dMinMoles = MINVAL(dMoles)
        lPass = lPass .AND. ALL(dMoles > 0D0)

        !=====================================================================================================
        ! SECTION 2: STRICT PRODUCTION DECODING
        !
        ! The old SUBG entry point must continue to reject this phase. The new
        ! SUBQ entry point then selects the different standalone formulation
        ! while reusing the source-faithful production topology translation.
        !=====================================================================================================
        call DecodeProductionSUBGPhase(iPhaseIndex,tRejectedModel,tRejectedInteraction,iInfo)
        lPass = lPass .AND. (iInfo == 2)

        call DecodeProductionSUBQPhase(iPhaseIndex,tModel,tInteraction,iInfo)
        lPass = lPass .AND. (iInfo == 0)
        if (iInfo == 0) then
            lPass = lPass .AND. (tModel%iModelType == MQMQA_MODEL_SUBQ)
            lPass = lPass .AND. (nQuad == 15)
            lPass = lPass .AND. (SIZE(tModel%iQuadruplet,1) == nQuad)
            lPass = lPass .AND. (MINVAL(tModel%dZeta) > 0D0)
            lPass = lPass .AND. (MAXVAL(DABS(tModel%dZeta-2.4D0)) <= 1D-12)

            nG = COUNT(tInteraction%iFamily == MQMQA_TERM_G)
            nQ = COUNT(tInteraction%iFamily == MQMQA_TERM_Q)
            lPass = lPass .AND. (nG == 6) .AND. (nQ == 8)

            call EvaluateS3Definitions(tModel,dMoles,dUniformLegacyS3,dUniformWeightedS3,iInfo, &
                dUniformPairFractionDifference)
            lPass = lPass .AND. (iInfo == 0)
            dUniformS3Difference = ScaledError(dUniformLegacyS3,dUniformWeightedS3)
            lPass = lPass .AND. (dUniformPairFractionDifference <= 1D-14)
            lPass = lPass .AND. (dUniformS3Difference <= 1D-14)

            !=================================================================================================
            ! SECTION 3: BLOCKWISE SCALAR-ENERGY PARITY
            !
            ! Each extensive standalone block is divided by total quadruplet
            ! moles before comparison because the production Euler sum is a
            ! molar phase energy at the imposed composition.
            !=================================================================================================
            call CompMQMQAGibbsEnergyUnconstrained(tModel,dMoles,1D0,tInteraction, &
                dG,dGReference,dGIdeal,dGExcess,iInfo)
            lPass = lPass .AND. (iInfo == 0)
            if (iInfo == 0) then
                call EvaluateProductionEnergy(iPhaseIndex,dMoles,dProductionReferenceIdeal, &
                    dProductionExcess,iInfo)
                lPass = lPass .AND. (iInfo == 0)
            end if

            if (iInfo == 0) then
                dErrorReferenceIdeal = ScaledError((dGReference+dGIdeal)/SUM(dMoles), &
                    dProductionReferenceIdeal)
                dErrorExcess = ScaledError(dGExcess/SUM(dMoles),dProductionExcess)
                dErrorTotal = ScaledError(dG/SUM(dMoles), &
                    dProductionReferenceIdeal+dProductionExcess)
                lPass = lPass .AND. (dErrorReferenceIdeal <= 1D-10)
                lPass = lPass .AND. (dErrorExcess <= 1D-10)
                lPass = lPass .AND. (dErrorTotal <= 1D-10)

                if (lReport) then
                    write(*,'(A)') 'MQ-2A Thermochimica-native SUBQ scalar-energy verification'
                    write(*,'(A)') 'native scope: nonmagnetic SUBQ reference/configurational/G/Q families'
                    write(*,'(A)') 'native exclusions: B, R, Hessian FD, constrained response, GEM integration'
                    write(*,'(A,A)') 'database = ',TRIM(cThermoFileName)
                    write(*,'(A,A)') 'phase = ',TRIM(cSolnPhaseName(iPhaseIndex))
                    write(*,'(A,I0)') 'quadruplet count = ',nQuad
                    write(*,'(A,I0)') 'active decoded G-family terms = ',nG
                    write(*,'(A,I0)') 'active decoded Q-family terms = ',nQ
                    write(*,'(A,ES14.6)') 'phase amount = ',dMolesPhase(iSlot)
                    write(*,'(A,ES14.6)') 'minimum quadruplet moles = ',dMinMoles
                    write(*,'(A,2ES18.8)') 'zeta range = ',MINVAL(tModel%dZeta),MAXVAL(tModel%dZeta)
                    write(*,'(A,ES14.6)') 'maximum weighted/ordinary pair-fraction difference = ', &
                        dUniformPairFractionDifference
                    write(*,'(A,ES14.6)') 'uniform-zeta weighted/ordinary S3 difference = ', &
                        dUniformS3Difference
                    write(*,'(A)') 'block                         standalone molar      production molar      scaled error'
                    write(*,'(A28,3ES22.11)') 'reference+configurational', &
                        (dGReference+dGIdeal)/SUM(dMoles),dProductionReferenceIdeal,dErrorReferenceIdeal
                    write(*,'(A28,3ES22.11)') 'G/Q excess',dGExcess/SUM(dMoles), &
                        dProductionExcess,dErrorExcess
                    write(*,'(A28,3ES22.11)') 'total',dG/SUM(dMoles), &
                        dProductionReferenceIdeal+dProductionExcess,dErrorTotal
                end if
            end if
        end if

        ! The assessed FeTiVO zeta values are uniform, so ordinary and weighted
        ! pair fractions coincide. A temporary nonuniform runtime row is needed
        ! to make this test sensitive to the historical SUBQ S3 array selection.
        call VerifyControlledNonuniformZeta(iPhaseIndex,dMoles,lPass,lReport)
    end if

    call ResetThermoAll
    if (lPass) then
        print *, 'TestMQMQASUBQNativeEnergyVerification: PASS'
        call EXIT(0)
    else
        print *, 'TestMQMQASUBQNativeEnergyVerification: FAIL <---'
        call EXIT(1)
    end if

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Exercise the production SUBQ S3 selection with nonuniform runtime zeta data.
    !>
    !> \details FeTiVO's assessed zeta values are uniform and therefore cannot
    !!          distinguish ordinary from zeta-weighted normalized pair fractions.
    !!          This helper saves the complete production row, installs a positive
    !!          deterministic nonuniform row for the active pair records, evaluates
    !!          production and independent standalone oracles, and restores the row
    !!          through one cleanup path. No database data or equilibrium state is
    !!          permanently changed.
    !---------------------------------------------------------------------------------------------------------
    subroutine VerifyControlledNonuniformZeta(iPhaseLocal,dMolesLocal,lAllPass,lVerbose)

        integer, intent(in) :: iPhaseLocal
        real(8), intent(in) :: dMolesLocal(:)
        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lVerbose

        integer, parameter :: nExpectedPairs = 5
        integer :: a, iGradientWorst, iInfoLocal, iPair, iSPI, x
        logical :: lCasePass, lDiagnosticsValid, lRestored
        real(8) :: dCorrectedMolar, dDeltaS3, dExcessLocal, dG, dGExcessLocal
        real(8) :: dGIdealLocal, dGLegacyMolar, dGReferenceLocal, dGradientAbsolute
        real(8) :: dGradientError, dGradientMaxAbsolute, dGradientMaxScaled
        real(8) :: dLegacyS3, dLegacySeparation, dParityError, dProductionReferenceIdeal
        real(8) :: dExcessParityError, dTotalParityError, dWeightedS3
        real(8), allocatable :: dGradient(:), dHessian(:,:), dMuProduction(:), dState(:), dZetaSave(:)
        type(MQMQAModelData) :: tModifiedModel
        type(MQMQAInteractionTerm), allocatable :: tModifiedInteraction(:)

        lCasePass = .TRUE.
        lDiagnosticsValid = .FALSE.
        iSPI = iPhaseSublattice(iPhaseLocal)
        if (iSPI < 1) then
            lAllPass = .FALSE.
            if (lVerbose) write(*,'(/,A,I0)') 'controlled zeta fixture invalid phase-storage index = ',iSPI
            return
        end if
        if (nPairsSRO(iSPI,1) /= nExpectedPairs) then
            lAllPass = .FALSE.
            if (lVerbose) write(*,'(/,A,I0)') 'controlled zeta fixture unexpected pair count = ', &
                nPairsSRO(iSPI,1)
            return
        end if

        allocate(dZetaSave(SIZE(dZetaSpecies,2)),dGradient(SIZE(dMolesLocal)), &
            dHessian(SIZE(dMolesLocal),SIZE(dMolesLocal)),dMuProduction(SIZE(dMolesLocal)), &
            dState(SIZE(dMolesLocal)))
        dZetaSave = dZetaSpecies(iSPI,:)
        dState = 0.8D0*dMolesLocal+0.2D0*SUM(dMolesLocal)/DFLOAT(SIZE(dMolesLocal))

        ! No early return is permitted after this mutation. Every outcome flows
        ! through the exact restoration and restoration assertion below.
        do iPair = 1, nExpectedPairs
            dZetaSpecies(iSPI,iPair) = 1.8D0+0.15D0*DFLOAT(iPair)
        end do

        dCorrectedMolar = 0D0
        dDeltaS3 = 0D0
        dExcessLocal = 0D0
        dG = 0D0
        dGExcessLocal = 0D0
        dGIdealLocal = 0D0
        dGLegacyMolar = 0D0
        dGReferenceLocal = 0D0
        dGradientAbsolute = HUGE(1D0)
        dGradientError = HUGE(1D0)
        dGradientMaxAbsolute = HUGE(1D0)
        dGradientMaxScaled = HUGE(1D0)
        dLegacyS3 = 0D0
        dLegacySeparation = 0D0
        dParityError = HUGE(1D0)
        dProductionReferenceIdeal = 0D0
        dExcessParityError = HUGE(1D0)
        dTotalParityError = HUGE(1D0)
        dWeightedS3 = 0D0
        dGradient = 0D0
        dHessian = 0D0
        dMuProduction = 0D0
        iGradientWorst = 0

        call DecodeProductionSUBQPhase(iPhaseLocal,tModifiedModel,tModifiedInteraction,iInfoLocal)
        lCasePass = lCasePass .AND. (iInfoLocal == 0)
        if (iInfoLocal == 0) then
            do iPair = 1, nExpectedPairs
                a = iConstituentSublattice(iSPI,1,iPair)
                x = iConstituentSublattice(iSPI,2,iPair)
                lCasePass = lCasePass .AND. &
                    (tModifiedModel%dZeta(a,x) == 1.8D0+0.15D0*DFLOAT(iPair))
            end do

            call CompMQMQAGibbsEnergyUnconstrained(tModifiedModel,dState,1D0, &
                tModifiedInteraction,dG,dGReferenceLocal,dGIdealLocal,dGExcessLocal,iInfoLocal)
            lCasePass = lCasePass .AND. (iInfoLocal == 0)
        end if
        if (iInfoLocal == 0) then
            call CompMQMQAHessianUnconstrained(tModifiedModel,dState,1D0, &
                tModifiedInteraction,dHessian,iInfoLocal,dGradient=dGradient)
            lCasePass = lCasePass .AND. (iInfoLocal == 0)
        end if
        if (iInfoLocal == 0) then
            call EvaluateProductionEnergy(iPhaseLocal,dState,dProductionReferenceIdeal, &
                dExcessLocal,iInfoLocal,dMuProduction)
            lCasePass = lCasePass .AND. (iInfoLocal == 0)
        end if
        if (iInfoLocal == 0) then
            call EvaluateS3Definitions(tModifiedModel,dState,dLegacyS3,dWeightedS3,iInfoLocal)
            lCasePass = lCasePass .AND. (iInfoLocal == 0)
        end if

        if (iInfoLocal == 0) then
            dCorrectedMolar = (dGReferenceLocal+dGIdealLocal)/SUM(dState)
            dDeltaS3 = dWeightedS3-dLegacyS3
            dGLegacyMolar = dCorrectedMolar-dDeltaS3/SUM(dState)
            dParityError = ScaledError(dProductionReferenceIdeal,dCorrectedMolar)
            dExcessParityError = ScaledError(dExcessLocal,dGExcessLocal/SUM(dState))
            dTotalParityError = ScaledError(dProductionReferenceIdeal+dExcessLocal,dG/SUM(dState))
            dLegacySeparation = ScaledError(dProductionReferenceIdeal,dGLegacyMolar)
            call ComputeVectorErrorMetrics(dGradient,dMuProduction,dGradientAbsolute,dGradientError, &
                dGradientMaxAbsolute,dGradientMaxScaled,iGradientWorst)

            lCasePass = lCasePass .AND. (dParityError <= 1D-10)
            lCasePass = lCasePass .AND. (dExcessParityError <= 1D-10)
            lCasePass = lCasePass .AND. (dTotalParityError <= 1D-10)
            lCasePass = lCasePass .AND. (dLegacySeparation >= 1D-8)
            lCasePass = lCasePass .AND. &
                (dLegacySeparation >= 100D0*DMAX1(dParityError,EPSILON(1D0)))
            lCasePass = lCasePass .AND. (dGradientError <= 1D-10)
            lCasePass = lCasePass .AND. (dGradientMaxScaled <= 1D-10)
            lDiagnosticsValid = .TRUE.
        end if

        dZetaSpecies(iSPI,:) = dZetaSave
        lRestored = ALL(dZetaSpecies(iSPI,:) == dZetaSave)
        lCasePass = lCasePass .AND. lRestored
        lAllPass = lAllPass .AND. lCasePass

        if (lVerbose) then
            write(*,'(/,A)') 'Controlled nonuniform-zeta SUBQ S3 production regression'
            write(*,'(A)') 'scope: modified runtime zeta row; assessed FeTiVO database remains unchanged'
            write(*,'(A,I0)') 'active pair records = ',nExpectedPairs
            write(*,'(A,2ES18.8)') 'controlled zeta range = ',1.95D0,2.55D0
            write(*,'(A,F7.3)') 'uniform interior blend fraction = ',0.2D0
            write(*,'(A,ES18.8)') 'minimum controlled quadruplet moles = ',MINVAL(dState)
            if (lDiagnosticsValid) then
                write(*,'(A,ES18.8)') 'legacy unweighted S3 = ',dLegacyS3
                write(*,'(A,ES18.8)') 'corrected weighted S3 = ',dWeightedS3
                write(*,'(A,ES18.8)') 'Delta S3 (weighted-legacy) = ',dDeltaS3
                write(*,'(A,ES18.8)') 'production reference/configurational molar = ', &
                    dProductionReferenceIdeal
                write(*,'(A,ES18.8)') 'corrected standalone reference/configurational molar = ', &
                    dCorrectedMolar
                write(*,'(A,ES18.8)') 'reconstructed legacy molar = ',dGLegacyMolar
                write(*,'(A,ES14.6)') 'corrected reference/configurational parity error = ',dParityError
                write(*,'(A,ES14.6)') 'complete excess-energy parity error = ',dExcessParityError
                write(*,'(A,ES14.6)') 'complete total-energy parity error = ',dTotalParityError
                write(*,'(A,ES14.6)') 'production/legacy scaled separation = ',dLegacySeparation
                write(*,'(A,ES14.6)') 'complete gradient normwise scaled error = ',dGradientError
                write(*,'(A,ES14.6,A,I0)') 'complete gradient worst component scaled error = ', &
                    dGradientMaxScaled,', quadruplet = ',iGradientWorst
                write(*,'(A)') 'quadruplet       analytic gradient       production partial molar       difference'
                do iPair = 1, SIZE(dGradient)
                    write(*,'(I6,3ES26.14)') iPair,dGradient(iPair),dMuProduction(iPair), &
                        dMuProduction(iPair)-dGradient(iPair)
                end do
            else
                write(*,'(A,I0)') 'controlled diagnostics unavailable; evaluator status = ',iInfoLocal
            end if
            write(*,'(A,L1)') 'runtime zeta row restored exactly = ',lRestored
            write(*,'(A,L1)') 'controlled regression pass = ',lCasePass
        end if

        deallocate(dGradient,dHessian,dMuProduction,dState,dZetaSave)

    end subroutine VerifyControlledNonuniformZeta


    !---------------------------------------------------------------------------------------------------------
    !> \brief Independently evaluate legacy-unweighted and corrected-weighted SUBQ S3.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateS3Definitions(tData,dState,dLegacy,dWeighted,iInfoLocal,dMaxPairDifference)

        type(MQMQAModelData), intent(in) :: tData
        real(8), intent(in) :: dState(:)
        real(8), intent(out) :: dLegacy, dWeighted
        integer, intent(out) :: iInfoLocal
        real(8), intent(out), optional :: dMaxPairDifference

        integer :: a, b, i, iPosition, iWeight, j, jPosition, n, nA, nX, q, x, y
        real(8) :: dEquivalentLog, dLegacyPairLog, dN, dOrdinarySum, dWeightedPairLog, dWeightedSum
        real(8), allocatable :: dEquivalent1(:), dEquivalent2(:), dFraction(:)
        real(8), allocatable :: dOrdinary(:,:), dWeightedAmount(:,:), dXOrdinary(:,:), dXWeighted(:,:)

        iInfoLocal = 0
        dLegacy = 0D0
        dWeighted = 0D0
        if (present(dMaxPairDifference)) dMaxPairDifference = HUGE(1D0)
        n = SIZE(dState)
        if ((tData%iModelType /= MQMQA_MODEL_SUBQ) .OR. &
            (n /= SIZE(tData%iQuadruplet,1)) .OR. ANY(dState <= 0D0)) then
            iInfoLocal = 1
            return
        end if
        allocate(dFraction(n),dEquivalent1(tData%nSublattice1), &
            dEquivalent2(tData%nSublattice2),dOrdinary(tData%nSublattice1,tData%nSublattice2), &
            dWeightedAmount(tData%nSublattice1,tData%nSublattice2), &
            dXOrdinary(tData%nSublattice1,tData%nSublattice2), &
            dXWeighted(tData%nSublattice1,tData%nSublattice2))

        dN = SUM(dState)
        dFraction = dState/dN
        dEquivalent1 = 0D0
        dEquivalent2 = 0D0
        dOrdinary = 0D0
        dWeightedAmount = 0D0
        do q = 1, n
            a = tData%iQuadruplet(q,1); b = tData%iQuadruplet(q,2)
            x = tData%iQuadruplet(q,3); y = tData%iQuadruplet(q,4)
            dEquivalent1(a) = dEquivalent1(a)+0.5D0*dFraction(q)
            dEquivalent1(b) = dEquivalent1(b)+0.5D0*dFraction(q)
            dEquivalent2(x) = dEquivalent2(x)+0.5D0*dFraction(q)
            dEquivalent2(y) = dEquivalent2(y)+0.5D0*dFraction(q)
            do i = 1, tData%nSublattice1
                nA = MERGE(1,0,a==i)+MERGE(1,0,b==i)
                do j = 1, tData%nSublattice2
                    nX = MERGE(1,0,x==j)+MERGE(1,0,y==j)
                    dOrdinary(i,j) = dOrdinary(i,j)+dState(q)*DFLOAT(nA*nX)
                    dWeightedAmount(i,j) = dWeightedAmount(i,j)+ &
                        dState(q)*DFLOAT(nA*nX)/tData%dZeta(i,j)
                end do
            end do
        end do
        dOrdinarySum = SUM(dOrdinary)
        dWeightedSum = SUM(dWeightedAmount)
        if ((dOrdinarySum <= 0D0) .OR. (dWeightedSum <= 0D0) .OR. &
            ANY(dEquivalent1 <= 0D0) .OR. ANY(dEquivalent2 <= 0D0)) then
            iInfoLocal = 2
            return
        end if
        dXOrdinary = dOrdinary/dOrdinarySum
        dXWeighted = dWeightedAmount/dWeightedSum
        if (present(dMaxPairDifference)) &
            dMaxPairDifference = MAXVAL(DABS(dXWeighted-dXOrdinary))

        do q = 1, n
            a = tData%iQuadruplet(q,1); b = tData%iQuadruplet(q,2)
            x = tData%iQuadruplet(q,3); y = tData%iQuadruplet(q,4)
            iWeight = 1
            if (a /= b) iWeight = 2*iWeight
            if (x /= y) iWeight = 2*iWeight
            dLegacyPairLog = 0D0
            dWeightedPairLog = 0D0
            do iPosition = 1, 2
                do jPosition = 3, 4
                    i = tData%iQuadruplet(q,iPosition)
                    j = tData%iQuadruplet(q,jPosition)
                    if ((dXOrdinary(i,j) <= 0D0) .OR. (dXWeighted(i,j) <= 0D0)) then
                        iInfoLocal = 3
                        return
                    end if
                    dLegacyPairLog = dLegacyPairLog+DLOG(dXOrdinary(i,j))
                    dWeightedPairLog = dWeightedPairLog+DLOG(dXWeighted(i,j))
                end do
            end do
            dEquivalentLog = DLOG(dEquivalent1(a))+DLOG(dEquivalent1(b))+ &
                DLOG(dEquivalent2(x))+DLOG(dEquivalent2(y))
            dLegacy = dLegacy+dState(q)*(DLOG(dFraction(q))-DLOG(DFLOAT(iWeight)) &
                -0.75D0*dLegacyPairLog+0.5D0*dEquivalentLog)
            dWeighted = dWeighted+dState(q)*(DLOG(dFraction(q))-DLOG(DFLOAT(iWeight)) &
                -0.75D0*dWeightedPairLog+0.5D0*dEquivalentLog)
        end do

    end subroutine EvaluateS3Definitions


    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate production SUBQ scalar blocks at a temporary positive mole state.
    !>
    !> \details The established SUBG production routine also evaluates SUBQ.
    !!          Every global vector it modifies is restored before returning, so
    !!          this local comparison cannot alter the converged calculation.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateProductionEnergy(iPhaseLocal,dMolesLocal,dReferenceIdeal,dExcess,iInfoLocal,dMu)

        integer, intent(in) :: iPhaseLocal
        real(8), intent(in) :: dMolesLocal(:)
        real(8), intent(out) :: dReferenceIdeal, dExcess
        integer, intent(out) :: iInfoLocal
        real(8), intent(out), optional :: dMu(:)

        integer :: iFirstLocal, iInfoSave, iLastLocal
        real(8), allocatable :: dChemicalSave(:), dFractionSave(:), dPartialSave(:), dX(:)

        iInfoLocal = 0
        iFirstLocal = nSpeciesPhase(iPhaseLocal-1)+1
        iLastLocal = nSpeciesPhase(iPhaseLocal)
        if ((SIZE(dMolesLocal) /= iLastLocal-iFirstLocal+1) .OR. &
            (SUM(dMolesLocal) <= 0D0) .OR. ANY(dMolesLocal <= 0D0)) then
            iInfoLocal = 1
            return
        end if
        if (PRESENT(dMu)) then
            if (SIZE(dMu) /= SIZE(dMolesLocal)) then
                iInfoLocal = 1
                return
            end if
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
        if (PRESENT(dMu)) dMu = dChemicalPotential(iFirstLocal:iLastLocal) + &
            dPartialExcessGibbs(iFirstLocal:iLastLocal)
        if (INFOThermo /= iInfoSave) iInfoLocal = 2

        dChemicalPotential = dChemicalSave
        dMolFraction = dFractionSave
        dPartialExcessGibbs = dPartialSave
        INFOThermo = iInfoSave
        deallocate(dChemicalSave,dFractionSave,dPartialSave,dX)

    end subroutine EvaluateProductionEnergy


    real(8) function ScaledError(dA,dB)

        real(8), intent(in) :: dA, dB
        ScaledError = DABS(dA-dB)/DMAX1(1D0,DABS(dA),DABS(dB))

    end function ScaledError

end program TestMQMQASUBQNativeEnergyVerification
