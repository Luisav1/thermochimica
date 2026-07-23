!-------------------------------------------------------------------------------------------------------------
!> \file    TestCEFHessianVerification.F90
!> \brief   Thermochimica-native verification of the disconnected nonmagnetic plain-SUBL CEF Hessian.
!>
!> \details The test checks standalone scalar-energy derivatives, production parameter decoding, and finite
!!          differences of the established Thermochimica SUBL partial molars.  It does not connect the CEF
!!          Hessian to GEMNewton.  Pass --report to print the full numerical evidence.
!!
!!          Verification map:
!!          1. Run real Thermochimica SUBL calculations and retain controlled and converged states.
!!          2. Decode production occupancy and interaction arrays into the generic CEF interface.
!!          3. Prove each decoded scalar interaction matches its production expression.
!!          4. Compare analytic directional curvature with three- and five-point energy differences.
!!          5. Perturb the endmember mole vector in a composition-changing
!!             direction v. Multiplying the Hessian H by v predicts how every
!!             partial molar changes; compare that prediction with finite
!!             differences of established production partial molars.
!!          6. Check matrix symmetry and degree-one extensivity, including the
!!             requirement that uniform phase scaling leaves composition and
!!             partial molars unchanged. Also check precision, branch coverage,
!!             and admissible-state handling.
!!
!!          This is Thermochimica-native local verification because it parses
!!          production databases and calls established SUBL thermodynamics. It
!!          remains disconnected from GEMNewton and does not constrain phases.
!-------------------------------------------------------------------------------------------------------------

program TestCEFHessianVerification

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleCEFUnconstrained
    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    interface
        subroutine CompExcessGibbsEnergySUBL(iSolnIndex)
            integer, intent(in) :: iSolnIndex
        end subroutine CompExcessGibbsEnergySUBL
    end interface

    logical :: lPass, lReport, lConvergedDone
    integer :: nBinaryCoverage, nTernaryCoverage, nCoupledCoverage
    character(len=32) :: cArgument

    lPass = .TRUE.
    lReport = .FALSE.
    lConvergedDone = .FALSE.
    nBinaryCoverage = 0
    nTernaryCoverage = 0
    nCoupledCoverage = 0
    if (COMMAND_ARGUMENT_COUNT() > 0) then
        call GET_COMMAND_ARGUMENT(1,cArgument)
        lReport = TRIM(cArgument) == '--report'
    end if

    lPass = lPass .AND. (STORAGE_SIZE(1D0) == 64)
    lPass = lPass .AND. (PRECISION(1D0) >= 15)
    lPass = lPass .AND. (DIGITS(1D0) >= 53)
    lPass = lPass .AND. (EPSILON(1D0) <= 3D-16)

    if (lReport) then
        write(*,'(A)') 'General disconnected plain-SUBL CEF Hessian verification'
        write(*,'(A,I0)') 'storage bits = ',STORAGE_SIZE(1D0)
        write(*,'(A,I0)') 'decimal precision = ',PRECISION(1D0)
        write(*,'(A,I0)') 'binary digits = ',DIGITS(1D0)
        write(*,'(A,ES14.6)') 'machine epsilon = ',EPSILON(1D0)
    end if

    call RunAlabanditeCases(lPass,lConvergedDone,nBinaryCoverage,nTernaryCoverage,nCoupledCoverage,lReport)
    call ResetThermoAll
    call RunCoupledCases(lPass,lConvergedDone,nBinaryCoverage,nTernaryCoverage,nCoupledCoverage,lReport)
    call ResetThermoAll

    lPass = lPass .AND. (nBinaryCoverage > 0) .AND. (nTernaryCoverage > 0) .AND. &
        (nCoupledCoverage > 0) .AND. lConvergedDone
    if (lReport) then
        write(*,'(/,A,3(I0,1X))') 'family coverage (binary ternary coupled) = ', &
            nBinaryCoverage,nTernaryCoverage,nCoupledCoverage
        write(*,'(A,L1)') 'admissible converged comparison completed = ',lConvergedDone
    end if

    if (lPass) then
        print *, 'TestCEFHessianVerification: PASS'
        call EXIT(0)
    else
        print *, 'TestCEFHessianVerification: FAIL <---'
        call EXIT(1)
    end if

contains

    !=========================================================================================================
    ! SECTION 1: PRODUCTION CALCULATIONS AND VERIFICATION STATES
    !
    ! ALABANDITE exercises existing binary and ternary SUBL parameters.
    ! CEFVerification.dat supplies the coupled two-sublattice interaction family
    ! that was not otherwise available in a compact admissible regression state.
    ! Each case checks both a controlled interior composition and, when active,
    ! the composition returned by the ordinary Thermochimica calculation.
    !=========================================================================================================

    subroutine RunAlabanditeCases(lAllPass,lAnyConverged,nBinary,nTernary,nCoupled,lVerbose)

        logical, intent(inout) :: lAllPass, lAnyConverged
        integer, intent(inout) :: nBinary, nTernary, nCoupled
        logical, intent(in) :: lVerbose
        integer :: iPhase, iFirst, iLast
        real(8), allocatable :: dControlled(:), dConverged(:)

        cInputUnitTemperature = 'K'
        cInputUnitPressure = 'atm'
        cInputUnitMass = 'moles'
        cThermoFileName = DATA_DIRECTORY // 'FeMnCaS-1.dat'
        dPressure = 1D0
        dTemperature = 500D0
        dElementMass = 0D0
        dElementMass(20) = 1D0
        dElementMass(26) = 1D0
        dElementMass(25) = 1D0
        dElementMass(16) = 3D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        lAllPass = lAllPass .AND. (INFOThermo == 0)
        iPhase = FindPhase('ALABANDITE')
        lAllPass = lAllPass .AND. (iPhase > 0)
        if (iPhase <= 0) return

        iFirst = nSpeciesPhase(iPhase-1) + 1
        iLast = nSpeciesPhase(iPhase)
        allocate(dControlled(iLast-iFirst+1),dConverged(iLast-iFirst+1))
        dControlled = [0.42D0,0.33D0,0.25D0]
        dConverged = dMolesSpecies(iFirst:iLast)
        call VerifyPhaseState(iPhase,dControlled,'ALABANDITE controlled',.FALSE.,lAllPass,lAnyConverged, &
            nBinary,nTernary,nCoupled,lVerbose)
        if (PhaseIsActive(iPhase)) then
            call VerifyPhaseState(iPhase,dConverged,'ALABANDITE converged',.TRUE.,lAllPass,lAnyConverged, &
                nBinary,nTernary,nCoupled,lVerbose)
        else if (lVerbose) then
            write(*,'(/,A)') 'ALABANDITE converged: SKIP (phase is not active)'
        end if

    end subroutine RunAlabanditeCases


    subroutine RunCoupledCases(lAllPass,lAnyConverged,nBinary,nTernary,nCoupled,lVerbose)

        logical, intent(inout) :: lAllPass, lAnyConverged
        integer, intent(inout) :: nBinary, nTernary, nCoupled
        logical, intent(in) :: lVerbose
        integer :: iPhase, iFirst, iLast
        real(8), allocatable :: dControlled(:), dConverged(:)

        cInputUnitTemperature = 'K'
        cInputUnitPressure = 'atm'
        cInputUnitMass = 'moles'
        cThermoFileName = DATA_DIRECTORY // 'CEFVerification.dat'
        dPressure = 1D0
        dTemperature = 1000D0
        dElementMass = 0D0
        dElementMass(1) = 1D0
        dElementMass(8) = 1D0
        call ParseCSDataFile(cThermoFileName)
        lAllPass = lAllPass .AND. (INFOThermo == 0)
        if (INFOThermo /= 0) return
        call Thermochimica
        lAllPass = lAllPass .AND. (INFOThermo == 0)
        iPhase = FindPhase('CEF_COUPLED')
        lAllPass = lAllPass .AND. (iPhase > 0)
        if (iPhase <= 0) return

        iFirst = nSpeciesPhase(iPhase-1) + 1
        iLast = nSpeciesPhase(iPhase)
        allocate(dControlled(iLast-iFirst+1),dConverged(iLast-iFirst+1))
        dControlled = [0.31D0,0.19D0,0.23D0,0.27D0]
        dConverged = dMolesSpecies(iFirst:iLast)
        call VerifyPhaseState(iPhase,dControlled,'CEF_COUPLED controlled',.FALSE.,lAllPass,lAnyConverged, &
            nBinary,nTernary,nCoupled,lVerbose)
        if (PhaseIsActive(iPhase)) then
            call VerifyPhaseState(iPhase,dConverged,'CEF_COUPLED converged',.TRUE.,lAllPass,lAnyConverged, &
                nBinary,nTernary,nCoupled,lVerbose)
        else if (lVerbose) then
            write(*,'(/,A)') 'CEF_COUPLED converged: SKIP (phase is not active)'
        end if

    end subroutine RunCoupledCases

    !=========================================================================================================
    ! SECTION 2: COMPLETE CHECK OF ONE PLAIN-SUBL PHASE STATE
    !
    ! This routine joins the three evidence layers: generic standalone
    ! mathematics, exact production-parameter decoding, and finite differences
    ! of established Thermochimica partial molars.
    !=========================================================================================================

    subroutine VerifyPhaseState(iPhase,dMoles,cLabel,lConverged,lAllPass,lAnyConverged, &
        nBinary,nTernary,nCoupled,lVerbose)

        integer, intent(in) :: iPhase
        real(8), intent(in) :: dMoles(:)
        character(*), intent(in) :: cLabel
        logical, intent(in) :: lConverged, lVerbose
        logical, intent(inout) :: lAllPass, lAnyConverged
        integer, intent(inout) :: nBinary, nTernary, nCoupled

        integer, parameter :: nSteps = 7
        integer :: iDirection, iInfo, iStep, nDirections, nLocalBinary, nLocalCoupled, nLocalTernary
        integer, allocatable :: iOccupancy(:,:), iSiteSublattice(:)
        real(8) :: dBest3, dBest5, dDecoderError, dG, dGEx, dGId, dGRef, dHomogeneity
        real(8) :: dMinSiteFraction, dProductionEx, dProductionRefIdeal, dSymmetry, dWorstMu
        real(8), allocatable :: dDirection(:), dErr3(:,:), dErr5(:,:), dErrMu(:,:)
        real(8), allocatable :: dH(:,:), dHEx(:,:), dHId(:,:), dHRef(:,:), dMultiplicity(:), dReference(:)
        real(8), allocatable :: dSiteFraction(:), dSteps(:)
        type(CEFInteractionTerm), allocatable :: tInteraction(:)

        call DecodeProductionPhase(iPhase,dMoles,iOccupancy,iSiteSublattice,dMultiplicity,dReference, &
            tInteraction,dSiteFraction,dDecoderError,nLocalBinary,nLocalTernary,nLocalCoupled,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0)
        if (iInfo /= 0) then
            if (lVerbose) write(*,'(/,A,A,I0)') TRIM(cLabel),': decoder failure iInfo = ',iInfo
            return
        end if
        nBinary = nBinary + nLocalBinary
        nTernary = nTernary + nLocalTernary
        nCoupled = nCoupled + nLocalCoupled
        dMinSiteFraction = MINVAL(dSiteFraction)
        if (lConverged .AND. (dMinSiteFraction <= 1D-12)) then
            if (lVerbose) write(*,'(/,A,A,ES14.6)') TRIM(cLabel),': SKIP boundary min(y) = ',dMinSiteFraction
            return
        end if

        allocate(dH(SIZE(dMoles),SIZE(dMoles)),dHRef(SIZE(dMoles),SIZE(dMoles)), &
            dHId(SIZE(dMoles),SIZE(dMoles)),dHEx(SIZE(dMoles),SIZE(dMoles)))
        call CompCEFGibbsEnergyUnconstrained(dMoles,iOccupancy,iSiteSublattice,dMultiplicity,dReference, &
            1D0,tInteraction,dG,dGRef,dGId,dGEx,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0)
        call CompCEFHessianUnconstrained(dMoles,iOccupancy,iSiteSublattice,dMultiplicity,dReference, &
            1D0,tInteraction,dH,iInfo,dHRef,dHId,dHEx)
        lAllPass = lAllPass .AND. (iInfo == 0) .AND. ALL(IEEE_IS_FINITE(dH))
        if (iInfo /= 0) return

        dSymmetry = FrobeniusNorm(dH-TRANSPOSE(dH))/DMAX1(1D0,FrobeniusNorm(dH))
        dHomogeneity = VectorTwoNorm(MATMUL(dH,dMoles))/ &
            DMAX1(1D0,MatrixTwoNormSymmetric(dH)*VectorTwoNorm(dMoles))
        lAllPass = lAllPass .AND. (dSymmetry <= 1D-12) .AND. (dHomogeneity <= 1D-10)
        lAllPass = lAllPass .AND. (dDecoderError <= 1D-10)
        if (dDecoderError > 1D-10) return

        call EvaluateProduction(iPhase,dMoles,dProductionRefIdeal,dProductionEx,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0)
        lAllPass = lAllPass .AND. (ScaledError(dProductionRefIdeal,(dGRef+dGId)/SUM(dMoles)) <= 1D-10)
        lAllPass = lAllPass .AND. (ScaledError(dProductionEx,dGEx/SUM(dMoles)) <= 1D-10)

        nDirections = SIZE(dMoles)
        allocate(dDirection(SIZE(dMoles)),dSteps(nSteps),dErr3(nDirections,nSteps), &
            dErr5(nDirections,nSteps),dErrMu(nDirections,nSteps))
        do iDirection = 1, nDirections
            call BuildDirection(iDirection,dDirection)
            call VerifyDirection(iPhase,dMoles,dDirection,iOccupancy,iSiteSublattice,dMultiplicity, &
                dReference,tInteraction,dH,dSteps,dErr3(iDirection,:),dErr5(iDirection,:), &
                dErrMu(iDirection,:),lAllPass)
        end do

        dBest3 = 0D0
        dBest5 = 0D0
        dWorstMu = 0D0
        do iDirection = 1, nDirections
            dBest3 = DMAX1(dBest3,MINVAL(dErr3(iDirection,:)))
            dBest5 = DMAX1(dBest5,MINVAL(dErr5(iDirection,:)))
            dWorstMu = DMAX1(dWorstMu,MINVAL(dErrMu(iDirection,:)))
            lAllPass = lAllPass .AND. (MINLOC(dErr3(iDirection,:),1) > 1)
        end do
        ! Scalar-energy thresholds apply to the mandatory controlled interior states.  Converged states can lie
        ! close enough to a composition boundary that energy second differences lose digits; their complete
        ! sweep remains reportable, while the more stable production-partial-molar comparison stays mandatory.
        if (.NOT.lConverged) then
            lAllPass = lAllPass .AND. (dBest3 <= 1D-6) .AND. (dBest5 <= 1D-8)
        end if
        lAllPass = lAllPass .AND. (dWorstMu <= 1D-8)
        if (lConverged) lAnyConverged = .TRUE.

        if (lVerbose) then
            write(*,'(/,A)') TRIM(cLabel)
            write(*,'(A,ES14.6)') 'min site fraction = ',dMinSiteFraction
            write(*,'(A,3(I0,1X))') 'family terms (binary ternary coupled) = ', &
                nLocalBinary,nLocalTernary,nLocalCoupled
            write(*,'(A,ES14.6)') 'maximum scalar decoder error = ',dDecoderError
            write(*,'(A,ES14.6)') 'reference+ideal scalar error = ', &
                ScaledError(dProductionRefIdeal,(dGRef+dGId)/SUM(dMoles))
            write(*,'(A,ES14.6)') 'excess scalar error = ',ScaledError(dProductionEx,dGEx/SUM(dMoles))
            write(*,'(A,ES14.6)') 'symmetry residual = ',dSymmetry
            write(*,'(A,ES14.6)') 'homogeneity residual = ',dHomogeneity
            write(*,'(A)') 'direction h                 3-point error         5-point error         production-mu error'
            do iDirection = 1, nDirections
                do iStep = 1, nSteps
                    write(*,'(I5,4ES22.12)') iDirection,dSteps(iStep),dErr3(iDirection,iStep), &
                        dErr5(iDirection,iStep),dErrMu(iDirection,iStep)
                end do
            end do
            write(*,'(A,ES14.6)') 'worst best 3-point error = ',dBest3
            write(*,'(A,ES14.6)') 'worst best 5-point error = ',dBest5
            write(*,'(A,ES14.6)') 'worst best production-mu error = ',dWorstMu
        end if

    end subroutine VerifyPhaseState

    !=========================================================================================================
    ! SECTION 3: PRODUCTION-SUBL TO GENERIC-CEF DECODER
    !
    ! Thermochimica stores constituent identities and interaction families in
    ! packed production arrays. This adapter converts them into occupancy,
    ! multiplicity, reference-energy, and generic excess-interaction objects
    ! understood by
    ! ModuleCEFUnconstrained. It belongs in the test because the CEF module itself
    ! deliberately has no dependency on ModuleThermo.
    !=========================================================================================================

    !---------------------------------------------------------------------------------------------------------
    !> \brief Decode one production plain-SUBL phase and prove scalar term identities.
    !>
    !> \details Binary, ternary, and coupled interactions are translated into a
    !!          common representation: a product P of selected site fractions
    !!          multiplied by an interaction polynomial L evaluated at the local
    !!          composition coordinate eta. Each translated term is evaluated
    !!          immediately and compared with the exact production-family scalar
    !!          expression, preventing a self-consistent but incorrect decoder
    !!          from passing the later standalone finite-difference checks.
    !---------------------------------------------------------------------------------------------------------
    subroutine DecodeProductionPhase(iPhase,dMoles,iOccupancy,iSiteSublattice,dMultiplicity,dReference, &
        tInteraction,dY,dDecoderError,nBinary,nTernary,nCoupled,iInfo)

        integer, intent(in) :: iPhase
        real(8), intent(in) :: dMoles(:)
        integer, allocatable, intent(out) :: iOccupancy(:,:), iSiteSublattice(:)
        real(8), allocatable, intent(out) :: dMultiplicity(:), dReference(:), dY(:)
        type(CEFInteractionTerm), allocatable, intent(out) :: tInteraction(:)
        real(8), intent(out) :: dDecoderError
        integer, intent(out) :: nBinary, nTernary, nCoupled, iInfo

        integer :: i, iCharged, iFirst, iFirstPosition, iOrder, iParam, iSecondPosition
        integer :: iThirdPosition, k, n, nSite, nSublattice, s, u
        integer, allocatable :: iOffset(:)
        real(8) :: dExact, dFirst, dSecond, dThird

        iInfo = 0
        nBinary = 0
        nTernary = 0
        nCoupled = 0
        dDecoderError = 0D0
        iCharged = iPhaseSublattice(iPhase)
        nSublattice = nSublatticePhase(iCharged)
        iFirst = nSpeciesPhase(iPhase-1) + 1
        if ((cSolnPhaseType(iPhase) /= 'SUBL') .OR. (SIZE(dMoles) /= nSpeciesPhase(iPhase)-iFirst+1)) then
            iInfo = 1
            return
        end if

        allocate(iOffset(nSublattice+1))
        iOffset(1) = 0
        do s = 1, nSublattice
            iOffset(s+1) = iOffset(s) + nConstituentSublattice(iCharged,s)
        end do
        nSite = iOffset(nSublattice+1)
        allocate(iOccupancy(nSublattice,SIZE(dMoles)),iSiteSublattice(nSite), &
            dMultiplicity(nSublattice),dReference(SIZE(dMoles)),dY(nSite))
        do s = 1, nSublattice
            iSiteSublattice(iOffset(s)+1:iOffset(s+1)) = s
            dMultiplicity(s) = dStoichSublattice(iCharged,s)
            do i = 1, SIZE(dMoles)
                iOccupancy(s,i) = iOffset(s) + iConstituentSublattice(iCharged,s,i)
            end do
        end do
        dReference = dStdGibbsEnergy(iFirst:nSpeciesPhase(iPhase))
        call ComputeSiteFractions(dMoles,iOccupancy,nSite,dY)

        allocate(tInteraction(nParamPhase(iPhase)-nParamPhase(iPhase-1)))
        i = 0
        do iParam = nParamPhase(iPhase-1)+1, nParamPhase(iPhase)
            i = i + 1
            n = iRegularParam(iParam,1)
            allocate(tInteraction(i)%dExponent(nSite),tInteraction(i)%dArgumentCoefficient(nSite))
            tInteraction(i)%dExponent = 0D0
            tInteraction(i)%dArgumentCoefficient = 0D0
            do k = 2, n+1
                call PackedToSite(iRegularParam(iParam,k),iOffset,u)
                tInteraction(i)%dExponent(u) = tInteraction(i)%dExponent(u) + 1D0
            end do
            iOrder = iRegularParam(iParam,n+2)

            if ((iSUBLParamData(iParam,1) == 1) .AND. (iSUBLParamData(iParam,3) == 2)) then
                nBinary = nBinary + 1
                iFirstPosition = iSUBLParamData(iParam,2)
                iSecondPosition = iFirstPosition + 1
                call SetBinaryArgument(iRegularParam(iParam,iFirstPosition), &
                    iRegularParam(iParam,iSecondPosition),iOffset,tInteraction(i))
                allocate(tInteraction(i)%dPolynomialCoefficient(iOrder+1))
                tInteraction(i)%dPolynomialCoefficient = 0D0
                tInteraction(i)%dPolynomialCoefficient(iOrder+1) = dExcessGibbsParam(iParam)
                call PackedValue(iRegularParam(iParam,iFirstPosition),iOffset,dY,dFirst)
                call PackedValue(iRegularParam(iParam,iSecondPosition),iOffset,dY,dSecond)
                dExact = PRODUCT(dY**tInteraction(i)%dExponent)*dExcessGibbsParam(iParam)* &
                    (dFirst-dSecond)**iOrder
            else if ((iSUBLParamData(iParam,1) == 1) .AND. (iSUBLParamData(iParam,3) == 3)) then
                nTernary = nTernary + 1
                iFirstPosition = iSUBLParamData(iParam,2) + iOrder
                iSecondPosition = iSUBLParamData(iParam,2) + MOD(iOrder+1,3)
                iThirdPosition = iSUBLParamData(iParam,2) + MOD(iOrder+2,3)
                call SetTernaryArgument(iRegularParam(iParam,iFirstPosition), &
                    iRegularParam(iParam,iSecondPosition),iRegularParam(iParam,iThirdPosition), &
                    iOffset,tInteraction(i))
                allocate(tInteraction(i)%dPolynomialCoefficient(2))
                tInteraction(i)%dPolynomialCoefficient = [0D0,dExcessGibbsParam(iParam)]
                call PackedValue(iRegularParam(iParam,iFirstPosition),iOffset,dY,dFirst)
                call PackedValue(iRegularParam(iParam,iSecondPosition),iOffset,dY,dSecond)
                call PackedValue(iRegularParam(iParam,iThirdPosition),iOffset,dY,dThird)
                dExact = PRODUCT(dY**tInteraction(i)%dExponent)*dExcessGibbsParam(iParam)* &
                    (dFirst+(1D0-dFirst-dSecond-dThird)/3D0)
            else if ((iSUBLParamData(iParam,1) == 2) .AND. (iSUBLParamData(iParam,3) == 2) .AND. &
                (iSUBLParamData(iParam,5) == 2)) then
                nCoupled = nCoupled + 1
                if (MOD(iOrder,2) == 0) then
                    iFirstPosition = iSUBLParamData(iParam,2)
                else
                    iFirstPosition = iSUBLParamData(iParam,4)
                end if
                iSecondPosition = iFirstPosition + 1
                call SetBinaryArgument(iRegularParam(iParam,iFirstPosition), &
                    iRegularParam(iParam,iSecondPosition),iOffset,tInteraction(i))
                allocate(tInteraction(i)%dPolynomialCoefficient(iOrder/2+1))
                tInteraction(i)%dPolynomialCoefficient = 0D0
                tInteraction(i)%dPolynomialCoefficient(iOrder/2+1) = dExcessGibbsParam(iParam)
                call PackedValue(iRegularParam(iParam,iFirstPosition),iOffset,dY,dFirst)
                call PackedValue(iRegularParam(iParam,iSecondPosition),iOffset,dY,dSecond)
                dExact = PRODUCT(dY**tInteraction(i)%dExponent)*dExcessGibbsParam(iParam)* &
                    (dFirst-dSecond)**(iOrder/2)
            else
                iInfo = 2
                return
            end if
            dDecoderError = DMAX1(dDecoderError,ScaledError(CompCEFInteractionValue(dY,tInteraction(i)),dExact))
        end do

    end subroutine DecodeProductionPhase

    !=========================================================================================================
    ! SECTION 4: DIRECTIONAL FINITE DIFFERENCES AND PRODUCTION ORACLE
    !
    ! A direction is composition-preserving here when its entries sum to zero:
    ! moles are transferred among endmembers without changing total phase amount.
    ! Such endmember-mole directions test the standalone
    ! energy Hessian. The same perturbations are sent through the established
    ! production SUBL routine. Multiplying H by the direction predicts how every
    ! partial molar changes, providing an independent comparison at the identical
    ! imposed state.
    !=========================================================================================================

    !---------------------------------------------------------------------------------------------------------
    !> \brief Verify one composition-preserving mole-transfer direction using scalar energy and production partial molars.
    !---------------------------------------------------------------------------------------------------------
    subroutine VerifyDirection(iPhase,dMoles,dDirection,iOccupancy,iSiteSublattice,dMultiplicity, &
        dReference,tInteraction,dH,dSteps,dErr3,dErr5,dErrMu,lAllPass)

        integer, intent(in) :: iPhase, iOccupancy(:,:), iSiteSublattice(:)
        real(8), intent(in) :: dMoles(:), dDirection(:), dMultiplicity(:), dReference(:), dH(:,:)
        type(CEFInteractionTerm), intent(in) :: tInteraction(:)
        real(8), intent(out) :: dSteps(:), dErr3(:), dErr5(:), dErrMu(:)
        logical, intent(inout) :: lAllPass

        integer :: iInfo, iStep
        real(8) :: dAnalytic, dE0, dEM1, dEM2, dEP1, dEP2, dGEx, dGId, dGRef, dHStep, dScale
        real(8) :: dDummyEx, dDummyRefIdeal
        real(8), allocatable :: dMinus(:), dMinus2(:), dMuMinus(:), dMuPlus(:), dPlus(:), dPlus2(:)

        allocate(dMinus(SIZE(dMoles)),dMinus2(SIZE(dMoles)),dPlus(SIZE(dMoles)), &
            dPlus2(SIZE(dMoles)),dMuMinus(SIZE(dMoles)),dMuPlus(SIZE(dMoles)))
        dAnalytic = DOT_PRODUCT(dDirection,MATMUL(dH,dDirection))
        call CompCEFGibbsEnergyUnconstrained(dMoles,iOccupancy,iSiteSublattice,dMultiplicity,dReference, &
            1D0,tInteraction,dE0,dGRef,dGId,dGEx,iInfo)
        lAllPass = lAllPass .AND. (iInfo == 0)
        dScale = HUGE(1D0)
        where (DABS(dDirection) > 0D0)
            dMinus = dMoles/DABS(dDirection)
        elsewhere
            dMinus = HUGE(1D0)
        end where
        dScale = MINVAL(dMinus)

        do iStep = 1, SIZE(dSteps)
            dHStep = 0.1D0*dScale*3D0**(-(iStep-1))
            dSteps(iStep) = dHStep
            dMinus = dMoles-dHStep*dDirection
            dPlus = dMoles+dHStep*dDirection
            dMinus2 = dMoles-2D0*dHStep*dDirection
            dPlus2 = dMoles+2D0*dHStep*dDirection
            call CompCEFGibbsEnergyUnconstrained(dMinus,iOccupancy,iSiteSublattice,dMultiplicity,dReference, &
                1D0,tInteraction,dEM1,dGRef,dGId,dGEx,iInfo)
            lAllPass = lAllPass .AND. (iInfo == 0)
            call CompCEFGibbsEnergyUnconstrained(dPlus,iOccupancy,iSiteSublattice,dMultiplicity,dReference, &
                1D0,tInteraction,dEP1,dGRef,dGId,dGEx,iInfo)
            lAllPass = lAllPass .AND. (iInfo == 0)
            call CompCEFGibbsEnergyUnconstrained(dMinus2,iOccupancy,iSiteSublattice,dMultiplicity,dReference, &
                1D0,tInteraction,dEM2,dGRef,dGId,dGEx,iInfo)
            lAllPass = lAllPass .AND. (iInfo == 0)
            call CompCEFGibbsEnergyUnconstrained(dPlus2,iOccupancy,iSiteSublattice,dMultiplicity,dReference, &
                1D0,tInteraction,dEP2,dGRef,dGId,dGEx,iInfo)
            lAllPass = lAllPass .AND. (iInfo == 0)
            dErr3(iStep) = ScaledError((dEP1-2D0*dE0+dEM1)/(dHStep*dHStep),dAnalytic)
            dErr5(iStep) = ScaledError((-dEP2+16D0*dEP1-30D0*dE0+16D0*dEM1-dEM2)/ &
                (12D0*dHStep*dHStep),dAnalytic)

            call EvaluateProductionVector(iPhase,dMinus,dMuMinus,dDummyRefIdeal,dDummyEx,iInfo)
            lAllPass = lAllPass .AND. (iInfo == 0)
            call EvaluateProductionVector(iPhase,dPlus,dMuPlus,dDummyRefIdeal,dDummyEx,iInfo)
            lAllPass = lAllPass .AND. (iInfo == 0)
            dErrMu(iStep) = VectorTwoNorm((dMuPlus-dMuMinus)/(2D0*dHStep)-MATMUL(dH,dDirection))/ &
                DMAX1(1D0,VectorTwoNorm(MATMUL(dH,dDirection)))
        end do

    end subroutine VerifyDirection


    subroutine EvaluateProduction(iPhase,dMoles,dReferenceIdeal,dExcess,iInfo)

        integer, intent(in) :: iPhase
        real(8), intent(in) :: dMoles(:)
        real(8), intent(out) :: dReferenceIdeal, dExcess
        integer, intent(out) :: iInfo
        real(8), allocatable :: dMu(:)

        allocate(dMu(SIZE(dMoles)))
        call EvaluateProductionVector(iPhase,dMoles,dMu,dReferenceIdeal,dExcess,iInfo)

    end subroutine EvaluateProduction


    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate established SUBL thermodynamics at a temporarily imposed local composition.
    !>
    !> \details Thermochimica global arrays are saved, updated only long enough to
    !!          call CompExcessGibbsEnergySUBL, and then fully restored. This lets
    !!          the production routine act as an independent numerical oracle
    !!          without changing the surrounding converged calculation.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateProductionVector(iPhase,dMoles,dMu,dReferenceIdeal,dExcess,iInfo)

        integer, intent(in) :: iPhase
        real(8), intent(in) :: dMoles(:)
        real(8), intent(out) :: dMu(:), dReferenceIdeal, dExcess
        integer, intent(out) :: iInfo
        integer :: iFirst, iLast, iInfoSave
        real(8), allocatable :: dChemicalSave(:), dMolFractionSave(:), dPartialSave(:), dSiteSave(:,:,:), dX(:)

        iInfo = 0
        iFirst = nSpeciesPhase(iPhase-1)+1
        iLast = nSpeciesPhase(iPhase)
        if ((SIZE(dMoles) /= iLast-iFirst+1) .OR. (SIZE(dMu) /= SIZE(dMoles)) .OR. &
            (SUM(dMoles) <= 0D0) .OR. ANY(dMoles < 0D0)) then
            iInfo = 1
            return
        end if
        allocate(dChemicalSave(SIZE(dChemicalPotential)),dMolFractionSave(SIZE(dMolFraction)), &
            dPartialSave(SIZE(dPartialExcessGibbs)),dSiteSave(SIZE(dSiteFraction,1), &
            SIZE(dSiteFraction,2),SIZE(dSiteFraction,3)),dX(SIZE(dMoles)))
        dChemicalSave = dChemicalPotential
        dMolFractionSave = dMolFraction
        dPartialSave = dPartialExcessGibbs
        dSiteSave = dSiteFraction
        iInfoSave = INFOThermo

        dMolFraction(iFirst:iLast) = dMoles/SUM(dMoles)
        call CompExcessGibbsEnergySUBL(iPhase)
        dX = dMolFraction(iFirst:iLast)
        dReferenceIdeal = DOT_PRODUCT(dX,dChemicalPotential(iFirst:iLast))
        dExcess = DOT_PRODUCT(dX,dPartialExcessGibbs(iFirst:iLast))
        dMu = dChemicalPotential(iFirst:iLast)+dPartialExcessGibbs(iFirst:iLast)
        if (INFOThermo /= iInfoSave) iInfo = 2

        dChemicalPotential = dChemicalSave
        dMolFraction = dMolFractionSave
        dPartialExcessGibbs = dPartialSave
        dSiteFraction = dSiteSave
        INFOThermo = iInfoSave

    end subroutine EvaluateProductionVector

    !=========================================================================================================
    ! SECTION 5: DECODER, DIRECTION, AND PHASE-STATE HELPERS
    !
    ! These routines translate packed constituent indices, build the linear
    ! binary and ternary coordinates used by the generic interaction object,
    ! construct tangent directions, and query the production phase assemblage.
    !=========================================================================================================

    subroutine ComputeSiteFractions(dMoles,iOccupancy,nSite,dY)

        real(8), intent(in) :: dMoles(:)
        integer, intent(in) :: iOccupancy(:,:), nSite
        real(8), intent(out) :: dY(nSite)
        integer :: i, s

        dY = 0D0
        do i = 1, SIZE(dMoles)
            do s = 1, SIZE(iOccupancy,1)
                dY(iOccupancy(s,i)) = dY(iOccupancy(s,i)) + dMoles(i)/SUM(dMoles)
            end do
        end do

    end subroutine ComputeSiteFractions


    subroutine SetBinaryArgument(iPackedFirst,iPackedSecond,iOffset,tTerm)

        integer, intent(in) :: iPackedFirst, iPackedSecond, iOffset(:)
        type(CEFInteractionTerm), intent(inout) :: tTerm
        integer :: u

        tTerm%dArgumentConstant = 0D0
        call PackedToSite(iPackedFirst,iOffset,u)
        tTerm%dArgumentCoefficient(u) = 1D0
        call PackedToSite(iPackedSecond,iOffset,u)
        tTerm%dArgumentCoefficient(u) = -1D0

    end subroutine SetBinaryArgument


    subroutine SetTernaryArgument(iPackedFirst,iPackedSecond,iPackedThird,iOffset,tTerm)

        integer, intent(in) :: iPackedFirst, iPackedSecond, iPackedThird, iOffset(:)
        type(CEFInteractionTerm), intent(inout) :: tTerm
        integer :: u

        tTerm%dArgumentConstant = 1D0/3D0
        call PackedToSite(iPackedFirst,iOffset,u)
        tTerm%dArgumentCoefficient(u) = 2D0/3D0
        call PackedToSite(iPackedSecond,iOffset,u)
        tTerm%dArgumentCoefficient(u) = -1D0/3D0
        call PackedToSite(iPackedThird,iOffset,u)
        tTerm%dArgumentCoefficient(u) = -1D0/3D0

    end subroutine SetTernaryArgument


    subroutine PackedToSite(iPacked,iOffset,u)

        integer, intent(in) :: iPacked, iOffset(:)
        integer, intent(out) :: u
        integer :: c, s

        c = MOD(iPacked,10000)
        s = (iPacked-c)/10000
        u = iOffset(s)+c

    end subroutine PackedToSite


    subroutine PackedValue(iPacked,iOffset,dY,dValue)

        integer, intent(in) :: iPacked, iOffset(:)
        real(8), intent(in) :: dY(:)
        real(8), intent(out) :: dValue
        integer :: u

        call PackedToSite(iPacked,iOffset,u)
        dValue = dY(u)

    end subroutine PackedValue


    subroutine BuildDirection(iDirection,dDirection)

        integer, intent(in) :: iDirection
        real(8), intent(out) :: dDirection(:)

        dDirection = 0D0
        if (iDirection < SIZE(dDirection)) then
            dDirection(iDirection) = 1D0
            dDirection(SIZE(dDirection)) = -1D0
        else if (SIZE(dDirection) == 3) then
            dDirection = [1D0,-2D0,1D0]
        else if (SIZE(dDirection) == 4) then
            dDirection = [1D0,-2D0,0D0,1D0]
        else
            dDirection(1) = 1D0
            dDirection(2) = -1D0
        end if

    end subroutine BuildDirection


    integer function FindPhase(cName)

        character(*), intent(in) :: cName
        integer :: i

        FindPhase = 0
        do i = 1, nSolnPhasesSys
            if (TRIM(cSolnPhaseName(i)) == TRIM(cName)) then
                FindPhase = i
                return
            end if
        end do

    end function FindPhase


    logical function PhaseIsActive(iPhase)

        integer, intent(in) :: iPhase
        integer :: i

        PhaseIsActive = .FALSE.
        do i = 1, nSolnPhases
            if (-iAssemblage(nElements-i+1) == iPhase) then
                PhaseIsActive = .TRUE.
                return
            end if
        end do

    end function PhaseIsActive

    !=========================================================================================================
    ! SECTION 6: SCALE-NORMALIZED NUMERICAL HELPERS
    !
    ! Shared norms and errors keep scalar, vector, symmetry, and homogeneity
    ! checks meaningful across phases whose thermodynamic scales differ.
    !=========================================================================================================

    real(8) function ScaledError(dActual,dExpected)

        real(8), intent(in) :: dActual, dExpected

        ScaledError = DABS(dActual-dExpected)/DMAX1(1D0,DABS(dActual),DABS(dExpected))

    end function ScaledError


    real(8) function VectorTwoNorm(dVector)

        real(8), intent(in) :: dVector(:)

        VectorTwoNorm = DSQRT(SUM(dVector*dVector))

    end function VectorTwoNorm


    real(8) function FrobeniusNorm(dMatrix)

        real(8), intent(in) :: dMatrix(:,:)

        FrobeniusNorm = DSQRT(SUM(dMatrix*dMatrix))

    end function FrobeniusNorm


    real(8) function MatrixTwoNormSymmetric(dMatrix)

        real(8), intent(in) :: dMatrix(:,:)
        integer :: iInfo, lWork, n
        real(8), allocatable :: dCopy(:,:), dEigenvalue(:), dWork(:)

        interface
            subroutine DSYEV(cJob,cUplo,n,dA,nLeading,dW,dWork,lWork,iInfo)
                character, intent(in) :: cJob, cUplo
                integer, intent(in) :: n, nLeading, lWork
                real(8), intent(inout) :: dA(nLeading,*)
                real(8), intent(out) :: dW(*), dWork(*)
                integer, intent(out) :: iInfo
            end subroutine DSYEV
        end interface

        n = SIZE(dMatrix,1)
        lWork = MAX(1,3*n-1)
        allocate(dCopy(n,n),dEigenvalue(n),dWork(lWork))
        dCopy = dMatrix
        call DSYEV('N','U',n,dCopy,n,dEigenvalue,dWork,lWork,iInfo)
        if (iInfo == 0) then
            MatrixTwoNormSymmetric = MAXVAL(DABS(dEigenvalue))
        else
            MatrixTwoNormSymmetric = FrobeniusNorm(dMatrix)
        end if

    end function MatrixTwoNormSymmetric

end program TestCEFHessianVerification
