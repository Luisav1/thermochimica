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
!!          native test of the nonuniform pair-specific-zeta distinction.
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
    !> \brief Evaluate production SUBQ scalar blocks at a temporary positive mole state.
    !>
    !> \details The established SUBG production routine also evaluates SUBQ.
    !!          Every global vector it modifies is restored before returning, so
    !!          this local comparison cannot alter the converged calculation.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateProductionEnergy(iPhaseLocal,dMolesLocal,dReferenceIdeal,dExcess,iInfoLocal)

        integer, intent(in) :: iPhaseLocal
        real(8), intent(in) :: dMolesLocal(:)
        real(8), intent(out) :: dReferenceIdeal, dExcess
        integer, intent(out) :: iInfoLocal

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
