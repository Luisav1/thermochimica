!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleMQMQAProductionAdapter.f90
!> \brief   Decode production plain-SUBG or SUBQ phases into the generic MQMQA interface.
!>
!> \details Thermochimica stores MQMQA topology, coefficients, and constituent
!!          grouping rules in filtered runtime arrays owned by ModuleThermo.
!!          This module translates those arrays into MQMQAModelData and
!!          MQMQAInteractionTerm objects without evaluating energy or changing
!!          the production state.
!!
!!          The adapter is deliberately narrower than the generic mathematics:
!!          - plain nonmagnetic SUBG and SUBQ local data only;
!!          - production G, Q, and B parameter labels;
!!          - the traced first-sublattice ternary orientation;
!!          - no reciprocal R, magnetic, or solver coupling behavior.
!!
!!          Keeping this translation separate lets native tests prove that the
!!          disconnected Hessian consumes the same topology and parameters as
!!          production Thermochimica. A later constrained-response mapper can
!!          reuse the adapter after native verification passes.
!-------------------------------------------------------------------------------------------------------------

module ModuleMQMQAProductionAdapter

    USE ModuleThermo, ONLY: cRegularParam, cSolnPhaseType, dCoordinationNumber, dExcessGibbsParam, &
        dStdGibbsEnergy, dZetaSpecies, iChemicalGroup, iConstituentSublattice, iInterpolationOverride, &
        iPairID, iPhaseSublattice, iRegularParam, nConstituentSublattice, nInterpolationOverride, &
        nPairsSRO, nParamPhase, nSpeciesPhase
    USE ModuleMQMQAUnconstrained, ONLY: MQMQAModelData, MQMQAInteractionTerm, &
        MQMQA_MODEL_SUBG, MQMQA_MODEL_SUBQ, MQMQA_TERM_G, MQMQA_TERM_Q, MQMQA_TERM_B

    implicit none
    private

    public :: DecodeProductionSUBGPhase, DecodeProductionSUBQPhase

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Translate one filtered production plain-SUBG phase.
    !>
    !> \param[in]  iPhase       Runtime solution-phase index.
    !> \param[out] tModel       Generic quadruplet topology and constant data.
    !> \param[out] tInteraction Active, supported excess-interaction records.
    !> \param[out] iInfo        Zero on success; nonzero when production data fall outside the traced scope.
    !>
    !> \details Parameters whose temperature-evaluated coefficient is exactly
    !!          zero are omitted, matching the early cycle in
    !!          CompExcessGibbsEnergySUBG. Second-sublattice constituent indices
    !!          are converted from production's combined index space to the local
    !!          indexing expected by ModuleMQMQAUnconstrained.
    !---------------------------------------------------------------------------------------------------------
    subroutine DecodeProductionSUBGPhase(iPhase,tModel,tInteraction,iInfo)

        integer, intent(in) :: iPhase
        type(MQMQAModelData), intent(out) :: tModel
        type(MQMQAInteractionTerm), allocatable, intent(out) :: tInteraction(:)
        integer, intent(out) :: iInfo

        call DecodeProductionPhase(iPhase,'SUBG',MQMQA_MODEL_SUBG,tModel,tInteraction,iInfo)

    end subroutine DecodeProductionSUBGPhase


    !---------------------------------------------------------------------------------------------------------
    !> \brief Translate one filtered production plain-SUBQ phase.
    !>
    !> \param[in]  iPhase       Runtime solution-phase index.
    !> \param[out] tModel       Generic quadruplet topology and constant data.
    !> \param[out] tInteraction Active, supported excess-interaction records.
    !> \param[out] iInfo        Zero on success; nonzero when production data fall outside the traced scope.
    !>
    !> \details SUBQ uses the same production topology and parameter records as
    !!          SUBG. The model selector activates the SUBQ configurational
    !!          exponents and environment weights and permits the decoded database
    !!          to supply pair-specific zeta values. Keeping this entry point
    !!          distinct prevents existing SUBG-only solver paths from becoming
    !!          SUBQ-enabled accidentally.
    !---------------------------------------------------------------------------------------------------------
    subroutine DecodeProductionSUBQPhase(iPhase,tModel,tInteraction,iInfo)

        integer, intent(in) :: iPhase
        type(MQMQAModelData), intent(out) :: tModel
        type(MQMQAInteractionTerm), allocatable, intent(out) :: tInteraction(:)
        integer, intent(out) :: iInfo

        call DecodeProductionPhase(iPhase,'SUBQ',MQMQA_MODEL_SUBQ,tModel,tInteraction,iInfo)

    end subroutine DecodeProductionSUBQPhase


    !---------------------------------------------------------------------------------------------------------
    !> \brief Shared production-array translation used by the strict SUBG and SUBQ entry points.
    !>
    !> \details The caller supplies both the required production phase label and
    !!          the corresponding generic formulation selector. All remaining
    !!          topology, zeta, reference-energy, and interaction decoding is
    !!          intentionally identical so the two public paths cannot drift.
    !---------------------------------------------------------------------------------------------------------
    subroutine DecodeProductionPhase(iPhase,cRequiredType,iModelType,tModel,tInteraction,iInfo)

        integer, intent(in) :: iPhase, iModelType
        character(len=*), intent(in) :: cRequiredType
        type(MQMQAModelData), intent(out) :: tModel
        type(MQMQAInteractionTerm), allocatable, intent(out) :: tInteraction(:)
        integer, intent(out) :: iInfo

        integer :: a, iFirst, iInteraction, iPair, iParam, iSPI
        integer :: nActive, nQuad, nSub1, nSub2, x
        integer, allocatable :: iZetaCount(:,:)

        iInfo = 0
        if ((iPhase < 1) .OR. (iPhase > SIZE(cSolnPhaseType))) then
            iInfo = 1
            return
        end if
        if (cSolnPhaseType(iPhase) /= cRequiredType) then
            iInfo = 2
            return
        end if

        iSPI = iPhaseSublattice(iPhase)
        nSub1 = nConstituentSublattice(iSPI,1)
        nSub2 = nConstituentSublattice(iSPI,2)
        nQuad = nPairsSRO(iSPI,2)
        iFirst = nSpeciesPhase(iPhase-1)+1
        if ((nSub1 <= 0) .OR. (nSub2 <= 0) .OR. (nQuad <= 0) .OR. &
            (nSpeciesPhase(iPhase)-iFirst+1 /= nQuad)) then
            iInfo = 3
            return
        end if

        tModel%iModelType = iModelType
        tModel%nSublattice1 = nSub1
        tModel%nSublattice2 = nSub2
        allocate(tModel%iQuadruplet(nQuad,4),tModel%dCoordination(nQuad,4), &
            tModel%dZeta(nSub1,nSub2),tModel%dReferenceEnergy(nQuad),iZetaCount(nSub1,nSub2))

        tModel%iQuadruplet(:,1:2) = iPairID(iSPI,1:nQuad,1:2)
        tModel%iQuadruplet(:,3:4) = iPairID(iSPI,1:nQuad,3:4)-nSub1
        tModel%dCoordination = dCoordinationNumber(iSPI,1:nQuad,1:4)
        tModel%dReferenceEnergy = dStdGibbsEnergy(iFirst:nSpeciesPhase(iPhase))

        ! Production stores zeta on retained A-X pair records. Reconstructing a
        ! matrix makes every pair lookup explicit and detects missing or duplicate
        ! records instead of silently relying on parser ordering.
        tModel%dZeta = 0D0
        iZetaCount = 0
        do iPair = 1, nPairsSRO(iSPI,1)
            a = iConstituentSublattice(iSPI,1,iPair)
            x = iConstituentSublattice(iSPI,2,iPair)
            if ((a < 1) .OR. (a > nSub1) .OR. (x < 1) .OR. (x > nSub2)) then
                iInfo = 4
                return
            end if
            tModel%dZeta(a,x) = dZetaSpecies(iSPI,iPair)
            iZetaCount(a,x) = iZetaCount(a,x)+1
        end do
        if (ANY(iZetaCount /= 1) .OR. ANY(tModel%dZeta <= 0D0)) then
            iInfo = 5
            return
        end if

        nActive = 0
        do iParam = nParamPhase(iPhase-1)+1, nParamPhase(iPhase)
            if (dExcessGibbsParam(iParam) /= 0D0) nActive = nActive+1
        end do
        allocate(tInteraction(nActive))

        iInteraction = 0
        do iParam = nParamPhase(iPhase-1)+1, nParamPhase(iPhase)
            if (dExcessGibbsParam(iParam) == 0D0) cycle
            iInteraction = iInteraction+1
            call DecodeInteraction(iPhase,iSPI,nSub1,nSub2,iParam,tInteraction(iInteraction),iInfo)
            if (iInfo /= 0) return
        end do

        ! Canonical ordering is part of the generic module's topology contract.
        if (ANY(tModel%iQuadruplet(:,1) > tModel%iQuadruplet(:,2)) .OR. &
            ANY(tModel%iQuadruplet(:,3) > tModel%iQuadruplet(:,4)) .OR. &
            ANY(tModel%dCoordination <= 0D0)) then
            iInfo = 6
            return
        end if

    end subroutine DecodeProductionPhase


    !---------------------------------------------------------------------------------------------------------
    !> \brief Decode one active G, Q, or B runtime parameter.
    !>
    !> \details G and Q terms carry asymmetric-group masks. These masks reproduce
    !!          the production distinction between constituents that interpolate
    !!          with endpoint A, endpoint B, or neither. B terms use weighted pair
    !!          fractions directly and therefore require no group masks.
    !---------------------------------------------------------------------------------------------------------
    subroutine DecodeInteraction(iPhase,iSPI,nSub1,nSub2,iParam,tTerm,iInfo)

        integer, intent(in) :: iPhase, iSPI, nSub1, nSub2, iParam
        type(MQMQAInteractionTerm), intent(out) :: tTerm
        integer, intent(out) :: iInfo

        iInfo = 0
        select case (cRegularParam(iParam))
        case ('G')
            tTerm%iFamily = MQMQA_TERM_G
        case ('Q')
            tTerm%iFamily = MQMQA_TERM_Q
        case ('B')
            tTerm%iFamily = MQMQA_TERM_B
        case default
            ! Reciprocal R is parser-recognized but currently has no implemented
            ! scalar contribution. Unknown labels are likewise never guessed.
            iInfo = 10
            return
        end select

        if ((iRegularParam(iParam,1) /= 3) .AND. (iRegularParam(iParam,1) /= 4)) then
            iInfo = 11
            return
        end if
        tTerm%iA = iRegularParam(iParam,2)
        tTerm%iB = iRegularParam(iParam,3)
        tTerm%iX = iRegularParam(iParam,4)-nSub1
        tTerm%iY = iRegularParam(iParam,5)-nSub1
        tTerm%iExponentP = iRegularParam(iParam,6)
        tTerm%iExponentQ = iRegularParam(iParam,7)
        tTerm%iExponentR = iRegularParam(iParam,8)
        tTerm%iTernaryConstituent = iRegularParam(iParam,10)
        tTerm%dCoefficient = dExcessGibbsParam(iParam)

        if ((iRegularParam(iParam,9) /= 0) .OR. (iRegularParam(iParam,11) /= 0) .OR. &
            (tTerm%iExponentP < 0) .OR. (tTerm%iExponentQ < 0) .OR. (tTerm%iExponentR < 0)) then
            iInfo = 12
            return
        end if
        if ((tTerm%iA < 1) .OR. (tTerm%iA > nSub1) .OR. &
            (tTerm%iB < 1) .OR. (tTerm%iB > nSub1) .OR. &
            (tTerm%iX < 1) .OR. (tTerm%iX > nSub2) .OR. &
            (tTerm%iY < 1) .OR. (tTerm%iY > nSub2)) then
            iInfo = 13
            return
        end if

        if (tTerm%iFamily == MQMQA_TERM_B) then
            if (tTerm%iTernaryConstituent /= 0) iInfo = 14
            return
        end if
        if ((tTerm%iA == tTerm%iB) .EQV. (tTerm%iX == tTerm%iY)) then
            iInfo = 15
            return
        end if
        if (tTerm%iTernaryConstituent > 0) then
            if ((tTerm%iX /= tTerm%iY) .OR. (tTerm%iTernaryConstituent > nSub1) .OR. &
                (tTerm%iExponentR < 1)) then
                iInfo = 16
                return
            end if
        end if

        call BuildAsymmetricGroups(iPhase,iSPI,nSub1,nSub2,tTerm,iInfo)

    end subroutine DecodeInteraction


    !---------------------------------------------------------------------------------------------------------
    !> \brief Reproduce production asymmetric interpolation groups for one G/Q term.
    !>
    !> \details For a fixed X-X environment, A and B define groups on the first
    !!          sublattice. Database interpolation overrides take precedence over
    !!          chemical-group defaults. For a fixed A-A environment, X and Y
    !!          define groups on the second sublattice; production uses chemical
    !!          groups directly in that orientation.
    !---------------------------------------------------------------------------------------------------------
    subroutine BuildAsymmetricGroups(iPhase,iSPI,nSub1,nSub2,tTerm,iInfo)

        integer, intent(in) :: iPhase, iSPI, nSub1, nSub2
        type(MQMQAInteractionTerm), intent(inout) :: tTerm
        integer, intent(out) :: iInfo

        integer :: i, k, l
        logical :: lIsException, lMatches

        iInfo = 0
        if (tTerm%iX == tTerm%iY) then
            allocate(tTerm%lGroup1(nSub1),tTerm%lGroup2(nSub1))
            tTerm%lGroup1 = .FALSE.
            tTerm%lGroup2 = .FALSE.
            tTerm%lGroup1(tTerm%iA) = .TRUE.
            tTerm%lGroup2(tTerm%iB) = .TRUE.

            do i = 1, nSub1
                lIsException = .FALSE.
                do k = 1, nInterpolationOverride(iPhase)
                    lMatches = .TRUE.
                    do l = 1, 3
                        if (.NOT. ((iInterpolationOverride(iPhase,k,l) == tTerm%iA) .OR. &
                            (iInterpolationOverride(iPhase,k,l) == tTerm%iB) .OR. &
                            (iInterpolationOverride(iPhase,k,l) == i))) lMatches = .FALSE.
                    end do
                    if (.NOT. lMatches) cycle
                    lIsException = .TRUE.
                    if (iInterpolationOverride(iPhase,k,5) == tTerm%iB) tTerm%lGroup1(i) = .TRUE.
                    if (iInterpolationOverride(iPhase,k,5) == tTerm%iA) tTerm%lGroup2(i) = .TRUE.
                    exit
                end do
                if (lIsException) cycle
                if (iChemicalGroup(iSPI,1,tTerm%iA) /= iChemicalGroup(iSPI,1,tTerm%iB)) then
                    if (iChemicalGroup(iSPI,1,i) == iChemicalGroup(iSPI,1,tTerm%iA)) then
                        tTerm%lGroup1(i) = .TRUE.
                    else if (iChemicalGroup(iSPI,1,i) == iChemicalGroup(iSPI,1,tTerm%iB)) then
                        tTerm%lGroup2(i) = .TRUE.
                    end if
                end if
            end do
        else
            allocate(tTerm%lGroup1(nSub2),tTerm%lGroup2(nSub2))
            tTerm%lGroup1 = .FALSE.
            tTerm%lGroup2 = .FALSE.
            tTerm%lGroup1(tTerm%iX) = .TRUE.
            tTerm%lGroup2(tTerm%iY) = .TRUE.
            if (iChemicalGroup(iSPI,2,tTerm%iX) /= iChemicalGroup(iSPI,2,tTerm%iY)) then
                do i = 1, nSub2
                    if (iChemicalGroup(iSPI,2,i) == iChemicalGroup(iSPI,2,tTerm%iX)) then
                        tTerm%lGroup1(i) = .TRUE.
                    else if (iChemicalGroup(iSPI,2,i) == iChemicalGroup(iSPI,2,tTerm%iY)) then
                        tTerm%lGroup2(i) = .TRUE.
                    end if
                end do
            end if
        end if

        if (ANY(tTerm%lGroup1 .AND. tTerm%lGroup2) .OR. &
            (.NOT. ANY(tTerm%lGroup1)) .OR. (.NOT. ANY(tTerm%lGroup2))) iInfo = 20

    end subroutine BuildAsymmetricGroups

end module ModuleMQMQAProductionAdapter
