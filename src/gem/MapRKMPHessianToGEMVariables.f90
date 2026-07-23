!-------------------------------------------------------------------------------------------------------------
!
!> \file    MapRKMPHessianToGEMVariables.f90
!> \brief   Map plain-RKMP local response curvature into GEMNewton trial systems.
!
!> \details This routine converts the verified RKMP local excess Hessian into a
!! constrained mole-fraction response, subtracts the corresponding ideal-response contribution, and applies the
!! resulting delta to the GEMNewton element-potential block and residual vector.  The mapper is intentionally
!! RKMP-only and is used both for accepted solver updates and for alpha-trust trial solves.
!!
!! "Condensing" a phase means solving for its local composition response and
!! then eliminating those local composition variables, leaving contributions
!! expressed only in the global GEM unknowns. The response correction has two
!! matching pieces. deltaA changes how
!! element-potential perturbations affect the element-balance equations. deltaB
!! changes the current residual using the same corrected local response:
!!   deltaA = corrected condensed matrix - ideal condensed matrix
!!   deltaB = corrected condensed residual - ideal condensed residual.
!!
!! This stage supplies an explicit trial alpha and controls whether accepted-run metrics are updated.  Trial calls
!! use the same thermodynamic correction as accepted calls, but leave global counters untouched so rejected
!! candidate alphas do not pollute the final audit summary.
!!
!! Conceptual pipeline for each active plain-RKMP phase:
!!   1. Read local species moles, mole fractions, stoichiometry C, and current
!!      chemical-potential forcing from the production GEM state. C(i,e) is the
!!      amount of element e carried by local species i.
!!   2. Obtain Hloc, which describes excess chemical-potential response to
!!      species-mole perturbations. At fixed total phase amount N, multiplying
!!      by N expresses that response per mole-fraction change. Add ideal mixing
!!      curvature, whose diagonal entry for species i is 1/x_i.
!!   3. Solve a bordered local equilibrium system for composition response. Its
!!      entries dx_i are infinitesimal mole-fraction changes, and they must sum
!!      to zero because all mole fractions remain normalized to one.
!!   4. Repeat the same condensation with ideal curvature alone.
!!   5. Apply only the difference between corrected and ideal responses:
!!         deltaA = A_response(corrected) - A_response(ideal)
!!         deltaB = B_response(corrected) - B_response(ideal).
!!
!! The subtraction does not remove RKMP curvature algebraically. It prevents
!! double-counting the ideal response already represented by GEMNewton. This
!! response condensation replaced an earlier diagnostic projection. That
!! projection correctly measured local curvature along selected element-driven
!! composition directions, but it did not solve how phase composition relaxes
!! under those perturbations and was not the reduced matrix required by GEMNewton.
!!
!! Additional notation:
!!   - A is the GEMNewton coefficient matrix and B is its right-hand-side residual.
!!   - Hloc is the local excess Hessian with respect to species moles.
!!   - A "response" is the solved composition change produced by a specified
!!     element-potential or chemical-potential perturbation.
!!   - alpha controls how much of a complete candidate correction is trusted in
!!     this nonlinear iteration; it does not weaken the derivative itself.
!
!> \param[in,out] A              GEMNewton matrix.  On return, receives alpha-scaled RKMP response deltas in the element block.
!> \param[in,out] B              GEMNewton right-hand side.  On return, receives the matching alpha-scaled RKMP residual delta.
!> \param[in]     nVar           Number of GEMNewton unknowns represented by A and B.
!> \param[in]     dAlphaInput.   Candidate correction-blend alpha selected by trust logic and internally clamped to
!!                               [0,1].  It scales mapped deltaA and deltaB only after the full RKMP curvature and
!!                               constrained response have been constructed.
!> \param[in]     lUpdateMetrics If true, update RKMP audit counters and emit detailed debug diagnostics when enabled.
!> \param[out]    lCorrectionOK  False when an RKMP mapped correction contains invalid floating-point values.
!> \param[out]    dTrialMaxRatio Maximum alpha-scaled A correction relative to the current element block.
!
!-------------------------------------------------------------------------------------------------------------

subroutine MapRKMPHessianToGEMVariables(A,B,nVar,dAlphaInput,lUpdateMetrics,lCorrectionOK,dTrialMaxRatio)

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    interface
        subroutine CompExcessGibbsEnergyRKMP_unconstrained(iSolnIndex,dHess)
            integer, intent(in)                  :: iSolnIndex
            real(8), intent(out), dimension(:,:) :: dHess
        end subroutine CompExcessGibbsEnergyRKMP_unconstrained

        subroutine DebugRKMPHessianFiniteDifference(iSolnIndex,dHess)
            integer, intent(in)                  :: iSolnIndex
            real(8), intent(in), dimension(:,:)  :: dHess
        end subroutine DebugRKMPHessianFiniteDifference
    end interface

    integer                              :: i, j, k, p, nVar, nPhaseSpecies, iSolnPhases
    integer                              :: iFirst, iLast, INFO, nDiagFlip, nDiagNegAfter, nDiagNegBefore
    real(8)                              :: dAlphaInput, dAlpha, dMaxA, dMaxDelta, dMaxApplied, dRatio, dTotalMoles
    real(8)                              :: dMaxDeltaB, dMaxAppliedB
    real(8)                              :: dTrialMaxRatio
    logical                              :: lBadDelta
    logical                              :: lUpdateMetrics, lCorrectionOK
    real(8), dimension(nVar,nVar)        :: A, ARKMP
    real(8), dimension(nVar)             :: B, BRKMP
    real(8), allocatable, dimension(:)   :: dX, dLocalMoles, dMu
    real(8), allocatable, dimension(:)   :: dDeltaB
    real(8), allocatable, dimension(:,:) :: dC, dHloc, dHx, dResponse, dIdealCandidate, dIdealResponse
    real(8), allocatable, dimension(:,:) :: dMuRHS, dMuResponse, dIdealMuResponse
    real(8), allocatable, dimension(:,:) :: dCandidate, dDelta

    ARKMP = 0D0
    BRKMP = 0D0
    lCorrectionOK = .TRUE.
    dTrialMaxRatio = 0D0

    if (nSolnPhases <= 0) return

    dAlpha = dAlphaInput
    dAlpha = DMAX1(0D0, DMIN1(1D0, dAlpha))

    !=========================================================================================================
    ! SECTION 1: BUILD ONE CORRECTION FROM EACH ACTIVE PLAIN-RKMP PHASE
    !
    ! Only phases in the current assemblage contribute. RKMPM and other solution
    ! models are skipped because their missing curvature terms would make this
    ! response model thermodynamically incomplete.
    !=========================================================================================================

    LOOP_SOLN: do iSolnPhases = 1, nSolnPhases

        k = -iAssemblage(nElements - iSolnPhases + 1)
        if (k <= 0) cycle LOOP_SOLN
        if (cSolnPhaseType(k) /= 'RKMP') cycle LOOP_SOLN

        iFirst = nSpeciesPhase(k-1) + 1
        iLast  = nSpeciesPhase(k)
        nPhaseSpecies = iLast - iFirst + 1
        if (nPhaseSpecies <= 1) cycle LOOP_SOLN

        allocate(dX(nPhaseSpecies), dLocalMoles(nPhaseSpecies), dMu(nPhaseSpecies), &
                 dC(nPhaseSpecies,nElements), &
                 dHloc(nPhaseSpecies,nPhaseSpecies), dHx(nPhaseSpecies,nPhaseSpecies), &
                 dResponse(nPhaseSpecies,nElements), dIdealResponse(nPhaseSpecies,nElements), &
                 dMuRHS(nPhaseSpecies,1), dMuResponse(nPhaseSpecies,1), &
                 dIdealMuResponse(nPhaseSpecies,1), dIdealCandidate(nElements,nElements), &
                 dCandidate(nElements,nElements), dDelta(nElements,nElements), dDeltaB(nElements))

        dLocalMoles = dMolesSpecies(iFirst:iLast)
        dTotalMoles = SUM(dLocalMoles)
        if (dTotalMoles <= 1D-30) then
            deallocate(dX, dLocalMoles, dMu, dC, dHloc, dHx, dResponse, dIdealResponse, &
                       dMuRHS, dMuResponse, dIdealMuResponse, dIdealCandidate, dCandidate, &
                       dDelta, dDeltaB)
            cycle LOOP_SOLN
        end if

        dX = DMAX1(dLocalMoles / dTotalMoles, 1D-30)
        dX = dX / SUM(dX)
        dMu = dChemicalPotential(iFirst:iLast)
        dMuRHS(:,1) = dMu

        do j = 1, nElements
            do i = 1, nPhaseSpecies
                p = iFirst + i - 1
                dC(i,j) = dStoichSpecies(p,j) / DFLOAT(iParticlesPerMole(p))
            end do
        end do

        !---------------------------------------------------------------------------------------------
        ! SECTION 1A: CORRECTED CONSTRAINED LOCAL RESPONSE
        !
        ! Hloc describes excess chemical-potential response to species-mole
        ! changes. Holding total phase amount fixed converts a mole-fraction
        ! change dx into a mole change N*dx. The diagonal 1/x_i terms add the
        ! corresponding ideal-mixing response for each species.
        !---------------------------------------------------------------------------------------------
        call CompExcessGibbsEnergyRKMP_unconstrained(k,dHloc)
        if (INFOThermo /= 0) then
            deallocate(dX, dLocalMoles, dMu, dC, dHloc, dHx, dResponse, dIdealResponse, &
                       dMuRHS, dMuResponse, dIdealMuResponse, dIdealCandidate, dCandidate, &
                       dDelta, dDeltaB)
            return
        end if

        if (lDebugRKMPHessianFD .AND. lUpdateMetrics) call DebugRKMPHessianFiniteDifference(k,dHloc)

        ! Convert mole-number excess curvature to fixed-phase-amount mole-fraction curvature, then add
        ! ideal mixing curvature to form the local response matrix.
        dHx = dTotalMoles * dHloc
        do i = 1, nPhaseSpecies
            dHx(i,i) = dHx(i,i) + 1D0 / dX(i)
        end do

        call SolveLocalResponse(nPhaseSpecies, nElements, dHx, dC, dResponse, INFO)
        if (INFO /= 0) then
            if (lDebugRKMPHessianFD) then
                write(*,'(A,1X,A,1X,I0,1X,A,1X,I0)') &
                    'RKMP_STAGE1D_MAP_DEBUG', TRIM(cSolnPhaseName(k)), k, 'response_solve_info=', INFO
            end if
            deallocate(dX, dLocalMoles, dMu, dC, dHloc, dHx, dResponse, dIdealResponse, &
                       dMuRHS, dMuResponse, dIdealMuResponse, dIdealCandidate, dCandidate, &
                       dDelta, dDeltaB)
            cycle LOOP_SOLN
        end if

        call SolveLocalResponse(nPhaseSpecies, 1, dHx, dMuRHS, dMuResponse, INFO)
        if (INFO /= 0) then
            deallocate(dX, dLocalMoles, dMu, dC, dHloc, dHx, dResponse, dIdealResponse, &
                       dMuRHS, dMuResponse, dIdealMuResponse, dIdealCandidate, dCandidate, &
                       dDelta, dDeltaB)
            cycle LOOP_SOLN
        end if

        !---------------------------------------------------------------------------------------------
        ! SECTION 1B: IDEAL BASELINE AND RESPONSE DELTA
        !
        ! GEMNewton already contains an ideal/simple phase response. Reconstruct
        ! that baseline with the identical constraint and forcing, then retain
        ! only the change caused by adding RKMP excess curvature.
        !---------------------------------------------------------------------------------------------
        dCandidate = MATMUL(TRANSPOSE(dC), dTotalMoles * dResponse)
        call ComputeIdealCandidate(nPhaseSpecies, nElements, dX, dC, dMuRHS, dTotalMoles, &
                                   dIdealCandidate, dIdealResponse, dIdealMuResponse, INFO)
        if (INFO /= 0) then
            if (lDebugRKMPHessianFD) then
                write(*,'(A,1X,A,1X,I0,1X,A,1X,I0)') &
                    'RKMP_STAGE1D_MAP_DEBUG', TRIM(cSolnPhaseName(k)), k, 'ideal_solve_info=', INFO
            end if
            deallocate(dX, dLocalMoles, dMu, dC, dHloc, dHx, dResponse, dIdealResponse, &
                       dMuRHS, dMuResponse, dIdealMuResponse, dIdealCandidate, dCandidate, &
                       dDelta, dDeltaB)
            cycle LOOP_SOLN
        end if

        dDelta = dCandidate - dIdealCandidate
        dDeltaB = MATMUL(TRANSPOSE(dC), dTotalMoles * &
            (dMuResponse(:,1) - dIdealMuResponse(:,1)))

        !---------------------------------------------------------------------------------------------
        ! SECTION 1C: FINITE-VALUE GUARDS AND AUDIT METRICS
        !
        ! Trial corrections are measured but not applied here. SolveRKMPAlphaTrust
        ! uses these metrics together with the solved update to choose alpha.
        !---------------------------------------------------------------------------------------------
        lBadDelta = .FALSE.
        do j = 1, nElements
            do i = 1, nElements
                if ((dDelta(i,j) /= dDelta(i,j)) .OR. (DABS(dDelta(i,j)) > 0.5D0 * HUGE(1D0))) then
                    lBadDelta = .TRUE.
                end if
            end do
        end do

        if (.NOT. lBadDelta) then
            dMaxA       = MAXVAL(DABS(A(1:nElements,1:nElements)))
            dMaxDelta   = MAXVAL(DABS(dDelta))
            dMaxDeltaB  = MAXVAL(DABS(dDeltaB))
            dMaxApplied = dAlpha * dMaxDelta
            dMaxAppliedB = dAlpha * dMaxDeltaB
            dRatio      = dMaxApplied / DMAX1(dMaxA, 1D-30)
            dTrialMaxRatio = DMAX1(dTrialMaxRatio, dRatio)
            if (lUpdateMetrics) then
                nRKMPHessianApplyCount = nRKMPHessianApplyCount + 1
                dRKMPHessianMaxAppliedA = DMAX1(dRKMPHessianMaxAppliedA, dMaxApplied)
                dRKMPHessianMaxAppliedB = DMAX1(dRKMPHessianMaxAppliedB, dMaxAppliedB)
                dRKMPHessianMaxAppliedRatio = DMAX1(dRKMPHessianMaxAppliedRatio, dRatio)
                dRKMPHessianMaxDeltaA = DMAX1(dRKMPHessianMaxDeltaA, dMaxDelta)
                dRKMPHessianMaxDeltaB = DMAX1(dRKMPHessianMaxDeltaB, dMaxDeltaB)
            end if
            nDiagFlip   = 0
            nDiagNegAfter  = 0
            nDiagNegBefore = 0
            do i = 1, nElements
                if (A(i,i) < 0D0) nDiagNegBefore = nDiagNegBefore + 1
                if ((A(i,i) + dAlpha*dDelta(i,i)) < 0D0) nDiagNegAfter = nDiagNegAfter + 1
                if ((A(i,i) * (A(i,i) + dAlpha*dDelta(i,i))) < 0D0) nDiagFlip = nDiagFlip + 1
            end do

            if (lDebugRKMPHessianFD .AND. lUpdateMetrics) then
                write(*,'(A,1X,A,1X,I0,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,ES14.6,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)') &
                    'RKMP_STAGE1D_MAP_DEBUG', TRIM(cSolnPhaseName(k)), k, &
                    'maxA=', dMaxA, 'maxDelta=', dMaxDelta, 'maxApplied=', dMaxApplied, &
                    'ratio=', dRatio, 'maxDeltaB=', dMaxDeltaB, 'maxAppliedB=', dMaxAppliedB, &
                    'diagNegBefore=', nDiagNegBefore, &
                    'diagNegAfter=', nDiagNegAfter, 'diagFlip=', nDiagFlip
            end if

            ARKMP(1:nElements,1:nElements) = ARKMP(1:nElements,1:nElements) + dDelta
            BRKMP(1:nElements) = BRKMP(1:nElements) + dDeltaB
        else
            lCorrectionOK = .FALSE.
            if (lDebugRKMPHessianFD .AND. lUpdateMetrics) then
                write(*,'(A,1X,A,1X,I0,1X,A)') &
                    'RKMP_STAGE1D_MAP_DEBUG', TRIM(cSolnPhaseName(k)), k, 'skipped_bad_delta'
            end if
        end if

        deallocate(dX, dLocalMoles, dMu, dC, dHloc, dHx, dResponse, dIdealResponse, &
                   dMuRHS, dMuResponse, dIdealMuResponse, dIdealCandidate, dCandidate, &
                   dDelta, dDeltaB)

    end do LOOP_SOLN

    !=========================================================================================================
    ! SECTION 2: APPLY THE AGGREGATED TRIAL CORRECTION
    !
    ! Alpha controls how much of the complete RKMP correction is used in this
    ! trial Newton system. This is part of globalization: the safeguards that
    ! keep a locally derived Newton step useful while the current state may
    ! still be far from equilibrium. Alpha does not rescale the underlying
    ! thermodynamic derivative.
    ! The element block is symmetrized as it is inserted into GEMNewton.
    !=========================================================================================================
    do j = 1, nElements
        do i = j, nElements
            A(i,j) = A(i,j) + 0.5D0 * dAlpha * (ARKMP(i,j) + ARKMP(j,i))
            A(j,i) = A(i,j)
        end do
    end do

    B(1:nElements) = B(1:nElements) + dAlpha * BRKMP(1:nElements)

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Solve the normalized local composition response to one or more forcing directions.
    !>
    !> \details The bordered matrix combines local chemical-potential curvature
    !!          with one normalization multiplier. Its final row requires the
    !!          infinitesimal mole-fraction changes to sum to zero, preserving
    !!          the definition that all mole fractions sum to one. The solve
    !!          changes composition within one already active phase; it neither
    !!          fixes nor constrains which phases belong to the global assemblage.
    !---------------------------------------------------------------------------------------------------------
    subroutine SolveLocalResponse(nSpecies, nElem, dHxLocal, dCLocal, dResp, INFO)

        integer, intent(in)                    :: nSpecies, nElem
        integer, intent(out)                   :: INFO
        real(8), intent(in), dimension(:,:)    :: dHxLocal, dCLocal
        real(8), intent(out), dimension(:,:)   :: dResp

        integer                                :: ii, jj, nEqn
        integer, dimension(:), allocatable     :: IPIV
        real(8), dimension(:,:), allocatable   :: dKKT, dRHS

        nEqn = nSpecies + 1
        allocate(dKKT(nEqn,nEqn), dRHS(nEqn,nElem), IPIV(nEqn))
        dKKT = 0D0
        dRHS = 0D0

        dKKT(1:nSpecies,1:nSpecies) = dHxLocal(1:nSpecies,1:nSpecies)
        do ii = 1, nSpecies
            dKKT(ii,nEqn) = -1D0
            dKKT(nEqn,ii) = 1D0
        end do
        dRHS(1:nSpecies,1:nElem) = dCLocal(1:nSpecies,1:nElem)

        call DGESV(nEqn, nElem, dKKT, nEqn, IPIV, dRHS, nEqn, INFO)
        if (INFO == 0) then
            do jj = 1, nElem
                dResp(1:nSpecies,jj) = dRHS(1:nSpecies,jj)
            end do
        else
            dResp = 0D0
        end if

        deallocate(dKKT, dRHS, IPIV)

    end subroutine SolveLocalResponse

    !---------------------------------------------------------------------------------------------------------
    !> \brief Reconstruct the ideal-only local response already represented by GEMNewton.
    !>
    !> \details Using the same normalization constraint and forcing as the
    !!          corrected solve isolates the excess-curvature effect as a
    !!          response difference rather than as a raw added stiffness.
    !---------------------------------------------------------------------------------------------------------
    subroutine ComputeIdealCandidate(nSpecies, nElem, dXLocal, dCLocal, dMuLocal, dPhaseMoles, &
                                     dIdeal, dRespIdeal, dMuRespIdeal, INFO)

        integer, intent(in)                    :: nSpecies, nElem
        integer, intent(out)                   :: INFO
        real(8), intent(in), dimension(:)      :: dXLocal
        real(8), intent(in), dimension(:,:)    :: dCLocal
        real(8), intent(in), dimension(:,:)    :: dMuLocal
        real(8), intent(in)                    :: dPhaseMoles
        real(8), intent(out), dimension(:,:)   :: dIdeal
        real(8), intent(out), dimension(:,:)   :: dRespIdeal, dMuRespIdeal

        integer                                :: ii
        real(8), dimension(:,:), allocatable   :: dHideal

        allocate(dHideal(nSpecies,nSpecies))
        dHideal = 0D0
        do ii = 1, nSpecies
            dHideal(ii,ii) = 1D0 / dXLocal(ii)
        end do

        call SolveLocalResponse(nSpecies, nElem, dHideal, dCLocal, dRespIdeal, INFO)
        if (INFO == 0) then
            dIdeal = MATMUL(TRANSPOSE(dCLocal), dPhaseMoles * dRespIdeal)
            call SolveLocalResponse(nSpecies, 1, dHideal, dMuLocal, dMuRespIdeal, INFO)
        else
            dIdeal = 0D0
            dMuRespIdeal = 0D0
        end if

        deallocate(dHideal)

    end subroutine ComputeIdealCandidate

end subroutine MapRKMPHessianToGEMVariables
