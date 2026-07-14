!-------------------------------------------------------------------------------------------------------------
!
!> \file    MapRKMPHessianToGEMVariables.f90
!> \brief   Map plain-RKMP local response curvature into GEMNewton trial systems.
!
!> \details This experimental Stage 1D/1E routine converts the validated RKMP local excess Hessian into a
!! constrained mole-fraction response, subtracts the corresponding ideal-response contribution, and applies the
!! resulting delta to GEMNewton's element-potential block and residual vector.  The mapper is intentionally
!! RKMP-only and is used both for accepted solver updates and for Stage 1E alpha-trust trial solves.
!!
!! Stage 1D supplies the response correction:
!!   deltaA = A_RKMP_response - A_ideal_response
!!   deltaB = N C^T (R_RKMP - R_ideal) mu
!!
!! Stage 1E supplies an explicit trial alpha and controls whether accepted-run metrics are updated.  Trial calls
!! use the same thermodynamic correction as accepted calls, but leave global counters untouched so rejected
!! candidate alphas do not pollute the final audit summary.
!
!> \param[in,out] A GEMNewton matrix.  On return, receives alpha-scaled RKMP response deltas in the element block.
!> \param[in,out] B GEMNewton right-hand side.  On return, receives the matching alpha-scaled RKMP residual delta.
!> \param[in] nVar Number of GEMNewton unknowns represented by A and B.
!> \param[in] dAlphaInput Candidate correction-blend alpha selected by trust logic and internally clamped to
!!                        [0,1].  It scales mapped deltaA and deltaB only after the full RKMP curvature and
!!                        constrained response have been constructed.
!> \param[in] lUpdateMetrics If true, update RKMP audit counters and emit detailed debug diagnostics when enabled.
!> \param[out] lCorrectionOK False when an RKMP mapped correction contains invalid floating-point values.
!> \param[out] dTrialMaxRatio Maximum alpha-scaled A correction relative to the current element block.
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

    do j = 1, nElements
        do i = j, nElements
            A(i,j) = A(i,j) + 0.5D0 * dAlpha * (ARKMP(i,j) + ARKMP(j,i))
            A(j,i) = A(i,j)
        end do
    end do

    B(1:nElements) = B(1:nElements) + dAlpha * BRKMP(1:nElements)

contains

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
