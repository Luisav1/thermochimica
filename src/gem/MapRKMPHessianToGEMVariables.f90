    !---------------------------------------------------------------------------------------------------------
    !
    ! Purpose:
    ! --------
    ! Blend mapped RKMP Hessian contributions into the GEMNewton constrained matrix A.
    !
    ! Current status:
    ! ---------------
    ! This routine is intentionally a scaffold.  It currently applies only the blend operation
    !
    !     A <- A + dAlpha * ARKMP
    !
    ! while ARKMP is still initialized to zero (i.e., mapping logic is not implemented yet).
    !
    ! Variable meaning:
    ! -----------------
    !   A(nVar,nVar)      In/out GEMNewton system matrix (constrained variable space).
    !   ARKMP(nVar,nVar)  Placeholder matrix for mapped RKMP curvature in GEM variable space.
    !   dAlpha            Blend factor used for conservative rollout of added curvature.
    !                     dAlpha is clamped to [0,1] to avoid accidental overshoot.
    !
    ! Notes:
    ! ------
    !   * This staged RKMP rollout intentionally targets RKMP first.
    !   * RKMPM magnetic second-order terms are not part of this scaffold.
    !
    ! TODO:
    !   1) Compute local unconstrained RKMP Hessian for each active RKMP phase.
    !   2) Map local phase/species curvature into GEM variable space.
    !   3) Accumulate mapped contributions into ARKMP before blending below.
    !
    !---------------------------------------------------------------------------------------------------------

subroutine MapRKMPHessianToGEMVariables(A,nVar)

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

    integer                              :: i, j, k, nVar
    integer                              :: iFirst, iLast, nPhaseSpecies, iSolnPhases
    real(8)                              :: dAlpha, dHnn
    real(8), dimension(nVar,nVar)        :: A, ARKMP
    real(8), allocatable, dimension(:,:) :: dHloc ! Local Hessian contribution for a single phase
    real(8), allocatable, dimension(:) :: dX

    ARKMP = 0D0

    if (nSolnPhases <= 0) return

    ! Loop over solution phases currently in the assemblage:
    LOOP_SOLN: do iSolnPhases = 1, nSolnPhases

        ! Absolute solution phase index for this active solution phase in assemblage: (from actual thermo data, not GEM variable indexing)
        k = -iAssemblage(nElements - iSolnPhases + 1)
        if (k <= 0) cycle LOOP_SOLN ! Skip non-solution phases

        ! Stage-1: only RKMP / RKMPM phases
        if (.NOT. (cSolnPhaseType(k) == 'RKMP')) cycle LOOP_SOLN

        ! Find number of species in this phase to determine local Hessian size.
        iFirst = nSpeciesPhase(k-1) + 1
        iLast  = nSpeciesPhase(k)
        nPhaseSpecies = iLast - iFirst + 1
        if (nPhaseSpecies <= 0) cycle LOOP_SOLN

        ! Allocate local Hessian and temporary vector for mapping.
        allocate(dHloc(nPhaseSpecies,nPhaseSpecies), dX(nPhaseSpecies))

        call CompExcessGibbsEnergyRKMP_unconstrained(k,dHloc)
        if (INFOThermo /= 0) then
            deallocate(dHloc, dX)
            return
        end if

        if (lDebugRKMPHessianFD) call DebugRKMPHessianFiniteDifference(k,dHloc)

        ! This is necessary because the unconstrained RKMP Hessian is in terms of species mole numbers, 
        ! while the GEMNewton system is in terms of phase mole numbers.  The mapping from species to phase 
        ! mole numbers is straightforward but must be applied carefully to ensure correct projection of
        ! curvature onto the solution-phase variables in the GEMNewton system.
         ! x_i = n_i / N_phase
        if (dMolesPhase(nElements - iSolnPhases + 1) <= 1D-30) then
            deallocate(dHloc,dX)
            cycle LOOP_SOLN
        end if

        do i = 1, nPhaseSpecies
            dX(i) = dMolesSpecies(iFirst + i - 1) / dMolesPhase(nElements - iSolnPhases + 1)
        end do

        ! Coarse J^T H J projection onto phase-mole variable N_p:
        ! J_s = dn_s/dN_p = x_s
        ! dHnn = x^T * Hloc * x
        dHnn = 0D0
        do i = 1, nPhaseSpecies
            do j = 1, nPhaseSpecies
                dHnn = dHnn + dX(i) * dHloc(i,j) * dX(j)
            end do
        end do

         ! Match GEMNewton indexing: solution-phase unknown lives at j = nElements + iSolnPhases
        j = nElements + iSolnPhases
         ! GEM variable index for this solution phase mole variable.
        if ((j >= 1) .AND. (j <= nVar)) then
            ARKMP(j,j) = ARKMP(j,j) + dHnn
        end if

        deallocate(dHloc, dX)

    end do LOOP_SOLN

    dAlpha = dRKMPHessianBlendAlpha
    dAlpha = DMAX1(0D0, DMIN1(1D0, dAlpha))

    do i = 1, nVar
        A(i,1:nVar) = A(i,1:nVar) + dAlpha * ARKMP(i,1:nVar)
    end do

end subroutine MapRKMPHessianToGEMVariables
