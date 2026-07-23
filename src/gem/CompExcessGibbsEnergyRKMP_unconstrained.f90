!-------------------------------------------------------------------------------------------------------------
!> \file    CompExcessGibbsEnergyRKMP_unconstrained.f90
!> \brief   Compute the local mole-space Hessian of supported plain-RKMP excess energy.
!>
!> \details This routine differentiates the extensive binary RKMP excess energy
!!          with respect to the species mole amounts of one production phase:
!!
!!              Hloc(i,j) = d2 Gex / (d n_i d n_j).
!!
!!          It reads the converged or current phase state from ModuleThermo. The
!!          result is excess-only local curvature: it does not include ideal
!!          mixing, preserve mole-fraction normalization during a composition
!!          response solve, map to element potentials, or modify
!!          GEMNewton. Those responsibilities belong to
!!          MapRKMPHessianToGEMVariables.
!!
!!          Conceptual pipeline:
!!          1. Select one active plain-RKMP phase and its species moles.
!!          2. For each supported binary parameter, write its extensive energy as
!!                G_lambda = Lambda * A(n) * L(Delta(n)),
!!             where Lambda is the database interaction coefficient, A is the
!!             amount factor n_a*n_b/N, Delta is the mole-fraction difference
!!             (n_a-n_b)/N, and L raises that difference to the parameter order.
!!          3. Apply the product and chain rules to A and Delta for every pair of
!!             local species-mole directions.
!!          4. Sum parameter contributions into the local Hessian.
!!
!!          Notation used below:
!!          - i and j label the local species whose mole amounts are being
!!            differentiated.
!!          - a and b label the two species named by one binary RKMP parameter.
!!          - n_i is the mole amount of local species i and N is total phase moles.
!!          - x_i=n_i/N is the species mole fraction.
!!          - Gex is the extensive excess Gibbs energy of this phase.
!!          - Hloc(i,j) measures how the excess chemical potential of species i
!!            changes when the mole amount of species j is perturbed.
!!          - A and Delta are temporary derivation symbols, not the GEMNewton
!!            coefficient matrix A or a solver update.
!!
!!          Scope is intentionally limited to nonmagnetic binary plain-RKMP
!!          parameters. Ternary, higher-component Muggiano, and RKMPM magnetic
!!          curvature are not implemented here.
!>
!> \param[in]  iSolnIndex Absolute production solution-phase index.
!> \param[out] dHess      Excess Hessian in local phase-species mole coordinates.
!-------------------------------------------------------------------------------------------------------------
    
subroutine CompExcessGibbsEnergyRKMP_unconstrained(iSolnIndex,dHess)

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer, intent(in)                  :: iSolnIndex
    real(8), intent(out), dimension(:,:) :: dHess

    integer :: i, j, iFirstSpecies, iLastSpecies, nPhaseSpecies
    integer :: iParam, iExponent, ia, ib
    real(8) :: dN, dNi, dNj, dB, dC, dA, dDelta
    real(8) :: dL0, dL1, dL2, dHij, dSym, dEpsN
    real(8) :: dDelta_ai, dDelta_aj, dDelta_bi, dDelta_bj
    real(8) :: dBa, dBb, dBab, dCa, dCb, dCab
    real(8) :: dAa, dAb, dAab, dDa, dDb, dDab

    !=========================================================================================================
    ! SECTION 1: PHASE DOMAIN AND LOCAL STATE
    !
    ! This is a production-linked local routine: phase topology, parameters, and
    ! species moles are obtained through the supplied ModuleThermo phase index.
    !=========================================================================================================

    ! Return if phase type does not match RKMP model:
    if (.NOT. (cSolnPhaseType(iSolnIndex) == 'RKMP')) then
        dHess = 0D0
        return
    end if

    iFirstSpecies = nSpeciesPhase(iSolnIndex-1) + 1
    iLastSpecies  = nSpeciesPhase(iSolnIndex)
    nPhaseSpecies = iLastSpecies - iFirstSpecies + 1

    ! Caller provides storage for dHess; fail early if dimensions are too small.
    if ((SIZE(dHess,1) < nPhaseSpecies) .OR. (SIZE(dHess,2) < nPhaseSpecies)) then
        ! Dedicated error code for Hessian buffer shape mismatch in RKMP unconstrained routine.
        INFOThermo = 46
        return
    end if

    ! Always rebuild the local Hessian from scratch for this phase.
    dHess(1:nPhaseSpecies,1:nPhaseSpecies) = 0D0

    ! dEpsN protects N-based denominators against near-zero phase totals.
    dEpsN = 1D-30
    dN = 0D0

    ! Sum up total moles in this phase by looping over local species indices. This is needed for all RKMP terms, so do it once here.
    do i = 1, nPhaseSpecies
        dN = dN + dMolesSpecies(iFirstSpecies + i - 1)
    end do
    ! dN appears in multiple denominators below. Clamp to avoid inf/NaN for tiny phases.
    dN = DMAX1(dN,dEpsN)

    !=========================================================================================================
    ! SECTION 2: BINARY RKMP PARAMETER CURVATURE
    !
    ! Each parameter contributes Lambda*A*L(Delta). Compact intermediates keep
    ! the implementation close to the analytic derivation while still allowing
    ! derivatives with respect to every local species, including species that
    ! enter only through the total phase amount N.
    !=========================================================================================================

    LOOP_Param: do iParam = nParamPhase(iSolnIndex-1)+1, nParamPhase(iSolnIndex)

        ! Implementation supports binary RKMP terms only.
        ! Ternary and higher-order cases are intentionally skipped.
        if (iRegularParam(iParam,1) /= 2) cycle LOOP_Param

        ! Local species IDs participating in this binary interaction parameter.
        ia = iRegularParam(iParam,2)
        ib = iRegularParam(iParam,3)

        if ((ia < 1) .OR. (ia > nPhaseSpecies) .OR. (ib < 1) .OR. (ib > nPhaseSpecies)) cycle LOOP_Param

        ! Polynomial exponent for the Redlich-Kister term of this parameter.
        iExponent = iRegularParam(iParam,4)
        if (iExponent < 0) cycle LOOP_Param

        dNi = dMolesSpecies(iFirstSpecies + ia - 1)
        dNj = dMolesSpecies(iFirstSpecies + ib - 1)

        ! Separate the amount and composition parts of this binary interaction.
        ! B is the product of the two participating species amounts. C is their
        ! difference. Dividing each by total phase moles N gives the extensive
        ! amount factor A=B/N and the dimensionless composition contrast
        ! Delta=C/N used by the Redlich-Kister polynomial.
        dB     = dNi * dNj
        dC     = dNi - dNj
        dA     = dB / dN
        dDelta = dC / dN

        ! Evaluate the interaction polynomial L and its first two derivatives
        ! with respect to the composition contrast Delta. These values are
        ! reused when differentiating with respect to every species-mole pair.
        dL0 = dDelta**iExponent
        if (iExponent == 0) then
            dL1 = 0D0
            dL2 = 0D0
        elseif (iExponent == 1) then
            dL1 = 1D0
            dL2 = 0D0
        else
            dL1 = DFLOAT(iExponent) * dDelta**(iExponent-1)
            if (iExponent == 2) then
                dL2 = 2D0
            else
                dL2 = DFLOAT(iExponent*(iExponent-1)) * dDelta**(iExponent-2)
            end if
        end if

        ! Expand the second derivative into every local (i,j) mole direction.
        ! Species outside the binary pair still affect A and Delta through N.
        do i = 1, nPhaseSpecies
            do j = 1, nPhaseSpecies

                ! These zero-or-one flags record whether derivative direction i
                ! or j changes either species named by this binary parameter.
                ! They let the same formulas cover participating species and
                ! all other species, which still affect the term through N.
                dDelta_ai = 0D0
                dDelta_aj = 0D0
                dDelta_bi = 0D0
                dDelta_bj = 0D0
                if (i == ia) dDelta_ai = 1D0
                if (i == ib) dDelta_aj = 1D0
                if (j == ia) dDelta_bi = 1D0
                if (j == ib) dDelta_bj = 1D0

                ! Differentiate the amount product B and amount difference C
                ! first, then propagate those changes through A=B/N and
                ! Delta=C/N. Suffix a means differentiation in species direction
                ! i, suffix b means direction j, and suffix ab means the mixed
                ! second derivative in directions i and j.
                dBa  = dDelta_ai*dNj + dDelta_aj*dNi
                dBb  = dDelta_bi*dNj + dDelta_bj*dNi
                dBab = dDelta_ai*dDelta_bj + dDelta_aj*dDelta_bi

                dCa  = dDelta_ai - dDelta_aj
                dCb  = dDelta_bi - dDelta_bj
                dCab = 0D0

                dAa  = dBa/dN - dB/(dN*dN)
                dAb  = dBb/dN - dB/(dN*dN)
                dAab = dBab/dN - dBa/(dN*dN) - dBb/(dN*dN) + 2D0*dB/(dN*dN*dN)

                ! Derivatives of the dimensionless composition contrast Delta.
                dDa  = dCa/dN - dC/(dN*dN)
                dDb  = dCb/dN - dC/(dN*dN)
                dDab = dCab/dN - dCa/(dN*dN) - dCb/(dN*dN) + 2D0*dC/(dN*dN*dN)

                ! Combine curvature of the amount factor and composition
                ! polynomial by the second-order product and chain rules.
                ! The database coefficient then gives the contribution of this parameter
                ! contribution to Hessian entry (i,j).
                dHij = dExcessGibbsParam(iParam) * ( dAab*dL0 + (dAa*dDb + dAb*dDa + dA*dDab)*dL1 + dA*dDa*dDb*dL2 )

                dHess(i,j) = dHess(i,j) + dHij

            end do
        end do

    end do LOOP_Param

    !=========================================================================================================
    ! SECTION 3: NUMERICAL SYMMETRY
    !
    ! Equality of mixed derivatives makes the theoretical Hessian symmetric.
    ! Average only roundoff-level evaluation-order differences before the matrix
    ! enters diagnostic or local-response linear algebra.
    !=========================================================================================================
    do i = 1, nPhaseSpecies
        do j = i + 1, nPhaseSpecies
            dSym = 0.5D0 * (dHess(i,j) + dHess(j,i))
            dHess(i,j) = dSym
            dHess(j,i) = dSym
        end do
    end do

end subroutine CompExcessGibbsEnergyRKMP_unconstrained
