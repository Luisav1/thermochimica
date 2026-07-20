    !---------------------------------------------------------------------------------------------------------
    !
    ! Purpose:
    ! --------
    ! Compute an unconstrained local Hessian for RKMP excess Gibbs energy terms in a single solution phase.
    !
    ! Scope / current limitations:
    !   * This routine currently handles binary RKMP parameters only (iRegularParam(:,1) == 2).
    !   * The contribution computed here is EXCESS-only curvature for staged integration.
    !   * Magnetic RKMPM second-order terms are not included here.
    !
    ! Inputs:
    !   iSolnIndex   Absolute solution-phase index.
    !
    ! Output:
    !   dHess(:,:)   Local phase Hessian in phase-species coordinates (1:nPhaseSpecies,1:nPhaseSpecies).
    !                The caller must provide sufficient storage.
    !
    ! Notation (matching derivation variables):
    !   dN      = total moles in this solution phase
    !   dB      = n_i * n_j
    !   dC      = n_i - n_j
    !   dA      = dB / dN
    !   dDelta  = dC / dN
    !   dL0,dL1,dL2 = L(Delta), L'(Delta), L''(Delta) for RK polynomial exponent
    !
    !---------------------------------------------------------------------------------------------------------
    
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

    LOOP_Param: do iParam = nParamPhase(iSolnIndex-1)+1, nParamPhase(iSolnIndex)

        ! Implementation supports binary RKMP terms only.
        ! Ternary and higher-order cases are intentionally skipped.
        if (iRegularParam(iParam,1) /= 2) cycle LOOP_Param

        ! Local species IDs participating in this binary interaction parameter.
        ia = iRegularParam(iParam,2)
        ib = iRegularParam(iParam,3)

        if ((ia < 1) .OR. (ia > nPhaseSpecies) .OR. (ib < 1) .OR. (ib > nPhaseSpecies)) cycle LOOP_Param

        ! Polynomial exponent for this parameter's Redlich-Kister term.
        iExponent = iRegularParam(iParam,4)
        if (iExponent < 0) cycle LOOP_Param

        dNi = dMolesSpecies(iFirstSpecies + ia - 1)
        dNj = dMolesSpecies(iFirstSpecies + ib - 1)

        ! Build compact intermediates once per parameter and reuse in all (i,j) entries.
        dB     = dNi * dNj
        dC     = dNi - dNj
        dA     = dB / dN
        dDelta = dC / dN

        ! Evaluate L(Delta), L'(Delta), L''(Delta) once; reused in every Hessian entry.
        ! Similar to RKMP excess Gibbs energy contribution, but with derivatives of L(Delta) included according to chain rule.
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

        do i = 1, nPhaseSpecies
            do j = 1, nPhaseSpecies

                ! Kronecker-delta flags capture whether loop indices i/j match parameter species ia/ib.
                ! This is how partial derivatives are expanded without branch-heavy symbolic code.
                dDelta_ai = 0D0
                dDelta_aj = 0D0
                dDelta_bi = 0D0
                dDelta_bj = 0D0
                if (i == ia) dDelta_ai = 1D0
                if (i == ib) dDelta_aj = 1D0
                if (j == ia) dDelta_bi = 1D0
                if (j == ib) dDelta_bj = 1D0

                ! First/second derivatives of compact terms B, C, A, Delta w.r.t local species moles.
                ! These derivatives are the building blocks of the final Hessian expression.
                dBa  = dDelta_ai*dNj + dDelta_aj*dNi
                dBb  = dDelta_bi*dNj + dDelta_bj*dNi
                dBab = dDelta_ai*dDelta_bj + dDelta_aj*dDelta_bi

                dCa  = dDelta_ai - dDelta_aj
                dCb  = dDelta_bi - dDelta_bj
                dCab = 0D0

                dAa  = dBa/dN - dB/(dN*dN)
                dAb  = dBb/dN - dB/(dN*dN)
                dAab = dBab/dN - dBa/(dN*dN) - dBb/(dN*dN) + 2D0*dB/(dN*dN*dN)

                ! The delta derivatives
                dDa  = dCa/dN - dC/(dN*dN)
                dDb  = dCb/dN - dC/(dN*dN)
                dDab = dCab/dN - dCa/(dN*dN) - dCb/(dN*dN) + 2D0*dC/(dN*dN*dN)

                ! Add this parameter's contribution to Hessian entry (i,j).
                ! Expression follows chain-rule expansion of the RKMP excess term.
                ! Multiplying by dExcessGibbsParam(iParam) applies the parameter value's specific contribution to the final contribution.
                dHij = dExcessGibbsParam(iParam) * ( dAab*dL0 + (dAa*dDb + dAb*dDa + dA*dDab)*dL1 + dA*dDa*dDb*dL2 )

                dHess(i,j) = dHess(i,j) + dHij

            end do
        end do

    end do LOOP_Param

    ! Numerical guard: theoretical Hessian is symmetric, but roundoff/order of operations can introduce
    ! tiny asymmetry. Average upper/lower entries so downstream linear algebra sees a symmetric matrix.
    do i = 1, nPhaseSpecies
        do j = i + 1, nPhaseSpecies
            dSym = 0.5D0 * (dHess(i,j) + dHess(j,i))
            dHess(i,j) = dSym
            dHess(j,i) = dSym
        end do
    end do

end subroutine CompExcessGibbsEnergyRKMP_unconstrained
