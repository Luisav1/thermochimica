!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleCEFUnconstrained.f90
!> \brief   Disconnected mole-space Hessian for the nonmagnetic plain-SUBL CEF forms in Thermochimica.
!>
!> \details This module evaluates a general CEF energy expressed through endmember reference terms,
!!          ideal sublattice mixing, and excess terms of the form P(y)L(c+q.y).  It depends only on generic
!!          arrays and does not access Thermochimica state, constrain a phase assemblage, or modify GEMNewton.
!!
!!          Conceptual pipeline:
!!          1. Endmember moles define the total phase amount and endmember fractions.
!!          2. Endmember occupancy maps those fractions to constituent site fractions y.
!!          3. Reference, ideal-mixing, and excess energies are evaluated in site-fraction space.
!!          4. Their site-space curvature B is contracted through the composition map D:
!!                H = transpose(D) B D / n.
!!          5. H is the unconstrained Hessian with respect to endmember moles. It is not a
!!             GEMNewton matrix contribution and does not constrain the phase assemblage.
!!
!!          Notation used below:
!!          - i labels an endmember, meaning one allowed constituent combination
!!            spanning all sublattices.
!!          - u labels one constituent on one sublattice after all sublattice
!!            constituent lists have been flattened into a single array.
!!          - y(u) is the fraction of sublattice sites occupied by constituent u.
!!          - n is the total number of endmember moles in the phase.
!!          - B is curvature with respect to site fractions y.
!!          - D(u,i) records how site fraction u responds when endmember i is
!!            perturbed, including the change in total phase moles.
!!          - H is the final curvature with respect to endmember mole amounts.
!!          - transpose(D) maps a site-space response back to endmember space.
!!
!!          In an excess term P(y)L(eta), P is a product of selected site
!!          fractions, eta is one local composition coordinate, and L is the
!!          interaction polynomial evaluated at that coordinate.
!!
!!          File map:
!!          1. Generic interaction description and public interface
!!          2. Public scalar-energy and mole-space Hessian evaluators
!!          3. Input checks and endmember-to-site composition mapping
!!          4. Generic site-space interaction and polynomial calculus
!-------------------------------------------------------------------------------------------------------------

module ModuleCEFUnconstrained

    implicit none
    private

    !=========================================================================================================
    ! SECTION 1: GENERIC CEF INTERACTION DESCRIPTION AND PUBLIC INTERFACE
    !
    ! The module consumes model-neutral arrays. The production SUBL decoder lives
    ! in the verification test so this thermodynamic core remains disconnected
    ! from ModuleThermo and reusable for controlled mathematical checks.
    !=========================================================================================================

    !> One generic excess interaction written as an occupancy product times a polynomial.
    !>
    !> For each site constituent u, nu(u) is its power in the occupancy product.
    !> The coefficients q(u) and constant c combine all site fractions into the
    !> scalar composition coordinate eta=c+sum(q(u)*y(u)). The stored polynomial
    !> coefficients then define L(eta).
    type, public :: CEFInteractionTerm
        real(8) :: dArgumentConstant = 0D0
        real(8), allocatable :: dExponent(:)
        real(8), allocatable :: dArgumentCoefficient(:)
        real(8), allocatable :: dPolynomialCoefficient(:)
    end type CEFInteractionTerm

    public :: CompCEFGibbsEnergyUnconstrained
    public :: CompCEFHessianUnconstrained
    public :: CompCEFInteractionValue

contains

    !=========================================================================================================
    ! SECTION 2: PUBLIC ENERGY AND HESSIAN EVALUATORS
    !
    ! The scalar routine evaluates the extensive reference, ideal, and excess
    ! energy blocks. The Hessian routine differentiates those blocks in
    ! site-fraction space and maps their curvature back to endmember-mole space.
    !=========================================================================================================

    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate the extensive disconnected CEF Gibbs energy.
    !>
    !> \param[in]  dMoles Endmember mole amounts.
    !> \param[in]  iOccupancy Site-variable index occupied by each (sublattice,endmember) pair.
    !> \param[in]  iSiteSublattice Sublattice owning each site variable.
    !> \param[in]  dSublatticeMultiplicity Number of sites represented by each sublattice.
    !> \param[in]  dReferenceEnergy Dimensionless reference energy of each endmember.
    !> \param[in]  dIdealScale Dimensionless ideal-mixing scale, normally one in Thermochimica units.
    !> \param[in]  tInteraction Excess interaction terms in generic P(y)L(eta) form.
    !> \param[out] dGibbs Total extensive Gibbs energy.
    !> \param[out] dGibbsReference Extensive reference contribution.
    !> \param[out] dGibbsIdeal Extensive ideal-mixing contribution.
    !> \param[out] dGibbsExcess Extensive excess contribution.
    !> \param[out] iInfo Zero on success; nonzero when dimensions or the local composition are invalid.
    !---------------------------------------------------------------------------------------------------------
    subroutine CompCEFGibbsEnergyUnconstrained(dMoles,iOccupancy,iSiteSublattice, &
        dSublatticeMultiplicity,dReferenceEnergy,dIdealScale,tInteraction,dGibbs, &
        dGibbsReference,dGibbsIdeal,dGibbsExcess,iInfo)

        real(8), intent(in) :: dMoles(:), dSublatticeMultiplicity(:), dReferenceEnergy(:), dIdealScale
        integer, intent(in) :: iOccupancy(:,:), iSiteSublattice(:)
        type(CEFInteractionTerm), intent(in) :: tInteraction(:)
        real(8), intent(out) :: dGibbs, dGibbsReference, dGibbsIdeal, dGibbsExcess
        integer, intent(out) :: iInfo

        integer :: i, s, u
        real(8) :: dMonomial, dN
        real(8), allocatable :: dSiteFraction(:)

        dGibbs = 0D0
        dGibbsReference = 0D0
        dGibbsIdeal = 0D0
        dGibbsExcess = 0D0

        call CheckCEFInputs(dMoles,iOccupancy,iSiteSublattice,dSublatticeMultiplicity, &
            dReferenceEnergy,tInteraction,dSiteFraction,dN,iInfo)
        if (iInfo /= 0) return

        ! A reference endmember represents one simultaneous constituent choice
        ! on every sublattice. Its probability is therefore the product of the
        ! corresponding site fractions.
        do i = 1, SIZE(dMoles)
            dMonomial = 1D0
            do s = 1, SIZE(dSublatticeMultiplicity)
                dMonomial = dMonomial*dSiteFraction(iOccupancy(s,i))
            end do
            dGibbsReference = dGibbsReference + dReferenceEnergy(i)*dMonomial
        end do

        ! Ideal CEF mixing occurs independently on each sublattice. Multiplicity
        ! scales the contribution by the number of sites of that sublattice.
        do u = 1, SIZE(dSiteFraction)
            s = iSiteSublattice(u)
            dGibbsIdeal = dGibbsIdeal + dIdealScale*dSublatticeMultiplicity(s)* &
                dSiteFraction(u)*DLOG(dSiteFraction(u))
        end do

        ! Each decoded excess term is an occupancy prefactor P(y) multiplied by
        ! a polynomial in one linear local-composition coordinate eta.
        do i = 1, SIZE(tInteraction)
            dGibbsExcess = dGibbsExcess + CompCEFInteractionValue(dSiteFraction,tInteraction(i))
        end do

        dGibbsReference = dN*dGibbsReference
        dGibbsIdeal = dN*dGibbsIdeal
        dGibbsExcess = dN*dGibbsExcess
        dGibbs = dGibbsReference + dGibbsIdeal + dGibbsExcess

    end subroutine CompCEFGibbsEnergyUnconstrained


    !---------------------------------------------------------------------------------------------------------
    !> \brief Compute the unconstrained endmember-mole Hessian of the disconnected CEF energy.
    !>
    !> \details B first describes curvature between pairs of constituent site
    !!          fractions. D maps a perturbation of endmember i into the resulting
    !!          changes of every site fraction u. The contraction
    !!          H=transpose(D)*B*D/n therefore expresses that same curvature in
    !!          endmember-mole coordinates. In
    !!          D(u,i)=occupancy_indicator(u,i)-y(u), the indicator is one when
    !!          endmember i contains constituent u and zero otherwise. Optional
    !!          outputs expose the separate reference, ideal, and excess blocks.
    !>
    !> \param[in]  dMoles Endmember mole amounts.
    !> \param[in]  iOccupancy Site-variable index occupied by each (sublattice,endmember) pair.
    !> \param[in]  iSiteSublattice Sublattice owning each site variable.
    !> \param[in]  dSublatticeMultiplicity Number of sites represented by each sublattice.
    !> \param[in]  dReferenceEnergy Dimensionless reference energy of each endmember.
    !> \param[in]  dIdealScale Dimensionless ideal-mixing scale, normally one in Thermochimica units.
    !> \param[in]  tInteraction Excess interaction terms in generic P(y)L(eta) form.
    !> \param[out] dHessian Total mole-space Hessian.
    !> \param[out] iInfo Zero on success; nonzero when dimensions or the local composition are invalid.
    !> \param[out] dHessianReference Optional reference-energy Hessian block.
    !> \param[out] dHessianIdeal Optional ideal-mixing Hessian block.
    !> \param[out] dHessianExcess Optional excess-energy Hessian block.
    !---------------------------------------------------------------------------------------------------------
    subroutine CompCEFHessianUnconstrained(dMoles,iOccupancy,iSiteSublattice, &
        dSublatticeMultiplicity,dReferenceEnergy,dIdealScale,tInteraction,dHessian,iInfo, &
        dHessianReference,dHessianIdeal,dHessianExcess)

        real(8), intent(in) :: dMoles(:), dSublatticeMultiplicity(:), dReferenceEnergy(:), dIdealScale
        integer, intent(in) :: iOccupancy(:,:), iSiteSublattice(:)
        type(CEFInteractionTerm), intent(in) :: tInteraction(:)
        real(8), intent(out) :: dHessian(:,:)
        integer, intent(out) :: iInfo
        real(8), intent(out), optional :: dHessianReference(:,:), dHessianIdeal(:,:), dHessianExcess(:,:)

        integer :: i, j, s, u
        real(8) :: dN
        real(8), allocatable :: dBExcess(:,:), dBIdeal(:,:), dBReference(:,:), dD(:,:), dExponent(:)
        real(8), allocatable :: dHExcess(:,:), dHIdeal(:,:), dHReference(:,:), dSiteFraction(:)

        dHessian = 0D0
        if (PRESENT(dHessianReference)) dHessianReference = 0D0
        if (PRESENT(dHessianIdeal)) dHessianIdeal = 0D0
        if (PRESENT(dHessianExcess)) dHessianExcess = 0D0

        call CheckCEFInputs(dMoles,iOccupancy,iSiteSublattice,dSublatticeMultiplicity, &
            dReferenceEnergy,tInteraction,dSiteFraction,dN,iInfo)
        if (iInfo /= 0) return
        if ((SIZE(dHessian,1) /= SIZE(dMoles)) .OR. (SIZE(dHessian,2) /= SIZE(dMoles))) then
            iInfo = 20
            return
        end if
        if (PRESENT(dHessianReference)) then
            if ((SIZE(dHessianReference,1) /= SIZE(dMoles)) .OR. &
                (SIZE(dHessianReference,2) /= SIZE(dMoles))) then
                iInfo = 21
                return
            end if
        end if
        if (PRESENT(dHessianIdeal)) then
            if ((SIZE(dHessianIdeal,1) /= SIZE(dMoles)) .OR. &
                (SIZE(dHessianIdeal,2) /= SIZE(dMoles))) then
                iInfo = 22
                return
            end if
        end if
        if (PRESENT(dHessianExcess)) then
            if ((SIZE(dHessianExcess,1) /= SIZE(dMoles)) .OR. &
                (SIZE(dHessianExcess,2) /= SIZE(dMoles))) then
                iInfo = 23
                return
            end if
        end if

        allocate(dBReference(SIZE(dSiteFraction),SIZE(dSiteFraction)), &
            dBIdeal(SIZE(dSiteFraction),SIZE(dSiteFraction)), &
            dBExcess(SIZE(dSiteFraction),SIZE(dSiteFraction)), &
            dD(SIZE(dSiteFraction),SIZE(dMoles)), dExponent(SIZE(dSiteFraction)), &
            dHReference(SIZE(dMoles),SIZE(dMoles)), dHIdeal(SIZE(dMoles),SIZE(dMoles)), &
            dHExcess(SIZE(dMoles),SIZE(dMoles)))
        dBReference = 0D0
        dBIdeal = 0D0
        dBExcess = 0D0

        ! Differentiate the endmember probability products in site space. The
        ! reference block is generally curved in y even though the final
        ! extensive reference energy is linear in independent endmember moles.
        do i = 1, SIZE(dMoles)
            dExponent = 0D0
            do s = 1, SIZE(dSublatticeMultiplicity)
                dExponent(iOccupancy(s,i)) = dExponent(iOccupancy(s,i)) + 1D0
            end do
            call AddMonomialHessian(dSiteFraction,dExponent,dReferenceEnergy(i),dBReference)
        end do

        ! The second derivative of y*log(y) is 1/y, so ideal site-space
        ! curvature is diagonal before the endmember-composition contraction.
        do u = 1, SIZE(dSiteFraction)
            s = iSiteSublattice(u)
            dBIdeal(u,u) = dIdealScale*dSublatticeMultiplicity(s)/dSiteFraction(u)
        end do

        do i = 1, SIZE(tInteraction)
            call AddInteractionHessian(dSiteFraction,tInteraction(i),dBExcess)
        end do

        ! D describes how every site fraction changes when one endmember mole is
        ! perturbed. Here u labels a site constituent and i labels an endmember.
        ! Start with -y(u), which accounts for dilution as total phase moles
        ! increase, then add one when endmember i actually contains constituent u.
        do i = 1, SIZE(dMoles)
            do u = 1, SIZE(dSiteFraction)
                dD(u,i) = -dSiteFraction(u)
            end do
            do s = 1, SIZE(dSublatticeMultiplicity)
                dD(iOccupancy(s,i),i) = dD(iOccupancy(s,i),i) + 1D0
            end do
        end do

        ! Contract site-space curvature into endmember-mole space. Division by
        ! total phase moles follows from differentiating normalized fractions.
        dHReference = MATMUL(TRANSPOSE(dD),MATMUL(dBReference,dD))/dN
        dHIdeal = MATMUL(TRANSPOSE(dD),MATMUL(dBIdeal,dD))/dN
        dHExcess = MATMUL(TRANSPOSE(dD),MATMUL(dBExcess,dD))/dN
        dHessian = dHReference + dHIdeal + dHExcess

        ! Remove only multiplication-order asymmetry; this does not alter the mathematical Hessian.
        do i = 1, SIZE(dMoles)
            do j = i + 1, SIZE(dMoles)
                dHessian(i,j) = 0.5D0*(dHessian(i,j)+dHessian(j,i))
                dHessian(j,i) = dHessian(i,j)
            end do
        end do

        if (PRESENT(dHessianReference)) dHessianReference = dHReference
        if (PRESENT(dHessianIdeal)) dHessianIdeal = dHIdeal
        if (PRESENT(dHessianExcess)) dHessianExcess = dHExcess

    end subroutine CompCEFHessianUnconstrained

    !=========================================================================================================
    ! SECTION 3: INPUT CHECKS AND ENDMEMBER-TO-SITE COMPOSITION MAPPING
    !
    ! CEF logarithms and 1/y curvature require a strictly positive interior
    ! state. Invalid topology and boundary compositions are rejected explicitly;
    ! the module never clips or silently renormalizes the supplied state.
    !=========================================================================================================

    !---------------------------------------------------------------------------------------------------------
    !> \brief Verify the generic CEF topology and construct constituent site fractions.
    !>
    !> \details Endmember mole fractions are projected through iOccupancy. Each
    !!          endmember contributes its fraction to exactly one constituent on
    !!          every sublattice, producing the local variables used by all energy
    !!          and curvature formulas.
    !---------------------------------------------------------------------------------------------------------

    subroutine CheckCEFInputs(dMoles,iOccupancy,iSiteSublattice,dSublatticeMultiplicity, &
        dReferenceEnergy,tInteraction,dSiteFraction,dN,iInfo)

        real(8), intent(in) :: dMoles(:), dSublatticeMultiplicity(:), dReferenceEnergy(:)
        integer, intent(in) :: iOccupancy(:,:), iSiteSublattice(:)
        type(CEFInteractionTerm), intent(in) :: tInteraction(:)
        real(8), allocatable, intent(out) :: dSiteFraction(:)
        real(8), intent(out) :: dN
        integer, intent(out) :: iInfo

        integer :: i, s, u
        real(8), allocatable :: dX(:)

        iInfo = 0
        dN = 0D0
        if ((SIZE(dMoles) <= 0) .OR. (SIZE(iSiteSublattice) <= 0) .OR. &
            (SIZE(dSublatticeMultiplicity) <= 0)) then
            iInfo = 1
            return
        end if
        if ((SIZE(dReferenceEnergy) /= SIZE(dMoles)) .OR. &
            (SIZE(iOccupancy,1) /= SIZE(dSublatticeMultiplicity)) .OR. &
            (SIZE(iOccupancy,2) /= SIZE(dMoles))) then
            iInfo = 2
            return
        end if
        if (ANY(dMoles < 0D0) .OR. ANY(dSublatticeMultiplicity <= 0D0)) then
            iInfo = 3
            return
        end if
        dN = SUM(dMoles)
        if (dN <= 0D0) then
            iInfo = 4
            return
        end if
        if (ANY(iSiteSublattice < 1) .OR. ANY(iSiteSublattice > SIZE(dSublatticeMultiplicity)) .OR. &
            ANY(iOccupancy < 1) .OR. ANY(iOccupancy > SIZE(iSiteSublattice))) then
            iInfo = 5
            return
        end if
        do i = 1, SIZE(dMoles)
            do s = 1, SIZE(dSublatticeMultiplicity)
                if (iSiteSublattice(iOccupancy(s,i)) /= s) then
                    iInfo = 6
                    return
                end if
            end do
        end do
        do s = 1, SIZE(dSublatticeMultiplicity)
            if (COUNT(iSiteSublattice == s) == 0) then
                iInfo = 7
                return
            end if
        end do

        allocate(dSiteFraction(SIZE(iSiteSublattice)),dX(SIZE(dMoles)))
        dX = dMoles/dN
        dSiteFraction = 0D0
        do i = 1, SIZE(dMoles)
            do s = 1, SIZE(dSublatticeMultiplicity)
                u = iOccupancy(s,i)
                dSiteFraction(u) = dSiteFraction(u) + dX(i)
            end do
        end do
        if (ANY(dSiteFraction <= 0D0)) then
            iInfo = 8
            return
        end if

        do i = 1, SIZE(tInteraction)
            if (.NOT.ALLOCATED(tInteraction(i)%dExponent) .OR. &
                .NOT.ALLOCATED(tInteraction(i)%dArgumentCoefficient) .OR. &
                .NOT.ALLOCATED(tInteraction(i)%dPolynomialCoefficient)) then
                iInfo = 9
                return
            end if
            if ((SIZE(tInteraction(i)%dExponent) /= SIZE(dSiteFraction)) .OR. &
                (SIZE(tInteraction(i)%dArgumentCoefficient) /= SIZE(dSiteFraction)) .OR. &
                (SIZE(tInteraction(i)%dPolynomialCoefficient) <= 0)) then
                iInfo = 10
                return
            end if
            if (ANY(tInteraction(i)%dExponent < 0D0)) then
                iInfo = 11
                return
            end if
        end do

    end subroutine CheckCEFInputs

    !=========================================================================================================
    ! SECTION 4: GENERIC SITE-SPACE INTERACTION CALCULUS
    !
    ! Production binary, ternary, and coupled SUBL parameters are decoded into
    ! the common form phi(y) = P(y)L(eta), eta = c + q.y. These helpers evaluate
    ! that form and its site-space curvature without knowing the production
    ! parameter family from which it came.
    !=========================================================================================================

    !> \brief Evaluate one generic P(y)L(eta) term for decoder identity checks.
    real(8) function CompCEFInteractionValue(dY,tTerm)

        real(8), intent(in) :: dY(:)
        type(CEFInteractionTerm), intent(in) :: tTerm
        real(8) :: dEta, dL, dLp, dLpp

        dEta = tTerm%dArgumentConstant + DOT_PRODUCT(tTerm%dArgumentCoefficient,dY)
        call CompPolynomial(tTerm%dPolynomialCoefficient,dEta,dL,dLp,dLpp)
        CompCEFInteractionValue = PRODUCT(dY**tTerm%dExponent)*dL

    end function CompCEFInteractionValue


    !---------------------------------------------------------------------------------------------------------
    !> \brief Add the site-space Hessian of a pure occupancy monomial.
    !>
    !> \details Reference endmember probabilities and excess prefactors are
    !!          products of site fractions. This helper supplies the curvature of
    !!          such a product before the D contraction into mole space.
    !---------------------------------------------------------------------------------------------------------
    subroutine AddMonomialHessian(dY,dExponent,dCoefficient,dB)

        real(8), intent(in) :: dY(:), dExponent(:), dCoefficient
        real(8), intent(inout) :: dB(:,:)
        integer :: u, v
        real(8) :: dP, dPuv

        dP = PRODUCT(dY**dExponent)
        do u = 1, SIZE(dY)
            do v = 1, SIZE(dY)
                dPuv = dP*dExponent(u)*dExponent(v)/(dY(u)*dY(v))
                if (u == v) dPuv = dPuv - dP*dExponent(u)/(dY(u)*dY(u))
                dB(u,v) = dB(u,v) + dCoefficient*dPuv
            end do
        end do

    end subroutine AddMonomialHessian


    !---------------------------------------------------------------------------------------------------------
    !> \brief Add the site-space Hessian of one P(y)L(c+q.y) excess term.
    !>
    !> \details The four contributions are the curvature of P, two cross terms
    !!          coupling the gradients of P and L, and the curvature of L along
    !!          its linear composition coordinate.
    !---------------------------------------------------------------------------------------------------------
    subroutine AddInteractionHessian(dY,tTerm,dB)

        real(8), intent(in) :: dY(:)
        type(CEFInteractionTerm), intent(in) :: tTerm
        real(8), intent(inout) :: dB(:,:)
        integer :: u, v
        real(8) :: dEta, dL, dLp, dLpp, dP, dPu, dPv, dPuv

        dP = PRODUCT(dY**tTerm%dExponent)
        dEta = tTerm%dArgumentConstant + DOT_PRODUCT(tTerm%dArgumentCoefficient,dY)
        call CompPolynomial(tTerm%dPolynomialCoefficient,dEta,dL,dLp,dLpp)
        do u = 1, SIZE(dY)
            dPu = dP*tTerm%dExponent(u)/dY(u)
            do v = 1, SIZE(dY)
                dPv = dP*tTerm%dExponent(v)/dY(v)
                dPuv = dP*tTerm%dExponent(u)*tTerm%dExponent(v)/(dY(u)*dY(v))
                if (u == v) dPuv = dPuv - dP*tTerm%dExponent(u)/(dY(u)*dY(u))
                dB(u,v) = dB(u,v) + dPuv*dL + &
                    dPu*dLp*tTerm%dArgumentCoefficient(v) + &
                    dPv*dLp*tTerm%dArgumentCoefficient(u) + &
                    dP*dLpp*tTerm%dArgumentCoefficient(u)*tTerm%dArgumentCoefficient(v)
            end do
        end do

    end subroutine AddInteractionHessian


    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate an interaction polynomial and its first two eta derivatives.
    !>
    !> \details Coefficients are stored in ascending power order. Horner-style
    !!          evaluation keeps the value and derivative calculations compact
    !!          without duplicating a production interaction formula.
    !---------------------------------------------------------------------------------------------------------
    subroutine CompPolynomial(dCoefficient,dEta,dValue,dFirst,dSecond)

        real(8), intent(in) :: dCoefficient(:), dEta
        real(8), intent(out) :: dValue, dFirst, dSecond
        integer :: k

        dValue = 0D0
        dFirst = 0D0
        dSecond = 0D0
        do k = SIZE(dCoefficient), 1, -1
            dValue = dValue*dEta + dCoefficient(k)
        end do
        do k = SIZE(dCoefficient), 2, -1
            dFirst = dFirst*dEta + DFLOAT(k-1)*dCoefficient(k)
        end do
        do k = SIZE(dCoefficient), 3, -1
            dSecond = dSecond*dEta + DFLOAT((k-1)*(k-2))*dCoefficient(k)
        end do

    end subroutine CompPolynomial

end module ModuleCEFUnconstrained
