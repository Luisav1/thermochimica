!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleCEFUnconstrained.f90
!> \brief   Disconnected mole-space Hessian for Thermochimica's nonmagnetic plain-SUBL CEF forms.
!>
!> \details This module evaluates a general CEF energy expressed through endmember reference terms,
!!          ideal sublattice mixing, and excess terms of the form P(y)L(c+q.y).  It depends only on generic
!!          arrays and does not access Thermochimica state, constrain a phase assemblage, or modify GEMNewton.
!-------------------------------------------------------------------------------------------------------------

module ModuleCEFUnconstrained

    implicit none
    private

    !> One excess term phi(y) = product(y(u)**nu(u)) * L(c + dot_product(q,y)).
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

        do i = 1, SIZE(dMoles)
            dMonomial = 1D0
            do s = 1, SIZE(dSublatticeMultiplicity)
                dMonomial = dMonomial*dSiteFraction(iOccupancy(s,i))
            end do
            dGibbsReference = dGibbsReference + dReferenceEnergy(i)*dMonomial
        end do

        do u = 1, SIZE(dSiteFraction)
            s = iSiteSublattice(u)
            dGibbsIdeal = dGibbsIdeal + dIdealScale*dSublatticeMultiplicity(s)* &
                dSiteFraction(u)*DLOG(dSiteFraction(u))
        end do

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
    !> \details The site-space block B is contracted as H = D^T B D / n, where
    !!          D(u,i) = delta(i occupies u) - y(u).  Optional outputs expose the reference, ideal, and excess
    !!          Hessian blocks for diagnostic comparisons without adding an analytic first-derivative API.
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

        ! Each reference endmember is a product containing its one constituent from every sublattice.
        do i = 1, SIZE(dMoles)
            dExponent = 0D0
            do s = 1, SIZE(dSublatticeMultiplicity)
                dExponent(iOccupancy(s,i)) = dExponent(iOccupancy(s,i)) + 1D0
            end do
            call AddMonomialHessian(dSiteFraction,dExponent,dReferenceEnergy(i),dBReference)
        end do

        do u = 1, SIZE(dSiteFraction)
            s = iSiteSublattice(u)
            dBIdeal(u,u) = dIdealScale*dSublatticeMultiplicity(s)/dSiteFraction(u)
        end do

        do i = 1, SIZE(tInteraction)
            call AddInteractionHessian(dSiteFraction,tInteraction(i),dBExcess)
        end do

        do i = 1, SIZE(dMoles)
            do u = 1, SIZE(dSiteFraction)
                dD(u,i) = -dSiteFraction(u)
            end do
            do s = 1, SIZE(dSublatticeMultiplicity)
                dD(iOccupancy(s,i),i) = dD(iOccupancy(s,i),i) + 1D0
            end do
        end do

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


    !> \brief Evaluate one generic P(y)L(eta) term for decoder identity checks.
    real(8) function CompCEFInteractionValue(dY,tTerm)

        real(8), intent(in) :: dY(:)
        type(CEFInteractionTerm), intent(in) :: tTerm
        real(8) :: dEta, dL, dLp, dLpp

        dEta = tTerm%dArgumentConstant + DOT_PRODUCT(tTerm%dArgumentCoefficient,dY)
        call CompPolynomial(tTerm%dPolynomialCoefficient,dEta,dL,dLp,dLpp)
        CompCEFInteractionValue = PRODUCT(dY**tTerm%dExponent)*dL

    end function CompCEFInteractionValue


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
