!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleFiniteDifferenceVerification.f90
!> \brief   Shared order-aware finite-difference assessment for derivative verification tests.
!>
!> \details A small error at one perturbation size is not, by itself, evidence
!!          that a finite-difference approximation is converging for the
!!          expected mathematical reason. This module evaluates a complete
!!          decreasing-step sweep and separates three questions:
!!
!!          1. Are all reported errors finite?
!!          2. Does a pre-roundoff region exhibit the expected algebraic order?
!!          3. Does at least one step satisfy the requested accuracy threshold?
!!
!!          The final row of a sweep has no observed order because no finer
!!          point exists. Orders are also unavailable when either error is zero,
!!          non-finite, or already at a conservative double-precision floor.
!-------------------------------------------------------------------------------------------------------------

module ModuleFiniteDifferenceVerification

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE

    implicit none

    private

    real(8), parameter, public :: FD_ORDER_SECOND_MIN = 1.7D0
    real(8), parameter, public :: FD_ORDER_SECOND_MAX = 2.3D0
    real(8), parameter, public :: FD_ORDER_FOURTH_MIN = 3.2D0
    real(8), parameter, public :: FD_ORDER_FOURTH_MAX = 4.8D0

    !> Summary of one decreasing-step finite-difference sweep.
    type, public :: FDSweepAssessment
        logical :: lFinite = .FALSE.
        logical :: lAccuracy = .FALSE.
        logical :: lOrderRegion = .FALSE.
        logical :: lRoundoffUpturn = .FALSE.
        logical :: lPassed = .FALSE.
        integer :: iBest = 0
        integer :: iOrderStart = 0
        real(8) :: dBestError = HUGE(1D0)
        real(8) :: dObservedOrderMin = HUGE(1D0)
        real(8) :: dObservedOrderMax = -HUGE(1D0)
    end type FDSweepAssessment

    public :: AssessFDSweep
    public :: ComputeObservedOrders
    public :: ComputeVectorErrorMetrics
    public :: VectorTwoNormFD

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Compute local observed order between every adjacent pair of sweep points.
    !>
    !> \param[in]  dStep       Strictly decreasing positive perturbation sizes.
    !> \param[in]  dError      Nonnegative errors corresponding to dStep.
    !> \param[out] dOrder      Local observed orders; the final value is zero and unavailable.
    !> \param[out] lAvailable  True where dOrder can be interpreted safely.
    !---------------------------------------------------------------------------------------------------------
    subroutine ComputeObservedOrders(dStep,dError,dOrder,lAvailable)

        real(8), intent(in) :: dStep(:), dError(:)
        real(8), intent(out) :: dOrder(:)
        logical, intent(out) :: lAvailable(:)

        integer :: i, n
        real(8) :: dRoundoffFloor

        n = SIZE(dStep)
        dOrder = 0D0
        lAvailable = .FALSE.
        if ((SIZE(dError) /= n) .OR. (SIZE(dOrder) /= n) .OR. (SIZE(lAvailable) /= n)) return

        ! Errors below this level are too close to double-precision arithmetic
        ! noise for a logarithmic ratio to carry useful convergence information.
        dRoundoffFloor = 100D0*EPSILON(1D0)
        do i = 1, n-1
            if ((dStep(i) <= dStep(i+1)) .OR. (dStep(i+1) <= 0D0)) cycle
            if ((dError(i) <= dRoundoffFloor) .OR. (dError(i+1) <= dRoundoffFloor)) cycle
            if ((.NOT.IEEE_IS_FINITE(dError(i))) .OR. (.NOT.IEEE_IS_FINITE(dError(i+1)))) cycle
            dOrder(i) = DLOG(dError(i)/dError(i+1))/DLOG(dStep(i)/dStep(i+1))
            lAvailable(i) = IEEE_IS_FINITE(dOrder(i))
        end do

    end subroutine ComputeObservedOrders


    !---------------------------------------------------------------------------------------------------------
    !> \brief Apply accuracy and consecutive-order requirements to one finite-difference sweep.
    !>
    !> \details The search stops at the minimum-error point, so a later
    !!          roundoff-driven increase cannot invalidate an otherwise valid
    !!          truncation region. Two consecutive accepted orders require
    !!          three successively smaller errors and therefore provide stronger
    !!          evidence than a single favourable ratio.
    !---------------------------------------------------------------------------------------------------------
    subroutine AssessFDSweep(dStep,dError,dOrderLower,dOrderUpper,dAccuracyTolerance,tResult, &
        dOrder,lOrderAvailable)

        real(8), intent(in) :: dStep(:), dError(:)
        real(8), intent(in) :: dOrderLower, dOrderUpper, dAccuracyTolerance
        type(FDSweepAssessment), intent(out) :: tResult
        real(8), intent(out), optional :: dOrder(:)
        logical, intent(out), optional :: lOrderAvailable(:)

        integer :: i, n
        logical, allocatable :: lAvailableLocal(:)
        real(8), allocatable :: dOrderLocal(:)

        tResult = FDSweepAssessment()
        n = SIZE(dStep)
        if ((n < 3) .OR. (SIZE(dError) /= n)) return

        allocate(dOrderLocal(n),lAvailableLocal(n))
        call ComputeObservedOrders(dStep,dError,dOrderLocal,lAvailableLocal)
        if (PRESENT(dOrder)) then
            if (SIZE(dOrder) == n) dOrder = dOrderLocal
        end if
        if (PRESENT(lOrderAvailable)) then
            if (SIZE(lOrderAvailable) == n) lOrderAvailable = lAvailableLocal
        end if

        tResult%lFinite = ALL(IEEE_IS_FINITE(dStep)) .AND. ALL(IEEE_IS_FINITE(dError)) .AND. &
            ALL(dStep > 0D0) .AND. ALL(dError >= 0D0)
        if (.NOT.tResult%lFinite) then
            deallocate(dOrderLocal,lAvailableLocal)
            return
        end if

        tResult%iBest = MINLOC(dError,1)
        tResult%dBestError = dError(tResult%iBest)
        tResult%lAccuracy = tResult%dBestError <= dAccuracyTolerance

        do i = 1, MIN(n-2,tResult%iBest-2)
            if ((.NOT.lAvailableLocal(i)) .OR. (.NOT.lAvailableLocal(i+1))) cycle
            if (.NOT.((dError(i+1) < dError(i)) .AND. (dError(i+2) < dError(i+1)))) cycle
            if ((dOrderLocal(i) < dOrderLower) .OR. (dOrderLocal(i) > dOrderUpper)) cycle
            if ((dOrderLocal(i+1) < dOrderLower) .OR. (dOrderLocal(i+1) > dOrderUpper)) cycle
            tResult%lOrderRegion = .TRUE.
            tResult%iOrderStart = i
            tResult%dObservedOrderMin = MIN(dOrderLocal(i),dOrderLocal(i+1))
            tResult%dObservedOrderMax = MAX(dOrderLocal(i),dOrderLocal(i+1))
            exit
        end do

        if (tResult%iBest < n) then
            do i = tResult%iBest+1, n
                if (dError(i) > dError(i-1)) then
                    tResult%lRoundoffUpturn = .TRUE.
                    exit
                end if
            end do
        end if

        tResult%lPassed = tResult%lFinite .AND. tResult%lAccuracy .AND. tResult%lOrderRegion
        deallocate(dOrderLocal,lAvailableLocal)

    end subroutine AssessFDSweep


    !---------------------------------------------------------------------------------------------------------
    !> \brief Compute normwise and componentwise errors for one vector derivative approximation.
    !>
    !> \details Normwise scaling uses the larger full-vector magnitude and one.
    !!          Componentwise scaling uses the larger local component magnitude
    !!          and the full-vector scale. Consequently, components whose exact
    !!          response is nearly zero are judged by an absolute tolerance tied
    !!          to the scale of the complete response rather than an unstable
    !!          relative error.
    !---------------------------------------------------------------------------------------------------------
    subroutine ComputeVectorErrorMetrics(dApproximate,dExact,dNormAbsolute,dNormScaled, &
        dMaxAbsolute,dMaxScaled,iWorst)

        real(8), intent(in) :: dApproximate(:), dExact(:)
        real(8), intent(out) :: dNormAbsolute, dNormScaled, dMaxAbsolute, dMaxScaled
        integer, intent(out) :: iWorst

        real(8) :: dFullScale
        real(8), allocatable :: dAbsolute(:), dComponentScale(:)

        if ((SIZE(dApproximate) == 0) .OR. (SIZE(dApproximate) /= SIZE(dExact))) then
            dNormAbsolute = HUGE(1D0)
            dNormScaled = HUGE(1D0)
            dMaxAbsolute = HUGE(1D0)
            dMaxScaled = HUGE(1D0)
            iWorst = 0
            return
        end if

        allocate(dAbsolute(SIZE(dExact)),dComponentScale(SIZE(dExact)))
        dAbsolute = DABS(dApproximate-dExact)
        dNormAbsolute = VectorTwoNormFD(dApproximate-dExact)
        dFullScale = DMAX1(1D0,VectorTwoNormFD(dApproximate),VectorTwoNormFD(dExact))
        dNormScaled = dNormAbsolute/dFullScale
        dComponentScale = MAX(dFullScale,DABS(dApproximate),DABS(dExact))
        dMaxAbsolute = MAXVAL(dAbsolute)
        iWorst = MAXLOC(dAbsolute/dComponentScale,1)
        dMaxScaled = MAXVAL(dAbsolute/dComponentScale)
        deallocate(dAbsolute,dComponentScale)

    end subroutine ComputeVectorErrorMetrics


    real(8) function VectorTwoNormFD(dVector)

        real(8), intent(in) :: dVector(:)
        VectorTwoNormFD = SQRT(SUM(dVector*dVector))

    end function VectorTwoNormFD

end module ModuleFiniteDifferenceVerification
