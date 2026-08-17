!-------------------------------------------------------------------------------------------------------------
!> \file    SetMQMQAHessianControls.f90
!> \brief   Configure the default-off fixed or adaptive MQMQA GEM correction experiment.
!>
!> \details These APIs choose how much of the completed reduced-GEM `(deltaA,deltaB)` correction may enter the
!!          Newton system.  Alpha does not scale the local SUBG or SUBQ Hessian.  Settings persist across
!!          Thermochimica initialization until explicitly reset, while per-calculation histories and readiness
!!          state are reset by `InitGEMSolver`.
!-------------------------------------------------------------------------------------------------------------

subroutine SetMQMQAHessianControls(lEnable,dAlpha,iInfo)

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleGEMSolver

    implicit none

    logical, intent(in) :: lEnable
    real(8), intent(in) :: dAlpha
    integer, intent(out) :: iInfo

    iInfo = 0
    if ((.NOT. IEEE_IS_FINITE(dAlpha)) .OR. (dAlpha < 0D0) .OR. (dAlpha > 1D0)) then
        iInfo = 1
        return
    end if

    ! Commit all requested controls only after every input has passed, so an
    ! invalid setter call cannot partially replace the previous configuration.
    lMQMQAHessianControlsConfigured = .TRUE.
    lMQMQAHessianRequestedEnable = lEnable
    dMQMQAHessianRequestedAlpha = dAlpha
    lMQMQAHessianAdaptiveControlsConfigured = .FALSE.
    lMQMQAHessianRequestedAdaptiveEnable = .FALSE.
    dMQMQAHessianRequestedAlphaMax = 0D0

end subroutine SetMQMQAHessianControls


!> \brief Restore the persistent MQMQA control request to its default-off state.
subroutine ResetMQMQAHessianControls

    USE ModuleGEMSolver

    implicit none

    lMQMQAHessianControlsConfigured = .FALSE.
    lMQMQAHessianRequestedEnable = .FALSE.
    dMQMQAHessianRequestedAlpha = 0D0

end subroutine ResetMQMQAHessianControls


!-------------------------------------------------------------------------------------------------------------
!> \brief Configure default-off adaptive MQMQA alpha trust without changing the fixed-alpha API contract.
!>
!> \details `dAlphaMax` bounds a descending trial list applied to the completed MQMQA reduced-GEM correction.
!!          The adaptive solver normally tries the largest allowed candidate first, but may return a reduced or
!!          zero-alpha solve when readiness or safety checks reject larger values.  Alpha zero is the untouched
!!          historical GEM Newton system.  A valid call atomically selects adaptive mode and disables the
!!          persistent fixed-alpha request; an invalid call preserves the previous configuration.
!-------------------------------------------------------------------------------------------------------------
subroutine SetMQMQAHessianAdaptiveControls(lEnable,dAlphaMax,iInfo)

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleGEMSolver

    implicit none

    logical, intent(in) :: lEnable
    real(8), intent(in) :: dAlphaMax
    integer, intent(out) :: iInfo

    iInfo = 0
    if ((.NOT. IEEE_IS_FINITE(dAlphaMax)) .OR. (dAlphaMax < 0D0) .OR. (dAlphaMax > 1D0)) then
        iInfo = 1
        return
    end if

    lMQMQAHessianAdaptiveControlsConfigured = .TRUE.
    lMQMQAHessianRequestedAdaptiveEnable = lEnable
    dMQMQAHessianRequestedAlphaMax = dAlphaMax
    lMQMQAHessianControlsConfigured = .FALSE.
    lMQMQAHessianRequestedEnable = .FALSE.
    dMQMQAHessianRequestedAlpha = 0D0

end subroutine SetMQMQAHessianAdaptiveControls


!> \brief Restore only the persistent adaptive MQMQA request to its default-off state.
subroutine ResetMQMQAHessianAdaptiveControls

    USE ModuleGEMSolver

    implicit none

    lMQMQAHessianAdaptiveControlsConfigured = .FALSE.
    lMQMQAHessianRequestedAdaptiveEnable = .FALSE.
    dMQMQAHessianRequestedAlphaMax = 0D0

end subroutine ResetMQMQAHessianAdaptiveControls
