!-------------------------------------------------------------------------------------------------------------
!> \file    SetMQMQAHessianControls.f90
!> \brief   Configure the default-off fixed-alpha MQMQA GEM correction experiment.
!>
!> \details Alpha weights the completed reduced-GEM matrix and residual correction. It does not scale the local
!!          SUBG or SUBQ Hessian. Settings persist across Thermochimica initialization until explicitly reset.
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

end subroutine SetMQMQAHessianControls


!> \brief Restore the persistent MQMQA control request to its default-off state.
subroutine ResetMQMQAHessianControls

    USE ModuleGEMSolver

    implicit none

    lMQMQAHessianControlsConfigured = .FALSE.
    lMQMQAHessianRequestedEnable = .FALSE.
    dMQMQAHessianRequestedAlpha = 0D0

end subroutine ResetMQMQAHessianControls
