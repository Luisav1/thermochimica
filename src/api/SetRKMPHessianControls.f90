!-------------------------------------------------------------------------------------------------------------
!> \file    SetRKMPHessianControls.f90
!> \brief   Configure the experimental plain-RKMP curvature path without editing solver source.
!>
!> \details The requested alpha is an upper trust bound on the completed mapped
!!          matrix and right-hand-side correction.  It does not scale the local
!!          RKMP derivatives.  The setting persists across Thermochimica
!!          initialization until ResetRKMPHessianControls is called.
!-------------------------------------------------------------------------------------------------------------

subroutine SetRKMPHessianControls(lEnable, dAlphaMax, lDebug, iInfo)

    USE ModuleGEMSolver

    implicit none

    logical, intent(in)  :: lEnable, lDebug
    real(8), intent(in)  :: dAlphaMax
    integer, intent(out) :: iInfo

    iInfo = 0
    if ((dAlphaMax /= dAlphaMax) .OR. (dAlphaMax < 0D0) .OR. (dAlphaMax > 1D0)) then
        iInfo = 1
        return
    end if

    lRKMPHessianControlsConfigured = .TRUE.
    lRKMPHessianRequestedEnable = lEnable
    lRKMPHessianRequestedDebug = lDebug
    dRKMPHessianRequestedAlphaMax = dAlphaMax

end subroutine SetRKMPHessianControls


subroutine ResetRKMPHessianControls

    USE ModuleGEMSolver

    implicit none

    lRKMPHessianControlsConfigured = .FALSE.
    lRKMPHessianRequestedEnable = .FALSE.
    lRKMPHessianRequestedDebug = .FALSE.
    dRKMPHessianRequestedAlphaMax = 0.10D0

end subroutine ResetRKMPHessianControls


!> \brief Configure the numerical trust gates used to globalize RKMP corrections.
subroutine SetRKMPHessianTrustThresholds(dUpdateRatioCap, dDirectionCosineMin, &
                                         dDirectionDifferenceCap, dProgressAllowance, iInfo)

    USE ModuleGEMSolver

    implicit none

    real(8), intent(in)  :: dUpdateRatioCap, dDirectionCosineMin
    real(8), intent(in)  :: dDirectionDifferenceCap, dProgressAllowance
    integer, intent(out) :: iInfo

    iInfo = 0
    if ((dUpdateRatioCap < 1D0) .OR. (dDirectionCosineMin < -1D0) .OR. &
        (dDirectionCosineMin > 1D0) .OR. (dDirectionDifferenceCap < 0D0) .OR. &
        (dProgressAllowance < 1D0)) then
        iInfo = 1
        return
    end if

    dRKMPTrustUpdateRatioCap = dUpdateRatioCap
    dRKMPTrustDirectionCosineMin = dDirectionCosineMin
    dRKMPTrustDirectionDifferenceCap = dDirectionDifferenceCap
    dRKMPTrustProgressAllowance = dProgressAllowance

end subroutine SetRKMPHessianTrustThresholds


subroutine ResetRKMPHessianTrustThresholds

    USE ModuleGEMSolver

    implicit none

    dRKMPTrustEmergencyRatioCap = 1D6
    dRKMPTrustUpdateRatioCap = 1.25D0
    dRKMPTrustDirectionCosineMin = 0.90D0
    dRKMPTrustDirectionDifferenceCap = 0.50D0
    dRKMPTrustLocalNormThreshold = 5D-2
    dRKMPTrustProgressAllowance = 1.05D0
    dRKMPTrustGibbsActivationTolerance = 1D-6
    dRKMPTrustGibbsRetentionTolerance = 1D-4

end subroutine ResetRKMPHessianTrustThresholds


!> \brief Enable or suppress the compact end-of-solve RKMP audit records.
subroutine SetRKMPHessianReportOutput(lEnable)

    USE ModuleGEMSolver

    implicit none

    logical, intent(in) :: lEnable

    lRKMPHessianReportSummary = lEnable

end subroutine SetRKMPHessianReportOutput
