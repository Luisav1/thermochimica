!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleGEMNewtonDiagnosticCapture.f90
!> \brief   Opt-in snapshot of the unreduced GEMNewton linear system for verification tests.
!>
!> \details MQMQA response mapping must be compared with the matrix and right-hand side that GEMNewton actually
!!          assembles, not only with a second copy of the source formulas in a test. This module provides a
!!          deliberately small diagnostic boundary: when capture is requested, GEMNewton copies its completed
!!          baseline system immediately before any experimental curvature correction or linear solve.
!!
!!          Capture is disabled by default, does not print, and never changes the supplied matrix or vector.
!-------------------------------------------------------------------------------------------------------------
module ModuleGEMNewtonDiagnosticCapture

    implicit none
    private

    logical, public :: lCaptureGEMNewtonSystem = .FALSE.
    logical, public :: lGEMNewtonSystemCaptured = .FALSE.
    integer, public :: nCapturedGEMNewtonVariables = 0
    real(8), allocatable, public :: dCapturedGEMNewtonA(:,:), dCapturedGEMNewtonB(:)

    public :: CaptureGEMNewtonSystem, ResetGEMNewtonDiagnosticCapture

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Copy one completed baseline Newton system when the opt-in request is active.
    !---------------------------------------------------------------------------------------------------------
    subroutine CaptureGEMNewtonSystem(dA,dB,nVar)

        integer, intent(in) :: nVar
        real(8), intent(in) :: dA(:,:), dB(:)

        if (.NOT. lCaptureGEMNewtonSystem) return
        if ((nVar <= 0) .OR. (SIZE(dA,1) < nVar) .OR. (SIZE(dA,2) < nVar) .OR. &
            (SIZE(dB) < nVar)) return

        if (allocated(dCapturedGEMNewtonA)) deallocate(dCapturedGEMNewtonA)
        if (allocated(dCapturedGEMNewtonB)) deallocate(dCapturedGEMNewtonB)
        allocate(dCapturedGEMNewtonA(nVar,nVar),dCapturedGEMNewtonB(nVar))
        dCapturedGEMNewtonA = dA(1:nVar,1:nVar)
        dCapturedGEMNewtonB = dB(1:nVar)
        nCapturedGEMNewtonVariables = nVar
        lGEMNewtonSystemCaptured = .TRUE.

    end subroutine CaptureGEMNewtonSystem


    !---------------------------------------------------------------------------------------------------------
    !> \brief Clear captured storage and restore the default-inactive diagnostic state.
    !---------------------------------------------------------------------------------------------------------
    subroutine ResetGEMNewtonDiagnosticCapture

        if (allocated(dCapturedGEMNewtonA)) deallocate(dCapturedGEMNewtonA)
        if (allocated(dCapturedGEMNewtonB)) deallocate(dCapturedGEMNewtonB)
        lCaptureGEMNewtonSystem = .FALSE.
        lGEMNewtonSystemCaptured = .FALSE.
        nCapturedGEMNewtonVariables = 0

    end subroutine ResetGEMNewtonDiagnosticCapture

end module ModuleGEMNewtonDiagnosticCapture
