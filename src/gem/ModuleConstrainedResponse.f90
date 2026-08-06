!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleConstrainedResponse.f90
!> \brief   Solve a local thermodynamic response subject to linear equality constraints.
!>
!> \details A phase composition cannot respond along directions that violate
!!          normalization or another local equality constraint. Given a local
!!          chemical-potential curvature H, constraint matrix C, and one or
!!          more forcing columns F, this module solves the symmetric bordered
!!          system
!!
!!              [ H  C^T ] [ X      ] = [ F ]
!!              [ C   0  ] [ Lambda ]   [ 0 ].
!!
!!          X contains the admissible composition responses. The module is
!!          deliberately model-independent: it does not decode thermodynamic
!!          data, alter global state, or assemble GEM matrix contributions.
!-------------------------------------------------------------------------------------------------------------

module ModuleConstrainedResponse

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE

    implicit none
    private

    public :: SolveConstrainedResponse

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Solve a symmetric constrained-response system for multiple forcing directions.
    !>
    !> \param[in]  dCurvature  Square local chemical-potential response matrix H.
    !> \param[in]  dConstraint Equality-constraint matrix C; each row is one constraint.
    !> \param[in]  dForcing    Forcing matrix F; each column is one perturbation direction.
    !> \param[out] dResponse   Composition response X with the same shape as dForcing.
    !> \param[out] iInfo       Zero on success; negative for invalid input; positive for a singular KKT solve.
    !> \param[out] dMultiplier Optional Lagrange-multiplier response for each constraint and forcing column.
    !>
    !> \details The sign of a Lagrange multiplier is conventional. Using the
    !!          same C in the lower-left block and C^T in the upper-right block
    !!          keeps the KKT matrix symmetric without changing the physical
    !!          composition response obtained from an equivalent mixed-sign
    !!          formulation.
    !---------------------------------------------------------------------------------------------------------
    subroutine SolveConstrainedResponse(dCurvature,dConstraint,dForcing,dResponse,iInfo,dMultiplier)

        real(8), intent(in) :: dCurvature(:,:), dConstraint(:,:), dForcing(:,:)
        real(8), intent(out) :: dResponse(:,:)
        integer, intent(out) :: iInfo
        real(8), intent(out), optional :: dMultiplier(:,:)

        integer :: nConstraint, nEquation, nForcing, nVariable
        integer, allocatable :: iPivot(:)
        real(8) :: dSymmetryScale
        real(8), allocatable :: dKKT(:,:), dRHS(:,:)

        iInfo = 0
        dResponse = 0D0
        if (PRESENT(dMultiplier)) dMultiplier = 0D0

        nVariable = SIZE(dCurvature,1)
        nConstraint = SIZE(dConstraint,1)
        nForcing = SIZE(dForcing,2)
        if ((nVariable < 2) .OR. (nConstraint < 1) .OR. (nForcing < 1)) then
            iInfo = -1
            return
        end if
        if ((SIZE(dCurvature,2) /= nVariable) .OR. &
            (SIZE(dConstraint,2) /= nVariable) .OR. &
            (SIZE(dForcing,1) /= nVariable) .OR. &
            (SIZE(dResponse,1) /= nVariable) .OR. &
            (SIZE(dResponse,2) /= nForcing)) then
            iInfo = -2
            return
        end if
        if (PRESENT(dMultiplier)) then
            if ((SIZE(dMultiplier,1) /= nConstraint) .OR. &
                (SIZE(dMultiplier,2) /= nForcing)) then
                iInfo = -3
                return
            end if
        end if
        if ((.NOT. ALL(IEEE_IS_FINITE(dCurvature))) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dConstraint))) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dForcing)))) then
            iInfo = -4
            return
        end if

        ! An asymmetric local matrix would make the advertised symmetric KKT
        ! interpretation false and usually indicates an upstream derivative or
        ! indexing defect. Do not repair it silently.
        dSymmetryScale = DMAX1(1D0,SQRT(SUM(dCurvature*dCurvature)))
        if (SQRT(SUM((dCurvature-TRANSPOSE(dCurvature))**2))/dSymmetryScale > 1D-12) then
            iInfo = -5
            return
        end if

        nEquation = nVariable+nConstraint
        allocate(dKKT(nEquation,nEquation),dRHS(nEquation,nForcing),iPivot(nEquation))
        dKKT = 0D0
        dRHS = 0D0
        dKKT(1:nVariable,1:nVariable) = dCurvature
        dKKT(1:nVariable,nVariable+1:nEquation) = TRANSPOSE(dConstraint)
        dKKT(nVariable+1:nEquation,1:nVariable) = dConstraint
        dRHS(1:nVariable,:) = dForcing

        call DGESV(nEquation,nForcing,dKKT,nEquation,iPivot,dRHS,nEquation,iInfo)
        if (iInfo == 0) then
            dResponse = dRHS(1:nVariable,:)
            if (PRESENT(dMultiplier)) dMultiplier = dRHS(nVariable+1:nEquation,:)
            if (.NOT. ALL(IEEE_IS_FINITE(dResponse))) iInfo = -6
            if (PRESENT(dMultiplier)) then
                if (.NOT. ALL(IEEE_IS_FINITE(dMultiplier))) iInfo = -6
            end if
        end if

        deallocate(dKKT,dRHS,iPivot)

    end subroutine SolveConstrainedResponse

end module ModuleConstrainedResponse
