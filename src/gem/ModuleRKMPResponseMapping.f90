!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleRKMPResponseMapping.f90
!> \brief   Shared constrained local-response solve used by the RKMP GEM mapper.
!>
!> \details The local composition fractions cannot vary independently because
!!          their infinitesimal changes must sum to zero.  This module forms the
!!          bordered linear system that combines the local chemical-potential
!!          curvature with that normalization condition.  Keeping the solve in
!!          one small module makes its singular-failure behavior directly
!!          testable without invoking the full nonlinear GEM solver.
!-------------------------------------------------------------------------------------------------------------

module ModuleRKMPResponseMapping

    implicit none

contains

    subroutine SolveRKMPConstrainedResponse(nSpecies, nForcing, dCurvature, dForcing, dResponse, iInfo)

        integer, intent(in)                    :: nSpecies, nForcing
        integer, intent(out)                   :: iInfo
        real(8), intent(in), dimension(:,:)    :: dCurvature, dForcing
        real(8), intent(out), dimension(:,:)   :: dResponse

        integer                                :: iSpecies, iDirection, nEquation
        integer, dimension(:), allocatable     :: iPivot
        real(8), dimension(:,:), allocatable   :: dKKT, dRHS

        iInfo = 0
        dResponse = 0D0
        if ((nSpecies <= 1) .OR. (nForcing <= 0)) then
            iInfo = -1
            return
        end if
        if ((SIZE(dCurvature,1) < nSpecies) .OR. (SIZE(dCurvature,2) < nSpecies) .OR. &
            (SIZE(dForcing,1) < nSpecies) .OR. (SIZE(dForcing,2) < nForcing) .OR. &
            (SIZE(dResponse,1) < nSpecies) .OR. (SIZE(dResponse,2) < nForcing)) then
            iInfo = -2
            return
        end if

        nEquation = nSpecies + 1
        allocate(dKKT(nEquation,nEquation), dRHS(nEquation,nForcing), iPivot(nEquation))
        dKKT = 0D0
        dRHS = 0D0

        dKKT(1:nSpecies,1:nSpecies) = dCurvature(1:nSpecies,1:nSpecies)
        do iSpecies = 1, nSpecies
            dKKT(iSpecies,nEquation) = -1D0
            dKKT(nEquation,iSpecies) = 1D0
        end do
        dRHS(1:nSpecies,1:nForcing) = dForcing(1:nSpecies,1:nForcing)

        call DGESV(nEquation, nForcing, dKKT, nEquation, iPivot, dRHS, nEquation, iInfo)
        if (iInfo == 0) then
            do iDirection = 1, nForcing
                dResponse(1:nSpecies,iDirection) = dRHS(1:nSpecies,iDirection)
            end do
        end if

        deallocate(dKKT, dRHS, iPivot)

    end subroutine SolveRKMPConstrainedResponse

end module ModuleRKMPResponseMapping
