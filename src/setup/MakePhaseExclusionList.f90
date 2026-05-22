
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    MakePhaseExclusionList.f90
    !> \brief   Builds or update the list of thermodynamic phases to exclude from or include to the calculation.
    !> \author  M. Poschmann
    !> \date    May 31, 2022
    ! TODO: sa
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/06/2022      M. Poschmann        Don't automatically exclude dummies from the calculation
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details This subroutine manages which thermodynamic phases should be excluded from calculations. If the user
    !! specifies an “excluded except” list — meaning only certain phases are to remain active, the routine scans all
    !! defined solution and pure condensed phases, and automatically adds every other phase to the exclusion list. This
    !! ensures the system only considers the allowed phases while all others are ignored. If no “excluded except” list
    !! is provided, the subroutine leaves the existing exclusion list unchanged, allowing users to manually specify 
    !! which individual phases to exclude.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nPhasesExcluded          Integer representing the number of phases to be excluded from the calculation.
    ! nPhasesExcludedExcept    Integer representing the number of phases to be included (exceptions to exclusion).
    ! cPhasesExcluded          Character array containing the names of phases to be excluded from
    ! cPhasesExcludedExcept    Character array containing the names of phases to be included (exceptions to exclusion)
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine MakePhaseExclusionList

    USE ModuleThermoIO
    USE ModuleParseCS

    implicit none

    integer :: i, j

    ! If there is an "excluded except" list, then add all other phases to exclusion list
    if (nPhasesExcludedExcept > 0) then
        ! Solution phases
        loop_checkExclusionSolution: do i = 1, nSolnPhasesSysCS
            ! Check if phase is on exception list
            do j = 1, nPhasesExcludedExcept
                if (cSolnPhaseNameCS(i) == cPhasesExcludedExcept(j)) cycle loop_checkExclusionSolution
            end do
            ! Check if phase is on exclusion list
            do j = 1, nPhasesExcluded
                if (cSolnPhaseNameCS(i) == cPhasesExcluded(j)) cycle loop_checkExclusionSolution
            end do
            ! If not, add to exclusion list
            nPhasesExcluded = nPhasesExcluded + 1
            cPhasesExcluded(nPhasesExcluded) = cSolnPhaseNameCS(i)
        end do loop_checkExclusionSolution

        ! Pure condensed phases
        loop_checkExclusionPureCondensed: do i = nSpeciesPhaseCS(nSolnPhasesSysCS) + 1, nSpeciesCS
            ! Check if dummy - don't exclude dummies automatically
            if (iPhaseCS(i) == -1) cycle loop_checkExclusionPureCondensed
            ! Check if phase is on exception list
            do j = 1, nPhasesExcludedExcept
                if (cSpeciesNameCS(i) == cPhasesExcludedExcept(j)) cycle loop_checkExclusionPureCondensed
            end do
            ! Check if phase is on exclusion list
            do j = 1, nPhasesExcluded
                if (cSpeciesNameCS(i) == cPhasesExcluded(j)) cycle loop_checkExclusionPureCondensed
            end do
            ! If not, add to exclusion list
            nPhasesExcluded = nPhasesExcluded + 1
            cPhasesExcluded(nPhasesExcluded) = cSpeciesNameCS(i)
        end do loop_checkExclusionPureCondensed
    end if


end subroutine MakePhaseExclusionList
