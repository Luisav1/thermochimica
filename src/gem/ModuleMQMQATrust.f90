!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleMQMQATrust.f90
!> \brief   State-free candidate construction and mixed-variable update checks for MQMQA globalization.
!>
!> \details This module contains the state-free part of MQ-4D.  It does not decide when the nonlinear solver is
!!          ready for curvature, assemble phase corrections, solve GEM, or alter Thermochimica state.  Instead,
!!          it answers two narrower questions for GEMNewton: which alpha values should be tried, and whether a
!!          trial correction/update is small and well aligned relative to the untouched historical solve.
!!
!!          MQMQA alpha weights a completed reduced-GEM `(deltaA,deltaB)` correction pair; it never weights the
!!          local analytic Hessian.  Alpha zero therefore means the historical GEM Newton system, not a
!!          first-order thermodynamic calculation.  Update comparisons keep element potentials,
!!          solution-phase logarithmic increments, and pure-phase amount displacements in separate groups so
!!          unlike units are not combined into one Euclidean trust metric.
!-------------------------------------------------------------------------------------------------------------
module ModuleMQMQATrust

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE

    implicit none
    private

    integer, parameter, public :: MQMQA_TRUST_ACCEPTED = 0
    integer, parameter, public :: MQMQA_TRUST_INVALID_INPUT = 1
    integer, parameter, public :: MQMQA_TRUST_UPDATE_REJECTED = 2
    integer, parameter, public :: MQMQA_TRUST_DIRECTION_REJECTED = 3

    public :: BuildMQMQAAlphaCandidateList, EvaluateMQMQACorrectionRatio, EvaluateMQMQAUpdateTrust

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Build descending alpha candidates bounded by the requested maximum.
    !---------------------------------------------------------------------------------------------------------
    subroutine BuildMQMQAAlphaCandidateList(dAlphaMaxInput,dCandidates,nCandidates)

        real(8), intent(in) :: dAlphaMaxInput
        real(8), intent(out) :: dCandidates(:)
        integer, intent(out) :: nCandidates

        integer :: i
        real(8) :: dAlphaMax
        real(8), parameter :: dStandard(4)=[1D0,1D-1,1D-2,1D-3]

        dCandidates = 0D0
        nCandidates = 0
        if ((.NOT. IEEE_IS_FINITE(dAlphaMaxInput)) .OR. (SIZE(dCandidates) < 5)) return

        dAlphaMax = DMAX1(0D0,DMIN1(1D0,dAlphaMaxInput))
        if (dAlphaMax > 0D0) then
            nCandidates = 1
            dCandidates(1) = dAlphaMax
            do i = 1,SIZE(dStandard)
                if (dStandard(i) >= dAlphaMax*(1D0-1D-12)) cycle
                nCandidates = nCandidates+1
                dCandidates(nCandidates) = dStandard(i)
            end do
        end if
        nCandidates = nCandidates+1
        dCandidates(nCandidates) = 0D0

    end subroutine BuildMQMQAAlphaCandidateList


    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate the matrix and residual correction sizes without applying either correction.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateMQMQACorrectionRatio(dABase,dBBase,nElement,dDeltaA,dDeltaB,dAlpha,dRatioCap, &
        lAccepted,dRatioA,dRatioB)

        real(8), intent(in) :: dABase(:,:), dBBase(:), dDeltaA(:,:), dDeltaB(:), dAlpha, dRatioCap
        integer, intent(in) :: nElement
        logical, intent(out) :: lAccepted
        real(8), intent(out) :: dRatioA, dRatioB

        real(8) :: dBaseNormA, dBaseNormB

        lAccepted = .FALSE.
        dRatioA = HUGE(1D0)
        dRatioB = HUGE(1D0)
        if ((nElement < 1) .OR. (SIZE(dABase,1) < nElement) .OR. (SIZE(dABase,2) < nElement) .OR. &
            (SIZE(dBBase) < nElement) .OR. (SIZE(dDeltaA,1) /= nElement) .OR. &
            (SIZE(dDeltaA,2) /= nElement) .OR. (SIZE(dDeltaB) /= nElement) .OR. &
            (.NOT. IEEE_IS_FINITE(dAlpha)) .OR. (dAlpha < 0D0) .OR. (dAlpha > 1D0) .OR. &
            (.NOT. IEEE_IS_FINITE(dRatioCap)) .OR. (dRatioCap < 0D0) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dABase))) .OR. (.NOT. ALL(IEEE_IS_FINITE(dBBase))) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dDeltaA))) .OR. (.NOT. ALL(IEEE_IS_FINITE(dDeltaB)))) return

        dBaseNormA = SQRT(SUM(dABase(1:nElement,1:nElement)**2))
        dBaseNormB = SQRT(SUM(dBBase(1:nElement)**2))
        dRatioA = dAlpha*SQRT(SUM(dDeltaA**2))/DMAX1(dBaseNormA,1D-30)
        dRatioB = dAlpha*SQRT(SUM(dDeltaB**2))/DMAX1(dBaseNormB,1D-30)
        lAccepted = DMAX1(dRatioA,dRatioB) <= dRatioCap

    end subroutine EvaluateMQMQACorrectionRatio


    !---------------------------------------------------------------------------------------------------------
    !> \brief Compare corrected and alpha-zero solved displacements in three physically distinct groups.
    !>
    !> \details Each group is compared only with its own alpha-zero displacement. Empty groups are neutral.
    !!          If both the baseline and trial norms lie below the group-specific resolved-step floor
    !!          (1e-8 for dimensionless element-potential and first-order logarithmic solution increments,
    !!          1e-12 mol for pure-phase amounts), their relative ratio and direction are treated as neutral.
    !!          This prevents ratios of two negligible steps from rejecting an otherwise resolved correction.
    !!          These floors and the remaining thresholds are numerical globalization heuristics rather than
    !!          thermodynamic identities.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateMQMQAUpdateTrust(dStepBase,dStepTrial,nElement,nSolution,nPure,dUpdateRatioCap, &
        dDirectionCosineMin,dDirectionDifferenceCap,lAccepted,dUpdateRatio,dDirectionCosine, &
        dDirectionDifference,dBaseNormGroup,dTrialNormGroup,iStatus)

        real(8), intent(in) :: dStepBase(:), dStepTrial(:)
        integer, intent(in) :: nElement, nSolution, nPure
        real(8), intent(in) :: dUpdateRatioCap, dDirectionCosineMin, dDirectionDifferenceCap
        logical, intent(out) :: lAccepted
        real(8), intent(out) :: dUpdateRatio(3), dDirectionCosine(3), dDirectionDifference(3)
        real(8), intent(out) :: dBaseNormGroup(3), dTrialNormGroup(3)
        integer, intent(out) :: iStatus

        integer :: iFirst(3), iLast(3), iGroup, nRequired
        real(8) :: dBaseNorm, dTrialNorm, dDifferenceNorm, dScale
        real(8), parameter :: dNormFloor(3)=[1D-8,1D-8,1D-12]

        lAccepted = .FALSE.
        iStatus = MQMQA_TRUST_INVALID_INPUT
        dUpdateRatio = HUGE(1D0)
        dDirectionCosine = -1D0
        dDirectionDifference = HUGE(1D0)
        dBaseNormGroup = HUGE(1D0)
        dTrialNormGroup = HUGE(1D0)

        nRequired = nElement+nSolution+nPure
        if ((nElement < 0) .OR. (nSolution < 0) .OR. (nPure < 0) .OR. &
            (SIZE(dStepBase) < nRequired) .OR. (SIZE(dStepTrial) < nRequired) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dStepBase(1:nRequired)))) .OR. &
            (.NOT. ALL(IEEE_IS_FINITE(dStepTrial(1:nRequired)))) .OR. &
            (.NOT. IEEE_IS_FINITE(dUpdateRatioCap)) .OR. (dUpdateRatioCap < 1D0) .OR. &
            (.NOT. IEEE_IS_FINITE(dDirectionCosineMin)) .OR. (dDirectionCosineMin < -1D0) .OR. &
            (dDirectionCosineMin > 1D0) .OR. (.NOT. IEEE_IS_FINITE(dDirectionDifferenceCap)) .OR. &
            (dDirectionDifferenceCap < 0D0)) return

        iFirst = [1,nElement+1,nElement+nSolution+1]
        iLast = [nElement,nElement+nSolution,nRequired]
        do iGroup = 1,3
            if (iLast(iGroup) < iFirst(iGroup)) then
                dUpdateRatio(iGroup) = 1D0
                dDirectionCosine(iGroup) = 1D0
                dDirectionDifference(iGroup) = 0D0
                dBaseNormGroup(iGroup) = 0D0
                dTrialNormGroup(iGroup) = 0D0
                cycle
            end if

            dBaseNorm = SQRT(SUM(dStepBase(iFirst(iGroup):iLast(iGroup))**2))
            dTrialNorm = SQRT(SUM(dStepTrial(iFirst(iGroup):iLast(iGroup))**2))
            dDifferenceNorm = SQRT(SUM((dStepTrial(iFirst(iGroup):iLast(iGroup))- &
                dStepBase(iFirst(iGroup):iLast(iGroup)))**2))
            dBaseNormGroup(iGroup) = dBaseNorm
            dTrialNormGroup(iGroup) = dTrialNorm
            dScale = DMAX1(dBaseNorm,dNormFloor(iGroup))
            dUpdateRatio(iGroup) = dTrialNorm/dScale
            dDirectionDifference(iGroup) = dDifferenceNorm/dScale
            if ((dBaseNorm <= dNormFloor(iGroup)) .AND. (dTrialNorm <= dNormFloor(iGroup))) then
                dUpdateRatio(iGroup) = 1D0
                dDirectionCosine(iGroup) = 1D0
                dDirectionDifference(iGroup) = 0D0
            else if ((dBaseNorm > dNormFloor(iGroup)) .AND. (dTrialNorm > dNormFloor(iGroup))) then
                dDirectionCosine(iGroup) = DOT_PRODUCT( &
                    dStepBase(iFirst(iGroup):iLast(iGroup)),dStepTrial(iFirst(iGroup):iLast(iGroup)))/ &
                    (dBaseNorm*dTrialNorm)
                dDirectionCosine(iGroup) = DMAX1(-1D0,DMIN1(1D0,dDirectionCosine(iGroup)))
            else
                dDirectionCosine(iGroup) = 0D0
            end if
        end do

        if (ANY(dUpdateRatio > dUpdateRatioCap)) then
            iStatus = MQMQA_TRUST_UPDATE_REJECTED
            return
        end if
        if (ANY(dDirectionCosine < dDirectionCosineMin) .OR. &
            ANY(dDirectionDifference > dDirectionDifferenceCap)) then
            iStatus = MQMQA_TRUST_DIRECTION_REJECTED
            return
        end if

        lAccepted = .TRUE.
        iStatus = MQMQA_TRUST_ACCEPTED

    end subroutine EvaluateMQMQAUpdateTrust

end module ModuleMQMQATrust
