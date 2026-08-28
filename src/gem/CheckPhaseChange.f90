
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckPhaseChange.f90
    !> \brief   Check whether a particular phase change is appropriate for further consideration.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      GEMNewton.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code
    !   07/31/2011      M.H.A. Piro         Clean up code: remove unnecessary variables, update variable names
    !   10/25/2011      M.H.A. Piro         Clean up code: Modules, simplify code
    !   04/26/2012      M.H.A. Piro         Convert to Gibbs energy minimization solver.
    !   05/08/2013      M.H.A. Piro         Specify a tolerance when the phase assemblage has never
    !                                        changed after a certain number of iterations (say 50).
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check whether a particular phase assemblage is a valid
    !! candidate.  It is possible for a particular combination of phases to yield non-real values when evaluating
    !! the Jacobian.  For example, suppose there is 1 mol of uranium in the system and the only uranium containing
    !! phase is removed from the assemblage.  Clearly, there must be at least one phase containing uranium for the
    !! system to be defined.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[out]  lPhasePass  A logical variable indicating whether the candidate phase assemblage is
    !!                           appropriate (i.e., .TRUE.) or not (i.e., .FALSE.).
    !> \param[out]  INFO        An integer scalar used to identify a successful exit or an error by the
    !!                           GEMNewton.f90 subroutine.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine CheckPhaseChange(lPhasePass,INFO)

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_DIVIDE_BY_ZERO, IEEE_GET_HALTING_MODE, &
        IEEE_INVALID, IEEE_IS_FINITE, IEEE_OVERFLOW, IEEE_SET_HALTING_MODE
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleGEMNewtonDiagnosticCapture, ONLY: CaptureMQMQAPhaseChangeCheck, &
        lMQMQADiagnosticMinimumNormStudy

    implicit none

    integer  ::  i, iActiveRank, j, INFO, nMiscPhases
    real(8)  ::  dMinActiveAmount, dTemp
    logical  ::  lPhasePass


    ! Initialize variables:
    INFO        = 0
    nMiscPhases = 0
    lPhasePass  = .TRUE.

    ! Count phases:
    j = nConPhases
    nConPhases  = 0
    CountCon: do i = 1, j
        if (iAssemblage(i) > 0) then
            nConPhases  = nConPhases  + 1
        else
            exit CountCon
        end if
    end do CountCon

    j = nSolnPhases
    nSolnPhases = 0
    CountSoln: do i = nElements, nElements + 1 - j, -1
        if (iAssemblage(i) < 0) then
            nSolnPhases = nSolnPhases + 1
        else
            exit CountSoln
        end if
    end do CountSoln

    ! Count the number of miscible phases:
    do i = 1, nSolnPhases
        j = -iAssemblage(nElements - i + 1)
        if (lMiscibility(j)) nMiscPhases = nMiscPhases + 1
    end do

    ! Determine tolerance:
    !if ((iterGlobal - iterLast < 10).OR.(iterLast == 0).OR. (nMiscPhases > 0)) then
    if ((iterGlobal > 50).AND.(iterLast == 0)) then
        dTemp = 1D50
    !elseif ((iterGlobal - iterLast < 20).OR.(iterLast == 0).OR. (nMiscPhases > 0)) then
    elseif ((iterGlobal - iterLast <= 30).OR.(iterLast == 0).OR. (nMiscPhases > 0)) then
        dTemp = 1D14
    !elseif (iterGlobal - iterLast <= 50) then
    elseif (iterGlobal - iterLast <= 100) then
        dTemp = 1D20
    else
        dTemp = 1D100
    end if

    ! Establish the Hessian matrix and compute the direction vector:
    call GEMNewton(INFO)

    ! Reinitialize variables:
    lRevertSystem = .FALSE.

    ! Check if this candidate phase assemblage is appropriate:
    if (INFO /= 0) lPhasePass = .FALSE.

    ! If the maximum value of the direction vector is above an arbitrarily large number, then the phase fails:
    if (MAXVAL(DABS(dUpdateVar)) >= dTemp) lPhasePass = .FALSE.

    dMinActiveAmount = HUGE(1D0)
    if (nConPhases > 0) dMinActiveAmount = MIN(dMinActiveAmount,MINVAL(dMolesPhase(1:nConPhases)))
    if (nSolnPhases > 0) dMinActiveAmount = MIN(dMinActiveAmount, &
        MINVAL(dMolesPhase(nElements-nSolnPhases+1:nElements)))
    if (nConPhases+nSolnPhases == 0) dMinActiveAmount = 0D0
    iActiveRank = 0
    if (lMQMQADiagnosticMinimumNormStudy) call AnalyzeActivePhaseRank(iActiveRank)
    call CaptureMQMQAPhaseChangeCheck(iterGlobal,nElements,nChargedConstraints,nSolnPhases,nConPhases, &
        INFO,lPhasePass,iActiveRank,MAXVAL(DABS(dUpdateVar)),dTemp,dMinActiveAmount,dTolerance(7),iAssemblage)

    return

contains

    !---------------------------------------------------------------------------------------------------------
    !> \brief Report the numerical rank of the active element-by-phase stoichiometry block.
    !---------------------------------------------------------------------------------------------------------
    subroutine AnalyzeActivePhaseRank(nRank)

        integer, intent(out) :: nRank
        integer :: iElement, iInfoSVD, iPhase, iSolution, lWork, nActive
        logical :: lHaltDivide, lHaltInvalid, lHaltOverflow
        real(8) :: dScale, dToleranceRank
        real(8), allocatable :: C(:,:), dSingular(:), dU(:,:), dVT(:,:), dWork(:)

        nRank = 0
        nActive = nConPhases+nSolnPhases
        if ((nElements <= 0) .OR. (nActive <= 0)) return
        allocate(C(nElements,nActive),dSingular(MIN(nElements,nActive)),dU(1,1),dVT(1,1))
        C = 0D0
        do iPhase = 1,nSolnPhases
            iSolution = -iAssemblage(nElements-iPhase+1)
            do iElement = 1,nElements
                C(iElement,iPhase) = dEffStoichSolnPhase(iSolution,iElement) * &
                    dMolesPhase(nElements-iPhase+1)
            end do
        end do
        do iPhase = 1,nConPhases
            C(:,nSolnPhases+iPhase) = dStoichSpecies(iAssemblage(iPhase),1:nElements)
        end do
        if (.NOT. ALL(IEEE_IS_FINITE(C))) return
        dScale = MAXVAL(DABS(C))
        if ((.NOT. IEEE_IS_FINITE(dScale)) .OR. (dScale <= 0D0)) return
        C = C/dScale
        lWork = MAX(1,5*MAX(nElements,nActive))
        allocate(dWork(lWork))
        call IEEE_GET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO,lHaltDivide)
        call IEEE_GET_HALTING_MODE(IEEE_INVALID,lHaltInvalid)
        call IEEE_GET_HALTING_MODE(IEEE_OVERFLOW,lHaltOverflow)
        call IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO,.FALSE.)
        call IEEE_SET_HALTING_MODE(IEEE_INVALID,.FALSE.)
        call IEEE_SET_HALTING_MODE(IEEE_OVERFLOW,.FALSE.)
        call DGESVD('N','N',nElements,nActive,C,nElements,dSingular,dU,1,dVT,1,dWork,lWork,iInfoSVD)
        call IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO,lHaltDivide)
        call IEEE_SET_HALTING_MODE(IEEE_INVALID,lHaltInvalid)
        call IEEE_SET_HALTING_MODE(IEEE_OVERFLOW,lHaltOverflow)
        if (iInfoSVD /= 0) return
        dToleranceRank = DBLE(MAX(nElements,nActive))*EPSILON(1D0)*dSingular(1)
        nRank = COUNT(dSingular > dToleranceRank)

    end subroutine AnalyzeActivePhaseRank

end subroutine CheckPhaseChange
