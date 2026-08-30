
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GEMSolver.f90
    !> \brief   Gibbs Energy Minimization solver.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      Thermochimica.f90
    !> \sa      InitGEMSolver.f90
    !> \sa      GEMNewton.f90
    !> \sa      GEMLineSearch.f90
    !> \sa      CheckPhaseAssemblage.f90
    !> \sa      CheckConvergence.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the quantities of species and phases at thermodynamic
    !! equilibrium using the Gibbs Energy Minimization (GEM) method.  This subroutine uses values of
    !! dMolesPhase, dChemicalPotential and iAssemblage from the Leveling and PostLeveling subroutines as initial
    !! estimates for computation.
    !!
    !! The main subroutines used by this solver are summarized below:
    !! <table border="1" width="800">
    !! <tr>
    !!    <td> <b> File name </td> <td> Description </b> </td>
    !! </tr>
    !! <tr>
    !!    <td> InitGEMSolver.f90 </td>
    !!    <td> Initialize the GEMSolver by establishing the initial phase assemblage and composition.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckSysOnlyPureConPhases.f90 </td>
    !!    <td> Check the system if there are only pure condensed phases. The system may already be converged.</td>
    !! </tr>
    !! <tr>
    !!    <td> GEMNewton.f90 </td>
    !!    <td> Compute the direction vector using Newton's method.  </td>
    !! </tr>
    !! <tr>
    !!    <td> GEMLineSearch.f90 </td>
    !!    <td> Perform a line search along the direction vector.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckPhaseAssemblage.f90 </td>
    !!    <td> Check if the phase assemblage needs to be adjusted.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckConvergence.f90 </td>
    !!    <td> Check if the system has converged.  </td>
    !! </tr>
    !! </table>
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nConPhases            The number of pure condensed phases in the assemblage
    ! nSolnPhases           The number of solution phases in the assemblage
    ! nSolnPhasesSys        The number of solution phases in the system
    ! iAssemblage           Integer vector containing the indices of phases in the assemblage
    !                        (1:nConphases represent pure condensed phases and (nElements-nSolnPhases:nSolnPhases)
    !                        represent solution phases.
    ! INFOThermo            An integer scalar identifying whether the program exits successfully or if
    !                        it encounters an error.
    ! INFO                  An integer scalar identifying an error from LAPACK.  This is used by the GEMNewton
    !                        subroutine to indicate whether there is a singularity in the Hessian matrix.
    ! lConverged            A logical variable indicating whether the code has convered (.TRUE.) or not (.FALSE.).
    ! lRevertSystem         A logical scalar indicating whether the system should be reverted to a previously
    !                        successful phase assemblage.
    ! dTolerance            A double real vector representing numerical tolerances (defined in InitThermo.f90).
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine GEMSolver

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleGEMNewtonDiagnosticCapture, ONLY: BeginMQMQATrajectoryAttempt, &
        CaptureMQMQATrajectoryPoint, CaptureMQMQAFixedPointState, &
        CaptureMQMQARecoveryPreStep, CaptureMQMQARecoveryNewton, &
        CaptureMQMQARecoveryAssemblage, &
        lCaptureMQMQAFixedPointAtConvergence, lMQMQAFixedPointCaptureAttempted, &
        iMQMQAFixedPointCaptureInfo, dMQMQAFixedPointRecomputeDifference, &
        dMQMQAFixedPointRestorationError, &
        dCapturedMQMQAFixedPointSolvedUpdate, &
        lCaptureGEMNewtonSystem, lCaptureGEMNewtonCorrectedSystem, &
        lCaptureFirstGEMNewtonCorrectionPair

    implicit none

    integer::   INFO, iCaptureInfo
    logical :: lAdaptiveSave, lRevertSave, lUseMQMQASave
    real(8) :: dAlphaSave, dGibbsSave
    integer, allocatable :: iAssemblageSave(:)
    real(8), allocatable :: dChemicalSave(:), dEffStoichSave(:,:), dElementSave(:), &
        dFractionSave(:), dMolesSave(:), dPhaseSave(:), dUpdateSave(:)


    ! Initialize the GEM solver:
    call InitGEMSolver
    call BeginMQMQATrajectoryAttempt(nElements,iterGlobalMax)

    !!!
    !!! CONSIDER MOVING THIS INTO THE InitGEMSolver SUBROUTINE:
    !!!
    ! The system may be converged if there aren't any solution phases:
    if ((nSolnPhases == 0).AND.(INFOThermo == 0).AND.(.NOT. lReinitLoaded)) then

        ! Check the system if only pure condensed phases are expected to appear:
        call CheckSysOnlyPureConPhases

        ! Report an error if this failed:
        if (.NOT.(lConverged)) then
            ! INFOThermo = 14
            ! return
        end if

    end if

    ! Begin the global iteration cycle:
    LOOP_GEMSolver: do iterGlobal = 1, iterGlobalMax

        call CaptureMQMQARecoveryPreStep(iterGlobal,iterLast,iterRevert,lRevertSystem, &
            iAssemblage(1:nElements),dGEMFunctionNorm)

        ! Ensures Newton and line-search paths know whether the current assemblage contains plain RKMP
        call UpdateRKMPHessianActivity

        ! If in debug mode, call the debugger:
        if (lDebugMode) call GEMDebug(1)

        ! Construct the Hessian matrix and compute the direction vector:
        call GEMNewton(INFO)
        call CaptureMQMQARecoveryNewton(iterGlobal,INFO,MAXVAL(DABS(dUpdateVar)), &
            dMQMQAHessianSelectedAlpha)

        ! Perform a line search using the direction vector:
        call GEMLineSearch

        ! Check if the estimated phase assemblage needs to be adjusted:
        call CheckPhaseAssemblage
        call CaptureMQMQARecoveryAssemblage(iterGlobal,iterLast,iterRevert,lRevertSystem)

        call CaptureMQMQATrajectoryPoint(iterGlobal,iterLast,iAssemblage(1:nElements), &
            dMolesPhase(1:nElements),dElementPotential(1:nElements), &
            nMQMQAHessianAcceptedSolveCount,nMQMQAHessianEligibleSolveCount, &
            nMQMQAHessianInteriorFallbackCount,dMQMQAHessianSelectedAlpha,dGEMFunctionNorm, &
            dMinGibbs,dMQMQAHessianMinimumRejectedFraction)

        ! Phase assemblage may have changed, so it needs to be refreshed before CheckConvergence and end-of-solve reporting.
        call UpdateRKMPHessianActivity

        ! Check convergence:
        ! if (iterGlobal /= iterLast) call CheckConvergence
        call CheckConvergence

        ! Opt-in fixed-point experiment: reconstruct the baseline and full
        ! MQMQA-corrected Newton systems at the exact live state that has just
        ! satisfied the production convergence test.  Every mutable solver
        ! array and runtime control is restored before normal exit.
        if (lConverged .AND. lCaptureMQMQAFixedPointAtConvergence .AND. &
            (.NOT. lMQMQAFixedPointCaptureAttempted)) then
            allocate(dChemicalSave(SIZE(dChemicalPotential)),dFractionSave(SIZE(dMolFraction)), &
                dMolesSave(SIZE(dMolesSpecies)),dPhaseSave(SIZE(dMolesPhase)), &
                dEffStoichSave(SIZE(dEffStoichSolnPhase,1),SIZE(dEffStoichSolnPhase,2)), &
                dElementSave(SIZE(dElementPotential)),dUpdateSave(SIZE(dUpdateVar)), &
                iAssemblageSave(SIZE(iAssemblage)))
            dChemicalSave = dChemicalPotential
            dFractionSave = dMolFraction
            dMolesSave = dMolesSpecies
            dPhaseSave = dMolesPhase
            dEffStoichSave = dEffStoichSolnPhase
            dElementSave = dElementPotential
            dUpdateSave = dUpdateVar
            iAssemblageSave = iAssemblage
            dGibbsSave = dGibbsEnergySys
            lRevertSave = lRevertSystem
            lUseMQMQASave = lUseMQMQAExactHessian
            lAdaptiveSave = lMQMQAHessianAdaptiveMode
            dAlphaSave = dMQMQAHessianAlpha

            call CaptureMQMQAFixedPointState(iAssemblage,dElementPotential,dChemicalPotential, &
                dMolesPhase,nElements,nSolnPhases,nConPhases)
            lUseMQMQAExactHessian = .TRUE.
            lMQMQAHessianAdaptiveMode = .FALSE.
            dMQMQAHessianAlpha = 1D0
            lCaptureGEMNewtonSystem = .TRUE.
            lCaptureGEMNewtonCorrectedSystem = .TRUE.
            lCaptureFirstGEMNewtonCorrectionPair = .TRUE.
            call GEMNewton(iCaptureInfo)
            iMQMQAFixedPointCaptureInfo = iCaptureInfo
            if (allocated(dCapturedMQMQAFixedPointSolvedUpdate)) &
                deallocate(dCapturedMQMQAFixedPointSolvedUpdate)
            allocate(dCapturedMQMQAFixedPointSolvedUpdate(SIZE(dUpdateVar)))
            dCapturedMQMQAFixedPointSolvedUpdate = dUpdateVar
            lMQMQAFixedPointCaptureAttempted = .TRUE.
            lCaptureGEMNewtonSystem = .FALSE.
            lCaptureGEMNewtonCorrectedSystem = .FALSE.
            lCaptureFirstGEMNewtonCorrectionPair = .FALSE.

            dMQMQAFixedPointRecomputeDifference = DMAX1(MAXVAL(ABS(dChemicalPotential-dChemicalSave)), &
                MAXVAL(ABS(dMolFraction-dFractionSave)),MAXVAL(ABS(dMolesSpecies-dMolesSave)), &
                MAXVAL(ABS(dMolesPhase-dPhaseSave)),MAXVAL(ABS(dEffStoichSolnPhase-dEffStoichSave)), &
                MAXVAL(ABS(dElementPotential-dElementSave)),ABS(dGibbsEnergySys-dGibbsSave))
            if (ANY(iAssemblage /= iAssemblageSave)) dMQMQAFixedPointRecomputeDifference = HUGE(1D0)

            dChemicalPotential = dChemicalSave
            dMolFraction = dFractionSave
            dMolesSpecies = dMolesSave
            dMolesPhase = dPhaseSave
            dEffStoichSolnPhase = dEffStoichSave
            dElementPotential = dElementSave
            dUpdateVar = dUpdateSave
            iAssemblage = iAssemblageSave
            dGibbsEnergySys = dGibbsSave
            lRevertSystem = lRevertSave
            lUseMQMQAExactHessian = lUseMQMQASave
            lMQMQAHessianAdaptiveMode = lAdaptiveSave
            dMQMQAHessianAlpha = dAlphaSave
            dMQMQAFixedPointRestorationError = DMAX1(MAXVAL(ABS(dChemicalPotential-dChemicalSave)), &
                MAXVAL(ABS(dMolFraction-dFractionSave)),MAXVAL(ABS(dMolesSpecies-dMolesSave)), &
                MAXVAL(ABS(dMolesPhase-dPhaseSave)),MAXVAL(ABS(dEffStoichSolnPhase-dEffStoichSave)), &
                MAXVAL(ABS(dElementPotential-dElementSave)),ABS(dGibbsEnergySys-dGibbsSave))
            if (ANY(iAssemblage /= iAssemblageSave)) dMQMQAFixedPointRestorationError = HUGE(1D0)
            deallocate(dChemicalSave,dFractionSave,dMolesSave,dPhaseSave,dEffStoichSave, &
                dElementSave,dUpdateSave,iAssemblageSave)
        end if

        ! If in debug mode, call the debugger:
        if (lDebugMode) call GEMDebug(9)

        ! Return control to the main program if an error has occured or if the solution has converged:
        if ((INFOThermo /= 0).OR.(lConverged)) exit LOOP_GEMSolver

    end do LOOP_GEMSolver

    ! Report an error if the GEMSolver did not converge but no other errors were encountered:
    if (.NOT.(lConverged).AND.(INFOThermo == 0)) INFOThermo = 12

    if (lDebugRKMPHessianFD .OR. &
        (lRKMPHessianReportSummary .AND. lUseRKMPExactHessian .AND. lRKMPHessianWasActive)) then
        ! Positional records keep the opt-in summary compact.  Fields follow the corresponding groups in
        ! ModuleGEMSolver: solve state, correction metrics/rejections, then nonlinear trust/direction metrics.
        write(*,*) 'RKMP_SOLVER_IMPACT', lUseRKMPExactHessian, dRKMPHessianBlendAlpha, &
            dRKMPHessianSelectedAlpha, dRKMPHessianUpdateNormRatio, iterGlobal, iterLast, iterRevert, &
            lConverged, lRevertSystem, INFOThermo, nSolnPhases, nConPhases, dGEMFunctionNorm
        write(*,*) 'RKMP_SOLVER_CORRECTION', dRKMPHessianMaxAppliedA, dRKMPHessianMaxAppliedB, &
            dRKMPHessianMaxAppliedRatio, dRKMPHessianMaxDeltaA, dRKMPHessianMaxDeltaB, &
            nRKMPHessianApplyCount, nRKMPHessianRejectDGESV, nRKMPHessianRejectBadDelta, &
            nRKMPHessianRejectRatio, nRKMPHessianRejectUpdate, nRKMPHessianRejectLocalResponse, &
            iRKMPHessianLastFailurePhase, iRKMPHessianLastFailureReason
        write(*,*) 'RKMP_NONLINEAR_TRUST', lRKMPHessianNonlinearReady, dRKMPHessianDirectionCosine, &
            dRKMPHessianDirectionDifference, nRKMPHessianRejectNonlinear, nRKMPHessianRejectDirection, &
            nRKMPHessianFullAlphaCount, dRKMPHessianMaxSelectedAlpha
    end if

    return

contains

    !> \brief Refresh whether the current assemblage contains a supported RKMP phase.
    !!
    !> \details Exact-curvature solver behavior is meaningful only while a plain RKMP phase is active.  The
    !! historical flag preserves end-of-solve diagnostics when that phase disappears before convergence.
    subroutine UpdateRKMPHessianActivity

        integer :: iPhase, kPhase

        lRKMPHessianActive = .FALSE.
        do iPhase = 1, nSolnPhases
            kPhase = -iAssemblage(nElements - iPhase + 1)
            if (kPhase <= 0) cycle
            if (cSolnPhaseType(kPhase) == 'RKMP') then
                lRKMPHessianActive = .TRUE.
                exit
            end if
        end do
        lRKMPHessianWasActive = lRKMPHessianWasActive .OR. lRKMPHessianActive

    end subroutine UpdateRKMPHessianActivity

end subroutine GEMSolver
