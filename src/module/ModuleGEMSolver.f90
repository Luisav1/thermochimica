
!-------------------------------------------------------------------------------------------------------------
    !
    !> \file        ModuleGEMSolver.f90
    !> \brief       Fortran module for input/output of the non-linear solver.
    !> \author      M.H.A. Piro
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !
    !> \param iterLast          The last global iteration that the phase assemblage was adjusted.
    !> \param iterHistory       An integer matrix representing all of the indices of phases that contribute
    !!                           to the equilibrium phase assemblage at each stage in the iteration history.
    !> \param dGEMFunctionNorm  A double real scalar representing the norm of the functional vector in the
    !!                          GEMSolver.
    !> \param dGEMFunctionNormLast A double real scalar representing the norm of the functional vector in the
    !!                              GEMSolver from the last iteration.
    !> \param dSumMolFractionSoln  A double real vector representing the sum of mole fractions in each solution
    !!                              phase.
    !> \param dUpdateVar        A double real vector representing the direction vector that updates the function
    !                            vector.
    !> \param dPartialExcessGibbs  A double real vector representing the partial molar excess Gibbs energy of
    !!                              mixing of each species in the system.
    !> \param dEffStoichSolnPhase   A double real matrix representing the effective stoichiometry of each
    !!                               solution phase.
    !> \param nRKMPHessianApplyCount Number of accepted RKMP response corrections applied during the current
    !!                               GEM solve.
    !> \param nRKMPHessianRejectDGESV Number of RKMP alpha-trust trial systems rejected because DGESV failed.
    !> \param nRKMPHessianRejectBadDelta Number of RKMP alpha-trust trial systems rejected because the mapped
    !!                                   correction contained invalid floating-point values.
    !> \param nRKMPHessianRejectRatio Number of RKMP alpha-trust trial systems rejected because the mapped
    !!                                correction was too large relative to the existing element block.
    !> \param nRKMPHessianRejectUpdate Number of RKMP alpha-trust trial systems rejected because the solved
    !!                                 update was too large relative to the alpha-zero update.
    !> \param nRKMPHessianRejectNonlinear Number of GEM iterations kept on the ideal direction because the
    !!                                    current state had not yet established local nonlinear trust.
    !> \param nRKMPHessianRejectDirection Number of RKMP alpha trials rejected because their solved direction
    !!                                    was insufficiently aligned with the alpha-zero direction.
    !> \param nRKMPHessianRejectLocalResponse Number of alpha trials rejected because a phase-local constrained
    !!                                        response or mapped residual could not be constructed.
    !> \param iRKMPHessianLastFailurePhase Absolute solution-phase index associated with the most recent mapper
    !!                                     failure, or zero when no failure has occurred.
    !> \param iRKMPHessianLastFailureReason Machine-readable RKMP mapper failure code.  Public constants below
    !!                                      distinguish Hessian, response, ideal-baseline, and finite-value errors.
    !> \param nRKMPHessianFullAlphaCount Number of GEM iterations that accepted the undamped RKMP correction.
    !> \param dRKMPHessianBlendAlpha Requested upper blend for the mapped RKMP GEM correction.  Trust logic
    !!                                selects the effective blend; neither value scales the RKMP derivatives or
    !!                                local excess Hessian.
    !> \param dRKMPHessianMaxAppliedA Maximum accepted alpha-scaled RKMP correction applied to the GEM A matrix.
    !> \param dRKMPHessianMaxAppliedB Maximum accepted alpha-scaled RKMP correction applied to the GEM B vector.
    !> \param dRKMPHessianMaxAppliedRatio Maximum accepted A correction relative to the current element block.
    !> \param dRKMPHessianMaxDeltaA Maximum unscaled RKMP correction candidate for the GEM A matrix.
    !> \param dRKMPHessianMaxDeltaB Maximum unscaled RKMP correction candidate for the GEM B vector.
    !> \param dRKMPHessianSelectedAlpha Largest RKMP alpha accepted by the Stage 1E alpha-trust trial solve.
    !> \param dRKMPHessianUpdateNormRatio Accepted update norm divided by the alpha-zero update norm.
    !> \param dRKMPHessianDirectionCosine Cosine between the accepted RKMP and alpha-zero update directions.
    !> \param dRKMPHessianDirectionDifference Relative two-norm difference between those update directions.
    !> \param lRKMPHessianNonlinearReady True after the current state is feasible, near the best Gibbs state,
    !!                                   locally settled, and making acceptable residual progress.
    !> \param lRKMPHessianActive True when the current assemblage contains a plain RKMP solution phase.
    !> \param lRKMPHessianWasActive True after a plain RKMP phase has appeared during the current GEM solve.
    !> \param lUseMQMQAExactHessian True only when the caller requests the default-off fixed-alpha MQMQA path.
    !> \param dMQMQAHessianAlpha Fixed weight applied to each completed aggregate `deltaA/deltaB` correction;
    !!                            it never scales the local SUBG or SUBQ Hessian.
    !> \param nMQMQAHessianApplyCount Number of valid MQMQA correction aggregates applied to trial systems.
    !> \param nMQMQAHessianAcceptedSolveCount Number of corrected trial systems accepted after a finite solve.
    !> \param nMQMQAHessianPhaseCorrectionCount Number of successful phase-local pairs included in aggregates.
    !> \param nMQMQAHessianChargedSkipCount Number of deliberately excluded charged MQMQA phases.
    !> \param nMQMQAHessianInteriorFallbackCount Number of routed states rejected only because at least one
    !!                                             local fraction is at or below the strict interior threshold.
    !> \param dMQMQAHessianMinimumRejectedFraction Smallest local fraction observed among those boundary states.
    !> \param dMQMQAHessianMaxRatioA Largest Frobenius norm of applied `alpha*deltaA` relative to the baseline
    !!                                GEM element block during the current calculation.
    !> \param dMQMQAHessianMaxRatioB Largest two-norm of applied `alpha*deltaB` relative to the baseline element
    !!                                residual during the current calculation.
    !> \param lMQMQAHessianFallbackUsed True after any routed build or corrected linear trial returns to the
    !!                                   untouched historical GEM system.
    !> \param lDebugMode        A logical variable used for debugging purposes.  When it is TRUE, a number
    !!                           of print statements are applied.
    !> \param lRevertSystem     A logical variable identifying whether the system should be reverted (TRUE)
    !!                           or not (FALSE).
    !> \param lConverged        A logical variable identifying whether the system has converged (TRUE) or not
    !!                           (FALSE).
    !> \param lSolnPhases       A logical vector indicating whether a particular solution phase is currently
    !!                           assumed to be stable (true) or not (false).
    !> \param lMiscibility      A logical vector indicating whether a particular solution phase has a
    !!                           miscibility gap (true) or not (false).
    !
    !
    ! CONSIDER REMOVING THE FOLLOWING VARIABLES:
    !
    ! iterLastCon           The last iteration that a pure condensed phase was either added to or removed from
    !                        the estimated phase assemblage.
    ! iterLastSoln          The last iteration that a pure condensed phase was either added to or removed from
    !                        the estimated phase assemblage.
    ! iConPhaseLast         The species index of the last pure condensed phase that was either added to, or
    !                        removed from, the estimated phase assemblage.
    ! iSolnPhaseLast        The species index of the last pure condensed phase that was either added to, or
    !                        removed from, the estimated phase assemblage.
    !
    !
!-------------------------------------------------------------------------------------------------------------

module ModuleGEMSolver

    implicit none

    SAVE

    integer                              ::  iterLast,      iterStep, iterRevert, iterGlobal
    integer                              ::  iterLastCon,   iterLastSoln,         iterSwap,   iterLastMiscGapCheck
    integer                              ::  nRKMPHessianApplyCount
    integer                              ::  nRKMPHessianRejectDGESV, nRKMPHessianRejectBadDelta
    integer                              ::  nRKMPHessianRejectRatio, nRKMPHessianRejectUpdate
    integer                              ::  nRKMPHessianRejectNonlinear, nRKMPHessianRejectDirection
    integer                              ::  nRKMPHessianRejectLocalResponse
    integer                              ::  nRKMPHessianFullAlphaCount
    integer                              ::  iRKMPHessianLastFailurePhase, iRKMPHessianLastFailureReason
    integer                              ::  nMQMQAHessianApplyCount, nMQMQAHessianAcceptedSolveCount
    integer                              ::  nMQMQAHessianPhaseCorrectionCount, nMQMQAHessianChargedSkipCount
    integer                              ::  nMQMQAHessianAggregateFailureCount, nMQMQAHessianApplicationFailureCount
    integer                              ::  nMQMQAHessianDGESVFallbackCount, nMQMQAHessianNonfiniteFallbackCount
    integer                              ::  nMQMQAHessianRKMPConflictCount, nMQMQAHessianInteriorFallbackCount
    integer                              ::  iMQMQAHessianLastFailurePhase, iMQMQAHessianLastFailureStatus
    integer, parameter                   ::  RKMP_MAP_SUCCESS = 0
    integer, parameter                   ::  RKMP_MAP_HESSIAN_FAILURE = 1
    integer, parameter                   ::  RKMP_MAP_ELEMENT_RESPONSE_FAILURE = 2
    integer, parameter                   ::  RKMP_MAP_RESIDUAL_RESPONSE_FAILURE = 3
    integer, parameter                   ::  RKMP_MAP_IDEAL_RESPONSE_FAILURE = 4
    integer, parameter                   ::  RKMP_MAP_INVALID_CORRECTION = 5
    integer, parameter                   ::  MQMQA_INTEGRATION_SUCCESS = 0
    integer, parameter                   ::  MQMQA_INTEGRATION_NO_APPLICABLE_PHASE = 1
    integer, parameter                   ::  MQMQA_INTEGRATION_AGGREGATE_FAILURE = 2
    integer, parameter                   ::  MQMQA_INTEGRATION_APPLICATION_FAILURE = 3
    integer, parameter                   ::  MQMQA_INTEGRATION_DGESV_FAILURE = 4
    integer, parameter                   ::  MQMQA_INTEGRATION_NONFINITE_UPDATE = 5
    integer, parameter                   ::  MQMQA_INTEGRATION_RKMP_CONFLICT = 6
    integer                              ::  iConPhaseLast, iSolnPhaseLast,       iSolnSwap,  iPureConSwap
    integer,                 parameter   ::  iterGlobalMax = 3000
    integer, dimension(:,:), allocatable ::  iterHistory

    real(8)                              ::  dGEMFunctionNorm,    dGEMFunctionNormLast, dMaxSpeciesChange, dMinGibbs
    real(8)                              ::  dRKMPHessianBlendAlpha
    real(8)                              ::  dRKMPHessianMaxAppliedA, dRKMPHessianMaxAppliedB
    real(8)                              ::  dRKMPHessianMaxAppliedRatio, dRKMPHessianMaxDeltaA
    real(8)                              ::  dRKMPHessianMaxDeltaB
    real(8)                              ::  dRKMPHessianSelectedAlpha, dRKMPHessianUpdateNormRatio
    real(8)                              ::  dRKMPHessianDirectionCosine, dRKMPHessianDirectionDifference
    real(8)                              ::  dRKMPHessianMaxSelectedAlpha
    real(8)                              ::  dMQMQAHessianAlpha
    real(8)                              ::  dMQMQAHessianMaxDeltaA, dMQMQAHessianMaxDeltaB
    real(8)                              ::  dMQMQAHessianMaxAppliedA, dMQMQAHessianMaxAppliedB
    real(8)                              ::  dMQMQAHessianMaxRatioA, dMQMQAHessianMaxRatioB
    real(8)                              ::  dMQMQAHessianMinimumRejectedFraction
    real(8), dimension(iterGlobalMax)     ::  dRKMPHessianAcceptedAlphaHistory
    real(8)                              ::  dRKMPTrustEmergencyRatioCap = 1D6
    real(8)                              ::  dRKMPTrustUpdateRatioCap = 1.25D0
    real(8)                              ::  dRKMPTrustDirectionCosineMin = 0.90D0
    real(8)                              ::  dRKMPTrustDirectionDifferenceCap = 0.50D0
    real(8)                              ::  dRKMPTrustLocalNormThreshold = 5D-2
    real(8)                              ::  dRKMPTrustProgressAllowance = 1.05D0
    real(8)                              ::  dRKMPTrustGibbsActivationTolerance = 1D-6
    real(8)                              ::  dRKMPTrustGibbsRetentionTolerance = 1D-4
    real(8), dimension(:),   allocatable ::  dSumMolFractionSoln, dMolesPhaseLast, dUpdateVar, dDrivingForceSoln
    real(8), dimension(:),   allocatable ::  dPartialExcessGibbs, dPartialExcessGibbsLast
    real(8), dimension(:,:), allocatable ::  dEffStoichSolnPhase

    logical                              ::  lDebugMode, lRevertSystem, lConverged
    logical                              ::  lUseRKMPExactHessian, lDebugRKMPHessianFD
    logical                              ::  lRKMPHessianNonlinearReady, lRKMPHessianActive, lRKMPHessianWasActive
    logical                              ::  lUseMQMQAExactHessian, lMQMQAHessianSupportedPhaseFound
    logical                              ::  lMQMQAHessianEligibleCorrectionBuilt, lMQMQAHessianAggregateBuilt
    logical                              ::  lMQMQAHessianCorrectionApplied, lMQMQAHessianCorrectedSolveAccepted
    logical                              ::  lMQMQAHessianFallbackUsed
    logical                              ::  lRKMPHessianControlsConfigured = .FALSE.
    logical                              ::  lRKMPHessianRequestedEnable = .FALSE.
    logical                              ::  lRKMPHessianRequestedDebug = .FALSE.
    logical                              ::  lRKMPHessianReportSummary = .FALSE.
    real(8)                              ::  dRKMPHessianRequestedAlphaMax = 0.10D0
    logical                              ::  lMQMQAHessianControlsConfigured = .FALSE.
    logical                              ::  lMQMQAHessianRequestedEnable = .FALSE.
    real(8)                              ::  dMQMQAHessianRequestedAlpha = 0D0
    logical, dimension(:),   allocatable ::  lSolnPhases, lMiscibility

end module ModuleGEMSolver
