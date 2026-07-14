
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
    !> \param nRKMPHessianFullAlphaCount Number of GEM iterations that accepted the undamped RKMP correction.
    !> \param dRKMPHessianBlendAlpha User-facing nominal blend factor for the experimental RKMP correction.
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
    integer                              ::  nRKMPHessianFullAlphaCount
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
    real(8), dimension(:),   allocatable ::  dSumMolFractionSoln, dMolesPhaseLast, dUpdateVar, dDrivingForceSoln
    real(8), dimension(:),   allocatable ::  dPartialExcessGibbs, dPartialExcessGibbsLast
    real(8), dimension(:,:), allocatable ::  dEffStoichSolnPhase

    logical                              ::  lDebugMode, lRevertSystem, lConverged
    logical                              ::  lUseRKMPExactHessian, lDebugRKMPHessianFD
    logical                              ::  lRKMPHessianNonlinearReady
    logical, dimension(:),   allocatable ::  lSolnPhases, lMiscibility

end module ModuleGEMSolver
