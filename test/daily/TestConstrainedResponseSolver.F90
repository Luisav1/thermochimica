!-------------------------------------------------------------------------------------------------------------
!> \file    TestConstrainedResponseSolver.F90
!> \brief   Unit verification of the model-independent constrained-response solve.
!>
!> \details Exercises multiple forcing columns, normalization, agreement with
!!          a null-space elimination, agreement with the established RKMP
!!          mixed-sign convention, and explicit singular/input failures. No
!!          Thermochimica phase state or GEM matrix is modified.
!-------------------------------------------------------------------------------------------------------------

program TestConstrainedResponseSolver

    USE ModuleConstrainedResponse, ONLY: SolveConstrainedResponse
    USE ModuleRKMPResponseMapping, ONLY: SolveRKMPConstrainedResponse

    implicit none

    integer :: iInfo, iInfoRKMP
    logical :: lPass
    real(8) :: dConstraintResidual, dKKTResidual, dNullDifference, dNullForcingDifference
    real(8) :: dRKMPDifference
    real(8) :: dConstraint(1,4), dCurvature(4,4), dForcing(4,2), dMultiplier(1,2)
    real(8) :: dResponse(4,2), dResponseNull(4,2), dResponseRKMP(4,2), dZ(4,3)
    real(8) :: dNullForcing(4,1), dNullForcingResponse(4,1)
    real(8) :: dReduced(3,3), dReducedRHS(3,2)
    real(8) :: dRedundant(2,4), dSingular(4,4), dAsymmetric(4,4)

    lPass = .TRUE.
    dCurvature = 0D0
    dCurvature(1,1) = 2D0
    dCurvature(2,2) = 3D0
    dCurvature(3,3) = 5D0
    dCurvature(4,4) = 7D0
    dConstraint = 1D0
    dForcing(:,1) = [1D0,-2D0,3D0,0.5D0]
    dForcing(:,2) = [-1D0,4D0,0.25D0,2D0]

    call SolveConstrainedResponse(dCurvature,dConstraint,dForcing,dResponse,iInfo,dMultiplier)
    lPass = lPass .AND. (iInfo == 0)
    dConstraintResidual = MAXVAL(DABS(MATMUL(dConstraint,dResponse)))
    dKKTResidual = MAXVAL(DABS(MATMUL(dCurvature,dResponse)+ &
        MATMUL(TRANSPOSE(dConstraint),dMultiplier)-dForcing))
    lPass = lPass .AND. (dConstraintResidual <= 1D-13) .AND. (dKKTResidual <= 1D-12)

    ! Columns of Z span all total-preserving composition changes. Solving in
    ! this smaller coordinate space must recover the bordered-system response.
    dZ = 0D0
    dZ(1,1) = 1D0
    dZ(2,2) = 1D0
    dZ(3,3) = 1D0
    dZ(4,:) = -1D0
    dReduced = MATMUL(TRANSPOSE(dZ),MATMUL(dCurvature,dZ))
    dReducedRHS = MATMUL(TRANSPOSE(dZ),dForcing)
    call SolveDense(dReduced,dReducedRHS,iInfo)
    dResponseNull = MATMUL(dZ,dReducedRHS)
    dNullDifference = MAXVAL(DABS(dResponse-dResponseNull))
    lPass = lPass .AND. (iInfo == 0) .AND. (dNullDifference <= 1D-12)

    call SolveRKMPConstrainedResponse(4,2,dCurvature,dForcing,dResponseRKMP,iInfoRKMP)
    dRKMPDifference = MAXVAL(DABS(dResponse-dResponseRKMP))
    lPass = lPass .AND. (iInfoRKMP == 0) .AND. (dRKMPDifference <= 1D-12)

    ! A forcing in the constraint-normal direction changes only the Lagrange
    ! multiplier. It cannot drive a normalized composition response.
    dNullForcing(:,1) = 2.5D0
    call SolveConstrainedResponse(dCurvature,dConstraint,dNullForcing,dNullForcingResponse,iInfo)
    dNullForcingDifference = MAXVAL(DABS(dNullForcingResponse))
    lPass = lPass .AND. (iInfo == 0) .AND. (dNullForcingDifference <= 1D-13)

    ! Redundant constraints and unconstrained tangent curvature must fail as
    ! singular systems instead of returning an apparently usable response.
    dRedundant(1,:) = 1D0
    dRedundant(2,:) = 2D0
    call SolveConstrainedResponse(dCurvature,dRedundant,dForcing,dResponse,iInfo)
    lPass = lPass .AND. (iInfo > 0)
    dSingular = 0D0
    call SolveConstrainedResponse(dSingular,dConstraint,dForcing,dResponse,iInfo)
    lPass = lPass .AND. (iInfo > 0)

    ! Strong scale separation is allowed when the solve remains finite and its
    ! physical and algebraic residuals remain bounded.
    dCurvature = 0D0
    dCurvature(1,1) = 1D-10
    dCurvature(2,2) = 1D0
    dCurvature(3,3) = 1D5
    dCurvature(4,4) = 1D10
    call SolveConstrainedResponse(dCurvature,dConstraint,dForcing,dResponse,iInfo,dMultiplier)
    dConstraintResidual = MAXVAL(DABS(MATMUL(dConstraint,dResponse)))
    dKKTResidual = MAXVAL(DABS(MATMUL(dCurvature,dResponse)+ &
        MATMUL(TRANSPOSE(dConstraint),dMultiplier)-dForcing)) / &
        DMAX1(1D0,MAXVAL(DABS(dForcing)))
    lPass = lPass .AND. (iInfo == 0) .AND. (dConstraintResidual <= 1D-8) .AND. &
        (dKKTResidual <= 1D-8)

    dAsymmetric = dCurvature
    dAsymmetric(1,2) = 1D0
    call SolveConstrainedResponse(dAsymmetric,dConstraint,dForcing,dResponse,iInfo)
    lPass = lPass .AND. (iInfo == -5)

    if (lPass) then
        print *, 'TestConstrainedResponseSolver: PASS'
        call EXIT(0)
    else
        print *, 'TestConstrainedResponseSolver: FAIL <---'
        print *, 'constraint residual = ',dConstraintResidual
        print *, 'KKT residual = ',dKKTResidual
        print *, 'null-space difference = ',dNullDifference
        print *, 'RKMP-sign difference = ',dRKMPDifference
        print *, 'null-forcing response = ',dNullForcingDifference
        call EXIT(1)
    end if

contains

    subroutine SolveDense(dMatrix,dRHSLocal,iInfoLocal)

        real(8), intent(in) :: dMatrix(:,:)
        real(8), intent(inout) :: dRHSLocal(:,:)
        integer, intent(out) :: iInfoLocal

        integer :: n
        integer, allocatable :: iPivot(:)
        real(8), allocatable :: dWork(:,:)

        n = SIZE(dMatrix,1)
        allocate(dWork(n,n),iPivot(n))
        dWork = dMatrix
        call DGESV(n,SIZE(dRHSLocal,2),dWork,n,iPivot,dRHSLocal,n,iInfoLocal)
        deallocate(dWork,iPivot)

    end subroutine SolveDense

end program TestConstrainedResponseSolver
