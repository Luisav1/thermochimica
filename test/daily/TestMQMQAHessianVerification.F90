!-------------------------------------------------------------------------------------------------------------
!> \file    TestMQMQAHessianVerification.F90
!> \brief   Standalone verification of SUBG and staged nonmagnetic SUBQ curvature.
!>
!> \details Independent ordinary-real and second-order-object evaluators are compared before finite
!!          differences test the analytic gradient and Hessian. The cases isolate ordinary
!!          configurational mixing, the G and Q binary parameter families, all three
!!          supported ternary-group branches, and the B weighted-pair family. Pass
!!          --report for the complete numerical evidence.
!!
!!          Verification map:
!!          1. Check the generic second-order calculus kernel.
!!          2. Build one complete positive synthetic topology and exercise SUBQ S3, chi, and pair zeta changes.
!!          3. Define isolated G, Q, ternary, and B interaction cases.
!!          4. Apply the scalar/gradient/Hessian verification ladder to each case.
!!          5. Quantify the known paper-versus-production SUBQ S3 distinction.
!!          6. Measure that distinction relative to complete controlled SUBQ
!!             energy, chemical-potential, and Hessian scales.
!!          7. Check the independently derived extensive B identity.
!!          8. Confirm unsupported and boundary inputs fail explicitly.
!-------------------------------------------------------------------------------------------------------------

program TestMQMQAHessianVerification

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleMQMQAUnconstrained
    USE ModuleFiniteDifferenceVerification

    implicit none

    type(MQMQAModelData) :: tModel, tSUBQModel
    type(MQMQAInteractionTerm), allocatable :: tAll(:), tOne(:), tNone(:), tSwapped(:)
    real(8), allocatable :: dMoles(:)
    logical :: lPass, lReport
    integer :: iKernelInfo
    character(len=32) :: cArgument

    lPass=.TRUE.
    lReport=.FALSE.
    if (COMMAND_ARGUMENT_COUNT()>0) then
        call GET_COMMAND_ARGUMENT(1,cArgument)
        lReport=TRIM(cArgument)=='--report'
    end if

    ! First verify the generic calculus engine without any MQMQA equations. A
    ! failure here means the chain-rule machinery is wrong, not the thermodynamics.
    call CheckMQMQADerivativeKernel(iKernelInfo)
    lPass=lPass.AND.(iKernelInfo==0)
    if (lReport) write(*,'(A,I0)') 'derivative-kernel iInfo = ',iKernelInfo

    call BuildModel(tModel,dMoles)
    tSUBQModel=tModel
    tSUBQModel%iModelType=MQMQA_MODEL_SUBQ
    ! SUBG databases provide one common zeta. This deliberately nonuniform
    ! matrix exercises the pair-specific SUBQ path and prevents a global-zeta
    ! implementation from passing by coincidence.
    tSUBQModel%dZeta=RESHAPE([3.8D0,4.0D0,4.2D0,4.4D0,4.1D0,4.3D0,4.5D0,4.7D0],[4,2])
    call BuildInteractions(tAll)
    allocate(tNone(0),tOne(1),tSwapped(1))
    call SetGQTerm(tSwapped(1),MQMQA_TERM_G,1,1,1,2,1,1,0,0, &
        [ .TRUE.,.FALSE. ],[ .FALSE.,.TRUE. ],48D0)

    if (lReport) then
        write(*,'(A)') 'Stage MQ-1 SUBG and staged SUBQ S3/chi/zeta Hessian verification'
        write(*,'(A,I0)') 'storage bits = ',STORAGE_SIZE(1D0)
        write(*,'(A,I0)') 'decimal precision = ',PRECISION(1D0)
        write(*,'(A,I0)') 'binary digits = ',DIGITS(1D0)
        write(*,'(A,ES14.6)') 'machine epsilon = ',EPSILON(1D0)
    end if
    lPass=lPass.AND.(STORAGE_SIZE(1D0)==64).AND.(PRECISION(1D0)>=15).AND.(DIGITS(1D0)>=53)

    ! Add one physical/mathematical layer at a time. This localizes a failure to
    ! configurational mixing, one production family, or one ternary branch.
    call VerifyCase('reference + ideal',tModel,dMoles,1D0,tNone,lPass,lReport)
    call VerifyCase('SUBQ nonuniform-zeta configurational',tSUBQModel,dMoles,1D0,tNone,lPass,lReport)
    tOne(1)=tAll(1); call VerifyCase('G binary',tModel,dMoles,1D0,tOne,lPass,lReport)
    tOne(1)=tAll(1); call VerifyCase('SUBQ G binary chi incidence',tSUBQModel,dMoles,1D0,tOne,lPass,lReport)
    call VerifyCase('SUBQ G binary swapped chi incidence',tSUBQModel,dMoles,1D0,tSwapped,lPass,lReport)
    tOne(1)=tAll(2); call VerifyCase('Q binary',tModel,dMoles,1D0,tOne,lPass,lReport)
    tOne(1)=tAll(3); call VerifyCase('ternary group 1',tModel,dMoles,1D0,tOne,lPass,lReport)
    tOne(1)=tAll(4); call VerifyCase('ternary group 2',tModel,dMoles,1D0,tOne,lPass,lReport)
    tOne(1)=tAll(5); call VerifyCase('ternary neither group',tModel,dMoles,1D0,tOne,lPass,lReport)
    tOne(1)=tAll(6); call VerifyCase('B family',tModel,dMoles,1D0,tOne,lPass,lReport)
    tOne(1)=tAll(6); call VerifyCase('SUBQ nonuniform-zeta configurational + B', &
        tSUBQModel,dMoles,1D0,tOne,lPass,lReport)
    call VerifyCase('integrated total',tModel,dMoles,1D0,tAll,lPass,lReport)
    call VerifyS3FormulationDifference(tSUBQModel,dMoles,lPass,lReport)
    call VerifyS3TotalSignificance(tSUBQModel,dMoles,tAll,lPass,lReport)
    call VerifyPairSpecificZeta(tSUBQModel,dMoles,tNone,tAll(6),lPass,lReport)
    call VerifyBProductionIdentity('SUBG common-zeta',tModel,dMoles,tAll(6),lPass,lReport)
    call VerifyBProductionIdentity('SUBQ pair-specific-zeta',tSUBQModel,dMoles,tAll(6),lPass,lReport)
    call VerifyFailures(tModel,dMoles,tAll(1),lPass,lReport)

    if (lPass) then
        print *, 'TestMQMQAHessianVerification: PASS'
        call EXIT(0)
    else
        print *, 'TestMQMQAHessianVerification: FAIL <---'
        call EXIT(1)
    end if

contains

    !=========================================================================================================
    ! SECTION 1: CONTROLLED MODEL AND INTERACTION FIXTURES
    !
    ! These routines construct a complete, nonuniform interior state and isolate
    ! each supported interaction branch. They provide controlled mathematical
    ! coverage; they do not decode a production Thermochimica SUBG phase.
    !=========================================================================================================

    !---------------------------------------------------------------------------------------------------------
    !> \brief Build a positive interior state containing every canonical quadruplet.
    !>
    !> \details The synthetic topology is intentionally complete rather than a toy
    !!          binary subset. Nonuniform moles, coordination numbers, and reference
    !!          energies prevent accidental cancellation from making an incorrect
    !!          derivative appear correct. This base SUBG model uses one common zeta;
    !!          the caller installs nonuniform pair-specific values for SUBQ.
    !---------------------------------------------------------------------------------------------------------
    subroutine BuildModel(tData,dState)

        type(MQMQAModelData), intent(out) :: tData
        real(8), allocatable, intent(out) :: dState(:)
        integer :: a,b,x,y,q,nQuad

        tData%nSublattice1=4
        tData%nSublattice2=2
        nQuad=(4*5/2)*(2*3/2)
        allocate(tData%iQuadruplet(nQuad,4),tData%dCoordination(nQuad,4), &
            tData%dZeta(4,2),tData%dReferenceEnergy(nQuad),dState(nQuad))
        q=0
        do x=1,2
            do y=x,2
                do a=1,4
                    do b=a,4
                        q=q+1
                        tData%iQuadruplet(q,:)=[a,b,x,y]
                        tData%dCoordination(q,:)=[4D0+0.05D0*a,4D0+0.05D0*b, &
                            4.2D0+0.04D0*x,4.2D0+0.04D0*y]
                        tData%dReferenceEnergy(q)=-2.5D0+0.11D0*q
                        dState(q)=0.35D0+0.017D0*q+0.003D0*MOD(q,4)
                    end do
                end do
            end do
        end do
        tData%dZeta=4.2D0

    end subroutine BuildModel


    !---------------------------------------------------------------------------------------------------------
    !> \brief Define isolated terms covering G, Q, every ternary branch, and B.
    !---------------------------------------------------------------------------------------------------------
    subroutine BuildInteractions(tTerm)

        type(MQMQAInteractionTerm), allocatable, intent(out) :: tTerm(:)

        allocate(tTerm(6))
        call SetGQTerm(tTerm(1),MQMQA_TERM_G,1,2,1,1,1,1,0,0, &
            [ .TRUE.,.FALSE.,.TRUE.,.FALSE. ],[ .FALSE.,.TRUE.,.FALSE.,.TRUE. ],120D0)
        call SetGQTerm(tTerm(2),MQMQA_TERM_Q,1,2,1,1,2,1,0,0, &
            [ .TRUE.,.FALSE.,.TRUE.,.FALSE. ],[ .FALSE.,.TRUE.,.FALSE.,.TRUE. ],-75D0)
        call SetGQTerm(tTerm(3),MQMQA_TERM_G,1,2,1,1,1,1,2,3, &
            [ .TRUE.,.FALSE.,.TRUE.,.FALSE. ],[ .FALSE.,.TRUE.,.FALSE.,.FALSE. ],32D0)
        call SetGQTerm(tTerm(4),MQMQA_TERM_Q,1,2,1,1,1,2,2,3, &
            [ .TRUE.,.FALSE.,.FALSE.,.FALSE. ],[ .FALSE.,.TRUE.,.TRUE.,.FALSE. ],-28D0)
        call SetGQTerm(tTerm(5),MQMQA_TERM_G,1,2,1,1,1,1,2,3, &
            [ .TRUE.,.FALSE.,.FALSE.,.FALSE. ],[ .FALSE.,.TRUE.,.FALSE.,.FALSE. ],18D0)
        tTerm(6)%iFamily=MQMQA_TERM_B
        tTerm(6)%iA=1; tTerm(6)%iB=2; tTerm(6)%iX=1; tTerm(6)%iY=2
        tTerm(6)%iExponentP=1; tTerm(6)%iExponentQ=0
        tTerm(6)%dCoefficient=55D0

    end subroutine BuildInteractions


    subroutine SetGQTerm(tTerm,iFamily,a,b,x,y,p,q,r,d,lGroup1,lGroup2,dCoefficient)

        type(MQMQAInteractionTerm), intent(out) :: tTerm
        integer, intent(in) :: iFamily,a,b,x,y,p,q,r,d
        logical, intent(in) :: lGroup1(:),lGroup2(:)
        real(8), intent(in) :: dCoefficient

        tTerm%iFamily=iFamily
        tTerm%iA=a; tTerm%iB=b; tTerm%iX=x; tTerm%iY=y
        tTerm%iExponentP=p; tTerm%iExponentQ=q; tTerm%iExponentR=r
        tTerm%iTernaryConstituent=d
        tTerm%dCoefficient=dCoefficient
        allocate(tTerm%lGroup1(SIZE(lGroup1)),tTerm%lGroup2(SIZE(lGroup2)))
        tTerm%lGroup1=lGroup1
        tTerm%lGroup2=lGroup2

    end subroutine SetGQTerm

    !=========================================================================================================
    ! SECTION 2: ENERGY, GRADIENT, AND HESSIAN VERIFICATION LADDER
    !
    ! Each case first compares the two independent energy paths, then checks the
    ! derivative path using finite differences and extensive-energy identities.
    !=========================================================================================================


    !---------------------------------------------------------------------------------------------------------
    !> \brief Apply the complete verification ladder to one selected energy model.
    !>
    !> \details The independent scalar and derivative paths must first agree in
    !!          value. Finite differences then check the gradient. Multiplying the
    !!          Hessian by a chosen mole-perturbation direction predicts the
    !!          corresponding change in every gradient component; directional
    !!          curvature predicts the scalar energy bending along that direction.
    !!          Symmetry, homogeneity, and
    !!          extensivity test structural properties that finite differences alone
    !!          would not establish.
    !---------------------------------------------------------------------------------------------------------
    subroutine VerifyCase(cName,tData,dState,dIdealScale,tTerm,lAllPass,lVerbose)

        character(*), intent(in) :: cName
        type(MQMQAModelData), intent(in) :: tData
        real(8), intent(in) :: dState(:),dIdealScale
        type(MQMQAInteractionTerm), intent(in) :: tTerm(:)
        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lVerbose

        integer, parameter :: nSteps=9
        integer :: i,j,k,n,iInfo
        real(8) :: dG,dGRef,dGIdeal,dGEx,dGDual,dGDualRef,dGDualIdeal,dGDualEx
        real(8) :: dScaleError,dSym,dHom,dNormH,dNormN,dHBase,dH,dGradientHBase
        real(8) :: dAnalytic,dFD3,dFD5,dErr3(nSteps),dErr5(nSteps)
        real(8) :: dAbsErr3(nSteps),dAbsErr5(nSteps),dStep(nSteps)
        real(8) :: dGradientError,dHvError,dExtError,dGPlus,dGMinus,dTmp1,dTmp2,dTmp3
        real(8) :: dGradNormAbs(nSteps),dGradNormScaled(nSteps),dGradMaxAbs(nSteps),dGradMaxScaled(nSteps)
        real(8) :: dHvNormAbs(nSteps),dHvNormScaled(nSteps),dHvMaxAbs(nSteps),dHvMaxScaled(nSteps)
        real(8) :: dOrder3(nSteps),dOrder5(nSteps),dOrderGradient(nSteps),dOrderHv(nSteps)
        integer :: iGradWorst(nSteps),iHvWorst(nSteps)
        real(8), allocatable :: dHessian(:,:),dHRef(:,:),dHIdeal(:,:),dHEx(:,:),dGradient(:),dDirection(:)
        real(8), allocatable :: dPlus(:),dMinus(:),dPlus2(:),dMinus2(:),dGradPlus(:),dGradMinus(:),dDummyH(:,:)
        real(8), allocatable :: dHv(:),dGradientFD(:),dHvFD(:)
        logical :: lCasePass
        logical :: lOrder3Available(nSteps),lOrder5Available(nSteps)
        logical :: lOrderGradientAvailable(nSteps),lOrderHvAvailable(nSteps)
        type(FDSweepAssessment) :: tAssess3,tAssess5,tAssessGradient,tAssessHv

        n=SIZE(dState)
        allocate(dHessian(n,n),dHRef(n,n),dHIdeal(n,n),dHEx(n,n),dGradient(n),dDirection(n), &
            dPlus(n),dMinus(n),dPlus2(n),dMinus2(n),dGradPlus(n),dGradMinus(n),dDummyH(n,n))
        allocate(dHv(n),dGradientFD(n),dHvFD(n))
        call CompMQMQAGibbsEnergyUnconstrained(tData,dState,dIdealScale,tTerm,dG,dGRef,dGIdeal,dGEx,iInfo)
        lCasePass=iInfo==0
        if (.NOT.lCasePass) then
            if (lVerbose) write(*,'(/,A,A,I0)') TRIM(cName),': scalar iInfo = ',iInfo
            lAllPass=.FALSE.
            return
        end if
        call CompMQMQAHessianUnconstrained(tData,dState,dIdealScale,tTerm,dHessian,iInfo,dGDual,dGradient, &
            dHRef,dHIdeal,dHEx,dGDualRef,dGDualIdeal,dGDualEx)
        lCasePass=lCasePass.AND.(iInfo==0)
        dScaleError=MAX(NormalizedDifference(dG,dGDual),NormalizedDifference(dGRef,dGDualRef), &
            NormalizedDifference(dGIdeal,dGDualIdeal),NormalizedDifference(dGEx,dGDualEx))
        lCasePass=lCasePass.AND.(dScaleError<=1D-12)
        lCasePass=lCasePass.AND.ALL(IEEE_IS_FINITE(dHessian)).AND.ALL(IEEE_IS_FINITE(dGradient))

        ! These are exact structural expectations for a smooth extensive energy.
        ! Symmetry means differentiating first with respect to species p and then
        ! q gives the same result in the opposite order. Homogeneity means scaling
        ! every quadruplet amount changes phase amount but not composition, so the
        ! Hessian applied to the current mole vector should be zero.
        dSym=FrobeniusNorm(dHessian-TRANSPOSE(dHessian))/MAX(1D0,FrobeniusNorm(dHessian))
        dNormH=MatrixTwoNorm(dHessian)
        dNormN=SQRT(DOT_PRODUCT(dState,dState))
        dHom=SQRT(DOT_PRODUCT(MATMUL(dHessian,dState),MATMUL(dHessian,dState))) / &
            MAX(1D0,dNormH*dNormN)
        lCasePass=lCasePass.AND.(dSym<=1D-12).AND.(dHom<=1D-10)

        ! A direction is a vector assigning a simultaneous mole perturbation to
        ! every quadruplet species. This mixed direction changes many quadruplets
        ! together and exercises
        ! off-diagonal Hessian entries that basis-only checks could miss.
        do i=1,n
            dDirection(i)=(-1D0)**i*(0.3D0+0.07D0*MOD(i,5))
        end do
        dDirection=dDirection/SQRT(DOT_PRODUCT(dDirection,dDirection))
        dAnalytic=DOT_PRODUCT(dDirection,MATMUL(dHessian,dDirection))
        dHBase=0.08D0*MINVAL(dState/MAX(ABS(dDirection),1D-12))
        ! Sweep h because truncation error dominates coarse perturbations while
        ! cancellation dominates very fine perturbations. Both three- and
        ! five-point energy stencils must show a genuine decreasing region.
        do k=1,nSteps
            dH=dHBase/(2D0**(k-1))
            dStep(k)=dH
            dPlus=dState+dH*dDirection; dMinus=dState-dH*dDirection
            dPlus2=dState+2D0*dH*dDirection; dMinus2=dState-2D0*dH*dDirection
            call ScalarValue(tData,dPlus,dIdealScale,tTerm,dGPlus,iInfo)
            call ScalarValue(tData,dMinus,dIdealScale,tTerm,dGMinus,iInfo)
            call ScalarValue(tData,dPlus2,dIdealScale,tTerm,dTmp1,iInfo)
            call ScalarValue(tData,dMinus2,dIdealScale,tTerm,dTmp2,iInfo)
            dFD3=(dGPlus-2D0*dG+dGMinus)/(dH*dH)
            dFD5=(-dTmp1+16D0*dGPlus-30D0*dG+16D0*dGMinus-dTmp2)/(12D0*dH*dH)
            dAbsErr3(k)=ABS(dAnalytic-dFD3)
            dAbsErr5(k)=ABS(dAnalytic-dFD5)
            dErr3(k)=NormalizedDifference(dAnalytic,dFD3)
            dErr5(k)=NormalizedDifference(dAnalytic,dFD5)
        end do
        call AssessFDSweep(dStep,dErr3,FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX,1D-6, &
            tAssess3,dOrder3,lOrder3Available)
        call AssessFDSweep(dStep,dErr5,FD_ORDER_FOURTH_MIN,FD_ORDER_FOURTH_MAX,1D-8, &
            tAssess5,dOrder5,lOrder5Available)
        lCasePass=lCasePass.AND.tAssess3%lPassed.AND.tAssess5%lPassed

        ! Check each analytic chemical-potential component against the independent
        ! scalar energy, one mole variable at a time. A common positivity-safe
        ! step gives one vector approximation at each refinement level.
        dGradientHBase=0.08D0*MINVAL(dState)
        do k=1,nSteps
            dH=dGradientHBase/(2D0**(k-1))
            do i=1,n
                dPlus=dState; dMinus=dState
                dPlus(i)=dPlus(i)+dH; dMinus(i)=dMinus(i)-dH
                call ScalarValue(tData,dPlus,dIdealScale,tTerm,dGPlus,iInfo)
                call ScalarValue(tData,dMinus,dIdealScale,tTerm,dGMinus,iInfo)
                dGradientFD(i)=(dGPlus-dGMinus)/(2D0*dH)
            end do
            call ComputeVectorErrorMetrics(dGradientFD,dGradient,dGradNormAbs(k),dGradNormScaled(k), &
                dGradMaxAbs(k),dGradMaxScaled(k),iGradWorst(k))
        end do
        call AssessFDSweep([(dGradientHBase/(2D0**(k-1)),k=1,nSteps)],dGradNormScaled, &
            FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX,1D-7,tAssessGradient, &
            dOrderGradient,lOrderGradientAvailable)
        dGradientError=tAssessGradient%dBestError
        lCasePass=lCasePass.AND.tAssessGradient%lPassed

        ! Multiplying H by direction v predicts the change of the complete
        ! gradient under that simultaneous mole perturbation. Finite-differencing
        ! analytic gradients checks this prediction and is numerically better
        ! conditioned than taking a second energy difference.
        ! Moving the state along dDirection changes each chemical potential.
        ! Multiplying the analytic Hessian by that mole direction predicts the
        ! complete vector of those first-order changes.
        dHv=MATMUL(dHessian,dDirection)
        do k=1,nSteps
            dH=dStep(k)
            dPlus=dState+dH*dDirection; dMinus=dState-dH*dDirection
            call CompMQMQAHessianUnconstrained(tData,dPlus,dIdealScale,tTerm,dDummyH,iInfo,dGradient=dGradPlus)
            call CompMQMQAHessianUnconstrained(tData,dMinus,dIdealScale,tTerm,dDummyH,iInfo,dGradient=dGradMinus)
            dHvFD=(dGradPlus-dGradMinus)/(2D0*dH)
            call ComputeVectorErrorMetrics(dHvFD,dHv,dHvNormAbs(k),dHvNormScaled(k), &
                dHvMaxAbs(k),dHvMaxScaled(k),iHvWorst(k))
        end do
        call AssessFDSweep(dStep,dHvNormScaled,FD_ORDER_SECOND_MIN,FD_ORDER_SECOND_MAX,1D-7, &
            tAssessHv,dOrderHv,lOrderHvAvailable)
        dHvError=tAssessHv%dBestError
        lCasePass=lCasePass.AND.tAssessHv%lPassed

        ! Scaling every quadruplet amount by the same factor changes only the
        ! amount of phase, not its composition. An extensive Gibbs energy must
        ! therefore scale by that same factor.
        dExtError=0D0
        do j=1,2
            if (j==1) then; dTmp3=0.5D0; else; dTmp3=2D0; end if
            call ScalarValue(tData,dTmp3*dState,dIdealScale,tTerm,dGPlus,iInfo)
            dExtError=MAX(dExtError,NormalizedDifference(dGPlus,dTmp3*dG))
        end do
        lCasePass=lCasePass.AND.(dExtError<=1D-12)
        lAllPass=lAllPass.AND.lCasePass

        if (lVerbose) then
            write(*,'(/,A)') TRIM(cName)
            write(*,'(A,ES12.4)') 'scalar-path block error = ',dScaleError
            write(*,'(A,ES12.4)') 'gradient FD error = ',dGradientError
            write(*,'(A,ES12.4)') 'Hv FD error = ',dHvError
            write(*,'(A,ES12.4)') 'raw symmetry residual = ',dSym
            write(*,'(A,ES12.4)') 'homogeneity residual = ',dHom
            write(*,'(A,ES12.4)') 'extensivity error = ',dExtError
            write(*,'(A)') 'controlled standalone inputs; expected scalar orders: 3pt=2, 5pt=4'
            write(*,'(A)') 'h          3pt abs       3pt scaled    3pt order   5pt abs       5pt scaled    5pt order'
            do k=1,nSteps
                write(*,'(3ES14.5,2X,A,2ES14.5,2X,A)') dStep(k),dAbsErr3(k),dErr3(k), &
                    TRIM(OrderLabel(dOrder3(k),lOrder3Available(k))),dAbsErr5(k),dErr5(k), &
                    TRIM(OrderLabel(dOrder5(k),lOrder5Available(k)))
            end do
            write(*,'(A,ES12.4,A,L1)') '3pt best scaled error = ',tAssess3%dBestError, &
                ' roundoff upturn = ',tAssess3%lRoundoffUpturn
            write(*,'(A,ES12.4,A,L1)') '5pt best scaled error = ',tAssess5%dBestError, &
                ' roundoff upturn = ',tAssess5%lRoundoffUpturn
            write(*,'(A)') 'scalar-to-gradient central derivative; expected order = 2'
            write(*,'(A)') 'h          norm abs      norm scaled   max abs       max scaled    worst  order'
            do k=1,nSteps
                write(*,'(5ES14.5,I7,2X,A)') dGradientHBase/(2D0**(k-1)),dGradNormAbs(k), &
                    dGradNormScaled(k),dGradMaxAbs(k),dGradMaxScaled(k),iGradWorst(k), &
                    TRIM(OrderLabel(dOrderGradient(k),lOrderGradientAvailable(k)))
            end do
            write(*,'(A)') 'gradient-to-Hessian-vector central derivative; expected order = 2'
            write(*,'(A)') 'h          norm abs      norm scaled   max abs       max scaled    worst  order'
            do k=1,nSteps
                write(*,'(5ES14.5,I7,2X,A)') dStep(k),dHvNormAbs(k),dHvNormScaled(k), &
                    dHvMaxAbs(k),dHvMaxScaled(k),iHvWorst(k), &
                    TRIM(OrderLabel(dOrderHv(k),lOrderHvAvailable(k)))
            end do
            write(*,'(A,L1)') 'case pass = ',lCasePass
        end if

    end subroutine VerifyCase


    character(len=16) function OrderLabel(dOrder,lAvailable)

        real(8), intent(in) :: dOrder
        logical, intent(in) :: lAvailable

        if (lAvailable) then
            write(OrderLabel,'(F10.4)') dOrder
        else
            OrderLabel='N/A'
        end if

    end function OrderLabel


    subroutine ScalarValue(tData,dState,dIdealScale,tTerm,dValue,iInfo)

        type(MQMQAModelData), intent(in) :: tData
        real(8), intent(in) :: dState(:),dIdealScale
        type(MQMQAInteractionTerm), intent(in) :: tTerm(:)
        real(8), intent(out) :: dValue
        integer, intent(out) :: iInfo
        real(8) :: dRef,dIdeal,dEx

        call CompMQMQAGibbsEnergyUnconstrained(tData,dState,dIdealScale,tTerm,dValue,dRef,dIdeal,dEx,iInfo)

    end subroutine ScalarValue

    !=========================================================================================================
    ! SECTION 3: MODEL-SPECIFIC IDENTITIES AND EXPECTED FAILURES
    !
    ! These checks cover requirements that a generic finite-difference sweep
    ! cannot establish by itself: the production-derived B prefactor and explicit
    ! rejection of unsupported or mathematically singular inputs.
    !=========================================================================================================


    !---------------------------------------------------------------------------------------------------------
    !> \brief Quantify the published-versus-production SUBQ S3 formulation difference.
    !>
    !> \details Current production evaluates the S3 quadruplet term with ordinary
    !!          pair fractions. Published Equations (5)--(6) define X_i/k from
    !!          pair amounts containing 1/zeta_i/k; Equations (16) and (29) use
    !!          those normalized zeta-weighted fractions in S3, while Equation
    !!          (31) retains the corresponding zeta-dependent derivative. This
    !!          diagnostic evaluates both definitions and finite-differences their
    !!          signed difference without deciding why production diverges.
    !---------------------------------------------------------------------------------------------------------
    subroutine VerifyS3FormulationDifference(tData,dState,lAllPass,lVerbose)

        type(MQMQAModelData), intent(in) :: tData
        real(8), intent(in) :: dState(:)
        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lVerbose

        real(8), allocatable :: dGradientCoarse(:),dGradientFine(:)
        real(8), allocatable :: dHessianCoarse(:,:),dHessianFine(:,:)
        real(8) :: dProduction,dPaper,dDelta,dEnergyScaled,dGradientScaled,dHessianScaled
        real(8) :: dGradientUncertainty,dHessianUncertainty,dStep,dGradientMax,dHessianMax
        integer :: n,iInfo,iGradientWorst,iHessianRow,iHessianColumn,i,j
        logical :: lCasePass

        n=SIZE(dState)
        allocate(dGradientCoarse(n),dGradientFine(n),dHessianCoarse(n,n),dHessianFine(n,n))
        call EvaluateS3Definitions(tData,dState,dProduction,dPaper,iInfo)
        lCasePass=iInfo==0
        if (.NOT.lCasePass) then
            if (lVerbose) write(*,'(/,A,I0)') 'SUBQ S3 formulation diagnostic iInfo = ',iInfo
            lAllPass=.FALSE.
            return
        end if

        dDelta=dPaper-dProduction
        dEnergyScaled=ABS(dDelta)/MAX(1D0,ABS(dPaper),ABS(dProduction))
        dStep=1D-3*MINVAL(dState)
        call NumericalS3DifferenceDerivatives(tData,dState,dStep,dGradientCoarse,dHessianCoarse,iInfo)
        lCasePass=lCasePass.AND.(iInfo==0)
        call NumericalS3DifferenceDerivatives(tData,dState,0.5D0*dStep,dGradientFine,dHessianFine,iInfo)
        lCasePass=lCasePass.AND.(iInfo==0)

        dGradientScaled=SQRT(DOT_PRODUCT(dGradientFine,dGradientFine))/MAX(1D0,ABS(dDelta))
        dHessianScaled=FrobeniusNorm(dHessianFine)/MAX(1D0,ABS(dDelta))
        dGradientUncertainty=SQRT(DOT_PRODUCT(dGradientFine-dGradientCoarse, &
            dGradientFine-dGradientCoarse))/MAX(1D0,SQRT(DOT_PRODUCT(dGradientFine,dGradientFine)))
        dHessianUncertainty=FrobeniusNorm(dHessianFine-dHessianCoarse)/ &
            MAX(1D0,FrobeniusNorm(dHessianFine))

        iGradientWorst=MAXLOC(ABS(dGradientFine),DIM=1)
        dGradientMax=ABS(dGradientFine(iGradientWorst))
        iHessianRow=1; iHessianColumn=1; dHessianMax=0D0
        do j=1,n
            do i=1,n
                if (ABS(dHessianFine(i,j))>dHessianMax) then
                    dHessianMax=ABS(dHessianFine(i,j))
                    iHessianRow=i; iHessianColumn=j
                end if
            end do
        end do

        ! The diagnostic must be non-vacuous and numerically resolved. These
        ! gates do not choose between the paper and production conventions.
        lCasePass=lCasePass.AND.ALL(IEEE_IS_FINITE(dGradientFine)) &
            .AND.ALL(IEEE_IS_FINITE(dHessianFine))
        lCasePass=lCasePass.AND.(dEnergyScaled>1D-10).AND.(dGradientScaled>1D-10) &
            .AND.(dHessianScaled>1D-10)
        lCasePass=lCasePass.AND.(dGradientUncertainty<=0.05D0*dGradientScaled) &
            .AND.(dHessianUncertainty<=0.10D0*dHessianScaled)
        lAllPass=lAllPass.AND.lCasePass

        if (lVerbose) then
            write(*,'(/,A)') 'SUBQ S3 paper-versus-production formulation diagnostic'
            write(*,'(A)') 'production S3 uses ordinary pair fractions; paper S3 uses zeta-weighted fractions'
            write(*,'(A,ES14.6)') 'production S3 = ',dProduction
            write(*,'(A,ES14.6)') 'paper S3 = ',dPaper
            write(*,'(A,ES14.6)') 'signed Delta S3 (paper-production) = ',dDelta
            write(*,'(A,ES12.4)') 'scaled energy difference = ',dEnergyScaled
            write(*,'(A,ES12.4)') 'gradient difference norm (scaled) = ',dGradientScaled
            write(*,'(A,ES12.4,A,I0)') 'maximum gradient difference = ',dGradientMax, &
                ' at quadruplet ',iGradientWorst
            write(*,'(A,ES12.4)') 'gradient refinement disagreement = ',dGradientUncertainty
            write(*,'(A,ES12.4)') 'Hessian difference norm (scaled) = ',dHessianScaled
            write(*,'(A,ES12.4,A,I0,A,I0,A)') 'maximum Hessian difference = ',dHessianMax, &
                ' at (',iHessianRow,',',iHessianColumn,')'
            write(*,'(A,ES12.4)') 'Hessian refinement disagreement = ',dHessianUncertainty
            write(*,'(A,L1)') 'formulations numerically distinct and resolved = ',lCasePass
        end if

    end subroutine VerifyS3FormulationDifference


    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate production-style and paper-style SUBQ S3 scalars independently.
    !>
    !> \details The two expressions differ only in the four pair fractions inside
    !!          each quadruplet logarithm. Both use SUBQ theta=3/4, psi=1/2,
    !!          identical equivalent fractions, and the same quadruplet multiplicity.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateS3Definitions(tData,dState,dProduction,dPaper,iInfo)

        type(MQMQAModelData), intent(in) :: tData
        real(8), intent(in) :: dState(:)
        real(8), intent(out) :: dProduction,dPaper
        integer, intent(out) :: iInfo

        integer :: n,q,i,j,a,b,x,y,nA,nX,iPosition,jPosition,iWeight
        real(8) :: dN,dOrdinarySum,dWeightedSum,dProdPairLog,dPaperPairLog,dEquivalentLog
        real(8), allocatable :: dFraction(:),dEquivalent1(:),dEquivalent2(:)
        real(8), allocatable :: dOrdinary(:,:),dWeighted(:,:),dXOrdinary(:,:),dXWeighted(:,:)

        iInfo=0; dProduction=0D0; dPaper=0D0
        n=SIZE(dState)
        if ((tData%iModelType/=MQMQA_MODEL_SUBQ).OR.(n/=SIZE(tData%iQuadruplet,1)) &
            .OR.ANY(dState<=0D0)) then
            iInfo=1
            return
        end if
        allocate(dFraction(n),dEquivalent1(tData%nSublattice1),dEquivalent2(tData%nSublattice2), &
            dOrdinary(tData%nSublattice1,tData%nSublattice2), &
            dWeighted(tData%nSublattice1,tData%nSublattice2), &
            dXOrdinary(tData%nSublattice1,tData%nSublattice2), &
            dXWeighted(tData%nSublattice1,tData%nSublattice2))

        dN=SUM(dState); dFraction=dState/dN
        dEquivalent1=0D0; dEquivalent2=0D0; dOrdinary=0D0; dWeighted=0D0
        do q=1,n
            a=tData%iQuadruplet(q,1); b=tData%iQuadruplet(q,2)
            x=tData%iQuadruplet(q,3); y=tData%iQuadruplet(q,4)
            dEquivalent1(a)=dEquivalent1(a)+0.5D0*dFraction(q)
            dEquivalent1(b)=dEquivalent1(b)+0.5D0*dFraction(q)
            dEquivalent2(x)=dEquivalent2(x)+0.5D0*dFraction(q)
            dEquivalent2(y)=dEquivalent2(y)+0.5D0*dFraction(q)
            do i=1,tData%nSublattice1
                nA=MERGE(1,0,a==i)+MERGE(1,0,b==i)
                do j=1,tData%nSublattice2
                    nX=MERGE(1,0,x==j)+MERGE(1,0,y==j)
                    dOrdinary(i,j)=dOrdinary(i,j)+dState(q)*DFLOAT(nA*nX)
                    dWeighted(i,j)=dWeighted(i,j)+dState(q)*DFLOAT(nA*nX)/tData%dZeta(i,j)
                end do
            end do
        end do
        dOrdinarySum=SUM(dOrdinary); dWeightedSum=SUM(dWeighted)
        if ((dOrdinarySum<=0D0).OR.(dWeightedSum<=0D0).OR.ANY(dEquivalent1<=0D0) &
            .OR.ANY(dEquivalent2<=0D0)) then
            iInfo=2
            return
        end if
        dXOrdinary=dOrdinary/dOrdinarySum
        dXWeighted=dWeighted/dWeightedSum

        do q=1,n
            a=tData%iQuadruplet(q,1); b=tData%iQuadruplet(q,2)
            x=tData%iQuadruplet(q,3); y=tData%iQuadruplet(q,4)
            iWeight=1
            if (a/=b) iWeight=2*iWeight
            if (x/=y) iWeight=2*iWeight
            dProdPairLog=0D0; dPaperPairLog=0D0
            do iPosition=1,2
                do jPosition=3,4
                    i=tData%iQuadruplet(q,iPosition)
                    j=tData%iQuadruplet(q,jPosition)
                    if ((dXOrdinary(i,j)<=0D0).OR.(dXWeighted(i,j)<=0D0)) then
                        iInfo=3
                        return
                    end if
                    dProdPairLog=dProdPairLog+DLOG(dXOrdinary(i,j))
                    dPaperPairLog=dPaperPairLog+DLOG(dXWeighted(i,j))
                end do
            end do
            dEquivalentLog=DLOG(dEquivalent1(a))+DLOG(dEquivalent1(b))+ &
                DLOG(dEquivalent2(x))+DLOG(dEquivalent2(y))
            dProduction=dProduction+dState(q)*(DLOG(dFraction(q))-DLOG(DFLOAT(iWeight)) &
                -0.75D0*dProdPairLog+0.5D0*dEquivalentLog)
            dPaper=dPaper+dState(q)*(DLOG(dFraction(q))-DLOG(DFLOAT(iWeight)) &
                -0.75D0*dPaperPairLog+0.5D0*dEquivalentLog)
        end do

    end subroutine EvaluateS3Definitions


    !---------------------------------------------------------------------------------------------------------
    !> \brief Finite-difference the signed paper-minus-production S3 discrepancy.
    !---------------------------------------------------------------------------------------------------------
    subroutine NumericalS3DifferenceDerivatives(tData,dState,dH,dGradient,dHessian,iInfo)

        type(MQMQAModelData), intent(in) :: tData
        real(8), intent(in) :: dState(:),dH
        real(8), intent(out) :: dGradient(:),dHessian(:,:)
        integer, intent(out) :: iInfo

        integer :: n,i,j,iLocalInfo
        real(8) :: dBase,dPlus,dMinus,dPP,dPM,dMP,dMM,dProd,dPaper
        real(8), allocatable :: dTrial(:)

        iInfo=0; n=SIZE(dState); allocate(dTrial(n))
        call EvaluateS3Definitions(tData,dState,dProd,dPaper,iLocalInfo)
        if (iLocalInfo/=0) then; iInfo=iLocalInfo; return; end if
        dBase=dPaper-dProd
        do i=1,n
            dTrial=dState; dTrial(i)=dTrial(i)+dH
            call EvaluateS3Definitions(tData,dTrial,dProd,dPaper,iLocalInfo); dPlus=dPaper-dProd
            dTrial=dState; dTrial(i)=dTrial(i)-dH
            call EvaluateS3Definitions(tData,dTrial,dProd,dPaper,iLocalInfo); dMinus=dPaper-dProd
            if (iLocalInfo/=0) then; iInfo=iLocalInfo; return; end if
            dGradient(i)=(dPlus-dMinus)/(2D0*dH)
            dHessian(i,i)=(dPlus-2D0*dBase+dMinus)/(dH*dH)
        end do
        do j=2,n
            do i=1,j-1
                dTrial=dState; dTrial(i)=dTrial(i)+dH; dTrial(j)=dTrial(j)+dH
                call EvaluateS3Definitions(tData,dTrial,dProd,dPaper,iLocalInfo); dPP=dPaper-dProd
                dTrial=dState; dTrial(i)=dTrial(i)+dH; dTrial(j)=dTrial(j)-dH
                call EvaluateS3Definitions(tData,dTrial,dProd,dPaper,iLocalInfo); dPM=dPaper-dProd
                dTrial=dState; dTrial(i)=dTrial(i)-dH; dTrial(j)=dTrial(j)+dH
                call EvaluateS3Definitions(tData,dTrial,dProd,dPaper,iLocalInfo); dMP=dPaper-dProd
                dTrial=dState; dTrial(i)=dTrial(i)-dH; dTrial(j)=dTrial(j)-dH
                call EvaluateS3Definitions(tData,dTrial,dProd,dPaper,iLocalInfo); dMM=dPaper-dProd
                if (iLocalInfo/=0) then; iInfo=iLocalInfo; return; end if
                dHessian(i,j)=(dPP-dPM-dMP+dMM)/(4D0*dH*dH)
                dHessian(j,i)=dHessian(i,j)
            end do
        end do

    end subroutine NumericalS3DifferenceDerivatives


    !---------------------------------------------------------------------------------------------------------
    !> \brief Compare the S3 formulation difference with complete controlled SUBQ scales.
    !>
    !> \details Three strictly positive states with the same total phase amount
    !!          probe the baseline composition and two deterministic composition
    !!          skews. The production-style total is evaluated by the complete
    !!          standalone SUBQ model. Because the paper-style interpretation
    !!          changes only S3, its complete energy, gradient, and Hessian equal
    !!          the production totals plus the independently finite-differenced
    !!          paper-minus-production S3 difference.
    !!
    !!          These ratios measure numerical significance within this controlled
    !!          model. They do not measure equilibrium or phase-stability effects
    !!          in an assessed database.
    !---------------------------------------------------------------------------------------------------------
    subroutine VerifyS3TotalSignificance(tData,dBaseState,tTerm,lAllPass,lVerbose)

        type(MQMQAModelData), intent(in) :: tData
        real(8), intent(in) :: dBaseState(:)
        type(MQMQAInteractionTerm), intent(in) :: tTerm(:)
        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lVerbose

        integer, parameter :: nCases=3
        character(len=24), parameter :: cCaseName(nCases)=[character(len=24) :: &
            'baseline composition','graded composition skew','alternating composition']
        real(8), allocatable :: dState(:,:),dGradientProduction(:),dGradientDeltaCoarse(:)
        real(8), allocatable :: dGradientDelta(:),dHessianProduction(:,:),dHessianDeltaCoarse(:,:)
        real(8), allocatable :: dHessianDelta(:,:),dGradientPaper(:),dHessianPaper(:,:)
        real(8) :: dGProduction,dGPaper,dS3Production,dS3Paper,dDeltaG,dTotalAmount,dStep
        real(8) :: dEnergyRatio,dGradientRatio,dHessianRatio,dGradientUncertainty,dHessianUncertainty
        real(8) :: dMinEnergyRatio,dMaxEnergyRatio,dMinGradientRatio,dMaxGradientRatio
        real(8) :: dMinHessianRatio,dMaxHessianRatio,dWeight
        integer :: n,i,iCase,iInfo,iLocalInfo
        logical :: lCasePass,lAllCasesPass

        n=SIZE(dBaseState)
        allocate(dState(n,nCases),dGradientProduction(n),dGradientDeltaCoarse(n), &
            dGradientDelta(n),dHessianProduction(n,n),dHessianDeltaCoarse(n,n), &
            dHessianDelta(n,n),dGradientPaper(n),dHessianPaper(n,n))

        dTotalAmount=SUM(dBaseState)
        dState(:,1)=dBaseState
        do i=1,n
            if (n>1) then
                dWeight=0.55D0+0.90D0*DFLOAT(i-1)/DFLOAT(n-1)
            else
                dWeight=1D0
            end if
            dState(i,2)=dBaseState(i)*dWeight
            dState(i,3)=dBaseState(i)*(0.55D0+0.18D0*DFLOAT(MOD(i,5)))
        end do
        do iCase=2,nCases
            dState(:,iCase)=dState(:,iCase)*dTotalAmount/SUM(dState(:,iCase))
        end do

        dMinEnergyRatio=HUGE(1D0); dMaxEnergyRatio=0D0
        dMinGradientRatio=HUGE(1D0); dMaxGradientRatio=0D0
        dMinHessianRatio=HUGE(1D0); dMaxHessianRatio=0D0
        lAllCasesPass=.TRUE.

        if (lVerbose) then
            write(*,'(/,A)') 'SUBQ S3 significance relative to the complete controlled model'
            write(*,'(A)') 'All states are synthetic, positive, nonuniform-zeta states with equal total phase amount.'
            write(*,'(A)') 'Ratios compare paper-minus-production S3 changes with complete production-style totals.'
            write(*,'(A)') 'The energy ratio depends on the chosen reference-energy zero; derivative ratios do not.'
            write(*,'(A)') 'state                       G production       Delta G   |Delta G|/|G|  '// &
                '||Delta mu||/||mu||  ||Delta H||/||H||'
        end if

        do iCase=1,nCases
            call CompMQMQAHessianUnconstrained(tData,dState(:,iCase),1D0,tTerm, &
                dHessianProduction,iInfo,dGibbs=dGProduction,dGradient=dGradientProduction)
            lCasePass=iInfo==0
            call EvaluateS3Definitions(tData,dState(:,iCase),dS3Production,dS3Paper,iLocalInfo)
            lCasePass=lCasePass.AND.(iLocalInfo==0)
            dDeltaG=dS3Paper-dS3Production

            dStep=1D-3*MINVAL(dState(:,iCase))
            call NumericalS3DifferenceDerivatives(tData,dState(:,iCase),dStep, &
                dGradientDeltaCoarse,dHessianDeltaCoarse,iLocalInfo)
            lCasePass=lCasePass.AND.(iLocalInfo==0)
            call NumericalS3DifferenceDerivatives(tData,dState(:,iCase),0.5D0*dStep, &
                dGradientDelta,dHessianDelta,iLocalInfo)
            lCasePass=lCasePass.AND.(iLocalInfo==0)

            dGPaper=dGProduction+dDeltaG
            dGradientPaper=dGradientProduction+dGradientDelta
            dHessianPaper=dHessianProduction+dHessianDelta
            dEnergyRatio=ABS(dDeltaG)/MAX(1D0,ABS(dGProduction),ABS(dGPaper))
            dGradientRatio=VectorNorm(dGradientDelta)/ &
                MAX(1D0,VectorNorm(dGradientProduction),VectorNorm(dGradientPaper))
            dHessianRatio=FrobeniusNorm(dHessianDelta)/ &
                MAX(1D0,FrobeniusNorm(dHessianProduction),FrobeniusNorm(dHessianPaper))
            dGradientUncertainty=VectorNorm(dGradientDelta-dGradientDeltaCoarse)/ &
                MAX(1D0,VectorNorm(dGradientDelta))
            dHessianUncertainty=FrobeniusNorm(dHessianDelta-dHessianDeltaCoarse)/ &
                MAX(1D0,FrobeniusNorm(dHessianDelta))

            lCasePass=lCasePass.AND.IEEE_IS_FINITE(dGProduction).AND.IEEE_IS_FINITE(dGPaper) &
                .AND.ALL(IEEE_IS_FINITE(dGradientPaper)).AND.ALL(IEEE_IS_FINITE(dHessianPaper))
            lCasePass=lCasePass.AND.(dEnergyRatio>1D-12).AND.(dGradientRatio>1D-12) &
                .AND.(dHessianRatio>1D-12)
            lCasePass=lCasePass.AND.(dGradientUncertainty<=0.05D0*dGradientRatio) &
                .AND.(dHessianUncertainty<=0.10D0*dHessianRatio)
            lAllCasesPass=lAllCasesPass.AND.lCasePass

            dMinEnergyRatio=MIN(dMinEnergyRatio,dEnergyRatio)
            dMaxEnergyRatio=MAX(dMaxEnergyRatio,dEnergyRatio)
            dMinGradientRatio=MIN(dMinGradientRatio,dGradientRatio)
            dMaxGradientRatio=MAX(dMaxGradientRatio,dGradientRatio)
            dMinHessianRatio=MIN(dMinHessianRatio,dHessianRatio)
            dMaxHessianRatio=MAX(dMaxHessianRatio,dHessianRatio)

            if (lVerbose) then
                write(*,'(A24,5ES19.6)') cCaseName(iCase),dGProduction,dDeltaG,dEnergyRatio, &
                    dGradientRatio,dHessianRatio
                write(*,'(A,2ES12.4,A,L1)') '  derivative refinement disagreements (mu,H) = ', &
                    dGradientUncertainty,dHessianUncertainty,' resolved = ',lCasePass
            end if
        end do

        lAllPass=lAllPass.AND.lAllCasesPass
        if (lVerbose) then
            write(*,'(A,2ES12.4)') 'energy-ratio range = ',dMinEnergyRatio,dMaxEnergyRatio
            write(*,'(A,2ES12.4)') 'chemical-potential-ratio range = ',dMinGradientRatio,dMaxGradientRatio
            write(*,'(A,2ES12.4)') 'Hessian-ratio range = ',dMinHessianRatio,dMaxHessianRatio
            write(*,'(A,L1)') 'complete-model significance diagnostic pass = ',lAllCasesPass
        end if

    end subroutine VerifyS3TotalSignificance


    !---------------------------------------------------------------------------------------------------------
    !> \brief Prove that nonuniform SUBQ zeta values affect the intended dependent quantities.
    !>
    !> \details A uniform-zeta SUBQ copy is compared with the nonuniform model at
    !!          the same positive mole state. The configurational comparison tests
    !!          the weighted pair distribution and its F marginals in S2. The
    !!          B-only comparison tests direct use of the same weighted pair
    !!          fractions. Requiring nonzero energy and Hessian differences prevents
    !!          pair-specific zeta support from being merely accepted as unused data.
    !---------------------------------------------------------------------------------------------------------
    subroutine VerifyPairSpecificZeta(tData,dState,tNone,tB,lAllPass,lVerbose)

        type(MQMQAModelData), intent(in) :: tData
        real(8), intent(in) :: dState(:)
        type(MQMQAInteractionTerm), intent(in) :: tNone(:),tB
        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lVerbose

        type(MQMQAModelData) :: tUniform
        type(MQMQAInteractionTerm) :: tOnly(1)
        real(8), allocatable :: dHNonuniform(:,:),dHUniform(:,:)
        real(8) :: dG,dRef,dIdealNonuniform,dIdealUniform,dEx,dEnergyDifference,dHessianDifference
        real(8) :: dBNonuniform,dBUniform,dBEnergyDifference,dBHessianDifference
        integer :: iInfo,n

        n=SIZE(dState)
        allocate(dHNonuniform(n,n),dHUniform(n,n))
        tUniform=tData
        tUniform%dZeta=SUM(tData%dZeta)/DFLOAT(SIZE(tData%dZeta))

        call CompMQMQAGibbsEnergyUnconstrained(tData,dState,1D0,tNone,dG,dRef,dIdealNonuniform,dEx,iInfo)
        lAllPass=lAllPass.AND.(iInfo==0)
        call CompMQMQAGibbsEnergyUnconstrained(tUniform,dState,1D0,tNone,dG,dRef,dIdealUniform,dEx,iInfo)
        lAllPass=lAllPass.AND.(iInfo==0)
        call CompMQMQAHessianUnconstrained(tData,dState,1D0,tNone,dHNonuniform,iInfo)
        lAllPass=lAllPass.AND.(iInfo==0)
        call CompMQMQAHessianUnconstrained(tUniform,dState,1D0,tNone,dHUniform,iInfo)
        lAllPass=lAllPass.AND.(iInfo==0)
        dEnergyDifference=NormalizedDifference(dIdealNonuniform,dIdealUniform)
        dHessianDifference=FrobeniusNorm(dHNonuniform-dHUniform)/ &
            MAX(1D0,FrobeniusNorm(dHNonuniform),FrobeniusNorm(dHUniform))

        tOnly(1)=tB
        call CompMQMQAGibbsEnergyUnconstrained(tData,dState,0D0,tOnly,dG,dRef,dIdealNonuniform,dBNonuniform,iInfo)
        lAllPass=lAllPass.AND.(iInfo==0)
        call CompMQMQAGibbsEnergyUnconstrained(tUniform,dState,0D0,tOnly,dG,dRef,dIdealUniform,dBUniform,iInfo)
        lAllPass=lAllPass.AND.(iInfo==0)
        call CompMQMQAHessianUnconstrained(tData,dState,0D0,tOnly,dHNonuniform,iInfo)
        lAllPass=lAllPass.AND.(iInfo==0)
        call CompMQMQAHessianUnconstrained(tUniform,dState,0D0,tOnly,dHUniform,iInfo)
        lAllPass=lAllPass.AND.(iInfo==0)
        dBEnergyDifference=NormalizedDifference(dBNonuniform,dBUniform)
        dBHessianDifference=FrobeniusNorm(dHNonuniform-dHUniform)/ &
            MAX(1D0,FrobeniusNorm(dHNonuniform),FrobeniusNorm(dHUniform))

        lAllPass=lAllPass.AND.(dEnergyDifference>1D-8).AND.(dHessianDifference>1D-8) &
            .AND.(dBEnergyDifference>1D-8).AND.(dBHessianDifference>1D-8)
        if (lVerbose) then
            write(*,'(/,A)') 'SUBQ pair-specific-zeta non-vacuity'
            write(*,'(A,ES12.4)') 'configurational energy difference = ',dEnergyDifference
            write(*,'(A,ES12.4)') 'configurational Hessian difference = ',dHessianDifference
            write(*,'(A,ES12.4)') 'B energy difference = ',dBEnergyDifference
            write(*,'(A,ES12.4)') 'B Hessian difference = ',dBHessianDifference
        end if

    end subroutine VerifyPairSpecificZeta


    !---------------------------------------------------------------------------------------------------------
    !> \brief Verify why B-family energy equals total quadruplet moles times its composition-only local modifier.
    !>
    !> \details N_Q denotes the sum of all quadruplet mole amounts. This check
    !!          independently differentiates the weighted-pair formula
    !!          in the same direct-plus-zeta structure used by production chemical
    !!          potentials. Agreement rules out choosing N_Q merely because it
    !!          makes the standalone finite differences self-consistent.
    !---------------------------------------------------------------------------------------------------------
    subroutine VerifyBProductionIdentity(cName,tData,dState,tB,lAllPass,lVerbose)

        character(*), intent(in) :: cName
        type(MQMQAModelData), intent(in) :: tData
        real(8), intent(in) :: dState(:)
        type(MQMQAInteractionTerm), intent(in) :: tB
        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lVerbose

        type(MQMQAModelData) :: tZeroReference
        type(MQMQAInteractionTerm) :: tOnly(1)
        integer :: n,q,i,j,a,b,x,y,nA,nX,iInfo
        real(8) :: dN,dS,dWA,dWB,dWAB,dF,dError
        real(8), allocatable :: dWeighted(:,:),dDerivativeW(:,:),dExpected(:),dGradient(:),dH(:,:)

        tZeroReference=tData
        tZeroReference%dReferenceEnergy=0D0
        tOnly(1)=tB
        n=SIZE(dState)
        allocate(dWeighted(tData%nSublattice1,tData%nSublattice2), &
            dDerivativeW(tData%nSublattice1,tData%nSublattice2),dExpected(n),dGradient(n),dH(n,n))
        ! Reconstruct the unnormalized zeta-weighted pair amounts. Their derivatives
        ! with respect to each quadruplet are simple incidence/zeta factors.
        dWeighted=0D0
        ! The B energy is total quadruplet moles N_Q multiplied by Delta g_B,
        ! a local modifier determined only by normalized weighted-pair
        ! composition. Differentiating this product produces one direct
        ! modifier contribution and additional contributions from how every
        ! weighted-pair amount changes.
        do q=1,n
            a=tData%iQuadruplet(q,1); b=tData%iQuadruplet(q,2)
            x=tData%iQuadruplet(q,3); y=tData%iQuadruplet(q,4)
            do i=1,tData%nSublattice1
                nA=MERGE(1,0,a==i)+MERGE(1,0,b==i)
                do j=1,tData%nSublattice2
                    nX=MERGE(1,0,x==j)+MERGE(1,0,y==j)
                    dWeighted(i,j)=dWeighted(i,j)+dState(q)*DFLOAT(nA*nX)/tData%dZeta(i,j)
                end do
            end do
        end do
        dN=SUM(dState); dS=SUM(dWeighted)
        dWA=dWeighted(tB%iA,tB%iX); dWB=dWeighted(tB%iB,tB%iY); dWAB=dWA+dWB
        dF=tB%dCoefficient*(dWA/dS)**(1+tB%iExponentP)*(dWB/dS)**(1+tB%iExponentQ) / &
            ((dWA+dWB)/dS)**(1+tB%iExponentP+tB%iExponentQ)
        dDerivativeW=-dF/dS
        dDerivativeW(tB%iA,tB%iX)=dDerivativeW(tB%iA,tB%iX)+ &
            dF*(DFLOAT(1+tB%iExponentP)/dWA-DFLOAT(1+tB%iExponentP+tB%iExponentQ)/dWAB)
        dDerivativeW(tB%iB,tB%iY)=dDerivativeW(tB%iB,tB%iY)+ &
            dF*(DFLOAT(1+tB%iExponentQ)/dWB-DFLOAT(1+tB%iExponentP+tB%iExponentQ)/dWAB)
        do q=1,n
            a=tData%iQuadruplet(q,1); b=tData%iQuadruplet(q,2)
            x=tData%iQuadruplet(q,3); y=tData%iQuadruplet(q,4)
            dExpected(q)=dF
            do i=1,tData%nSublattice1
                nA=MERGE(1,0,a==i)+MERGE(1,0,b==i)
                do j=1,tData%nSublattice2
                    nX=MERGE(1,0,x==j)+MERGE(1,0,y==j)
                    dExpected(q)=dExpected(q)+dN*dDerivativeW(i,j)*DFLOAT(nA*nX)/tData%dZeta(i,j)
                end do
            end do
        end do
        call CompMQMQAHessianUnconstrained(tZeroReference,dState,0D0,tOnly,dH,iInfo,dGradient=dGradient)
        dError=0D0
        do q=1,n
            dError=MAX(dError,NormalizedDifference(dExpected(q),dGradient(q)))
        end do
        lAllPass=lAllPass.AND.(iInfo==0).AND.(dError<=1D-12)
        if (lVerbose) write(*,'(/,A,A,ES12.4)') TRIM(cName), &
            ' B production-gradient identity error = ',dError

    end subroutine VerifyBProductionIdentity


    !> Confirm that unresolved R terms and boundary mole states fail explicitly.
    subroutine VerifyFailures(tData,dState,tTemplate,lAllPass,lVerbose)

        type(MQMQAModelData), intent(in) :: tData
        real(8), intent(in) :: dState(:)
        type(MQMQAInteractionTerm), intent(in) :: tTemplate
        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lVerbose

        type(MQMQAModelData) :: tBadModel
        type(MQMQAInteractionTerm) :: tBad(1)
        real(8), allocatable :: dBadState(:)
        real(8) :: dG,dRef,dIdeal,dEx
        integer :: iInfo

        tBad(1)=tTemplate
        tBad(1)%iFamily=MQMQA_TERM_R
        call CompMQMQAGibbsEnergyUnconstrained(tData,dState,1D0,tBad,dG,dRef,dIdeal,dEx,iInfo)
        lAllPass=lAllPass.AND.(iInfo==11)
        dBadState=dState
        dBadState(1)=0D0
        tBad(1)=tTemplate
        call CompMQMQAGibbsEnergyUnconstrained(tData,dBadState,1D0,tBad,dG,dRef,dIdeal,dEx,iInfo)
        lAllPass=lAllPass.AND.(iInfo==4)
        tBadModel=tData
        tBadModel%dZeta(1,1)=1.25D0*tBadModel%dZeta(1,1)
        call CompMQMQAGibbsEnergyUnconstrained(tBadModel,dState,1D0,tBad,dG,dRef,dIdeal,dEx,iInfo)
        lAllPass=lAllPass.AND.(iInfo==7)
        tBadModel=tData
        tBadModel%iModelType=0
        call CompMQMQAGibbsEnergyUnconstrained(tBadModel,dState,1D0,tBad,dG,dRef,dIdeal,dEx,iInfo)
        lAllPass=lAllPass.AND.(iInfo==6)
        if (lVerbose) write(*,'(/,A)') &
            'failure checks: R, boundary state, nonuniform SUBG zeta, and unknown formulation rejected'

    end subroutine VerifyFailures

    !=========================================================================================================
    ! SECTION 4: SCALE-NORMALIZED NUMERICAL HELPERS
    !
    ! Centralize norms and convergence-trend checks so every thermodynamic case
    ! is judged using the same scale-aware criteria.
    !=========================================================================================================


    real(8) function NormalizedDifference(dA,dB)
        real(8), intent(in) :: dA,dB
        NormalizedDifference=ABS(dA-dB)/MAX(1D0,ABS(dA),ABS(dB))
    end function NormalizedDifference


    real(8) function FrobeniusNorm(dA)
        real(8), intent(in) :: dA(:,:)
        FrobeniusNorm=SQRT(SUM(dA*dA))
    end function FrobeniusNorm


    real(8) function VectorNorm(dA)
        real(8), intent(in) :: dA(:)
        VectorNorm=SQRT(DOT_PRODUCT(dA,dA))
    end function VectorNorm


    real(8) function MatrixTwoNorm(dA)
        real(8), intent(in) :: dA(:,:)
        integer :: i
        real(8) :: dNorm
        real(8), allocatable :: dVector(:),dNext(:)

        allocate(dVector(SIZE(dA,2)),dNext(SIZE(dA,2)))
        dVector=1D0/SQRT(DFLOAT(SIZE(dVector)))
        do i=1,100
            dNext=MATMUL(TRANSPOSE(dA),MATMUL(dA,dVector))
            dNorm=SQRT(DOT_PRODUCT(dNext,dNext))
            if (dNorm<=TINY(1D0)) then
                MatrixTwoNorm=0D0
                return
            end if
            dVector=dNext/dNorm
        end do
        dNext=MATMUL(dA,dVector)
        MatrixTwoNorm=SQRT(DOT_PRODUCT(dNext,dNext))
    end function MatrixTwoNorm


end program TestMQMQAHessianVerification
