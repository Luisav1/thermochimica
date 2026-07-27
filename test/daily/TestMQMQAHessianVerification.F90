!-------------------------------------------------------------------------------------------------------------
!> \file    TestMQMQAHessianVerification.F90
!> \brief   Standalone verification of the supported nonmagnetic SUBG Hessian.
!>
!> \details Independent ordinary-real and second-order-object evaluators are compared before finite
!!          differences test the analytic gradient and Hessian. The cases isolate ordinary
!!          configurational mixing, the G and Q binary parameter families, all three
!!          supported ternary-group branches, and the B weighted-pair family. Pass
!!          --report for the complete numerical evidence.
!!
!!          Verification map:
!!          1. Check the generic second-order calculus kernel.
!!          2. Build one complete positive synthetic SUBG topology.
!!          3. Define isolated G, Q, ternary, and B interaction cases.
!!          4. Apply the scalar/gradient/Hessian verification ladder to each case.
!!          5. Check the independently derived extensive B identity.
!!          6. Confirm unsupported and boundary inputs fail explicitly.
!-------------------------------------------------------------------------------------------------------------

program TestMQMQAHessianVerification

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
    USE ModuleMQMQAUnconstrained
    USE ModuleFiniteDifferenceVerification

    implicit none

    type(MQMQAModelData) :: tModel
    type(MQMQAInteractionTerm), allocatable :: tAll(:), tOne(:), tNone(:)
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
    call BuildInteractions(tAll)
    allocate(tNone(0),tOne(1))

    if (lReport) then
        write(*,'(A)') 'Stage MQ-1 supported nonmagnetic SUBG Hessian verification'
        write(*,'(A,I0)') 'storage bits = ',STORAGE_SIZE(1D0)
        write(*,'(A,I0)') 'decimal precision = ',PRECISION(1D0)
        write(*,'(A,I0)') 'binary digits = ',DIGITS(1D0)
        write(*,'(A,ES14.6)') 'machine epsilon = ',EPSILON(1D0)
    end if
    lPass=lPass.AND.(STORAGE_SIZE(1D0)==64).AND.(PRECISION(1D0)>=15).AND.(DIGITS(1D0)>=53)

    ! Add one physical/mathematical layer at a time. This localizes a failure to
    ! configurational mixing, one production family, or one ternary branch.
    call VerifyCase('reference + ideal',tModel,dMoles,1D0,tNone,lPass,lReport)
    tOne(1)=tAll(1); call VerifyCase('G binary',tModel,dMoles,1D0,tOne,lPass,lReport)
    tOne(1)=tAll(2); call VerifyCase('Q binary',tModel,dMoles,1D0,tOne,lPass,lReport)
    tOne(1)=tAll(3); call VerifyCase('ternary group 1',tModel,dMoles,1D0,tOne,lPass,lReport)
    tOne(1)=tAll(4); call VerifyCase('ternary group 2',tModel,dMoles,1D0,tOne,lPass,lReport)
    tOne(1)=tAll(5); call VerifyCase('ternary neither group',tModel,dMoles,1D0,tOne,lPass,lReport)
    tOne(1)=tAll(6); call VerifyCase('B family',tModel,dMoles,1D0,tOne,lPass,lReport)
    call VerifyCase('integrated total',tModel,dMoles,1D0,tAll,lPass,lReport)
    call VerifyBProductionIdentity(tModel,dMoles,tAll(6),lPass,lReport)
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
    !!          binary subset. Nonuniform moles, coordination numbers, zeta values,
    !!          and reference energies prevent accidental cancellation from making
    !!          an incorrect derivative appear correct.
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
        tData%dZeta=RESHAPE([3.8D0,4.0D0,4.2D0,4.4D0,4.1D0,4.3D0,4.5D0,4.7D0],[4,2])

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
    !> \brief Verify why B-family energy equals total quadruplet moles times its composition-only local modifier.
    !>
    !> \details N_Q denotes the sum of all quadruplet mole amounts. This check
    !!          independently differentiates the weighted-pair formula
    !!          in the same direct-plus-zeta structure used by production chemical
    !!          potentials. Agreement rules out choosing N_Q merely because it
    !!          makes the standalone finite differences self-consistent.
    !---------------------------------------------------------------------------------------------------------
    subroutine VerifyBProductionIdentity(tData,dState,tB,lAllPass,lVerbose)

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
        if (lVerbose) write(*,'(/,A,ES12.4)') 'B production-gradient identity error = ',dError

    end subroutine VerifyBProductionIdentity


    !> Confirm that unresolved R terms and boundary mole states fail explicitly.
    subroutine VerifyFailures(tData,dState,tTemplate,lAllPass,lVerbose)

        type(MQMQAModelData), intent(in) :: tData
        real(8), intent(in) :: dState(:)
        type(MQMQAInteractionTerm), intent(in) :: tTemplate
        logical, intent(inout) :: lAllPass
        logical, intent(in) :: lVerbose

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
        if (lVerbose) write(*,'(/,A)') 'failure checks: R rejected and boundary state rejected'

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
