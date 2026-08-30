!-------------------------------------------------------------------------------------------------------------
!> \file    ModuleMQMQAUnconstrained.f90
!> \brief   Disconnected Hessian for SUBG and staged nonmagnetic SUBQ scalar forms.
!>
!> \details The module evaluates reference, traced configurational, G-family, Q-family, B-family,
!!          and supported ternary scalar terms from generic arrays. It neither reads Thermochimica global
!!          state nor maps curvature into GEMNewton. The ordinary-real energy evaluator is intentionally
!!          independent of the private value/gradient/Hessian evaluator used for analytic derivatives.
!!          
!!          Conceptual map of the energy:
!!          - Reference: the standard-state energy assigned to each quadruplet.
!!          - S1: entropy from distributing individual constituents over the two sublattices.
!!          - S2: correction for how unlike cation-anion pairs depart from independent random pairing.
!!          - S3: correction for how complete A-B-X-Y quadruplets depart from the distribution implied
!!            by their pair and constituent fractions.
!!          - G and Q: Thermochimica names for two excess-interaction parameter families. G uses a
!!            local binary coordinate (chi); Q uses a broader normalized projected coordinate (xi).
!!          - B: a third production interaction family written directly in terms of two
!!            zeta-weighted pair fractions.
!!          
!!          The labels G, Q, and B are database-family names. They do not mean gradient,
!!          heat, or the total Gibbs energy returned by this module.
!!
!!          Notation used below:
!!          - A and B label constituents on the first sublattice; X and Y label
!!            constituents on the second. Together A-B-X-Y identifies one
!!            quadruplet species, not a matrix expression.
!!          - n(p) is the mole amount of quadruplet species p, and N_Q is the sum
!!            of all quadruplet mole amounts in the phase.
!!          - An ordinary pair fraction describes how frequently one first-sublattice
!!            constituent is paired with one second-sublattice constituent.
!!          - zeta is fixed model data with zero mole derivatives. SUBG supplies one
!!            common zeta for every pair. SUBQ may supply a different zeta for each
!!            pair, creating a second, normalized zeta-weighted pair distribution.
!!            S2 and B use that weighted distribution. S3 uses ordinary pair
!!            fractions for SUBG and the weighted distribution for SUBQ, together
!!            with the model-specific theta/psi exponents.
!!          - chi and xi are normalized composition coordinates used by the G
!!            and Q parameter families; they are calculated from quadruplet
!!            populations and are not additional independent solver variables.
!!          - The gradient contains first derivatives of energy with respect to
!!            every n(p). The Hessian contains all corresponding second derivatives.
!!
!!          In MQMQAInteractionTerm, the fields named exponent P, Q, and R are
!!          integer powers from the production parameter definition. They are
!!          distinct from the G/Q/B family label and from the total N_Q.
!!
!!          File map:
!!          1. Public model data and interaction descriptions
!!          2. Verification of the private second-order calculus kernel
!!          3. Public scalar-energy and Hessian entry points
!!          4. Input, topology, and mathematical-domain checks
!!          5. Independent ordinary-real scalar evaluator
!!             - dependent MQMQA composition measures
!!             - configurational S1/S2/S3 terms
!!             - G, Q, B, and supported ternary excess terms
!!          6. Derivative-object evaluator of the same thermodynamic model
!!          7. Private second-order arithmetic and chain-rule primitives
!!
!!          Sections 5 and 6 intentionally implement the energy through separate
!!          expression paths. Finite differences of the ordinary-real path can
!!          therefore test the analytic derivatives without differentiating the
!!          same coded expression that produced them.
!-------------------------------------------------------------------------------------------------------------

module ModuleMQMQAUnconstrained

    implicit none
    private

    !=========================================================================================================
    ! SECTION 1: PUBLIC MODEL AND INTERACTION DESCRIPTIONS
    !
    ! These types describe a local MQMQA phase without reading ModuleThermo. They
    ! are the boundary between a future production-data adapter and the
    ! disconnected thermodynamic mathematics implemented below.
    !=========================================================================================================

    integer, parameter, public :: MQMQA_MODEL_UNSET = 0
    integer, parameter, public :: MQMQA_MODEL_SUBG = 1
    integer, parameter, public :: MQMQA_MODEL_SUBQ = 2

    integer, parameter, public :: MQMQA_TERM_G = 1
    integer, parameter, public :: MQMQA_TERM_Q = 2
    integer, parameter, public :: MQMQA_TERM_B = 3
    integer, parameter, public :: MQMQA_TERM_R = 4

    !> Test-only ablation of the corrected SUBQ S3 pair distribution.
    !>
    !> The default is the production zeta-weighted definition.  Private solver
    !! diagnostics may temporarily select the former ordinary-pair definition
    !! to determine whether the S3 correction changed MQMQA globalization.
    !! CompExcessGibbsEnergySUBG reads the same flag so production partial
    !! molars and the analytical scalar/Hessian always use one formulation.
    logical, public :: lMQMQADiagnosticLegacyS3 = .FALSE.

    !> Generic topology and constant data for one disconnected MQMQA phase.
    !>
    !> Each row of iQuadruplet identifies [A,B,X,Y], where A/B belong to the
    !> first sublattice and X/Y belong to the second. Coordination numbers
    !> convert quadruplet amounts into site amounts. Zeta values provide the
    !> pair weighting used by the modified quasichemical composition variables.
    type, public :: MQMQAModelData
        !> Configurational formulation. Callers must explicitly select SUBG or SUBQ.
        integer :: iModelType = MQMQA_MODEL_UNSET
        !> Number of constituents available on the cation-like sublattice.
        integer :: nSublattice1 = 0
        !> Number of constituents available on the anion-like sublattice.
        integer :: nSublattice2 = 0
        !> Canonically ordered quadruplet identities [A,B,X,Y].
        integer, allocatable :: iQuadruplet(:,:)
        !> Coordination number associated with each of the four quadruplet positions.
        real(8), allocatable :: dCoordination(:,:)
        !> Pair zeta values indexed by first- and second-sublattice constituent.
        !> All entries must be equal for SUBG; SUBQ permits pair-specific values.
        real(8), allocatable :: dZeta(:,:)
        !> Standard Gibbs energy assigned to each quadruplet.
        real(8), allocatable :: dReferenceEnergy(:)
    end type MQMQAModelData

    !> One traced SUBG excess interaction. Sublattice-2 indices are local, not offset.
    !>
    !> A/B/X/Y identify the binary interaction axis. P and Q are the binary
    !> composition powers. A positive ternary constituent activates the traced
    !> ternary modifier, whose branch is selected by the fixed group masks and
    !> whose order is R. B-family terms use P/Q but do not accept ternary data.
    type, public :: MQMQAInteractionTerm
        !> Production family: G, Q, B, or explicitly unsupported R.
        integer :: iFamily = 0
        !> First-sublattice binary endpoints.
        integer :: iA = 0
        integer :: iB = 0
        !> Second-sublattice binary endpoints.
        integer :: iX = 0
        integer :: iY = 0
        !> Nonnegative binary powers and the supported ternary order.
        integer :: iExponentP = 0
        integer :: iExponentQ = 0
        integer :: iExponentR = 0
        !> Zero for a binary term; otherwise the third constituent in a traced ternary branch.
        integer :: iTernaryConstituent = 0
        !> Temperature-evaluated interaction coefficient supplied by the caller.
        real(8) :: dCoefficient = 0D0
        !> Fixed asymmetric-group membership on the active binary sublattice.
        logical, allocatable :: lGroup1(:)
        logical, allocatable :: lGroup2(:)
    end type MQMQAInteractionTerm

    !> A scalar value carrying all first and second derivatives with respect to quadruplet moles.
    !>
    !> This private type acts like a small second-order automatic-differentiation
    !> object. Building the MQMQA equations from these objects applies the chain
    !> rule through every normalized site, pair, binary, and ternary quantity.
    type :: SecondOrderScalar
        real(8) :: dValue = 0D0
        real(8), allocatable :: dGradient(:)
        real(8), allocatable :: dHessian(:,:)
    end type SecondOrderScalar

    public :: CompMQMQAGibbsEnergyUnconstrained
    public :: CompMQMQAHessianUnconstrained
    public :: CheckMQMQADerivativeKernel

contains

    !=========================================================================================================
    ! SECTION 2: SECOND-ORDER CALCULUS KERNEL VERIFICATION
    !
    ! Verify the generic differentiation before it is asked to carry
    ! the much larger MQMQA expression. A failure here is an arithmetic or chain-
    ! rule failure, not evidence against the thermodynamic model.
    !=========================================================================================================

    !---------------------------------------------------------------------------------------------------------
    !> \brief Check every private second-order primitive against analytic and finite-difference controls.
    !> \param[out] iInfo Zero when every primitive passes; otherwise the failing operation number.
    !---------------------------------------------------------------------------------------------------------
    subroutine CheckMQMQADerivativeKernel(iInfo)

        integer, intent(out) :: iInfo
        real(8), parameter :: dX=1.7D0, dY=0.8D0
        real(8) :: dGradient(2), dHessian(2,2)
        type(SecondOrderScalar) :: tX, tY, tResult

        iInfo=0
        ! The kernel is checked independently of MQMQA. This separates errors in
        ! product/chain-rule propagation from errors in the thermodynamic equations.
        tX=VariableSO(dX,1,2)
        tY=VariableSO(dY,2,2)

        dGradient=0D0; dHessian=0D0
        call CheckPrimitive(ConstantSO(2.3D0,2),2.3D0,dGradient,dHessian,0,iInfo)
        if (iInfo/=0) return
        dGradient=[1D0,0D0]
        call CheckPrimitive(tX,dX,dGradient,dHessian,1,iInfo)
        if (iInfo/=0) return

        dGradient=[1D0,1D0]
        tResult=AddSO(tX,tY)
        call CheckPrimitive(tResult,dX+dY,dGradient,dHessian,2,iInfo)
        call CheckPrimitiveFD(tResult,dX,dY,2,iInfo)
        if (iInfo/=0) return

        dGradient=[1D0,-1D0]
        tResult=SubtractSO(tX,tY)
        call CheckPrimitive(tResult,dX-dY,dGradient,dHessian,3,iInfo)
        call CheckPrimitiveFD(tResult,dX,dY,3,iInfo)
        if (iInfo/=0) return

        dGradient=[0D0,-1D0]
        tResult=NegateSO(tY)
        call CheckPrimitive(tResult,-dY,dGradient,dHessian,4,iInfo)
        call CheckPrimitiveFD(tResult,dX,dY,4,iInfo)
        if (iInfo/=0) return

        dGradient=[2.5D0,0D0]
        tResult=ScaleSO(tX,2.5D0)
        call CheckPrimitive(tResult,2.5D0*dX,dGradient,dHessian,5,iInfo)
        call CheckPrimitiveFD(tResult,dX,dY,5,iInfo)
        if (iInfo/=0) return

        dGradient=[dY,dX]; dHessian=0D0
        dHessian(1,2)=1D0; dHessian(2,1)=1D0
        tResult=MultiplySO(tX,tY)
        call CheckPrimitive(tResult,dX*dY,dGradient,dHessian,6,iInfo)
        call CheckPrimitiveFD(tResult,dX,dY,6,iInfo)
        if (iInfo/=0) return

        dGradient=[-1D0/dX**2,0D0]; dHessian=0D0; dHessian(1,1)=2D0/dX**3
        tResult=ReciprocalSO(tX)
        call CheckPrimitive(tResult,1D0/dX,dGradient,dHessian,7,iInfo)
        call CheckPrimitiveFD(tResult,dX,dY,7,iInfo)
        if (iInfo/=0) return

        dGradient=[1D0/dY,-dX/dY**2]; dHessian=0D0
        dHessian(1,2)=-1D0/dY**2; dHessian(2,1)=dHessian(1,2); dHessian(2,2)=2D0*dX/dY**3
        tResult=DivideSO(tX,tY)
        call CheckPrimitive(tResult,dX/dY,dGradient,dHessian,8,iInfo)
        call CheckPrimitiveFD(tResult,dX,dY,8,iInfo)
        if (iInfo/=0) return

        dGradient=[1D0/dX,0D0]; dHessian=0D0; dHessian(1,1)=-1D0/dX**2
        tResult=LogSO(tX)
        call CheckPrimitive(tResult,LOG(dX),dGradient,dHessian,9,iInfo)
        call CheckPrimitiveFD(tResult,dX,dY,9,iInfo)
        if (iInfo/=0) return

        dGradient=[3D0*dX**2,0D0]; dHessian=0D0; dHessian(1,1)=6D0*dX
        tResult=PowerSO(tX,3)
        call CheckPrimitive(tResult,dX**3,dGradient,dHessian,10,iInfo)
        call CheckPrimitiveFD(tResult,dX,dY,10,iInfo)
        if (iInfo/=0) return

        dGradient=0D0; dHessian=0D0
        call CheckPrimitive(PowerSO(tX,0),1D0,dGradient,dHessian,11,iInfo)

    end subroutine CheckMQMQADerivativeKernel


    subroutine CheckPrimitive(tActual,dValue,dGradient,dHessian,iOperation,iInfo)

        type(SecondOrderScalar), intent(in) :: tActual
        real(8), intent(in) :: dValue, dGradient(:), dHessian(:,:)
        integer, intent(in) :: iOperation
        integer, intent(inout) :: iInfo

        if (iInfo/=0) return
        if ((ABS(tActual%dValue-dValue)>1D-12) .OR. &
            (MAXVAL(ABS(tActual%dGradient-dGradient))>1D-12) .OR. &
            (MAXVAL(ABS(tActual%dHessian-dHessian))>1D-12)) iInfo=100+iOperation

    end subroutine CheckPrimitive


    subroutine CheckPrimitiveFD(tActual,dX,dY,iOperation,iInfo)

        type(SecondOrderScalar), intent(in) :: tActual
        real(8), intent(in) :: dX,dY
        integer, intent(in) :: iOperation
        integer, intent(inout) :: iInfo
        real(8), parameter :: dDirection(2)=[0.6D0,-0.4D0], dH=1D-4
        real(8) :: dBase,dPlus,dMinus,dFirst,dSecond,dAnalyticFirst,dAnalyticSecond

        if (iInfo/=0) return
        dBase=PrimitiveScalar(dX,dY,iOperation)
        dPlus=PrimitiveScalar(dX+dH*dDirection(1),dY+dH*dDirection(2),iOperation)
        dMinus=PrimitiveScalar(dX-dH*dDirection(1),dY-dH*dDirection(2),iOperation)
        dFirst=(dPlus-dMinus)/(2D0*dH)
        dSecond=(dPlus-2D0*dBase+dMinus)/(dH*dH)
        dAnalyticFirst=DOT_PRODUCT(tActual%dGradient,dDirection)
        dAnalyticSecond=DOT_PRODUCT(dDirection,MATMUL(tActual%dHessian,dDirection))
        if ((ABS(dFirst-dAnalyticFirst)>1D-8*MAX(1D0,ABS(dAnalyticFirst))) .OR. &
            (ABS(dSecond-dAnalyticSecond)>1D-6*MAX(1D0,ABS(dAnalyticSecond)))) iInfo=200+iOperation

    end subroutine CheckPrimitiveFD


    real(8) function PrimitiveScalar(dX,dY,iOperation)

        real(8), intent(in) :: dX,dY
        integer, intent(in) :: iOperation

        select case(iOperation)
        case(2);  PrimitiveScalar=dX+dY
        case(3);  PrimitiveScalar=dX-dY
        case(4);  PrimitiveScalar=-dY
        case(5);  PrimitiveScalar=2.5D0*dX
        case(6);  PrimitiveScalar=dX*dY
        case(7);  PrimitiveScalar=1D0/dX
        case(8);  PrimitiveScalar=dX/dY
        case(9);  PrimitiveScalar=LOG(dX)
        case(10); PrimitiveScalar=dX**3
        case default; PrimitiveScalar=0D0
        end select

    end function PrimitiveScalar

    !=========================================================================================================
    ! SECTION 3: PUBLIC THERMODYNAMIC ENTRY POINTS
    !
    ! The scalar entry point is the independent numerical reference. The Hessian
    ! entry point evaluates the same physical model with derivative-carrying
    ! objects and optionally exposes its energy and gradient for verification.
    ! Neither entry point reads or changes Thermochimica global state.
    !=========================================================================================================

    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate the selected disconnected MQMQA scalar energy using ordinary real arithmetic.
    !> \param[in]  tModel          Generic MQMQA topology, formulation, and model constants.
    !> \param[in]  dMoles          Strictly positive quadruplet mole amounts.
    !> \param[in]  dIdealScale     Scale multiplying the traced configurational-energy expression.
    !> \param[in]  tInteraction    Supported nonmagnetic G, Q, and B interaction terms.
    !> \param[out] dGibbs          Total extensive Gibbs energy.
    !> \param[out] dGibbsReference Reference-energy contribution.
    !> \param[out] dGibbsIdeal     Configurational-energy contribution.
    !> \param[out] dGibbsExcess    Sum of supported excess contributions.
    !> \param[out] iInfo           Zero on success; nonzero for invalid or unsupported input.
    !---------------------------------------------------------------------------------------------------------
    subroutine CompMQMQAGibbsEnergyUnconstrained(tModel,dMoles,dIdealScale,tInteraction, &
        dGibbs,dGibbsReference,dGibbsIdeal,dGibbsExcess,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        real(8), intent(in) :: dMoles(:), dIdealScale
        type(MQMQAInteractionTerm), intent(in) :: tInteraction(:)
        real(8), intent(out) :: dGibbs, dGibbsReference, dGibbsIdeal, dGibbsExcess
        integer, intent(out) :: iInfo

        call EvaluateScalarEnergy(tModel,dMoles,dIdealScale,tInteraction,dGibbsReference, &
            dGibbsIdeal,dGibbsExcess,iInfo)
        dGibbs = dGibbsReference + dGibbsIdeal + dGibbsExcess

    end subroutine CompMQMQAGibbsEnergyUnconstrained


    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate the analytic quadruplet-mole Hessian with private second-order derivative objects.
    !>
    !> \details Optional energy and gradient outputs expose the derivative path for verification against the
    !!          independent ordinary-real evaluator. The raw Hessian is returned without symmetrization.
    !> \param[in]  tModel             Generic MQMQA topology, formulation, and model constants.
    !> \param[in]  dMoles             Strictly positive quadruplet mole amounts.
    !> \param[in]  dIdealScale        Scale multiplying the traced configurational-energy expression.
    !> \param[in]  tInteraction       Supported nonmagnetic G, Q, and B interaction terms.
    !> \param[out] dHessian           Raw total quadruplet-mole Hessian.
    !> \param[out] iInfo              Zero on success; nonzero for invalid or unsupported input.
    !> \param[out] dGibbs             Optional derivative-path total energy.
    !> \param[out] dGradient          Optional analytic quadruplet-mole gradient.
    !> \param[out] dHessianReference  Optional reference block Hessian.
    !> \param[out] dHessianIdeal      Optional configurational block Hessian.
    !> \param[out] dHessianExcess     Optional excess block Hessian.
    !> \param[out] dGibbsReference    Optional derivative-path reference energy.
    !> \param[out] dGibbsIdeal        Optional derivative-path configurational energy.
    !> \param[out] dGibbsExcess       Optional derivative-path excess energy.
    !---------------------------------------------------------------------------------------------------------
    subroutine CompMQMQAHessianUnconstrained(tModel,dMoles,dIdealScale,tInteraction,dHessian,iInfo, &
        dGibbs,dGradient,dHessianReference,dHessianIdeal,dHessianExcess, &
        dGibbsReference,dGibbsIdeal,dGibbsExcess)

        type(MQMQAModelData), intent(in) :: tModel
        real(8), intent(in) :: dMoles(:), dIdealScale
        type(MQMQAInteractionTerm), intent(in) :: tInteraction(:)
        real(8), intent(out) :: dHessian(:,:)
        integer, intent(out) :: iInfo
        real(8), intent(out), optional :: dGibbs, dGradient(:)
        real(8), intent(out), optional :: dHessianReference(:,:), dHessianIdeal(:,:), dHessianExcess(:,:)
        real(8), intent(out), optional :: dGibbsReference, dGibbsIdeal, dGibbsExcess

        integer :: nQuad
        type(SecondOrderScalar) :: tReference, tIdeal, tExcess, tTotal

        nQuad = SIZE(dMoles)
        dHessian = 0D0
        if (PRESENT(dGibbs)) dGibbs = 0D0
        if (PRESENT(dGradient)) dGradient = 0D0
        if (PRESENT(dHessianReference)) dHessianReference = 0D0
        if (PRESENT(dHessianIdeal)) dHessianIdeal = 0D0
        if (PRESENT(dHessianExcess)) dHessianExcess = 0D0
        if (PRESENT(dGibbsReference)) dGibbsReference = 0D0
        if (PRESENT(dGibbsIdeal)) dGibbsIdeal = 0D0
        if (PRESENT(dGibbsExcess)) dGibbsExcess = 0D0

        if ((SIZE(dHessian,1) /= nQuad) .OR. (SIZE(dHessian,2) /= nQuad)) then
            iInfo = 30
            return
        end if
        if (PRESENT(dGradient)) then
            if (SIZE(dGradient) /= nQuad) then
                iInfo = 31
                return
            end if
        end if
        if (PRESENT(dHessianReference)) then
            if ((SIZE(dHessianReference,1) /= nQuad) .OR. (SIZE(dHessianReference,2) /= nQuad)) then
                iInfo = 32
                return
            end if
        end if
        if (PRESENT(dHessianIdeal)) then
            if ((SIZE(dHessianIdeal,1) /= nQuad) .OR. (SIZE(dHessianIdeal,2) /= nQuad)) then
                iInfo = 33
                return
            end if
        end if
        if (PRESENT(dHessianExcess)) then
            if ((SIZE(dHessianExcess,1) /= nQuad) .OR. (SIZE(dHessianExcess,2) /= nQuad)) then
                iInfo = 34
                return
            end if
        end if

        call EvaluateDerivativeEnergy(tModel,dMoles,dIdealScale,tInteraction,tReference,tIdeal,tExcess,iInfo)
        if (iInfo /= 0) return
        tTotal = AddSO(AddSO(tReference,tIdeal),tExcess)

        dHessian = tTotal%dHessian
        if (PRESENT(dGibbs)) dGibbs = tTotal%dValue
        if (PRESENT(dGradient)) dGradient = tTotal%dGradient
        if (PRESENT(dHessianReference)) dHessianReference = tReference%dHessian
        if (PRESENT(dHessianIdeal)) dHessianIdeal = tIdeal%dHessian
        if (PRESENT(dHessianExcess)) dHessianExcess = tExcess%dHessian
        if (PRESENT(dGibbsReference)) dGibbsReference = tReference%dValue
        if (PRESENT(dGibbsIdeal)) dGibbsIdeal = tIdeal%dValue
        if (PRESENT(dGibbsExcess)) dGibbsExcess = tExcess%dValue

    end subroutine CompMQMQAHessianUnconstrained

    !=========================================================================================================
    ! SECTION 4: INPUT, TOPOLOGY, AND DOMAIN CHECKS
    !
    ! The local Hessian exists only for a valid canonical quadruplet topology and
    ! a strictly positive interior composition. Unsupported production families
    ! and untraced ternary orientations are rejected here rather than approximated.
    !=========================================================================================================


    !---------------------------------------------------------------------------------------------------------
    !> \brief Enforce the mathematical domain assumed by the disconnected MQMQA equations.
    !>
    !> \details The analytic Hessian is an interior-state object. Zero constituent
    !!          populations would make logarithms or normalized ratios singular, so
    !!          invalid states are rejected rather than clipped or renormalized.
    !---------------------------------------------------------------------------------------------------------
    subroutine CheckInputs(tModel,dMoles,tInteraction,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        real(8), intent(in) :: dMoles(:)
        type(MQMQAInteractionTerm), intent(in) :: tInteraction(:)
        integer, intent(out) :: iInfo

        integer :: i, nQuad, nMask

        iInfo = 0
        nQuad = SIZE(dMoles)
        ! The model must contain both sublattices and at least one quadruplet species.
        if ((tModel%nSublattice1 <= 0) .OR. (tModel%nSublattice2 <= 0) .OR. (nQuad <= 0)) then
            iInfo = 1
            return
        end if
        if ((tModel%iModelType /= MQMQA_MODEL_SUBG) .AND. &
            (tModel%iModelType /= MQMQA_MODEL_SUBQ)) then
            iInfo = 6
            return
        end if
        if (.NOT. ALLOCATED(tModel%iQuadruplet) .OR. .NOT. ALLOCATED(tModel%dCoordination) .OR. &
            .NOT. ALLOCATED(tModel%dZeta) .OR. .NOT. ALLOCATED(tModel%dReferenceEnergy)) then
            iInfo = 2
            return
        end if
        if ((SIZE(tModel%iQuadruplet,1) /= nQuad) .OR. (SIZE(tModel%iQuadruplet,2) /= 4) .OR. &
            (SIZE(tModel%dCoordination,1) /= nQuad) .OR. (SIZE(tModel%dCoordination,2) /= 4) .OR. &
            (SIZE(tModel%dZeta,1) /= tModel%nSublattice1) .OR. &
            (SIZE(tModel%dZeta,2) /= tModel%nSublattice2) .OR. &
            (SIZE(tModel%dReferenceEnergy) /= nQuad)) then
            iInfo = 3
            return
        end if
        ! Positive values keep every logarithm, reciprocal, and normalized pair
        ! measure inside the smooth domain where a finite Hessian exists.
        if (ANY(dMoles <= 0D0) .OR. ANY(tModel%dCoordination <= 0D0) .OR. ANY(tModel%dZeta <= 0D0)) then
            iInfo = 4
            return
        end if
        ! Thermochimica parses one common FNN/SNN ratio for SUBG but one value
        ! per constituent pair for SUBQ. Enforcing that distinction prevents a
        ! nominal SUBG fixture from silently exercising SUBQ weighting.
        if (tModel%iModelType == MQMQA_MODEL_SUBG) then
            if (MAXVAL(ABS(tModel%dZeta-tModel%dZeta(1,1))) > &
                1D-12*MAX(1D0,ABS(tModel%dZeta(1,1)))) then
                iInfo = 7
                return
            end if
        end if
        if (ANY(tModel%iQuadruplet(:,1:2) < 1) .OR. &
            ANY(tModel%iQuadruplet(:,1:2) > tModel%nSublattice1) .OR. &
            ANY(tModel%iQuadruplet(:,3:4) < 1) .OR. &
            ANY(tModel%iQuadruplet(:,3:4) > tModel%nSublattice2) .OR. &
            ANY(tModel%iQuadruplet(:,1) > tModel%iQuadruplet(:,2)) .OR. &
            ANY(tModel%iQuadruplet(:,3) > tModel%iQuadruplet(:,4))) then
            iInfo = 5
            return
        end if
        ! MQ-1 accepts only the production interaction branches whose underlying
        ! scalar energies have been traced. Unknown families are not guessed.
        do i = 1, SIZE(tInteraction)
            if ((tInteraction(i)%iFamily < MQMQA_TERM_G) .OR. &
                (tInteraction(i)%iFamily > MQMQA_TERM_R)) then
                iInfo = 10
                return
            end if
            if (tInteraction(i)%iFamily == MQMQA_TERM_R) then
                ! R is a production label whose scalar form remains outside MQ-1.
                iInfo = 11
                return
            end if
            if ((tInteraction(i)%iExponentP < 0) .OR. (tInteraction(i)%iExponentQ < 0) .OR. &
                (tInteraction(i)%iExponentR < 0)) then
                iInfo = 12
                return
            end if
            if ((tInteraction(i)%iA < 1) .OR. (tInteraction(i)%iA > tModel%nSublattice1) .OR. &
                (tInteraction(i)%iB < 1) .OR. (tInteraction(i)%iB > tModel%nSublattice1) .OR. &
                (tInteraction(i)%iX < 1) .OR. (tInteraction(i)%iX > tModel%nSublattice2) .OR. &
                (tInteraction(i)%iY < 1) .OR. (tInteraction(i)%iY > tModel%nSublattice2) .OR. &
                (tInteraction(i)%iA > tInteraction(i)%iB) .OR. &
                (tInteraction(i)%iX > tInteraction(i)%iY)) then
                iInfo = 13
                return
            end if
            if (tInteraction(i)%iFamily == MQMQA_TERM_B) then
                ! No production ternary extension has been traced for B terms.
                if (tInteraction(i)%iTernaryConstituent /= 0) then
                    iInfo = 14
                    return
                end if
                cycle
            end if
            ! A G/Q interaction varies on exactly one sublattice. The other
            ! sublattice is diagonal and supplies the fixed local environment.
            if ((tInteraction(i)%iA == tInteraction(i)%iB) .EQV. &
                (tInteraction(i)%iX == tInteraction(i)%iY)) then
                iInfo = 15
                return
            end if
            if (.NOT. ALLOCATED(tInteraction(i)%lGroup1) .OR. &
                .NOT. ALLOCATED(tInteraction(i)%lGroup2)) then
                iInfo = 16
                return
            end if
            if (tInteraction(i)%iX == tInteraction(i)%iY) then
                nMask = tModel%nSublattice1
            else
                nMask = tModel%nSublattice2
            end if
            if ((SIZE(tInteraction(i)%lGroup1) /= nMask) .OR. &
                (SIZE(tInteraction(i)%lGroup2) /= nMask) .OR. &
                ANY(tInteraction(i)%lGroup1 .AND. tInteraction(i)%lGroup2)) then
                iInfo = 17
                return
            end if
            if ((.NOT. ANY(tInteraction(i)%lGroup1)) .OR. (.NOT. ANY(tInteraction(i)%lGroup2))) then
                iInfo = 18
                return
            end if
            ! The supported ternary orientation adds a third first-sublattice
            ! constituent while the second-sublattice environment remains diagonal.
            if (tInteraction(i)%iTernaryConstituent > 0) then
                if ((tInteraction(i)%iX /= tInteraction(i)%iY) .OR. &
                    (tInteraction(i)%iTernaryConstituent > tModel%nSublattice1) .OR. &
                    (tInteraction(i)%iExponentR < 1)) then
                    iInfo = 19
                    return
                end if
            end if
        end do

    end subroutine CheckInputs

    !=========================================================================================================
    ! SECTION 5: INDEPENDENT ORDINARY-REAL SCALAR EVALUATOR
    !
    ! This path constructs the MQMQA state and energy using ordinary real numbers.
    ! It is deliberately independent of the derivative-object path in Section 6,
    ! making it suitable as the scalar oracle for finite-difference verification.
    !=========================================================================================================


    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate the independent ordinary-real MQMQA energy path.
    !>
    !> \details This path deliberately does not carry derivatives. Finite differences
    !!          of this independently coded scalar energy therefore test the analytic
    !!          derivative-object path instead of merely repeating it.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateScalarEnergy(tModel,dMoles,dIdealScale,tInteraction,dReference,dIdeal,dExcess,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        real(8), intent(in) :: dMoles(:), dIdealScale
        type(MQMQAInteractionTerm), intent(in) :: tInteraction(:)
        real(8), intent(out) :: dReference, dIdeal, dExcess
        integer, intent(out) :: iInfo

        integer :: a, x, q, iPosition, jPosition, iWeight
        real(8) :: dN, dS1, dS2, dS3, dDen, dTheta, dPsi
        real(8) :: dPairLogBlock, dEquivalentLogBlock, dLogRatio
        real(8), allocatable :: dQuadFraction(:), dSiteAmount1(:), dSiteAmount2(:), dSiteFraction1(:), dSiteFraction2(:)
        real(8), allocatable :: dEquivalent1(:), dEquivalent2(:), dPairAmount(:,:), dPairWeightedAmount(:,:)
        real(8), allocatable :: dPairFraction(:,:), dPairWeightedFraction(:,:), dF1(:), dF2(:)

        dReference = 0D0
        dIdeal = 0D0
        dExcess = 0D0
        call CheckInputs(tModel,dMoles,tInteraction,iInfo)
        if (iInfo /= 0) return

        call AllocateScalarState(tModel,dMoles,dQuadFraction,dSiteAmount1,dSiteAmount2,dSiteFraction1, &
            dSiteFraction2,dEquivalent1,dEquivalent2,dPairAmount,dPairWeightedAmount,dPairFraction, &
            dPairWeightedFraction,dF1,dF2,dN)

        ! Reference energy: each quadruplet amount multiplies its standard-state
        ! energy. This block is linear in moles and therefore has zero curvature.
        dReference = DOT_PRODUCT(dMoles,tModel%dReferenceEnergy)

        ! S1 is ordinary sublattice mixing entropy. It measures the combinatorial
        ! cost of distributing individual constituents over each sublattice.
        dS1 = 0D0
        do a = 1, tModel%nSublattice1
            if (dSiteFraction1(a) <= 0D0) then
                iInfo = 40
                return
            end if
            dS1 = dS1 + dSiteAmount1(a)*DLOG(dSiteFraction1(a))
        end do
        do x = 1, tModel%nSublattice2
            if (dSiteFraction2(x) <= 0D0) then
                iInfo = 40
                return
            end if
            dS1 = dS1 + dSiteAmount2(x)*DLOG(dSiteFraction2(x))
        end do

        ! S2 measures pair ordering. A random pair distribution would equal the
        ! product of its cation and anion marginals, F_A*F_X. The logarithmic ratio
        ! records how the actual zeta-weighted A-X pairs depart from that limit.
        dS2 = 0D0
        do a = 1, tModel%nSublattice1
            do x = 1, tModel%nSublattice2
                dDen = dF1(a)*dF2(x)
                if ((dPairWeightedFraction(a,x) <= 0D0) .OR. (dDen <= 0D0)) then
                    iInfo = 41
                    return
                end if
                dS2 = dS2 + dPairWeightedAmount(a,x)*DLOG(dPairWeightedFraction(a,x)/dDen)
            end do
        end do

        ! S3 measures ordering at the complete quadruplet level. The pair-log block
        ! contains the four A-X, A-Y, B-X, and B-Y pair fractions. The equivalent-
        ! fraction block contains the two first-sublattice and two second-sublattice
        ! constituent fractions. SUBG gives both blocks unit weight. SUBQ retains the
        ! same logarithmic derivative structure but uses theta=3/4 for the pair block
        ! and psi=1/2 for the equivalent-fraction block.
        dTheta = 0D0
        dPsi = 0D0
        select case(tModel%iModelType)
        case(MQMQA_MODEL_SUBG)
            dTheta = 1D0
            dPsi = 1D0
        case(MQMQA_MODEL_SUBQ)
            dTheta = 3D0/4D0
            dPsi = 1D0/2D0
        end select
        dS3 = 0D0
        do q = 1, SIZE(dMoles)
            iWeight = 1
            if (tModel%iQuadruplet(q,1) /= tModel%iQuadruplet(q,2)) iWeight = 2*iWeight
            if (tModel%iQuadruplet(q,3) /= tModel%iQuadruplet(q,4)) iWeight = 2*iWeight
            dPairLogBlock = 0D0
            do iPosition = 1, 2
                do jPosition = 3, 4
                    if ((tModel%iModelType == MQMQA_MODEL_SUBG) .OR. &
                        ((tModel%iModelType == MQMQA_MODEL_SUBQ) .AND. lMQMQADiagnosticLegacyS3)) then
                        dDen = dPairFraction(tModel%iQuadruplet(q,iPosition), &
                            tModel%iQuadruplet(q,jPosition))
                    else
                        ! SUBQ S3 uses the normalized zeta-weighted X_i/k
                        ! distribution defined by the updated-MQMQA equations.
                        dDen = dPairWeightedFraction(tModel%iQuadruplet(q,iPosition), &
                            tModel%iQuadruplet(q,jPosition))
                    end if
                    if (dDen <= 0D0) then
                        iInfo = 42
                        return
                    end if
                    dPairLogBlock = dPairLogBlock+DLOG(dDen)
                end do
            end do
            if ((dQuadFraction(q) <= 0D0) .OR. &
                (dEquivalent1(tModel%iQuadruplet(q,1)) <= 0D0) .OR. &
                (dEquivalent1(tModel%iQuadruplet(q,2)) <= 0D0) .OR. &
                (dEquivalent2(tModel%iQuadruplet(q,3)) <= 0D0) .OR. &
                (dEquivalent2(tModel%iQuadruplet(q,4)) <= 0D0)) then
                iInfo = 42
                return
            end if
            dEquivalentLogBlock = DLOG(dEquivalent1(tModel%iQuadruplet(q,1)))+ &
                DLOG(dEquivalent1(tModel%iQuadruplet(q,2)))+ &
                DLOG(dEquivalent2(tModel%iQuadruplet(q,3)))+ &
                DLOG(dEquivalent2(tModel%iQuadruplet(q,4)))
            dLogRatio = DLOG(dQuadFraction(q))-DLOG(DFLOAT(iWeight))- &
                dTheta*dPairLogBlock+dPsi*dEquivalentLogBlock
            dS3 = dS3+dMoles(q)*dLogRatio
        end do
        dIdeal = dIdealScale*(dS1+dS2+dS3)

        ! The nonideal interaction families are added after the shared composition
        ! measures have been built. These terms change curvature beyond S1-S3.
        call EvaluateScalarExcess(tModel,dMoles,dN,dQuadFraction,dPairWeightedFraction, &
            tInteraction,dExcess,iInfo)

    end subroutine EvaluateScalarEnergy


    !---------------------------------------------------------------------------------------------------------
    !> \brief Convert independent quadruplet moles into all dependent MQMQA composition measures.
    !>
    !> \details A quadruplet amount says how much A-B-X-Y local structure is present.
    !!          The energy equations also need constituent site fractions, pair
    !!          probabilities, zeta-weighted pair probabilities, and their marginal
    !!          sums. None of these are independent variables; they all change when
    !!          any quadruplet mole amount changes.
    !---------------------------------------------------------------------------------------------------------
    subroutine AllocateScalarState(tModel,dMoles,dQuadFraction,dSiteAmount1,dSiteAmount2,dSiteFraction1, &
        dSiteFraction2,dEquivalent1,dEquivalent2,dPairAmount,dPairWeightedAmount,dPairFraction, &
        dPairWeightedFraction,dF1,dF2,dN)

        type(MQMQAModelData), intent(in) :: tModel
        real(8), intent(in) :: dMoles(:)
        real(8), allocatable, intent(out) :: dQuadFraction(:), dSiteAmount1(:), dSiteAmount2(:)
        real(8), allocatable, intent(out) :: dSiteFraction1(:), dSiteFraction2(:), dEquivalent1(:), dEquivalent2(:)
        real(8), allocatable, intent(out) :: dPairAmount(:,:), dPairWeightedAmount(:,:), dPairFraction(:,:)
        real(8), allocatable, intent(out) :: dPairWeightedFraction(:,:), dF1(:), dF2(:)
        real(8), intent(out) :: dN

        integer :: q, a, b, x, y, i, j, nA, nX
        real(8) :: dSum1, dSum2, dPairSum, dWeightedSum

        allocate(dQuadFraction(SIZE(dMoles)),dSiteAmount1(tModel%nSublattice1), &
            dSiteAmount2(tModel%nSublattice2),dSiteFraction1(tModel%nSublattice1), &
            dSiteFraction2(tModel%nSublattice2),dEquivalent1(tModel%nSublattice1), &
            dEquivalent2(tModel%nSublattice2),dPairAmount(tModel%nSublattice1,tModel%nSublattice2), &
            dPairWeightedAmount(tModel%nSublattice1,tModel%nSublattice2), &
            dPairFraction(tModel%nSublattice1,tModel%nSublattice2), &
            dPairWeightedFraction(tModel%nSublattice1,tModel%nSublattice2), &
            dF1(tModel%nSublattice1),dF2(tModel%nSublattice2))
        dSiteAmount1 = 0D0; dSiteAmount2 = 0D0; dEquivalent1 = 0D0; dEquivalent2 = 0D0
        dPairAmount = 0D0; dPairWeightedAmount = 0D0; dF1 = 0D0; dF2 = 0D0
        ! N_Q is the sum of all quadruplet mole amounts and therefore carries the
        ! overall amount of the phase. Dividing each amount by N_Q produces
        ! composition fractions that remain unchanged if the whole phase is scaled.
        dN = SUM(dMoles)
        dQuadFraction = dMoles/dN

        ! Project each quadruplet onto three complementary descriptions:
        ! 1. constituent site amounts, adjusted by coordination numbers;
        ! 2. equivalent constituent fractions, counting the two positions equally;
        ! 3. cross-sublattice pair populations, counting every A-X pairing.
        !
        ! For pair i-j and quadruplet q, nA and nX are the numbers of times i
        ! and j occur in q. The ordinary incidence is nA*nX. The weighted
        ! incidence is the production Eq. 31 quantity (nA*nX)/zeta(i,j).
        ! Zeta is fixed model data, so only the quadruplet mole amount carries
        ! derivatives; pair-specific SUBQ values change the resulting normalized
        ! weighted distribution without changing the ordinary pair fractions.
        do q = 1, SIZE(dMoles)
            a = tModel%iQuadruplet(q,1); b = tModel%iQuadruplet(q,2)
            x = tModel%iQuadruplet(q,3); y = tModel%iQuadruplet(q,4)
            dSiteAmount1(a) = dSiteAmount1(a) + dMoles(q)/tModel%dCoordination(q,1)
            dSiteAmount1(b) = dSiteAmount1(b) + dMoles(q)/tModel%dCoordination(q,2)
            dSiteAmount2(x) = dSiteAmount2(x) + dMoles(q)/tModel%dCoordination(q,3)
            dSiteAmount2(y) = dSiteAmount2(y) + dMoles(q)/tModel%dCoordination(q,4)
            dEquivalent1(a) = dEquivalent1(a) + 0.5D0*dQuadFraction(q)
            dEquivalent1(b) = dEquivalent1(b) + 0.5D0*dQuadFraction(q)
            dEquivalent2(x) = dEquivalent2(x) + 0.5D0*dQuadFraction(q)
            dEquivalent2(y) = dEquivalent2(y) + 0.5D0*dQuadFraction(q)
            do i = 1, tModel%nSublattice1
                nA = MERGE(1,0,a==i)+MERGE(1,0,b==i)
                do j = 1, tModel%nSublattice2
                    nX = MERGE(1,0,x==j)+MERGE(1,0,y==j)
                    dPairAmount(i,j) = dPairAmount(i,j) + dMoles(q)*DFLOAT(nA*nX)
                    ! SUBG receives the same zeta in every entry, whereas SUBQ can
                    ! assign a distinct fixed value to this exact constituent pair.
                    dPairWeightedAmount(i,j) = dPairWeightedAmount(i,j) + &
                        dMoles(q)*DFLOAT(nA*nX)/tModel%dZeta(i,j)
                end do
            end do
        end do
        ! Normalize the ordinary and zeta-weighted pair amounts separately.
        ! F1 and F2 are marginals of the weighted distribution, so their values
        ! also change when SUBQ pair-specific zeta values differ. Their quotient
        ! and logarithm derivative formulas do not need to be re-derived.
        dSum1 = SUM(dSiteAmount1); dSum2 = SUM(dSiteAmount2)
        dSiteFraction1 = dSiteAmount1/dSum1
        dSiteFraction2 = dSiteAmount2/dSum2
        dPairSum = SUM(dPairAmount); dWeightedSum = SUM(dPairWeightedAmount)
        dPairFraction = dPairAmount/dPairSum
        dPairWeightedFraction = dPairWeightedAmount/dWeightedSum
        do i = 1, tModel%nSublattice1
            dF1(i) = SUM(dPairWeightedFraction(i,:))
        end do
        do j = 1, tModel%nSublattice2
            dF2(j) = SUM(dPairWeightedFraction(:,j))
        end do

    end subroutine AllocateScalarState

    !---------------------------------------------------------------------------------------------------------
    ! SECTION 5A: SCALAR EXCESS-ENERGY FAMILIES
    !
    ! G and Q combine a local composition modifier with a topology-dependent
    ! extensive amount. B uses its proven N_Q prefactor and weighted-pair
    ! composition. Supported ternary factors modify only the traced G/Q branches.
    !---------------------------------------------------------------------------------------------------------


    !---------------------------------------------------------------------------------------------------------
    !> \brief Assemble the nonideal energy contributed by the production interaction families.
    !>
    !> \details Every family supplies an intensive composition-dependent modifier.
    !!          Multiplying by the correct phase or topology amount turns that
    !!          modifier into an extensive Gibbs-energy contribution.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateScalarExcess(tModel,dMoles,dN,dQuadFraction,dPairWeightedFraction, &
        tInteraction,dExcess,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        real(8), intent(in) :: dMoles(:), dN, dQuadFraction(:), dPairWeightedFraction(:,:)
        type(MQMQAInteractionTerm), intent(in) :: tInteraction(:)
        real(8), intent(out) :: dExcess
        integer, intent(inout) :: iInfo

        integer :: i
        real(8) :: dModifier, dOuter

        dExcess = 0D0
        do i = 1, SIZE(tInteraction)
            if (tInteraction(i)%iFamily == MQMQA_TERM_B) then
                ! B is a direct competition between two selected weighted pair
                ! fractions. Its proven extensive prefactor is the total N_Q.
                call ScalarBModifier(dPairWeightedFraction,tInteraction(i),dModifier,iInfo)
                if (iInfo /= 0) return
                ! Multiplying the composition-only modifier by total quadruplet
                ! moles N_Q makes this an extensive energy. Differentiating that
                ! product produces both a direct term and the response of the
                ! zeta-weighted pair fractions.
                dExcess = dExcess + dN*dModifier
            else
                ! G and Q separate into a composition polynomial and a topology
                ! amount representing the quadruplets affected by that interaction.
                call ScalarGQModifier(tModel,dQuadFraction,tInteraction(i),dModifier,iInfo)
                if (iInfo /= 0) return
                call ScalarOuterAmount(tModel,dMoles,tInteraction(i),dOuter,iInfo)
                if (iInfo /= 0) return
                dExcess = dExcess + dOuter*dModifier
            end if
        end do

    end subroutine EvaluateScalarExcess


    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate the B-family competition between two zeta-weighted pair types.
    !>
    !> \details The numerator rewards simultaneous presence of the selected A-X
    !!          and B-Y pairs. The denominator normalizes their combined amount,
    !!          making the modifier composition dependent but scale independent.
    !---------------------------------------------------------------------------------------------------------
    subroutine ScalarBModifier(dPairWeightedFraction,tTerm,dModifier,iInfo)

        real(8), intent(in) :: dPairWeightedFraction(:,:)
        type(MQMQAInteractionTerm), intent(in) :: tTerm
        real(8), intent(out) :: dModifier
        integer, intent(inout) :: iInfo

        real(8) :: dFirst, dSecond, dDen

        dFirst = dPairWeightedFraction(tTerm%iA,tTerm%iX)
        dSecond = dPairWeightedFraction(tTerm%iB,tTerm%iY)
        dDen = dFirst+dSecond
        if ((dFirst <= 0D0) .OR. (dSecond <= 0D0) .OR. (dDen <= 0D0)) then
            iInfo = 50
            dModifier = 0D0
            return
        end if
        dModifier = tTerm%dCoefficient*dFirst**(1+tTerm%iExponentP)*dSecond**(1+tTerm%iExponentQ) / &
            dDen**(1+tTerm%iExponentP+tTerm%iExponentQ)

    end subroutine ScalarBModifier


    !---------------------------------------------------------------------------------------------------------
    !> \brief Evaluate the binary G/Q composition polynomial and optional ternary modifier.
    !>
    !> \details Chi is a local binary coordinate weighted by how many of the two
    !!          fixed-environment positions contain the selected constituent. SUBG
    !!          admits only a fully matching diagonal environment, while SUBQ uses
    !!          half the matching-position count and therefore permits weights 1,
    !!          1/2, and 0. Xi remains a broader projection of asymmetric-group
    !!          occurrence. G-family parameters use chi; Q-family parameters use xi.
    !---------------------------------------------------------------------------------------------------------
    subroutine ScalarGQModifier(tModel,dQuadFraction,tTerm,dModifier,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        real(8), intent(in) :: dQuadFraction(:)
        type(MQMQAInteractionTerm), intent(in) :: tTerm
        real(8), intent(out) :: dModifier
        integer, intent(inout) :: iInfo

        integer :: q, a, b, x, y, nFixed
        real(8) :: dChi1, dChi2, dChiDen, dXi1, dXi2, dXiDen, dTernary
        real(8) :: dEnvironmentWeight

        dChi1 = 0D0; dChi2 = 0D0; dChiDen = 0D0; dXi1 = 0D0; dXi2 = 0D0
        do q = 1, SIZE(dQuadFraction)
            a = tModel%iQuadruplet(q,1); b = tModel%iQuadruplet(q,2)
            x = tModel%iQuadruplet(q,3); y = tModel%iQuadruplet(q,4)
            if (tTerm%iX == tTerm%iY) then
                ! Binary branch with fixed X-X environment: A/B carry the
                ! changing binary composition on the first sublattice. This mirrors
                ! production SUBG/SUBQ selection directly from the two constituent
                ! indices: both matches give one, a single SUBQ match gives one half.
                dEnvironmentWeight = 0D0
                if ((tTerm%iX == x) .AND. (tTerm%iX == y)) then
                    dEnvironmentWeight = 1D0
                else if ((tModel%iModelType == MQMQA_MODEL_SUBQ) .AND. &
                    ((tTerm%iX == x) .OR. (tTerm%iX == y))) then
                    dEnvironmentWeight = 0.5D0
                end if
                if (dEnvironmentWeight > 0D0) then
                    if (tTerm%lGroup1(a) .AND. tTerm%lGroup1(b)) &
                        dChi1 = dChi1+dEnvironmentWeight*dQuadFraction(q)
                    if (tTerm%lGroup2(a) .AND. tTerm%lGroup2(b)) &
                        dChi2 = dChi2+dEnvironmentWeight*dQuadFraction(q)
                    if ((tTerm%lGroup1(a).OR.tTerm%lGroup2(a)) .AND. &
                        (tTerm%lGroup1(b).OR.tTerm%lGroup2(b))) &
                        dChiDen = dChiDen+dEnvironmentWeight*dQuadFraction(q)
                end if
                nFixed = MERGE(1,0,x==tTerm%iX)+MERGE(1,0,y==tTerm%iX)
                if (tTerm%lGroup1(a)) dXi1=dXi1+0.25D0*dQuadFraction(q)*nFixed
                if (tTerm%lGroup1(b)) dXi1=dXi1+0.25D0*dQuadFraction(q)*nFixed
                if (tTerm%lGroup2(a)) dXi2=dXi2+0.25D0*dQuadFraction(q)*nFixed
                if (tTerm%lGroup2(b)) dXi2=dXi2+0.25D0*dQuadFraction(q)*nFixed
            else
                ! Sublattice-swapped binary branch with fixed A-A environment:
                ! X/Y carry the changing composition on the second sublattice.
                dEnvironmentWeight = 0D0
                if ((tTerm%iA == a) .AND. (tTerm%iA == b)) then
                    dEnvironmentWeight = 1D0
                else if ((tModel%iModelType == MQMQA_MODEL_SUBQ) .AND. &
                    ((tTerm%iA == a) .OR. (tTerm%iA == b))) then
                    dEnvironmentWeight = 0.5D0
                end if
                if (dEnvironmentWeight > 0D0) then
                    if (tTerm%lGroup1(x) .AND. tTerm%lGroup1(y)) &
                        dChi1=dChi1+dEnvironmentWeight*dQuadFraction(q)
                    if (tTerm%lGroup2(x) .AND. tTerm%lGroup2(y)) &
                        dChi2=dChi2+dEnvironmentWeight*dQuadFraction(q)
                    if ((tTerm%lGroup1(x).OR.tTerm%lGroup2(x)) .AND. &
                        (tTerm%lGroup1(y).OR.tTerm%lGroup2(y))) &
                        dChiDen=dChiDen+dEnvironmentWeight*dQuadFraction(q)
                end if
                nFixed = MERGE(1,0,a==tTerm%iA)+MERGE(1,0,b==tTerm%iA)
                if (tTerm%lGroup1(x)) dXi1=dXi1+0.25D0*dQuadFraction(q)*nFixed
                if (tTerm%lGroup1(y)) dXi1=dXi1+0.25D0*dQuadFraction(q)*nFixed
                if (tTerm%lGroup2(x)) dXi2=dXi2+0.25D0*dQuadFraction(q)*nFixed
                if (tTerm%lGroup2(y)) dXi2=dXi2+0.25D0*dQuadFraction(q)*nFixed
            end if
        end do
        dXiDen = dXi1+dXi2
        if ((dChi1 <= 0D0) .OR. (dChi2 <= 0D0) .OR. (dChiDen <= 0D0) .OR. &
            (dXi1 <= 0D0) .OR. (dXi2 <= 0D0) .OR. (dXiDen <= 0D0)) then
            iInfo = 51
            dModifier = 0D0
            return
        end if
        dChi1=dChi1/dChiDen; dChi2=dChi2/dChiDen
        call ScalarTernaryFactor(tModel,dQuadFraction,tTerm,dXi1,dXi2,dTernary,iInfo)
        if (iInfo /= 0) return
        if (tTerm%iFamily == MQMQA_TERM_G) then
            ! G parameters use the two local binary shares directly. The p/q
            ! powers determine how strongly the interaction depends on each side.
            dModifier=tTerm%dCoefficient*dChi1**tTerm%iExponentP*dChi2**tTerm%iExponentQ*dTernary
        else
            ! Q parameters use projected xi amounts normalized by their combined
            ! amount, so the polynomial responds to composition rather than scale.
            dModifier=tTerm%dCoefficient*dXi1**tTerm%iExponentP*dXi2**tTerm%iExponentQ / &
                dXiDen**(tTerm%iExponentP+tTerm%iExponentQ)*dTernary
        end if

    end subroutine ScalarGQModifier


    !---------------------------------------------------------------------------------------------------------
    !> \brief Extend a binary G/Q interaction to one supported third constituent.
    !>
    !> \details The third constituent D can belong to asymmetric group 1, group 2,
    !!          or neither. Production SUBG uses a different normalization in each
    !!          case. The factor is one when no ternary constituent is requested,
    !!          so the binary expression is recovered exactly.
    !---------------------------------------------------------------------------------------------------------
    subroutine ScalarTernaryFactor(tModel,dQuadFraction,tTerm,dXi1,dXi2,dFactor,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        real(8), intent(in) :: dQuadFraction(:), dXi1, dXi2
        type(MQMQAInteractionTerm), intent(in) :: tTerm
        real(8), intent(out) :: dFactor
        integer, intent(inout) :: iInfo

        integer :: q, a, b, x, y, nD, nA, nB, nX, d
        real(8) :: dYdk, dYik, dYjk, dBase

        dFactor = 1D0
        d = tTerm%iTernaryConstituent
        if (d <= 0) return
        dYdk=0D0; dYik=0D0; dYjk=0D0
        do q=1,SIZE(dQuadFraction)
            a=tModel%iQuadruplet(q,1); b=tModel%iQuadruplet(q,2)
            x=tModel%iQuadruplet(q,3); y=tModel%iQuadruplet(q,4)
            nX=MERGE(1,0,x==tTerm%iX)+MERGE(1,0,y==tTerm%iX)
            nD=MERGE(1,0,a==d)+MERGE(1,0,b==d)
            nA=MERGE(1,0,a==tTerm%iA)+MERGE(1,0,b==tTerm%iA)
            nB=MERGE(1,0,a==tTerm%iB)+MERGE(1,0,b==tTerm%iB)
            dYdk=dYdk+0.25D0*dQuadFraction(q)*nD*nX
            dYik=dYik+0.25D0*dQuadFraction(q)*nA*nX
            dYjk=dYjk+0.25D0*dQuadFraction(q)*nB*nX
        end do
        if (tTerm%lGroup2(d)) then
            ! D belongs with the second binary group: normalize D by xi2 and
            ! remove the fraction already occupied by the B-like endpoint.
            dBase=1D0-dYjk/dXi2
            if ((dYdk<=0D0).OR.(dBase<=0D0)) then; iInfo=52; return; end if
            dFactor=dYdk/dXi2*dBase**(tTerm%iExponentR-1)
        else if (tTerm%lGroup1(d)) then
            ! D belongs with the first binary group: the same construction uses
            ! xi1 and the A-like endpoint.
            dBase=1D0-dYik/dXi1
            if ((dYdk<=0D0).OR.(dBase<=0D0)) then; iInfo=52; return; end if
            dFactor=dYdk/dXi1*dBase**(tTerm%iExponentR-1)
        else
            ! D belongs to neither asymmetric group, so its available background
            ! is the composition left outside both xi1 and xi2.
            dBase=1D0-dXi1-dXi2
            if ((dYdk<=0D0).OR.(dBase<=0D0)) then; iInfo=52; return; end if
            dFactor=dYdk*dBase**(tTerm%iExponentR-1)
        end if

    end subroutine ScalarTernaryFactor


    !---------------------------------------------------------------------------------------------------------
    !> \brief Convert an intensive G/Q modifier into its affected extensive amount.
    !>
    !> \details The base interaction contributes half of its defining quadruplet.
    !!          Quadruplets containing one substituted constituent also contribute
    !!          through coordination-number ratios. This is the scalar counterpart
    !!          of the production direct and diagonal-cation/anion T1-T3 assembly.
    !---------------------------------------------------------------------------------------------------------
    subroutine ScalarOuterAmount(tModel,dMoles,tTerm,dOuter,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        real(8), intent(in) :: dMoles(:)
        type(MQMQAInteractionTerm), intent(in) :: tTerm
        real(8), intent(out) :: dOuter
        integer, intent(inout) :: iInfo

        integer :: qBase, qOther, c, z, iPosition

        qBase=FindQuadruplet(tModel,tTerm%iA,tTerm%iB,tTerm%iX,tTerm%iY)
        if (qBase<=0) then; iInfo=53; dOuter=0D0; return; end if
        dOuter=0.5D0*dMoles(qBase)
        if ((tTerm%iA==tTerm%iB).AND.(tTerm%iX/=tTerm%iY)) then
            ! A=A is fixed: gather all first-sublattice substitutions coupled to X-Y.
            do c=1,tModel%nSublattice1
                if (c==tTerm%iA) cycle
                qOther=FindQuadruplet(tModel,MIN(tTerm%iA,c),MAX(tTerm%iA,c),tTerm%iX,tTerm%iY)
                if (qOther<=0) then; iInfo=53; return; end if
                iPosition=MERGE(1,2,tTerm%iA<c)
                dOuter=dOuter+0.25D0*dMoles(qOther)*tModel%dCoordination(qBase,1)/ &
                    tModel%dCoordination(qOther,iPosition)
            end do
        else if ((tTerm%iA/=tTerm%iB).AND.(tTerm%iX==tTerm%iY)) then
            ! X=X is fixed: gather the second-sublattice substitutions coupled to A-B.
            do z=1,tModel%nSublattice2
                if (z==tTerm%iX) cycle
                qOther=FindQuadruplet(tModel,tTerm%iA,tTerm%iB,MIN(tTerm%iX,z),MAX(tTerm%iX,z))
                if (qOther<=0) then; iInfo=53; return; end if
                iPosition=MERGE(3,4,tTerm%iX<z)
                dOuter=dOuter+0.25D0*dMoles(qOther)*tModel%dCoordination(qBase,3)/ &
                    tModel%dCoordination(qOther,iPosition)
            end do
        end if

    end subroutine ScalarOuterAmount


    !> Locate the canonical [A,B,X,Y] row needed by an interaction prefactor.
    integer function FindQuadruplet(tModel,a,b,x,y)

        type(MQMQAModelData), intent(in) :: tModel
        integer, intent(in) :: a,b,x,y
        integer :: q

        FindQuadruplet=0
        do q=1,SIZE(tModel%iQuadruplet,1)
            if (ALL(tModel%iQuadruplet(q,:)==[a,b,x,y])) then
                FindQuadruplet=q
                return
            end if
        end do

    end function FindQuadruplet

    !=========================================================================================================
    ! SECTION 6: ANALYTIC DERIVATIVE-OBJECT EVALUATOR
    !
    ! Rebuild every dependent composition measure and every energy block with
    ! SecondOrderScalar objects. Values, gradients, and Hessians then propagate
    ! together from independent quadruplet moles through all normalization,
    ! ordering, and interaction formulas.
    !=========================================================================================================


    !---------------------------------------------------------------------------------------------------------
    !> \brief Rebuild the complete scalar model with derivative-carrying quantities.
    !>
    !> \details This follows the same thermodynamic decomposition as the ordinary-real
    !!          path, but every intermediate quantity carries its value, gradient,
    !!          and Hessian. The chain rule therefore passes from quadruplet moles
    !!          through all normalizations into the final energy automatically.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateDerivativeEnergy(tModel,dMoles,dIdealScale,tInteraction,tReference,tIdeal,tExcess,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        real(8), intent(in) :: dMoles(:), dIdealScale
        type(MQMQAInteractionTerm), intent(in) :: tInteraction(:)
        type(SecondOrderScalar), intent(out) :: tReference, tIdeal, tExcess
        integer, intent(out) :: iInfo

        integer :: nQuad, q, a, b, x, y, i, j, nA, nX, iPosition, jPosition, iWeight
        real(8) :: dTheta, dPsi
        type(SecondOrderScalar) :: tN, tS1, tS2, tS3, tDen, tTermValue
        type(SecondOrderScalar) :: tPairLogBlock, tEquivalentLogBlock, tLogRatio
        type(SecondOrderScalar) :: tMoles(SIZE(dMoles)), tQuadFraction(SIZE(dMoles))
        type(SecondOrderScalar) :: tSiteAmount1(tModel%nSublattice1), tSiteAmount2(tModel%nSublattice2)
        type(SecondOrderScalar) :: tSiteFraction1(tModel%nSublattice1), tSiteFraction2(tModel%nSublattice2)
        type(SecondOrderScalar) :: tEquivalent1(tModel%nSublattice1)
        type(SecondOrderScalar) :: tEquivalent2(tModel%nSublattice2)
        type(SecondOrderScalar) :: tPairAmount(tModel%nSublattice1,tModel%nSublattice2)
        type(SecondOrderScalar) :: tPairWeightedAmount(tModel%nSublattice1,tModel%nSublattice2)
        type(SecondOrderScalar) :: tPairFraction(tModel%nSublattice1,tModel%nSublattice2)
        type(SecondOrderScalar) :: tPairWeightedFraction(tModel%nSublattice1,tModel%nSublattice2)
        type(SecondOrderScalar) :: tF1(tModel%nSublattice1), tF2(tModel%nSublattice2)
        type(SecondOrderScalar) :: tSiteSum1, tSiteSum2, tPairSum, tWeightedSum

        call CheckInputs(tModel,dMoles,tInteraction,iInfo)
        if (iInfo /= 0) return
        nQuad=SIZE(dMoles)

        ! Seed one independent derivative variable for every quadruplet mole.
        ! Reference energy is linear, while division by total N_Q creates the
        ! composition derivatives of every quadruplet fraction.
        tN=ConstantSO(0D0,nQuad)
        tReference=ConstantSO(0D0,nQuad)
        do q=1,nQuad
            tMoles(q)=VariableSO(dMoles(q),q,nQuad)
            tN=AddSO(tN,tMoles(q))
            tReference=AddSO(tReference,ScaleSO(tMoles(q),tModel%dReferenceEnergy(q)))
        end do
        do q=1,nQuad
            tQuadFraction(q)=DivideSO(tMoles(q),tN)
        end do
        ! Initialize the dependent site, equivalent-constituent, pair, and
        ! pair-marginal quantities before accumulating quadruplet contributions.
        do a=1,tModel%nSublattice1
            tSiteAmount1(a)=ConstantSO(0D0,nQuad)
            tEquivalent1(a)=ConstantSO(0D0,nQuad)
            tF1(a)=ConstantSO(0D0,nQuad)
        end do
        do x=1,tModel%nSublattice2
            tSiteAmount2(x)=ConstantSO(0D0,nQuad)
            tEquivalent2(x)=ConstantSO(0D0,nQuad)
            tF2(x)=ConstantSO(0D0,nQuad)
        end do
        do a=1,tModel%nSublattice1
            do x=1,tModel%nSublattice2
                tPairAmount(a,x)=ConstantSO(0D0,nQuad)
                tPairWeightedAmount(a,x)=ConstantSO(0D0,nQuad)
            end do
        end do

        ! This is the derivative-carrying counterpart of AllocateScalarState.
        ! Because each tMoles(q) is seeded independently, these projections also
        ! construct their complete first- and second-order composition response.
        ! The coefficient (nA*nX)/zeta(i,j) is constant with respect to moles:
        ! SUBQ changes its pair-specific numerical value, not the outer product,
        ! quotient, logarithm, or Hessian propagation rules.
        do q=1,nQuad
            a=tModel%iQuadruplet(q,1); b=tModel%iQuadruplet(q,2)
            x=tModel%iQuadruplet(q,3); y=tModel%iQuadruplet(q,4)
            tSiteAmount1(a)=AddSO(tSiteAmount1(a),ScaleSO(tMoles(q),1D0/tModel%dCoordination(q,1)))
            tSiteAmount1(b)=AddSO(tSiteAmount1(b),ScaleSO(tMoles(q),1D0/tModel%dCoordination(q,2)))
            tSiteAmount2(x)=AddSO(tSiteAmount2(x),ScaleSO(tMoles(q),1D0/tModel%dCoordination(q,3)))
            tSiteAmount2(y)=AddSO(tSiteAmount2(y),ScaleSO(tMoles(q),1D0/tModel%dCoordination(q,4)))
            tEquivalent1(a)=AddSO(tEquivalent1(a),ScaleSO(tQuadFraction(q),0.5D0))
            tEquivalent1(b)=AddSO(tEquivalent1(b),ScaleSO(tQuadFraction(q),0.5D0))
            tEquivalent2(x)=AddSO(tEquivalent2(x),ScaleSO(tQuadFraction(q),0.5D0))
            tEquivalent2(y)=AddSO(tEquivalent2(y),ScaleSO(tQuadFraction(q),0.5D0))
            do i=1,tModel%nSublattice1
                nA=MERGE(1,0,a==i)+MERGE(1,0,b==i)
                do j=1,tModel%nSublattice2
                    nX=MERGE(1,0,x==j)+MERGE(1,0,y==j)
                    tPairAmount(i,j)=AddSO(tPairAmount(i,j),ScaleSO(tMoles(q),DFLOAT(nA*nX)))
                    tPairWeightedAmount(i,j)=AddSO(tPairWeightedAmount(i,j), &
                        ScaleSO(tMoles(q),DFLOAT(nA*nX)/tModel%dZeta(i,j)))
                end do
            end do
        end do

        tSiteSum1=ConstantSO(0D0,nQuad); tSiteSum2=ConstantSO(0D0,nQuad)
        ! Normalize the weighted amounts and form their first- and second-
        ! sublattice marginals. These F values feed S2; the complete weighted
        ! pair distribution also feeds B and the SUBQ form of S3. SUBG S3 uses
        ! the separately normalized ordinary pair distribution.
        do a=1,tModel%nSublattice1
            tSiteSum1=AddSO(tSiteSum1,tSiteAmount1(a))
        end do
        do x=1,tModel%nSublattice2
            tSiteSum2=AddSO(tSiteSum2,tSiteAmount2(x))
        end do
        do a=1,tModel%nSublattice1
            tSiteFraction1(a)=DivideSO(tSiteAmount1(a),tSiteSum1)
        end do
        do x=1,tModel%nSublattice2
            tSiteFraction2(x)=DivideSO(tSiteAmount2(x),tSiteSum2)
        end do

        tPairSum=ConstantSO(0D0,nQuad); tWeightedSum=ConstantSO(0D0,nQuad)
        do a=1,tModel%nSublattice1
            do x=1,tModel%nSublattice2
                tPairSum=AddSO(tPairSum,tPairAmount(a,x))
                tWeightedSum=AddSO(tWeightedSum,tPairWeightedAmount(a,x))
            end do
        end do
        do a=1,tModel%nSublattice1
            do x=1,tModel%nSublattice2
                tPairFraction(a,x)=DivideSO(tPairAmount(a,x),tPairSum)
                tPairWeightedFraction(a,x)=DivideSO(tPairWeightedAmount(a,x),tWeightedSum)
                tF1(a)=AddSO(tF1(a),tPairWeightedFraction(a,x))
                tF2(x)=AddSO(tF2(x),tPairWeightedFraction(a,x))
            end do
        end do

        ! S1: individual constituent mixing on both sublattices.
        tS1=ConstantSO(0D0,nQuad)
        do a=1,tModel%nSublattice1
            if (tSiteFraction1(a)%dValue<=0D0) then; iInfo=40; return; end if
            tS1=AddSO(tS1,MultiplySO(tSiteAmount1(a),LogSO(tSiteFraction1(a))))
        end do
        do x=1,tModel%nSublattice2
            if (tSiteFraction2(x)%dValue<=0D0) then; iInfo=40; return; end if
            tS1=AddSO(tS1,MultiplySO(tSiteAmount2(x),LogSO(tSiteFraction2(x))))
        end do

        ! S2: pair-ordering correction relative to independent pair marginals.
        tS2=ConstantSO(0D0,nQuad)
        do a=1,tModel%nSublattice1
            do x=1,tModel%nSublattice2
                tDen=MultiplySO(tF1(a),tF2(x))
                if ((tPairWeightedFraction(a,x)%dValue<=0D0).OR.(tDen%dValue<=0D0)) then
                    iInfo=41
                    return
                end if
                tTermValue=MultiplySO(tPairWeightedAmount(a,x), &
                    LogSO(DivideSO(tPairWeightedFraction(a,x),tDen)))
                tS2=AddSO(tS2,tTermValue)
            end do
        end do

        ! S3: complete quadruplet-ordering correction relative to the pair model.
        ! The derivative objects below retain every logarithm explicitly. SUBQ changes
        ! only the fixed weights multiplying the pair and equivalent-fraction blocks.
        dTheta=0D0
        dPsi=0D0
        select case(tModel%iModelType)
        case(MQMQA_MODEL_SUBG)
            dTheta=1D0
            dPsi=1D0
        case(MQMQA_MODEL_SUBQ)
            dTheta=3D0/4D0
            dPsi=1D0/2D0
        end select
        tS3=ConstantSO(0D0,nQuad)
        do q=1,nQuad
            iWeight=1
            if (tModel%iQuadruplet(q,1)/=tModel%iQuadruplet(q,2)) iWeight=2*iWeight
            if (tModel%iQuadruplet(q,3)/=tModel%iQuadruplet(q,4)) iWeight=2*iWeight
            tPairLogBlock=ConstantSO(0D0,nQuad)
            do iPosition=1,2
                do jPosition=3,4
                    if ((tModel%iModelType==MQMQA_MODEL_SUBG).OR. &
                        ((tModel%iModelType==MQMQA_MODEL_SUBQ).AND.lMQMQADiagnosticLegacyS3)) then
                        tDen=tPairFraction(tModel%iQuadruplet(q,iPosition), &
                            tModel%iQuadruplet(q,jPosition))
                    else
                        ! Keep the derivative structure unchanged while selecting
                        ! the weighted SUBQ pair distribution used by the scalar path.
                        tDen=tPairWeightedFraction(tModel%iQuadruplet(q,iPosition), &
                            tModel%iQuadruplet(q,jPosition))
                    end if
                    if (tDen%dValue<=0D0) then; iInfo=42; return; end if
                    tPairLogBlock=AddSO(tPairLogBlock,LogSO(tDen))
                end do
            end do
            if ((tQuadFraction(q)%dValue<=0D0).OR. &
                (tEquivalent1(tModel%iQuadruplet(q,1))%dValue<=0D0).OR. &
                (tEquivalent1(tModel%iQuadruplet(q,2))%dValue<=0D0).OR. &
                (tEquivalent2(tModel%iQuadruplet(q,3))%dValue<=0D0).OR. &
                (tEquivalent2(tModel%iQuadruplet(q,4))%dValue<=0D0)) then
                iInfo=42
                return
            end if
            tEquivalentLogBlock=AddSO(AddSO(LogSO(tEquivalent1(tModel%iQuadruplet(q,1))), &
                LogSO(tEquivalent1(tModel%iQuadruplet(q,2)))), &
                AddSO(LogSO(tEquivalent2(tModel%iQuadruplet(q,3))), &
                LogSO(tEquivalent2(tModel%iQuadruplet(q,4)))))
            tLogRatio=SubtractSO(LogSO(tQuadFraction(q)),ConstantSO(DLOG(DFLOAT(iWeight)),nQuad))
            tLogRatio=SubtractSO(tLogRatio,ScaleSO(tPairLogBlock,dTheta))
            tLogRatio=AddSO(tLogRatio,ScaleSO(tEquivalentLogBlock,dPsi))
            tS3=AddSO(tS3,MultiplySO(tMoles(q),tLogRatio))
        end do
        tIdeal=ScaleSO(AddSO(AddSO(tS1,tS2),tS3),dIdealScale)

        ! Add nonideal G, Q, B, and supported ternary curvature after the
        ! configurational terms have established the common composition state.
        call EvaluateDerivativeExcess(tModel,tMoles,tN,tQuadFraction,tPairWeightedFraction, &
            tInteraction,tExcess,iInfo)

    end subroutine EvaluateDerivativeEnergy

    !---------------------------------------------------------------------------------------------------------
    ! SECTION 6A: DERIVATIVE-CARRYING EXCESS-ENERGY FAMILIES
    !
    ! These routines mirror the physical decomposition of Section 5A while using
    ! second-order objects. They are not called by the independent scalar path.
    !---------------------------------------------------------------------------------------------------------


    !---------------------------------------------------------------------------------------------------------
    !> \brief Assemble derivative-carrying excess terms using the proven extensive prefactors.
    !---------------------------------------------------------------------------------------------------------
    subroutine EvaluateDerivativeExcess(tModel,tMoles,tN,tQuadFraction,tPairWeightedFraction, &
        tInteraction,tExcess,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        type(SecondOrderScalar), intent(in) :: tMoles(:), tN, tQuadFraction(:), tPairWeightedFraction(:,:)
        type(MQMQAInteractionTerm), intent(in) :: tInteraction(:)
        type(SecondOrderScalar), intent(out) :: tExcess
        integer, intent(inout) :: iInfo

        integer :: i, nQuad
        type(SecondOrderScalar) :: tModifier, tOuter

        nQuad=SIZE(tMoles)
        tExcess=ConstantSO(0D0,nQuad)
        do i=1,SIZE(tInteraction)
            if (tInteraction(i)%iFamily==MQMQA_TERM_B) then
                ! Differentiating N_Q*Delta g_B generates both the direct B
                ! contribution and the zeta-weighted composition-response terms.
                call DerivativeBModifier(tPairWeightedFraction,tInteraction(i),tModifier,iInfo)
                if (iInfo/=0) return
                tExcess=AddSO(tExcess,MultiplySO(tN,tModifier))
            else
                ! G/Q curvature comes from both factors: the intensive modifier
                ! and the extensive amount of topology affected by the parameter.
                call DerivativeGQModifier(tModel,tQuadFraction,tInteraction(i),tModifier,iInfo)
                if (iInfo/=0) return
                call DerivativeOuterAmount(tModel,tMoles,tInteraction(i),tOuter,iInfo)
                if (iInfo/=0) return
                tExcess=AddSO(tExcess,MultiplySO(tOuter,tModifier))
            end if
        end do

    end subroutine EvaluateDerivativeExcess


    !---------------------------------------------------------------------------------------------------------
    !> \brief Derivative-object form of the scale-independent B pair competition.
    !---------------------------------------------------------------------------------------------------------
    subroutine DerivativeBModifier(tPairWeightedFraction,tTerm,tModifier,iInfo)

        type(SecondOrderScalar), intent(in) :: tPairWeightedFraction(:,:)
        type(MQMQAInteractionTerm), intent(in) :: tTerm
        type(SecondOrderScalar), intent(out) :: tModifier
        integer, intent(inout) :: iInfo

        type(SecondOrderScalar) :: tFirst, tSecond, tDen

        tFirst=tPairWeightedFraction(tTerm%iA,tTerm%iX)
        tSecond=tPairWeightedFraction(tTerm%iB,tTerm%iY)
        tDen=AddSO(tFirst,tSecond)
        if ((tFirst%dValue<=0D0).OR.(tSecond%dValue<=0D0).OR.(tDen%dValue<=0D0)) then
            iInfo=50
            return
        end if
        tModifier=ScaleSO(DivideSO(MultiplySO(PowerSO(tFirst,1+tTerm%iExponentP), &
            PowerSO(tSecond,1+tTerm%iExponentQ)),PowerSO(tDen,1+tTerm%iExponentP+tTerm%iExponentQ)), &
            tTerm%dCoefficient)

    end subroutine DerivativeBModifier


    !---------------------------------------------------------------------------------------------------------
    !> \brief Derivative-object form of the binary G/Q and optional ternary modifier.
    !>
    !> \details The branch logic intentionally mirrors ScalarGQModifier, while the
    !!          arithmetic is independently expressed through SecondOrderScalar.
    !---------------------------------------------------------------------------------------------------------
    subroutine DerivativeGQModifier(tModel,tQuadFraction,tTerm,tModifier,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        type(SecondOrderScalar), intent(in) :: tQuadFraction(:)
        type(MQMQAInteractionTerm), intent(in) :: tTerm
        type(SecondOrderScalar), intent(out) :: tModifier
        integer, intent(inout) :: iInfo

        integer :: q,a,b,x,y,nFixed,nQuad
        real(8) :: dEnvironmentWeight
        type(SecondOrderScalar) :: tChi1,tChi2,tChiDen,tXi1,tXi2,tXiDen,tTernary

        nQuad=SIZE(tQuadFraction)
        tChi1=ConstantSO(0D0,nQuad); tChi2=ConstantSO(0D0,nQuad); tChiDen=ConstantSO(0D0,nQuad)
        tXi1=ConstantSO(0D0,nQuad); tXi2=ConstantSO(0D0,nQuad)
        do q=1,nQuad
            a=tModel%iQuadruplet(q,1); b=tModel%iQuadruplet(q,2)
            x=tModel%iQuadruplet(q,3); y=tModel%iQuadruplet(q,4)
            if (tTerm%iX==tTerm%iY) then
                ! Fixed X-X environment: only the incidence coefficient differs
                ! between SUBG and SUBQ; DivideSO retains the quotient derivatives.
                dEnvironmentWeight=0D0
                if ((tTerm%iX==x).AND.(tTerm%iX==y)) then
                    dEnvironmentWeight=1D0
                else if ((tModel%iModelType==MQMQA_MODEL_SUBQ).AND. &
                    ((tTerm%iX==x).OR.(tTerm%iX==y))) then
                    dEnvironmentWeight=0.5D0
                end if
                if (dEnvironmentWeight>0D0) then
                    if (tTerm%lGroup1(a).AND.tTerm%lGroup1(b)) &
                        tChi1=AddSO(tChi1,ScaleSO(tQuadFraction(q),dEnvironmentWeight))
                    if (tTerm%lGroup2(a).AND.tTerm%lGroup2(b)) &
                        tChi2=AddSO(tChi2,ScaleSO(tQuadFraction(q),dEnvironmentWeight))
                    if ((tTerm%lGroup1(a).OR.tTerm%lGroup2(a)).AND. &
                        (tTerm%lGroup1(b).OR.tTerm%lGroup2(b))) &
                        tChiDen=AddSO(tChiDen,ScaleSO(tQuadFraction(q),dEnvironmentWeight))
                end if
                nFixed=MERGE(1,0,x==tTerm%iX)+MERGE(1,0,y==tTerm%iX)
                if (tTerm%lGroup1(a)) tXi1=AddSO(tXi1,ScaleSO(tQuadFraction(q),0.25D0*nFixed))
                if (tTerm%lGroup1(b)) tXi1=AddSO(tXi1,ScaleSO(tQuadFraction(q),0.25D0*nFixed))
                if (tTerm%lGroup2(a)) tXi2=AddSO(tXi2,ScaleSO(tQuadFraction(q),0.25D0*nFixed))
                if (tTerm%lGroup2(b)) tXi2=AddSO(tXi2,ScaleSO(tQuadFraction(q),0.25D0*nFixed))
            else
                ! Fixed A-A environment: differentiate the sublattice-swapped X/Y coordinates.
                dEnvironmentWeight=0D0
                if ((tTerm%iA==a).AND.(tTerm%iA==b)) then
                    dEnvironmentWeight=1D0
                else if ((tModel%iModelType==MQMQA_MODEL_SUBQ).AND. &
                    ((tTerm%iA==a).OR.(tTerm%iA==b))) then
                    dEnvironmentWeight=0.5D0
                end if
                if (dEnvironmentWeight>0D0) then
                    if (tTerm%lGroup1(x).AND.tTerm%lGroup1(y)) &
                        tChi1=AddSO(tChi1,ScaleSO(tQuadFraction(q),dEnvironmentWeight))
                    if (tTerm%lGroup2(x).AND.tTerm%lGroup2(y)) &
                        tChi2=AddSO(tChi2,ScaleSO(tQuadFraction(q),dEnvironmentWeight))
                    if ((tTerm%lGroup1(x).OR.tTerm%lGroup2(x)).AND. &
                        (tTerm%lGroup1(y).OR.tTerm%lGroup2(y))) &
                        tChiDen=AddSO(tChiDen,ScaleSO(tQuadFraction(q),dEnvironmentWeight))
                end if
                nFixed=MERGE(1,0,a==tTerm%iA)+MERGE(1,0,b==tTerm%iA)
                if (tTerm%lGroup1(x)) tXi1=AddSO(tXi1,ScaleSO(tQuadFraction(q),0.25D0*nFixed))
                if (tTerm%lGroup1(y)) tXi1=AddSO(tXi1,ScaleSO(tQuadFraction(q),0.25D0*nFixed))
                if (tTerm%lGroup2(x)) tXi2=AddSO(tXi2,ScaleSO(tQuadFraction(q),0.25D0*nFixed))
                if (tTerm%lGroup2(y)) tXi2=AddSO(tXi2,ScaleSO(tQuadFraction(q),0.25D0*nFixed))
            end if
        end do
        tXiDen=AddSO(tXi1,tXi2)
        if ((tChi1%dValue<=0D0).OR.(tChi2%dValue<=0D0).OR.(tChiDen%dValue<=0D0).OR. &
            (tXi1%dValue<=0D0).OR.(tXi2%dValue<=0D0).OR.(tXiDen%dValue<=0D0)) then
            iInfo=51
            return
        end if
        tChi1=DivideSO(tChi1,tChiDen); tChi2=DivideSO(tChi2,tChiDen)
        call DerivativeTernaryFactor(tModel,tQuadFraction,tTerm,tXi1,tXi2,tTernary,iInfo)
        if (iInfo/=0) return
        if (tTerm%iFamily==MQMQA_TERM_G) then
            tModifier=ScaleSO(MultiplySO(MultiplySO(PowerSO(tChi1,tTerm%iExponentP), &
                PowerSO(tChi2,tTerm%iExponentQ)),tTernary),tTerm%dCoefficient)
        else
            tModifier=ScaleSO(MultiplySO(DivideSO(MultiplySO(PowerSO(tXi1,tTerm%iExponentP), &
                PowerSO(tXi2,tTerm%iExponentQ)),PowerSO(tXiDen,tTerm%iExponentP+tTerm%iExponentQ)), &
                tTernary),tTerm%dCoefficient)
        end if

    end subroutine DerivativeGQModifier


    !---------------------------------------------------------------------------------------------------------
    !> \brief Propagate ternary curvature through the group-1, group-2, or neither branch.
    !---------------------------------------------------------------------------------------------------------
    subroutine DerivativeTernaryFactor(tModel,tQuadFraction,tTerm,tXi1,tXi2,tFactor,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        type(SecondOrderScalar), intent(in) :: tQuadFraction(:),tXi1,tXi2
        type(MQMQAInteractionTerm), intent(in) :: tTerm
        type(SecondOrderScalar), intent(out) :: tFactor
        integer, intent(inout) :: iInfo

        integer :: q,a,b,x,y,nD,nA,nB,nX,d,nQuad
        type(SecondOrderScalar) :: tYdk,tYik,tYjk,tBase

        nQuad=SIZE(tQuadFraction)
        tFactor=ConstantSO(1D0,nQuad)
        d=tTerm%iTernaryConstituent
        if (d<=0) return
        tYdk=ConstantSO(0D0,nQuad); tYik=ConstantSO(0D0,nQuad); tYjk=ConstantSO(0D0,nQuad)
        do q=1,nQuad
            a=tModel%iQuadruplet(q,1); b=tModel%iQuadruplet(q,2)
            x=tModel%iQuadruplet(q,3); y=tModel%iQuadruplet(q,4)
            nX=MERGE(1,0,x==tTerm%iX)+MERGE(1,0,y==tTerm%iX)
            nD=MERGE(1,0,a==d)+MERGE(1,0,b==d)
            nA=MERGE(1,0,a==tTerm%iA)+MERGE(1,0,b==tTerm%iA)
            nB=MERGE(1,0,a==tTerm%iB)+MERGE(1,0,b==tTerm%iB)
            tYdk=AddSO(tYdk,ScaleSO(tQuadFraction(q),0.25D0*nD*nX))
            tYik=AddSO(tYik,ScaleSO(tQuadFraction(q),0.25D0*nA*nX))
            tYjk=AddSO(tYjk,ScaleSO(tQuadFraction(q),0.25D0*nB*nX))
        end do
        if (tTerm%lGroup2(d)) then
            ! Third constituent shares the second asymmetric group.
            tBase=SubtractSO(ConstantSO(1D0,nQuad),DivideSO(tYjk,tXi2))
            if ((tYdk%dValue<=0D0).OR.(tBase%dValue<=0D0)) then; iInfo=52; return; end if
            tFactor=MultiplySO(DivideSO(tYdk,tXi2),PowerSO(tBase,tTerm%iExponentR-1))
        else if (tTerm%lGroup1(d)) then
            ! Third constituent shares the first asymmetric group.
            tBase=SubtractSO(ConstantSO(1D0,nQuad),DivideSO(tYik,tXi1))
            if ((tYdk%dValue<=0D0).OR.(tBase%dValue<=0D0)) then; iInfo=52; return; end if
            tFactor=MultiplySO(DivideSO(tYdk,tXi1),PowerSO(tBase,tTerm%iExponentR-1))
        else
            ! Third constituent lies outside both groups.
            tBase=SubtractSO(SubtractSO(ConstantSO(1D0,nQuad),tXi1),tXi2)
            if ((tYdk%dValue<=0D0).OR.(tBase%dValue<=0D0)) then; iInfo=52; return; end if
            tFactor=MultiplySO(tYdk,PowerSO(tBase,tTerm%iExponentR-1))
        end if

    end subroutine DerivativeTernaryFactor


    !---------------------------------------------------------------------------------------------------------
    !> \brief Differentiate the topology amount that makes each G/Q term extensive.
    !---------------------------------------------------------------------------------------------------------
    subroutine DerivativeOuterAmount(tModel,tMoles,tTerm,tOuter,iInfo)

        type(MQMQAModelData), intent(in) :: tModel
        type(SecondOrderScalar), intent(in) :: tMoles(:)
        type(MQMQAInteractionTerm), intent(in) :: tTerm
        type(SecondOrderScalar), intent(out) :: tOuter
        integer, intent(inout) :: iInfo

        integer :: qBase,qOther,c,z,iPosition

        qBase=FindQuadruplet(tModel,tTerm%iA,tTerm%iB,tTerm%iX,tTerm%iY)
        if (qBase<=0) then; iInfo=53; return; end if
        tOuter=ScaleSO(tMoles(qBase),0.5D0)
        if ((tTerm%iA==tTerm%iB).AND.(tTerm%iX/=tTerm%iY)) then
            do c=1,tModel%nSublattice1
                if (c==tTerm%iA) cycle
                qOther=FindQuadruplet(tModel,MIN(tTerm%iA,c),MAX(tTerm%iA,c),tTerm%iX,tTerm%iY)
                if (qOther<=0) then; iInfo=53; return; end if
                iPosition=MERGE(1,2,tTerm%iA<c)
                tOuter=AddSO(tOuter,ScaleSO(tMoles(qOther),0.25D0*tModel%dCoordination(qBase,1)/ &
                    tModel%dCoordination(qOther,iPosition)))
            end do
        else if ((tTerm%iA/=tTerm%iB).AND.(tTerm%iX==tTerm%iY)) then
            do z=1,tModel%nSublattice2
                if (z==tTerm%iX) cycle
                qOther=FindQuadruplet(tModel,tTerm%iA,tTerm%iB,MIN(tTerm%iX,z),MAX(tTerm%iX,z))
                if (qOther<=0) then; iInfo=53; return; end if
                iPosition=MERGE(3,4,tTerm%iX<z)
                tOuter=AddSO(tOuter,ScaleSO(tMoles(qOther),0.25D0*tModel%dCoordination(qBase,3)/ &
                    tModel%dCoordination(qOther,iPosition)))
            end do
        end if

    end subroutine DerivativeOuterAmount


    !=========================================================================================================
    ! SECTION 7: PRIVATE SECOND-ORDER ARITHMETIC
    !
    ! These primitives are the chain-rule engine used only by Section 6. Each
    ! operation returns the value of an expression together with how that
    ! expression changes under every independent quadruplet-mole perturbation.
    ! Their independent verification lives in Section 2.
    !=========================================================================================================

    !> Create a number that has no dependence on any quadruplet mole.
    function ConstantSO(dValue,n) result(tResult)

        real(8), intent(in) :: dValue
        integer, intent(in) :: n
        type(SecondOrderScalar) :: tResult

        tResult%dValue=dValue
        allocate(tResult%dGradient(n),tResult%dHessian(n,n))
        tResult%dGradient=0D0
        tResult%dHessian=0D0

    end function ConstantSO


    !> Seed one independent quadruplet mole with unit first derivative.
    function VariableSO(dValue,iVariable,n) result(tResult)

        real(8), intent(in) :: dValue
        integer, intent(in) :: iVariable,n
        type(SecondOrderScalar) :: tResult

        tResult=ConstantSO(dValue,n)
        tResult%dGradient(iVariable)=1D0

    end function VariableSO


    !> Add two expressions; their sensitivities and curvatures add componentwise.
    function AddSO(tA,tB) result(tResult)

        type(SecondOrderScalar), intent(in) :: tA,tB
        type(SecondOrderScalar) :: tResult

        tResult=ConstantSO(tA%dValue+tB%dValue,SIZE(tA%dGradient))
        tResult%dGradient=tA%dGradient+tB%dGradient
        tResult%dHessian=tA%dHessian+tB%dHessian

    end function AddSO


    function SubtractSO(tA,tB) result(tResult)

        type(SecondOrderScalar), intent(in) :: tA,tB
        type(SecondOrderScalar) :: tResult

        tResult=AddSO(tA,ScaleSO(tB,-1D0))

    end function SubtractSO


    function NegateSO(tA) result(tResult)

        type(SecondOrderScalar), intent(in) :: tA
        type(SecondOrderScalar) :: tResult

        tResult=ScaleSO(tA,-1D0)

    end function NegateSO


    function ScaleSO(tA,dScale) result(tResult)

        type(SecondOrderScalar), intent(in) :: tA
        real(8), intent(in) :: dScale
        type(SecondOrderScalar) :: tResult

        tResult=ConstantSO(dScale*tA%dValue,SIZE(tA%dGradient))
        tResult%dGradient=dScale*tA%dGradient
        tResult%dHessian=dScale*tA%dHessian

    end function ScaleSO


    !> Apply the product rule through second order.
    !>
    !> The two gradient outer products are the cross-curvature created when both
    !> factors change under the same mole perturbations.
    function MultiplySO(tA,tB) result(tResult)

        type(SecondOrderScalar), intent(in) :: tA,tB
        type(SecondOrderScalar) :: tResult
        integer :: i,j,n

        n=SIZE(tA%dGradient)
        tResult=ConstantSO(tA%dValue*tB%dValue,n)
        tResult%dGradient=tA%dGradient*tB%dValue+tB%dGradient*tA%dValue
        tResult%dHessian=tA%dHessian*tB%dValue+tB%dHessian*tA%dValue
        do i=1,n
            do j=1,n
                tResult%dHessian(i,j)=tResult%dHessian(i,j)+ &
                    tA%dGradient(i)*tB%dGradient(j)+tB%dGradient(i)*tA%dGradient(j)
            end do
        end do

    end function MultiplySO


    !> Apply the chain rule to 1/A, including curvature from the changing denominator.
    function ReciprocalSO(tA) result(tResult)

        type(SecondOrderScalar), intent(in) :: tA
        type(SecondOrderScalar) :: tResult
        integer :: i,j,n
        real(8) :: dFirst,dSecond

        n=SIZE(tA%dGradient)
        dFirst=-1D0/(tA%dValue*tA%dValue)
        dSecond=2D0/(tA%dValue*tA%dValue*tA%dValue)
        tResult=ConstantSO(1D0/tA%dValue,n)
        tResult%dGradient=dFirst*tA%dGradient
        tResult%dHessian=dFirst*tA%dHessian
        do i=1,n
            do j=1,n
                tResult%dHessian(i,j)=tResult%dHessian(i,j)+dSecond*tA%dGradient(i)*tA%dGradient(j)
            end do
        end do

    end function ReciprocalSO


    !> Express division as multiplication by a reciprocal so one rule is authoritative.
    function DivideSO(tA,tB) result(tResult)

        type(SecondOrderScalar), intent(in) :: tA,tB
        type(SecondOrderScalar) :: tResult

        tResult=MultiplySO(tA,ReciprocalSO(tB))

    end function DivideSO


    !> Apply the logarithm chain rule used throughout configurational entropy.
    function LogSO(tA) result(tResult)

        type(SecondOrderScalar), intent(in) :: tA
        type(SecondOrderScalar) :: tResult
        integer :: i,j,n
        real(8) :: dFirst,dSecond

        n=SIZE(tA%dGradient)
        dFirst=1D0/tA%dValue
        dSecond=-1D0/(tA%dValue*tA%dValue)
        tResult=ConstantSO(DLOG(tA%dValue),n)
        tResult%dGradient=dFirst*tA%dGradient
        tResult%dHessian=dFirst*tA%dHessian
        do i=1,n
            do j=1,n
                tResult%dHessian(i,j)=tResult%dHessian(i,j)+dSecond*tA%dGradient(i)*tA%dGradient(j)
            end do
        end do

    end function LogSO


    !> Apply the integer-power chain rule used by G, Q, B, and ternary polynomials.
    function PowerSO(tA,iExponent) result(tResult)

        type(SecondOrderScalar), intent(in) :: tA
        integer, intent(in) :: iExponent
        type(SecondOrderScalar) :: tResult
        integer :: i,j,n
        real(8) :: dFirst,dSecond

        n=SIZE(tA%dGradient)
        if (iExponent==0) then
            tResult=ConstantSO(1D0,n)
            return
        end if
        dFirst=DFLOAT(iExponent)*tA%dValue**(iExponent-1)
        if (iExponent==1) then
            dSecond=0D0
        else
            dSecond=DFLOAT(iExponent*(iExponent-1))*tA%dValue**(iExponent-2)
        end if
        tResult=ConstantSO(tA%dValue**iExponent,n)
        tResult%dGradient=dFirst*tA%dGradient
        tResult%dHessian=dFirst*tA%dHessian
        do i=1,n
            do j=1,n
                tResult%dHessian(i,j)=tResult%dHessian(i,j)+dSecond*tA%dGradient(i)*tA%dGradient(j)
            end do
        end do

    end function PowerSO

end module ModuleMQMQAUnconstrained
