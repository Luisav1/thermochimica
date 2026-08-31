# MQMQA Local Response and GEM Mapping Audit

This document records two related development tracks. The first sections trace
the completed plain-`SUBG` path from MQ-3A through the reusable MQ-4B correction
builder. The final section records the separate SUBQ MQ-3A audit that decides
which response and mapping concepts may be reused and which require new native
evidence.

## Scope

MQ-3A identifies the local variables, equality constraints, derived MQMQA
quantities, element-potential forcing, and candidate existing GEM baseline that
a future plain-`SUBG` response mapper must verify. It is a source-and-equation audit only.
It does not modify `GEMNewton`, inject a Hessian, add solver controls, or claim
native database coverage for production-family `Q` or `B` parameters.

The traced production path is:

1. `CompMolFraction` sends every non-`IDMX` phase to `Subminimization`.
2. `SubMinInit` treats the phase's solution species as local composition
   variables and constructs their element-potential chemical potentials.
3. `CompExcessGibbsEnergySUBG` evaluates the reference, configurational, and
   excess partial molars from the current quadruplet fractions.
4. `SubMinNewton` uses a simple diagonal local response together with
   normalization and, when present, charge neutrality.
5. `GEMNewton` uses the same simple species response to assemble the global
   element and phase-amount blocks.

## Independent Local Variables

For a plain-`SUBG` phase with `nQuad` quadruplet species, MQ-3B must construct
one authoritative local state from the active assemblage slot:

```text
x_q = dMolFraction(iFirst+q-1)
      / sum_r dMolFraction(iFirst+r-1),
N   = dMolesPhase(activeSlot),
n_q = N x_q.
```

The explicit normalization protects the response calculation from small drift
in the stored production fractions. `N` comes from the active phase amount, not
from `SUM(dMolesSpecies)`. `GEMNewton` floors each `dMolesSpecies` entry while
assembling its system, so the stored species amounts may not sum exactly to
`dMolesPhase(activeSlot)` and may differ from `N*x`. MQ-3B must report that
difference and must not silently treat the floored species amounts as the
authoritative physical state.

`ModuleMQMQAUnconstrained` differentiates with respect to `n_q`. At fixed phase
amount, an infinitesimal composition change satisfies

```text
dn = N dx.
```

Therefore, if `Hn = d(mu)/d(n)` is the verified extensive-mole Hessian, the
fixed-amount composition-space curvature required by the local response is

```text
Hx = N Hn.
```

Reference energy is linear in `n_q`, so its Hessian block is zero. `Hx`
contains the traced configurational and excess curvature.

## Units and Held-Fixed Quantities

Thermochimica normalizes Gibbs energy by `R*T`. The normalized extensive energy
`Ghat=G/(R*T)` therefore has amount units, while its mole derivative
`mu=d(Ghat)/d(n)` is dimensionless. Under this convention:

- `Hn = d(mu)/d(n)` has inverse-amount units;
- `Hx = N*Hn` is the dimensionless chemical-potential response with respect to
  phase composition;
- the element potentials `Gamma`, production partial molars `mu`, the
  stoichiometric forcing matrix `S`, and the composition `x` are dimensionless;
- the constrained response `R=dx/dGamma` is dimensionless;
- GEM contributions such as `N*S^T*R`, `Abase`, and the corresponding residual
  contribution have amount units, normally moles.

The local derivative and finite-difference response hold temperature `T`,
pressure `P`, phase amount `N`, phase topology, decoded parameter records, and
the selected thermodynamic-model branch fixed. Only the independent
element-potential forcing and the resulting interior phase composition vary.

## Independent Equality Constraints

### Normalization

Production `Subminimization` explicitly imposes

```text
sum_q x_q = 1,
```

so every admissible fixed-amount response obeys

```text
sum_q dx_q = 0.
```

This is the only local equality constraint for the uncharged plain-`SUBG`
prototype represented by the converged `CuFeC-Kang.dat` Liquid.

### Optional charge neutrality

When `iPhaseElectronID(iPhase) /= 0`, `SubMinNewton` adds a second bordered
row using the electron-component stoichiometry. Production does not divide
this row by `iParticlesPerMole`. Its exact entries and infinitesimal form are

```text
C_charge(q) = dStoichSpecies(iFirst+q-1, iPhaseElectronID(iPhase)),
sum_q C_charge(q) dx_q = 0.
```

This row is distinct from the per-particle element-forcing matrix `S`; a future
charged implementation must not reuse an `S` column as the charge constraint.

The MQ-3B response solver should accept a general constraint matrix so this
row can be represented. The first native response test remains deliberately
restricted to the uncharged plain-`SUBG` Liquid. A later production mapper
must either implement the charge row or explicitly reject charged `SUBG`
phases; it must not silently apply an uncharged response.

## Derived Identities, Not Additional Constraints

For fixed topology and model parameters, the following are deterministic
functions of the quadruplet amounts:

- quadruplet fractions;
- constituent site amounts and site fractions;
- equivalent constituent fractions;
- ordinary and zeta-weighted pair populations;
- pair fractions and weighted-pair fractions;
- the weighted-pair marginals `F1` and `F2`;
- binary `chi` and `xi` coordinates;
- supported ternary composition factors;
- the `S1`, `S2`, and `S3` configurational-energy quantities.

They are reconstructed in both the production evaluator and
`ModuleMQMQAUnconstrained`. They do not receive independent Lagrange
multipliers. Their first- and second-order responses already enter through the
chain rule in the verified gradient and Hessian.

The quadruplet identities, coordination numbers, zeta values, and asymmetric
group masks are fixed model data, not composition constraints. Positivity is
an inequality/domain requirement handled by admissible perturbations and later
globalization; it is not another row in the local equality-constraint matrix.
Global element mass balance belongs to GEM and supplies forcing/coupling. It is
not imposed again as a local phase constraint. The phase assemblage is also not
constrained by this work.

## Element-Potential Forcing

Define the per-particle quadruplet stoichiometry matrix

```text
S(q,e) = dStoichSpecies(iFirst+q-1,e)
         / iParticlesPerMole(iFirst+q-1).
```

`SubMinInit` constructs the chemical potential supplied by the element
potentials as

```text
mu_star(q) = sum_e S(q,e) Gamma(e).
```

Consequently, an element-potential perturbation supplies the local forcing

```text
d(mu_star) = S dGamma.
```

For a constraint matrix `C`, MQ-3B will use the symmetric bordered convention

```text
[ Hx    C^T ] [ dx      ] = [ S dGamma ]
[ C      0  ] [ dlambda ]   [     0    ].
```

The multiplier sign is conventional; choosing matching `C^T` and `C` blocks
keeps the KKT matrix symmetric. The resulting `dx` must be compared with the
existing mixed-sign RKMP response solver to demonstrate that changing the
multiplier convention does not change the physical response.

For the uncharged prototype, `C` is one row of ones. This is structurally the
same normalized local-response problem used by the final RKMP mapper; the
difference lies in the MQMQA curvature, not in an invented collection of site
or pair constraints.

### Tangent-space stability

A successful linear solve is necessary but is not evidence that the local
thermodynamic response is stable. MQ-3B must construct a basis `Z` for the
admissible tangent space,

```text
C Z = 0,
Ktangent = Z^T Hx Z,
```

and report:

- the numerical rank of `C`;
- the smallest and largest singular values of the bordered KKT matrix;
- a condition estimate;
- the eigenvalues or inertia of `Ktangent`;
- whether the converged state is positive definite on the tangent space;
- the constraint residual `||C dx||`;
- the complete KKT residual.

For the first native prototype, which is deliberately a stable, interior,
uncharged plain-`SUBG` state, positive definiteness of `Ktangent` is an
acceptance gate, not merely a reported diagnostic. MQ-3B must document the
eigenvalue tolerance used for this decision. If `Ktangent` is not positive
definite within that tolerance, the prototype must stop and classify the cause
as one or more of:

- an unstable or metastable thermodynamic state;
- a decoding or mole/composition scaling defect;
- a rank or constraint defect;
- numerical conditioning.

Treatment of genuinely indefinite states belongs to later rejection or
globalization work; the first response-verification case must not proceed
silently with an unstable tangent-space matrix.

The raw extensive Hessian has the expected radial null direction. Normalization
removes that direction, but MQ-3B must detect any additional tangent-space
degeneracy or instability. The null-space response

```text
dx = Z (Z^T Hx Z)^(-1) Z^T f
```

must provide an algebraic cross-check of the bordered solve.

### Independent forcing directions

Element-potential columns are not automatically independent or nonzero after
normalization. MQ-3B must project the forcing into the tangent space, determine
the numerical rank of `Z^T S`, and construct a forcing-basis matrix `P` whose
columns make

```text
Z^T S P
```

linearly independent. The nonlinear finite-difference oracle must perturb one
supported forcing direction at a time:

```text
Gamma -> Gamma +/- h P(:,j).
```

These comparisons verify the response on the supported forcing subspace. They
do not create independent information in null forcing directions. The full
element-coordinate response may be represented from the supported basis, but
true null-response columns must be checked separately with absolute residuals
rather than meaningless relative errors against zero.

## Existing Thermochimica Baseline

`SubMinNewton` does not use the model-specific `SUBG` Hessian. Its local matrix
uses

```text
Hbase(q,q) = 1/x_q
```

with the normalization border and optional charge row. The chemical-potential
residual contains the complete production `SUBG` partial molars, but the local
response matrix is the same ideal/simple diagonal approximation used for other
nonideal solution models.

`GEMNewton` assembles the matching global baseline. For one active solution
phase, its element block and element-to-phase-amount column are

```text
Aee_base = N S^T diag(x) S,
Aep_base = N S^T x.
```

The separate phase-amount unknown carries the radial phase-scaling direction.
Eliminating that direction gives the fixed-amount normalized response

```text
N S^T [diag(x) - x x^T] S,
```

which is algebraically the response obtained by condensing `diag(1/x)` with
`sum(dx)=0` for an uncharged normalized interior state. This source audit
therefore identifies the candidate plain-`SUBG` subtraction baseline. MQ-3B
must establish its numerical equivalence to the current Thermochimica
implementation before relying on that subtraction.

The MQMQA correction must therefore follow the established response-delta
architecture:

```text
corrected response = constrained solve with Hx = N Hn,
baseline response  = constrained solve with diag(1/x),
applied response   = corrected response - baseline response.
```

It must not add the total MQMQA Hessian directly to the GEM element block.

### Baseline reconstruction gate

MQ-3B must require three-way agreement among:

1. the symmetric bordered solve with `Hbase=diag(1/x)`;
2. the closed-form normalized response `diag(x)-x*x^T`;
3. the corresponding phase contribution reconstructed directly from the
   current `GEMNewton` arrays.

The third MQ-3B comparison is a **source-faithful phase-local reconstruction
from the arrays consumed by `GEMNewton`**. It reproduces the selected phase's
element block, element-to-phase-amount column, and residual formula without
mixing in contributions from other phases. It does not capture the internally
assembled `A` matrix directly. Live reduced-system verification in MQ-4 must
strengthen this evidence with a before-and-after snapshot around the selected
phase's production assembly or an equivalent internal diagnostic capture.

The comparison must also cover the raw element-block symmetry residual before
any explicit symmetrization, the element-to-phase-amount column, scaling with
phase amount and composition, and consistent amount scaling and units for the
element block, element-to-phase-amount column, and residual contribution. This
converts a source-derived identity into tested Thermochimica evidence.

## Candidate Reduced Contributions

Let `Rcorr` and `Rbase` be the constrained responses to all element-potential
columns, and let `rMuCorr` and `rMuBase` be the corresponding constrained
responses to the current production partial-molar vector. The RKMP derivation
suggests the candidate fixed-amount deltas

```text
deltaA = N S^T (Rcorr - Rbase),
deltaB = N S^T (rMuCorr - rMuBase).
```

`deltaA` is the leading MQ-4 candidate. `deltaB` remains explicitly provisional
until the GEM residual convention is derived from production source. Before
MQ-4, the derivation must begin with the current `GEMNewton` element residual,
including its `mu-1` terms, and establish:

- the residual forcing vector and its sign;
- why constant terms vanish under the normalization response;
- every phase-amount factor;
- why no additional solution-phase row correction is required, if that is the
  final result;
- how an off-equilibrium local stationarity residual enters the condensation.

Both candidates must first reproduce the baseline and then pass an independent
off-equilibrium reduced-residual or response finite difference. Formula
agreement inside a mapper implementation is not sufficient evidence by itself.

## MQ-3B Verification Requirements

The next stage should proceed in this order:

1. Implement a model-independent constrained-response solver accepting `H`,
   `C`, and multiple forcing columns. Unit-test synthetic positive-definite,
   singular, redundant-constraint, and ill-conditioned cases.
2. Compare the symmetric KKT convention against the existing RKMP response
   solver and the null-space solve. Implement this solver independently; do not
   refactor or replace the validated RKMP production path during MQ-3B.
3. Numerically pass the three-way diagonal-baseline gate using a source-faithful
   phase-local reconstruction from the arrays consumed by `GEMNewton`. Document
   that this is not yet a capture of the internally assembled matrix; that live
   reduced-system check belongs to MQ-4.
4. Decode the converged uncharged `CuFeC-Kang.dat` Liquid using normalized
   production `dMolFraction`, `N=dMolesPhase(activeSlot)`, and `n=N*x`. Report
   the discrepancy from the floored `dMolesSpecies` values.
5. Obtain `Hn`, the production gradient, `S`, `C`, `x`, and `N`; audit tangent
   stability and the independent rank of `Z^T*S`. Require tangent-space
   positive definiteness for this first stable prototype or stop and classify
   the failure.
6. Solve the corrected and baseline responses for an independent forcing
   basis `P`, verify `Z^T*S*P` has full column rank, represent the response on
   that supported forcing subspace, check null directions separately, and
   verify both parts of the complete KKT residual,
   `Hx*R+C^T*Lambda-F` and `C*R`. Retain the multiplier-eliminated projected
   residual only as an additional algebraic cross-check. Use an orthonormal
   tangent basis when reporting eigenvalues and tangent condition estimates.
7. Independently perturb each forcing-basis direction and reconverge the local
   composition using the established production `CompExcessGibbsEnergySUBG`
   partial molars. Build the test-only nonlinear oracle from finite differences
   of those production partial molars; do not use the new analytic MQMQA Hessian
   in the oracle Jacobian. Rebuild that finite-difference Jacobian at each final
   perturbed root and use it, rather than the analytic Hessian, to convert the
   remaining stationarity residual into an oracle-uncertainty estimate.
8. Keep phase amount, phase type, topology, active parameter records, and
   production branch fixed. Warm-start from the unperturbed composition, save
   and restore every touched global array, and reject branch changes, failed
   positivity, boundary clipping, or failed convergence.
9. Classify a finite-difference point as oracle-resolved only when the estimated
   stationarity uncertainty satisfies `u(h) <= 0.5*e(h)`. Require the accepted
   second-order region to contain resolved points, apply the accuracy gate to
   the best resolved point, and report both the raw minimum and best resolved
   error. Use a positivity-safe sweep with observed-order reporting,
   componentwise errors, worst-quadruplet indices, and vector norms rather than
   accepting one favourable perturbation.
10. Keep the complete MQ-3B path diagnostic-only and leave `GEMNewton`
    unchanged. Stop for review after the native local response passes. MQ-3B
    does not authorize assembly of `deltaA`, assembly of `deltaB`, activation
    in `GEMNewton`, or refactoring of the validated RKMP production mapper.

The public production `Subminimization` path remains relevant source evidence,
but its historical update-size stopping rule is intentionally looser than a
derivative finite-difference oracle requires. MQ-3B therefore solves the same
production partial-molar stationarity equations to a documented tighter
tolerance using a test-only finite-difference Jacobian. This oracle remains
independent of the new analytic MQMQA Hessian and does not modify production
`Subminimization` behavior.

Only after this response agrees should MQ-4 construct and independently verify
the reduced `A` and `B` deltas.

## Database-Coverage Boundary

MQ-3A does not require an assessed plain-`SUBG` database containing `Q` or
`B`. The variable space, normalization, forcing, baseline response, and
condensation architecture are common to the supported plain-`SUBG` families.
The available assessed `G` case is sufficient for the first native response
and mapping tests.

An authoritative `Q`/`B` database becomes fundamental only before claiming or
releasing **database-backed production support for those parameter families**.
Without it, the project may still complete:

- standalone mathematical verification of the traced `Q` and `B` scalar forms;
- native `G` decoding and local-response verification;
- a plain-`SUBG` mapper exercised by assessed `G` cases.

It may not claim native thermodynamic verification or regression coverage for
plain-`SUBG` `Q` or `B`. If no assessed database is available, that limitation
must remain explicit; a synthetic parser fixture could test software mechanics
but would not replace thermodynamic evidence.

## MQ-3A Exit Decision

The state, unit, constraint, charge-row, and symmetric KKT sign conventions are
now fixed sufficiently for MQ-3B to begin. The uncharged plain-`SUBG` local
response uses explicitly normalized production quadruplet fractions and the
active `dMolesPhase` amount, with one normalization constraint. MQMQA site,
pair, coordination, zeta, and ternary quantities are derived functions, not
additional constrained unknowns. The verified total quadruplet-mole Hessian
converts to fixed-amount composition curvature as `Hx=N*Hn`.

The source audit identifies the bordered `diag(1/x)` response as the candidate
Thermochimica baseline and establishes its algebraic equivalence only for the
uncharged normalized interior case. MQ-3B must still demonstrate numerical
baseline equivalence, tangent-space stability, independent forcing coverage,
source-faithful phase-local GEM baseline reconstruction, and agreement with the
production nonlinear oracle on the supported forcing basis. Direct capture of
the live assembled reduced system, response-delta mapping, and GEM integration
remain MQ-4 work; no production solver behavior has changed.

## MQ-4A Diagnostic Reduced Mapping

MQ-4A keeps the correction outside `GEMNewton` and first derives it from the
actual species update used by `GEMLineSearch`. For species `q`, that update is
proportional to

```text
1 + lambdaPhase + S(q,:)*Gamma - mu(q).
```

Here `lambdaPhase` is the solution-phase radial update variable. It changes the
total amount of the phase without being an additional composition coordinate.
The local composition response remains tangent to normalization.

Let `C` be the row of ones, `S` the quadruplet-to-element stoichiometric
forcing matrix, and `Hx=N*Hn`. The corrected element-potential response and
the response to the current species residual are defined by the candidate KKT
systems

```text
[ Hx  C^T ] [ Rcorr   ] = [ S ]
[ C    0  ] [ LambdaS ]   [ 0 ]

[ Hx  C^T ] [ rMuCorr ] = [ fMu ]
[ C    0  ] [ lambdaMu ]   [  0  ].
```

The corresponding baseline objects use `Hbase=diag(1/x)`. The production
element residual contains `mu-1`, so the source-faithful choice is
`fMu=mu-1`. Solving with `mu` gives the same composition response: the two
right-hand sides differ by `C^T`, and that constant forcing is absorbed by the
normalization multiplier because `C*rMu=0`. MQ-4A verifies this cancellation
numerically rather than assuming it.

Expanding the species update in the element-balance equation gives the
phase-local candidate corrections

```text
deltaA = N*S^T*(Rcorr-Rbase),
deltaB = N*S^T*(rMuCorr-rMuBase).
```

The two terms must also satisfy one combined affine-response identity. For a
trial element-potential update `dGamma`, the corrected and baseline local
composition changes are

```text
dxCorr = Rcorr*dGamma-rMuCorr,
dxBase = Rbase*dGamma-rMuBase.
```

Their mapped difference is therefore

```text
N*S^T*(dxCorr-dxBase) = deltaA*dGamma-deltaB.
```

The minus sign follows from the species forcing
`S*dGamma-(mu-1)`. MQ-4A verifies this complete identity with a single
nonlinear finite-difference experiment, so separate sign errors in `deltaA`
or `deltaB` cannot be hidden by using two unrelated tests.

Both corrections have amount units and target only the element equations.
The current element-to-solution-phase column remains `N*S^T*x`; its transpose
and the current solution-phase residual `SUM(n*mu)` are unchanged because the
new curvature changes the tangent response, not the current state or radial
direction. A later builder must therefore return zero corrections for those
objects unless a broader derivation proves otherwise.

Two production chemical-potential conventions must not be mixed:

- after `CompChemicalPotential`, `dChemicalPotential` already contains the
  complete SUBG partial molar;
- after a direct low-level call to `CompExcessGibbsEnergySUBG`, the complete
  partial molar is `dChemicalPotential+dPartialExcessGibbs`.

Adding the excess array to the first form would double count it.

`GEMNewton` floors each live species amount before assembling its baseline.
The candidate response instead uses the authoritative thermodynamic state
`n=N*x`. MQ-4A reports their discrepancy and accepts the candidate only when
it is immaterial. The test temporarily installs the same positive,
off-equilibrium composition used by the nonlinear response oracle, recomputes
production chemical potentials, and requests the opt-in snapshot in
`ModuleGEMNewtonDiagnosticCapture`. The snapshot copies the completed baseline
immediately before any experimental correction or solve. Contributions from
the other active solution phases are removed from that captured system to
expose the selected Liquid contribution. Every modified production array is
then restored to its converged value.

`TestMQMQAGEMMappingVerification` verifies:

1. live phase-local element block, element residual, solution-phase column,
   and solution-phase residual reconstruction;
2. zero correction when corrected and baseline responses are identical;
3. finite values, normalization residuals, `deltaA` symmetry, constant-forcing
   cancellation, and non-vacuous `deltaB`;
4. the complete `deltaA`, `deltaB`, and combined affine identity against
   independent central finite differences of nonlinear states converged with
   production SUBG partial molars;
5. second-order convergence, oracle-resolution, normwise error, and maximum
   scaled component error;
6. every `deltaA` column using its own resolved step and uncertainty after both
   are mapped into the element-residual space through `N*S^T`.

An oracle point is resolved only when its mapped, scale-normalized uncertainty
is no more than half of the measured mapped error. This prevents a small
composition-space residual or an excellent aggregate matrix norm from hiding
an unresolved element direction.

MQ-4A itself did not provide a reusable correction builder or mutate a GEM
matrix. MQ-4B packages that verified result as described below. Controls,
solver activation, trust selection, and globalization remain MQ-4C concerns.

## MQ-4B Reusable Correction Builder

`ModuleMQMQAResponseMapping` turns the MQ-4A derivation into two deliberately
separate software operations:

```text
BuildMQMQAGEMCorrection
    reads one live eligible plain-SUBG phase
    returns its complete unscaled deltaA and deltaB

ApplyMQMQAGEMCorrection
    receives caller-owned A and B arrays
    adds a caller-selected fraction only to the element equations
```

Neither routine is called by normal `GEMNewton` execution in MQ-4B. The
existing diagnostic capture hook is unchanged, and no corrected Newton solve
occurs.

### Builder layers

The live builder checks that the requested phase is active, plain `SUBG`,
uncharged, positive, interior, finite, and dimensionally consistent. It then
uses `ModuleMQMQAProductionAdapter` and
`CompMQMQAHessianUnconstrained` to obtain the production-decoded mole Hessian.
It does not duplicate database parameters or thermodynamic formulas.

The response algebra is isolated in `BuildMQMQAReducedCorrection`. This
state-based kernel accepts

```text
N, x, mu, S, Hx
```

and constructs the baseline curvature, normalization constraint, corrected
and baseline element responses, and corrected and baseline residual responses.
Separating this kernel from production decoding has two purposes:

1. decoding failures and response-solve failures receive different statuses;
2. singular response systems can be tested directly without a test-only switch
   in the live builder.

Both entry points clear `deltaA` and `deltaB` before checking inputs. A failed
stage therefore cannot expose a partially assembled correction.

### Applicability and status

`lApplicable` is false for a phase outside the supported model scope, including
a charged phase. It becomes true after a live phase passes model and domain
eligibility. A later decoding, Hessian, response, or correction failure keeps
`lApplicable` true but returns a nonzero stage-specific status. This separates
ordinary inapplicability from a numerical or implementation failure in a phase
that should have been mappable.

The public status constants distinguish:

- success;
- not applicable;
- unsupported charged phase;
- invalid input;
- decoding or Hessian failure;
- corrected or baseline element-response failure;
- corrected or baseline residual-response failure;
- invalid correction or application.

### Structural checks

Every successful build requires:

- complete corrected and baseline KKT residuals below tolerance;
- normalization residuals below tolerance;
- finite `deltaA` and `deltaB`;
- raw `deltaA` symmetric within numerical tolerance, followed by removal of
  roundoff-level skew so the returned matrix is exactly symmetric;
- correct element-space dimensions.

The builder reads complete production `dChemicalPotential` values after
`CompChemicalPotential`. It does not add `dPartialExcessGibbs`, and it does not
write any Thermochimica global array.

The applicator validates dimensions, finite inputs, correction symmetry, and a
caller-owned weight in `[0,1]`. A materially nonsymmetric correction is rejected
before either caller-owned array is changed. The weight is a future
globalization choice applied to the fully formed correction. It does not scale
the local MQMQA Hessian. The applicator can change only

```text
A(1:nElements,1:nElements)
B(1:nElements).
```

Phase rows and columns, phase residuals, pure-phase equations, and production
thermodynamic state remain untouched.

### MQ-4B verification boundary

`TestMQMQAGEMMappingVerification` retains the independent MQ-4A construction
and nonlinear production-partial-molar oracle. It additionally verifies that:

1. builder `deltaA` and `deltaB` match the independent MQ-4A reference;
2. the builder changes no production state array;
3. inapplicable, invalid-input, charged, and singular-response cases return zero
   corrections with distinct statuses;
4. full and zero-weight application modify copied arrays exactly as expected;
5. sequential supplied correction pairs add linearly;
6. invalid alpha or materially nonsymmetric application leaves copied arrays
   unchanged;
7. all independent `deltaA`, `deltaB`, affine, and per-column oracle checks
   continue to pass.

MQ-4C may call this builder from an experimental solver path and decide when a
correction is trustworthy. That future stage must add controls, phase-by-phase
aggregation in the live iteration, alpha-zero reference solves, and nonlinear
globalization tests. MQ-4B does none of those things.

## SUBQ MQ-3A Local Response Reuse Audit

### Scope

SUBQ MQ-3A determines whether the completed plain-`SUBG` local-response
derivation applies to production `SUBQ`, now that the SUBQ scalar energy,
gradient, Hessian, production decoder, and native FeTiVO derivative comparison
have passed. This is a source-and-equation audit only. It does not broaden
`BuildMQMQAGEMCorrection`, modify `GEMNewton`, or claim constrained-response
verification for SUBQ.

The audit distinguishes three questions that must not be collapsed:

1. Does production `SUBQ` use the same local unknowns and constraints as
   production `SUBG`?
2. Does the verified SUBQ mole Hessian enter that response with the same
   mole-to-composition scaling?
3. Has the resulting constrained SUBQ response been independently verified?

The source answers the first two questions affirmatively. The third remains
SUBQ MQ-3B work.

### Shared production path

There is no separate production `CompExcessGibbsEnergySUBQ` routine.
`CompExcessGibbsEnergy` dispatches both model labels through
`CompExcessGibbsEnergySUBG`, which selects the SUBG or SUBQ thermodynamic
formula internally. Outside that model evaluator, the production response path
does not branch on `SUBG` versus `SUBQ`:

```text
CompMolFraction / Subminimization
    -> solution-species mole fractions
    -> per-particle element-potential forcing
    -> generic normalization and optional charge rows

GEMNewton
    -> generic solution-species element block and residual
    -> generic solution-phase amount row and column
```

`SubMinNewton` constructs the same diagonal local approximation
`diag(1/x)` for every nonideal solution model. `GEMNewton` similarly assembles
its baseline from `dMolFraction`, `dMolesPhase`, `dStoichSpecies`, and
`iParticlesPerMole` without a SUBG/SUBQ branch. SUBQ therefore does not require
a new baseline formula merely because its thermodynamic curvature differs.

### Local variables and amount scaling

For an active SUBQ phase, the independent local thermodynamic variables remain
the moles of its quadruplet solution species:

```text
x_q = normalized dMolFraction for quadruplet q,
N   = dMolesPhase(activeSlot),
n_q = N*x_q.
```

The standalone SUBQ routine and the native MQ-2B comparison use the same
unconstrained quadruplet-mole coordinates. The complete production partial
molar is the derivative with respect to `n_q`; no reference-species or gauge
transformation is required. Consequently,

```text
Hn = d(mu)/d(n),
dn = N*dx at fixed N,
Hx = N*Hn.
```

These units and chain-rule relations are model-independent. Pair-specific
zeta changes values inside `mu`, `Hn`, and `Hx`; it does not change the
independent mole coordinate or introduce a second phase-amount scale.

As in the SUBG audit, `N*x` is the authoritative interior thermodynamic state.
The species amounts floored by `GEMNewton` are an assembly safeguard and must
be compared with, rather than silently substituted for, `N*x`.

### Equality constraints and derived quantities

The mandatory fixed-amount constraint remains quadruplet normalization:

```text
C_normalization = [1, 1, ..., 1],
C_normalization*dx = sum_q dx_q = 0.
```

`SubMinNewton` has generic machinery to add a charge-neutrality row when
`iPhaseElectronID(iPhase) /= 0`. Current `CompThermoData`, however, assigns a
nonzero phase electron ID only to `SUBL` and `SUBLM`. The current production
implementation therefore treats SUBQ phases as normalization-only; this is an
implementation fact, not a claim that charged SUBQ is mathematically
impossible. The assessed FeTiVO response test must still assert
`iPhaseElectronID(SlagBsoln)==0` so the current behavior is recorded as
executable evidence rather than silently inherited.

SUBQ changes several derived thermodynamic quantities:

- the `S3` pair-log and equivalent-fraction exponents;
- the fixed environment incidence weights used by `chi`;
- the pair-specific zeta-weighted populations and their propagated marginals;
- the assessed production-family `Q` contributions in FeTiVO;
- potentially other supported decoded terms in future databases.

All remain deterministic functions of the quadruplet moles and fixed model
data. They enter the verified gradient and Hessian through the chain rule. They
are not independent response unknowns and do not receive separate Lagrange
multipliers. Site fractions, pair fractions, zeta-weighted pair distributions,
equivalent fractions, `chi`, `xi`, and `F` quantities therefore add no rows to
the KKT constraint matrix.

Positivity remains an inequality/domain requirement. Global element balance
remains the GEM equation that supplies forcing. The phase assemblage remains
unconstrained by this work.

### Element-potential forcing and baseline response

The SUBQ forcing matrix is the same production per-particle stoichiometric
matrix used for SUBG:

```text
S(q,e) = dStoichSpecies(iFirst+q-1,e)
         / iParticlesPerMole(iFirst+q-1).
```

Thus the candidate corrected response continues to satisfy

```text
[ Hx  C^T ] [ Rcorr   ] = [ S ]
[ C    0  ] [ LambdaS ]   [ 0 ],
```

while the production baseline uses `Hbase=diag(1/x)`. Define the closed-form
normalization-only species-response operator

```text
Mbase = diag(x) - x*x^T.
```

For the element-potential right-hand side used above, the corresponding
composition response is

```text
Rbase = Mbase*S.
```

After mapping through the per-particle stoichiometry, the candidate phase-local
response delta retains the completed SUBG form

```text
deltaA = N*S^T*(Rcorr-Rbase).
```

For the element residual, define the current off-stationary species forcing

```text
fMu = mu - 1.
```

The corrected residual response is obtained from

```text
[ Hx  C^T ] [ rMuCorr ] = [ fMu ]
[ C    0  ] [ lambdaMu ]   [  0  ],
```

and the analogous baseline solve replaces `Hx` by `Hbase=diag(1/x)` to obtain
`rMuBase`. The candidate residual correction is therefore

```text
deltaB = N*S^T*(rMuCorr-rMuBase).
```

This freezes the forcing, sign, phase-amount factor, and `mu-1` convention.
Adding a constant-one vector to the forcing changes only the normalization
multiplier, not the tangent composition response, but the subsequent SUBQ
mapping test must reverify that cancellation and the combined affine identity
before the live builder is broadened.

The ordinary positive-interior baseline is the only baseline path in the first
SUBQ response test. `SubMinNewton` also clamps or bounds diagonal entries in a
fallback solve after its arrow solver fails or compositions become extremely
small. MQ-3B deliberately excludes that fallback behavior: every stationary
and perturbed composition must remain strictly positive and the ordinary
`1/x` solve must succeed. Boundary and fallback support belongs to later
eligibility or globalization work.

### Production chemical-potential convention

The same convention applies to both production model labels:

- after the normal `CompChemicalPotential` path, `dChemicalPotential` already
  contains the complete SUBQ partial molar because `CompExcessGibbsEnergy`
  adds `dPartialExcessGibbs`;
- immediately after a direct low-level call to
  `CompExcessGibbsEnergySUBG`, the complete SUBQ partial molar is
  `dChemicalPotential+dPartialExcessGibbs`.

MQ-2B verified the second convention against the standalone analytic gradient.
A future live builder must use the first convention and must not add the excess
array again.

### Evidence already available

The assessed FeTiVO `SlagBsoln` case provides the following SUBQ evidence
before constrained response is attempted:

- 15 quadruplet species and a 14-dimensional normalization tangent space;
- strict production SUBQ decoding with six `G` and eight `Q` records;
- machine-level scalar and direct-gradient parity;
- raw Hessian symmetry and radial homogeneity;
- second-order finite differences of production partial molars in all 14
  independent total-preserving directions;
- a positive interior seed obtained by blending the converged composition
  without changing total phase amount.

This evidence establishes the local derivative object and production coordinate
contract. It does not establish tangent-space positive definiteness, response
to element-potential forcing, the rank of the supported forcing subspace, or
nonlinear response agreement. Those are MQ-3B gates.

FeTiVO has uniform zeta `2.4` and no `B` or `R` records. The first response test
is therefore an assessed SUBQ `G/Q` test, not native response verification of
every SUBQ feature. Production and the analytic implementation now both use the
published zeta-weighted SUBQ `S3` definition. A separate controlled modified-
runtime regression distinguishes that selection from the legacy unweighted form
and verifies the corresponding weighted-normalization derivative in production
partial molars;
FeTiVO itself remains unable to distinguish them because its zeta is uniform.

### Reuse decision

| Response component | SUBQ MQ-3A decision |
| --- | --- |
| Quadruplet-mole variables and `n=N*x` state | Reuse |
| Mole-to-composition scaling `Hx=N*Hn` | Reuse |
| Normalization constraint and generic charge machinery | Reuse; current SUBQ is normalization-only and must assert that state |
| General bordered KKT solver | Reuse unchanged |
| Null-space and forcing-basis construction | Reuse unchanged |
| Per-particle stoichiometric forcing `S` | Reuse construction, rebuild from FeTiVO arrays |
| Diagonal `diag(1/x)` GEM baseline and `Mbase=diag(x)-x*x^T` operator | Reuse formulas, reverify numerically for SUBQ |
| Nonlinear finite-difference oracle architecture | Reuse structure; converge production SUBQ stationarity at every forcing |
| Production decoder | Use strict `DecodeProductionSUBQPhase` |
| Analytic curvature | Use verified SUBQ `CompMQMQAHessianUnconstrained` result |
| Live MQ-4B builder | Not reused during MQ-3A; a strict SUBQ entry point was added only after MQ-4A passed |
| `GEMNewton` activation | Out of scope until SUBQ response and mapping pass |

### SUBQ MQ-3B requirements

The next stage should add a dedicated native SUBQ response test rather than
silently routing SUBQ through the completed SUBG test:

1. Reproduce the assessed FeTiVO state and construct the same strictly positive
   20%-blended `SlagBsoln` composition used by MQ-2B. Treat this composition as
   a nonlinear-solve seed, not automatically as the response-verification
   state.
2. Assert the model label, active phase, zero electron ID, 15-quadruplet
   topology, six `G` records, eight `Q` records, uniform zeta, and preserved
   phase amount.
3. Decode with `DecodeProductionSUBQPhase`, reconstruct `S` from production
   stoichiometry and particle counts, and retain the converged production
   element potentials `Gamma`.
4. At fixed `T`, `P`, `N`, `Gamma`, model, topology, parameters, and phase
   identity, solve the production SUBQ partial-molar stationarity equations
   from the blended seed to a nearby stationary composition `x0`. Report the
   seed-to-root change, final stationarity residual, and `MINVAL(x0)`.
5. Recompute the complete production partial molars and analytic SUBQ Hessian
   at `n0=N*x0`, then form `Hx=N*Hn`. Build an orthonormal
   normalization-tangent basis `Z` and report the inertia/eigenvalues of
   `Z^T*Hx*Z`. Positive definiteness is required to claim that `x0` is a stable
   local minimum. Otherwise classify the stationary point and do not complete
   the first stable-response prototype with that state.
6. Determine the independent rank of `Z^T*S` with a rank-revealing,
   column-pivoted QR construction and construct the supported element-potential
   forcing basis `P`. Report the selected element names and columns, rank, and
   null forcing directions.
7. Pass the three-way ordinary-interior baseline comparison among the bordered
   `diag(1/x0)` solve, `[diag(x0)-x0*x0^T]*S`, and a source-faithful phase-local
   reconstruction from the FeTiVO arrays consumed by `GEMNewton`. Do not invoke
   or claim coverage of the clamped fallback baseline.
8. Solve corrected and baseline responses on the independent forcing basis.
   Gate the KKT top-block residual, constraint residual, projected residual,
   KKT-versus-null-space agreement, and normwise and componentwise response
   residuals. Include a synthetic null forcing in `range(C^T)` and require zero
   composition response.
9. For every supported direction `P(:,j)` and perturbation size `h`, hold the
   quantities listed in step 4 fixed except the forcing, set
   `GammaPlus=Gamma+h*P(:,j)` and `GammaMinus=Gamma-h*P(:,j)`, and re-solve the
   same nonlinear production SUBQ stationarity equations from `x0` to obtain
   `xPlus` and `xMinus`. Do not globally re-equilibrate or permit a phase,
   topology, parameter, or model-branch change.
10. Use `(xPlus-xMinus)/(2*h)` as the independent nonlinear response oracle.
    Require every nonlinear root and every intermediate accepted state to
    remain strictly positive; clipping is not permitted.
11. Rebuild the production-partial-molar finite-difference tangent Jacobian at
    each converged `xPlus` and `xMinus` root. Estimate the oracle uncertainty
    from the nonlinear stationarity residuals using those independent
    Jacobians, not the analytic SUBQ Hessian under test. Classify a point as
    resolved only when `u(h) <= 0.5*e(h)`. Apply convergence-order gates only
    to resolved points. Require and report the best resolved normwise error,
    maximum scaled component error and worst quadruplet, observed order, and
    the raw minimum separately.
12. Keep MQ-3B diagnostic-only. Do not broaden the MQ-4B builder, construct live
    GEM corrections, or alter `GEMNewton` in the same stage.

### SUBQ MQ-3A exit decision

SUBQ uses the same production response coordinates, amount scaling, equality
constraints, element-potential forcing, and historical GEM baseline as plain
SUBG. Its pair-specific zeta, revised configurational powers, chi incidence,
and `Q` terms change the thermodynamic gradient and curvature but do not create
new independent constraints. The verified unconstrained SUBQ Hessian therefore
has the correct variable space for the existing constrained-response
architecture.

SUBQ MQ-3A is complete as a derivation and source audit. SUBQ MQ-3B may now use
the assessed FeTiVO `G/Q` case to converge and classify a positive stationary
root, then verify its native constrained response. The existing live correction
builder must remain SUBG-only until that response and the subsequent
reduced-mapping evidence pass.

### SUBQ MQ-3B exit result

`TestMQMQASUBQResponseVerification.F90` now performs the dedicated native
response experiment required above. The 20%-blended FeTiVO composition is used
only as a positive nonlinear-solve seed. At fixed temperature, pressure, phase
amount, production element potentials, model branch, topology, and decoded
parameters, the independent production-partial-molar solver moves by
`1.133449E-01` in composition two-norm to a stationary root. The root remains
strictly positive (`MINVAL(x0)=1.128405E-09`) and its tangent stationarity
residual is `1.348160E-14`.

The decoded root contains the assessed 15-quadruplet, six-`G`, eight-`Q`,
uniform-zeta (`2.4`) SUBQ model. Its orthonormal tangent Hessian is positive
definite: the minimum and maximum tangent eigenvalues are `1.352461E+00` and
`8.271186E+08`. The resulting condition estimate, `6.115654E+08`, is reported
because forward-solution comparisons must be interpreted separately from
backward equation residuals.

All four independent element-potential forcing directions are supported. A
rank-revealing column-pivoted QR construction selects them in the order O, Fe,
Ti, and V. The bordered KKT solve has a top-equation residual of `1.776357E-15`, a
normalization residual of `9.150666E-17`, and a projected residual of
`1.554312E-15`. A separate null-space solve agrees with the bordered response
within `3.686123E-09`, below its condition-aware `1.357948E-07` forward-error
bound; its scaled reduced-equation residual is `7.578160E-17`. A synthetic
normalization-only forcing produces exactly zero composition response.

Three ordinary-interior GEM response/block formulations plus the residual
reconstruction agree independently:

- bordered ideal solve versus the closed normalized response: `1.292369E-16`;
- phase-local centered element block reconstruction: `3.383117E-16`;
- element-to-phase reconstruction: `1.110223E-16`;
- phase-local residual reconstruction: `2.842327E-16`.

For every supported forcing direction, independently reconverged production
SUBQ roots at `Gamma +/- h*P(:,j)` exhibit a resolved second-order region. The
uncertainty classifier rebuilds the production-partial-molar finite-difference
tangent Jacobian at each converged plus/minus root; it does not use the analytic
SUBQ Hessian being tested. The worst raw minimum response error is
`2.049665E-11`; after enforcing the oracle
resolution rule `u(h)<=0.5*e(h)`, the worst best resolved normwise error is
`4.981605E-10`, the worst corresponding scaled component error is
`3.755623E-10`, and the worst accepted stationarity uncertainty is
`1.232034E-11`.

The stationary root is interior under the test contract but close to the
composition boundary, and its tangent condition estimate is approximately
`6.12E+08`. MQ-3B therefore establishes the correctness of the analytic local
response in this difficult state. It does not establish that a future live
mapper should activate at every similarly small minimum fraction or condition
number; MQ-4 eligibility and globalization must define and test that policy.

This completes native constrained-response verification for the assessed
uniform-zeta FeTiVO SUBQ `G/Q` case. It does not verify SUBQ `B`, `R`, an
assessed database-native nonuniform-zeta response, reduced `deltaA/deltaB`
mapping, the live correction builder, or `GEMNewton` activation. The separate
controlled runtime regression covers the corrected nonuniform-zeta local
thermodynamics, not constrained response. The next SUBQ stage is a
diagnostic reduced-mapping verification analogous to MQ-4A; the existing
production correction builder remains strict plain-`SUBG` until that gate
passes.

### SUBQ MQ-4A exit result

`TestMQMQASUBQGEMMappingVerification.F90` now verifies the complete
phase-local reduced mapping for the assessed FeTiVO `SlagBsoln` SUBQ `G/Q`
case. The test uses a strictly positive, deliberately off-equilibrium
composition while holding temperature, pressure, phase amount, topology,
decoded parameters, and model branch fixed. It does not call the existing
plain-`SUBG` correction builder and does not modify `GEMNewton` behavior.

At that same off-equilibrium state, the diagnostic capture around
`GEMNewton` agrees with a source-faithful reconstruction of the selected
phase's historical baseline contribution. The scaled differences are
`3.477985E-16` for the element block, `1.464734E-15` for the element residual,
and exactly zero for both the element-to-phase column and solution-phase
residual. The authoritative `N*x` species amounts and GEM's floored species
amounts are identical in this state.

The corrected and historical constrained responses use the same normalization
constraint and forcing conventions. The candidate corrections are

```text
deltaA = N*S^T*(R_corrected - R_baseline)
deltaB = N*S^T*(r_mu,corrected - r_mu,baseline),
```

where the off-equilibrium residual forcing is `mu-1`. A constant-one forcing
produces no normalized composition response. The constrained corrected and
baseline responses to `mu` and `mu-1` agree to `8.869314E-16` and
`1.751567E-15`, respectively.

The complete KKT equations are checked independently for corrected and
baseline element forcing and corrected and baseline `mu-1` forcing. The four
scaled top-equation residuals are `1.468099E-15`, `3.058540E-17`,
`8.581426E-16`, and exactly zero. Their normalization-constraint residuals are
`9.714451E-17`, `1.266348E-16`, `1.665335E-15`, and `1.104672E-14`. The
candidate `deltaA` symmetry residual is `2.143770E-16`, and `deltaB` is
non-vacuous with maximum magnitude `7.104417E+00`.

Independent nonlinear production-partial-molar oracles verify the full
reduced mapping. Every perturbed state is reconverged without using the
analytic response as the nonlinear solver Jacobian. Oracle uncertainty is
estimated from a production-partial-molar finite-difference tangent rebuilt at
each converged plus/minus state. The results are:

- complete columnwise-reconstructed `deltaA` matrix error: `2.597968E-09`;
- complete `deltaA` matrix maximum scaled component error: `1.509257E-09`;
- worst oracle-resolved `deltaA` column error: `2.062968E-09`, with uncertainty
  `1.100866E-10` for that column;
- best resolved mixed-forcing `deltaA*dGamma` error: `3.195113E-10`, with
  corresponding maximum scaled component error `2.110043E-10`;
- complete off-equilibrium `deltaB` error: `4.188294E-09` at its best resolved
  step;
- combined affine identity error for `deltaA*dGamma-deltaB`: `4.520129E-09`
  at its best resolved step;
- best-resolved maximum scaled component errors of `4.786684E-09` and
  `5.142649E-09` for `deltaB` and the affine identity, respectively.

All three sweeps contain the expected resolved second-order region. Smaller
steps that fall below the independent oracle's resolving power are reported
but excluded from acceptance. Each reconstructed `deltaA` column is gated
against its own stored uncertainty, preventing an unresolved column from being
hidden by the full-matrix norm.

SUBQ MQ-4A therefore establishes the signs, amount scaling, `mu-1` convention,
affine reduced-GEM interpretation, and live phase-local baseline for this
assessed SUBQ `G/Q` case. It does not establish a reusable SUBQ correction
builder, correction application, multiple active SUBQ phase aggregation,
trust/globalization, or live solver activation. The next stage is SUBQ MQ-4B:
extend the reusable correction-builder contract only after preserving the
strict plain-`SUBG` path and all current rejection behavior.

### SUBQ MQ-4B exit result

`ModuleMQMQAResponseMapping.f90` now exposes two strict production entry
points. `BuildMQMQAGEMCorrection` continues to accept only plain `SUBG`, while
`BuildMQMQASUBQGEMCorrection` accepts only `SUBQ`. Both use the same verified
reduced-response kernel after dispatching through their own production decoder;
neither entry point can silently reinterpret the other MQMQA model type.

The SUBQ builder packages the MQ-4A equations as a reusable, state-preserving
operation. It reads the active FeTiVO `SlagBsoln` state, decodes its assessed
`G/Q` model, evaluates the complete local Hessian, constructs corrected and
historical constrained responses, and returns the unscaled phase-local
`deltaA` and `deltaB`. It does not mutate global thermodynamic state or any live
GEM array. The existing applicator remains model-independent and changes only
caller-owned copies of the element block and element residual.

The extended `TestMQMQASUBQGEMMappingVerification.F90` compares this reusable
builder against the independent MQ-4A construction at the same positive,
off-equilibrium state. The measured results are:

- builder status `MQMQA_MAP_SUCCESS`, with the phase marked applicable;
- builder versus independent `deltaA` error: `1.095489E-16`;
- builder versus independent `deltaB` error: exactly zero;
- production-state mutation metric: exactly zero;
- copied-array full-application error at alpha one: exactly zero;
- copied-array alpha-zero no-op error: exactly zero;
- sequential correction-pair additivity error: exactly zero.

The test also verifies safe zero outputs and distinct status behavior for
invalid input, charged phases, singular corrected responses, nonsymmetric
application data, and alpha outside `[0,1]`. Both model-specific entry points
are tested against the opposite phase type and return `NOT_APPLICABLE` without
producing a correction. The builder returns an exactly symmetric `deltaA` only
after the raw matrix passes its symmetry gate.

This completes the diagnostic-only SUBQ MQ-4B software contract for the
assessed uniform-zeta FeTiVO `G/Q` case. It does not activate the correction in
`GEMNewton`, select alpha, define trust/globalization policy, prove physical
aggregation of multiple active SUBQ phases, or extend native evidence to `B`,
`R`, magnetism, or nonuniform-zeta assessed response data. Those are MQ-4C or later
questions.

## MQ-4C Default-Off GEM Integration

MQ-4C connects both strict MQ-4B entry points to `GEMNewton` behind persistent,
default-off controls. `SetMQMQAHessianControls` accepts only finite alpha in
`[0,1]`; invalid calls leave the previous request unchanged. Reset restores
`enable=false` and `alpha=0`. A generic `MQMQAModelData` remains intentionally
invalid until its caller explicitly selects `SUBG` or `SUBQ`.

The atomic setter contract is exercised separately for negative alpha, alpha
greater than one, NaN, positive infinity, and negative infinity. Every invalid
call preserves the previously accepted enable flag and alpha exactly.

For every assembled Newton system, the live router inspects active solution
phases and dispatches exact model types only:

```text
SUBG -> BuildMQMQAGEMCorrection
SUBQ -> BuildMQMQASUBQGEMCorrection
other model types -> ignored
```

Charged phases are deliberate exclusions. `MQMQA_MAP_OUTSIDE_INTERIOR`
separates a well-formed state with `min(x)<=1E-12` from malformed or nonfinite
`MQMQA_MAP_INVALID_INPUT` data; the threshold itself is unchanged. A successful uncharged phase adds
one complete `(deltaA,deltaB)` pair to temporary aggregate arrays. Any other
failure after strict routing erases the entire aggregate, so a partially
corrected global system cannot survive. Controlled sequential correction-pair
additivity verifies this aggregation algebra; it is not native evidence for an
assessed equilibrium containing multiple simultaneous MQMQA phases.

### Historical baseline and linear transaction

The historical baseline is the standard, fully assembled Thermochimica `A/B`
after existing charged-phase safeguards and before any MQMQA correction. Since
`DGESV` overwrites its arguments, the corrected solve receives private trial
copies only:

```text
Atrial = Abase + alpha*deltaA
Btrial = Bbase + alpha*deltaB.
```

The original baseline arrays remain untouched until the trial has passed
application validation, `DGESV`, and a finite-update check. Application,
corrected-solve, or nonfinite-update failure solves fresh copies of the
historical baseline. Alpha zero bypasses construction and application entirely,
which gives bit-for-bit baseline `A/B`, solve, and final-output behavior.

Alpha is therefore a fixed correction weight in MQ-4C. It multiplies the
completed reduced matrix and residual corrections together; it is not a scale
factor on the analytic Hessian. Adaptive alpha selection, nonlinear merit
checks, and recovery from poor full-curvature steps are handled separately by
MQ-4D below.

### RKMP ownership

Enabling both experimental controls is not itself a conflict. A conflict exists
only when an eligible active RKMP phase and a successfully built nonzero-alpha
MQMQA aggregate occur in the same solve. The untouched historical `A/B` is then
passed directly to `SolveRKMPAlphaTrust`; the MQMQA aggregate is discarded and
no MQMQA trial reaches RKMP trust. The four dual-control cases are checked:
only RKMP eligible, only MQMQA eligible, neither eligible, and both eligible.
The conflict case uses controlled RKMP ownership around a real SUBQ state
because no assessed simultaneous RKMP/MQMQA assemblage is available. Its
resulting update agrees exactly with the equivalent RKMP-owned solve, and no
MQMQA trial is captured. This does not claim native multiphase coverage.

### Diagnostics and verification boundary

The integration metrics distinguish control request, supported phases,
successful phase corrections, completed aggregates, applications, accepted
corrected solves, and fallbacks. They separately count aggregate, application,
`DGESV`, nonfinite-update, charged-exclusion, strict-interior, and RKMP-conflict
outcomes. The interior metric also retains the minimum rejected fraction. With
the element block/residual denoted by subscripts `ee/e`, the reported ratios are

```text
rhoA = ||alpha*deltaA||F / max(||Abase,ee||F,1E-30)
rhoB = ||alpha*deltaB||2 / max(||Bbase,e||2,1E-30).
```

`TestMQMQAGEMIntegration.F90` verifies control persistence and reset, exact
alpha-zero identity, strict live SUBG and SUBQ routing, state-free pair
additivity and aggregate erasure, charged exclusion, application/`DGESV`/
nonfinite fallback, and the RKMP ownership rule. At a controlled positive
interior FeTiVO state, the captured corrected pre-solve system differs from the
captured historical baseline only by the requested element-block and element-
residual corrections. The maximum application discrepancy is `9.7700E-15`,
and the alpha-one corrected linear solve returns a finite update with `INFO=0`.
Independent solves of the captured baseline and corrected systems also require
a numerically resolved change in the pre-line-search GEM linear-system solution
vector, rather than accepting array insertion alone as integration evidence.
The relative `L2` separation between the corrected and historical pre-line-
search GEM solution vectors, normalized by the larger solution-vector norm, is
`2.6485E-01`. This demonstrates a nontrivial effect on the solved vector; it is
not an error measure or a claim of a 26.5 percent physical-output change. The
report also separates magnitude and orientation by printing historical/
corrected norms of `4.3514E+01` and `4.5485E+01`, a corrected-to-historical norm
ratio of `1.0453E+00`, cosine similarity of `9.6432E-01`, and an angle of
`1.5352E+01` degrees in the current unscaled coordinate representation. Thus
this state contains both a modest solution-magnitude change and a resolved
coordinate-space rotation. `GEMLineSearch` may subsequently limit the applied
thermodynamic-state displacement, so this angle does not describe that final
applied displacement. The corrected replay reproduces the solution vector
accepted by the live linear path exactly within the reported double-precision
comparison.
Just-below and just-above threshold controls freeze the dedicated interior
status at the unchanged `1E-12` boundary.

The separate full FeTiVO fixed-alpha-one run is evidence rather than an MQ-4C
gate. The observed run converges in 95 iterations with no reversions, 82
accepted corrected linear solves, and no application, `DGESV`, or nonfinite
fallback. All 16 early fallbacks are classified solely as
`MQMQA_MAP_OUTSIDE_INTERIOR`, with zero malformed aggregate failures; the
smallest rejected fraction is `1.7504E-19`. Relative to the historical
run, the scaled final differences are `3.4540E-13` for total Gibbs energy,
`5.0564E-09` for the complete mole-fraction array, `2.0179E-09` for species
moles, and `8.3206E-10` for phase amounts. This is encouraging evidence, but
the small final-state differences do not establish that the correction is
physically important in this particular FeTiVO calculation. Instead, the
non-negligible linear-system ratios (`rhoA=3.5183E-02` and `rhoB=7.2634E-03`)
show that the correction was genuinely active while both paths converged to
nearly the same equilibrium. The 16 strict-interior outcomes are intentional
eligibility exclusions, not silent numerical or decoding failures.

The pre-line-search solution-vector comparison uses an ordinary Euclidean norm
over the complete GEM unknown vector. That vector combines element-potential and phase-related
unknowns with different physical meanings and numerical scales. The result is
therefore strong non-vacuity evidence, but it is not a solver-safety or trust
metric. Variable-group scaling and actual nonlinear residual or merit-function
progress are addressed by the MQ-4D globalization layer below.

The MQ-4C claim is therefore **live default-off integration demonstrated and
transactionally verified**, not finished MQMQA solver integration. MQ-4C does
not promise alpha-one robustness across calculations. Native
evidence remains limited to assessed plain-SUBG cases and the uniform-zeta
FeTiVO SUBQ `G/Q` case, plus a controlled modified-runtime production regression
of the corrected nonuniform-zeta S3 selection. SUBQ `B`, `R`,
magnetism, assessed database-native nonuniform-zeta data, native simultaneous multi-MQMQA
assemblages, and combined RKMP/MQMQA corrections remain outside this claim.

## MQ-4D: default-off adaptive globalization

### Reader orientation: what MQ-4D changes

MQ-4D changes how the already verified MQMQA correction is admitted into a
live Newton solve.  It does not derive another Hessian and it does not replace
Thermochimica's historical GEM solver.  At an eligible solver state, the
historical code has assembled a linear system

```text
A_base u_base = B_base,
```

where `A_base` is the historical GEM Newton matrix, `B_base` is its
right-hand-side/residual vector, and `u_base` is the resulting Newton update.
The MQMQA mapper independently supplies a completed correction pair
`(deltaA,deltaB)`.  An alpha candidate therefore solves

```text
(A_base + alpha*deltaA) u_alpha = B_base + alpha*deltaB.
```

The important meanings are:

- `alpha=0` is the untouched historical GEM Newton solve.  It is **not** a
  first-order method and it does not mean that Thermochimica has stopped using
  all thermodynamic derivatives.
- `alpha=1` applies the complete mapped MQMQA excess-curvature correction and
  is the desired candidate whenever it is safe.
- `0<alpha<1` applies the same completed correction pair more cautiously.  It
  does not scale or alter the local analytic Hessian itself.
- The accepted Newton direction is still passed to the existing GEM line
  search.  MQ-4D adds a safety decision before that line search; it does not
  replace it.

The objective is reliable nonlinear convergence with sustained full-alpha use
near the solution, not the largest possible count of `alpha=1` selections.
More aggressive full-alpha use can be counterproductive if it repeatedly
pushes the nonlinear iteration away from a settled region.

### One eligible solve in plain language

For each GEM solve at which an active, supported MQMQA phase is eligible,
`SolveMQMQAAlphaTrust` performs the following sequence:

1. Build the complete phase-aggregated `(deltaA,deltaB)` pair once.  If any
   applicable phase fails, discard the whole pair rather than applying a
   partial correction.
2. Solve an untouched copy of the historical system.  This provides both the
   exact fallback and the reference update used by the trust checks.
3. Ask whether the current phase assemblage and recent nonlinear progress are
   settled enough to try a nonzero correction.  If not, return the historical
   solution exactly.
4. If ready, try alpha candidates from largest to smallest, beginning with
   `alpha_max` (normally one).
5. For each candidate, apply both corrections together, solve the trial linear
   system, and compare its update with the historical update in separately
   scaled variable groups.
6. Accept the first safe candidate.  If every positive candidate is rejected,
   return the already solved `alpha=0` result.
7. Let the existing GEM line search process the selected update.  On the next
   nonlinear iteration, use the observed residual and Gibbs-energy progress to
   retain or revoke readiness.

Readiness and candidate safety answer different questions.  **Readiness asks
whether this point in the nonlinear trajectory is mature enough to attempt the
correction.**  **Candidate safety asks whether one particular corrected linear
solve is acceptably close to the historical solve.**  A solve can therefore be
ready while rejecting `alpha=1` and accepting a smaller positive alpha.

### Code-to-concept map

| Location | Responsibility |
|---|---|
| `ModuleMQMQAResponseMapping.f90` | Builds the phase-local and aggregated `(deltaA,deltaB)` correction without changing the live GEM arrays. |
| `ModuleMQMQATrust.f90` | Constructs descending alpha candidates and performs state-free correction-size and grouped-update checks. |
| `GEMNewton.f90::SolveMQMQAAlphaTrust` | Orchestrates the baseline solve, readiness decision, candidate trials, exact fallback, and diagnostic history. |
| `GEMNewton.f90::BuildMQMQAGroupedDisplacement` | Converts the mixed GEM solution vector into three meaningful update groups before norms and directions are compared. |
| `SetMQMQAHessianControls.f90` | Exposes the narrow default-off fixed and adaptive control APIs. |
| `ModuleGEMSolver.f90` | Stores controls, heuristic limits, readiness state, counters, and alpha/rejection histories. |
| `InitGEMSolver.f90` | Resets adaptive runtime state so one calculation cannot leak history into the next. |
| `TestMQMQAAdaptiveTrust.F90` | Checks the utilities, controls, failure behavior, live FeTiVO path, and bounded sensitivity matrix. |

### Terminology used in the diagnostics

| Term | Meaning |
|---|---|
| Historical or baseline system | The GEM matrix and right-hand side that Thermochimica would solve with MQMQA curvature integration disabled. |
| Correction pair | The completed `deltaA` and `deltaB`; they are constructed and applied together. |
| Eligible solve | A call for which an active supported MQMQA phase is present and the correction builder can be considered.  It is not a configured maximum or a convergence limit. |
| Readiness activation/reset | The adaptive state begins/stops allowing positive-alpha trials based on assemblage stability and measured nonlinear progress. |
| Full/reduced/zero alpha | Selected `alpha=1`, `0<alpha<1`, or `alpha=0`, respectively. |
| Fallback | Returning the untouched historical system because the correction was inapplicable, invalid, unsafe, or could not be solved. |
| Final full-alpha window | Consecutive eligible solves ending at the last eligible solve before convergence, all using `alpha=1` without a readiness reset. |
| `First+` / `First1` | First global iteration selecting positive alpha / first iteration of the final sustained full-alpha window. |
| `F/R/Z` | Counts of eligible solves selecting full, reduced, and zero alpha. |
| Function norm | The existing GEM nonlinear residual measure used to judge whether the state is locally resolved. |
| Gibbs gap | A scaled comparison with the best Gibbs energy observed for the current settled assemblage. |
| Transactional | Either the complete correction pair is accepted, or the original arrays are retained exactly; no partial phase correction survives a failure. |

An eligible-solve count need not equal the printed main Newton iteration count.
For the FeTiVO evidence below, the 77 eligible solves consist of one
phase-assemblage/initialization solve recorded at global iteration zero plus 76
main-loop solves.  This is a property of where `GEMNewton` is called, not a
77-iteration limit; the normal Thermochimica iteration limit remains unchanged.

### Implementation and controls

MQ-4D preserves the complete MQ-4C correction and transaction.  It builds one
aggregate `(deltaA,deltaB)` pair, solves an untouched alpha-zero copy, and then
tests completed correction pairs in descending order:

```text
[alpha_max, 0.1, 0.01, 0.001, 0]
```

Duplicates and candidates above `alpha_max` are omitted.  Alpha always weights
the completed reduced GEM matrix and residual corrections together; it never
weights the local SUBG or SUBQ Hessian.  If no positive candidate passes, the
already solved alpha-zero arrays are returned exactly.  The historical GEM line
search is unchanged.

The persistent adaptive interface is separate from the fixed-alpha interface:

```text
SetMQMQAHessianAdaptiveControls(enable,alpha_max,info)
ResetMQMQAHessianAdaptiveControls()
```

Both modes remain default off.  A valid setter call selects its mode and clears
the other request.  Nonfinite or out-of-range `alpha_max` values are rejected
before any persistent setting changes.  Thus fixed and adaptive ownership is
unambiguous while the MQ-4C API and evidence remain available unchanged.

### Readiness and candidate safety

The additional MQMQA correction is withheld while the phase assemblage or
nonlinear state is not locally settled; the historical GEM Newton solve remains
active.  Readiness uses an established Gibbs minimum, at least five
iterations since the last assemblage change, a local function norm, and actual
next-iteration residual and Gibbs progress.  The residual-ratio check uses the
existing `1.05` allowance while the function norm exceeds `1E-6`; below the
existing GEM line-search small-norm threshold of `1E-6`, the absolute norm is
treated as resolved and a ratio of two sub-threshold residuals does not revoke
readiness.  Readiness activations, resets, function norms, norm ratios, and
relative Gibbs gaps are retained per iteration.

Each positive candidate must also pass finite correction/application checks,
`DGESV`, an emergency correction-to-baseline cap of `1E6`, and groupwise update
checks relative to the alpha-zero solve.  The three groups are:

1. dimensionless element-potential displacements;
2. constituent-level first-order logarithmic mole increments for active
   solution phases;
3. pure-phase amount displacements in moles.

The second group is reconstructed from the exact expression subsequently used
by `GEMLineSearch`; the solution-phase scalar GEM unknown is not incorrectly
treated as a composition increment.  The ratio, cosine, and relative-direction
limits are `1.25`, `0.90`, and `0.50`, respectively.  Ratios and directions of
two negligible steps are neutral below group-specific resolved-step floors:
`1E-8` for the two dimensionless groups and `1E-12 mol` for pure-phase amounts.
This prevents a large ratio of two numerically immaterial steps from becoming a
false safety rejection without relaxing the stated ratio limit.

These limits and floors are numerical globalization heuristics, not
thermodynamic identities.  Their use must be accompanied by sensitivity and
broader assessed-database evidence before any general robustness claim.

### FeTiVO exit evidence

`TestMQMQAAdaptiveTrust.F90` verifies candidate ordering, emergency-ratio,
groupwise magnitude and direction rejection, negligible-step handling,
mutually exclusive persistent controls, atomic invalid-call rejection, exact
adaptive-alpha-zero identity, and the live FeTiVO nonlinear path.  The existing
MQ-4C test separately retains no-applicable-phase, boundary-state, application,
singular-solve, and nonfinite fallback coverage against untouched baseline
arrays.

For FeTiVO with `alpha_max=1`, the adaptive run converges in 76 iterations.  It
records 77 eligible solves: 13 at full alpha, 50 at reduced alpha, and 14 at
alpha zero.  The final 13 eligible corrected solves all select `alpha=1`, with
no readiness reset inside that final window.  Early reductions are accompanied
by rejection masks for every larger candidate.  The final scaled differences
from the historical result are approximately `1.34E-15` for Gibbs energy,
`4.81E-10` for mole fractions, `5.88E-10` for species moles, and `8.33E-10` for
phase amounts. Candidate-by-candidate alpha and rejection-mask histories make
each larger rejected candidate individually identifiable; the reduced-alpha
gate does not rely only on an undifferentiated per-iteration mask.

### FeTiVO trust-setting sensitivity

The sensitivity driver reruns the same FeTiVO inputs against one historical
default-off reference.  It changes one internal trust setting at a time, then
restores every production setting exactly before the next run.  The two final
rows are exploratory combined policies; neither changes the production
defaults.  `First+` is the first global iteration selecting positive alpha,
`First1` begins the final sustained alpha-one window, and `F/R/Z` counts full,
reduced, and zero-alpha eligible solves.

The classifications are intentionally behavioral rather than thermodynamic:

- **ROBUST:** converged, agreed with the historical final state within the
  established `1E-8` scaled tolerance, retained at least three final alpha-one
  solves, and did not differ materially from the default trajectory;
- **SAFE BUT SENSITIVE:** preserved all safety and final-state requirements but
  changed iterations or first-positive activation by at least five, changed
  the final alpha-one window by at least three, changed readiness-reset count,
  or did not retain the three-solve final alpha-one window;
- **FAILED:** violated finiteness, alpha bounds, candidate-specific rejection
  evidence, fallback safety, setting restoration, or historical agreement.

| Case | Changed value | Class | Iter | Eligible | First+ | First1 | F/R/Z | Final window | Act/reset |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| Production default | defaults | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Local residual lower | 0.025 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Local residual upper | 0.10 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Settled period conservative | 8 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Settled period permissive | 3 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Gibbs activation lower | 1E-7 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Gibbs activation upper | 1E-5 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Gibbs retention lower | 1E-5 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Gibbs retention upper | 1E-3 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Progress allowance lower | 1.02 | SAFE BUT SENSITIVE | 64 | 65 | 10 | 59 | 12/31/22 | 6 | 13/12 |
| Progress allowance upper | 1.10 | SAFE BUT SENSITIVE | 74 | 75 | 10 | 63 | 12/51/12 | 12 | 3/2 |
| Update ratio lower | 1.10 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Update ratio upper | 1.50 | SAFE BUT SENSITIVE | 103 | 104 | 10 | 85 | 69/3/32 | 19 | 23/22 |
| Direction cosine conservative | 0.95 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Direction cosine permissive | 0.85 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Relative direction lower | 0.35 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Relative direction upper | 0.75 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Earlier activation policy | 0.10 / 3 / 1E-5 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |
| Later activation policy | 0.025 / 8 / 1E-7 | ROBUST | 76 | 77 | 10 | 64 | 13/50/14 | 13 | 5/4 |

All cases converged with `INFOThermo=0`, selected alpha only in `[0,1]`,
reached maximum alpha one, documented every larger rejected candidate whenever
a reduced alpha was selected, restored the production settings, and remained
within the historical final-state tolerance.  No correction, ratio,
linear-solve, or nonfinite failure occurred.  Across the matrix, scaled final
differences were at most `2.30E-15` in Gibbs energy, `1.23E-9` in constituent
mole fractions, `6.26E-10` in species moles, and `8.42E-10` in phase amounts.

The exploratory earlier policy did **not** activate curvature earlier: it was
identical to the default in activation iteration, total iterations, alpha
counts, readiness resets, final alpha-one window, and final differences.  The
later policy was also identical.  On this calculation the changed local,
settling, and activation limits are therefore not the controlling readiness
conditions.

The progress allowance is behaviorally active in both neighboring directions.
At `1.02`, convergence shortened to 64 iterations but readiness cycled 13/12
times and the final alpha-one window shortened to six.  At `1.10`, convergence
required 74 iterations with 3/2 readiness activation/reset events.

With the update-ratio limit increased to `1.50`, full alpha was accepted at
iteration 10, but the following residual norm increased by a factor of `2.92`
and readiness was revoked.  Similar full-alpha/recovery cycles produced 23
activations and 22 resets, increasing convergence from 76 to 103 iterations.
The calculation remained safe and reached the historical equilibrium,
demonstrating successful fallback, but the permissive threshold produced
inferior nonlinear progress.  The production value of `1.25` was therefore
retained.

The current defaults are retained.  They lie within a tested safe and
convergent neighborhood, while the default progress and update-ratio values
avoid the material trajectory changes observed at neighboring permissive or
stricter settings.  This is evidence for retaining the current FeTiVO policy,
not evidence that the defaults are globally optimal.

> The study evaluates local sensitivity of the MQ-4D heuristics for FeTiVO. It
> does not establish universal convergence or robustness for arbitrary MQMQA
> assessments.

The fixed-alpha-one MQ-4C evidence remains unchanged: it converges in 95
iterations with 82 accepted corrected solves and only 16 intentional
strict-interior exclusions.  MQ-4D therefore establishes default-off adaptive
connection and a sustained settled full-alpha window for this FeTiVO case.  It
does not establish universal alpha-one robustness, broad model/database
coverage, constrained-KKT compatibility, or assessed molten-salt application
evidence.

Private local MSD-TC FLiBe verification now supplies database-native
nonuniform-zeta `G`-family scalar, gradient, and Hessian evidence over three
documented compositions.  That closes the local derivative-evidence item, but
it does not extend the current MQ-4D solver claim: the FLiBe cases have not yet
been exercised through live GEM mapping or adaptive globalization.  MQ-4E must
reuse the documented database and application range to record eligibility,
strict-interior and boundary exclusions, accepted corrections, alpha histories,
fallback frequency, convergence, and final-equilibrium agreement.  `B` remains
controlled standalone coverage rather than database-native solver evidence.

## MQ-4E-A: portable public-database solver coverage

### Purpose and evidence boundary

MQ-4E-A broadens the live solver experiment without using the private MSD-TC
files reserved for MQ-4E-B.  `TestMQMQASolverCoverage.F90` runs thirteen
public-data states: the repository's `CuFeC-Kang.dat` plain-`SUBG` case at
1400 K and a twelve-state `FeTiVO.dat` `SUBQ` sweep.  The FeTiVO sweep changes
temperature from 1900 to 2100 K while Fe/Ti change linearly from 0.55/0.45 to
0.45/0.55 mol; V=0.5 mol, O=2 mol, and pressure=1 atm remain fixed.  It is a
controlled parameter trajectory, not a physical time integration.

Each state is solved three times from a fresh Thermochimica state:

1. historical/default-off GEM, which supplies the case-local reference;
2. fixed `alpha=1`, retained as a stress comparison with the MQ-4C interface;
3. adaptive `alpha_max=1`, which is the MQ-4D production candidate.

The test records convergence, iterations, eligible/application/accepted-solve
counts, full/reduced/zero-alpha selections, final full-alpha window, boundary
exclusions, correction ratios, readiness changes, rejection reasons, and
scaled final differences.  It also reports fixed/historical and
adaptive/historical iteration ratios so safety evidence is not mistaken for a
performance improvement.  Every reduced-alpha selection must retain a
candidate-specific rejection mask for each larger rejected candidate.  The
default-off run must remain silent: no correction, eligibility, or boundary
diagnostic is allowed to activate.

The classifications deliberately separate what was demonstrated:

- `FULL_CURVATURE_EVIDENCE`: at least one full-alpha correction was accepted;
- `REDUCED_CURVATURE_EVIDENCE`: positive corrections were accepted, but none
  at full alpha;
- `SAFE_BOUNDARY_FALLBACK`: only documented strict-interior exclusions were
  observed;
- `SAFE_TRUST_FALLBACK`: the adaptive path considered eligible states but
  retained only the historical update;
- `REFERENCE_STATE_DIFFERENCE`: the run converged finitely with accepted
  corrections, but its raw final state arrays did not agree with the
  case-local historical arrays within the stated comparison tolerance;
- `FAILURE`: convergence, finiteness, transactional safety, or documented
  reduced-alpha evidence failed.

Only the first two classifications are curvature-effect evidence.  The two
`SAFE_*` classes are fallback evidence.  `REFERENCE_STATE_DIFFERENCE` is a
flag for further interpretation, not a successful final-state-equivalence
claim.  Historical comparison uses scaled tolerances of `1E-8` for total Gibbs
energy and `1E-6` for mole fractions, species moles, and phase amounts.  The
state tolerance is still substantially tighter than the `1E-3` relative
application checks in Tests 56 and 57, while allowing the broader trajectory's
solver-level numerical variation.

### Measured public-data results

All thirteen historical runs and all twenty-six corrected runs converged with
finite outputs.  No aggregate-construction, application, corrected-linear-
solve, or nonfinite-update failure occurred.  Every reduced-alpha decision was
candidate-specifically documented.

For adaptive mode, all thirteen states agreed with their historical references
within the stated tolerances.  Eleven states selected `alpha=1` at least once.
FeTiVO steps 4 and 11 selected only reduced positive corrections, with maximum
alpha `0.1`; they are therefore reduced-curvature evidence rather than full-
curvature evidence.  Three fixed FeTiVO runs encountered strict-interior
exclusions, but still accepted positive corrections elsewhere and converged.

| Public case | Historical iterations | Adaptive iterations | Adaptive/historical | Full/reduced/zero | Maximum alpha | Final full-alpha window | Adaptive class |
|---|---:|---:|---:|---:|---:|---:|---|
| CuFeC 1400 K, SUBG | 448 | 657 | 1.467 | 26/271/306 | 1.0 | 0 | full-curvature evidence |
| FeTiVO step 1 | 45 | 69 | 1.533 | 2/42/28 | 1.0 | 1 | full-curvature evidence |
| FeTiVO step 2 | 66 | 116 | 1.758 | 18/62/27 | 1.0 | 18 | full-curvature evidence |
| FeTiVO step 3 | 45 | 67 | 1.489 | 1/40/29 | 1.0 | 0 | full-curvature evidence |
| FeTiVO step 4 | 46 | 66 | 1.435 | 0/44/25 | 0.1 | 0 | reduced-curvature evidence |
| FeTiVO step 5 | 87 | 107 | 1.230 | 19/31/60 | 1.0 | 18 | full-curvature evidence |
| FeTiVO step 6 | 79 | 95 | 1.203 | 18/27/53 | 1.0 | 15 | full-curvature evidence |
| FeTiVO step 7 | 38 | 40 | 1.053 | 7/24/10 | 1.0 | 7 | full-curvature evidence |
| FeTiVO step 8 | 41 | 67 | 1.634 | 30/30/8 | 1.0 | 30 | full-curvature evidence |
| FeTiVO step 9 | 33 | 73 | 2.212 | 1/58/15 | 1.0 | 0 | full-curvature evidence |
| FeTiVO step 10 | 32 | 81 | 2.531 | 1/64/17 | 1.0 | 0 | full-curvature evidence |
| FeTiVO step 11 | 32 | 84 | 2.625 | 0/66/19 | 0.1 | 0 | reduced-curvature evidence |
| FeTiVO step 12 | 31 | 52 | 1.677 | 2/39/12 | 1.0 | 0 | full-curvature evidence |

The public plain-SUBG case is a clear performance warning.  Adaptive mode
preserved the historical final state, but required 657 iterations versus 448
historically, an iteration ratio of `1.467`.  It recorded 159
readiness activations and 158 resets.  Its rejections were dominated by update
magnitude (`551`), followed by nonlinear readiness (`305`) and direction
(`11`).  This is safe guarded behavior, not evidence that adaptive curvature
improves convergence for CuFeC.  It should remain visible in MQ-4E-C's final
performance interpretation rather than prompting FeTiVO-tuned threshold
changes during MQ-4E-A.

The performance issue is not limited to CuFeC.  Adaptive mode required more
iterations than the historical solver in every public case, with ratios from
`1.053` to `2.625`.  This matrix therefore demonstrates safety and coverage,
not acceleration.  MQ-4E-C must use the public and private MQ-4E evidence to
decide whether the overhead is an acceptable cost of guarded experimental
curvature, whether readiness/trust needs another bounded revision, or whether
adaptive curvature should remain an opt-in diagnostic capability.  Thresholds
are not changed in MQ-4E-A solely to improve this table.

Fixed `alpha=1` on CuFeC converged in 466 iterations and applied 383 corrected
solves, but the raw species/phase state arrays differed materially from the
historical arrays even though the scaled Gibbs-energy difference was only
`5.66E-14`.  MQ-4E-A therefore labels this result
`REFERENCE_STATE_DIFFERENCE`; it does not call the fixed run equivalent or use
it as the adaptive acceptance gate.  Whether the raw difference represents an
alternative near-degenerate assemblage or a meaningful fixed-curvature solver
difference requires phase-identity-aware interpretation in MQ-4E-C.

MQ-4E-A establishes portable public-data execution through both supported
plain-SUBG and SUBQ routes, with exact default-off behavior, transactional
safety, and adaptive historical agreement over the tested states.  It does not
establish that full alpha is selected near every solution, that adaptive mode
is faster, or that the FeTiVO trust heuristics are optimal for plain SUBG.
Private assessed molten-salt solver evidence remains MQ-4E-B.

## MQ-4E-B: private assessed FLiBe solver gate

### Purpose and private-evidence boundary

MQ-4E-B reused the local MSD-TC V4.1 fluoride assessment at 1000 K, 1 atm,
and LiF/BeF2 feed ratios 45/55, 50/50, and 55/45.  The database, optional
driver, runner, and detailed numerical record remain ignored local assets
because the external assessment is not distributed with Thermochimica.  The
tracked audit records the result and claim boundary without publishing the
database or database-dependent test source.

Unlike the earlier derivative fixture, the live solver experiment did not
blend quadruplet amounts toward an artificial interior state.  Each actual
assessed equilibrium was run as an untouched historical reference, a second
untouched reproducibility check, fixed `alpha=1`, and adaptive
`alpha_max=1`.  The repeated historical calculations were bitwise identical
in the reported Gibbs and state-array comparisons and reproduced the same
`MSFL` assemblage and iteration count.  The reference is therefore stable.

### Measured blocking result

All corrected calculations returned `INFOThermo=0`, but none satisfied the
established historical final-state gate.  The important results are:

| LiF/BeF2 | Historical | Fixed `alpha=1` | Adaptive `alpha_max=1` |
|---|---|---|---|
| 45/55 | 1267 iterations; `MSFL` | 1117; solid + liquid; scaled `dG=5.120E-3` | 2884; `MSFL`; scaled `dG=2.271E-7`, `dx=4.339E-3`; F/R/Z=30/79/2124 |
| 50/50 | 1586 iterations; `MSFL` | 5510; `gas_ideal`; scaled `dG=4.619E-2` | 4720; `gas_ideal`; scaled `dG=4.619E-2`; F/R/Z=37/72/1134 |
| 55/45 | 3511 iterations; `MSFL` | 5459; solid + liquid; scaled `dG=1.174E-2` | 1002; solid + liquid; scaled `dG=4.668E-3`; F/R/Z=16/14/455 |

The alternative assemblages have less negative, and therefore higher, Gibbs
energies than the reproducible historical `MSFL` states.  These are not merely
different representations of a near-degenerate equilibrium.  Fixed mode also
recorded 18, 10, and 27 corrected-linear-solve fallbacks, respectively.  Those
fallbacks were transactional, but accepted corrections earlier in the
trajectory had already changed the subsequent phase search.

Adaptive mode recorded no aggregate-construction, application,
accepted-path linear-solve, or nonfinite fallback, and every reduced-alpha
selection retained a rejection reason for all larger candidates.  It still
ended with no final full-alpha window and frequent readiness cycling:
55/54, 77/77, and 22/22 activation/reset events.  Local step acceptance and
next-iteration progress checks therefore did not guarantee recovery of the
historical global phase assemblage for these assessed states.

The separate database-native derivative gate remains successful.  It still
shows nonuniform-zeta `G`-family scalar, gradient, Hessian, and production-mu
`H*v` consistency for all three compositions.  MQ-4E-B therefore distinguishes
a global phase-path/globalization defect from a failure of the local analytic
derivatives or mapper algebra.

The hardening investigation added two further localization checks.  First, the
actual converged `MSFL` tangent Hessians were positive definite at all three
states.  Their minimum/maximum tangent eigenvalues were
`1.608/1.164E5`, `2.428/5.883E4`, and `2.194/2.658E4`, respectively.
Second, the independent MQ-4A reduced-mapping finite-difference experiment was
repeated without the earlier 20 percent interior blend.  All three unmodified
assessed states passed live baseline capture, builder-versus-independent
`deltaA/deltaB` comparison, and second-order finite-difference checks of
`deltaA`, `deltaB`, and their combined affine action.  The worst reported
scaled mapping error remained below `2E-9`.  Thus the evidence now localizes
the unresolved behavior beyond the local Hessian, constrained response, and
reduced GEM mapping layers.

Two bounded globalization experiments were also rejected rather than retained
in production.  Comparing the corrected and historical line-search outcomes
from the same pre-step state did not preserve the eventual historical phase
path: a locally preferable step can still enter a worse basin several phase
search decisions later.  Delaying curvature until Thermochimica's established
post-1000-iteration convergence-eligibility regime likewise failed at 55/45
and left small but gate-breaking state differences at 45/55 and 50/50.  These
experiments show that neither a one-step merit comparison nor a fixed global
iteration delay is an adequate phase-path safeguard.  Both experimental code
changes were removed.

The three established RKMP phase-path safeguards were then ported directly to
MQMQA as one further bounded experiment: live SUBG/SUBQ activity tracking,
suppression of the residual-only convergence shortcuts while MQMQA curvature
was active, and the larger mixed solution/pure-phase Wolfe-search budget.  The
same three unchanged FLiBe states still failed the historical-state exit gate.
At 45/55 and 50/50, fixed `alpha=1` reached the reported 6001-iteration limit;
at 55/45 it converged to `gas_ideal` rather than the historical `MSFL` state.
Adaptive mode ended in alternative solid/liquid or liquid/gas assemblages at
all three compositions.  Its readiness rejection counts remained substantial,
and no run recovered the historical final state.  The stricter convergence
path usefully exposed that some earlier alternative states were not globally
settled, but it did not restore the correct phase path.  Therefore, the missing
RKMP safeguards were relevant omissions but were not the complete MQMQA
solution.  The experimental production changes were removed after this gate.

The remaining failure was then localized with attempt-aware trajectory and
ablation diagnostics.  This distinction matters because Thermochimica can
restart the GEM solve after the nominal 3000-iteration limit; comparing only
the final attempt can otherwise assign a phase-path change to the wrong
iteration.  Four controls all reproduced the historical result exactly:

- a second unmodified historical solve;
- adaptive mode with `alpha_max=0`;
- entry into the complete adaptive correction-building and baseline-solve path
  while positive candidates were never tested; and
- normal candidate evaluation followed by diagnostic rejection of every
  otherwise accepted positive candidate.

Consequently, the divergence is not random solver variability, corruption from
constructing the correction, or a side effect of the adaptive bookkeeping.  It
requires an accepted positive curvature correction.

Aligned historical/corrected traces placed the first phase-assemblage
divergence at, or after, an accepted correction in every composition.  Skipping
only that first accepted event delayed some divergences but did not restore the
historical result: a later accepted correction produced the same qualitative
redirection.  This rules out a single anomalous Newton iteration and instead
shows a repeatable sensitivity of the discrete phase search to the corrected
continuous update.

One final causal probe separated ordinary GEMNewton calls from the speculative
assemblage solves issued by `CheckPhaseChange`.  The normal adaptive runs made
1179, 919, and 380 such probe calls at 45/55, 50/50, and 55/45, but only 5, 5,
and 2 of those calls accepted positive curvature.  Suppressing MQMQA curvature
in every speculative probe reduced those accepted-probe counts to zero without
changing the final failure classification.  The first accepted main-solve
corrections under that suppression were only `0.01`, `0.01`, and `0.001`, yet
the three calculations still left the historical path.  Curvature inside the
speculative probes is therefore not the primary cause, and the behavior is not
an `alpha=1` blow-up.

The assessed FLiBe states also contain no active `R`, `Q`, or `B` terms: their
active excess contribution is the verified nonuniform-zeta `G` family.  The
absence of an `R` branch therefore cannot make the target historical solution
unrepresentable.  Together with the positive-definite tangent Hessians and the
successful unmodified-state mapping finite differences, the evidence supports
the following narrower diagnosis: a locally admissible and algebraically
correct continuous Newton correction can be large enough in phase-selection
coordinates to change a later discrete assemblage decision in this difficult,
multi-basin FLiBe solve.  The present local trust gates do not measure that
longer-horizon phase-path consequence.

That diagnosis was then tested against the independently supplied legacy
MSTDB-TC V3.1 fluoride assessment.  The same three 1000 K LiF/BeF2 inputs used
one complete database version per fresh process.  V3.1 changed the wider
database topology relative to V4.1, but the active binary `MSFL` state still
decoded as a six-quadruplet SUBQ phase with five `G` terms, no active `R`, `Q`,
or `B` terms, and zeta values spanning 2.4 to 6.0.  Historical and repeated
historical V3.1 solves agreed exactly at all three compositions, with final
`MSFL` assemblages after 1526, 902, and 2700 iterations.  Force-zero and
evaluate-then-reject controls also reproduced those references exactly.

Retained positive curvature nevertheless failed the same reference-state gate
in V3.1.  Adaptive runs required 2657, 1637, and 1259 iterations and departed
from the historical phase path or final state.  At 45/55, the first retained
positive correction was followed immediately by the first recorded divergence:
the corrected path added `Be_S1(s)` while the historical path retained only
`MSFL`.  The V3.1 production-linked nonuniform-zeta derivative evidence still
passed, with a worst best `H*v` error of `4.87E-10` and the expected
second-order region.  Thus the V4.1 reassessment is not the sole cause, and the
common failure is not explained by missing `R`, `Q`, or `B` families or by an
incorrect local Hessian.

The next diagnostic captured the exact normalized pure- and solution-phase
driving forces compared by `CheckPhaseAssemblage`, together with the phase-set
difference at the first historical/corrected divergence.  A negative force is
eligible to enter the assemblage; the more negative pure/solution value is
considered first.  The measured adaptive results were:

| Database and LiF/BeF2 | First phase-set difference | Historical leading pure / solution | Corrected leading pure / solution | Corrected pure-minus-solution |
|---|---|---:|---:|---:|
| V4.1 45/55 | add `Be_S1(s)` | `0 / 0` | `-45.613 / -16.572` | `-29.041` |
| V4.1 50/50 | remove `gas_ideal` | `0 / 0` | `0 / 0` at the nearest add check | `0` |
| V4.1 55/45 | add `Be_S1(s)` | `0 / 0` | `-13.019 / 0` | `-13.019` |
| V3.1 45/55 | add `Be_S1(s)` | `0 / -30.643` | `-18.480 / 0` | `-18.480` |
| V3.1 50/50 | add `Be_S1(s)` | `0 / -127.197` | `-20.602 / 0` | `-20.602` |
| V3.1 55/45 | remove `Be_S1(s)` | `-10.966 / 0` | `0 / -6.963` | `+6.963` |

The V4.1 45/55 and 55/45 additions occur at the same iteration as a retained
full-alpha correction.  V3.1 45/55 diverges on a retained reduced correction
(`alpha=0.001` at the recorded decision), demonstrating that this is not only
an `alpha=1` instability.  V3.1 50/50 already differs during an outer
zero-alpha iteration because earlier internal accepted curvature activity has
changed the state delivered to that iteration.  The V4.1 50/50 event is a
solution-phase removal, whose amount-based removal criterion is not represented
by the phase-addition driving-force capture; its nearest add check is therefore
reported as inconclusive rather than interpreted as a tie.

Those aligned-trajectory values identify what the two already-separated paths
eventually presented to `CheckPhaseAssemblage`; they do **not** by themselves
measure the immediate effect of one corrected Newton candidate.  A subsequent
same-pre-step replay therefore solved the untouched `alpha=0` and selected
positive-alpha systems from one identical Thermochimica state, without
committing either result.  It recorded the element-potential targets, rebuilt
the pure-phase force component by component, and reran the production
inactive-solution ranking.  The componentwise pure-force reconstruction agreed
with the direct production formula to roundoff.

This causal replay separates three behaviors that the earlier trajectory-only
table combined:

- In V4.1 45/55, the candidate at the later `Be_S1(s)` addition changed the
  element potentials by at most `1.61E-3`, while the immediate pure and
  solution rankings remained at roundoff.  The order-10 phase forces seen six
  iterations later were therefore accumulated nonlinear/phase-search effects,
  not the instantaneous response of the same candidate.
- In V4.1 50/50, the selected `alpha=0.01` candidate changed the Be/Li element
  potentials by only about `2.66E-6`.  The corrected trajectory nevertheless
  reached the exact solution-removal test with `gas_ideal=5.121E-12`, below the
  `1E-11` threshold, while the historical trajectory still contained
  `gas_ideal=2.232E-1` at the same global iteration.  The removal routine thus
  behaved as coded; the continuous trajectories had separated before the
  phase-set change.
- Other nearby recorded gas removals also occurred on the historical path and
  are explicitly labelled non-causal by the diagnostic.  A removal event is
  treated as explanatory only when the matched historical trace still retains
  that same phase at the same global iteration.
- In V4.1 55/45, the first retained full-alpha candidate coincided with the
  `Be_S1(s)` addition, but both replayed pure-phase forces were negative only
  at roundoff (`-4.44E-16` and `-8.88E-16`).  This is a discrete boundary tie,
  not an order-10 immediate force displacement.
- V3.1 contains both regimes.  Its 45/55 same-state pure-force change was only
  `4.00E-15`; at 50/50 the full-alpha candidate changed the element-potential
  targets by as much as `35.51` and the inactive-solution force from `-76.56`
  to `-49.92`; at 55/45 the `alpha=0.1` candidate changed the `Be_S1(s)` force
  from `-8.30` to `-15.16` and the leading inactive-solution composition by at
  most `5.18E-5`.

Post-assemblage phase-amount and element-potential traces place the first
resolved continuous-state separation before, or at, the first discrete phase
split.  This evidence rules out a single universal explanation such as an
incorrect V4.1 assessment, a missing parameter family, or a defective removal
condition.  Depending on the state, a locally accepted correction can either
move the phase-selection coordinates directly or introduce a small continuous
change that is amplified by subsequent nonlinear iterations and a nearly
degenerate add/remove boundary.

### Bounded same-state phase-path gate experiment

A bounded experiment next evaluated every candidate that had already passed
the existing correction-ratio, linear-solve, finiteness, and grouped-update
trust checks.  From the same pre-step state it measured:

- changes in the leading inactive pure or solution phase;
- crossings of the production phase-addition eligibility threshold;
- resolved reversals of the pure-versus-solution ordering;
- crossings of the production active-phase removal amount threshold; and
- the largest scaled force and active-phase target-amount displacement.

These are local candidate diagnostics only.  They do not consult a known final
assemblage, a historical equilibrium result, a database name, or a prescribed
phase identity.  Candidate evaluation snapshots and restores all production
arrays and therefore does not seed the live calculation.

The portable successful set consisted of one public SUBG state and twelve
public SUBQ FeTiVO states.  Across all twelve FeTiVO adaptive runs, the
diagnostic observed **zero** leading-identity, eligibility, ordering, and
removal crossings.  Cu-Fe-C had 130 leading-identity changes in 284 measured
candidates, but zero eligibility, ordering, or removal crossings.  This shows
why leading identity by itself is not an acceptable gate: a numerically leading
inactive phase can change without changing any production decision boundary.

The six diverting FLiBe runs produced the following ungated measurements:

| Database | LiF/BeF2 | Candidates | Identity | Eligibility | Ordering | Removal |
|---|---:|---:|---:|---:|---:|---:|
| V4.1 | 45/55 | 109 | 8 | 6 | 4 | 0 |
| V4.1 | 50/50 | 294 | 50 | 34 | 5 | 0 |
| V4.1 | 55/45 | 30 | 8 | 6 | 3 | 0 |
| V3.1 | 45/55 | 87 | 6 | 5 | 0 | 0 |
| V3.1 | 50/50 | 85 | 8 | 3 | 1 | 0 |
| V3.1 | 55/45 | 88 | 16 | 9 | 4 | 0 |

The separation is real but not sufficient for a local acceptance rule.  Three
test-only policies were evaluated and then removed:

1. an amount-crossing gate rejected no candidate and reproduced the ungated
   reference-state differences;
2. an eligibility/ordering gate rejected between 9 and 45 candidates per run,
   but all six calculations still ended outside the established reference-state
   tolerances; and
3. the combined gate was equivalent to the ranking gate because no immediate
   removal crossing occurred.

The negative result is important.  The V4.1 50/50 gas removal is preceded by a
continuous multi-iteration separation, so the phase amount is still locally
safe at the earlier accepted candidates.  Rejecting only the later local
ranking crossings also changes the path without proving that the replacement
path reaches the desired minimum.  Accordingly, no phase-path acceptance gate
is retained in the solver.  Only the opt-in read-only diagnostic remains.

The next remedy must therefore remain model-independent and must not hard-code
the historical phase identity, database, composition, or iteration.  Before a
multi-iteration rollback is selected, the next bounded experiment should test
whether candidate acceptance can monitor the immediate phase-driving-force and
active-phase amount consequences without rejecting benign near-degenerate
steps.  If no local phase-aware criterion can distinguish the successful and
diverting candidates, checkpoint recovery becomes justified as a genuinely
longer-horizon globalization mechanism rather than a substitute for diagnosis.
Neither outcome implies that the local Hessian or reduced mapper should be
altered.  MQ-4E-B remains blocked until the assessed live-solver gate passes.

The MQ-4E-B **assessment and localization work is complete, but its live-solver
exit gate is blocked**.  Before MQ-4E-B can pass, the point at which accepted
curvature diverts the assemblage search must be addressed.  A phase-identity/
Gibbs-aware multi-iteration recovery mechanism is one candidate, but it is not
yet selected; other phase-path safeguards may be evaluated first.  The three
assessed states must then be rerun against the unchanged historical tolerances.
Do not widen the tolerances, weaken trust thresholds, or count an alternative
converged assemblage as a pass.  MQ-4E-C remains evidence synthesis after this
solver gate is resolved.  `B` remains controlled standalone coverage, and
constrained-KKT execution remains outside this unconstrained assessed gate.

### Compound-component and live matrix-rank diagnostic

A subsequent bounded experiment tested the proposed explanation that elemental
Li-Be-F components leave the live GEM system rank-deficient and that replacing
them by LiF and BeF2 components would restore uniqueness.  The diagnostic used
Thermochimica's existing compound-component input path and captured the exact
live GEM matrix.  Numerical rank was measured with an SVD after scaling each
matrix by its largest absolute entry, using

```text
rank tolerance = matrix dimension * machine epsilon * largest singular value.
```

For both the MSD-TC V4.1 and MSTDB-TC V3.1 fluoride assessments, at all three
45/55, 50/50, and 55/45 compositions:

- the final historical elemental-component matrix had rank `3/4`;
- the first paired pre-correction elemental matrix had rank `4/6`;
- the corresponding fixed-`alpha=1` corrected matrix also had rank `4/6`;
- the historical LiF-BeF2 compound-component matrix had rank `3/3`; and
- the compound-basis condition estimates ranged from `9.20` to `16.45`.

Thus the live elemental GEM system is genuinely rank-deficient, but the MQMQA
correction does **not create** that deficiency: the paired baseline and
corrected matrices have the same numerical rank.  Rank deficiency alone also
does not prove that the returned element potentials are physically corrupt,
because historical Thermochimica converges reproducibly while using the same
deficient component representation.  It instead establishes a concrete
mechanism by which two otherwise acceptable linear solves can choose different
representatives of a non-unique potential space and therefore alter downstream
phase-driving-force coordinates.

The LiF-BeF2 result is diagnostic evidence, not an acceptable production fix.
Several internal MQMQA pair/quadruplet states in the assessed `MSFL` model
cannot individually be represented as nonnegative combinations of LiF and
BeF2.  Thermochimica's compound conversion consequently assigns those entries
zero compound stoichiometry and removes them.  After setup filtering was made
consistent for solution species, MQM pairs, interpolation overrides, and
temporary array capacity, the historical compound-component calculation
converged and its reduced matrix was full rank.  Attempting to continue into
the corrected MQMQA path then encountered missing pair topology in the SUBQ
energy evaluation.  The apparent improvement in rank therefore comes with a
changed and incomplete thermodynamic model, rather than a coordinate-only
reparameterization of the original assessment.

The temporary setup-filtering edits used to complete this bounded experiment
were removed afterward.  They are not retained as production changes because
the compound representation is not a valid remedy and initially disturbed an
existing elemental SUBQ regression.

Accordingly, compound components are rejected as the MQ-4E-B remedy.  The next
bounded experiment should retain the elemental formulation and compare the
current linear solve with a rank-revealing minimum-norm solve on the *same*
paired baseline and corrected matrices.  It must report null-space residuals,
element-potential and update differences, and the resulting inactive-phase
driving-force rankings before any pseudoinverse or regularization is considered
for production.  The experiment must also establish whether the relevant
rank deficiency is structural and shared by historical Thermochimica, rather
than describing it as an MQMQA-Hessian defect.

### Same-matrix minimum-norm replay

That comparison was then implemented as an opt-in, read-only diagnostic.  At
the first live corrected trial, the exact elemental-component baseline system
and its fixed-`alpha=1` correction were each solved twice:

1. with production `DGESV`; and
2. with an SVD pseudoinverse using the same numerical-rank tolerance as the
   preceding rank experiment.

Neither replay replaces the Newton update returned by `GEMNewton`.  The study
records direct right-hand-side residuals, normwise backward errors, solution
norms, matrix-relative difference residuals, element-potential differences,
and inactive-phase rankings.  It was repeated for the V4.1 and V3.1 fluoride
assessments at 45/55, 50/50, and 55/45 LiF/BeF2.

All twelve paired matrices again had numerical rank `4/6`.  Their smallest
retained scaled singular values were `0.112`--`0.245`, while the largest
discarded values were only `1.33e-17`--`1.83e-16`, below tolerances of
`1.67e-15`--`2.70e-15`.  Production LU returned very large raw solution norms,
`8.27e17`--`5.77e18`.  Its normwise backward errors remained near machine
precision.  Eleven of the twelve direct residuals relative to the right-hand
side ranged from `0.441` to `2.51`; the remaining corrected V3.1 45/55 replay
was `6.32e-16`.  This combination is evidence of severe forward sensitivity
and cancellation: `DGESV` is backward stable, yet the raw solution is not a
reliable unique representative of the nearly singular system.

The SVD replay separated the three compositions:

| LiF/BeF2 | SVD result on baseline/corrected systems | Baseline-to-corrected element-potential difference |
|---|---|---:|
| 45/55 | consistent minimum-norm solution; direct residual `3.13e-16`--`6.35e-16` | `<= 2.12e-15` |
| 50/50 | least-squares only; direct residual `0.340`/`0.478` | `<= 2.66e-15` |
| 55/45 | least-squares only; direct residual `0.341`/`0.471` | `<= 2.12e-15` |

The V4.1 and V3.1 matrices gave the same SVD classification and essentially
the same minimum-norm phase rankings.  At 45/55, the minimum-norm baseline and
corrected systems produced the same `MSFL` solution-phase driving force
(`-0.07236`); their leading pure-phase forces were numerically tied near zero.
At 50/50 and 55/45, the baseline and corrected minimum-norm element potentials
were also identical within approximately `3e-15` and gave the same leading
inactive-phase rankings.  In contrast, the LU representatives differed
strongly between baseline and corrected systems and produced driving forces as
large as approximately `1e18`.

This supports the narrower causal mechanism: the elemental GEM system contains
near-null directions, and the MQMQA correction can change which enormous LU
representative is selected even though it does not create the deficiency.
Those representative changes can move inactive-phase rankings.  It does not,
however, justify replacing `DGESV` with the tested pseudoinverse.  At 50/50 and
55/45 the standard rank cutoff discards directions needed to reproduce the
right-hand side, so the bounded SVD solution changes the linear problem into a
least-squares approximation.  The full solution vector also changes by about
`0.285`--`0.305` between the baseline and corrected least-squares systems even
though their element-potential blocks agree, showing that a raw unscaled
minimum norm mixes physically different GEM variable groups.

Therefore the minimum-norm replay is positive diagnostic evidence but a
negative production-fix result.  A production remedy, if pursued, must define
the structural gauge or component constraints explicitly and use physically
meaningful variable scaling or a constrained solve.  It must preserve the
original GEM equations and verify phase-driving-force invariance; simply
discarding small singular values is not acceptable.  No production solve,
phase decision, default, or MQMQA Hessian formula was changed by this study.

## Null-mode, gauge, and exact-scaling diagnostic

The two discarded right singular vectors were next tested as possible physical
gauges. For each paired baseline/corrected matrix, the diagnostic records
`||A v||/||A||`, the discarded left-vector projection of the right-hand side,
the corresponding singular coefficient, the element/solution/pure block
content, and the change in the production leading phase forces produced by a
unit maximum element-potential perturbation along the mode. It then performs
two read-only solves:

1. a null-space selection minimizing the physically grouped element-potential,
   constituent-logarithm, and pure-amount displacements; and
2. a row/column-equilibrated `DGESV` solve of the algebraically unchanged
   square equations.

The 45/55 systems exhibit a genuine gauge-like regime. Their discarded-mode
right-hand-side projections are `3.82e-17`--`1.91e-16`, their phase-force
changes are at most `2.39e-15`, and the physically selected representative has
norm approximately `169` while satisfying the original equations to
`3.13e-16`--`6.39e-16` relative residual. Baseline and corrected rankings are
unchanged. This confirms that bounded, physically selected representatives can
exist when the right-hand side is compatible with the numerical range.

The 50/50 and 55/45 systems are qualitatively different. Their discarded-mode
right-hand-side projections are `1.71e-2`--`4.73e-1`; division by singular
values near `1e-17` requires coefficients of `2.72e16`--`2.18e18`. Their
phase-force changes under the normalized mode perturbation are approximately
`0.1375` and `0.1177`, respectively. The leading phase identities do not
change in this local unit probe, but the forces themselves are not gauge
invariant. Physically grouped null-space selection remains bounded
(`127`--`149`) only by retaining the same least-squares residuals as the
truncated SVD: approximately `0.340`/`0.478` at 50/50 and
`0.341`/`0.471` at 55/45 for baseline/corrected systems.

Exact algebraic equilibration does not remove the obstruction. Depending on
roundoff and database version, `DGESV` either reports a singular pivot or
returns a solution of order `1e18`; its small normwise backward error coexists
with direct right-hand-side residuals as large as `3.95` because the solve is
again cancellation dominated. Scaling therefore improves neither uniqueness
nor forward reliability.

This rejects a blanket gauge-fixing or scaling-only production change. The
structural deficiency is real, but only the compatible 45/55 state admits the
desired bounded exact representative under the captured equations. At the
other compositions, an exact bounded update would require changing or reducing
the assembled equations based on a separately justified thermodynamic
constraint; choosing a different numerical norm cannot create it. No
production linear algebra, phase search, or adaptive policy was changed by
this diagnostic.

### Structural origin of the incompatible FLiBe forcing

The discarded left-null projections were then decomposed by exact GEM equation
identity, and the active phase stoichiometry block

\[
C=A_{1:n_e,\,n_e+1:n_e+n_p}
\]

was analyzed independently.  In all six V4.1/V3.1 captures, the three active
phase columns had numerical rank `2/3`.  The element-inventory forcing lay in
`range(C)` to `1e-17`--`2e-16`, so element balance is not the source of the
incompatibility.  The obstruction instead lies in the phase-equation forcing,
which must belong to `range(C^T)` for all active stationarity equations to be
satisfied simultaneously.

| LiF/BeF2 | Dependent active-phase combination (sign is arbitrary) | Raw phase-energy closure | Interpretation |
|---|---|---:|---|
| 45/55 | `-1.000 MSFL + O(1e-16) Li2BeF4 + O(1e-16) BeF2` | `9.08e-14` | `MSFL` stoichiometry-column norm is `1.38e-15` and its phase RHS is zero; this is an effectively absent phase coordinate |
| 50/50 | `-0.8341 MSFL + 0.4970 Li2BeF4 + 0.2393 BeF2` | `180.9817` | substantive dependent three-phase combination does not satisfy energy closure |
| 55/45 | `-0.8351 MSFL + 0.4956 Li2BeF4 + 0.2386 BeF2` | `180.4419` | substantive dependent three-phase combination does not satisfy energy closure |

The V4.1 and V3.1 assessments gave the same ranks, phase combinations, and
closure magnitudes.  SVD-vector signs differed in some captures, as expected,
without changing the result.  The normalized phase-energy residual changed
between the baseline and corrected systems because the norm of the complete
right-hand side changed; the **raw closure was identical** before and after the
MQMQA correction.  The baseline-to-corrected left and right null-subspace
rotations were zero to approximately `3e-8`, and the direct correction forcing
projected onto the baseline left null space at approximately `1e-16`.

Therefore the MQMQA Hessian correction neither creates the dependent phase
combination nor injects its incompatible forcing.  The historical GEM assembly
already contains an inconsistent singular saddle system at the captured
50/50 and 55/45 transient assemblages.  Production `DGESV` can still return an
enormous cancellation-dominated representative of that system, and a small
curvature perturbation in its well-determined range can select a materially
different representative and phase-search path.  This explains why local
derivative, mapping, and transactional checks can all pass while the global
assemblage trajectory changes.

This diagnosis rejects raw pseudoinverse substitution, gauge selection,
scaling alone, and phase-specific rollback as root-cause remedies.  A future
production remedy must either (a) detect an overcomplete active phase set and
reduce it through model-independent thermodynamic phase-selection logic, or
(b) reformulate the GEM equations on an independently justified component
basis while preserving the complete SUBQ topology.  No such production change
is made in this checkpoint.

### Why the rank-deficient three-phase set survives the active-set controls

A diagnostic trace of every `CheckPhaseChange` trial localized when the
overcomplete V4.1 50/50 assemblage enters the live path.  The three-phase set
is **not** accepted during initialization.  At global iteration zero,
`MSFL + Li2BeF4(s) + BeF2(l)` is tried three times and rejected each time:
the active stoichiometry block has rank `2/3`, and the largest Newton update is
`6.897383e17`, above the then-active `1e14` trial threshold.  At iteration two,
the two-pure-phase set `Li2BeF4(s) + BeF2(l)` is accepted with rank `2/2`.
At iteration seven, `gas_ideal` is added and the resulting three-phase set has
full active-stoichiometry rank `3/3`.  Thus three phases in a Li-Be-F system
are not intrinsically invalid; the defect is the later acceptance of a
linearly dependent active set.

The decisive transition occurs at global iteration 49.  The trial
`MSFL + Li2BeF4(s) + BeF2(l)` again has rank `2/3`, but `CheckPhaseChange`
accepts it because `GEMNewton` returns `INFO=0` and the maximum update,
`9.659123e15`, is below the time-dependent threshold, which has relaxed to
`1e20`.  The trial check contains no active-stoichiometry rank or dependent
phase-energy closure criterion.

Three existing controls then fail to remove the structural redundancy for
distinct reasons:

1. `CheckPureConPhaseAdd` and `CheckSolnPhaseAdd` use the nominal limit
   `nElements - nChargedConstraints`.  They permit direct addition below that
   count and require swapping only at equality.  They do not reduce this limit
   when the active phase stoichiometry spans fewer independent directions.
2. `CorrectPhaseRule` is called only when the number of active phases is
   **greater than** that same nominal limit.  Here the count is `3`, the limit
   is `3`, and the measured rank is only `2`; the count-based correction is
   therefore never triggered.
3. The ordinary pure- and solution-phase removal paths are amount based.  At
   the accepted iteration-49 trial the minimum active phase amount is `1.0`,
   while the removal tolerance is `1e-11`, so no member of the dependent set
   is considered small enough to remove.

The causal chain is therefore: a later global phase-change trial creates a
dependent three-phase candidate; a relaxed update-size gate accepts it; the
nominal phase-rule count does not recognize it as overcomplete; and the
amount-based removal logic retains it.  This is a pre-existing, model-
independent active-set deficiency exposed by the assessed FLiBe trajectory,
not an error in the MQMQA Hessian formulas.  The curvature correction can
change which cancellation-dominated solution and later phase path is selected,
but it does not create the rank-two stoichiometry or the raw closure defect.
The reported update magnitudes are from the accepted Linux/Docker verification
environment.  The macOS Accelerate build returned different order-`1e15`--
`1e18` representatives while preserving every rank, threshold, and pass/fail
classification, as expected for the diagnosed inconsistent singular system.

The next bounded investigation should test a model-independent structural gate
at the phase-addition/swap acceptance boundary.  A rank-deficient candidate
must not be rejected merely because its columns are dependent: a redundant
coordinate may be harmless when its dependent phase-energy closure is also
compatible, as in the effectively absent 45/55 `MSFL` coordinate.  The gate
must distinguish that case from the substantive 50/50 and 55/45 closure
defects, and any rejection or phase reduction must preserve a feasible lower-
Gibbs assemblage.  This checkpoint records the cause only; it does not alter
production phase selection, removal, thresholds, or linear algebra.

### Bounded structural-gate experiment: rejected

The proposed gate was implemented temporarily behind a private, default-off
diagnostic control.  At every `CheckPhaseChange` trial it formed the active
element-by-phase stoichiometry block `C`, computed its numerical right null
space, and rejected an otherwise historically acceptable candidate when

\[
\frac{\lVert Y^T g_{\rm phase}\rVert_2}
     {\max(\lVert g_{\rm phase}\rVert_2,1)} > 10^{-10},
\]

where the columns of `Y` span the dependent active-phase combinations.  The
complete private matrix covered LiF mole fractions 0.45, 0.50, and 0.55 in
both MSD-TC V4.1 and MSTDB-TC V3.1.  Historical, fixed-alpha-one, adaptive,
forced-zero, rejected-candidate, and one-step-ablation modes remained
separately classified.

The experiment failed its intended discrimination test.  In all six assessed
states, the first structural rejection occurred at the same decisive early
phase-set check (global iteration 47--55, depending on database and
composition), with rank `2/3`, raw closure `0.4157`, and normalized closure
`1.226e-3`.  The gate then rejected hundreds to thousands of dependent
candidates per corrected run.  Nevertheless, every fixed or adaptive
positive-curvature result remained different from the historical reference;
the V3.1 50/50 fixed-alpha calculation failed rather than recovering the
reference path.  The 45/55 case was especially decisive: although its later
captured corrected system admits a compatible bounded minimum-norm
representative, the trial-time gate still rejected intermediate dependent
sets and did not preserve the reference trajectory.

The reason is mathematical rather than a tolerance-tuning accident.  The
condition `Y^T g_phase = 0` is required when the dependent phase stationarity
equations are simultaneously satisfied at the state being characterized.  A
candidate passed to `CheckPhaseChange` is generally an off-equilibrium Newton
trial, so its current phase-energy residual need not already satisfy that
closure.  Enforcing near-zero closure at every trial confuses an equilibrium
compatibility condition with a nonlinear-iteration acceptance condition.  A
looser fixed threshold would merely move this arbitrary boundary and would
not establish that the retained candidate is feasible, lower in Gibbs energy,
or on the correct basin.

Accordingly, the gate, its counters, and its private command-line control were
removed after the experiment.  The existing read-only rank and null-space
diagnostics were retained.  No public API, production default, phase-search
decision, or MQMQA Hessian formula was changed.  This negative result rules
out instantaneous dependent-energy closure as a standalone phase-set gate;
any future rank-aware phase reduction must evaluate a thermodynamically
defined reduced candidate or nonlinear progress rather than demand
equilibrium closure before the candidate has converged.

### Bounded reduced-phase-set convergence experiment

The next diagnostic replaced the invalid off-equilibrium closure question with
a thermodynamic comparison of independently converged reduced phase sets.  At
the first adaptive FLiBe phase-addition trial whose active stoichiometry block
had rank `2/3`, the diagnostic recorded:

- the incumbent assemblage immediately before the trial;
- the rank-dependent three-phase trial assemblage; and
- the global iteration and measured rank.

It then restored the same database, temperature, pressure, and elemental
inventory for each independent calculation and ran ordinary Thermochimica
restricted equilibrium for (a) the incumbent set and (b) each set formed by
omitting one member of the dependent trial.  The historical final assemblage
was not used to select a result, no phase name was hard-coded into the
selection rule, and the MQMQA curvature correction was disabled during these
comparison solves.  Feasibility required ordinary convergence, finite Gibbs
energy, mass balance within the production tolerance, nonnegative phase
amounts within tolerance, and a final active assemblage contained in the
candidate set.

Thermochimica's setup requires a pure-element reference species for every
system element.  The diagnostic therefore constructs the smallest greedy set
of reference-carrier phases from parsed species stoichiometry.  Both FLiBe
databases selected `gas_ideal`.  A carrier may participate in setup, but a
candidate is rejected if a carrier that was not requested becomes active in
the converged result.  In the reported matrix the carrier remained inactive
for every reduced candidate.

| Database | LiF/BeF2 | Captured iteration | Rank | Lower-Gibbs distinct assemblage | Pure-set minus liquid Gibbs (J) |
|---|---:|---:|---:|---|---:|
| MSD-TC V4.1 | 45/55 | 54 | 2/3 | `MSFL` | 4747.8 |
| MSD-TC V4.1 | 50/50 | 46 | 2/3 | `MSFL` | 4587.8 |
| MSD-TC V4.1 | 55/45 | 46 | 2/3 | `MSFL` | 4127.0 |
| MSTDB-TC V3.1 | 45/55 | 54 | 2/3 | `MSFL` | 4747.8 |
| MSTDB-TC V3.1 | 50/50 | 46 | 2/3 | `MSFL` | 4587.8 |
| MSTDB-TC V3.1 | 55/45 | 46 | 2/3 | `MSFL` | 4127.0 |

All `24/24` restricted candidate calculations converged and satisfied the
feasibility checks.  In every state, omitting either one of the redundant pure
phases allowed the calculation to converge to the same distinct `MSFL`
assemblage.  Omitting `MSFL`, or retaining the incumbent set, converged to
`Li2BeF4(s) + BeF2(l)` at higher Gibbs energy.  The two routes that converged
to `MSFL` differed by at most approximately `0.051 J`, or about `6e-8`
relative to the system Gibbs magnitude.  That small route-dependent spread is
reported as numerical solver evidence, not as a thermodynamically distinct
winner between the two omission labels.

This result supplies the missing positive evidence from the rejected
instantaneous gate experiment: for all six captured decisions, a feasible
independent reduced set exists, preserves the liquid `MSFL` phase, and has
lower Gibbs energy than the pure-phase alternative.  It supports a future
model-independent strategy that detects an overcomplete active trial,
converges feasible reduced candidates, and selects among **distinct converged
assemblages** by Gibbs energy.  It does not yet establish production behavior:
candidate enumeration cost, state restoration, duplicate-assemblage handling,
failure recovery, and the exact trigger boundary still require design and live
solver verification.  Accordingly, no phase-selection logic or production
default was changed by this checkpoint.

### Rank frequency, candidate cost, and restoration hardening

The one-phase-omitted construction above is complete only for the observed
`rank 2/3` trials.  In general, if an active trial contains `k` phase columns
with numerical rank `r`, its nullity is `d=k-r`, at least `d` columns must be
removed to form an `r`-column basis, and a brute-force search can expose up to

\[
\binom{k}{r}=\binom{k}{d}
\]

candidate subsets.  Not every such subset is guaranteed to be independent.
Therefore, the diagnostic now records the complete number of phase-change
checks, rank-deficient checks, rank-deficient checks that passed the existing
local `CheckPhaseChange` test, maximum nullity, and maximum combinatorial basis
count.  A bounded detailed history retains at most 512 individual records,
while separate aggregate counters continue to the end of the calculation and
report any discarded detailed records.  Here, *locally passing* means only
that `CheckPhaseChange` returned `lPhasePass=.TRUE.`; it does not by itself
mean that the outer active-set algorithm permanently accepted the phase set.

The 13-state portable MQ-4E-A matrix was rerun in fixed-alpha-one and adaptive
modes with rank capture enabled.  Across all 26 corrected calculations, every
recorded phase-change trial was full column rank: there were zero
rank-deficient checks, zero locally passing rank-deficient checks, zero
dropped detailed records, and maximum nullity zero.  The public CuFeC and
FeTiVO cases therefore do not show that dependent candidate phase sets are a
generic consequence of enabling MQMQA curvature.

The two private FLiBe database versions show a sharply different pattern:

| Database | LiF/BeF2 | All checks | Rank deficient | Locally passing deficient | Maximum nullity | Maximum `C(k,r)` | Forward/reverse CPU time (s) | Restoration mismatches |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| MSD-TC V4.1 | 45/55 | 1179 | 1098 | 8 | 1 | 3 | 0.341/0.337 | 0 |
| MSD-TC V4.1 | 50/50 | 919 | 795 | 17 | 1 | 3 | 0.245/0.238 | 0 |
| MSD-TC V4.1 | 55/45 | 380 | 354 | 6 | 1 | 3 | 0.238/0.229 | 0 |
| MSTDB-TC V3.1 | 45/55 | 889 | 817 | 11 | 1 | 3 | 0.224/0.221 | 0 |
| MSTDB-TC V3.1 | 50/50 | 297 | 242 | 6 | 1 | 3 | 0.215/0.208 | 0 |
| MSTDB-TC V3.1 | 55/45 | 138 | 78 | 4 | 1 | 3 | 0.147/0.151 | 0 |

Thus `3384/3802` observed FLiBe phase-change trials were rank deficient, but
only `52` of those deficient trials passed the existing local phase-change
test.  Every observed deficiency was `rank 2/3`, so the three one-phase-
omitted subsets are the complete `C(3,2)=3` reduced-basis search for this
matrix.  No current evidence establishes that future MQMQA systems must also
have nullity one; higher-nullity systems could create combinatorial growth.

For each of the six captured decisions, the diagnostic solved four candidates
in forward order--the incumbent plus all three one-phase-omitted subsets--and
then repeated the same four candidates in reverse order.  This produced
`24+24` restricted equilibrium calculations.  All forward and reverse solves
remained feasible, the same candidate labels returned the same final
assemblages and Gibbs energies within the stated numerical comparison, and
all six cases reported zero restoration mismatches.  The summed CPU times
were approximately `1.41 s` forward and `1.38 s` in reverse order on the test container.  These
times characterize the small external diagnostic only; they are not a
production performance estimate.

This establishes order-independent restoration for the conservative oracle
used here: each candidate calls `ResetThermoAll`, reparses the same database,
and reapplies the same temperature, pressure, elemental inventory, and phase
restriction.  It does **not** establish that an in-place nested solve can save
and restore every mutable GEM, phase-history, line-search, and adaptive-trust
field.  A production implementation must either retain this expensive fresh-
solve isolation or separately prove an explicit solver-state transaction.

The measured trigger frequency also rules out blindly launching a reduced-set
search at every rank-deficient FLiBe check.  The next production design must
combine rank exposure with a later, thermodynamically meaningful trigger,
screen dependent subsets before nonlinear solution, group duplicate converged
assemblages, bound candidate count and failure handling, and demonstrate that
candidate evaluation cannot recursively invoke itself.  None of those
production mechanisms is implemented by this bounded checkpoint.

### Default-off rank-rejection prototype: negative result

MQ-4E-B.1 tested the cheapest possible use of the rank information before
implementing reduced-candidate solves.  An opt-in prototype was inserted after
the ordinary `CheckPhaseChange` Newton trial and its existing pass/fail tests.
It normalized each active phase's structural stoichiometry column, measured
the column rank by SVD, and changed an otherwise passing result to
`lPhasePass=.FALSE.` whenever the rank was smaller than the number of active
phases.  Rank-analysis failure preserved the ordinary result.  The prototype
also counted total rejections, repeated rejection of the same assemblage, and
the longest consecutive rejection sequence.  It did not use phase names,
database identities, a historical result, phase amounts as column scales, or
the MQMQA correction mode in its decision.

The public 13-state MQ-4E-A matrix remained a clean control.  Its 26 fixed and
adaptive calculations passed, and the guard observed zero dependent trials
and made zero interventions.  The six FLiBe states, however, showed that
structural dependence is not by itself a valid reason to reject a temporary
phase-search trial:

| Database | LiF/BeF2 | Guarded historical rejects | Repeated | Maximum consecutive | Guarded/reference iterations | Fixed result | Adaptive result |
|---|---:|---:|---:|---:|---:|---|---|
| MSD-TC V4.1 | 45/55 | 447 | 444 | 381 | 2002/1267 | iteration-limit failure | reference-state difference |
| MSD-TC V4.1 | 50/50 | 946 | 939 | 377 | 4841/1586 | reference-state difference | reference-state difference |
| MSD-TC V4.1 | 55/45 | 41 | 40 | 41 | 318/3511 | reference-state difference | reference-state difference |
| MSTDB-TC V3.1 | 45/55 | 175 | 160 | 45 | 1333/1526 | reference-state difference | reference-state difference |
| MSTDB-TC V3.1 | 50/50 | 731 | 588 | 398 | 5113/902 | reference-state difference | reference-state difference |
| MSTDB-TC V3.1 | 55/45 | 171 | 157 | 50 | 1188/2700 | reference-state difference | reference-state difference |

The guarded historical calculation itself no longer reproduced the ordinary
historical result.  Of the twelve corrected fixed/adaptive calculations,
eleven ended in a different state and the remaining fixed calculation reached
its iteration limit.  A smaller iteration count in an individual guarded row is
not an improvement because that calculation failed the reference-state gate.
The hundreds of repeated rejections show that the outer search can legitimately
revisit temporary dependent assemblages while attempting to find a feasible
independent equilibrium set.  Rank deficiency describes the geometry of the
trial; it does not determine which dependent phase should be removed or which
reduced set has the lowest converged Gibbs energy.

Accordingly, the rank-rejection implementation and its command-line control
were removed after the experiment.  No production behavior, public API, or
default was retained.  The result strengthens the next design boundary:
rank deficiency may trigger a bounded reduced-candidate comparison, but it
cannot be used as an unconditional phase-rejection rule.  The positive reduced-
set evidence above remains the thermodynamic basis for MQ-4E-B's next
prototype.

### MQ-4E-B.2: default-off RRQR-screened reduced-candidate comparison

The next prototype retained rank information as a **screen and ordering tool**
rather than turning it into an accept/reject rule.  It activates only in the
private diagnostic when the ordinary `CheckPhaseChange` result already passes
and the active phase columns are dependent.  The production default and phase-
selection result remain unchanged.

For the captured active structural stoichiometry matrix, the prototype:

1. normalizes each phase column so pivot order is not determined merely by its
   current phase amount;
2. applies rank-revealing QR with column pivoting (`DGEQP3`);
3. retains the pivot order as a numerical independence ordering, not as a
   thermodynamic phase choice;
4. constructs the complete three one-phase-omitted subsets for the observed
   `rank 2/3` case;
5. rejects a subset before nonlinear solution unless its SVD rank is two;
6. converges each retained subset independently from the same database,
   temperature, pressure, and elemental inventory; and
7. groups duplicate converged assemblages and compares the distinct results by
   Gibbs energy.

The search is explicitly capped at eight candidates.  All currently observed
FLiBe trials require only three, but the cap is a diagnostic safety boundary,
not evidence that future databases cannot exhibit larger nullity or candidate
counts.

| Database | LiF/BeF2 | Captured iteration | Rank | RRQR primary basis contains `MSFL` | Lowest-Gibbs distinct result | Forward/reverse CPU time (s) |
|---|---:|---:|---:|:---:|---|---:|
| MSD-TC V4.1 | 45/55 | 54 | 2/3 | no | `MSFL` | 0.232/0.224 |
| MSD-TC V4.1 | 50/50 | 46 | 2/3 | yes | `MSFL` | 0.162/0.163 |
| MSD-TC V4.1 | 55/45 | 46 | 2/3 | no | `MSFL` | 0.166/0.171 |
| MSTDB-TC V3.1 | 45/55 | 54 | 2/3 | yes | `MSFL` | 0.175/0.178 |
| MSTDB-TC V3.1 | 50/50 | 46 | 2/3 | yes | `MSFL` | 0.173/0.174 |
| MSTDB-TC V3.1 | 55/45 | 46 | 2/3 | yes | `MSFL` | 0.120/0.113 |

All six runs passed the candidate-cap, nullity, forward-feasibility,
reverse-feasibility, distinct-result, and restoration gates.  Every candidate
was full rank after its indicated omission, every forward and reverse solve
was feasible, and all six reverse-order comparisons reported zero restoration
mismatches.  Each state produced two distinct equilibria: liquid `MSFL` and
the two-pure-phase alternative.  `MSFL` was lower in Gibbs energy in every
case, consistent with the earlier unscreened reduced-set experiment.

The pivot evidence is deliberately non-vacuous.  In the V4.1 45/55 and 55/45
states, RRQR selected the two pure phases as its primary numerical basis; that
basis converged approximately 4.75 kJ and 4.13 kJ above `MSFL`, respectively.
Thus RRQR is useful for exposing independent bases cheaply, but **cannot** be
used alone to decide which phases are thermodynamically stable.  The converged
Gibbs comparison is the step that correctly overrules the numerical pivot
choice.

A bounded attempt was also made to hand the selected restricted equilibrium
to an unrestricted curvature-enabled solve through Thermochimica's existing
reinitialization API.  That route is not structurally valid: the restricted
parse contained 27 solution variables whereas restoring full phase
availability produced 38, and the bounds-checked replay detected the `38/27`
dimension mismatch in `SwapSolnPhase`.  The replay code was removed.  This is
not a failure of the RRQR/Gibbs candidate oracle; it confirms that a production
handoff needs either an explicit full-space state mapping or a proven
transaction that preserves one common parsed phase space.  Ordinary reinit
data cannot be copied between differently restricted phase spaces.

MQ-4E-B.2 therefore establishes a successful, default-off **selection
diagnostic**, not production solver integration.  Before this mechanism can
become live logic, the implementation must still define a later trigger that
avoids launching it at every dependent FLiBe trial, evaluate candidates
without recursive invocation, restore one full-space solver state safely,
and demonstrate corrected fixed/adaptive convergence after the selected basis
is installed.  No production phase behavior, public control, or default has
been added at this checkpoint.

### Null-direction thermodynamic basis pivot: positive bounded diagnostic

The reduced-candidate results suggested a narrower alternative to enumerating
and converging every one-phase-omitted set.  For the first locally passing,
rank-dependent active set in each FLiBe run, a private default-off diagnostic
asks whether the dependent phase can be selected directly from the
thermodynamics of that same state.  It does not alter live phase selection.

The mixed GEM columns used by the existing rank capture are not appropriate
for this question because each solution-phase column includes its current
phase amount.  The diagnostic therefore also captures the physical
amount-basis stoichiometry

\[
  C_{\mathrm{amt}} = [c_1\;c_2\;\ldots\;c_k]
\]

in the natural active-phase order.  When the active set has nullity one,
`DGESVD` supplies a right-null vector \(y\) satisfying

\[
  C_{\mathrm{amt}}y \simeq 0.
\]

Consequently, \(n(t)=n+t y\) preserves the elemental inventory to first
order.  The nonnegativity constraints on the phase amounts define a feasible
interval in \(t\); each finite endpoint makes at least one phase amount zero.
With the current solution-phase compositions frozen, the diagnostic forms the
phase-Gibbs vector from `dGibbsSolnPhase/dMolesPhase` for solution phases and
`dStdGibbsEnergy` for pure phases.  It selects the endpoint for which

\[
  \Delta G_{\mathrm{frozen}} = t\,y^T g < 0.
\]

This decision uses only the current active state.  The independently
converged reduced candidates are evaluated afterward as an oracle; their
answers are not supplied to the pivot rule.

| Database | LiF fraction | Active rank | Relative null residual | Predicted \(\Delta G_{\mathrm{frozen}}\) | Predicted leaving phase | Reached lowest-Gibbs branch |
|---|---:|---:|---:|---:|---|---|
| MSD-TC V4.1 | 0.45 | 2/3 | `8.15e-17` | `-2.22e2` | BeF2 liquid | yes |
| MSD-TC V4.1 | 0.50 | 2/3 | `2.02e-17` | `-1.74e2` | BeF2 liquid | yes |
| MSD-TC V4.1 | 0.55 | 2/3 | `2.52e-17` | `-1.24e2` | BeF2 liquid | yes |
| MSTDB-TC V3.1 | 0.45 | 2/3 | `2.86e-17` | `-2.22e2` | BeF2 liquid | yes |
| MSTDB-TC V3.1 | 0.50 | 2/3 | `3.23e-17` | `-1.74e2` | BeF2 liquid | yes |
| MSTDB-TC V3.1 | 0.55 | 2/3 | `3.53e-17` | `-1.24e2` | BeF2 liquid | yes |

The sign and scale of \(y\), and therefore the separate signs and magnitudes
of \(t\) and \(y^Tg\), depend on the arbitrary SVD orientation and
normalization.  Their product is the meaningful quantity.  In all six runs,
the feasible downhill endpoint removed the BeF2 liquid phase and predicted
the same MSFL branch later identified as the lowest-Gibbs independently
converged reduced candidate.

This is stronger evidence than choosing the RRQR-dependent column alone:
RRQR identifies numerical independence, whereas the null-direction pivot
uses the local thermodynamic slope to choose between feasible independent
bases.  It also addresses the cause of the singular trial more directly than
rejecting every dependent set without selecting an alternative.

The claim remains bounded.  The evidence covers only the first captured
locally passing nullity-one decision in six binary FLiBe calculations, uses a
frozen-composition phase-Gibbs slope, and relies on reduced solves only as a
post-decision oracle.  A live prototype must install the selected endpoint
inside one unchanged parsed phase space, update solution-species and phase
amounts consistently, preserve existing behavior for full-rank and
higher-nullity sets, and demonstrate both fixed and adaptive convergence plus
clean regressions.  No production default, public control, or solver behavior
was changed by this diagnostic.

### Live null-direction pure-phase replacement: negative bounded prototype

The positive frozen-state result above was next exercised inside the existing
`AddSolnPhase`/`CheckPhaseChange` transaction.  The default-off prototype acted
only when one newly added solution phase created a nullity-one active set and
the downhill feasible endpoint selected an incumbent pure phase to leave.  It
shifted physical phase amounts along the right-null vector, updated the
solution-species amounts at fixed composition, and used the existing `INFO`
convention to remove the selected pure phase and retest the reduced set.  The
historical equilibrium was used only as a final-state oracle; it was not used
by the live selection rule.

The local calculations were numerically sound.  Every applied pivot preserved
the elemental inventory to approximately `4e-16` or better, and no nonfinite
amount or failed linear solve was observed.  The complete six-state result was
nevertheless negative:

| Database | LiF | Fixed pivot applied/attempted | Fixed final assemblage | Adaptive pivot applied/attempted | Adaptive final assemblage |
|---|---:|---:|---|---:|---|
| MSD-TC V4.1 | 0.45 | 22/31 | `MSFL` | 7/14 | `MSFL` |
| MSD-TC V4.1 | 0.50 | 11/22 | two pure phases | 18/24 | `MSFL` |
| MSD-TC V4.1 | 0.55 | 13/23 | gas | 0/6 | two pure phases |
| MSTDB-TC V3.1 | 0.45 | 17/24 | two pure phases | 2/10 | `MSFL` |
| MSTDB-TC V3.1 | 0.50 | 15/21 | two pure phases plus gas | 4/9 | two pure phases |
| MSTDB-TC V3.1 | 0.55 | 12/24 | two pure phases plus gas | 0/4 | two pure phases |

Fixed alpha one recovered the historical `MSFL` assemblage in only one of six
states; adaptive curvature recovered that named assemblage in three of six.
None of the twelve corrected runs passed the established strict final Gibbs,
composition, species-mole, phase-amount, and assemblage gates.  The failures
were therefore not hidden by accepting name-only agreement.

The decisive limitation was structural rather than a bad null vector.  Several
later dependent trials selected a solution phase, rather than an incumbent
pure phase, as the downhill endpoint.  Other trials had a numerically flat
frozen Gibbs slope, and two adaptive 55/45 runs never encountered a supported
pure-leaving transaction.  Thus the first captured `new solution / old pure`
event was representative enough to motivate the experiment but not sufficient
to describe the live phase-search topology across a complete calculation.

The live code was removed after this negative result.  The experiment rules out
promoting a specialized one-pure-phase replacement as the MQ-4E-B solution.  A
subsequent design must either support solution-phase endpoints and the other
observed transaction classes through a model-independent mechanism, or explain
why those events can safely remain on the historical path.  It must not infer
success merely from a conservative local pivot or from final phase names.

### Historical phase search followed by curvature refinement: negative tradeoff experiment

A bounded private experiment tested whether MQMQA curvature could be removed
from the global phase search and enabled only after the historical solver had
converged.  This is intentionally different from reduced-candidate
enumeration: no phase subsets were generated and no phase choice was imposed.
Each run used one unchanged full parsed phase space, saved the converged
historical state through Thermochimica's existing reinitialization API, and
then repeated the same state with either no curvature, fixed alpha one, or
adaptive curvature.  Direct-from-scratch fixed and adaptive runs were retained
as comparisons.

The alpha-zero restart is an essential control.  It measures the cost and
state change caused by reinitialization itself rather than attributing every
difference to the Hessian.  The table reports total search-plus-refinement
iterations relative to one historical solve.  CPU ratios were also measured,
but the individual calculations are short and single-run timing is
cache-sensitive; iteration ratios are the more reproducible comparison.

| Database | LiF | alpha-zero restart total/base | fixed total/base | fixed result | adaptive total/base | adaptive positive corrections | adaptive strict agreement |
|---|---:|---:|---:|---|---:|---:|:---:|
| MSD-TC V4.1 | 0.45 | 1.213 | 5.736 | iteration limit | 2.684 | 265 | no |
| MSD-TC V4.1 | 0.50 | 1.359 | 4.784 | iteration limit | 1.566 | 131 | no |
| MSD-TC V4.1 | 0.55 | 1.001 | 1.805 | different assemblage | 1.001 | 0 | yes |
| MSTDB-TC V3.1 | 0.45 | 1.230 | 4.933 | iteration limit | 2.170 | 171 | no |
| MSTDB-TC V3.1 | 0.50 | 1.001 | 6.296 | different assemblage | 1.001 | 0 | yes |
| MSTDB-TC V3.1 | 0.55 | 1.095 | 3.223 | iteration limit | 2.207 | 53 | no |

All alpha-zero restarts and all adaptive refinements converged to an `MSFL`
assemblage.  However, the alpha-zero restart reproduced the strict historical
Gibbs/composition/species/phase state in only two of six cases.  In the other
four, the existing reinitialization path itself returned a different internal
MSFL representation despite retaining the same named assemblage.  It added
9.5% to 35.9% more iterations in three of those cases.  The current reinit API
is therefore not a transparent boundary for this proposed two-stage algorithm.

Fixed alpha one preserved the strict historical result in zero of six cases:
four refinements reached the 6001-iteration limit and two converged to
different assemblages.  Adaptive refinement converged in all six, but the only
two strict-agreement cases selected alpha zero throughout.  In every state
where adaptive refinement accepted positive curvature, strict state agreement
failed and total iterations increased by 56.6% to 168.4%.

This experiment rejects the tested **reinitialize-after-convergence** design.
It does not prove that all separation of phase search and curvature is
impossible.  An in-loop transition within one continuous solver state could
avoid the reinitialization drift, but activation after complete convergence
may be vacuous, whereas activation before convergence reintroduces the phase-
path safety question.  That tradeoff must be addressed explicitly before any
such mechanism is promoted to production logic.  No production switch,
rollback, phase freezing, or altered default was added by this experiment.

### Same-state fixed-point preservation and GEM gauge diagnostic

The preceding global-path experiments could not distinguish an incorrect
local curvature integration from a correct local response interacting badly
with Thermochimica's phase search.  A default-off diagnostic therefore captures
one historical FLiBe state at the instant `CheckConvergence` accepts it, before
`PostProcessThermo` changes the internal representation.  Without reparsing,
reinitializing, running a line search, or invoking phase-change logic, it builds
the baseline and fixed-alpha-one GEM systems at that exact live state.  All
mutable solver arrays and controls are restored afterward; the measured
restoration error was exactly zero in all six cases.

Four linear targets are replayed from the captured pair:

1. the historical system, (A_0x_0=B_0);
2. the complete correction,
   ((A_0+\Delta A)x_{AB}=B_0+\Delta B);
3. the matrix-only hybrid,
   ((A_0+\Delta A)x_A=B_0); and
4. the right-hand-side-only hybrid,
   (A_0x_B=B_0+\Delta B).

This is a diagnostic decomposition, not four proposed production algorithms.
The solved targets are converted to the same element-potential,
solution-logarithm, and pure-phase-amount displacement groups used by the trust
layer.  SVD rank estimates and normwise backward errors are reported alongside
the ordinary `DGESV` replays.

| Database | LiF | Historical/corrected rank | Corrected-system backward error at (x_0) | Historical element-potential displacement | Full-correction element-potential displacement | Full-correction solution-log displacement |
|---|---:|---:|---:|---:|---:|---:|
| MSD-TC V4.1 | 0.45 | 3/4, 3/4 | `5.20e-17` | `1.72e18` | `2.71e2` | `2.96e-2` |
| MSD-TC V4.1 | 0.50 | 3/4, 3/4 | `2.03e-17` | `1.55e18` | `3.72e3` | `6.85e-3` |
| MSD-TC V4.1 | 0.55 | 3/4, 3/4 | `9.41e-10` | `2.45e2` | `1.01e19` | `1.38e3` |
| MSTDB-TC V3.1 | 0.45 | 3/4, 3/4 | `3.71e-6` | `1.74e2` | `1.31e2` | `1.13e-1` |
| MSTDB-TC V3.1 | 0.50 | 3/4, 3/4 | `2.50e-10` | `4.94e2` | `1.85e2` | `1.05e-5` |
| MSTDB-TC V3.1 | 0.55 | 3/4, 3/4 | `1.91e-7` | `7.14e2` | `2.03e2` | `6.76e-3` |

Both the baseline and corrected systems are numerically rank three in four
unknowns at every captured state.  Their scaled smallest singular values lie
between approximately `1e-17` and `1e-16`, below the case-local rank tolerances
of approximately `1.5e-15`.  Thus the analytical MQMQA correction neither
creates the rank deficiency nor restores a unique four-variable solve.

The ordinary LU solution is consequently a gauge-dependent representative of
an underdetermined system.  A tiny normwise backward error proves that a
returned vector solves the captured equations; it does **not** prove that the
representative is unique, bounded, or safe for the later inactive-phase
driving-force calculation.  This explains the apparently contradictory
behavior across neighboring V4.1 states: at LiF fractions 0.45 and 0.50 the
complete correction replaces cancellation-dominated element-potential targets
of order `1e18` with bounded targets, while at 0.55 the corrected LU target
instead grows to order `1e19` and contaminates the physical solution-log
group.  All of these linear solves have small backward errors.

The hybrid replays also rule out a simple claim that `deltaB` alone moves a
well-defined fixed point.  In the two V4.1 states where the full correction is
bounded, the `deltaA`-only system is also bounded while the `deltaB`-only
system retains the large baseline scale.  At 0.55, however, no ungauged hybrid
provides a generally safe interpretation.  The matrix correction changes which
representative LU selects, but the response condensation is still embedded in
a singular global coordinate system.

This evidence narrows the missing production mechanism.  The next candidate
must define an explicit, physically scaled gauge or solve in an independent
reduced coordinate space, then prove that phase-driving forces are invariant
to the eliminated null coordinate.  Replacing `DGESV` with an unconstrained
minimum-norm routine is not sufficient by itself: a minimum depends on variable
scaling, and the global phase search must use potentials consistent with the
chosen gauge.  Likewise, phase-specific rollback would treat the downstream
symptom without defining the missing global coordinate convention.  No
production solver behavior or default is changed by this diagnostic.

### Rank-deficient recovery timeline: historical behavior identified

The fixed-point rank evidence established that both the historical and
corrected GEM systems can be rank deficient, but it did not show how the live
historical calculation nevertheless leaves temporary dependent phase sets.
A default-off recovery trace therefore correlated every passing dependent
`CheckPhaseChange` trial with four stages of the same live calculation:

1. the assemblage and functional norm before the main Newton solve;
2. the returned Newton update and selected MQMQA alpha;
3. the final line-search step and functional norm; and
4. the post-`CheckPhaseAssemblage` phase set, phase-history index, and reversion
   state.

A trial was classified as *installed* only when its recorded phase identities
matched the post-assemblage state at that same global iteration.  The trace
then followed that exact phase set until the first later iteration with a
different assemblage.  A changed reversion counter distinguished
`RevertSystem` from ordinary phase addition, removal, or swapping.  Capture is
inactive by default and does not alter any solver decision.

The three MSD-TC V4.1 FLiBe states gave:

| LiF fraction | Historical installed / exited / retained / reverted | Adaptive installed / exited / retained / reverted | Historical iterations | Adaptive iterations | Adaptive final result |
|---:|---:|---:|---:|---:|---|
| 0.45 | 3 / 3 / 0 / 0 | 8 / 8 / 0 / 0 | 1267 | 2884 | same named `MSFL` phase, different internal state |
| 0.50 | 8 / 8 / 0 / 0 | 17 / 17 / 0 / 0 | 1586 | 4720 | `gas_ideal`, higher Gibbs energy |
| 0.55 | 12 / 12 / 0 / 0 | 6 / 6 / 0 / 0 | 3511 | 1002 | two condensed phases, higher Gibbs energy |

Historical Thermochimica did **not** invoke a special rank-aware solve or
`RevertSystem` for any of the 23 installed dependent sets.  It allowed the
ordinary LU update, applied the existing line search, and subsequently left
every dependent set through normal active-set operations.  Some historical
dependent iterations were numerically severe: examples included maximum
updates of approximately `8.2e38` and `1.6e18`, paired with line-search steps
of approximately `1.2e-39` and `2.0`.  Thus historical success does not mean
that rank deficiency was removed or that LU returned a unique physical
representative.  The existing damping and phase-management sequence happened
to continue to a lower-Gibbs final state.

The adaptive calculation used the same recovery machinery.  Almost every
dependent set was installed on an iteration that selected alpha zero; one
0.55-LiF installation retained alpha `0.1`.  Nevertheless, positive curvature
accepted on earlier iterations changed the continuous state delivered to later
phase decisions.  This changed both the identities and frequency of dependent
sets even where the installation iteration itself used the historical GEM
matrix.  For example, at 0.45 LiF the historical path installed three
dependent sets while the adaptive path installed eight; at 0.50 the counts
were eight and seventeen.

This result rules out the hypothesis that a hidden historical rank-recovery
algorithm merely needs to be copied into the corrected path.  There is no such
special mechanism in the observed calculations.  It also shows why rejecting
only the dependent trial is downstream of the first cause: the curvature has
already changed the state and phase-driving-force history before that trial is
proposed.  The next bounded design question is therefore whether adaptive
readiness can recognize *phase-search instability before accepting positive
curvature* while still allowing the exact Hessian to remain active within a
settled assemblage.  Any proposed gate must be based on live, model-independent
phase-search observables and must be tested against both successful FeTiVO
activation and these FLiBe paths.  No such gate was added in this diagnostic.

The existing same-state candidate diagnostic was then used to test the most
local version of that proposal.  At the first observed assemblage divergence,
the alpha-zero and accepted corrected targets were evaluated from the same
pre-step GEM system using the production pure- and solution-phase driving-force
calculations.  The result was:

| LiF fraction | Candidate iteration | Accepted alpha | Alpha-zero versus corrected leading phases | Immediate result |
|---:|---:|---:|---|---|
| 0.45 | 366 | 1 | same pure and solution identities | pure-force change `6.7e-15`; solution composition unchanged |
| 0.50 | 145 | 0.01 | same pure and solution identities | both leading forces exactly zero at report precision; solution composition unchanged |
| 0.55 | 240 | 1 | same pure and solution identities | pure-force change `4.4e-16`; solution composition unchanged |

Thus the first divergent phase decision was **not** preceded by an immediate
identity, eligibility, ordering, or active-amount boundary crossing that could
uniquely reject the accepted candidate.  The candidate can be locally benign
according to the quantities used by the next phase-addition check while its
repeated accepted updates gradually move the continuous state toward a later
phase removal or alternative addition.  This explains why the earlier
instantaneous phase-path gates changed the trajectory but did not recover the
historical lower-Gibbs result.

Accordingly, no immediate phase-ranking acceptance rule is promoted.  A useful
next experiment must use a short, model-independent history of live phase-search
stability--for example, recent assemblage changes, repeated near-degenerate
driving-force margins, active phase amounts approaching removal tolerances, and
recent dependent proposals--rather than comparing only the two instantaneous
candidate rankings.  This is still an adaptive-readiness question, not a change
to the analytical Hessian or its GEM mapping.

### Matched SUBQ S3 ablation: weighted correction is not the FLiBe cause

The FLiBe investigation was performed after correcting the SUBQ S3 term to use
the normalized zeta-weighted pair distribution.  Consequently, both the
alpha-zero reference and adaptive-curvature calculations in the preceding
sections used the corrected production Gibbs-energy formulation.  Alpha zero
disabled the MQMQA GEM curvature correction; it did **not** restore the former
ordinary-pair S3 expression.

A matched four-mode ablation was therefore added to distinguish these two
questions.  A default-false diagnostic switch selected either the corrected
zeta-weighted S3 expression or the former ordinary-pair expression.  The same
selection was applied simultaneously to production partial molars and to the
independent analytical scalar, gradient, and Hessian path.  For each
formulation, an adaptive calculation was compared only with its own alpha-zero
reference:

1. corrected weighted S3 with alpha zero;
2. corrected weighted S3 with adaptive MQMQA curvature;
3. former ordinary-pair S3 with alpha zero; and
4. former ordinary-pair S3 with adaptive MQMQA curvature.

This matched design is necessary because changing S3 changes the thermodynamic
model itself.  The corrected and former alpha-zero results may therefore differ
without identifying a curvature-integration error.  The causal question is
whether restoring the former S3 expression also restores agreement between the
adaptive and alpha-zero paths for that same formulation.

The MSD-TC V4.1 results were:

| LiF fraction | Corrected alpha-zero | Corrected adaptive | Former alpha-zero | Former adaptive | Interpretation |
|---:|---|---|---|---|---|
| 0.45 | `MSFL`, 1267 iterations | `MSFL`, 2884 iterations; scaled dG `2.27e-7` | `MSFL`, 1467 iterations | two condensed phases, 1194 iterations; scaled dG `8.29e-3` | both matched pairs disagree |
| 0.50 | `MSFL`, 1586 iterations | `gas_ideal`, 4720 iterations; scaled dG `4.62e-2` | `MSFL`, 1003 iterations | two condensed phases, 1774 iterations; scaled dG `8.56e-3` | both matched pairs disagree |
| 0.55 | `MSFL`, 3511 iterations | two condensed phases, 1002 iterations; scaled dG `4.67e-3` | `MSFL`, 797 iterations | two condensed phases, 1936 iterations; scaled dG `8.51e-3` | both matched pairs disagree |

Here, scaled dG is the absolute Gibbs-energy difference between each adaptive
result and its formulation-matched alpha-zero reference, divided by
`max(1,abs(G_reference))`.  All twelve calculations were finite and converged,
and both adaptive formulations applied positive curvature.  However, neither
formulation recovered matched adaptive/reference agreement in any of the three
states.  At 0.45 the corrected adaptive calculation retained the same named
phase but still differed in its internal composition and phase amount; this is
not counted as agreement.

The ablation therefore rules out the weighted-S3 correction as the specific
cause of the assessed FLiBe globalization failure.  Restoring the former S3
expression changes the model, iteration histories, and final adaptive basins,
but it does not restore the formulation-matched alpha-zero result.  This result
supports retaining the corrected weighted expression and returning the solver
investigation to phase-search readiness and accumulated trajectory effects.
The legacy switch is diagnostic-only, defaults to false, and changes no normal
Thermochimica behavior.

### Cross-database chloride audit: the FLiBe failure is not generic SUBQ behavior

A second private assessed-database audit tested whether the FLiBe result was a
generic consequence of activating analytical curvature in a nonuniform-zeta
SUBQ phase.  The MSD-TC V4.1 chloride database was exercised at 1000 K and
1 atm for LiCl/MgCl2 mole ratios 45/55, 50/50, and 55/45.  All three untouched
alpha-zero calculations reproducibly converged to the single liquid `MSCL`
phase.  The decoded active model had zeta values from 2.4 to 4.0 and three
G-family interaction records; it had no Q-, B-, or reciprocal-R-family records.
This is therefore cross-database nonuniform-zeta G-family evidence, not B- or
R-family coverage.

| LiCl fraction | Historical iterations | Fixed-alpha iterations and result | Adaptive iterations | Adaptive full/reduced/zero | Active-state agreement | Rank-deficient checks |
|---:|---:|---|---:|---:|---|---:|
| 0.45 | 75 | 1095; `MSCL`, slightly outside the active-state tolerance | 75 | 6 / 0 / 73 | yes | 0 |
| 0.50 | 78 | 2242; `FM3M`, scaled Gibbs difference `2.705e-3` | 74 | 10 / 1 / 67 | yes | 0 |
| 0.55 | 56 | 1119; `FM3M + gas_ideal`, scaled Gibbs difference `2.674e-3` | 56 | 1 / 1 / 61 | yes | 0 |

Here, *active-state agreement* requires the same named assemblage, the
established Gibbs- and phase-amount tolerances, and componentwise mole-fraction
and species-mole agreement over the active phase only.  The distinction matters
for this large database: the largest whole-array differences at 0.45 and 0.55
were stored values belonging to inactive phases, while the active `MSCL` state
agreed.  The untouched repeat calculations remained exactly reproducible.

The result separates three effects.  First, fixed alpha one is not a robust
production policy: it was 14.6 to 28.7 times more expensive than the historical
solve and selected higher-Gibbs nonliquid assemblages in two states.  Second,
the current adaptive policy safely retained the historical `MSCL` solution in
all three states and did not increase the iteration count.  Third, none of the
chloride paths encountered the rank-deficient phase-change trials observed in
FLiBe.  The assessed FLiBe failure is consequently not explained by
nonuniform zeta or SUBQ curvature alone; its active-set topology and transient
phase-search history remain material parts of the mechanism.

This is useful but not a complete MQ-4E exit.  Positive curvature was genuinely
accepted in every chloride state, including full-alpha candidates, but the
final stable full-alpha window was zero in all three runs.  The chloride result
therefore strengthens the safety and cross-database evidence while also showing
that the FeTiVO sustained-full-alpha criterion has not yet generalized to this
assessed salt system.  The private database and driver remain outside the
registered public suite, and no default solver behavior is changed by this
audit.

The zero final-window count should not be read as an absence of sustained
full-curvature behavior.  The complete eligible-solve histories showed:

| LiCl fraction | Longest consecutive full-alpha window | What ended it before convergence |
|---:|---|---|
| 0.45 | 5 solves, global iterations 70-74 | a second eligible solve at iteration 74 failed the residual-magnitude/progress readiness tests; the final solve at 75 also followed a recent assemblage change |
| 0.50 | 9 solves, global iterations 64-72 | iteration 73 accepted alpha 0.1, then readiness was revoked for poor residual progress; the final solve at 74 also followed a recent assemblage change |
| 0.55 | 1 solve at global iteration 55 | the final solve at 56 followed a recent assemblage change and failed the residual-progress test |

Thus 0.45 and 0.50 nearly met the intent of the settled-full-alpha gate, but not
its required endpoint condition: every solve from the final full-alpha window
through the last eligible solve must retain alpha one.  The 0.55 case was not
close by that measure.  Across these runs there were no correction-construction,
correction-ratio, linear-solve, or nonfinite rejections.  The terminal resets
were driven by the current nonlinear-readiness policy and, at some earlier
candidate tests, grouped update/direction safeguards.  This localizes the open
question to readiness retention near convergence rather than to a failed
chloride Hessian or mapper.

The complete late-iteration trace resolved why the endpoint condition was not
met.  In every composition, the last transient phase was removed immediately
before convergence to the single `MSCL` liquid:

| LiCl fraction | Late full-alpha behavior | Final phase event | Convergence |
|---:|---|---:|---:|
| 0.45 | five full-alpha solves at iterations 70--74 | `gas_ideal` removed at 74 | iteration 75 |
| 0.50 | nine full-alpha solves at iterations 64--72; alpha 0.1 at 73 | `gas_ideal` removed at 73 | iteration 74 |
| 0.55 | full alpha selected at iteration 55 | `Mg_L1(liq)` removed at 55 | iteration 56 |

The zero terminal-window metric therefore does not mean that full curvature
was absent late in the solve.  It means that the final phase event correctly
reset readiness, after which the historical convergence test was satisfied
before the new one-phase assemblage could accumulate another settled window.
For assessed cross-database evidence, the event-aware gate consequently
requires either (a) at least three consecutive late full-alpha solves before
the final phase event, or (b) a full-alpha solve that produces the final
assemblage followed by convergence on the next iteration.  All three chloride
states satisfy this gate and agree with the alpha-zero active `MSCL` state.
The original stricter terminal-window requirement remains unchanged for the
FeTiVO MQ-4D checkpoint; this refinement records a distinct late-phase-event
case rather than weakening any trust threshold.

Two tempting ways to manufacture a terminal window were tested and rejected.
First, ordinary convergence was held for at most twelve additional convergence
events while all production readiness, line-search, and phase-search logic
remained active.  The three calculations expanded from 75, 74, and 56
iterations to 304, 515, and 456 iterations, respectively; all exhausted the
hold budget and still ended with zero terminal full-alpha window.  Deferring
the exit re-entered the nonlinear and phase machinery rather than providing a
quiet certification tail.  No such hold is retained in the implementation.

Second, a full-alpha Newton system was reconstructed at the converged state
without applying its update.  The solves were finite and the captured state
was restored exactly, but the maximum raw mixed-variable updates were about
`79.9`, `212.2`, and `164.5`.  These unscaled values cannot support the claim
that an extra full-alpha step would be negligible, and forcing such a step
would bypass the purpose of adaptive readiness.  This replay is retained only
as negative diagnostic reasoning, not as an MQ-4E-B acceptance condition.

Accordingly, the chloride audit completes the bounded MQ-4E-B assessed
globalization evidence: the default-off adaptive path converges efficiently,
uses genuine full curvature, preserves the active assessed liquid state, and
does not encounter the rank-deficient transient topology that blocked FLiBe.
FLiBe remains documented as a separate active-set/globalization limitation;
it is not silently converted into passing evidence.  The private assessed
database and driver remain outside the public registered suite.

### Future evidence and defensible claim framework

The remaining work must build an eventual claim in explicit layers rather than
turn one successful FeTiVO calculation into a universal robustness statement.
The following items are future evidence requirements and are not claims of the
current MQ-4D checkpoint:

1. **Mathematical scope:** show that the implementation follows the supported
   SUBQ equations and variable definitions without depending on the layout or
   parameter values of one database.
2. **Software scope:** retain transactional fallback evidence for invalid,
   singular, nonfinite, boundary, and locally unsafe corrected systems so that
   a rejected curvature correction cannot corrupt the historical GEM solve.
3. **Verification scope:** identify every controlled fixture, production-native
   calculation, and physically assessed database case used as evidence.  Here,
   *assessed* means that the thermodynamic parameters were developed to
   represent a material system using experimental, published, and/or
   first-principles evidence; merely parsing a file through Thermochimica makes
   a case native, not necessarily assessed.
4. **Empirical robustness scope:** extend the completed bounded FeTiVO MQ-4D
   threshold-sensitivity study to representative supported assessed systems.
   Establish that the selected heuristics are not a FeTiVO-specific knife edge
   without claiming global optimality.
5. **Application scope:** demonstrate the selected molten-salt problem over its
   stated database, temperature, composition, phase, and fraction ranges, and
   keep those bounds attached to the resulting application claim.

Before the MQMQA solver work is presented as complete, the audit must contain a
traceable evidence table mapping each layer to its equations, tests, databases,
input ranges, results, and remaining exclusions.  In particular, controlled,
native, and assessed are not mutually exclusive labels: controlled/assessed
describe data provenance, whereas native describes execution through
Thermochimica's production path.

The intended final claim, subject to completion of the listed evidence, is:

> The adaptive MQMQA curvature integration was demonstrated across the tested
> supported systems and remained stable under the reported threshold-sensitivity
> study. Unsafe corrections reverted transactionally to the historical GEM
> system. Universal convergence for arbitrary thermodynamic assessments is not
> claimed.

This framework leaves broad database coverage, the assessed molten-salt case,
and application-range evidence as explicit future tasks rather than implicit
assumptions.
