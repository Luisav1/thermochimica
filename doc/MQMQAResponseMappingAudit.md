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
