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
every SUBQ feature. The documented paper-versus-production `S3` discrepancy is
also not adjudicated here: the response path must follow current production
Thermochimica until that model-definition question is resolved.

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
| Live MQ-4B builder | Do not reuse yet; it explicitly accepts only plain `SUBG` |
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
uniform-zeta FeTiVO SUBQ `G/Q` case. It does not verify SUBQ `B`, `R`,
nonuniform-zeta production data, reduced `deltaA/deltaB` mapping, the live
correction builder, or `GEMNewton` activation. The next SUBQ stage is a
diagnostic reduced-mapping verification analogous to MQ-4A; the existing
production correction builder remains strict plain-`SUBG` until that gate
passes.
