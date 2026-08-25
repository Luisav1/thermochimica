# Hessian Verification and RKMP Integration Evidence

This document records the verification scope and representative numerical
evidence for the local RKMP, CEF, and MQMQA Hessians and for the experimental
plain-RKMP GEM response mapper. MQMQA GEM integration is not part of this work.

## Finite-Difference Standard

All mandatory interior-state sweeps require both an accuracy threshold and an
asymptotic convergence region. Observed order is calculated from adjacent
errors,

```text
p = log(error_coarse/error_fine) / log(h_coarse/h_fine).
```

Scalar comparisons use the scale-normalized error

```text
e_s = |a-b| / max(1, |a|, |b|),
```

and vector comparisons use

```text
e_v = ||a-b||_2 / max(1, ||a||_2, ||b||_2).
```

These definitions behave like relative errors when the compared quantities
are appreciably nonzero, while avoiding division by a nearly zero reference.

The accepted windows are `[0.7, 1.3]` for one-sided forward/backward scalar
second differences, `[1.7, 2.3]` for three-point scalar second
differences and central gradient/partial-molar differences, and `[3.2, 4.8]`
for five-point scalar second differences. At least two consecutive orders must
pass before or at the error minimum. Zero, non-finite, and finite-precision-level
errors are reported as `N/A`.

The automated Fortran assessment owns this pass/fail decision. Figure labels
called `in-range observed p` summarize all individual pre-minimum slopes that
fall inside the relevant window; they are descriptive and do not claim that
the plotting script independently reproduced the stricter consecutive-order
gate.

The plotting code first finds the longest measured in-range convergence region
and anchors each dashed guide directly to one point `(h_a,e_a)` in that region:

```text
e_guide(h) = e_a (h/h_a)^p.
```

There is no arbitrary vertical multiplier. The resulting theoretical
`C h^p` trend is then extrapolated across the complete displayed perturbation
range, including the shaded small-`h` region. This shows how truncation error
would continue decreasing if finite-precision effects did not take over.
Measured curves are solid with markers; extrapolated theoretical trends are
dashed and markerless. Plot limits are computed from measured errors only. A
high-order guide may therefore leave through the bottom of the plot rather than
expanding the axis and compressing the measured evidence.

Gray shading is empirical: it begins at the first observed post-minimum error
increase as `h` decreases. The plot does not uniquely prove the cause of that
increase. Its behavior is consistent with finite-precision cancellation because
the stencil subtracts nearly equal Gibbs-energy values and then divides the
remaining numerical error by `h^2`. Machine epsilon describes the underlying
floating-point resolution; cancellation and division by `h^2` can amplify its
effect until the plotted error is much larger than epsilon. The scientifically
bounded interpretation is therefore **an observed small-`h` error-upturn region
consistent with finite-precision cancellation**.

## Representative Results

| Verification | Best or worst-best scaled error | Structural evidence |
| --- | ---: | ---: |
| Controlled RKMP scalar, forward | `3.2087E-9` | order approximately `1.00`; small-`h` error upturn at `h=1.1290E-6` |
| Controlled RKMP scalar, backward | `3.5040E-9` | order approximately `1.00`; small-`h` error upturn at `h=3.7634E-7` |
| Controlled RKMP scalar, three point | `1.5883E-12` | symmetry `0`; radial residual `2.6531E-17` |
| Controlled RKMP scalar, five point | `3.9900E-15` | observed pre-upturn orders approximately `4.00` |
| Exponent-eight RKMP scalar, seven point | `2.2456E-19` | report-only; two visible orders `6.07`, `6.02`; small-`h` upturn at `h=1.5625E-3` |
| Exponent-eight RKMP scalar, nine point | `2.9911E-21` | report-only; one visible order `8.00`; small-`h` upturn at `h=3.1250E-3` |
| Native RKMP partial molars | `4.6425E-15` worst-best | order unavailable: native low-order polynomial is finite-precision-limited from the coarsest step |
| Controlled CEF scalar, three point | `3.0903E-7` worst-best | symmetry `0`; homogeneity at or below `1.25E-16` |
| Controlled CEF scalar, five point | `1.8582E-9` worst-best | binary, ternary, and coupled families covered |
| Native CEF partial molars | `3.2268E-10` worst-best | controlled and admissible converged states covered |
| Standalone MQMQA scalar, three point | below `2.59E-8` in reported cases | raw symmetry below `5.13E-16` |
| Standalone MQMQA scalar, five point | below `1.03E-10` in reported cases | homogeneity below `6.93E-16` |
| Controlled SUBQ scalar, five point | `3.99E-11` to `1.42E-10` across plotted cases | observed orders `4.00` to `4.09`; nonuniform zeta, chi incidence, and B covered |
| Native MQMQA partial molars | `1.1640E-10` worst-best | symmetry `7.5047E-20`; homogeneity `7.3714E-24` |
| Native SUBQ partial molars | `7.2228E-10` worst-best | direct gradient `2.9756E-17`; all 14 tangent directions show second-order convergence |

The RKMP scalar order fixture uses TestThermo30's parsed plain-RKMP model but
temporarily evaluates its binary interaction at exponent four. This supplies
the sixth directional derivative needed to visibly demonstrate the five-point
fourth-order truncation region. The original database exponent is then restored
and the analytic Hessian is recomputed before the native TestThermo30
partial-molar comparison. The native low-order polynomial already agrees to
approximately floating-point precision at the coarsest retained step, so that separate production-data check
enforces accuracy but does not claim an observed truncation order.

The exponent-four sweep now continues four additional factor-of-three
refinements. Both one-sided formulas reach a minimum and turn upward. Their
signed scale-normalized errors are

```text
s_+ = (D_+ - a) / max(1, |D_+|, |D_-|, |a|),
s_- = (D_- - a) / max(1, |D_+|, |D_-|, |a|),
```

where `a = v^T H v`. In the shared first-order region, the signs are opposite
and the magnitude ratio approaches one, as predicted by the retained
`+h f'''` and `-h f'''` terms. A separate symmetric-log plot shows the later
sign changes inside the observed small-`h` upturn region rather than hiding
them with absolute values.

For an exponent `m`, the controlled constant-total-moles RKMP energy has
directional polynomial degree `m+2`: the binary mixing factor is quadratic in
`h`, while the composition contrast contributes degree `m`. Exponent four is
therefore degree six, so seven- and nine-point centered formulas would be exact
apart from finite-precision effects. A separate exponent-eight fixture is degree ten and has
the nonzero eighth- and tenth-derivative terms needed to expose sixth- and
eighth-order truncation. Its three-, five-, seven-, and nine-point tables report
orders, best errors, and small-`h` error-upturn points. Seven- and nine-point results remain
report-only because the seven-point curve has only two visible sixth-order
intervals and the nine-point curve reaches the conservative double-precision
floor after one visible eighth-order interval.

The RKMP/CEF truncation-coefficient figure plots `e(h)/h^p` only over each
measured expected-order region. A near-horizontal segment demonstrates
`e(h) approximately C h^p`. Plateau heights are not expected to match: `C`
depends on each model's higher derivatives, thermodynamic state, and error
normalization. The diagnostic therefore tests the power-law behavior without
misreading vertical curve separation as Hessian disagreement.

For presentation, the controlled RKMP scalar-energy figure is the primary
stencil-order result because it displays the expected first-, second-, and
fourth-order regions. The native RKMP production partial-molar figure is backup
accuracy evidence: it begins near the finite-precision limit and therefore rises as the plotted
step moves left toward smaller `h`.

The remaining figures use the same visual language without forcing every model
into the RKMP stencil experiment:

- The controlled CEF figure shows three- and five-point scalar-energy checks,
  including the small-`h` region where a post-minimum error upturn is measured.
- The standalone MQMQA figure shows controlled nonmagnetic `SUBG` `G`, `Q`, and
  `B` cases. It verifies the standalone scalar forms and derivative propagation,
  not native database decoding.
- The standalone SUBQ figure shows controlled nonuniform-zeta configurational,
  mixed-environment chi-incidence, and `B` cases. It verifies the staged SUBQ
  scalar forms and derivative propagation independently of the existing SUBG
  figure. It is not yet a native `FeTiVO.dat` comparison.
- The native MQMQA figure compares analytic Hessian-vector products with
  central finite differences of production partial molars along five
  total-preserving directions at a converged `CuFeC-Kang.dat` state. It is not
  shaded because no small-`h` upturn is observed in the retained sweep.
- The native SUBQ figure performs the production-partial-molar comparison along
  all 14 independent total-preserving directions at a fixed interior
  composition of the assessed `FeTiVO.dat` `SlagBsoln` phase. Its scope is
  uniform-zeta reference/configurational/G/Q behavior, not every SUBQ feature.
- Every figure reports descriptive in-range observed orders and best scaled
  errors in the plot, while the Fortran executables remain the sole owners of
  automated acceptance.

The converged ALABANDITE CEF state has minimum site fraction
`1.382909E-4`. It is reported as boundary-adjacent, so scalar-order enforcement
is skipped for that state without clipping. Its complete sweep is retained and
the production partial-molar comparison remains mandatory. Controlled
interior CEF states enforce scalar order.

## RKMP GEM Mapping

`TestRKMPGEMMappingVerification` contains two deliberately separate layers.
The structural layer reconstructs the normalized ideal and ideal-plus-excess
responses with the same local linear-response mathematics as the mapper; this
catches implementation, indexing, and unit errors. The independent layer
perturbs species-level thermodynamic driving forces, re-solves the nonlinear
normalized phase composition with the established production RKMP
partial-molar routine and a finite-difference Newton Jacobian, and then
condenses the observed response into the reduced GEM objects.

| Mapping check | Scaled residual or error |
| --- | ---: |
| Ideal GEM response reconstruction | `2.220446E-17` |
| Mole-fraction normalization | `8.881784E-16` |
| Local stationarity response | `1.082467E-15` |
| Independent mixed-direction `dx/dGamma` FD | `1.819768E-9` |
| Independent full `dx/dGamma` FD | `1.192622E-8` |
| Independent reduced matrix response, delta A | `1.283343E-8` |
| Independent reduced RHS response, delta B | `1.740565E-9` |
| Maximum absolute mapped delta B | `1.284228E2` |
| Structural mapper reconstruction, delta A | `1.069096E-16` |
| Structural mapper reconstruction, delta B | `0.000000E+00` |

The nonlinear response sweep spans `h=4.0E-2` to `4.8828125E-6` and decreases
monotonically from `8.137236E-2` to `1.819768E-9`, providing a visible
second-order region rather than relying on one favorable perturbation. Its
report prints the observed order beside every adjacent `h`/error pair; after
the coarsest nonlinear step, the retained orders approach `2.00`. The nonzero
delta-B magnitude confirms that the right-hand-side comparison is not a
vacuous zero-vector test.

The test also verifies alpha-zero and no-active-phase no-op behavior, direct
rejection of a singular bordered local system, and phase-specific propagation
of a controlled non-finite local residual response. The mapper never applies a
partially assembled correction after one active phase fails.

## Curvature-Enabled Regression

The supported Fortran control is:

```fortran
call SetRKMPHessianControls(enable, alpha_max, debug, info)
```

For process-level testing, the equivalent inputs are:

```text
THERMOCHIMICA_RKMP_EXACT=1
THERMOCHIMICA_RKMP_ALPHA_MAX=1
```

`alpha_max` is an upper trust bound on the completed `delta A`, `delta B`
correction, not derivative scaling. For example, `alpha_max=0.1` tries
`0.1`, `0.01`, `0.001`, then zero.

The exact-on regression records:

| Case | Iterations | Positive / zero-alpha iterations | Maximum alpha | Mass residual |
| --- | ---: | ---: | ---: | ---: |
| TestThermo30, alpha max 1 | 21 | 15 / 6 | 1 | `1.9762E-14` |
| TestThermo33, alpha max 1 | 488 | 35 / 453 | 1 | `9.4357E-9` |
| TestThermo90, alpha max 1 | 438 | 67 / 256 | 1 | `1.7319E-14` |
| TestThermo30, alpha max 0.1 | 21 | 15 / 6 | 0.1 | `1.1102E-15` |
| TestThermo30, tighter trust gates | 21 | 15 / 6 | 1 | `1.9762E-14` |
| TestThermo30, looser trust gates | 21 | 15 / 6 | 1 | `1.9762E-14` |

TestThermo90 contains transient RKMP activity but converges to its expected
pure-phase assemblage; iterations without an active RKMP phase are not included
in either alpha count. Every listed case preserves its original thermodynamic
output checks and reports zero local mapper failures.

The complete suite passes `70/70` with defaults off and `70/70` with requested
alpha one. A final default-off run is required after any exact-on experiment.
These regressions demonstrate that the curvature-enabled path produces stable,
converged Thermochimica solutions under the tested configurations. They do not
establish that exact curvature reduces iteration count or wall-clock time.

## Coverage Boundaries

- RKMP native evidence covers converged plain RKMP and the established
  production partial-molar implementation. It excludes RKMPM magnetism and
  other solution models.
- Standalone CEF evidence covers nonmagnetic plain-SUBL binary, ternary, and
  coupled interaction families. CEF is not mapped into GEMNewton.
- Standalone MQMQA evidence covers controlled nonmagnetic SUBG reference,
  configurational, G, Q, B, and traced ternary formulas.
- Standalone SUBQ evidence covers controlled nonmagnetic SUBQ configurational
  exponents, mixed-environment chi weights, pair-specific zeta propagation
  through weighted-pair quantities, and B-family curvature. The corrected
  weighted-versus-legacy-unweighted S3 diagnostic follows published Eqs. (5),
  (6), (16), (29), and (31). Production and both standalone evaluators now use
  the weighted definition for SUBQ. The resulting complete-model derivative
  differences from the legacy definition are modest but resolved and
  state-dependent in the controlled states tested.
- Native MQMQA Hessian evidence covers nonmagnetic plain SUBG reference,
  configurational, and G-family behavior from `CuFeC-Kang.dat`.
- Native SUBQ scalar-energy evidence covers reference/configurational and G/Q
  behavior for the uniform-zeta `SlagBsoln` phase from `FeTiVO.dat`. Native SUBQ
  MQ-2B evidence compares the standalone unconstrained gradient directly with
  complete production partial molars, then compares the analytic Hessian with
  order-aware finite differences of those partial molars at a fixed interior
  composition. The direct gradient error is `2.9756E-17`; the worst-best
  Hessian-vector error is `7.2228E-10`; the worst componentwise scaled error at
  the normwise-best steps is `5.5123E-10`, with second-order convergence in all
  14 independent total-preserving directions. This is a native assessed SUBQ
  G/Q case, not native verification of every SUBQ feature.
- A controlled modified-runtime production regression verifies the corrected
  nonuniform-zeta SUBQ S3 selection, including the derivative of the weighted
  normalization, and restores FeTiVO's parsed zeta row exactly. At its positive
  interior state, the corrected reference/configurational, complete excess, and
  complete total parity errors are `1.47E-16`, `8.33E-17`, and `1.47E-16`.
  The complete-gradient normwise error is `6.66E-17`, and the worst componentwise
  scaled gradient error is `5.95E-17`. The corrected production block remains
  resolved from the reconstructed legacy block by `1.42E-4`. In the unchanged
  assessed uniform-zeta state, the maximum direct weighted/ordinary pair-fraction
  difference is `1.11E-16`, making explicit why that database cannot distinguish
  the two S3 definitions.
  Private local MSD-TC FLiBe evidence now provides database-native
  nonuniform-zeta `G`-family scalar, gradient, and Hessian verification across
  three compositions. Its database-dependent assets remain excluded from the
  public regression suite. Native `B`-family parameters remain controlled
  standalone coverage; reciprocal `R` terms and magnetism remain outside the
  demonstrated scope.

## Reproduction

```bash
docker exec thermochimica bash -lc \
  "cd /work && make clean && make test && ./run_tests"

docker exec thermochimica bash -lc \
  "cd /work && env THERMOCHIMICA_RKMP_EXACT=1 \
   THERMOCHIMICA_RKMP_ALPHA_MAX=1 ./run_tests"

# Plotting is optional and requires the Debian Matplotlib package:
docker exec -u root thermochimica bash -lc \
  "apt-get update && apt-get install -y python3-matplotlib"

docker exec thermochimica bash -lc \
  "cd /work && python3 scripts/plot_hessian_verification.py"
```

The plotting utility writes complete reports, vector SVG figures, and 300-DPI
PNG figures to `outputs/hessian_verification/`. Matplotlib is required only for
this optional presentation layer; the automated Fortran tests do not depend on
plotting software.
