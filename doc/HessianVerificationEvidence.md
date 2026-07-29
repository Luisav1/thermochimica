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

The accepted windows are `[1.7, 2.3]` for three-point scalar second
differences and central gradient/partial-molar differences, and `[3.2, 4.8]`
for five-point scalar second differences. At least two consecutive orders must
pass before or at the error minimum. Zero, non-finite, and roundoff-level
errors are reported as `N/A`.

## Representative Results

| Verification | Best or worst-best scaled error | Structural evidence |
| --- | ---: | ---: |
| Controlled RKMP scalar, three point | `1.5883E-12` | symmetry `0`; radial residual `2.6531E-17` |
| Controlled RKMP scalar, five point | `3.9900E-15` | observed pre-roundoff orders approximately `4.00` |
| Native RKMP partial molars | `4.6425E-15` worst-best | order unavailable: native low-order polynomial begins at roundoff |
| Controlled CEF scalar, three point | `3.0903E-7` worst-best | symmetry `0`; homogeneity at or below `1.25E-16` |
| Controlled CEF scalar, five point | `1.8582E-9` worst-best | binary, ternary, and coupled families covered |
| Native CEF partial molars | `3.2268E-10` worst-best | controlled and admissible converged states covered |
| Standalone MQMQA scalar, three point | below `2.59E-8` in reported cases | raw symmetry below `5.13E-16` |
| Standalone MQMQA scalar, five point | below `1.03E-10` in reported cases | homogeneity below `6.93E-16` |
| Native MQMQA partial molars | `1.1640E-10` worst-best | symmetry `7.5047E-20`; homogeneity `7.3714E-24` |

The RKMP scalar order fixture uses TestThermo30's parsed plain-RKMP model but
temporarily evaluates its binary interaction at exponent four. This supplies
the sixth directional derivative needed to visibly demonstrate the five-point
fourth-order truncation region. The original database exponent is then restored
and the analytic Hessian is recomputed before the native TestThermo30
partial-molar comparison. The native low-order polynomial is already exact to
roundoff at the coarsest retained step, so that separate production-data check
enforces accuracy but does not claim an observed truncation order.

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
- Native MQMQA evidence covers nonmagnetic plain SUBG reference,
  configurational, and G-family behavior from `CuFeC-Kang.dat`.
- Native SUBQ, Q-family parameters, B-family parameters, reciprocal R terms,
  magnetism, constrained MQMQA response, and MQMQA GEM integration remain
  outside the demonstrated scope.

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
