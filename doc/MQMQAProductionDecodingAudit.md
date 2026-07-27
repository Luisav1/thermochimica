# MQ-2A: MQMQA Production Decoding Audit

## Scope

This document freezes the translation from Thermochimica's current plain
`SUBG` production state into the generic inputs accepted by
`ModuleMQMQAUnconstrained`.

MQ-2A is a source audit only. It does not add a production adapter, call the
new Hessian from a Thermochimica calculation, modify `GEMNewton`, or support
`SUBQ`, reciprocal `R` terms, or magnetic energy.

The production data flow is:

1. `ParseCSDataBlock.f90` reads the plain-`SUBG` zeta value and pair counts.
2. `ParseCSDataBlockSUBG.f90` reads constituents, chemical groups,
   quadruplet topology, coordination numbers, and excess parameters.
3. `CheckSystemExcess.f90` removes constituents and quadruplets that are not
   present in the selected chemical system and renumbers the survivors.
4. `CompThermoData.f90` evaluates temperature-dependent reference and excess
   coefficients and converts them to Thermochimica's dimensionless units.
5. `CompExcessGibbsEnergySUBG.f90` consumes the resulting runtime arrays.

A future adapter must consume the **runtime arrays after steps 3 and 4**. It
must not mix original parser indices with filtered runtime indices.

## What Zeta Means

Zeta is part of the MQMQA model, not an adapter invented by Thermochimica.
Thermochimica describes it as the first-nearest-neighbour to
second-nearest-neighbour ratio. In the energy equations it changes the
statistical weight of an A-X pair:

```text
weighted A-X amount = ordinary A-X amount / zeta(A,X)
```

The weighted pair distribution is used by the configurational `S2` term and
the production `B` interaction family.

For plain `SUBG`, the current database format reads one zeta value and applies
it to every A-X pair. The generic module stores `dZeta(A,X)` as a matrix so the
mathematics remains explicit. The production adapter will fill every matrix
entry from the corresponding runtime A-X pair record and will verify that all
entries are present and positive. Pair-dependent zeta values occur in
Thermochimica's `SUBQ` format, which is outside the current scope.

## Phase and Array Boundaries

For a candidate production phase `iSolnIndex`, the adapter must first require:

```fortran
cSolnPhaseType(iSolnIndex) == 'SUBG'
```

The relevant runtime ranges are:

```text
iSPI   = iPhaseSublattice(iSolnIndex)
iFirst = nSpeciesPhase(iSolnIndex - 1) + 1
iLast  = nSpeciesPhase(iSolnIndex)
nQuad  = nPairsSRO(iSPI, 2)
```

The quadruplet species are stored in `iFirst:iLast`, and the audit requires
`iLast - iFirst + 1 == nQuad`.

## `MQMQAModelData` Mapping

| Generic field | Thermochimica runtime source | Required translation |
| --- | --- | --- |
| `nSublattice1` | `nConstituentSublattice(iSPI,1)` | Copy directly. |
| `nSublattice2` | `nConstituentSublattice(iSPI,2)` | Copy directly. |
| `iQuadruplet(q,1:2)` | `iPairID(iSPI,q,1:2)` | Copy the first-sublattice A and B indices. |
| `iQuadruplet(q,3:4)` | `iPairID(iSPI,q,3:4)` | Subtract `nSublattice1`; production stores second-sublattice X and Y indices after the first-sublattice index range. |
| `dCoordination(q,:)` | `dCoordinationNumber(iSPI,q,:)` | Copy the finalized A, B, X, and Y coordination numbers. |
| `dZeta(a,x)` | `dZetaSpecies(iSPI,m)` | Find `m` for which `iConstituentSublattice(iSPI,1,m)==a` and `iConstituentSublattice(iSPI,2,m)==x`. |
| `dReferenceEnergy(q)` | `dStdGibbsEnergy(iFirst+q-1)` | Copy the temperature-evaluated, dimensionless quadruplet reference energy. |

The zeta lookup covers `m=1:nPairsSRO(iSPI,1)`, the retained A-X endmember
pairs. Every `(a,x)` cell must be assigned exactly once.

The generic module does not need the following arrays for local energy
evaluation:

- `dSublatticeCharge`: used while the parser constructs missing coordination
  numbers and quadruplet stoichiometry;
- `dConstituentCoefficients`: used while pair reference energies are converted
  into quadruplet reference energies;
- `dStoichSpecies`: needed later for element coupling and constrained GEM
  response, but not for MQ-2A local energy decoding.

By the time the adapter reads `dCoordinationNumber` and `dStdGibbsEnergy`,
those setup calculations have already been completed.

## Local Composition and Phase Amount

`CompExcessGibbsEnergySUBG` evaluates the local model using
`dMolFraction(iFirst:iLast)`. The generic module instead accepts extensive
quadruplet amounts and constructs its own normalized fractions.

For an active phase in assemblage slot `k`, the production translation is:

```text
x_q = dMolFraction(iFirst+q-1) / sum(dMolFraction(iFirst:iLast))
n_q = dMolesPhase(k) * x_q
```

where `iAssemblage(k) == -iSolnIndex`.

For a scale-independent local scalar comparison, `n_q=x_q` is sufficient.
For a Hessian with respect to physical quadruplet mole amounts, the phase
amount must be retained because an extensive-energy Hessian scales as the
inverse phase amount.

Only strictly positive interior states are admissible to the current generic
Hessian. A future native test must report and skip a boundary state rather
than clipping its composition.

## Units

Thermochimica evaluates both the standard quadruplet energies and excess
coefficients at the current temperature, then multiplies them by
`1/(R*T)`. The local production routine consequently works in dimensionless
energy units.

The native adapter must therefore use:

```text
dIdealScale = 1
dReferenceEnergy(q) = dStdGibbsEnergy(iFirst+q-1)
dCoefficient        = dExcessGibbsParam(parameterIndex)
```

The generic outputs will then be dimensionless and directly comparable with
the nonmagnetic outputs of `CompExcessGibbsEnergySUBG`. Multiplying the total
energy, gradient, and Hessian by `R*T` converts them to physical energy units;
that conversion must not be applied to only one energy block.

## `MQMQAInteractionTerm` Mapping

For each retained parameter

```text
parameterIndex =
    nParamPhase(iSolnIndex-1)+1 : nParamPhase(iSolnIndex)
```

the runtime mapping is:

| Generic field | Thermochimica runtime source | Translation |
| --- | --- | --- |
| `iFamily` | `cRegularParam(parameterIndex)` | `G`, `Q`, and `B` map to their public module constants. Reject `R` and every other label. |
| `iA` | `iRegularParam(parameterIndex,2)` | Copy the filtered first-sublattice index. |
| `iB` | `iRegularParam(parameterIndex,3)` | Copy the filtered first-sublattice index. |
| `iX` | `iRegularParam(parameterIndex,4)` | Subtract `nSublattice1`. |
| `iY` | `iRegularParam(parameterIndex,5)` | Subtract `nSublattice1`. |
| `iExponentP` | `iRegularParam(parameterIndex,6)` | Copy; require a nonnegative integer. |
| `iExponentQ` | `iRegularParam(parameterIndex,7)` | Copy; require a nonnegative integer. |
| `iExponentR` | `iRegularParam(parameterIndex,8)` | Copy the supported first-sublattice ternary order. |
| unused fourth exponent | `iRegularParam(parameterIndex,9)` | Production reads this value but does not use it. MQ-2 support requires it to be zero instead of silently discarding a nonzero value. |
| `iTernaryConstituent` | `iRegularParam(parameterIndex,10)` | Copy the filtered first-sublattice constituent index. Zero means binary. |
| unsupported second-sublattice ternary | `iRegularParam(parameterIndex,11)` | Production currently reads this index but does not use it in the scalar/partial-molar implementation. Require zero. |
| `dCoefficient` | `dExcessGibbsParam(parameterIndex)` | Copy the temperature-evaluated dimensionless coefficient. |

`iRegularParam(:,1)` records the database parameter arity. It is useful for
adapter validation but does not replace the explicit topology and ternary
checks above.

The parser accepts only `G`, `Q`, `R`, and `B`. Although the production
evaluator contains an `H` alias in the `B` branch, the current `SUBG` parser
rejects `H`; MQ-2 must not expose an unreachable parser branch.

## Asymmetric Group Masks

The generic `lGroup1` and `lGroup2` masks are fixed model data for one
interaction. "Fixed" means they are derived from database metadata and do not
change with composition.

For a first-sublattice A-B interaction in a fixed X-X environment:

1. Set A in group 1 and B in group 2.
2. If a matching interpolation override exists for the A-B-third-constituent
   triad, use the override's constant constituent to assign the third
   constituent.
3. Otherwise, when A and B have different `iChemicalGroup` values, assign
   constituents sharing A's chemical group to group 1 and constituents
   sharing B's chemical group to group 2.
4. Constituents matching neither group remain in neither mask.

For a second-sublattice X-Y interaction in a fixed A-A environment, the
current production routine constructs the masks from
`iChemicalGroup(iSPI,2,:)`. It does not apply interpolation overrides in that
branch.

Supported ternary decoding is currently narrower than the parser storage:

- only `iRegularParam(:,10)>0` is traced;
- it must be a first-sublattice ternary with X equal to Y;
- the third constituent can be in group 1, group 2, or neither;
- `iRegularParam(:,11)>0` is rejected as untraced.

## Explicit Scope Rejections

The future MQ-2 adapter must reject, not approximate:

- `cSolnPhaseType == 'SUBQ'`, whose configurational exponents and pair
  treatment differ from plain `SUBG`;
- production-family `R`, whose extensive scalar energy remains untraced;
- parser-unreachable `H`;
- second-sublattice ternary index slot 11;
- a nonzero fourth exponent in slot 9;
- negative exponents;
- missing or duplicate A-X zeta records;
- nonpositive coordination or zeta values;
- boundary compositions for the analytic Hessian;
- magnetic contributions, which are assembled outside
  `CompExcessGibbsEnergySUBG`.

## Native Verification Fixture Audit

An exact search of all tracked database files finds these MQMQA parameter
records:

| Database | Phase model | `G` records | `Q` records | `B` records | `R` records |
| --- | --- | ---: | ---: | ---: | ---: |
| `data/CuFeC-Kang.dat` | two `SUBG` phases | 22 | 0 | 0 | 0 |
| `data/ClAlNa.dat` | one `SUBQ` phase | 7 | 0 | 0 | 0 |
| `data/FeTiVO.dat` | one `SUBQ` phase | 15 | 8 | 0 | 0 |

No other tracked database contains a `G`, `Q`, `B`, or `R` MQMQA parameter
record. The locally added `data/bergeron_Th-U-Pu-O_.dat` likewise contains no
`SUBG`, `SUBQ`, or MQMQA parameter records.

Therefore, `Q` is not absent from Thermochimica's available databases:
`FeTiVO.dat` supplies eight real `Q` records, and Tests 57--60 exercise its
stable `SlagBsoln` `SUBQ` phase. Those records do not provide a fixture for
the present adapter because MQ-2 is intentionally restricted to plain
`SUBG`; `SUBQ` uses different configurational exponents, pair treatment, and
zeta storage. This is a scope distinction, not evidence that the `Q` family
is unavailable in Thermochimica.

The parser accepts `R`, and the production evaluator contains an `R` branch,
but the database search finds zero `R` records. The evaluator resets `dGex`
to zero at the start of each parameter and does not assign an `R` expression
before adding `dGex`; consequently, the current `R` branch is a no-op. Thus
`R` is recognized by the file format but is neither thermodynamically
implemented by this routine nor database-covered in this checkout. It remains
explicitly rejected by the disconnected Hessian. This project will not try to
resurrect or infer an intended `R` formulation.

Within the narrower plain-`SUBG` scope, the repository contains one database:

```text
data/CuFeC-Kang.dat
```

`TestThermo56` runs this database at 1400 K and confirms an active Liquid
phase. That Liquid is a suitable first native fixture for:

- parser-to-runtime topology;
- finalized coordination numbers;
- the uniform plain-`SUBG` zeta value (`2.4` in this database);
- reference and configurational energy;
- binary first-sublattice `G` parameters;
- first-sublattice ternary `G` parameters whose third constituent belongs to
  neither asymmetric group.

Its plain-`SUBG` phases do **not** cover:

- `Q` parameters;
- `B` parameters;
- second-sublattice G/Q interactions;
- ternary group-1 and group-2 branches;
- a nonuniform zeta matrix, which belongs to `SUBQ`, not plain `SUBG`.

The added assessed Bergeron database,
`data/bergeron_Th-U-Pu-O_.dat`, does not contain `SUBG` or `SUBQ`. Its
solution phases use `IDMX`, `SUBI`, `RKMP`, and `SUBL`, so it cannot provide
native MQMQA parameter coverage.

MQ-2B should use `CuFeC-Kang.dat` for the first converged plain-`SUBG`
production comparison. It must not copy `SUBQ` parameters into a `SUBG`
fixture or invent `Q` or `B` coefficients and present them as physical or
thermodynamic evidence. Until an authoritative assessed plain-`SUBG` database
containing those parameter families is available, the evidence must be
reported in two distinct categories:

- `Q` and `B` mathematical implementation: covered by the independent
  standalone MQ-1 tests;
- plain-`SUBG` production parser/data decoding for `Q` and `B`: not yet
  demonstrated with an assessed database. Existing `SUBQ` `Q` coverage is a
  separate future adapter target.

A search for additional native plain-`SUBG` `Q` or `B` cases is not a gate for
MQ-2B. Reassess whether that evidence is needed only after the standalone
SUBG Hessian and the available native `G` comparison both pass.

A synthetic parser fixture could test software mechanics only if it were
clearly labelled as nonphysical test data. It is not required for MQ-2B and
must never be described as thermodynamic verification. The current project
decision is to defer such a fixture and preserve the production-data coverage
gap explicitly.

## MQ-2A Exit Decision

All fields required by `MQMQAModelData` and `MQMQAInteractionTerm` have a
specific runtime source and conversion rule. The current plain-`SUBG` units
and filtered-index boundary are identified, and unsupported production
branches have explicit rejection rules.

## MQ-2B Exit Result

MQ-2B implements the audited conversion in
`ModuleMQMQAProductionAdapter.f90` and exercises it in
`TestMQMQANativeHessianVerification.F90`. The test converges the TestThermo56
state from `CuFeC-Kang.dat`, decodes its active plain-`SUBG` Liquid, and
compares the disconnected MQ-1 implementation with production Thermochimica.

The native `G` case passes all intended gates:

- reference plus configurational scalar energy agrees exactly at the reported
  precision;
- excess scalar energy agrees to `3.47E-18` normalized error;
- the analytic gradient agrees with production partial molars to `1.22E-15`;
- finite differences of production partial molars agree with Hessian-vector
  products to a worst best normalized error of `1.16E-10`;
- raw symmetry and homogeneity residuals are `7.50E-20` and `7.37E-24`.

These results verify the available assessed plain-`SUBG` `G` path from parsed
data through runtime decoding and local derivatives. They do not create
database-backed evidence for `Q` or `B`, extend the result to `SUBQ`, or
connect the Hessian to `GEMNewton`.

The next mathematical stage is the constrained local MQMQA composition
response: determine how quadruplet amounts respond to perturbations while
preserving the model's local normalization and topology identities. That
response should be verified independently before deriving any reduced GEM
matrix contribution.
