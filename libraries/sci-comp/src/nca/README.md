# NCA — Non-Compartmental Analysis core

Pure-math, stateless computation core for Non-Compartmental Analysis.
Operates on plain `Float64Array` / `Uint8Array` inputs.

| Capability | Scope |
|---|---|
| **Parameters** | Cmax, Tmax, AUClast, AUCinf, %AUCextrap, λz, t½, CL (or CL/F), Vz (or Vz/F), AUMClast, AUMCinf, MRT, Vss (IV only), Tlag (extravascular only), %AUMCextrap |
| **Routes**     | IV bolus (with `c0` back-extrapolation), IV infusion (with the `T_inf/2` MRT correction), extravascular (PO/SC/IM treated uniformly) |
| **AUC / AUMC methods** | linear, log-linear, linear-up/log-down — each as both naive Float64 and Neumaier-compensated summation |
| **BLQ rules**  | `set-zero`, `set-half-lloq`, `exclude`, `missing` — per-phase configurable (`preFirstMeasurable`, `embedded`, `afterLast`, `consecutiveAfterLast`) |
| **λz**         | auto best-fit (subset search; PKNCA/WinNonlin flat-tolerance tie-break — most points within `adjRSquaredFactor` of the global-max adjusted R²) + manual point selection |
| **Sparse sampling** | composite AUClast with a Holder standard error, Nedelman-Jia df and Student-t CI; stratified bootstrap with a BCa interval for the nonlinear parameters |

## Quickstart

```typescript
import {nca} from '@datagrok-libraries/sci-comp';

const inputs: nca.ProfileInputs = {
  time:    new Float64Array([0, 0.25, 0.5, 1, 2, 4, 8, 12]),
  conc:    new Float64Array([0, 1.5,  2.4, 3.0, 1.8, 0.9, 0.3, 0.1]),
  blqMask: new Uint8Array(8),
  lloq:    0.05,
  dose:               2.5,
  doseUnits:          'mg',
  concentrationUnits: 'mg/L',
  timeUnits:          'h',
  route:              nca.ROUTE_PO,         // 0=IV-bolus, 1=IV-infusion, 2=PO, 3=SC, 4=IM, 5=other
  infusionDuration:   null,
  bodyWeight:         null,
};

const rules: nca.NcaRules = {
  aucMethod: 'linear-up-log-down',
  blq: {
    preFirstMeasurable: 'set-zero',
    embedded:           'set-zero',
    afterLast:          'set-zero',
    consecutiveAfterLast: 'set-zero',
  },
  lambdaZ: {
    mode:              'auto-best-fit',
    minPoints:         3,
    minRSquared:       0.85,
    excludeCmax:       true,
    adjRSquaredFactor: 1e-4,
    // minSpanRatio:   2,   // optional; diagnostic only — see below
  },
  extrapWarnPct:        20,
  extrapErrorPct:       50,
  extrapWarnPctAumc:    20,
  compensatedSummation: false,
};

const result = nca.computeNca(inputs, rules);

// result.values.{cmax, tmax, aucLast, aucInf, pctExtrap, lambdaZ, halfLife, cl, vz,
//                aumcLast, aumcInf, mrt, vss, tlag, pctExtrapAumc}
// result.provenance.lambdaZ      → LambdaZResult: pointsUsed, R², tStart, tEnd, intercept, spanRatio
// result.provenance.blqApplied   → which points were modified by BLQ rules
// result.provenance.warnings     → AUC_EXTRAP_HIGH / AUMC_EXTRAP_HIGH / LAMBDAZ_FEW_POINTS /
//                                  LAMBDAZ_LOW_SPAN / BLQ_HIGH_FRACTION
// result.status                  → 'ok' | 'partial' | 'failed'
```

### Terminal-phase span ratio

`provenance.lambdaZ.spanRatio` is `(tEnd − tStart) / halfLife` — how many half-lives
the fitted terminal window actually covers. PKNCA reports the same quantity as
`span.ratio`, and its convention flags values below 2.

**It is a diagnostic, never a gate.** Setting `lambdaZ.minSpanRatio` only makes
`computeNca` emit a `LAMBDAZ_LOW_SPAN` warning; the selected fit is identical with the
threshold set and unset. The full contract — and why PKNCA behaves the same way — is on
the `LambdaZResult.spanRatio` TSDoc.

`minSpanRatio` deliberately defaults to `undefined` rather than to PKNCA's 2: at 2, four
of the twenty-seven profiles in the reference corpus would emit a warning with no consumer
opt-in.

Why it earns a place next to adjusted R²: adj-R² cannot discriminate at small `n`. With
three points and one residual degree of freedom, a slope resting on 0.69 half-lives still
reports adj-R² ≈ 0.994 (`02_indometh` subject 1 — the corpus worst case, asserted in
`reference-suite.test.ts`). Span ratio is the only statistic in the PKNCA output set that
separates that fit from a well-characterised one.

## Sparse / destructive sampling

When each animal contributes one sample (destructive) or a few (batch), per-subject NCA
is undefined and `computeNca` does not apply. `sparseAuc` instead integrates a composite
mean profile and reports an honest interval for it — Bailer's variance generalised to
Holder's covariance estimator, with a Nedelman-Jia correlated Satterthwaite df.

```typescript
const sparse: nca.SparseInput = {
  //                              t=0.5          t=1          t=2          t=4
  nominalTime: new Float64Array([0.5, 0.5, 0.5,  1, 1, 1,    2, 2, 2,    4, 4, 4]),
  conc:        new Float64Array([1.8, 2.1, 1.6,  2.9, 3.2, 2.7,  1.7, 1.9, 1.5,  0.6, 0.8, 0.5]),
  blqMask:     new Uint8Array(12),
  lloq:        0.05,
  // Destructive design: every animal appears exactly once.
  // Pass `null` to assume destructive (Bailer, r_ij = 0) — you get a warning saying so.
  animalId:    new Int32Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]),
};

const r = nca.sparseAuc(sparse, {ciLevel: 0.95});
// r.auc       → 6.3      composite AUClast (linear trapezoidal, anchored at t=0, C=0)
// r.se        → 0.2344   Holder standard error
// r.df        → 5.343    Nedelman-Jia correlated Satterthwaite df
// r.ci        → [5.709, 6.891]   Student-t interval at ciLevel
// r.topology  → 'destructive'    derived from the animal × time table, not the label
// r.composite → mean / SD / %CV / n / nBlq per nominal timepoint
// r.warnings  → SPARSE_TOPOLOGY_MISMATCH / SPARSE_TERMINAL_OVEREST /
//               SPARSE_VARIANCE_MODELED / SPARSE_DF_UNAVAILABLE /
//               SPARSE_ANIMAL_ID_ABSENT_ASSUMED_DESTRUCTIVE
```

`sparseAuc` is closed-form and applies only to AUClast. The **nonlinear** parameters
(Cmax, t½, CL) have no closed-form sparse variance — bootstrap those instead:

```typescript
const boot = nca.summarizeBootstrap(sparse, (s) => nca.sparseAuc(s).auc,
  {iterations: 1000, masterSeed: 42});
// boot.ci → [5.95, 6.70]   boot.ciMethod → 'BCa'   boot.suppressed → false
```

The bootstrap **suppresses itself** rather than returning a fake interval when the design
is too small to resample honestly — `h = Π_i C(2n_i − 1, n_i) ≤ 360` distinct stratified
resamples, or any timepoint below `minNPerTimepoint` (default 3). Check `suppressed` and
fall back to the closed-form CI, which is never gated; `suppressReason` says which gate
fired. At a fixed `masterSeed` the resampling is deterministic.

`buildCompositeProfile(input, blqRule?)` returns just the per-timepoint table when you
want the composite profile without the interval.

## Direct kernel access

Every kernel behind `computeNca` and `sparseAuc` is also exported individually — pure
functions, each tested in isolation — for finer-grained control or a custom pipeline:
BLQ (`applyBlqStrategy`), peak (`findCmax`), back-extrapolation (`estimateC0`,
`insertC0`), terminal slope (`lambdaZBestFit`, `lambdaZManual`), integration (`auc*`,
`aumc*`, `neumaierSum`), derived parameters (`halfLifeFromLambdaZ`, `clearance`,
`volumeTerminal`, `pctExtrapolated`, `meanResidenceTime`, `volumeSteadyState`,
`pctExtrapolatedAumc`, `tlag`) and seeding (`mulberry32`, `deriveWorkerSeeds`).
`core/index.ts` is the authoritative list.

## What this module does NOT do

- **No DataFrame I/O.** Convert columns to `Float64Array` / `Uint8Array`
  yourself.
- **No hierarchical rule resolution.** Pass the already-resolved
  `NcaRules` for one profile.
- **No batch orchestration / workers.** `computeNca` is per-profile and
  synchronous; batch over profiles in the caller.
