# NCA — Non-Compartmental Analysis core

Pure-math, stateless computation core for Non-Compartmental Analysis.
Operates on plain `Float64Array` / `Uint8Array` inputs.

| Capability | MVP scope |
|---|---|
| **Parameters** | Cmax, Tmax, AUClast, AUCinf, %AUCextrap, λz, t½, CL (or CL/F), Vz (or Vz/F) |
| **Routes**     | IV bolus (with `c0` back-extrapolation), IV infusion, extravascular (PO/SC/IM treated uniformly) |
| **AUC methods** | linear, log-linear, linear-up/log-down — each as both naive Float64 and Neumaier-compensated summation |
| **BLQ rules**  | `set-zero`, `set-half-lloq`, `exclude`, `missing` — per-phase configurable (`preFirstMeasurable`, `embedded`, `afterLast`, `consecutiveAfterLast`) |
| **λz**         | auto best-fit (subset search; PKNCA/WinNonlin flat-tolerance tie-break — most points within `adjRSquaredFactor` of the global-max adjusted R²) + manual point selection |

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
    // minSpanRatio:   2,   // optional; DIAGNOSTIC ONLY — see below
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

**It is a diagnostic, never a gate.** Setting `lambdaZ.minSpanRatio` makes `computeNca`
emit a `LAMBDAZ_LOW_SPAN` warning (severity `warning`); it does **not** discard candidate
windows, does **not** make `lambdaZBestFit` return `null`, does **not** change which window
is selected, and does **not** affect `status`. The selected fit is byte-identical with the
threshold set and unset — `lambda-z.test.ts` and `compute-nca.test.ts` both assert that
identity, because it *is* the contract. This matches PKNCA, where `min.span.ratio` is not
read by `pk.nca()` at all and is applied post-hoc by `exclude_nca_*()`.

`minSpanRatio` deliberately defaults to `undefined` rather than to PKNCA's 2: at 2, four of
the eighteen profiles in the reference corpus would emit a warning with no consumer opt-in.

Why it earns a place next to adjusted R²: adj-R² cannot discriminate at small `n`. With
three points and one residual degree of freedom, a slope resting on 0.69 half-lives still
reports adj-R² ≈ 0.994 (`02_indometh` subject 1 — the corpus worst case, asserted in
`reference-suite.test.ts`). Span ratio is the only statistic in the PKNCA output set that
separates that fit from a well-characterised one.

The exposed namespace also gives direct access to every kernel
(`applyBlqStrategy`, `findCmax`, `lambdaZBestFit`, `lambdaZManual`,
`estimateC0`, `insertC0`, `aucLinearUpLogDownNaive`,
`aucExtrapolateToInfinity`, `halfLifeFromLambdaZ`, `clearance`,
`volumeTerminal`, `pctExtrapolated`, `mulberry32`, `deriveWorkerSeeds`)
when you need finer-grained control or want to plug them into a custom
pipeline.

## What this module does NOT do

- **No DataFrame I/O.** Convert columns to `Float64Array` / `Uint8Array`
  yourself.
- **No hierarchical rule resolution.** Pass the already-resolved
  `NcaRules` for one profile.
- **No batch orchestration / workers.** `computeNca` is per-profile and
  synchronous; batch over profiles in the caller.
