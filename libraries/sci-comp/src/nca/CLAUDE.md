# CLAUDE.md — nca

Non-Compartmental Analysis (PK).

## Architecture

```
src/nca/
  index.ts                          # Public namespace entry — re-exports core
  README.md                         # Module overview + quickstart
  core/    
    types.ts                        # ProfileInputs, NcaRules, ComputeResult, BlqStrategy, LambdaZStrategy, ParameterValues, ProfileProvenance, RouteCode
    prng.ts                         # mulberry32 + deriveWorkerSeeds
    blq.ts                          # applyBlqStrategy — 4 BLQ rules × 4 phases
    auc.ts                          # AUC: 3 methods × {naive, Neumaier-compensated} = 6 + neumaierSum + aucExtrapolateToInfinity
    cmax.ts                         # findCmax — first-occurrence Cmax/Tmax
    lambda-z.ts                     # lambdaZBestFit (auto subset, adj-R² + tie-break) + lambdaZManual; centered-sum OLS
    c0.ts                           # estimateC0 + insertC0 — IV bolus back-extrapolation;
    derived.ts                      # halfLifeFromLambdaZ, clearance, volumeTerminal, pctExtrapolated
    compute-nca.ts                  # computeNca orchestrator — full pipeline
    __tests__/                      # Per-module tests + reference-suite vs fixtures
  __tests__/                        # Cross-module assets
    datasets/                       # CSV inputs (committed)
    fixtures/                       # Reference values (committed JSON)
```

## Key design patterns

- **One orchestrator, isolated kernels**: `computeNca` is the only entry point that fuses the steps. Each kernel (`applyBlqStrategy`, `findCmax`, `aucLinearUpLogDownNaive`, `lambdaZBestFit`, `estimateC0`, `halfLifeFromLambdaZ`, …) is a pure function tested in isolation and re-exported from the namespace for direct use. Keep new logic in kernels; only add to the orchestrator when it genuinely fuses kernel results.

- **t=0 augmentation lives in the orchestrator only**: for IV bolus profiles without a t=0 sample the orchestrator inserts `(0, c0)` (estimated by `insertC0`); for extravascular profiles it inserts `(0, 0)` by convention. Kernels stay route-agnostic and never know about augmentation.

- **Observed vs. computed Cmax**: reported Cmax/Tmax are the OBSERVED peak from the original (non-augmented) profile, even when the orchestrator inserts a t=0 point. Internal lambda_z fit and AUC integration use the AUGMENTED profile. Don't conflate.

- **Status flag separates degeneracy modes**: `'failed'` (no measurable point), `'partial'` (Cmax/AUClast computed but lambda_z not estimable → no AUCinf, t½, CL, Vz), `'ok'` (all parameters). All numeric fields default to `NaN` when not computed.

- **Reference data lives with tests**: CSV inputs in `__tests__/datasets/` and JSON fixtures in `__tests__/fixtures/` are committed source artifacts.
