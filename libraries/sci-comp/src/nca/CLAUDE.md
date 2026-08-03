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
    aumc.ts                         # AUMC first-moment: 3 methods × {naive, compensated} = 6 + two-term aumcExtrapolateToInfinity
    cmax.ts                         # findCmax — first-occurrence Cmax/Tmax
    lambda-z.ts                     # lambdaZBestFit (auto subset; PKNCA/WinNonlin flat-tolerance adj-R² tie-break — most points within adjRSquaredFactor of the global-max adj-R²) + lambdaZManual; centered-sum OLS
    c0.ts                           # estimateC0 + insertC0 — IV bolus back-extrapolation;
    derived.ts                      # halfLifeFromLambdaZ, clearance, volumeTerminal, pctExtrapolated, meanResidenceTime, volumeSteadyState, pctExtrapolatedAumc, tlag
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

- **Reference data lives with tests**: CSV inputs in `__tests__/datasets/` and JSON fixtures in `__tests__/fixtures/` are committed source artifacts. Regenerate via `__tests__/regen-fixtures.R` + `merge-fixtures.mjs` (PKNCA 0.12.1 oracle) — see `__tests__/REGEN.md`.

- **Route gates live in the orchestrator, not the kernels**: `vss` is `NaN` for non-IV routes (an extravascular Vss would be `Vss/F` confounded by absorption); `tlag` is `NaN` for IV routes (no absorption phase). `meanResidenceTime`/`volumeSteadyState`/`tlag` stay pure and route-agnostic; `computeNca` applies the gate and copies the `NaN` sentinel. The writer in nca-studio copies the sentinel — it never re-derives the gate.

- **AUMC has its OWN moment kernels** (`aumc.ts`), not AUC of a `t·C` array — the log-linear interval has a distinct closed form. The infinite tail is **two-term** (`(tLast·cLast)/λz + cLast/λz²`); the one-term form silently under-reports AUMC/MRT/Vss. MRT is a single unified column for all routes (`aumcInf/aucInf − T_inf/2`, `T_inf = 0` for bolus/EV).
