# CLAUDE.md — nca

Non-Compartmental Analysis (PK).

## Architecture

```
src/nca/
  index.ts                          # Public namespace entry — re-exports core
  README.md                         # Module overview + quickstart
  core/
    index.ts                        # Authoritative export list for the namespace
    types.ts                        # ProfileInputs, NcaRules, ComputeResult, BlqStrategy, LambdaZStrategy, ParameterValues, ProfileProvenance, RouteCode
    prng.ts                         # mulberry32 + deriveWorkerSeeds
    blq.ts                          # applyBlqStrategy — 4 BLQ rules × 4 phases
    auc.ts                          # AUC: 3 methods × {naive, Neumaier-compensated} = 6 + neumaierSum + aucExtrapolateToInfinity
    aumc.ts                         # AUMC first-moment: 3 methods × {naive, compensated} = 6 + two-term aumcExtrapolateToInfinity
    cmax.ts                         # findCmax — first-occurrence Cmax/Tmax
    lambda-z.ts                     # lambdaZBestFit (auto subset; PKNCA/WinNonlin flat-tolerance adj-R² tie-break — most points within adjRSquaredFactor of the global-max adj-R²) + lambdaZManual; centered-sum OLS; reports spanRatio (diagnostic only, never gates)
    c0.ts                           # estimateC0 / estimateC0Detailed + insertC0 — IV bolus back-extrapolation (PKNCA c0/logslope/c1/cmin/set0 method chain; the detailed form reports WHICH method answered)
    augment.ts                      # augmentProfile — pipeline Steps 1–3 as ONE exported kernel: BLQ → observed Cmax → route-aware dose-time augmentation; returns the effective BLQ mask, the drop set, sourceIndex (augmented ↔ input), and the c0 estimate
    derived.ts                      # halfLifeFromLambdaZ, clearance, volumeTerminal, pctExtrapolated, meanResidenceTime, volumeSteadyState, pctExtrapolatedAumc, tlag
    compute-nca.ts                  # computeNca orchestrator — augmentProfile + Steps 4–7 (integration, λz, derived, warnings)
    sparse.ts                       # sparseAuc + buildCompositeProfile — composite AUClast for destructive/batch designs; Holder covariance SE, Nedelman-Jia Satterthwaite df, Student-t CI
    bootstrap.ts                    # summarizeBootstrap — stratified-by-timepoint resampling + BCa interval for the nonlinear parameters, with a self-suppression gate
    __tests__/                      # Per-module tests + reference-suite vs fixtures
  __tests__/                        # Cross-module assets
    datasets/                       # CSV inputs (committed)
    fixtures/                       # Reference values (committed JSON)
```

## Key design patterns

- **One orchestrator, isolated kernels**: `computeNca` is the only entry point that fuses the steps. Each kernel (`applyBlqStrategy`, `findCmax`, `aucLinearUpLogDownNaive`, `lambdaZBestFit`, `estimateC0`, `halfLifeFromLambdaZ`, …) is a pure function tested in isolation and re-exported from the namespace for direct use. Keep new logic in kernels; only add to the orchestrator when it genuinely fuses kernel results.

- **t=0 augmentation lives in ONE exported kernel (`augmentProfile`)**, consumed by `computeNca` and by callers that need the engine's augmented index space (nca-studio's manual-λz editor maps `pointsUsed` / `manualPoints` through `sourceIndex`). The kernel is route-aware: for IV bolus a positive measured t=0 sample is the c0 (`method: 'observed'`); a MISSING t=0 row is prepended as `(0, c0)`, and a PRESENT-but-BLQ / non-positive t=0 row is REMOVED and replaced by `(0, c0)` (`replacedDoseTimeRow`) — a pre-dose sample can never anchor the back-extrapolation, whatever the substitution rule (PKNCA's own `c0` agrees, measured). For extravascular / IV infusion a kept t=0 row (even `conc = 0`) counts; otherwise `(0, 0)` is prepended by convention. Only index 0 is ever the dose-time row. The numeric kernels (`applyBlqStrategy`, `findCmax`, `insertC0`, `auc*`, `lambdaZ*`) stay route-agnostic and never know about augmentation. Do NOT re-orchestrate Steps 1–3 anywhere else — a positional mirror cannot track REPLACE (index 0 becomes synthetic and input row 0 disappears), which is exactly why `sourceIndex` is reported.

- **Two masks, on purpose**: `AugmentedProfile.blqMask` is the EFFECTIVE BLQ mask (input mask ∪ `exclude`d) that the OBSERVED quantities read — Cmax/Tmax skip it, Tlag treats it as 0, `BLQ_HIGH_FRACTION` counts it. `AugmentedProfile.dropMask` is what integration and the λz regression SKIP. They are separate because a substituted BLQ value (`set-zero`, `set-half-lloq`) must reach the integrator as the value the rule wrote while Cmax keeps skipping it — one mask cannot express both. Never collapse them.

- **IV-bolus AUC integrates FROM the back-extrapolated c0** (Phoenix WinNonlin convention), so the dose-time → first-sample area is part of AUClast whether or not a pre-dose sample was drawn. Stock PKNCA does not do this (its raw-profile `auclast` is NA without a t=0 datum and integrates from an observed `(0, 0)` when present — measured 2026-09-22); the reference fixtures feed PKNCA the augmented profile. A deliberate, documented divergence — see `__tests__/REGEN.md`. `provenance.c0.pctAucBackExtrap` reports how much of AUCinf that segment is (Phoenix `AUC_%Back_Ext`), diagnostic only, never a gate.

- **Observed vs. computed Cmax**: reported Cmax/Tmax are the OBSERVED peak from the original (non-augmented) profile, even when the kernel inserts a t=0 point. Internal lambda_z fit and AUC integration use the AUGMENTED profile. Don't conflate.

- **Status flag separates degeneracy modes**: `'failed'` (no measurable point), `'partial'` (Cmax/AUClast computed but lambda_z not estimable → no AUCinf, t½, CL, Vz), `'ok'` (all parameters). All numeric fields default to `NaN` when not computed.

- **Span ratio is a diagnostic, NOT a gate**: `LambdaZResult.spanRatio` = `(tEnd − tStart)/halfLife` is always reported; `LambdaZStrategy.minSpanRatio` only decides whether `computeNca` emits a `LAMBDAZ_LOW_SPAN` warning. It must never discard a candidate window, never make `lambdaZBestFit` return `null`, never touch `status`. The selected fit is identical with the threshold set and unset — asserted in both `lambda-z.test.ts` and `compute-nca.test.ts`, and that identity IS the contract. Defaults to `undefined`, deliberately not to PKNCA's 2: gating at 2 would flip 4 of the 27 reference profiles off the lambda_z PKNCA returns, destroying the parity the fixtures exist to protect. Full rationale is on the `LambdaZResult.spanRatio` TSDoc — don't restate it in a third place.

- **Reference data lives with tests**: CSV inputs in `__tests__/datasets/` and JSON fixtures in `__tests__/fixtures/` are committed source artifacts. Regenerate via `__tests__/regen-fixtures.R` + `merge-fixtures.mjs` (PKNCA 0.12.1 oracle) — see `__tests__/REGEN.md`. `06_blq_rules.json` is the per-rule BLQ oracle (five `BlqStrategy` blocks × five subjects, each block a real PKNCA `conc.blq` run); three of its cells are DOCUMENTED divergences from PKNCA (PKNCA `drop` leaves AUC NA without a t=0 datum; PKNCA reports λz fits below our adj-R² floor — `pk.calc.half.life` never consults `min.hl.r.squared`; PKNCA extrapolates AUCinf from `clast.obs` while integrating `auclast` through LLOQ/2 substitutes) and the suite asserts them as such — never widen a tolerance or skip a cell to make a divergence disappear.

- **The 27-profile snapshot is a permanent byte-identity gate** (`reference-suite.test.ts`, `__snapshots__/`): every `ParameterValues` field + λz window + status on 01/02/03/04 in both summation modes. The tolerance assertions prove PKNCA parity; the snapshot proves "unchanged". A deliberate numbers-changing fix updates it (`jest -u`) in its own commit with the CHANGELOG naming the moved profiles. An unexplained snapshot diff is a defect.

- **Route gates live in the orchestrator, not the kernels**: `vss` is `NaN` for non-IV routes (an extravascular Vss would be `Vss/F` confounded by absorption); `tlag` is `NaN` for IV routes (no absorption phase). `meanResidenceTime`/`volumeSteadyState`/`tlag` stay pure and route-agnostic; `computeNca` applies the gate and copies the `NaN` sentinel. The writer in nca-studio copies the sentinel — it never re-derives the gate.

- **AUMC has its OWN moment kernels** (`aumc.ts`), not AUC of a `t·C` array — the log-linear interval has a distinct closed form. The infinite tail is **two-term** (`(tLast·cLast)/λz + cLast/λz²`); the one-term form silently under-reports AUMC/MRT/Vss. MRT is a single unified column for all routes (`aumcInf/aucInf − T_inf/2`, `T_inf = 0` for bolus/EV).

- **Sparse is a separate entry point, not a mode of `computeNca`**: destructive/batch designs make per-subject NCA undefined, so `sparseAuc` takes a different input shape (`SparseInput` — flat columnar over animal × nominal time) and returns a different result. There is no fallback path between the two, and `computeNca` must not grow one.

- **Sparse topology comes from the data, never from the label**: the r_ij overlap matrix decides destructive / batch / serial; a caller-supplied `declaredTopology` is only cross-checked (`SPARSE_TOPOLOGY_MISMATCH`). One code path covers all three because Holder's covariance estimator reduces to Bailer's when every r_ij = 0 — don't branch on topology.

- **`blq.ts` is deliberately NOT reused by `sparse.ts`**: the phase model (preFirstMeasurable / embedded / afterLast) is a per-profile time-ordering concept, and the composite pools many animals at ONE nominal timepoint where no phase exists. Only the `BlqRule` type and its per-rule semantics are shared, and BLQ is imputed **before** averaging. Routing sparse through `applyBlqStrategy` looks like de-duplication but forces an ill-fitting model.

- **The bootstrap suppresses itself instead of returning a fake interval**: `summarizeBootstrap` returns `suppressed: true` + `suppressReason` when the design is combinatorially degenerate (`h = Π_i C(2n_i − 1, n_i) ≤ 360`, an unconditional floor) or any timepoint falls below the calibrated `minNPerTimepoint`. Callers fall back to the closed-form CI, which is **never** gated. Don't add a force/override flag — the suppression IS the result.
