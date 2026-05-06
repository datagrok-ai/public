# NCA core — internal documentation

Internal notes for maintainers of the NCA core. Public API consumers should
read [`../README.md`](../README.md) instead.

## Phase 1 status — complete

| File             | Purpose                                                                           | Phase 1 task |
|------------------|-----------------------------------------------------------------------------------|--------------|
| `types.ts`       | Public type contracts + `RouTE_*` constants                                       | 1.1 ✅       |
| `prng.ts`        | Mulberry32 + `deriveWorkerSeeds`                                                  | 1.2 ✅       |
| `blq.ts`         | BLQ pre-processing (4 strategies × 4 phases)                                      | 1.3 ✅       |
| `auc.ts`         | AUC: 3 methods × {naive, compensated} = 6 functions + `neumaierSum` + extrap     | 1.4 + 1.5 ✅ |
| `cmax.ts`        | `findCmax` (PKNCA first-occurrence convention)                                    | 1.6 ✅       |
| `lambda-z.ts`    | `lambdaZBestFit` (auto subset search) + `lambdaZManual` + stable centered-sum OLS | 1.7 ✅       |
| `c0.ts`          | `estimateC0` + `insertC0` — IV bolus back-extrapolation                          | 1.7 ✅       |
| `derived.ts`     | `halfLifeFromLambdaZ`, `clearance`, `volumeTerminal`, `pctExtrapolated`           | 1.8 ✅       |
| `compute-nca.ts` | `computeNca` orchestrator — full pipeline                                         | 1.9 ✅       |

Phase 2 will add `sparse.ts` (composite profile assembly) and
`bootstrap.ts` (single bootstrap iteration). Worker-pool orchestration
stays in `packages/NCA/`, not here.

## Conventions

- **Vectors are `Float64Array`**. `number[]` is not used for any
  time-series data; index lists use `Int32Array`; boolean masks use
  `Uint8Array`.
- **Pure functions, no module-level state.** Same input, same output.
  PRNG state, when needed, is passed as a parameter.
- **No throw on numerically-degenerate inputs.** Kernel functions return
  `null` (e.g. `findCmax`, `lambdaZBestFit`) or `NaN`/`±Infinity` for the
  numeric helpers; `computeNca` encodes the failure mode in
  `ComputeResult.status` (`'failed'` / `'partial'` / `'ok'`).
- **PKNCA-equivalent semantics.** Where there is a published convention,
  this module ports it 1-for-1 from PKNCA 0.12 source rather than
  reinventing it. Specific places where this matters are flagged below.

## Formulas and references

| Parameter             | Formula                                                                            | Source                       |
|-----------------------|------------------------------------------------------------------------------------|------------------------------|
| AUC linear            | `Σ (t_{i+1} − t_i) · (c_i + c_{i+1}) / 2`                                          | Gabrielsson & Weiner Ch. 2   |
| AUC log-linear        | `Σ (t_{i+1} − t_i) · (c_i − c_{i+1}) / ln(c_i / c_{i+1})`                          | Gabrielsson & Weiner Ch. 2   |
| AUC linear-up/log-down| linear if `c_{i+1} ≥ c_i`; log-linear if `c_{i+1} < c_i ∧ c_{i+1} ≠ 0`             | PKNCA `lin up/log down`      |
| AUCinf                | `AUClast + Clast / λz`                                                             | FDA guidance                 |
| %AUCextrap            | `(AUCinf − AUClast) / AUCinf · 100`                                                | FDA guidance                 |
| λz                    | `−slope` of OLS `ln(c) = α + slope · t` over best-fit terminal subset              | PKNCA auto best-fit          |
| t½                    | `ln 2 / λz`                                                                        | Standard                     |
| CL (or CL/F)          | `Dose / AUCinf`                                                                    | Standard                     |
| Vz (or Vz/F)          | `Dose / (λz · AUCinf)`                                                             | Standard                     |
| c0 (IV bolus)         | `c1 · exp(−slope·(t1 − tDose))` where `slope = (ln c2 − ln c1)/(t2 − t1)`          | PKNCA `pk.calc.c0` logslope  |

## PKNCA-equivalent design notes

These conventions matter for PKNCA reproducibility and were not obvious
from a textbook reading:

- **BLQ phasing is by first / last MEASURABLE point, not by Cmax.** Our
  `applyBlqStrategy` ignores `cmaxIdx` (kept as a documented unused
  parameter for forward compatibility with the v2.1 interface).
- **`lin up/log down` interval classification.** PKNCA tags each interval
  as `linear` / `log` / `zero`:
  - `c1 == 0 ∧ c2 == 0`           → `zero` (skipped, contributes 0).
  - `c2 < c1 ∧ c2 ≠ 0`            → `log`.
  - else                          → `linear`.
  Our `aucLinearUpLogDownNaive` matches this on every branch (linear
  formula naturally returns 0 for the (0, 0) interval).
- **λz best-fit tie-breaking.** PKNCA's `adj.r.squared.factor = 1e-4`
  effectively scores `score = adjR² + factor · n`, so larger subsets win
  on ties. Exposed as `LambdaZStrategy.adjRSquaredFactor` in TS.
- **IV bolus c0 chain** (`pk.calc.c0` default): try `c0` (pick existing
  observation at dose time), then `logslope` (back-extrapolate first 2
  decreasing post-dose points), then `c1` (use first post-dose value),
  then `cmin`, then `set0`. Our `estimateC0` ports the chain literally,
  including the `c2 < c1 ∧ c2 ≠ 0` guard in `logslope`.
- **t=0 augmentation in the orchestrator.** When the input has no t=0
  observation, `computeNca` prepends one before AUC and λz: `(0, c0)` for
  IV bolus, `(0, 0)` for extravascular. This is what PKNCA's R pipeline
  does internally; without it AUClast under-estimates the back-extrapolated
  area and λz windowing on IV bolus diverges from the reference fixtures.
- **Observed Cmax is preserved** even when the orchestrator inserts a
  point. The reported `cmax` / `tmax` come from the original profile, so
  for IV bolus they are NOT the inserted `c0`.

## Validation

End-to-end validation lives in
[`__tests__/reference-suite.test.ts`](__tests__/reference-suite.test.ts):
26 profiles across 3 datasets (Theophylline 12, Indomethacin 6, synthetic
rat 8), every Phase-1 parameter compared against PKNCA 0.12.1 fixtures
within the §9.2 tolerances of the development plan.

The `lambda-z-validation.test.ts` keeps a focused gate on `lambdaZBestFit`
(plus `insertC0` for IV bolus) — useful for fast diagnosis when only the
λz fit changes.

## Sources

- PKNCA documentation — https://billdenney.github.io/pknca/
- PKNCA source (browsed via `body(...)` in R 4.6 + PKNCA 0.12.1) for the
  exact algorithms of `pk.calc.c0`, `choose_interval_method`,
  `aucintegrate_linear`, `aucintegrate_log`.
- FDA Guidance for Industry: Bioavailability and Bioequivalence Studies.
- EMA Guideline on the Investigation of Bioequivalence.
- Gabrielsson J., Weiner D. — *Pharmacokinetic and Pharmacodynamic Data
  Analysis*.
- Neumaier, A. (1974). *Rundungsfehleranalyse einiger Verfahren zur
  Summation endlicher Summen.* ZAMM 54: 39–51 (compensated AUC).
