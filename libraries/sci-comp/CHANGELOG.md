# sci-comp changelog

## 0.11.0 (2026-09-22 — GROK-20960)

Non-Compartmental Analysis — two bug fixes in `computeNca` that **change computed
values on affected profiles**. Every one of the 27 committed reference profiles is
byte-identical (pinned by a jest snapshot committed before either fix); the fixes
move only the regimes named below. A MINOR bump on purpose: under the `0.x` caret
rule a consumer takes these numbers only by widening its range.

* **BLQ substitutions now reach the integrator and the λz filter.** `applyBlqStrategy`
  wrote the substituted values, but `computeNca` integrated and fitted λz on the
  effective BLQ mask, so every rule was numerically `exclude` and a declared
  `set-zero` / `set-half-lloq` was a no-op. AUC/AUMC now integrate the value the rule
  wrote (`set-zero` → 0, `set-half-lloq` → LLOQ/2), and positive substitutes are
  λz-eligible (the fit's own `conc ≤ 0` guard still keeps zeros out). Cmax/Tmax and
  Tlag are observed quantities and keep reading the BLQ mask. **This changes AUC /
  AUMC / λz-derived values on every profile with an EMBEDDED or dose-time BLQ under
  `set-zero` / `set-half-lloq`, and on every profile with trailing BLQ under
  `set-half-lloq`** (fixture P1: AUClast `exclude` 15.8010 → `set-zero` 13.4729).
  Trailing `set-zero` substitutes are still trimmed (BUG-05), so trailing-BLQ-only
  profiles under `set-zero` are unchanged. `missing` ≡ `exclude` on every parameter.
  An extravascular profile whose t=0 sample is a substituted BLQ no longer gets a
  second `(0, 0)` prepended. Per-rule contract on the `BlqRule` TSDoc.
* New fixture `06_blq_rules` (5 subjects × 5 rule blocks, each a real PKNCA 0.12.1
  `conc.blq` run) — the first reference assertion for nca-studio's shipped default.
  Three cells are documented divergences from PKNCA, asserted as such (PKNCA `drop`
  leaves AUC NA without a t=0 datum; PKNCA reports λz fits below our adj-R² floor —
  its `pk.calc.half.life` never consults `min.hl.r.squared`; PKNCA extrapolates AUCinf
  from `clast.obs` while integrating `auclast` through LLOQ/2 substitutes). See
  `__tests__/REGEN.md`.
* **Route-aware IV-bolus dose-time gate.** An IV-bolus profile whose `t = 0` row is
  BLQ-flagged or non-positive / non-finite (a pre-dose sample encoded as `0`) was
  taken as a measured dose-time value and integrated from `(0, 0)`. It is now
  treated as ABSENT: the row is removed and `(0, c0)` from the PKNCA c0 chain
  takes its slot, so the profile gives the same parameters as one with no
  pre-dose row. **This changes AUClast / AUCinf / AUMC / CL / Vz / Vss / MRT on
  every IV-bolus profile with such a row** (indometh subject 1 + `(0, 0)`:
  AUClast 1.7194 → 2.0099, the no-row value). A positive measured `t = 0` value is
  used as-is; extravascular and IV-infusion profiles are unchanged. PKNCA's own
  `c0` ignores a substituted dose-time BLQ the same way (measured, fixture
  `c0_pknca`).
* `nca.augmentProfile(inputs, blq)` — pipeline Steps 1–3 (BLQ → observed Cmax →
  dose-time augmentation) as one exported, stateless kernel returning the
  augmented arrays, the effective BLQ mask, the drop set, `sourceIndex`
  (augmented ↔ input index map; `-1` = synthetic point) and the c0 estimate.
  `computeNca` now calls it; consumers that need the engine's augmented index
  space (manual λz point selection) should too — a positional mirror cannot
  track the replace case.
* `provenance.c0?: C0Provenance | null` — `{value, method, replacedDoseTimeRow,
  pctAucBackExtrap}`; `null` for non-IV-bolus routes and on `'failed'`.
  `pctAucBackExtrap` is the back-extrapolated share of AUCinf (Phoenix
  `AUC_%Back_Ext`), diagnostic only — on the indometh corpus it is 16–28 %.
  **Compile-time impact: none** — the field is OPTIONAL, so code that constructs a
  `ProfileProvenance` keeps compiling. `C0Method` moved to `types.ts` (still
  exported from the namespace).
* `nca.estimateC0Detailed` — `estimateC0` plus the chain method that answered;
  `insertC0` now returns `method` too. `estimateC0` is unchanged.
* New warning `C0_FALLBACK` (severity `warning`) when c0 fell back to `c1` /
  `cmin` / `set0` — the log-slope was not estimable. The union stays open.
* IV-bolus AUC convention documented (`__tests__/REGEN.md`): sci-comp integrates
  from the back-extrapolated c0 (WinNonlin); stock PKNCA does not (raw-profile
  `auclast` NA; 1.719365 with an observed `(0, 0)` row). Pre-existing behaviour,
  now stated.
* Fixtures: `02_indometh.json` gains provenance `c0_pknca` (PKNCA's own `c0`,
  equal to the committed `c0_extrapolated` to ≥ 10 digits on 6/6 subjects) and
  `pct_auc_back_extrap`; no previously committed value changed.

## 0.10.0 (2026-08-11)

Non-Compartmental Analysis — terminal-phase span ratio (PKNCA `span.ratio`):

* `LambdaZResult.spanRatio` — `(tEnd − tStart) / halfLife`, always reported.
* `LambdaZStrategy.minSpanRatio?` — opt-in threshold for a new `LAMBDAZ_LOW_SPAN`
  warning. Defaults to `undefined`, not PKNCA's 2.
* **Diagnostic only — no computed value changes.** Window selection, λz and every
  derived parameter are unaffected; the threshold never discards a candidate.
* Fixtures regenerated against PKNCA 0.12.1 to carry its own `span.ratio` as the
  oracle. No previously committed value changed.

## 0.9.1 (2026-07-16)

Non-Compartmental Analysis — λz best-fit window-selection fix (VAL-01-LZ-R019):

* `lambdaZBestFit` now uses PKNCA / WinNonlin's best-fit rule: among windows whose
  adjusted R² is within `adjRSquaredFactor` of the global-maximum adjusted R², keep the
  one with the most points. The prior additive score (`adjRSquared + adjRSquaredFactor·n`)
  over-selected long windows. **This changes computed λz / t½ / Vz / %AUCextrap on affected
  profiles.** Verified against PKNCA 0.12.1 (rat-IV R019: n=8→n=4, λz 0.23494→0.24047); all
  reference-suite fixtures unchanged. `factor=0` now breaks exact ties toward more points.

## 0.9.0 (2026-06-15)

Non-Compartmental Analysis — sparse / destructive-sampling NCA (UC-04 / FR-301..306):

* `nca.sparseAuc(input, options)` — design-aware closed-form composite AUClast
  with an honest standard error and degrees of freedom. Holder 2001 eq (A1)
  covariance variance + eq (A3) **unbiased** estimator, with the Nedelman-Jia
  1998 **correlated Satterthwaite df** (matrix form; reduces to the scalar
  independence df, i.e. Bailer, when no animal is sampled at two timepoints).
  Sampling topology (destructive / batch / serial) is derived from the animal ×
  nominal-time overlap matrix and cross-checked against a declared label. Honesty
  guards: linear-trapezoidal-only SE with a `SPARSE_TERMINAL_OVEREST` flag
  (Jia-Nedelman 1996), Nedelman variance-borrowing for n=1 timepoints
  (`SPARSE_VARIANCE_MODELED`), and an explicit destructive-on-absent-animal-ID
  warning. Validated against PKNCA 0.12.1 `pk.calc.sparse_auclast` to
  floating-point round-off (destructive) and a hand-derived Holder oracle (batch,
  where PKNCA returns `df = NA`). Fixture: `__tests__/fixtures/05_mouse_sparse.json`
  (generator `regen-sparse-fixture.R`).
* `nca.buildCompositeProfile(input, blqRule)` — arithmetic mean / SD / %CV / n /
  %BLQ per nominal timepoint, BLQ imputed **before** averaging, with PKNCA's
  `arithmetic mean, <=50% BLQ` rule for the AUClast endpoint.
* `nca.summarizeBootstrap(input, statistic, options)` — stratified-by-timepoint
  bootstrap with a **BCa** interval for nonlinear parameters that have no
  closed-form sparse variance. Min-n gated: an unconditional hard floor at ≤ 360
  distinct stratified resamples (Bonate 1998 eq 8 — e.g. a 5×2 destructive
  design) plus a calibrated per-timepoint minimum. Deterministic at a fixed
  master seed (`mulberry32`).

## 0.8.0 (2026-06-11)

Statistics — simple linear regression primitive (NCA dose-proportionality / UC-03):

* `stats.linearFit(x, y, {ciLevel})` — OLS fit of `y = intercept + slope·x` with
  a Student-t slope confidence interval. Returns `{slope, intercept, slopeSe,
  slopeCI, rSquared, df, n}`; default `ciLevel = 0.90` (Smith 1−2α). NaN pairs
  are dropped; degenerate inputs return (not throw): zero x-spread → NaN slope,
  `n = 2` (df 0) → slope defined but SE/CI NaN. Thin wrapper over the existing
  OLS engine + `studentTInv` — no new numerical math.
* `stats.oneWayAnova(values, groups)` — fixed-effects one-way ANOVA F-test
  (`values ~ C(groups)`), the no-covariate complement to `runAncova`. Returns
  `{fStatistic, dfBetween, dfWithin, pValue, groups, n}`; `null` on insufficient
  data (`< 2` groups, `N ≤ k`, or zero within-group variation). NaN responses
  dropped. Used for the secondary dose-normalized-AUC comparison (FR-412).
* Internal: the normal-equations `fitOls` is lifted from `tests/ancova.ts` to a
  shared `stats/internal/ols.ts` so ANCOVA and `linearFit` share one
  implementation. ANCOVA output is byte-identical (regression suite unchanged).

## 0.7.1 (2026-06-11)

Non-Compartmental Analysis — bug fix (GROK-20219 / BUG-05):

* `computeNca` now drops the **trailing run of non-positive (≤ 0) concentrations**
  from the measurable profile. A trailing `conc = 0` is an unflagged below-LLOQ
  washout sample (datasets that encode BLQ as `0` with no LLOQ/BLQ-flag column
  yield an all-zeros `blqMask`). Previously such a point anchored `cLast = 0`,
  which (a) forced `status = 'partial'` — withholding AUCinf/t½/CL/Vz — even
  though λz was perfectly well-formed, and (b) inflated `AUClast` with a spurious
  tail-to-zero trapezoid. Now matches PKNCA 0.12.1 `conc.blq` trailing-exclude
  behaviour (verified on rat-IV R005/R013). **Embedded** zeros are unchanged
  (PKNCA set-zero semantics; `lambdaZBestFit` already excludes `conc ≤ 0` from
  the regression). Profiles without a trailing zero are byte-identical.

## 0.7.0 (2026-06-09)

Non-Compartmental Analysis — FR-200 derived (moment) parameters:

* AUMClast / AUMCinf — first-moment area, own log-linear moment kernels
  (3 methods × {naive, Neumaier-compensated}) + two-term extrapolation tail
* MRT — mean residence time, unified all-route column (−T_inf/2 for IV infusion)
* Vss — steady-state volume (IV-only route gate)
* Tlag — absorption lag time (extravascular-only route gate)
* %AUMCextrap + `AUMC_EXTRAP_HIGH` warning (`NcaRules.extrapWarnPctAumc`)
* IV-infusion compute branch wires `infusionDuration` (Perrier & Mayersohn 1982)
* New IV-infusion PKNCA reference fixture; existing fixtures regenerated with
  moment/lag columns (PKNCA 0.12.1)

Breaking: `NcaRules` gains the required `extrapWarnPctAumc` field;
`ParameterValues` grows 9 → 15 fields (additive for readers).

## 0.6.0 (2026-05-08)

Non-Compartmental Analysis:

* Cmax, Tmax, AUClast, AUCinf, AUCextrap, λz, t_half, CL, Vz
* 3 AUC integration methods
* 4 BLQ-handling rules
* λz auto best-fit
* IV bolus c0 back-extrapolation
* Mulberry32 PRNG

## 0.5.1 (2026-05-05)

Update docs

## 0.5.0 (2026-05-05)

Statistics:

* Welch's t-test
* Mann-Whitney U test
* Hedges' g effect size
* Spearman rank correlation (with severity-trend helper)
* Fisher's exact test (2×2)
* Welch pairwise comparisons
* Dunnett's test
* Cochran-Armitage trend test
* Williams' test
* ANCOVA
* Jonckheere-Terpstra trend test (approximate / permutation / exact)
* Bonferroni multiple-comparison correction
* Boschloo's exact test

## 0.4.1 (2026-04-28)

Fix build

## 0.4.0 (2026-04-27)

Optimization:

* The L-BFGS-B method

## 0.3.0 (2026-04-17)

Optimization:

* The L-BFGS method

## 0.2.0 (2026-04-14)

Time series

* Calculator of 45 basic features

## 0.1.0 (2026-03-23)

Initial release

* Single-objective optimization (Nelder-Mead, PSO, GD, Adam)
* Benchmarks
* Constraint handling via penalty functions
* Optimizer registry pattern
* Multi-objective optimization (MOEA/D)
