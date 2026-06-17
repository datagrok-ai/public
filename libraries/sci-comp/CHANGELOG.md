# sci-comp changelog

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
