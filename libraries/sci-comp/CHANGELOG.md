# sci-comp changelog

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
