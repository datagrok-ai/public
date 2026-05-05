# sci-comp changelog

## Unreleased

* Dunnett's test: doubled the 2D Simpson quadrature density from 80×80 to
  160×160 panels. Improves agreement with `scipy.stats.dunnett` from
  ~3e-4 worst-case to ~1e-4 — at scipy's own QMC precision floor for the
  multivariate-t CDF (Genz QMC tolerance is 1e-5; scipy's spread across
  20 RNG seeds on small-p cases is comparable to our remaining
  disagreement). Doubling further does not help.
* Boschloo's exact test (unconditional 2×2). Uses Fisher's one-sided
  p-value as the test statistic and maximises the rejection probability
  over the nuisance parameter via a grid + golden-section refinement.
  Agreement with `scipy.stats.boschloo_exact` is ~1e-11 absolute on
  preclinical-sized tables. Convenience wrapper `incidenceExactBoth`
  returns Boschloo (primary) + Fisher (alternative) together with the
  sample odds ratio, mirroring the SEND `incidence_exact_both` shape.
* `hypgeomCdf` exposed in the distributions wrapper module.
* Jonckheere-Terpstra: `continuity` default flipped from `true` to `false`
  to match `clinfun::jonckheere.test`, `PMCMRplus::jonckheereTest`, SAS
  PROC FREQ, and the SEND `jonckheere_test` Python library. Pass
  `{continuity: true}` explicitly to restore the previous Wilcoxon-style
  ±0.5 shift.
* Jonckheere-Terpstra: tightened input validation. The `k < 3` error now
  points callers to `mannWhitneyU` for two-sample comparisons. The pooled
  sample-size threshold raised from `N ≥ 3` to `N ≥ 4` (matches SEND), and
  `null` is now returned when fewer than two groups have any finite
  observations.

## 0.5.0 (2026-05-04)
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
