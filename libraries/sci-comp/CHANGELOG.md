# sci-comp changelog

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
