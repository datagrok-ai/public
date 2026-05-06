# CLAUDE.md

## Overview

`@datagrok-libraries/sci-comp` is a pure TypeScript library of numerical methods for the Datagrok platform. It has **no dependency on datagrok-api**.

## Architecture

The library is organized by numerical domain: **optimization**, **stats**, and **time-series**.

```
index.ts                          # Entry point: re-exports {singleObjective, multiObjective, timeSeries} namespaces
src/optimization/
  single-objective/
    types.ts                      # ObjectiveFunction, AsyncObjectiveFunction, OptimizationResult, Constraint, CommonSettings
    optimizer.ts                  # Abstract Optimizer<S> base class (validation, penalty wiring, minimize/maximize + async variants)
    penalty.ts                    # applyPenalty(), applyPenaltyAsync(), boxConstraints() — constraint → penalized objective
    registry.ts                   # registerOptimizer/getOptimizer/listOptimizers — name-based lookup
    optimizers/
      nelder-mead.ts              # NelderMead extends Optimizer<NelderMeadSettings>
      pso.ts                      # PSO extends Optimizer<PSOSettings>
      gradient-descent.ts         # GradientDescent extends Optimizer<GradientDescentSettings>
      adam.ts                     # Adam extends Optimizer<AdamSettings>
      lbfgs.ts                    # LBFGS extends Optimizer<LBFGSSettings> (quasi-Newton with two-loop recursion + Armijo line search)
      lbfgs-b/                    # L-BFGS-B (limited-memory BFGS with native box constraints)
        index.ts                  # LBFGSB class (public surface) + withDefaults validation
        types.ts                  # LBFGSBSettings, LBFGSBLineSearchSettings, LBFGSBBounds, BOUND_* constants
        driver.ts                 # runSync / runAsync outer loop (Cauchy → subspace → line search → memory update)
        line-search.ts            # dcsrch + dcstep — Moré–Thuente strong-Wolfe line search + NaN-bisection wrapper
        bfgs-mat.ts               # BFGSMat compact representation (ring buffer, W-products, block Cholesky, solveM)
        bounds.ts                 # normalizeBounds, classifyBounds, project, projectedGradient, maxFeasibleStep
        cauchy.ts                 # Generalized Cauchy point + binary min-heap
        subspace.ts               # Subspace minimisation (SMW direct primal) + Morales–Nocedal 2011 endpoint selection (directional-derivative test) + xc-anchored 1997 truncation with bound snap
    __tests__/                    # Jest tests per optimizer + registry (sync & async)
      helpers.ts                  # Test utilities: rosenbrock, sphere, gaussian, quadratic3d, etc.
      lbfgs-b.test.ts             # Scaffold + settings validation
      lbfgs-b-line-search.test.ts # dcsrch/dcstep coverage + Moré–Thuente §5 functions
      lbfgs-b-bfgs-mat.test.ts    # Compact rep: single-pair vs closed form, secant equation, K·M·z round-trip
      lbfgs-b-bounds.test.ts      # Bounds helpers unit tests
      lbfgs-b-cauchy.test.ts      # Cauchy sweep: unconstrained, 1-D snap, feasibility, c = Wᵀ(xc-x) invariant
      lbfgs-b-subspace.test.ts    # Subspace min: early exits, Newton-like property, M-N 2011 endpoint selection (directional-derivative test, xc-anchored truncation, bound snap)
      lbfgs-b-integration.test.ts # End-to-end: Rosenbrock / Sphere / bounded / fixed / half-bounded / maximize
    examples/                     # Runnable examples (npx tsx src/optimization/single-objective/examples/*.ts)
      unconstrained.ts            # Rosenbrock, Sphere, Gaussian examples
      constrained.ts              # Box constraints example
      async-and-callbacks.ts      # Async objective functions + iteration callbacks
      gradient-descent.ts         # Gradient descent specific example
      adam.ts                     # Adam specific example
      lbfgs.ts                    # L-BFGS specific example
      lbfgs-b.ts                  # L-BFGS-B specific example (native box constraints)
      khai-unconstraint.ts        # NelderMead + PSO on a 3D quadratic + product-surface problem
      khai-constraint.ts          # NelderMead + PSO with equality constraint via applyPenalty
      registry.ts                 # Registry lookup example
    benchmarks/
      test-functions.ts           # Shared suite of classical test functions + BOUNDED_PROBLEMS + HIMMELBLAU_MINIMA
      unconstrained-benchmarks.ts # 15 standard test functions, single x₀ per problem, comparison runner
      multistart-benchmarks.ts    # Same 15 problems × 3 x₀ per problem (baseline + adversarial + near-optimum) — exposes x₀-sensitivity
      bounded-benchmarks.ts       # 7 bounded problems: L-BFGS-B native bounds vs others via penalty layer
  multi-objectives/
    moead/                        # MOEA/D multi-objective optimizer (defs.ts, moead.ts, utils.ts)
src/stats/
  index.ts                        # Public entry — re-exports tests, types, distributions, multiple-comparison helpers
  README.md                       # Per-domain method tables, code snippets, run instructions for all 14 stats methods
  types.ts                        # NumericInput, Alternative, TestResult, SpearmanResult, FisherResult
  distributions.ts                # Typed wrappers around jstat (normal, t, F, chi², hypergeom, special functions)
  internal/
    normalize.ts                  # NumericInput → Float64Array, NaN stripping, mean / variance / std / sum
    rank.ts                       # Rank assignment with average-rank tie handling
    matrix.ts                     # Linear algebra helpers (matrix inverse via jstat)
    random.ts                     # mulberry32 seedable PRNG + shuffleInPlace (Fisher–Yates)
  tests/                          # One file per statistical test
    welch-t.ts, mann-whitney.ts, hedges-g.ts, spearman.ts, fisher-exact.ts
    welch-pairwise.ts, dunnett.ts, cochran-armitage.ts, ancova.ts
    williams.ts, williams-tables.ts
    jonckheere.ts                 # Full JT: approximate / permutation / exact, ±continuity, tie-corrected variance
    boschloo-exact.ts             # Unconditional exact test: Fisher one-sided p as test stat, sup over nuisance π via grid + golden-section refinement; `incidenceExactBoth` returns Boschloo + Fisher together
  multiple-comparison/
    bonferroni.ts                 # bonferroniCorrect — multiplicity adjustment
  __tests__/                      # Jest tests, one per method, all driven by JSON fixtures in __tests__/fixtures/
    helpers.ts                    # loadFixture, expectClose (note: uses `diff >= tol`, so pick tol > 0)
    fixtures/                     # 18 JSON fixtures (179 cases) generated by scripts/generate-fixtures.py from scipy/numpy
      ancova.json, bonferroni.json, cochran-armitage.json, cochran-armitage-modified.json
      dunnett.json, fisher-exact.json, hedges-g.json, mann-whitney.json
      severity-trend.json, spearman.json, welch-pairwise.json, welch-t.json, williams.json
      jonckheere-rp-approximate.json  # 24 cases vs Python regressionpack (tight)
      jonckheere-rp-exact.json    # 5 cases vs Python regressionpack exact method
      jonckheere-clinfun.json     # 144 cases vs R clinfun::jonckheere.test
      jonckheere-pmcmr.json       # 144 cases vs R PMCMRplus permutation
      boschloo-exact.json         # 31 cases vs scipy.stats.boschloo_exact (hand-picked + 60 randomised)
  examples/                       # 9 runnable npx tsx examples + _helpers.ts (ancova, boschloo-exact, compare-two-groups, correlation, fisher-exact, pairwise-vs-control, threshold-test, trend-tests, williams-step-down)
src/time-series/
  feature-extraction/
    types.ts                      # NumericArray, TimeSeriesDataFrame, TimeSeriesColumn, FeatureColumn, FeatureMatrix, ExtractOptions
    extract.ts                    # extractFeatures() — main entry point, orchestrates passes + median + linear trend
    calculators.ts                # Multi-pass feature calculators: pass1 (basic stats), pass2 (diffs, uniqueness, strikes), pass3 (threshold-based), computeMedian
    linear-trend.ts               # linearTrend() — slope, intercept, rvalue, pvalue, stderr via least-squares
    dataframe.ts                  # buildIndex(), validate() — sample grouping and input validation
    stats-utils.ts                # Statistical utility functions
    __tests__/                    # Jest tests for calculators, extract, linear-trend, naming, types, validate
      helpers.ts                  # Test utilities for time-series tests
    examples/                     # Runnable examples
      single-sample.ts            # Single sample feature extraction
      multiple-samples.ts         # Multi-sample feature extraction
      multiple-columns.ts         # Multi-column feature extraction
src/nca/
  index.ts                        # Public namespace entry — re-exports core
  README.md                       # Module overview + quickstart
  core/
    README.md                     # Internal documentation: formulas, sources, pipeline
    types.ts                      # ProfileInputs, NcaRules, ComputeResult, BlqStrategy, LambdaZStrategy, ParameterValues, ProfileProvenance, RouteCode constants
    prng.ts                       # mulberry32 + deriveWorkerSeeds (PKNCA-Variant-A sub-stream derivation for bootstrap)
    blq.ts                        # applyBlqStrategy — 4 BLQ rules × 4 phases (PKNCA convention)
    auc.ts                        # AUC: 3 methods × {naive, Neumaier-compensated} = 6 functions + neumaierSum + aucExtrapolateToInfinity
    cmax.ts                       # findCmax — first-occurrence Cmax/Tmax (PKNCA convention)
    lambda-z.ts                   # lambdaZBestFit (auto subset search by adj-R² + PKNCA tie-breaking factor) + lambdaZManual; numerically stable centered-sum OLS
    c0.ts                         # estimateC0 + insertC0 — IV bolus back-extrapolation; faithful port of PKNCA pk.calc.c0 (chain: c0 → logslope → c1 → cmin → set0)
    derived.ts                    # halfLifeFromLambdaZ, clearance, volumeTerminal, pctExtrapolated
    compute-nca.ts                # computeNca orchestrator — full pipeline: BLQ → t=0 augmentation (IV bolus c0 / extravascular zero) → observed Cmax → AUClast → lambda_z → AUCinf + derived → quality warnings
    __tests__/                    # Jest tests, one per module + reference-suite + lambda-z-validation
      prng.test.ts, blq.test.ts, auc.test.ts, cmax.test.ts
      lambda-z.test.ts, lambda-z-validation.test.ts
      c0.test.ts, derived.test.ts, compute-nca.test.ts
      reference-suite.test.ts     # Full computeNca pipeline vs PKNCA fixtures (26 profiles, all 8 parameters within §9.2 tolerances)
  __tests__/                      # Cross-module assets
    datasets/                     # CSV input datasets (committed source artifacts)
      01_theoph.csv, 02_indometh.csv, 03_rat_simple.csv, README.md
    fixtures/                     # PKNCA reference values (generated by scripts/generate-nca-fixtures.R)
      01_theoph.json, 02_indometh.json, 03_rat_simple.json
```

### Key design patterns

- **Optimizer base class** (`optimizer.ts`): All solvers extend `Optimizer<S>`. Subclasses implement `runInternal()`, `runInternalAsync()`, and `withDefaults()`. The base class handles input validation, constraint penalty wrapping, and the minimize/maximize inversion.
- **Sync + async API**: Each optimizer exposes `minimize`/`maximize` (sync) and `minimizeAsync`/`maximizeAsync` (for async objective functions). Penalty wrappers have sync (`applyPenalty`) and async (`applyPenaltyAsync`) variants.
- **Namespace re-exports**: The public API uses namespace re-exports (`singleObjective`, `multiObjective`, `timeSeries`) to avoid name collisions between submodules. Consumers import as `import {singleObjective, timeSeries} from '@datagrok-libraries/sci-comp'`.
- **Float64Array everywhere**: All point vectors use `Float64Array`, not `number[]`.
- **Registry pattern**: Optimizers self-register at import time via side-effect imports in `single-objective/index.ts`.
- **Iteration callbacks**: `onIteration` callback in settings allows progress monitoring and early stopping (return `true` to stop).

### NCA (Non-Compartmental Analysis)

- **PKNCA-equivalent semantics**: Where the algorithm has a published convention, this module ports it 1-for-1 from PKNCA 0.12 source — including the BLQ phase-detection rules, the `lin up/log down` AUC method, the lambda_z best-fit subset search with `adj.r.squared.factor` tie-breaking, and the `c0` back-extrapolation chain (`c0 → logslope → c1 → cmin → set0`). Validated against PKNCA on 26 profiles to within §9.2 tolerances of `docs/nca_development_plan_v2.md`.
- **One orchestrator, isolated kernels**: `computeNca` is the only entry point that fuses the steps. Each kernel (`applyBlqStrategy`, `findCmax`, `aucLinearUpLogDownNaive`, `lambdaZBestFit`, `estimateC0`, `halfLifeFromLambdaZ`, …) is a pure function tested in isolation and re-exported from the namespace for direct use.
- **t=0 augmentation in the orchestrator**: For IV bolus profiles without a t=0 sample the orchestrator inserts `(0, c0)` (estimated by `insertC0`); for extravascular profiles it inserts `(0, 0)` by convention. The kernels themselves don't know about this — keeps them route-agnostic.
- **Observed vs. computed Cmax**: The reported Cmax/Tmax are the OBSERVED peak from the original (non-augmented) profile, even when the orchestrator inserts a t=0 point. Internal lambda_z fit and AUC integration use the augmented profile.
- **Status flag separates degeneracy modes**: `'failed'` (no measurable point), `'partial'` (Cmax/AUClast computed but lambda_z not estimable → no AUCinf, t½, CL, Vz), `'ok'` (all parameters). All numeric fields default to `NaN` when not computed.
- **No DataGrok dependency**: Same invariant as the rest of sci-comp. The DataFrame adapter, validation, controller, and worker pool live in the separate `packages/NCA/` package (planned, not in this library).
- **Reference data lives with tests**: CSV inputs in `src/nca/__tests__/datasets/`, PKNCA JSON fixtures in `src/nca/__tests__/fixtures/`. Generators in `scripts/`: `generate-dataset-03-rat-simple.R` (synthetic rat data) and `generate-nca-fixtures.R` (PKNCA → JSON).

### Time-series feature extraction

- **tsfresh-compatible**: Produces ~45 features per value column following tsfresh naming convention (`{col}__{feature}` or `{col}__{feature}__{param}_{value}`).
- **Multi-pass architecture** (`calculators.ts`): Computations are split into three passes to minimize iterations over each sample slice — pass1 (basic stats), pass2 (diffs, uniqueness, crossings, strikes), pass3 (threshold-based counts), plus separate median and linear trend.
- **Columnar I/O**: Input is a `TimeSeriesDataFrame` with `ids` (sample grouping), `time`, and value `columns`. Output is a `FeatureMatrix` with one row per sample and one `FeatureColumn` per feature.
- **NumericArray**: Accepts `Int32Array`, `Uint32Array`, `Float32Array`, or `Float64Array` as input data types.
- **Contiguous id grouping**: All rows for the same sample id must be contiguous in the input (sorted by id).

### Benchmarks

The single-objective benchmark suite is split into two complementary runners that
share objective functions through `benchmarks/test-functions.ts`:

- `unconstrained-benchmarks.ts` — one x₀ per problem, direct head-to-head table.
- `multistart-benchmarks.ts` — three x₀ per problem (baseline + adversarial
  perturbation + near-optimum) and a success-rate summary across all 45 runs.
- `bounded-benchmarks.ts` — 7 box-constrained problems. **L-BFGS-B** uses its
  native `settings.bounds`; all other optimizers route box constraints through
  `boxConstraints()` + the quadratic-penalty layer. Success requires both
  numerical accuracy AND exact feasibility (`feas_vio < 1e-6`) — the latter is
  where penalty-based methods typically fail (~1e-3 violations).

When adding a new optimizer, register it in **all three** runners and
regenerate the `.md` reports. When adding a new test function, export it from
`test-functions.ts` (do not duplicate the body in the runner files); bounded
problems are recorded in `BOUNDED_PROBLEMS` there.

**x₀-sensitivity caveat — important when interpreting results.** The single-start
tables can make local optimizers look stronger than they are on multimodal
problems. Example: L-BFGS solves Rastrigin / Lévi N.13 in one iteration because
the baseline x₀ is integer-aligned and zeroes out the `sin(kπxᵢ)` gradient terms;
a 0.1 perturbation of x₀ destroys that effect (visible in the multi-start tables).
Treat impressive single-start results on multimodal objectives as hypotheses to
verify against `multistart-benchmarks.md`.

### Test structure

Optimization tests are grouped by problem, each containing `sync` and `async` variants:

```
describe('minimize Rosenbrock 2D → min ≈ 0 at (1, 1)', () => {
  it('sync', () => { ... });
  it('async', async () => { ... });
});
```

Time-series tests are organized by module (calculators, extract, linear-trend, naming, types, validate).

## Skills

- `/add-optimizer <algorithm name>` — scaffold a new single-objective optimizer (class, tests, registration, exports)
