# CLAUDE.md

## Overview

`@datagrok-libraries/sci-comp` is a pure TypeScript library of numerical methods for the Datagrok platform. It has **no dependency on datagrok-api**.

## Architecture

The library is organized by numerical domain: **optimization** and **time-series**.

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
        subspace.ts               # Subspace minimisation (SMW direct primal) + Morales–Nocedal 2011 project + angle-test endpoint, 1997 truncation as fallback
    __tests__/                    # Jest tests per optimizer + registry (sync & async)
      helpers.ts                  # Test utilities: rosenbrock, sphere, gaussian, quadratic3d, etc.
      lbfgs-b.test.ts             # Scaffold + settings validation
      lbfgs-b-line-search.test.ts # dcsrch/dcstep coverage + Moré–Thuente §5 functions
      lbfgs-b-bfgs-mat.test.ts    # Compact rep: single-pair vs closed form, secant equation, K·M·z round-trip
      lbfgs-b-bounds.test.ts      # Bounds helpers unit tests
      lbfgs-b-cauchy.test.ts      # Cauchy sweep: unconstrained, 1-D snap, feasibility, c = Wᵀ(xc-x) invariant
      lbfgs-b-subspace.test.ts    # Subspace min: early exits, Newton-like property, M-N 2011 endpoint selection (angle test + 1997 truncation)
      lbfgs-b-integration.test.ts # End-to-end: Rosenbrock / Sphere / bounded / fixed / half-bounded / maximize
    examples/                     # Runnable examples (npx tsx src/optimization/single-objective/examples/*.ts)
      unconstrained.ts            # Rosenbrock, Sphere, Gaussian examples
      constrained.ts              # Box constraints example
      async-and-callbacks.ts      # Async objective functions + iteration callbacks
      gradient-descent.ts         # Gradient descent specific example
      adam.ts                     # Adam specific example
      lbfgs.ts                    # L-BFGS specific example
      registry.ts                 # Registry lookup example
    benchmarks/
      test-functions.ts           # Shared suite of classical test functions + BOUNDED_PROBLEMS + HIMMELBLAU_MINIMA
      unconstrained-benchmarks.ts # 15 standard test functions, single x₀ per problem, comparison runner
      multistart-benchmarks.ts    # Same 15 problems × 3 x₀ per problem (baseline + adversarial + near-optimum) — exposes x₀-sensitivity
      bounded-benchmarks.ts       # 7 bounded problems: L-BFGS-B native bounds vs others via penalty layer
  multi-objectives/
    moead/                        # MOEA/D multi-objective optimizer (defs.ts, moead.ts, utils.ts)
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
```

### Key design patterns

- **Optimizer base class** (`optimizer.ts`): All solvers extend `Optimizer<S>`. Subclasses implement `runInternal()`, `runInternalAsync()`, and `withDefaults()`. The base class handles input validation, constraint penalty wrapping, and the minimize/maximize inversion.
- **Sync + async API**: Each optimizer exposes `minimize`/`maximize` (sync) and `minimizeAsync`/`maximizeAsync` (for async objective functions). Penalty wrappers have sync (`applyPenalty`) and async (`applyPenaltyAsync`) variants.
- **Namespace re-exports**: The public API uses namespace re-exports (`singleObjective`, `multiObjective`, `timeSeries`) to avoid name collisions between submodules. Consumers import as `import {singleObjective, timeSeries} from '@datagrok-libraries/sci-comp'`.
- **Float64Array everywhere**: All point vectors use `Float64Array`, not `number[]`.
- **Registry pattern**: Optimizers self-register at import time via side-effect imports in `single-objective/index.ts`.
- **Iteration callbacks**: `onIteration` callback in settings allows progress monitoring and early stopping (return `true` to stop).

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
