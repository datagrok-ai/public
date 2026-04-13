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
    __tests__/                    # Jest tests per optimizer + registry (sync & async)
      helpers.ts                  # Test utilities: rosenbrock, sphere, gaussian, quadratic3d, etc.
    examples/                     # Runnable examples (npx tsx src/optimization/single-objective/examples/*.ts)
      unconstrained.ts            # Rosenbrock, Sphere, Gaussian examples
      constrained.ts              # Box constraints example
      async-and-callbacks.ts      # Async objective functions + iteration callbacks
      gradient-descent.ts         # Gradient descent specific example
      adam.ts                     # Adam specific example
      registry.ts                 # Registry lookup example
    benchmarks/
      unconstrained-benchmarks.ts # 15 standard test functions (Sphere, Rosenbrock, Ackley, etc.) with comparison runner
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
