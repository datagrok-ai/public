# Sci Comp

Pure TypeScript library of numerical methods for the [Datagrok](https://datagrok.ai) platform.

## Modules

* **Optimization:**
  * single-objective ([detailed docs](./src/optimization/single-objective/README.md)):
    - [Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) - derivative-free simplex method
    - [PSO (Particle Swarm Optimization)](https://en.wikipedia.org/wiki/Particle_swarm_optimization) - stochastic population-based method
    - [Gradient Descent](https://en.wikipedia.org/wiki/Gradient_descent) - first-order method with numerical gradients, momentum, and learning rate decay
    - [Adam](https://arxiv.org/abs/1412.6980) - adaptive moment estimation with per-parameter learning rates
    - [L-BFGS](https://en.wikipedia.org/wiki/Limited-memory_BFGS) - limited-memory quasi-Newton method with two-loop recursion and Armijo line search
    - [L-BFGS-B](https://doi.org/10.1137/0916069) - limited-memory BFGS with native box constraints, Moré–Thuente strong-Wolfe line search, compact representation, and Morales–Nocedal 2011 subspace refinement
  * multi-objective:
    - [MOEA/D](https://ieeexplore.ieee.org/document/4358754)

* **[Time Series:](https://en.wikipedia.org/wiki/Time_series)**
  * extraction of 45 features per value column (statistics, trend, complexity, threshold-based metrics)

## Installation

```bash
npm install @datagrok-libraries/sci-comp
```

## Quick start

### Optimization

Minimize the [Rosenbrock function](https://en.wikipedia.org/wiki/Rosenbrock_function) with Nelder-Mead:

```typescript
import {singleObjective} from '@datagrok-libraries/sci-comp';

const rosenbrock = (x: Float64Array): number => {
  let sum = 0;
  for (let i = 0; i < x.length - 1; i++)
    sum += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2;
  return sum;
};

const nm = new singleObjective.NelderMead();
const result = nm.minimize(rosenbrock, new Float64Array([-1.2, 1.0]), {
  maxIterations: 5_000,
  tolerance: 1e-12,
});
// result.point ≈ [1, 1], result.value ≈ 0
```

See [single-objective optimization docs](./src/optimization/single-objective/README.md) for

* constrained optimization
* implemented methods
* sync and async objectives
* callbacks
* registry
* benchmarks
* examples

### Time-series feature extraction

Extract tsfresh-compatible features from time-series data:

```typescript
import {timeSeries} from '@datagrok-libraries/sci-comp';

const df: timeSeries.featureExtraction.TimeSeriesDataFrame = {
  ids: new Uint32Array([1, 1, 1, 1, 1]),
  time: new Float64Array([0, 1, 2, 3, 4]),
  rowCount: 5,
  columns: [{
    name: 'temperature',
    data: new Float64Array([20.1, 20.5, 21.0, 20.8, 21.3]),
  }],
};

const result = timeSeries.featureExtraction.extractFeatures(df);
// result.nSamples = 1, result.columns.length = 45
// e.g. temperature__mean, temperature__standard_deviation, temperature__linear_trend__attr_slope, ...
```

See [feature extraction docs](./src/time-series/feature-extraction/README.md) for

* full feature list (45 features across 11 categories)
* input/output format
* multiple samples and multiple columns
* validation
* naming convention
* examples
