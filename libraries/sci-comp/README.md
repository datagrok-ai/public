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

* **Statistics** ([detailed docs](./src/stats/README.md)):
  * two-group comparison: [Welch's t-test](https://doi.org/10.1093/biomet/34.1-2.28), [Mann-Whitney U](https://doi.org/10.1214/aoms/1177730491) (exact + asymptotic), [Hedges' g](https://doi.org/10.3102/10769986006002107) effect size
  * rank correlation: [Spearman ρ](https://doi.org/10.2307/1412159), severity-trend wrapper
  * multi-group vs control: Welch pairwise + [Bonferroni](https://en.wikipedia.org/wiki/Bonferroni_correction), [Dunnett's many-to-one](https://doi.org/10.1080/01621459.1955.10501294)
  * trend tests: [Jonckheere-Terpstra](https://doi.org/10.1093/biomet/41.1-2.133) (continuous), [Cochran-Armitage](https://doi.org/10.2307/3001775) (incidence, with Buonaccorsi-modified variant)
  * sequential threshold test for proportions ([Young 1985](https://www.sra.org/))
  * categorical: [Fisher 2×2 exact](https://en.wikipedia.org/wiki/Fisher%27s_exact_test) (two-sided minlike + one-sided + odds ratio)
  * dose-response: [PAVA isotonic regression](https://en.wikipedia.org/wiki/Isotonic_regression), [Williams step-down](https://doi.org/10.2307/2528930) with 1971/1972 critical-value tables
  * covariate-adjusted analysis: ANCOVA (LS means, slope homogeneity, effect decomposition)

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

### Statistics

Run Welch's t-test on two samples with NaN handling:

```typescript
import {stats} from '@datagrok-libraries/sci-comp';

const a = [69, 70, 66, 63, 68, 70, 69, 67, 62, 63];
const b = [68, 62, 67, 68, 69, 67, 61, 59, 62, 61];

const r = stats.welchTTest(a, b);
// r.statistic ≈ 1.511, r.pValue ≈ 0.149
```

Or compare each treated group to a control with Dunnett's family-wise-error-controlled test:

```typescript
const control = [7.40, 8.50, 7.20, 8.24, 9.84, 8.32];
const treated = [
  {doseLevel: 1, values: [9.76, 8.80, 7.68, 9.36]},
  {doseLevel: 2, values: [12.80, 9.68, 12.16, 9.20, 10.55]},
];

const dun = stats.dunnettPairwise(control, treated);
// dun[1].pValueAdj ≈ 0.0058  → drug B differs from control at α = 0.05
```

See [statistics docs](./src/stats/README.md) for

* full method list (14 tests across 8 domains)
* input types (`number[]`, `Float32Array`, `Float64Array`, `Int32Array`, …)
* NaN handling (NaN as missing-value sentinel, stripped per-method)
* worked examples reproducing published references (Dunnett 1955, NIST Iris, Williams 1971/1972, Young 1985, Montgomery 15.10 vs SAS PROC GLM)
* validation against scipy via JSON fixtures (179 cases, 14 fixture files)
