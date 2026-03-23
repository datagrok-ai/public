# Sci Comp

Pure TypeScript library of numerical methods for the [Datagrok](https://datagrok.ai) platform.

* Optimization:
  * single-objective:
    - [Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) - derivative-free simplex method
    - [PSO (Particle Swarm Optimization)](https://en.wikipedia.org/wiki/Particle_swarm_optimization) - stochastic population-based method
    - [Gradient Descent](https://en.wikipedia.org/wiki/Gradient_descent) - first-order method with numerical gradients, momentum, and learning rate decay
    - [Adam](https://arxiv.org/abs/1412.6980) - adaptive moment estimation with per-parameter learning rates
  * multi-objective:
    - [MOEA/D](https://ieeexplore.ieee.org/document/4358754)

## Installation

```bash
npm install @datagrok-libraries/sci-comp
```

## Optimization

```typescript
import {singleObjective} from '@datagrok-libraries/sci-comp';

const {NelderMead, PSO, applyPenalty, applyPenaltyAsync, boxConstraints, getOptimizer, listOptimizers} = singleObjective;
```

### Single-objective

Four built-in solvers, each supporting synchronous and asynchronous objective functions:

- [Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) - derivative-free simplex method
- [PSO (Particle Swarm Optimization)](https://en.wikipedia.org/wiki/Particle_swarm_optimization) - stochastic population-based method
- [Gradient Descent](https://en.wikipedia.org/wiki/Gradient_descent) - first-order method with numerical gradients, momentum, and learning rate decay
- [Adam](https://arxiv.org/abs/1412.6980) - adaptive moment estimation with per-parameter learning rates

| Method | Sync | Async |
|--------|------|-------|
| Minimize | `minimize(fn, x0, settings)` | `minimizeAsync(fn, x0, settings)` |
| Maximize | `maximize(fn, x0, settings)` | `maximizeAsync(fn, x0, settings)` |

#### Unconstrained minimize

Minimize the [Rosenbrock function](https://en.wikipedia.org/wiki/Rosenbrock_function):

$$f(x, y) = 100\,(y - x^2)^2 + (1 - x)^2$$

- **Goal:** minimize
- **Start:** $x_0 = (-1.2,\; 1.0)$
- **Expected:** $\min = 0$ at $(1, 1)$

```typescript
const rosenbrock = (x: Float64Array): number => {
  let sum = 0;
  for (let i = 0; i < x.length - 1; i++)
    sum += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2;
  return sum;
};

const nm = new NelderMead();
const result = nm.minimize(rosenbrock, new Float64Array([-1.2, 1.0]), {
  maxIterations: 5_000,
  tolerance: 1e-12,
});
// result.point ≈ [1, 1], result.value ≈ 0
```

#### Unconstrained maximize

Maximize the Gaussian function:

$$f(x, y) = e^{-(x^2 + y^2)}$$

- **Goal:** maximize
- **Start:** $x_0 = (2,\; -3)$
- **Expected:** $\max = 1$ at $(0, 0)$

```typescript
const gaussian = (x: Float64Array): number =>
  Math.exp(-(x[0] ** 2 + x[1] ** 2));

const result = nm.maximize(gaussian, new Float64Array([2, -3]), {
  maxIterations: 5_000,
});
// result.point ≈ [0, 0], result.value ≈ 1
```

#### Constrained optimization (penalty method)

Minimize the sphere function with box constraints:

$$f(x, y) = x^2 + y^2$$

- **Goal:** minimize
- **Subject to:** $2 \le x \le 5,\; 2 \le y \le 5$
- **Method:** quadratic penalty, $\mu = 10\,000$
- **Start:** $x_0 = (3, 3)$
- **Expected:** $\min = 8$ at $(2, 2)$

Constraints can be passed directly via settings - the penalty is applied
with the correct sign for both `minimize` and `maximize`:

```typescript
const sphere = (x: Float64Array): number => {
  let sum = 0;
  for (let i = 0; i < x.length; i++) sum += x[i] ** 2;
  return sum;
};

const nm = new NelderMead();
const result = nm.minimize(sphere, new Float64Array([3, 3]), {
  maxIterations: 5_000,
  constraints: boxConstraints(
    new Float64Array([2, 2]),  // lower bounds
    new Float64Array([5, 5]),  // upper bounds
  ),
  penaltyOptions: {mu: 10_000},
});
// result.point ≈ [2, 2], result.value ≈ 8
```

#### Custom constraints (inequality + equality)

$$f(x, y) = (x - 3)^2 + (y - 3)^2$$

- **Goal:** minimize
- **Subject to:** $x + y \le 4$ (inequality), $x = y$ (equality)
- **Method:** quadratic penalty, $\mu = 100\,000$
- **Start:** $x_0 = (0, 0)$
- **Expected:** $\min = 2$ at $(2, 2)$

```typescript
const fn = (x: Float64Array) => (x[0] - 3) ** 2 + (x[1] - 3) ** 2;

const constraints: singleObjective.Constraint[] = [
  {type: 'ineq', fn: (x) => x[0] + x[1] - 4},  // x + y <= 4
  {type: 'eq',   fn: (x) => x[0] - x[1]},        // x = y
];

const constrained = applyPenalty(fn, constraints, {mu: 100_000});

const nm = new NelderMead();
const result = nm.minimize(constrained, new Float64Array([0, 0]), {
  maxIterations: 10_000,
});
// result.point ≈ [2, 2], result.value ≈ 2
```

#### PSO with reproducible seed

Minimize the Rosenbrock function with PSO using a fixed seed for deterministic results:

$$f(x, y) = 100\,(y - x^2)^2 + (1 - x)^2$$

- **Goal:** minimize
- **Start:** $x_0 = (-1.2,\; 1.0)$
- **Expected:** $\min = 0$ at $(1, 1)$

```typescript
const pso = new PSO();
const result = pso.minimize(rosenbrock, new Float64Array([-1.2, 1.0]), {
  swarmSize: 40,
  maxIterations: 3_000,
  seed: 42,  // deterministic results
});
```

#### Async objective function

When the objective function involves asynchronous work (API calls, simulations, file I/O),
use `minimizeAsync` / `maximizeAsync`:

```typescript
const asyncObjective = async (x: Float64Array): Promise<number> => {
  // e.g. call a remote simulation service
  return rosenbrock(x);
};

const nm = new NelderMead();
const result = await nm.minimizeAsync(asyncObjective, new Float64Array([-1.2, 1.0]), {
  maxIterations: 5_000,
  tolerance: 1e-12,
});
```

Async penalty wrappers are also available via `applyPenaltyAsync`.

#### Iteration callback

Use `onIteration` to monitor progress or stop early:

```typescript
const result = nm.minimize(sphere, new Float64Array([5, -3, 7]), {
  maxIterations: 5_000,
  onIteration: (state) => {
    // Log every 100th iteration
    if (state.iteration % 100 === 0)
      console.log(`iter ${state.iteration}: best = ${state.bestValue}`);

    // Return true to stop early
    if (state.bestValue < 1e-6) return true;
  },
});
```

#### Registry

Optimizers can be looked up by name:

```typescript
console.log(listOptimizers()); // ['nelder-mead', 'pso']

const optimizer = getOptimizer('nelder-mead');
const result = optimizer.minimize(sphere, new Float64Array([5, -3, 7]), {
  maxIterations: 5_000,
});
```

### Benchmarks

Run all 15 standard test functions (Sphere, Rosenbrock, Beale, Booth, Matyas, Himmelblau, Three-Hump Camel,
Rastrigin, Ackley, Lévi N.13, Griewank, Styblinski-Tang, Easom, Goldstein-Price, McCormick) across all
4 optimizers with default hyperparameters:

```bash
npx tsx src/optimization/single-objective/benchmarks/unconstrained-benchmarks.ts
```

The output is a per-problem comparison table showing found value, distance to optimum, convergence status,
iteration count, function evaluations, and wall-clock time.

### Running examples

Runnable examples are located in `src/optimization/single-objective/examples/`:

```bash
# Unconstrained minimize & maximize (Rosenbrock, Sphere, Gaussian)
npx tsx src/optimization/single-objective/examples/unconstrained.ts

# Constrained optimization (box constraints, custom ineq/eq, log-barrier)
npx tsx src/optimization/single-objective/examples/constrained.ts

# Registry usage
npx tsx src/optimization/single-objective/examples/registry.ts

# Async objective functions and onIteration callbacks
npx tsx src/optimization/single-objective/examples/async-and-callbacks.ts
```
