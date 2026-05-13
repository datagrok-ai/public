# Single-Objective Optimization

Six built-in solvers, each supporting synchronous and asynchronous objective functions:

- [Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) - derivative-free simplex method
- [PSO (Particle Swarm Optimization)](https://en.wikipedia.org/wiki/Particle_swarm_optimization) - stochastic population-based method
- [Gradient Descent](https://en.wikipedia.org/wiki/Gradient_descent) - first-order method with numerical gradients, momentum, and learning rate decay
- [Adam](https://arxiv.org/abs/1412.6980) - adaptive moment estimation with per-parameter learning rates
- [L-BFGS](https://en.wikipedia.org/wiki/Limited-memory_BFGS) - limited-memory quasi-Newton method with two-loop recursion and Armijo line search
- [L-BFGS-B](https://doi.org/10.1137/0916069) - limited-memory BFGS with **native box constraints**, Moré–Thuente strong-Wolfe line search, compact representation (Byrd–Nocedal–Schnabel 1994), and the Morales–Nocedal (2011) subspace refinement

```typescript
import {singleObjective} from '@datagrok-libraries/sci-comp';

const {NelderMead, PSO, applyPenalty, applyPenaltyAsync, boxConstraints, getOptimizer, listOptimizers} = singleObjective;
```

| Method | Sync | Async |
|--------|------|-------|
| Minimize | `minimize(fn, x0, settings)` | `minimizeAsync(fn, x0, settings)` |
| Maximize | `maximize(fn, x0, settings)` | `maximizeAsync(fn, x0, settings)` |

## Unconstrained minimize

Minimize the [Rosenbrock function](https://en.wikipedia.org/wiki/Rosenbrock_function):

$$f(x, y) = 100(y - x^2)^2 + (1 - x)^2$$

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

## Unconstrained maximize

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

## Native box constraints with L-BFGS-B

L-BFGS-B handles box constraints `l ≤ x ≤ u` geometrically (no penalty). The
optimum sits exactly on active bounds, with no trade-off against objective
accuracy:

$$f(x, y) = 100(y - x^2)^2 + (1 - x)^2,\quad -1.5 \le x, y \le 0.5$$

- **Goal:** minimize
- **Subject to:** both coords in `[-1.5, 0.5]`
- **Expected:** $\min = 0.25$ at $(0.5,\; 0.25)$ (upper bound on `x` active)

```typescript
const {LBFGSB} = singleObjective;

const opt = new LBFGSB();
const result = opt.minimize(rosenbrock, new Float64Array([0, 0]), {
  bounds: {lower: -1.5, upper: 0.5},
  // Scalar lower/upper applies to every coordinate; arrays pick per-dim.
  // Omitted side defaults to ±Infinity.
  gradTolerance: 1e-8,
});
// result.point ≈ [0.5, 0.25], result.value ≈ 0.25
// All bounds satisfied exactly.
```

Per-dimension bounds, half-bounded variables, and fixed variables (`l == u`):

```typescript
const result = opt.minimize(fn, x0, {
  bounds: {
    lower: [0, -Infinity, 2, -1],    // variable 1 unbounded below, variable 2 fixed
    upper: [10, 5, 2, Infinity],     // variable 2 fixed (l == u), variable 3 unbounded above
  },
});
```

Supplying an analytic gradient (`gradFn`) activates the fast path — otherwise
central finite differences are used:

```typescript
const result = opt.minimize(fn, x0, {
  bounds: {lower: 0},
  gradFn: (x, gOut) => {
    for (let i = 0; i < x.length; i++) gOut[i] = 2 * x[i];
  },
});
```

Compared to running other optimizers with `boxConstraints()` + quadratic
penalty, L-BFGS-B on 7 bounded benchmarks ([bounded-benchmarks.md](./benchmarks/bounded-benchmarks.md)):

| Optimizer | Mode | Success | Feasibility |
|-----------|------|--------:|-------------|
| **L-BFGS-B** | native | **7/7** | exact |
| L-BFGS / PSO | penalty | 3/7 | violates by ~1e-3 |
| Adam / GradientDescent / Nelder-Mead | penalty | 2/7 | violates by ~1e-3 |

## Constrained optimization (penalty method)

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

## Custom constraints (inequality + equality)

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

## PSO with reproducible seed

Minimize the Rosenbrock function with PSO using a fixed seed for deterministic results:

$$f(x, y) = 100(y - x^2)^2 + (1 - x)^2$$

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

## Async objective function

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

## Iteration callback

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

## Registry

Optimizers can be looked up by name:

```typescript
console.log(listOptimizers()); // ['nelder-mead', 'pso']

const optimizer = getOptimizer('nelder-mead');
const result = optimizer.minimize(sphere, new Float64Array([5, -3, 7]), {
  maxIterations: 5_000,
});
```

## Benchmarks

Run 15 [standard test functions](./benchmarks/unconstrained-benchmark-functions.md) (Sphere, Rosenbrock, Beale, Booth, Matyas, Himmelblau, Three-Hump Camel,
Rastrigin, Ackley, Lévi N.13, Griewank, Styblinski-Tang, Easom, Goldstein-Price, McCormick) across all
6 optimizers with default hyperparameters:

```bash
npx tsx src/optimization/single-objective/benchmarks/unconstrained-benchmarks.ts
```

The output is a per-problem comparison table showing found value, distance to optimum, convergence status,
iteration count, function evaluations, and wall-clock time.

Find the results in [this summary](./benchmarks/unconstrained-benchmarks.md).

A complementary [**multi-start benchmark**](./benchmarks/multistart-benchmarks.md) runs every
optimizer from three carefully-chosen starting points per problem (baseline + adversarial
perturbation + near-optimum) to expose x₀-sensitivity and local-trap failure modes:

```bash
npx tsx src/optimization/single-objective/benchmarks/multistart-benchmarks.ts
```

A third runner — [**bounded benchmarks**](./benchmarks/bounded-benchmarks.md) —
compares L-BFGS-B's native box-constraint handling against the penalty layer
used by the other five optimizers across 7 bounded problems (active bounds,
fixed variables, half-bounded, and bounds-inactive sanity checks):

```bash
npx tsx src/optimization/single-objective/benchmarks/bounded-benchmarks.ts
```

> **Why multi-start matters.** A single-start comparison can make a local optimizer
> look better than it really is. For example, L-BFGS reaches the Rastrigin / Lévi N.13
> global minimum in one iteration in the single-start benchmark — but only because the
> chosen x₀ happens to be integer-aligned (so the `sin(kπxᵢ)` gradient contributions
> zero out and the landscape reduces to a quadratic bowl). A 0.1 perturbation of x₀
> destroys that effect and L-BFGS traps in a nearby local well. When picking a method
> for a **multimodal** objective, consult the multi-start tables or wrap a local
> optimizer in a multi-start strategy.

## Running examples

Runnable examples are located in `examples/`:

```bash
# Unconstrained minimize & maximize (Rosenbrock, Sphere, Gaussian)
npx tsx src/optimization/single-objective/examples/unconstrained.ts

# Constrained optimization (box constraints, custom ineq/eq, log-barrier)
npx tsx src/optimization/single-objective/examples/constrained.ts

# Registry usage
npx tsx src/optimization/single-objective/examples/registry.ts

# Async objective functions and onIteration callbacks
npx tsx src/optimization/single-objective/examples/async-and-callbacks.ts

# L-BFGS specific example (historySize tuning, comparison vs GD/Adam)
npx tsx src/optimization/single-objective/examples/lbfgs.ts
```
