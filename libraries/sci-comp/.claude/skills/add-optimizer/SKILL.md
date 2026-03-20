---
name: add-optimizer
description: Add a new single-objective optimization algorithm to the sci-comp library. Use when the user asks to implement a new optimizer (e.g. gradient descent, BFGS, simulated annealing, differential evolution).
user-invocable: true
---

# Add a new single-objective optimizer

The user wants to add a new optimizer: **$ARGUMENTS**

Follow these steps exactly. Do not skip any step.

## Step 1: Understand the algorithm

Before writing code, research the algorithm:
- What are its tunable hyperparameters?
- Is it population-based or single-point?
- Does it require gradients?
- What is the iteration logic?

## Step 2: Create the optimizer file

Create `src/optimization/single-objective/optimizers/<name>.ts` following this structure:

```typescript
import {Optimizer} from '../optimizer';
import type {
  ObjectiveFunction,
  OptimizationResult,
  CommonSettings,
} from '../types';

// 1. Settings interface — extend CommonSettings with algorithm-specific params
export interface <Name>Settings extends CommonSettings {
  // All fields optional with defaults documented in JSDoc
  /** Description (default: value). */
  param?: number;
}

// 2. Class — extend Optimizer<Settings>
export class <Name> extends Optimizer<<Name>Settings> {
  constructor() {
    super('<kebab-name>');  // registry key, lowercase with hyphens
  }

  protected withDefaults(s: <Name>Settings): <Name>Settings {
    return {
      maxIterations: 1000,
      tolerance: 1e-8,
      // ...algorithm defaults...
      ...s,
    };
  }

  protected runInternal(
    fn: ObjectiveFunction,
    x0: Float64Array,
    s: <Name>Settings,
  ): OptimizationResult {
    const n = x0.length;
    const costHistory: number[] = [];
    let converged = false;

    // ... algorithm logic ...
    // Use Float64Array for all point vectors
    // Call this.notify(s.onIteration, {iteration, bestValue, bestPoint}) each iteration
    // If notify returns true, stop early

    return {
      point: bestPoint,
      value: bestValue,
      iterations: iter,
      converged,
      costHistory: new Float64Array(costHistory),
    };
  }
}
```

### Rules

- **Float64Array** for all point vectors — never `number[]`
- **All settings fields optional** — `withDefaults()` fills them in
- Constructor calls `super('<kebab-name>')` — this becomes the registry key
- `runInternal` receives the already-penalized objective — do NOT handle constraints
- `runInternal` receives settings with defaults already applied via `withDefaults`
- Call `this.notify()` every iteration for callback support
- Track `costHistory` as best value per iteration
- Set `converged = true` only if a convergence criterion was met, not just maxIterations

## Step 3: Export from index.ts

Edit `src/optimization/single-objective/index.ts`:

1. Add export for the class and settings type (in the "Built-in optimizers" section):
```typescript
export {<Name>} from './optimizers/<name>';
export type {<Name>Settings} from './optimizers/<name>';
```

2. Add auto-registration (in the "Auto-register" section at the bottom):
```typescript
import {<Name>} from './optimizers/<name>';
registerOptimizer('<kebab-name>', () => new <Name>());
```

## Step 4: Write tests

Create `src/optimization/single-objective/__tests__/<name>.test.ts`:

```typescript
import {<Name>, boxConstraints} from '..';
import type {Constraint} from '..';
import {
  rosenbrock, sphere, gaussian, expectPointClose,
} from './helpers';

describe('<Name>', () => {
  const opt = new <Name>();

  // Minimum set of tests:
  it('minimize Rosenbrock 2D -> min = 0 at (1, 1)', () => { ... });
  it('minimize Sphere 3D -> min = 0 at (0, 0, 0)', () => { ... });
  it('maximize Gaussian 2D -> max = 1 at (0, 0)', () => { ... });
  it('minimize with box constraints', () => { ... });
  it('throws on empty x0', () => {
    expect(() => opt.minimize(sphere, new Float64Array([]), {})).toThrow();
  });
});
```

Use test functions from `./helpers.ts` (rosenbrock, sphere, gaussian, etc.).
Use `expectPointClose(result, expected, tolerance)` for point assertions.
Use `toBeCloseTo(value, decimalDigits)` for value assertions.

## Step 5: Verify

Run in order:
1. `npm run lint-fix` — fix style issues
2. `npm run build` — must compile without errors
3. `npm test` — all tests must pass (old + new)

If any step fails, fix and re-run.

## Step 6: Update CLAUDE.md

In `CLAUDE.md`, update the architecture tree — add the new optimizer file under `optimizers/` with a short comment (class name and description).

## Step 7: Add to README.md

In `README.md`, add the new optimizer to the list under "Single-objective" section, with a Wikipedia or reference link.
