---
name: add-optimizer
description: Add a new single-objective optimization algorithm to the sci-comp library. Use when the user asks to implement a new optimizer (e.g. gradient descent, BFGS, simulated annealing, differential evolution).
user-invocable: true
---

# Add a new single-objective optimizer

The user wants to add a new optimizer: **$ARGUMENTS**

Follow these steps exactly. Do not skip any step.

> **MANDATORY RULE — test failures require user approval before ANY changes:**
> After running tests (Step 7), if ANY tests fail, you MUST:
> 1. Compile the full list of failed tests (name + reason).
> 2. Present the list to the user and ask what to do next.
> 3. **Wait for the user's response** — do NOT proceed with fixes, adjustments, or any further steps until the user explicitly approves a course of action.
> This rule applies to every test run, including re-runs after fixes.

## Step 1: Understand the algorithm

Before writing code, research the algorithm:
- What are its tunable hyperparameters?
- Is it population-based or single-point?
- Does it require gradients?
- What is the iteration logic?

## Step 2: Read the reference implementation

Read the Nelder-Mead optimizer and its tests to understand the codebase patterns:
- `src/optimization/single-objective/optimizers/nelder-mead.ts` — reference optimizer implementation
- `src/optimization/single-objective/__tests__/nelder-mead.test.ts` — reference test file (the new optimizer MUST include all the same test cases)
- `src/optimization/single-objective/__tests__/helpers.ts` — test functions and helpers
- `src/optimization/single-objective/examples/unconstrained.ts` — reference example file

## Step 3: Create the optimizer file

Create `src/optimization/single-objective/optimizers/<name>.ts` following this structure:

```typescript
import {Optimizer} from '../optimizer';
import type {
  ObjectiveFunction,
  AsyncObjectiveFunction,
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

  // --- Synchronous path ---

  protected runInternal(
    fn: ObjectiveFunction,
    x0: Float64Array,
    s: <Name>Settings,
  ): OptimizationResult {
    const n = x0.length;
    const maxIter = s.maxIterations!;

    // Pre-allocate all buffers before the loop (see Array Operations rules below)
    const costHistory = new Float64Array(maxIter);
    let costLen = 0;
    const bestPoint = new Float64Array(n);
    // ... other scratch buffers ...
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
      costHistory: costHistory.subarray(0, costLen),
    };
  }

  // --- Asynchronous path ---

  protected async runInternalAsync(
    fn: AsyncObjectiveFunction,
    x0: Float64Array,
    s: <Name>Settings,
  ): Promise<OptimizationResult> {
    // Same algorithm logic as runInternal, but with `await fn(x)` instead of `fn(x)`.
    // Duplicate the loop structure; the only difference is `await` at function evaluation points.
    // Helper methods that call fn should have async variants (e.g. initSwarmAsync).
  }
}
```

### Rules

- **Float64Array** for all point vectors — never `number[]`
- **All settings fields optional** — `withDefaults()` fills them in
- Constructor calls `super('<kebab-name>')` — this becomes the registry key
- `runInternal` receives the already-penalized objective — do NOT handle constraints
- `runInternal` receives settings with defaults already applied via `withDefaults`
- **Both `runInternal` and `runInternalAsync` must be implemented** — the base class declares both as abstract
- `runInternalAsync` mirrors the sync logic but uses `await fn(x)` for evaluations
- Call `this.notify()` every iteration for callback support
- Track `costHistory` as best value per iteration
- Set `converged = true` only if a convergence criterion was met, not just maxIterations

### Array operations

Follow the patterns from `ARRAY-OPERATIONS.md` (see `implement-interactive-scientific-application-from-spec/references/reference/ARRAY-OPERATIONS.md`):

- **Pre-allocate and reuse:** All buffers (`costHistory`, `bestPoint`, scratch arrays) must be allocated once before the loop — zero allocations inside the iteration body
- **Out-parameter pattern:** Helper methods that compute into an array (e.g. gradient, centroid) must write into a caller-provided buffer, not allocate and return a new one
- **Logical length:** `costHistory` must be a pre-allocated `Float64Array(maxIter)` with a separate `costLen` counter; return `costHistory.subarray(0, costLen)` (zero-copy view)
- **Bulk copy with `set()`:** Use `bestPoint.set(x)` and `Float64Array.from(x0)` instead of manual loops
- **In-place transforms:** Gradient clipping, velocity updates, and similar transforms should modify arrays in-place when the previous values are no longer needed

## Step 4: Export from index.ts

Edit `src/optimization/single-objective/index.ts`:

1. Add export for the class and settings type (in the "Built-in optimizers" section):
```typescript
export {<Name>} from './optimizers/<name>';
export type {<Name>Settings} from './optimizers/<name>';
```

1. Add auto-registration (in the "Auto-register" section at the bottom):
```typescript
import {<Name>} from './optimizers/<name>';
registerOptimizer('<kebab-name>', () => new <Name>());
```

## Step 5: Write tests

Create `src/optimization/single-objective/__tests__/<name>.test.ts`.

**IMPORTANT:** The new test file MUST include ALL the same test cases as `nelder-mead.test.ts`. Read it first and replicate every `describe` block. The test structure is:

```typescript
import {<Name>, applyPenalty, applyPenaltyAsync, boxConstraints} from '..';
import type {Constraint} from '..';
import {
  rosenbrock, sphere, gaussian, quadratic3d,
  productSurface, quadraticMixed, negQuadraticMixed,
  ellipseConstraint, gaussianBump, expectPointClose, toAsync,
} from './helpers';

describe('<Name>', () => {
  const opt = new <Name>();

  // === ALL of the following test groups are REQUIRED ===
  // Each group must have both 'sync' and 'async' variants.
  // Settings (maxIterations, tolerance, learningRate, etc.) may be tuned
  // for the specific algorithm, but the test structure and assertions
  // must match nelder-mead.test.ts.

  // 1. minimize Rosenbrock 2D → min ≈ 0 at (1, 1)
  // 2. minimize Sphere 3D → min ≈ 0 at (0, 0, 0)
  // 3. maximize Gaussian 2D → max ≈ 1 at (0, 0)
  // 4. maximize Product Surface → max ≈ 1/27 at (1/3, 1/3)
  // 5. minimize Gaussian Bump → min ≈ 0 at (0, 0)
  // 6. maximize Gaussian Bump (x0 = [2, 2]) → max ≈ 2/e at (1, 0)
  // 7. maximize Gaussian Bump (x0 = [-3, -2]) → max ≈ 2/e at (-1, 0)
  // 8. minimize Sphere 2D with box constraints [2,5]² → min ≈ 8 at (2, 2)
  // 9. maximize Gaussian with box constraints [1,5]² via settings → max ≈ 0.1353 at (1, 1)
  // 10. minimize with custom ineq + eq constraints → min ≈ 2 at (2, 2)
  // 11. minimize Sphere 2D with log-barrier [2,5]² → min ≈ 8 at (2, 2)
  // 12. minimize Quadratic 3D → min ≈ -9 at (2, 1, 1)
  // 13. minimize x²+12xy+2y² s.t. 4x²+y²=25 → min ≈ -50 at (2, -3)
  // 14. maximize x²+12xy+2y² s.t. 4x²+y²=25 → max ≈ 106.25 at (1.5, 4)
  // 15. throws on empty x0

  // Refer to nelder-mead.test.ts for exact assertions and structure of each test.
  // You may adjust:
  //   - x0 starting points (if the algorithm needs a closer start)
  //   - maxIterations, tolerance, and algorithm-specific settings
  //   - Precision in toBeCloseTo / expectPointClose (if the algorithm is less precise)
  // You must NOT:
  //   - Remove any test group
  //   - Change expected values or expected points
});
```

Use test functions from `./helpers.ts` (rosenbrock, sphere, gaussian, quadratic3d, productSurface, quadraticMixed, negQuadraticMixed, ellipseConstraint, gaussianBump, etc.).
Use `toAsync(fn)` from `./helpers.ts` to wrap sync functions for async tests.
Use `expectPointClose(result, expected, tolerance)` for point assertions.
Use `toBeCloseTo(value, decimalDigits)` for value assertions.
For constrained async tests, use `applyPenaltyAsync(toAsync(fn), constraints, options)`.

## Step 6: Create example file

Create `src/optimization/single-objective/examples/<name>.ts` — a runnable example demonstrating the new optimizer.

Use `src/optimization/single-objective/examples/unconstrained.ts` as a reference for structure. The example file should:

1. Import the new optimizer from `..`
2. Define at least 2-3 test functions inline (Rosenbrock, Sphere, Gaussian — same as in the reference)
3. Show both `minimize` and `maximize` usage
4. Show at least one constrained example using `boxConstraints` + `applyPenalty`
5. Print results with `console.log`
6. Be runnable with `npx tsx src/optimization/single-objective/examples/<name>.ts`

## Step 7: Verify

Run in order:
1. `npm run lint-fix` — fix style issues
2. `npm run build` — must compile without errors
3. `npm test` — run all tests

**CRITICAL:** If any tests fail, do NOT silently fix or skip them. Instead:
1. Collect the full list of failing test names
2. Present the list to the user
3. Ask the user how to proceed (adjust settings, relax tolerances, skip specific tests, etc.)
4. Only continue after user confirmation

## Step 8: Update CLAUDE.md

In `CLAUDE.md`, update the architecture tree — add the new optimizer file under `optimizers/` with a short comment (class name and description).

## Step 9: Add to README.md

In `README.md`, add the new optimizer to the list under "Single-objective" section, with a Wikipedia or reference link.