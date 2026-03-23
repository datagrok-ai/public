# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`@datagrok-libraries/sci-comp` is a pure TypeScript library of numerical methods for the Datagrok platform. It has **no dependency on datagrok-api**.

## Commands

```bash
npm run build        # tsc → dist/
npm run lint         # ESLint (Google style, 120 char max-len, 2-space indent)
npm run lint-fix     # ESLint with --fix
npm test             # Jest (ts-jest, all __tests__/*.test.ts files)
npx jest --testPathPattern nelder-mead   # Run a single test file
```

## Architecture

The library is organized by numerical domain, currently **optimization**:

```
index.ts                          # Entry point: re-exports {singleObjective, multiObjective} namespaces
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
```

### Key design patterns

- **Optimizer base class** (`optimizer.ts`): All solvers extend `Optimizer<S>`. Subclasses implement `runInternal()`, `runInternalAsync()`, and `withDefaults()`. The base class handles input validation, constraint penalty wrapping, and the minimize/maximize inversion.
- **Sync + async API**: Each optimizer exposes `minimize`/`maximize` (sync) and `minimizeAsync`/`maximizeAsync` (for async objective functions). Penalty wrappers have sync (`applyPenalty`) and async (`applyPenaltyAsync`) variants.
- **Namespace re-exports**: The public API uses namespace re-exports (`singleObjective`, `multiObjective`) to avoid name collisions between submodules. Consumers import as `import {singleObjective} from '@datagrok-libraries/sci-comp'`.
- **Float64Array everywhere**: All point vectors use `Float64Array`, not `number[]`.
- **Registry pattern**: Optimizers self-register at import time via side-effect imports in `single-objective/index.ts`.
- **Iteration callbacks**: `onIteration` callback in settings allows progress monitoring and early stopping (return `true` to stop).

### Test structure

Tests are grouped by problem, each containing `sync` and `async` variants:

```
describe('minimize Rosenbrock 2D → min ≈ 0 at (1, 1)', () => {
  it('sync', () => { ... });
  it('async', async () => { ... });
});
```

## Skills

- `/add-optimizer <algorithm name>` — scaffold a new single-objective optimizer (class, tests, registration, exports)

## Code Style

- ESLint extends `google` config with: 2-space indent, 120-char max line length, `curly: multi-or-nest`
- Single quotes, semicolons required
- TypeScript strict mode
