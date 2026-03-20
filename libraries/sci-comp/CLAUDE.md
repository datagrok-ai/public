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
    types.ts                      # ObjectiveFunction, OptimizationResult, Constraint, CommonSettings
    optimizer.ts                  # Abstract Optimizer<S> base class (validation, penalty wiring, minimize/maximize)
    penalty.ts                    # applyPenalty(), boxConstraints() — constraint → penalized objective
    registry.ts                   # registerOptimizer/getOptimizer/listOptimizers — name-based lookup
    optimizers/
      nelder-mead.ts              # NelderMead extends Optimizer<NelderMeadSettings>
      pso.ts                      # PSO extends Optimizer<PSOSettings>
    __tests__/                    # Jest tests per optimizer + registry
    examples/                     # Runnable examples (npx tsx src/optimization/single-objective/examples/*.ts)
  multi-objectives/
    moead/                        # MOEA/D multi-objective optimizer (defs.ts, moead.ts, utils.ts)
```

### Key design patterns

- **Optimizer base class** (`optimizer.ts`): All solvers extend `Optimizer<S>`. Subclasses implement `runInternal()` and `withDefaults()`. The base class handles input validation, constraint penalty wrapping, and the minimize/maximize inversion.
- **Namespace re-exports**: The public API uses namespace re-exports (`singleObjective`, `multiObjective`) to avoid name collisions between submodules. Consumers import as `import {singleObjective} from '@datagrok-libraries/sci-comp'`.
- **Float64Array everywhere**: All point vectors use `Float64Array`, not `number[]`.
- **Registry pattern**: Optimizers self-register at import time via side-effect imports in `single-objective/index.ts`.

## Code Style

- ESLint extends `google` config with: 2-space indent, 120-char max line length, `curly: multi-or-nest`
- Single quotes, semicolons required
- TypeScript strict mode
