# CLAUDE.md — optimization

Single-objective and multi-objective solvers. No dependency on `datagrok-api`.

## Architecture

```
src/optimization/
  index.ts                          # Namespace re-exports + side-effect imports that self-register optimizers
  single-objective/
    types.ts                        # ObjectiveFunction, AsyncObjectiveFunction, OptimizationResult, Constraint, CommonSettings
    optimizer.ts                    # Abstract Optimizer<S> base class — validation, penalty wiring, minimize/maximize ± async
    penalty.ts                      # applyPenalty(), applyPenaltyAsync(), boxConstraints() — constraint → penalized objective
    registry.ts                     # registerOptimizer / getOptimizer / listOptimizers — name-based lookup
    optimizers/
      nelder-mead.ts                # NelderMead extends Optimizer<NelderMeadSettings>
      pso.ts                        # PSO extends Optimizer<PSOSettings>
      gradient-descent.ts           # GradientDescent extends Optimizer<GradientDescentSettings>
      adam.ts                       # Adam extends Optimizer<AdamSettings>
      lbfgs.ts                      # LBFGS extends Optimizer<LBFGSSettings> (quasi-Newton, two-loop recursion, Armijo)
      lbfgs-b/                      # L-BFGS-B (limited-memory BFGS with native box constraints)
        index.ts                    # LBFGSB class (public surface) + withDefaults validation
        types.ts                    # LBFGSBSettings, LBFGSBLineSearchSettings, LBFGSBBounds, BOUND_* constants
        driver.ts                   # runSync / runAsync outer loop (Cauchy → subspace → line search → memory update)
        line-search.ts              # dcsrch + dcstep — Moré–Thuente strong-Wolfe + NaN-bisection wrapper
        bfgs-mat.ts                 # BFGSMat compact representation (ring buffer, W-products, block Cholesky, solveM)
        bounds.ts                   # normalizeBounds, classifyBounds, project, projectedGradient, maxFeasibleStep
        cauchy.ts                   # Generalized Cauchy point + binary min-heap
        subspace.ts                 # Subspace minimisation (SMW direct primal) + Morales–Nocedal 2011 + 1997 truncation
    __tests__/                      # Jest tests per optimizer + registry (sync & async); helpers.ts hosts shared problems
    examples/                       # Runnable: npx tsx src/optimization/single-objective/examples/*.ts
    benchmarks/
      test-functions.ts             # Shared classical test functions + BOUNDED_PROBLEMS + HIMMELBLAU_MINIMA
      unconstrained-benchmarks.ts   # 15 problems × 1 x₀
      multistart-benchmarks.ts      # 15 problems × 3 x₀ (baseline + adversarial + near-optimum)
      bounded-benchmarks.ts         # 7 box-constrained problems: native bounds vs penalty
  multi-objectives/
    moead/                          # MOEA/D (defs.ts, moead.ts, utils.ts)
```

## Key design patterns

- **Optimizer base class** (`optimizer.ts`): subclasses implement `runInternal()`, `runInternalAsync()`, and `withDefaults()`. The base class handles input validation, constraint penalty wrapping, and the minimize→maximize sign inversion. Don't reimplement those concerns inside subclasses.
- **Sync + async API**: each optimizer exposes `minimize` / `maximize` (sync) and `minimizeAsync` / `maximizeAsync`. Penalty wrappers come in both flavours (`applyPenalty`, `applyPenaltyAsync`).
- **Registry pattern**: optimizers self-register at import time via side-effect imports in `single-objective/index.ts`. To add a solver, add it there or it won't be discoverable through `getOptimizer` / `listOptimizers`.
- **Iteration callbacks**: `onIteration` in settings supports progress monitoring and early stopping (return `true` to stop). Honour it in both sync and async run loops.
- **`Float64Array` vectors**: all point vectors use `Float64Array`, never `number[]`.

## Test structure

Group tests by problem; each `describe` has a `sync` and an `async` `it`:

```ts
describe('minimize Rosenbrock 2D → min ≈ 0 at (1, 1)', () => {
  it('sync', () => { ... });
  it('async', async () => { ... });
});
```

Shared problems (rosenbrock, sphere, gaussian, quadratic3d, …) live in `__tests__/helpers.ts`.

## Benchmarks

Three runners share objective functions through `benchmarks/test-functions.ts`:

- `unconstrained-benchmarks.ts` — one x₀ per problem, head-to-head table.
- `multistart-benchmarks.ts` — three x₀ per problem, success-rate summary across 45 runs.
- `bounded-benchmarks.ts` — 7 box-constrained problems. **L-BFGS-B** uses its native `settings.bounds`; everyone else routes through `boxConstraints()` + the quadratic-penalty layer. Success requires both numerical accuracy AND exact feasibility (`feas_vio < 1e-6`) — penalty methods typically violate by ~1e-3.

When adding a new optimizer, register it in **all three** runners and regenerate the `.md` reports. When adding a new test function, export it from `test-functions.ts` (don't duplicate the body in runner files); bounded problems go in `BOUNDED_PROBLEMS` there.

**x₀-sensitivity caveat.** Single-start tables can flatter local optimizers on multimodal problems. Example: L-BFGS solves Rastrigin / Lévi N.13 in one iteration because the baseline x₀ is integer-aligned and zeros out the `sin(kπxᵢ)` gradient terms; a 0.1 perturbation destroys that effect (visible in multi-start). Treat impressive single-start results on multimodal objectives as hypotheses to verify against `multistart-benchmarks.md`.

## Skills

- `/add-optimizer <algorithm name>` — scaffold a new single-objective optimizer (class, tests, registration, exports).
