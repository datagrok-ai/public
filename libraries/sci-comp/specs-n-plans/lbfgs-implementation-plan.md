# L-BFGS Optimizer — Implementation Plan

## Overview

L-BFGS (Limited-memory Broyden–Fletcher–Goldfarb–Shanno) is a quasi-Newton optimization method that approximates the inverse Hessian using a fixed-size history of past gradient and position changes. Unlike full BFGS which stores an n×n matrix, L-BFGS stores only `m` vector pairs (typically m=5–20), making it practical for large-scale problems while retaining superlinear convergence near the optimum.

This is the third gradient-based optimizer in sci-comp (after GradientDescent and Adam), and the natural extraction point for shared gradient utilities.

**Reference implementation:** SciPy `scipy.optimize._lbfgsb_py._minimize_lbfgsb` — Python wrapper around Fortran L-BFGS-B v3.0 by Byrd, Lu, Nocedal & Morales.

**References:**
- J. Nocedal, "Updating Quasi-Newton Matrices with Limited Storage" (1980), Mathematics of Computation, 35, 151, pp. 773–782.
- R. H. Byrd, P. Lu, J. Nocedal, "A Limited Memory Algorithm for Bound Constrained Optimization" (1995), SIAM J. Sci. Comput., 16, 5, pp. 1190–1208.
- D. C. Liu, J. Nocedal, "On the limited memory BFGS method for large scale optimization" (1989), Mathematical Programming, 45, pp. 503–528.

---

## Scope of the Port

### What to port (pure algorithm in TypeScript)

| SciPy component | Port target |
|---|---|
| Two-loop recursion (inverse Hessian–gradient product) | `twoLoopRecursion()` private method |
| Curvature history ring buffer (s, y pairs) | Circular buffer of `Float64Array` pairs |
| Wolfe line search (strong or Armijo backtracking) | `lineSearch()` private method |
| Projected gradient convergence check (`pgtol`, max-norm) | Max-norm gradient convergence gate — `max\|∂f/∂xᵢ\| < gradTolerance` (matches SciPy) |
| Finite-difference gradient (no user-supplied gradient) | Reuse existing `computeGradient` / `computeGradientAsync` pattern |

### What NOT to port

| SciPy component | Reason |
|---|---|
| Fortran `_lbfgsb` BLAS kernel | Replaced by pure TypeScript loops on `Float64Array` |
| Bound-constrained (L-BFGS-**B**) projected gradient / Cauchy point | sci-comp handles constraints via the penalty layer; native bounds can be a follow-up |
| `LbfgsInvHessProduct` LinearOperator | Not needed — we don't expose the Hessian approximation |
| `maxfun` (separate function evaluation cap) | sci-comp uses `maxIterations` uniformly |
| `factr` convergence factor | Map to `tolerance` (equivalent: `factr * eps ≈ tolerance`) |
| `disp` / `iprint` logging | sci-comp uses `onIteration` callback |

---

## Core Algorithm — Two-Loop Recursion

The key differentiator of L-BFGS is the two-loop recursion that computes the product H·g without ever forming the matrix H. Given history vectors s_k, y_k (k = 0..count-1, **chronological order** — index 0 is the oldest stored pair, index count-1 is the most recent):

```
q = g (current gradient)

// --- First loop (backward over history, newest → oldest) ---
for i = count-1 downto 0:
    ρ[i] = 1 / (y[i] · s[i])
    α[i] = ρ[i] · (s[i] · q)
    q = q − α[i] · y[i]

// --- Scaling (using the most recent pair) ---
γ = (s[count-1] · y[count-1]) / (y[count-1] · y[count-1])    // initial Hessian scaling
r = γ · q

// --- Second loop (forward over history, oldest → newest) ---
for i = 0 to count-1:
    β = ρ[i] · (y[i] · r)
    r = r + (α[i] − β) · s[i]

// r is now H·g. Search direction: d = -r
```

**Implementation note:** the actual storage is a ring buffer (see below) — the chronological iteration above is obtained by walking the buffer from `(histHead - 1 + m) % m` backward for `histCount` steps in the first loop, and by reversing that traversal in the second loop. `ρ[i]` and `α[i]` in the pseudocode refer to logical history positions, not ring-buffer slots; reuse the same `alpha` scratch `Float64Array` indexed by logical position.

The scaling factor γ adapts the initial Hessian approximation H₀ = γ·I based on the most recent curvature information, which is critical for good step sizes without manual learning rate tuning.

---

## Line Search

SciPy uses a Moré–Thuente line search satisfying the strong Wolfe conditions. For the initial port, implement **backtracking line search with Armijo condition** (simpler, well-tested, sufficient for most use cases):

```
Given: f, x, d (search direction), grad, c1=1e-4, c2=0.9, maxSteps=20

α = 1.0  (initial step size)
f0 = f(x)
slope = grad · d  (directional derivative, must be negative)

for step = 0 to maxSteps:
    x_new = x + α · d
    f_new = f(x_new)

    // Guard against NaN/Inf — objective may blow up for large α
    // (e.g. log/pow of out-of-domain values). NaN fails Armijo comparison
    // silently (all comparisons with NaN are false), so check explicitly.
    if !isFinite(f_new):
        α *= 0.5
        continue

    if f_new <= f0 + c1 · α · slope:     // Armijo (sufficient decrease)
        return (α, x_new, f_new)          // accept

    α *= 0.5  (backtrack)

return (α, x_new, f_new)  // return best found
```

**Follow-up (optional):** Replace with strong Wolfe line search (Moré–Thuente) for better performance on ill-conditioned problems. Can be added as a `lineSearchMethod: 'backtracking' | 'wolfe'` setting.

---

## Settings Interface

```typescript
export interface LBFGSSettings extends CommonSettings {
  /** Number of corrections to store (history size m). Higher = better Hessian
   *  approximation but more memory. Typical range: 3–20 (default: 10). */
  historySize?: number;
  /** Step size h for central finite-difference gradient (default: 1e-7). */
  finiteDiffStep?: number;
  /** Gradient convergence threshold — stop when `max|∂f/∂xᵢ| < gradTolerance`
   *  (infinity norm, matches SciPy's `pgtol`/`gtol`). Default: 1e-5. */
  gradTolerance?: number;
  /** Armijo sufficient decrease parameter c₁ ∈ (0, 1) (default: 1e-4). */
  c1?: number;
  /** Maximum number of line search steps per iteration (default: 20). */
  maxLineSearchSteps?: number;
  /** Initial step size for line search (default: 1.0). */
  initialStepSize?: number;
}
```

Note: no `learningRate`. L-BFGS determines step sizes automatically via the line search and Hessian scaling — this is the main advantage over GD/Adam.

---

## Implementation Constraints (from `add-optimizer` skill)

These rules come from `.claude/skills/add-optimizer/SKILL.md` and are binding — they apply in addition to everything described below.

### Zero allocations in the iteration body

All buffers must be allocated **once**, before the main loop, and reused across iterations:

- Core state: `x`, `grad`, `prevGrad`, `direction`, `xNew`, `gradNew` — each `Float64Array(n)`
- History ring buffer: `S`, `Y` as `Float64Array[m]`, each inner array `Float64Array(n)` pre-allocated
- `rho`, `alpha` — each `Float64Array(m)`, pre-allocated
- `costHistory` — `Float64Array(maxIterations + 1)` with a separate `costLen: number` counter; on return, slice with `costHistory.subarray(0, costLen)`

Private helpers (`twoLoopRecursion`, `lineSearch`, `lineSearchAsync`, `computeGradient`, `computeGradientAsync`) must NOT allocate — they receive pre-allocated buffers via parameters and write into them (out-parameter pattern).

### Line search return shape

`LineSearchResult` (containing `stepSize`, `accepted`, etc. — no `xNew` field) is a single reusable scratch object created before the loop. The line search writes the new point into a caller-provided `xNew` buffer and the new function value into a caller-provided slot (mutable struct or a single-element scratch `Float64Array(1)`). This avoids per-iteration object allocation.

### `converged` semantics

Set `result.converged = true` **only** when one of the two real convergence gates fires (gradient max-norm or relative function change). Hitting `maxIterations` is NOT convergence — in that case `converged = false` and the final `bestPoint` / `bestValue` is still returned.

---

## File Structure

```
src/optimization/single-objective/optimizers/lbfgs.ts    — LBFGSSettings + LBFGS class
src/optimization/single-objective/__tests__/lbfgs.test.ts — Jest tests
src/optimization/single-objective/examples/lbfgs.ts       — Runnable examples
```

Changes to existing files:

- `src/optimization/single-objective/index.ts` — add export and registration
- `src/optimization/single-objective/README.md` — add to solver list and examples
- `src/optimization/single-objective/benchmarks/unconstrained-benchmarks.ts` — add LBFGS to benchmark runner

---

## Class Skeleton

```typescript
export class LBFGS extends Optimizer<LBFGSSettings> {
  constructor() { super('L-BFGS'); }

  protected withDefaults(s: LBFGSSettings): LBFGSSettings { ... }
  protected runInternal(fn, x0, s): OptimizationResult { ... }
  protected async runInternalAsync(fn, x0, s): Promise<OptimizationResult> { ... }

  // --- Private helpers ---
  private computeGradient(fn, x, grad, h, n): void { ... }
  private async computeGradientAsync(fn, x, grad, h, n): Promise<void> { ... }
  private twoLoopRecursion(grad, S, Y, rho, start, count, n, direction): void { ... }
  private lineSearch(fn, x, direction, grad, fx, c1, maxSteps, n): LineSearchResult { ... }
  private async lineSearchAsync(fn, x, direction, grad, fx, c1, maxSteps, n): Promise<LineSearchResult> { ... }
}
```

---

## Algorithm — Sync Path (`runInternal`)

### State

| Variable | Type | Description |
|---|---|---|
| `x` | `Float64Array(n)` | Current parameter vector |
| `grad` | `Float64Array(n)` | Current gradient |
| `prevGrad` | `Float64Array(n)` | Previous gradient (for computing y_k) |
| `direction` | `Float64Array(n)` | Search direction (output of two-loop recursion) |
| `S` | `Float64Array[m]` | Ring buffer of s_k = x_{k+1} - x_k vectors |
| `Y` | `Float64Array[m]` | Ring buffer of y_k = grad_{k+1} - grad_k vectors |
| `rho` | `Float64Array(m)` | Cached 1 / (y_k · s_k) values |
| `alpha` | `Float64Array(m)` | Scratch for two-loop recursion |
| `bestPoint` | `Float64Array(n)` | Best point found so far |
| `bestValue` | `number` | Best objective value |
| `histCount` | `number` | Number of stored history pairs (0..m) |
| `histHead` | `number` | Ring buffer insertion index |
| `xNew` | `Float64Array(n)` | Scratch: trial point from line search |
| `gradNew` | `Float64Array(n)` | Scratch: gradient at trial point |
| `costHistory` | `Float64Array(maxIter+1)` | Pre-allocated cost history buffer |
| `costLen` | `number` | Number of filled entries in `costHistory` |
| `lineSearchResult` | object | Reused scratch struct for line-search status (no `xNew` field) |

### Per-Iteration Steps

**Important:** f(x) and ∇f(x) are carried across iteration boundaries. On iteration 0 we evaluate both once up front; on subsequent iterations we reuse `fNew` and `gradNew` produced by the line search and post-line-search gradient of the *previous* iteration. Each finite-difference gradient costs 2n objective evaluations — re-computing it at the top of every iteration would double the cost unnecessarily.

```
Setup (before loop):
    fCur  = f(x)
    gradCur = ∇f(x)  (2n evaluations)

Loop (iteration = 0, 1, ...):
 1. Check gradient convergence (max-norm, matches SciPy):
      gNormInf = max_i |gradCur[i]|
      if gNormInf < gradTolerance → converged
 2. Update bestPoint / bestValue if fCur improved
 3. Record costHistory: `costHistory[costLen++] = fCur` (pre-allocated buffer, no push)
 4. Fire onIteration callback (may request early stop)
 5. Compute search direction:
      if histCount == 0:
        direction = -gradCur                 // steepest descent for first iteration
      else:
        twoLoopRecursion(gradCur, S, Y, rho, histHead, histCount, n, direction)
        negate direction                     // d = -H·g
 6. Verify descent: gradCur · direction < 0
      if not: reset history (histCount = 0), direction = -gradCur
 7. Line search along direction → (stepSize, xNew, fNew)
 8. Compute gradNew at xNew (2n evaluations)
 9. Update history:
      s_k = xNew - x
      y_k = gradNew - gradCur
      ys  = y_k · s_k
      yy  = y_k · y_k
      if ys > 1e-10 · max(1, yy):            // relative curvature condition
        store s_k, y_k in ring buffer
        rho[head] = 1 / ys
        advance histHead, cap histCount at m
10. Check function value convergence:
      |fNew - fCur| / max(|fCur|, |fNew|, 1) < tolerance → converged
11. x = xNew, fCur = fNew, gradCur = gradNew
12. iteration++
```

**Return value:** `costHistory` in the final `OptimizationResult` must be `costHistory.subarray(0, costLen)`. Set `converged = true` only if gate 1 or gate 2 fired (or `onIteration` returned `true`); exiting via `maxIterations` yields `converged = false`.

### Ring Buffer Design

The history uses a circular buffer to avoid array shifts:

```typescript
const m = s.historySize!;
const S: Float64Array[] = Array.from({length: m}, () => new Float64Array(n));
const Y: Float64Array[] = Array.from({length: m}, () => new Float64Array(n));
const rho = new Float64Array(m);
let histHead = 0;   // next write position
let histCount = 0;   // filled slots (0..m)

// Store new pair:
S[histHead].set(s_k);
Y[histHead].set(y_k);
rho[histHead] = 1 / ys;
histHead = (histHead + 1) % m;
if (histCount < m) histCount++;
```

The two-loop recursion iterates backward from `(histHead - 1 + m) % m` for `histCount` steps.

---

## Algorithm — Async Path (`runInternalAsync`)

Identical logic with `await fn(x)` and `await this.computeGradientAsync(...)` and `await this.lineSearchAsync(...)`.

---

## Gradient Utility Extraction (Refactoring)

This is the third gradient-based optimizer. Per the convention established in the Adam plan, extract shared gradient code now:

```typescript
// New file: src/optimization/single-objective/gradient-utils.ts

/** Central finite-difference gradient: ∂f/∂xᵢ ≈ (f(x+hᵢ) − f(x−hᵢ)) / 2h */
export function computeGradient(
  fn: ObjectiveFunction, x: Float64Array, grad: Float64Array, h: number, n: number
): void { ... }

/** Async variant. */
export async function computeGradientAsync(
  fn: AsyncObjectiveFunction, x: Float64Array, grad: Float64Array, h: number, n: number
): Promise<void> { ... }
```

**Decision:** This extraction is *optional* for the initial PR. The Adam plan deferred extraction ("extract when the third gradient optimizer arrives"), so this is the trigger point. However, it can be done as a separate refactoring commit to keep the L-BFGS PR focused. Recommend: do the extraction in the same PR, separate commit.

**Critical:** when extracting, update *all three* optimizers (`GradientDescent`, `Adam`, `LBFGS`) to import from `gradient-utils.ts` in the same commit. Leaving two copies of `computeGradient` (one in Adam, one imported) is worse than not extracting at all — it forks the definition and defeats the purpose. The extraction commit is a pure refactoring: no behavior changes, existing tests must continue to pass unchanged.

---

## Convergence Criteria

Two gates (any one triggers `converged = true`):

1. **Gradient max-norm:** `max_i |∂f/∂xᵢ| < gradTolerance` (default 1e-5, matches SciPy's `pgtol`/`gtol` — SciPy compares the infinity norm of the *projected* gradient; for unconstrained problems the projected gradient equals the raw gradient, so this is numerically compatible with SciPy out of the box)
2. **Function value:** `|f_{k+1} - f_k| / max(|f_k|, |f_{k+1}|, 1) < tolerance` (matches SciPy's `ftol`/`factr` criterion)

**Why max-norm (L∞), not Euclidean (L₂):** SciPy's `projgr` routine ([src/lbfgsb.c:2751](https://github.com/scipy/scipy/blob/main/scipy/optimize/src/lbfgsb.c#L2751)) computes `max_i |proj g_i|`. Using L₂ here would produce convergence thresholds that scale as √n relative to SciPy — the same tolerance value would mean different things at different dimensionalities, breaking cross-implementation comparison. Max-norm is also scale-invariant per coordinate, which is the property users expect from a gradient-based stopping criterion.

**Deliberate departure from GD/Adam:** L-BFGS does **not** use a `noImprovementLimit` streak counter. The function-value relative change criterion above is more appropriate for quasi-Newton methods, which take aggressive steps and can legitimately see no improvement for one iteration without being stuck. Document this in the `LBFGSSettings` TSDoc so users do not expect `noImprovementLimit` to be a valid option.

---

## Default Values Summary

| Parameter | Default | Rationale |
|---|---|---|
| `historySize` | `10` | SciPy default `maxcor=10` |
| `finiteDiffStep` | `1e-7` | Consistent with GD/Adam |
| `gradTolerance` | `1e-5` | SciPy default `gtol=1e-5` — applied as `max_i \|gᵢ\| < gradTolerance` (infinity norm, SciPy-compatible) |
| `c1` | `1e-4` | Standard Armijo parameter |
| `maxLineSearchSteps` | `20` | SciPy default `maxls=20` |
| `initialStepSize` | `1.0` | Standard for quasi-Newton |
| `maxIterations` | `1_000` | L-BFGS converges much faster than GD; SciPy uses 15000 but that's generous |
| `tolerance` | `1e-8` | Consistent with other optimizers |

---

## Registration

```typescript
// index.ts additions:
export {LBFGS} from './optimizers/lbfgs';
export type {LBFGSSettings} from './optimizers/lbfgs';

import {LBFGS} from './optimizers/lbfgs';
registerOptimizer('l-bfgs', () => new LBFGS());
```

---

## Edge Cases and Robustness

| Case | Handling |
|---|---|
| First iteration (no history) | Fall back to steepest descent: d = −∇f |
| Curvature condition fails (y·s ≤ 1e-10·max(1, y·y)) | Skip history update for this iteration |
| Line search fails (no Armijo step found) | Accept smallest step, do not update history |
| Line search produces NaN/Inf f(x_new) | Continue backtracking (α *= 0.5); do not treat NaN as accepted |
| Non-descent direction after two-loop | Reset history, use steepest descent for this step |
| Penalty boundary (box-constraint etc.) | Central FD through kinked penalty can produce spiky gradients — expect ys to fail curvature check near the boundary, which correctly skips the update |
| 1D problem (n=1) | Works normally — two-loop recursion is dimension-agnostic |
| Very large gradient (numerical issues) | Clip gradient norm before two-loop (same as GD/Adam) |
| NaN in gradient | Propagate as non-convergence — max-norm computation must check `isNaN` and treat it as "not converged" (do not return NaN from the compare); SciPy's `projgr` propagates NaN explicitly ([src/lbfgsb.c:2781](https://github.com/scipy/scipy/blob/main/scipy/optimize/src/lbfgsb.c#L2781)) |

---

## Implementation Steps (ordered)

### Step 1 — Core optimizer file

Create `src/optimization/single-objective/optimizers/lbfgs.ts`:
- `LBFGSSettings` interface
- `LBFGS` class with `withDefaults`, `runInternal`, `runInternalAsync`
- Private: `twoLoopRecursion`, `lineSearch`/`lineSearchAsync`, `computeGradient`/`computeGradientAsync`

### Step 2 — Registration and exports

Update `src/optimization/single-objective/index.ts`:
- Add `export {LBFGS}` and `export type {LBFGSSettings}`
- Add `registerOptimizer('l-bfgs', () => new LBFGS())`

### Step 3 — Tests

Create `src/optimization/single-objective/__tests__/lbfgs.test.ts`.

**The file MUST replicate every `describe` block from `nelder-mead.test.ts` (14 groups)** — the `add-optimizer` skill requires this verbatim. Each group has `it('sync')` + `it('async')`. Permitted adjustments: `x0`, `maxIterations`, algorithm-specific settings, and `toBeCloseTo` precision. **Forbidden**: removing a group, changing expected optimum values or expected points.

The required 14 groups (matching `nelder-mead.test.ts` at the time of writing):

1. minimize Rosenbrock 2D → min ≈ 0 at (1, 1)
2. minimize Sphere 3D → min ≈ 0 at (0, 0, 0)
3. maximize Gaussian 2D → max ≈ 1 at (0, 0)
4. maximize Product Surface → max ≈ 1/27 at (1/3, 1/3)
5. minimize Gaussian Bump → min ≈ 0 at (0, 0)
6. maximize Gaussian Bump (x0 = [2, 2]) → max ≈ 2/e at (1, 0)
7. maximize Gaussian Bump (x0 = [-3, -2]) → max ≈ 2/e at (-1, 0)
8. minimize Sphere 2D with box constraints [2,5]² → min ≈ 8 at (2, 2)
9. maximize Gaussian with box constraints [1,5]² via settings → max ≈ 0.1353 at (1, 1)
10. minimize with custom ineq + eq constraints → min ≈ 2 at (2, 2)
11. minimize Sphere 2D with log-barrier [2,5]² → min ≈ 8 at (2, 2)
12. minimize Quadratic 3D → min ≈ -9 at (2, 1, 1)
13. minimize x²+12xy+2y² s.t. 4x²+y²=25 → min ≈ -50 at (2, -3)
14. maximize x²+12xy+2y² s.t. 4x²+y²=25 → max ≈ 106.25 at (1.5, 4)

Re-verify the list against the current state of `nelder-mead.test.ts` before starting — if new groups have been added there, include them too.

In addition to the verbatim replication, add L-BFGS-specific coverage:
- Empty x0 throws
- Verify convergence in significantly fewer iterations than GD
- History reset via non-descent path (staircase function — see Testing Checklist)
- NaN safety in line search (log of out-of-domain argument)
- Max-norm convergence parity with SciPy (see Testing Checklist)

### Step 4 — Example

Create `src/optimization/single-objective/examples/lbfgs.ts`:
- Unconstrained Rosenbrock, Sphere, Gaussian
- Constrained with box constraints
- Show `historySize` tuning
- Compare iteration count vs GD/Adam on same problem

### Step 5 — Benchmarks

Add LBFGS to `benchmarks/unconstrained-benchmarks.ts`:
- Import `LBFGS`
- Add an entry to the `optimizers` array:
  ```typescript
  {
    name: 'L-BFGS',
    optimizer: new LBFGS(),
    settings: {maxIterations: 1_000, historySize: 10},
  },
  ```
- Run `npx tsx src/optimization/single-objective/benchmarks/unconstrained-benchmarks.ts`
- Regenerate `src/optimization/single-objective/benchmarks/unconstrained-benchmarks.md` as markdown tables (same format as existing file)
- Update the problem-count banner in the runner / `.md` file if it has changed

### Step 6 — Documentation

Update `README.md`:
- Add L-BFGS to solver list
- Add usage example

### Step 7 (optional) — Gradient utility extraction

Create `gradient-utils.ts`, update GD/Adam/LBFGS to import from it.
Separate commit for clean diff.

---

## Testing Checklist

- [ ] **Sphere 2D/3D** — should converge in <50 iterations (vs ~thousands for GD)
- [ ] **Rosenbrock 2D** — verify it handles the narrow valley; expect <200 iterations
- [ ] **Gaussian maximize** — verify via `maximize()` / `maximizeAsync()`
- [ ] **Quadratic 3D** — verify on non-trivially coupled variables
- [ ] **ProductSurface maximize** — verify maximize path
- [ ] **GaussianBump** — verify multiple local optima (from different starting points)
- [ ] **Box constraints via penalty** — Sphere with [2,5]² → min at (2,2). Also assert that curvature-skip near the boundary does not prevent eventual convergence.
- [ ] **Custom ineq+eq constraints** — same test case as GD
- [ ] **Async path** — verify identical results with `toAsync()` wrapper
- [ ] **1D problem** — simple parabola, verify it works
- [ ] **High-dimensional** — Sphere 10D or 20D, verify it still converges fast
- [ ] **History reset** — verify non-descent direction triggers reset. Concrete pathological target: `f(x) = Math.floor(10 * (x[0]*x[0] + x[1]*x[1])) / 10` — piecewise-constant staircase where central FD gives near-zero gradient on flat terraces and spikes at step boundaries, reliably producing non-descent directions after a few history updates.
- [ ] **NaN/Inf safety** — line search must not accept a step that produces NaN (e.g. `f(x) = Math.log(x[0])` started from x=1 with a large initial direction)
- [ ] **Iteration callback** — verify `onIteration` fires, early stop via `return true`
- [ ] **Empty x0** — throws error
- [ ] **Benchmark comparison** — on the 15 standard functions, L-BFGS should outperform GD/Adam on smooth unimodal problems, may struggle on multimodal (Rastrigin, Ackley)
- [ ] **Max-norm convergence parity with SciPy** — on Sphere/Rosenbrock with identical `x0`, `historySize=10`, `gradTolerance=1e-5`, `tolerance=factr*eps`, iteration count and final `max|gᵢ|` should be in the same order of magnitude as `scipy.optimize.minimize(method='L-BFGS-B')`. Small differences from line-search method (Armijo vs Moré–Thuente) are expected; large (≫10×) divergences indicate a bug in the stopping criterion or Hessian scaling.

---

## Expected Performance Characteristics

| Problem type | L-BFGS vs GD/Adam |
|---|---|
| Smooth, convex (Sphere, Quadratic) | Much faster (10–100× fewer iterations) |
| Smooth, non-convex unimodal (Rosenbrock) | Faster (5–20× fewer iterations) |
| Ill-conditioned (elongated valleys) | Much better than GD, comparable to Adam |
| Multimodal (Rastrigin, Ackley) | Worse than PSO — L-BFGS is local, will find nearest local minimum |
| Noisy objectives | Worse than Adam — curvature estimates are corrupted by noise |
| Very high-dimensional (n > 1000) | Sweet spot for L-BFGS — memory scales as O(m·n) vs O(n²) for full BFGS |

---

## Future Extensions (out of scope for initial PR)

1. **Native bound constraints (L-BFGS-B):** Projected gradient + generalized Cauchy point computation. Would avoid the penalty method overhead for simple box constraints.
2. **Strong Wolfe line search (Moré–Thuente):** Better theoretical convergence guarantees. Add as `lineSearchMethod: 'backtracking' | 'wolfe'` setting.
3. **User-supplied gradient:** Accept optional `gradient: (x) => Float64Array` in settings to skip finite differences. Major performance win for problems where analytic gradients are available.
4. **Damped L-BFGS:** Modified update rule that ensures positive definiteness even when curvature condition fails (Powell's damping).
