# L-BFGS-B in TypeScript: Implementation Specification

This document is a self-contained, implementation-ready specification for **L-BFGS-B** — the limited-memory BFGS method with box constraints — in pure TypeScript. It combines the mathematics of Byrd, Lu, Nocedal & Zhu (1995), the ACM TOMS 778 algorithm of Zhu et al. (1997), the Morales–Nocedal (2011) correction, the Moré–Thuente line search of MINPACK-2, and the implementation patterns of SciPy, L-BFGS-B-C, and LBFGSpp. An engineer (or Claude Code) should be able to produce a correct, tested, scipy-faithful implementation from this spec alone.

---

## 0. Integration with `@datagrok-libraries/sci-comp`

This section is the **delta** between the generic specification (§1–§13, Appendices A–C) and the form in which L-BFGS-B is delivered inside this library. Where §0 and §1–§13 disagree, §0 wins. Where §0 is silent, the generic spec applies verbatim.

### 0.1 Decision log

| # | Decision | Choice |
|---|---|---|
| 1 | Class / file path | `class LBFGSB` in `src/optimization/single-objective/optimizers/lbfgs-b/` (subfolder, not single file) |
| 2 | Integration pattern | `extends Optimizer<LBFGSBSettings>`; returns `OptimizationResult`; both sync + async paths required |
| 3 | Box-constraints API | Field `bounds?: { lower?: number \| ArrayLike<number>; upper?: number \| ArrayLike<number> }` in `LBFGSBSettings`; orthogonal to `settings.constraints` (ineq/eq via penalty layer) |
| 4 | Gradient source | Optional `gradFn?` / `gradFnAsync?` in settings; fallback to central finite differences using `finiteDiffStep` (default `1e-7`) |
| 5 | Settings & defaults | See §0.3 table below |
| 6 | Convergence semantics | `converged = true` only on T1 (projected-gradient) or T2 (ΔfRel); all other exits → `converged = false`; `point` / `value` = best-so-far + `project()` (feasibility invariant for users) |
| 7 | Line-search strategy | Full Moré–Thuente from commit 1; no Armijo intermediate version |
| 8 | File layout | Subfolder `lbfgs-b/` with `index.ts`, `types.ts`, `driver.ts`, `cauchy.ts`, `subspace.ts`, `line-search.ts`, `bfgs-mat.ts`, `bounds.ts` |
| 9 | Linalg organisation | Vector ops inlined at call sites; `potrf`/`trsv`/W-specialised products live private inside `bfgs-mat.ts`; no separate `linalg.ts` |
| 10 | Benchmarks | Register in existing `unconstrained-benchmarks.ts` + `multistart-benchmarks.ts`; add new `bounded-benchmarks.ts` (6–10 problems) comparing L-BFGS-B (native bounds) vs others (bounds via penalty) |
| 11 | Phasing | One PR, 9 logical commits (each CI-green independently); commit order in §0.7 |

### 0.2 Public API

```ts
import { LBFGSB } from '@datagrok-libraries/sci-comp/.../optimizers/lbfgs-b';

const opt = new LBFGSB();
const result = opt.minimize(fn, x0, {
  bounds: { lower: 0, upper: [1, 5, 10] },
  gradFn: (x, gOut) => { /* write ∇f(x) into gOut */ },
  historySize: 10,
  gradTolerance: 1e-5,
  tolerance: 1e-8,
  maxIterations: 1000,
  maxFunctionEvaluations: 15000,
  lineSearch: { ftol: 1e-3, gtol: 0.9, xtol: 0.1, maxSteps: 20 },
  onIteration: ({ iteration, bestValue, extra }) => {
    // extra: { projGradInfNorm, functionEvaluations, stepSize,
    //          lineSearchSteps, historyCount, activeBounds }
    return false;
  },
});
// result: { point, value, iterations, converged, costHistory }
```

No `minimize(fg, x0, bounds, options)` top-level function. No `FgFn`, `LBFGSBOptions`, `LBFGSBResult`, `ExitStatus`, `IterationInfo` — these are replaced by the library's `ObjectiveFunction`, `LBFGSBSettings` (extends `CommonSettings`), `OptimizationResult`, and `IterationState`. The `maximize`, `minimizeAsync`, `maximizeAsync` entry points are inherited from `Optimizer` for free.

### 0.3 Settings type and defaults

```ts
export interface LBFGSBSettings extends CommonSettings {
  bounds?: { lower?: number | ArrayLike<number>; upper?: number | ArrayLike<number> };
  gradFn?: (x: Float64Array, gOut: Float64Array) => void;
  gradFnAsync?: (x: Float64Array, gOut: Float64Array) => Promise<void>;
  historySize?: number;                    // m, default 10
  gradTolerance?: number;                   // pgtol, default 1e-5
  maxFunctionEvaluations?: number;          // maxfun, default 15000
  finiteDiffStep?: number;                  // FD h, default 1e-7
  lineSearch?: {
    ftol?: number;                          // Armijo c₁, default 1e-3
    gtol?: number;                          // curvature c₂, default 0.9
    xtol?: number;                          // interval tol, default 0.1
    maxSteps?: number;                      // maxls, default 20
  };
}
```

Parameter cross-reference (supersedes Appendix C):

| SciPy | Generic spec §2.2 | This library |
|---|---|---|
| `maxcor` | `m` | `historySize` |
| `factr · εmach` | `factr` | `tolerance` (direct threshold, not factr-scaled) |
| `gtol` (pgtol) | `pgtol` | `gradTolerance` |
| `maxiter` | `maxIterations` | `maxIterations` (**default 1000**, not 15000) |
| `maxfun` | `maxFunctionEvaluations` | `maxFunctionEvaluations` |
| `maxls` | `maxLineSearch` | `lineSearch.maxSteps` |
| — | `ftol` (line search) | `lineSearch.ftol` |
| — | `gtol` (curvature) | `lineSearch.gtol` |
| — | `xtol` | `lineSearch.xtol` |
| `bounds` (positional) | `Bounds` arg | `settings.bounds` |

Key divergence from SciPy: `maxIterations` default is `1000` (consistent with `GradientDescent` / `Adam` / `LBFGS` in this library), not `15000`. Benchmarks and SciPy-parity tests override explicitly.

### 0.4 Convergence tests (supersedes §11)

| # | Test | Condition | Sets `converged`? |
|---|---|---|---|
| T1 | Projected-gradient | `‖P(x − g, l, u) − x‖∞ ≤ gradTolerance` | **true** |
| T2 | Relative Δf | `\|f_k − f_{k+1}\| / max(\|f_k\|, \|f_{k+1}\|, 1) ≤ tolerance` | **true** |
| T3 | Iteration cap | `iter ≥ maxIterations` | false |
| T4 | Evaluation cap | `nfev ≥ maxFunctionEvaluations` | false |
| T5 | Line-search failure + retry failure | §9 recovery exhausted | false |
| T6 | NaN/Inf in `f` or `g` mid-run | `throw` at x0, return at runtime | false |
| T7 | User callback returned `true` | — | false |

Per-iteration order: (0) T1 on initial point; (after step) T2 → T1 (on new g) → T3/T4 → callback → T7.

### 0.5 File layout (supersedes §2.1)

```
src/optimization/single-objective/optimizers/lbfgs-b/
  index.ts          // exports class LBFGSB + LBFGSBSettings only
  types.ts          // LBFGSBSettings, internal BoundType enum
  driver.ts         // runInternal / runInternalAsync (outer loop, §4)
  cauchy.ts         // generalised Cauchy point (§5) + binary min-heap
  subspace.ts       // subspace minimisation + Morales–Nocedal backtrack (§6)
  line-search.ts    // dcsrch + dcstep (§7)
  bfgs-mat.ts       // BFGSMat + compact-form products + private potrf/trsv (§3)
  bounds.ts         // classifyBounds, project, projectedGradient, maxFeasibleStep (§8)
```

Tests sit at the existing `optimizers/../__tests__/` level (not inside `lbfgs-b/`):

```
__tests__/
  lbfgs-b.test.ts                 // integration: Rosenbrock/Beale/bounded/...
  lbfgs-b-cauchy.test.ts
  lbfgs-b-subspace.test.ts
  lbfgs-b-line-search.test.ts     // regression vs Fortran dcstep tabulated outputs
  lbfgs-b-bfgs-mat.test.ts
  lbfgs-b-bounds.test.ts
```

Registration in `single-objective/index.ts`:

```ts
import { LBFGSB } from './optimizers/lbfgs-b';
registerOptimizer('L-BFGS-B', () => new LBFGSB());
```

### 0.6 Benchmarks (supersedes §12.2 defaults)

L-BFGS-B is registered in **three** runners:

1. `unconstrained-benchmarks.ts` — existing. Sanity check: L-BFGS-B with no bounds ≈ L-BFGS quality or better.
2. `multistart-benchmarks.ts` — existing. Guards against x₀-sensitivity masquerading as performance (per the CLAUDE.md caveat).
3. `bounded-benchmarks.ts` — **new**. 6–10 problems where bounds are active (bounded Rosenbrock, half-bounded quadratic, fixed variables, Hock–Schittkowski HS3, etc.). Compares:
   - L-BFGS-B — native bounds via `settings.bounds`.
   - `LBFGS` / `GradientDescent` / `Adam` / `PSO` / `NelderMead` — bounds translated to `Constraint[]` via `boxConstraints()` and handled by the penalty layer.
   - Success = ‖x − x*‖∞ ≤ 1e-4 **and** all bounds respected exactly. A violation of any bound is failure, regardless of `f`.

Test-function problem records are added to `benchmarks/test-functions.ts`; runners only orchestrate.

### 0.7 Commit phasing (supersedes §12 test-plan sequencing)

| # | Commit | Content | ~ LoC (code + tests) |
|---|---|---|---|
| 1 | scaffold | Folder structure, `types.ts`, empty `LBFGSB` class, `withDefaults` with validation, registration | 150 + 100 |
| 2 | line-search | `dcsrch` + `dcstep` + unit tests vs Fortran tabulated outputs | 400 + 200 |
| 3 | bfgs-mat | Compact form, ring buffer, W-products, `solveM`, private `potrf`/`trsv` + unit tests (compact `B_k v` vs explicit `H v` on a quadratic) | 300 + 150 |
| 4 | bounds | `classifyBounds`, `project`, `projectedGradient`, `maxFeasibleStep` + unit tests | 150 + 100 |
| 5 | cauchy | Generalised Cauchy point + binary min-heap + unit tests (2D breakpoint sequence to 14 digits) | 300 + 150 |
| 6 | subspace | Subspace minimisation + Morales–Nocedal backtrack + unit tests (unconstrained-equivalent ≡ L-BFGS to 1e-12) | 250 + 150 |
| 7 | driver | `runInternal` + `runInternalAsync`, line-search retry, integration tests | 250 + 200 |
| 8 | benchmarks | Register in three runners; regenerate `.md` reports | 200 + 0 |
| 9 | docs | TSDoc, section README, top-level README, `CLAUDE.md` tree update (checklist §1–5) | docs |

Each commit: compiles, passes `npm run lint-fix && npm run build && npm test`, and is independently reviewable.

### 0.8 Override table — what §1–§13 + Appendices claim vs what we do

| Generic-spec section | Status in this library |
|---|---|
| §2.1 Module surface (root-level `lbfgsb.ts`, `cauchy.ts`, …) | **Replaced** → §0.5 (subfolder under `optimizers/`) |
| §2.2 Types (`FgFn`, `LBFGSBOptions`, `LBFGSBResult`, `ExitStatus`, `IterationInfo`) | **Replaced** → §0.2, §0.3 (library types) |
| §2.3 `minimize(fg, x0, bounds, options)` | **Replaced** → inherited `Optimizer.minimize(fn, x0, settings)` |
| §2.4 Defaults | **Replaced** → §0.3 table |
| §3 Compact representation | **Accepted** verbatim |
| §4 Outer loop | **Accepted** with §0.4 clarification on test order and best-so-far tracking |
| §5 Generalised Cauchy point | **Accepted** verbatim |
| §6 Subspace minimisation + Morales–Nocedal | **Accepted** verbatim |
| §7 Moré–Thuente line search | **Accepted** verbatim (no Armijo intermediate) |
| §8 Bounds API | **Narrowed**: internal `nbd`/`classifyBounds` kept; public surface is `{lower, upper}` object |
| §9 Numerical edge cases | **Accepted** + extended with cases C1–C13 (n=1, all-fixed, start-at-optimum, etc.) |
| §10 Linalg primitives | **Narrowed** → §0.5, §9 of decisions: only matrix ops, private inside `bfgs-mat.ts` |
| §11 Convergence criteria | **Replaced** → §0.4 (7-test table with `converged` mapping) |
| §12 Test plan | **Accepted** + extended with `bounded-benchmarks.ts` (§0.6) |
| §13.1–§13.12 TS notes | **Mostly accepted**; §13.2 "avoid `class`" overridden — public class required by `Optimizer<S>` |
| Appendix A Fortran workspace sizing | **Ignored** — not TS-relevant |
| Appendix B Exit-message strings | **Ignored** — `OptimizationResult` has no `message` field |
| Appendix C SciPy ↔ Fortran cross-reference | **Replaced** → §0.3 table (our names) |

---

## 1. Overview and background

### 1.1 Problem

L-BFGS-B minimises a smooth nonlinear objective subject to simple bounds:

$$\min_{x\in\mathbb{R}^n}\; f(x)\quad\text{s.t.}\quad l \le x \le u,\qquad l_i\in\mathbb{R}\cup\{-\infty\},\;\; u_i\in\mathbb{R}\cup\{+\infty\}.$$

`f : R^n → R` is assumed continuously differentiable; the gradient `g(x) = ∇f(x)` is supplied by the caller (analytic derivatives) or approximated by finite differences. Any variable may be unbounded (`±Infinity`), bounded on one side only, or fixed (`l_i = u_i`).

### 1.2 Origin and lineage

- **Byrd, Lu, Nocedal & Zhu (1995)**, *A Limited Memory Algorithm for Bound Constrained Optimization*, SIAM J. Sci. Comput. 16(5):1190–1208 — the original algorithm, introducing the generalized Cauchy point on the projected-gradient path and subspace minimisation over the free set.
- **Zhu, Byrd, Lu & Nocedal (1997)**, *Algorithm 778: L-BFGS-B*, ACM TOMS 23(4):550–560 — the reference Fortran 77 implementation with reverse communication.
- **Morales & Nocedal (2011)**, *Remark on Algorithm 778*, ACM TOMS 38(1):7 — two critical fixes: (i) a refinement of the subspace minimisation step that projects-and-backtracks rather than truncating along a straight line to the Cauchy point; (ii) replacement of the home-grown `dpmeps` machine-ε routine, which returned ≈10⁻¹⁴ on some Linux/gcc builds, with the Fortran intrinsic `epsilon(1.0d0)` (true 2.22×10⁻¹⁶). The two changes produce materially faster convergence on CUTEst problems and restore the stop-on-function-decrease criterion. The resulting source is distributed as **v3.0** (dated Feb 2011).
- **Byrd, Nocedal & Schnabel (1994)**, *Representations of Quasi-Newton Matrices*, Math. Programming 63:129–156 — supplies the **compact representation** `B_k = θI − W M Wᵀ` that makes the $O(mn)$ storage and $O(m^2)$-per-segment Cauchy search possible.
- **Moré & Thuente (1994)**, *Line Search Algorithms with Guaranteed Sufficient Decrease*, ACM TOMS 20(3):286–307 — the `dcsrch`/`dcstep` line search enforcing strong Wolfe conditions.

### 1.3 Relationship to BFGS and L-BFGS

- **BFGS** stores a dense $n\times n$ Hessian approximation ($O(n^2)$ memory). Suitable for $n$ up to a few thousand.
- **L-BFGS** stores only the last $m$ pairs $(s_i,y_i)$ and applies $H_k g$ implicitly via the **two-loop recursion** in $O(mn)$ flops; no bounds.
- **L-BFGS-B** retains $O(mn)$ storage but must (i) detect which variables are active (at a bound) and (ii) minimise the quadratic model over the free variables. Instead of the two-loop recursion it uses the **compact representation** of the same limited-memory matrix because the compact form allows $O(m^2)$ updates of directional derivatives along the piecewise-linear path inside the Cauchy-point routine — something the two-loop form cannot do.

Default memory size `m = 10`; common range 3–20. Memory cost: **2mn + 11m² + 5n + 8m** doubles (v3.0 workspace).

---

## 2. Public TypeScript API

### 2.1 Module surface

File layout (recommended):

```
src/
  index.ts          // re-exports
  lbfgsb.ts         // minimize() driver (mainlb analogue)
  cauchy.ts         // generalized Cauchy point
  subspace.ts       // subspace minimisation (subsm)
  lineSearch.ts     // dcsrch + dcstep
  bfgsMat.ts        // compact representation (W, M, ring buffer, formk, formt)
  linalg.ts         // dot, axpy, nrm2, scal, copy, trsv, potrf, trsm
  bounds.ts         // project(), projectedGradient(), maxFeasibleStep(), classifyBounds()
  types.ts          // public types
```

All modules use ES module syntax, named exports only (no `default` exports). The package is isomorphic (Node.js ≥ 18 and modern browsers; no `fs`/`process` imports in hot code).

### 2.2 Types (`types.ts`)

```ts
export type Vec = Float64Array;                        // dense vector
export type FgFn = (x: Vec, gOut: Vec) => number;      // writes gradient into gOut, returns f(x)

/** Bound-encoding matches Nocedal's nbd[i]: 0=free, 1=lower, 2=both, 3=upper. */
export const enum BoundType { Free = 0, Lower = 1, Both = 2, Upper = 3 }

export interface LBFGSBOptions {
  /** Memory size (number of correction pairs). Default 10. */
  m?: number;
  /** Relative function-value tolerance: stop when Δf ≤ factr · εmach. Default 1e7. */
  factr?: number;
  /** Projected gradient ∞-norm tolerance. Default 1e-5. */
  pgtol?: number;
  /** Maximum outer iterations. Default 15000. */
  maxIterations?: number;
  /** Maximum f/g evaluations (across line searches). Default 15000. */
  maxFunctionEvaluations?: number;
  /** Maximum Moré–Thuente iterations per outer step. Default 20. */
  maxLineSearch?: number;
  /** Line-search Armijo parameter (ftol). Default 1e-3 (matches L-BFGS-B v3.0). */
  ftol?: number;
  /** Line-search curvature parameter (gtol). Default 0.9. */
  gtol?: number;
  /** Line-search relative interval tolerance (xtol). Default 0.1. */
  xtol?: number;
  /** Optional verbose callback invoked after each accepted step. Return true to abort. */
  callback?: (info: IterationInfo) => boolean | void;
  /** Verbosity: 0 silent, 1 summary, ≥2 per-iteration detail. Default 0. */
  verbose?: 0 | 1 | 2;
}

export interface IterationInfo {
  readonly iter: number;
  readonly nfev: number;
  readonly f: number;
  readonly projGradInfNorm: number;
  readonly x: Vec;         // view into solver state; copy if retained
  readonly g: Vec;
}

export type Bounds =
  | { lower?: number | ArrayLike<number>; upper?: number | ArrayLike<number> }
  | ReadonlyArray<readonly [number | null, number | null]>;

export const enum ExitStatus {
  Converged = 0,           // projected-gradient or Δf test satisfied
  IterationLimit = 1,      // max iter / maxfun hit
  Abnormal = 2,            // line-search failure after retry
  InvalidInput = 3,        // thrown, never returned
}

export interface LBFGSBResult {
  x: Vec;                  // solution (guaranteed feasible)
  f: number;
  g: Vec;                  // gradient at solution
  projGradInfNorm: number;
  iterations: number;
  functionEvaluations: number;
  status: ExitStatus;
  converged: boolean;
  message: string;
  history?: ReadonlyArray<IterationInfo>; // only populated if verbose ≥ 2
}
```

### 2.3 Entry point (`lbfgsb.ts`)

```ts
export function minimize(
  fg: FgFn,
  x0: ArrayLike<number>,
  bounds?: Bounds,
  options?: LBFGSBOptions,
): LBFGSBResult;
```

`x0` may be a plain `number[]` or a `Float64Array`; it is copied internally. The solution in the result is always a fresh `Float64Array`. Unspecified bounds default to `(-Infinity, +Infinity)` per dimension.

**Why `Float64Array` everywhere internally.** It is contiguous, avoids boxed-number boxing in V8, gives deterministic bit-for-bit behaviour across platforms, and keeps hidden-class stability inside hot loops. Public API accepts `ArrayLike<number>` for ergonomics and converts once at the entry.

**Combined `fg` vs separate `f`/`g`.** The API is combined. This mirrors SciPy's `jac=True` mode, minimises redundant work when `f` and `g` share computation (the common case), and removes ambiguity about whether `x` may be mutated. A convenience overload that accepts `(f, grad)` separately can be layered on top.

### 2.4 Defaults (scipy-compatible)

| Option | Default | Source |
|---|---|---|
| `m` | 10 | SciPy `maxcor=10` |
| `factr` | 1e7 | classic; maps to `ftol = 1e7 · 2.22e-16 ≈ 2.22e-9` |
| `pgtol` | 1e-5 | SciPy `gtol=1e-5` |
| `maxIterations` | 15000 | SciPy |
| `maxFunctionEvaluations` | 15000 | SciPy |
| `maxLineSearch` | 20 | SciPy `maxls=20`; also Fortran `iback` |
| `ftol` (line search Armijo $c_1$) | 1e-3 | L-BFGS-B v3.0 `lnsrlb` |
| `gtol` (line search curvature $c_2$) | 0.9 | quasi-Newton standard |
| `xtol` | 0.1 | L-BFGS-B v3.0 `lnsrlb` |

---

## 3. Compact limited-memory BFGS representation

### 3.1 Stored data

For the $m$ most recent accepted pairs,

$$s_i = x_{i+1}-x_i,\qquad y_i = g_{i+1}-g_i,\qquad i=k-m,\ldots,k-1,$$

arranged as columns of

$$S_k=[s_{k-m},\ldots,s_{k-1}]\in\mathbb{R}^{n\times m},\qquad Y_k=[y_{k-m},\ldots,y_{k-1}]\in\mathbb{R}^{n\times m}.$$

### 3.2 Scalar $\theta_k$ (Shanno–Phua)

$$\boxed{\;\theta_k \;=\; \dfrac{y_{k-1}^{T}\,y_{k-1}}{s_{k-1}^{T}\,y_{k-1}}.\;}$$

$\theta_k I$ is the "initial" matrix to which the $m$ BFGS updates are applied. Recomputed every iteration after an accepted update. (Morales–Nocedal emphasise that consistently updating $\theta$ is critical.)

### 3.3 Compact form of $B_k$

$$\boxed{\;B_k \;=\; \theta_k I \;-\; W_k\,M_k\,W_k^{T},\;}\qquad W_k = [\,Y_k\;\;\theta_k S_k\,]\in\mathbb{R}^{n\times 2m},$$

where the $2m\times 2m$ matrix $M_k$ and its easily-assembled inverse $K_k := M_k^{-1}$ are

$$K_k \;=\; \begin{bmatrix}-D_k & L_k^{T}\\ L_k & \theta_k\,S_k^{T}S_k\end{bmatrix},\qquad M_k \;=\; K_k^{-1}.$$

The $m\times m$ blocks are derived element-wise from $s_i^T y_j$:

$$(L_k)_{ij}=\begin{cases}s_{k-m-1+i}^{T}y_{k-m-1+j}, & i>j,\\ 0, & i\le j,\end{cases}\qquad D_k=\mathrm{diag}\bigl(s_{k-m}^{T}y_{k-m},\ldots,s_{k-1}^{T}y_{k-1}\bigr).$$

In code, pre-compute the full Gram matrix $G := S_k^T Y_k \in \mathbb{R}^{m\times m}$; then $D_k = \mathrm{diag}(G)$ and $L_k = \mathrm{strictLowerTri}(G)$.

### 3.4 Useful auxiliary factorisation (`formt`)

Define

$$T_k \;:=\; \theta_k\,S_k^{T}S_k + L_k\,D_k^{-1}\,L_k^{T}\quad\in\mathbb{R}^{m\times m},$$

which is **symmetric positive-definite** whenever the stored pairs satisfy the curvature condition. Compute its upper-triangular Cholesky factor $T_k = J_k J_k^{T}$ (LINPACK `dpofa` in the reference). Together with $D_k^{1/2}$ this gives the block factorisation

$$K_k \;=\; \begin{bmatrix}D_k^{1/2} & 0\\ -L_k D_k^{-1/2} & J_k\end{bmatrix}\begin{bmatrix}-D_k^{1/2} & D_k^{-1/2}L_k^{T}\\ 0 & J_k^{T}\end{bmatrix},$$

which is what `cauchy` and `subsm` use to apply $M_k = K_k^{-1}$ to a $2m$-vector in $O(m^2)$.

### 3.5 Ring-buffer update (`bfgsMat.ts`, `matupd` analogue)

Maintain:

```ts
class BFGSMat {
  readonly n: number;
  readonly m: number;
  S: Float64Array;   // length m*n, column-major: column i at offset i*n
  Y: Float64Array;   // same layout
  G: Float64Array;   // m*m, G[i,j] = s_i^T y_j (full, not just lower)
  SS: Float64Array;  // m*m, SS[i,j] = s_i^T s_j (symmetric)
  Wt: Float64Array;  // m*m upper-triangular Cholesky factor J_k
  theta: number;
  col: number;       // current number of stored pairs 0..m
  head: number;      // ring-buffer index of oldest column 0..m-1
}
```

On accepted pair $(s_\text{new}, y_\text{new})$:

1. If `col < m`: new column is appended at `(head + col) mod m`; `col++`.
2. Else: overwrite column `head`; `head = (head + 1) mod m`.
3. Update `G`, `SS` in $O(mn)$ (just one new row/column against all stored vectors).
4. $\theta \leftarrow y_\text{new}^T y_\text{new} / s_\text{new}^T y_\text{new}$.
5. Recompute $T_k$ and Cholesky factor it into `Wt`.

**Never shift columns.** Indexing with `(head + j) mod m` for logical column $j$ (oldest) through `col-1` (newest) is $O(1)$ and allocation-free.

### 3.6 Product `B_k · v` and `Mv`

- $W_k^T v = \begin{bmatrix} Y_k^T v \\ \theta_k S_k^T v\end{bmatrix}$, $O(mn)$.
- Solving $K_k x = z$ (i.e. $x = M_k z$) via the block factorisation in §3.4, $O(m^2)$.
- $B_k v = \theta_k v - W_k (M_k (W_k^T v))$, total $O(mn + m^2)$.

---

## 4. Core L-BFGS-B outer loop

Numbered steps of `mainlb` (`lbfgsb.ts`):

```
1.  Input x0 (possibly infeasible), fg, bounds, options.
2.  Classify bounds → nbd[i] ∈ {0,1,2,3}; validate l ≤ u element-wise (throw on l>u).
3.  x ← project(x0, l, u).
4.  (f, g) ← fg(x, gbuf);  nfev ← 1.
5.  Initialise BFGSMat with col=0, theta=1, head=0.
6.  Compute pg := projectedGradient(x, g, l, u); if ||pg||_∞ ≤ pgtol → return Converged.
7.  Main loop, iter = 0, 1, …:
    a. (xc, c) ← cauchy(x, g, l, u, θ, S, Y, W-products, M-factors).
    b. Identify free set F = { i : l_i < xc_i < u_i }.
    c. If col > 0 and |F| > 0:
         d̂ ← subspaceMin(x, g, l, u, xc, c, θ, S, Y, F, M-factors)  // full-n direction
       else:
         d̂ ← xc − x.
    d. d ← d̂ (search direction from x).
    e. stpmax ← maxFeasibleStep(x, d, l, u).
    f. For iter == 0 only: α0 ← min(1, 1/||d||_∞); else α0 ← 1.
       Clip α0 ≤ stpmax.
    g. (stp, f_new, g_new, nLS) ← lineSearch(fg, x, d, f, g, α0, stpmax, …).
       nfev += nLS.
    h. x_new ← x + stp · d.   // same as x + stp · d̂; but for robustness, after line search compute x_new = project(x_new, l, u) to absorb float rounding at active bounds.
    i. Relative Δf convergence test:
         if |f − f_new| ≤ factr·εmach · max(|f|,|f_new|,1): return Converged.
    j. s = x_new − x;  y = g_new − g.
    k. dr = sᵀy;  if dr > εmach · (−gᵀd)·stp (curvature):
         push (s, y) into BFGSMat, update G, SS, θ, T_k = J J^T.
       else: nskip++.
    l. (x, f, g) = (x_new, f_new, g_new);  iter++.
    m. Projected-gradient convergence test on new g.
    n. Iteration/evaluation-limit check.
    o. Invoke callback(info); if callback returns true → return Abnormal/UserAbort.
8.  Return result with status, x, f, g, iterations, nfev, message.
```

**Line-search failure recovery (step g).** If `dcsrch` exits with an "abnormal termination" status:
1. If `col == 0` (memory already empty) and this is the first retry: return status = `Abnormal`.
2. Otherwise reset the memory (`col = 0; theta = 1`), reuse the current `(x, f, g)`, and restart from step 7a. One retry is allowed per iteration; if it fails again, return `Abnormal`.

---

## 5. Generalized Cauchy Point (`cauchy.ts`)

The Cauchy point $x^c$ is the first local minimiser of the quadratic model $m_k$ along the piecewise-linear projected-gradient path $x(t) = P(x_k - t g_k, l, u)$.

### 5.1 Breakpoints

For each variable $i$:

$$t_i \;=\; \begin{cases}(x_i - u_i)/g_i, & g_i<0\text{ and }u_i<+\infty,\\ (x_i - l_i)/g_i, & g_i>0\text{ and }l_i>-\infty,\\ +\infty, & \text{otherwise.}\end{cases}$$

Initial direction: `d_i = -g_i` for `t_i > 0` else `d_i = 0`.

### 5.2 Invariants along the path

On segment $[t^{(j-1)}, t^{(j)}]$, with $\Delta t = t - t^{(j-1)}$, let $z^{(j-1)} = x^{(j-1)} - x_k$. Maintain $2m$-vectors

$$p := W^{T} d,\qquad c := W^{T} z.$$

Then (eqs. (4.9)–(4.10) of BLNZ '95):

$$\boxed{\;f' \;=\; g^{T}d \;+\; \theta\,d^{T}z \;-\; p^{T}M\,c,\qquad f'' \;=\; \theta\,d^{T}d \;-\; p^{T}M\,p.\;}$$

The minimiser on the segment is $\Delta t^* = -f'/f''$.

### 5.3 Update when variable $b$ hits bound at $t^{(j)}$

Let $w_b \in \mathbb{R}^{2m}$ be the $b$-th row of $W$ (i.e. `[Y[b,:], θ S[b,:]]`), $g_b = g_b$, $z_b = x^c_b - x_{k,b}$ ($= l_b - x_{k,b}$ or $u_b - x_{k,b}$).

Advance to breakpoint:
$$c \leftarrow c + \Delta t^{(j-1)}\,p,\qquad z \leftarrow z + \Delta t^{(j-1)}\,d.$$

Then update (eqs. (4.11)–(4.12) of BLNZ '95):

$$\boxed{\begin{aligned}f'_{\rm new} &= f'_{\rm old} + \Delta t^{(j-1)}\,f''_{\rm old} + g_b^{2} + \theta\,g_b\,z_b - g_b\,w_b^{T}M\,c_{\rm new},\\ f''_{\rm new} &= f''_{\rm old} - \theta\,g_b^{2} - 2\,g_b\,w_b^{T}M\,p_{\rm old} - g_b^{2}\,w_b^{T}M\,w_b,\\ p_{\rm new} &= p_{\rm old} + g_b\,w_b,\\ d_b &\leftarrow 0. \end{aligned}}$$

Each inner step costs $O(m^2)$ (the three $w_b^T M v$ products).

### 5.4 Breakpoint ordering

Use a **binary min-heap** over $\{t_i\}$. Construction $O(n)$; each extraction $O(\log n)$. Only `n_int` breakpoints (≤ number of segments actually visited) are extracted, which is typically ≪ n.

### 5.5 Pseudocode (Algorithm CP)

```
function cauchy(x, g, l, u, nbd, θ, S, Y, G, SS, Wt):
  n ← length(x)
  d  ← new Float64Array(n)
  xc ← copy(x)
  tArr ← new Float64Array(n)         // breakpoints
  heapIdx ← []

  // (1) breakpoints, initial d, fill heap
  fPrime ← 0
  dTd    ← 0
  for i in 0..n-1:
    if nbd[i] ≠ Free and g[i] ≠ 0:
      if      g[i] < 0 and nbd[i] ∈ {Upper,Both}: tArr[i] = (x[i]-u[i])/g[i]
      else if g[i] > 0 and nbd[i] ∈ {Lower,Both}: tArr[i] = (x[i]-l[i])/g[i]
      else:                                       tArr[i] = +∞
    else: tArr[i] = +∞
    if tArr[i] > 0:
      d[i] = -g[i]
      fPrime -= g[i]*g[i]             // g^T d = -||d_free||^2
      dTd    += g[i]*g[i]
      heap.push(tArr[i], i)

  // (2) compute p = W^T d  (O(mn))
  p ← WT_times_vec(Y, θ, S, d, col, head)  // length 2m

  // (3) initial fPrime already = g^T d;  add θ·d^T z (z=0 here) and subtract p^T M p via formula
  // Initial z=0, c=0, so fPrime = g^T d, fDoublePrime = θ d^T d - p^T M p
  Mp ← solveM(Mfactors, p)             // Mp = M p
  fDoublePrime ← θ*dTd - dot(p, Mp)
  fDoublePrime ← max(fDoublePrime, εmach·|fDoublePrime|)   // guard tiny +/- 0

  c ← zeros(2m)
  dtMin ← -fPrime / fDoublePrime       // candidate interior minimum
  tOld  ← 0
  if heap empty: goto FINALIZE
  (tCurrent, b) ← heap.pop()
  Δt ← tCurrent - tOld

  // (4) main sweep
  while dtMin ≥ Δt and heap has elements:
    xc[b] ← (d[b] > 0) ? u[b] : l[b]    // snap to bound
    zb    ← xc[b] - x[b]
    tOld  ← tCurrent
    // advance c, then fPrime, fDoublePrime, then p
    axpy(c, Δt, p)                       // c += Δt · p
    wb ← row_b_of_W(Y, θ, S, b)          // 2m-vector
    Mc ← solveM(Mfactors, c)
    Mp_old ← Mp                          // keep old
    Mwb ← solveM(Mfactors, wb)
    fPrime       += Δt*fDoublePrime + g[b]*g[b] + θ*g[b]*zb - g[b]*dot(wb, Mc)
    fDoublePrime += -θ*g[b]*g[b] - 2*g[b]*dot(wb, Mp_old) - g[b]*g[b]*dot(wb, Mwb)
    fDoublePrime = max(fDoublePrime, εmach*|fDoublePrime|)
    axpy(p, g[b], wb)                     // p += g[b] · wb
    Mp = solveM(Mfactors, p)              // refresh for next iter
    d[b] = 0
    dtMin = -fPrime / fDoublePrime
    if heap empty: break
    (tCurrent, b) = heap.pop()
    Δt = tCurrent - tOld
    if fPrime ≥ 0:                        // derivative turned positive → done
      dtMin = 0
      break

  // (5) FINALIZE
  dtMin = max(dtMin, 0)
  tOld += dtMin
  for i in 0..n-1:
    if tArr[i] ≥ tCurrent and nbd[i] ≠ Free:   // variable still free at end
      xc[i] = x[i] + tOld * d[i]
    // else xc[i] already set to bound or remains x[i] if d[i]==0
  axpy(c, dtMin, p)                        // final c = W^T (xc - x)
  return (xc, c, freeSet = { i : l_i < xc[i] < u_i })
```

Notes:
- `solveM(Mfactors, v)` solves $K_k\,x = v$ using the block factorisation of §3.4; cost $O(m^2)$.
- Treat `ties` in breakpoints by processing all equal-$t$ indices together (extract all heap items with the same $t$ before updating $p$); this preserves numerical equivalence with the Fortran for degenerate bounds.
- Guard for `g[i] == 0 exactly` and for $\pm\infty$ bounds — no breakpoint is generated.

---

## 6. Subspace minimisation (`subspace.ts`)

### 6.1 Formulation

Let $F = \{i : l_i < x^c_i < u_i\}$, $t = |F|$, and $Z \in \mathbb{R}^{n\times t}$ select the free coordinates (columns are unit vectors $e_i$ for $i \in F$). Solve approximately

$$\min_{\bar d\in\mathbb{R}^t} \bar r^{\,c\,T}\bar d + \tfrac12 \bar d^{T} \bar B \bar d,\qquad \bar B = \theta I_t - N^{T} M_k N,\;\; N := W^{T} Z\in\mathbb{R}^{2m\times t},$$

with the **reduced gradient at $x^c$** (BLNZ eq. (5.4)):

$$\boxed{\;\bar r^{\,c} = Z^{T}\bigl[g_k + \theta(x^c-x_k) - W_k M_k c\bigr].\;}$$

### 6.2 Direct primal method via Sherman–Morrison–Woodbury

$$\boxed{\;\bar B^{-1} = \tfrac{1}{\theta}I_t + \tfrac{1}{\theta^2} N^{T}\bigl(I_{2m} - \tfrac{1}{\theta} M_k N N^{T}\bigr)^{-1} M_k N.\;}$$

Unconstrained subspace minimiser:

$$\boxed{\;\bar d^{\,u} = -\tfrac{1}{\theta}\bar r^{c} - \tfrac{1}{\theta^2} N^{T}\bigl(I - \tfrac{1}{\theta} M N N^{T}\bigr)^{-1} M N \bar r^{c}.\;}$$

Form $N\bar r^c$ ($2m$-vector, $O(mn)$ since $N = W^T Z$), compute $M N \bar r^c$ ($O(m^2)$), assemble the $2m\times 2m$ matrix $I - (1/\theta) M N N^T$, and Cholesky-solve it (fall back to LDLᵀ via Bunch–Kaufman if indefinite).

### 6.3 Morales–Nocedal 2011 refinement

After computing $\bar d^u$:

1. Let $\hat x = x^c + Z\bar d^u$ (extended to full $n$ by zeros on non-free).
2. If $l \le \hat x \le u$: accept $\hat x$ directly.
3. Otherwise **project + backtrack**:
   ```
   λ = 1
   for step in 0..maxBacktrack (≈ 20):
     x_trial = project(x^c + λ·Z·bar_du, l, u)
     if m_k(x_trial) < m_k(x^c): accept and break
     λ /= 2
   ```
   where $m_k(x) = f_k + g_k^T(x-x_k) + \tfrac12(x-x_k)^T B_k (x-x_k)$.
4. In the worst case (no decrease after backtracks) set $\hat x = x^c$ (fall back to Cauchy).

This is the `c-jlm-jn` fix in `subsm.f`: it replaces the 1997 code's "truncate along the straight line from $x^c$ to $\hat x_u$ at the first bound hit" behaviour, which produced poor directions when many bounds were crossed simultaneously.

### 6.4 Pseudocode

```
function subspaceMin(x, g, l, u, xc, c, θ, S, Y, G, SS, Wt, freeSet):
  t = |freeSet|
  if t == 0 or col == 0: return xc   // nothing to do

  // (1) reduced gradient r̄^c
  Mc   = solveM(Mfactors, c)              // 2m
  W_Mc = W_times_vec(Y, θ, S, Mc)         // full-n
  r    = new Float64Array(t)
  for k,i in freeSet.enumerate():
    r[k] = g[i] + θ*(xc[i] - x[i]) - W_Mc[i]

  // (2) N = W^T Z (only rows of W corresponding to freeSet)
  //     in practice compute N·r̄^c directly: a 2m-vector
  //     Nr[2m] = sum over freeSet of W[i,:] * r[k]
  Nr = new Float64Array(2*col)
  for k,i in freeSet.enumerate():
    for j in 0..col-1:
      Nr[j]       += Y[i, j] * r[k]        // Y block
      Nr[col + j] += θ * S[i, j] * r[k]    // θS block

  // (3) form (I − M N N^T / θ) — need M N N^T as 2m×2m
  //     A_ij = sum_k W[i,k] W[j,k] = (W^T W) at free rows
  A = zeros(2m, 2m)                         // = N N^T
  for k,i in freeSet.enumerate():
    for a in 0..2m-1:
      for b in 0..2m-1:
        A[a,b] += W[i,a] * W[i,b]           // O(t · 4m^2) naive; can be done incrementally
  // Better: maintain Y^T Z Z^T Y, S^T Z Z^T Y, S^T Z Z^T S via complement trick.

  MA = solveM_batch(Mfactors, A)            // column-wise: each solve O(m^2)
  MN_r = solveM(Mfactors, Nr)               // 2m
  Kmat = I_{2m} - (1/θ) * MA                 // 2m × 2m dense
  v    = Kmat \ MN_r                         // Cholesky or LDL^T (2m)^3/3

  // (4) unconstrained subspace step
  du = new Float64Array(t)
  // du = -r/θ - (1/θ^2) N^T v
  //    (N^T v)[k] = sum_j W[i_k, j] v[j]
  for k,i in freeSet.enumerate():
    NTv_k = 0
    for j in 0..2*col-1:
      NTv_k += (j < col ? Y[i,j] : θ*S[i, j-col]) * v[j]
    du[k] = -r[k]/θ - NTv_k/(θ*θ)

  // (5) Morales-Nocedal refinement
  x_hat = copy(xc)
  feasible = true
  for k,i in freeSet.enumerate():
    xi = xc[i] + du[k]
    if xi < l[i] or xi > u[i]: feasible = false
    x_hat[i] = xi

  if feasible:
    return x_hat

  // project + backtrack on m_k
  mk_xc = 0   // m_k(xc) = (1/2) c^T M c − θ·0  relative to xc; compare differences
  λ = 1
  for step in 0..20:
    for k,i in freeSet.enumerate():
      v_i = xc[i] + λ*du[k]
      x_hat[i] = clamp(v_i, l[i], u[i])
    Δm = evaluateModelChange(x_hat, xc, g, θ, W, M, c)
    if Δm < 0: return x_hat
    λ /= 2
  return xc    // fall back
```

`evaluateModelChange(x_hat, xc, g, …)` computes $m_k(x_\text{hat}) - m_k(x^c)$ using the compact form ($O(mn)$); only the difference is required.

---

## 7. Line search — Moré–Thuente (`lineSearch.ts`)

Implement the MINPACK-2 reverse-communication routine `dcsrch` with its safeguarded cubic/quadratic step routine `dcstep`. The function is encapsulated so the caller does not see reverse communication; internally a `LineSearchState` object plays the role of `isave`/`dsave`.

### 7.1 Strong Wolfe conditions

Directional: $\phi(\alpha) = f(x + \alpha d)$, $\phi'(\alpha) = \nabla f(x + \alpha d)^T d$. Given $0 < \texttt{ftol} < \texttt{gtol} < 1$:

$$\phi(\alpha) \le \phi(0) + \texttt{ftol}\cdot\alpha\cdot\phi'(0),\qquad |\phi'(\alpha)| \le \texttt{gtol}\cdot|\phi'(0)|.$$

Auxiliary function used in stage 1:
$$\psi(\alpha) = \phi(\alpha) - \phi(0) - \texttt{ftol}\cdot\alpha\cdot\phi'(0),\qquad \psi'(\alpha) = \phi'(\alpha) - \texttt{ftol}\cdot\phi'(0).$$
$\psi(0)=0$, $\psi'(0)<0$, sufficient decrease $\Leftrightarrow \psi(\alpha)\le 0$.

### 7.2 `LineSearchState` persistent fields

```ts
interface LineSearchState {
  brackt: boolean;
  stage: 1 | 2;
  ginit: number;   // φ'(0)
  gtest: number;   // ftol · ginit (< 0)
  finit: number;   // φ(0)
  fx: number;      // best f so far (stage-2 uses φ, stage-1 uses ψ but stored as φ)
  gx: number;
  fy: number;
  gy: number;
  stx: number;     // best step
  sty: number;     // other endpoint of uncertainty interval
  stmin: number;
  stmax: number;
  width: number;
  width1: number;
}
```

Constants: `p5 = 0.5`, `p66 = 0.66`, `xtrapl = 1.1`, `xtrapu = 4.0`.

### 7.3 `dcsrch` pseudocode (driver)

```
function dcsrch(stp, f, g, state, ftol, gtol, xtol, stpmin, stpmax, task):
  if task == 'START':
    validate (stp in [stpmin, stpmax], g < 0, ftol,gtol,xtol ≥ 0, stpmin ≥ 0, stpmax ≥ stpmin)
    state.brackt = false; state.stage = 1
    state.finit = f; state.ginit = g; state.gtest = ftol * g
    state.width  = stpmax - stpmin
    state.width1 = 2 * state.width
    state.stx = 0; state.fx = f; state.gx = g
    state.sty = 0; state.fy = f; state.gy = g
    state.stmin = 0; state.stmax = stp + xtrapu*stp
    return { stp, task: 'FG' }

  ftest = state.finit + stp * state.gtest
  if state.stage == 1 and f ≤ ftest and g ≥ 0: state.stage = 2

  // --- convergence & warning tests ---
  if state.brackt and (stp ≤ state.stmin or stp ≥ state.stmax):
      return { stp: state.stx, task: 'WARNING_ROUNDING' }
  if state.brackt and state.stmax - state.stmin ≤ xtol * state.stmax:
      return { stp: state.stx, task: 'WARNING_XTOL' }
  if stp == stpmax and f ≤ ftest and g ≤ state.gtest:
      return { stp, task: 'WARNING_STP_EQ_STPMAX' }
  if stp == stpmin and (f > ftest or g ≥ state.gtest):
      return { stp, task: 'WARNING_STP_EQ_STPMIN' }
  if f ≤ ftest and |g| ≤ gtol * (-state.ginit):
      return { stp, task: 'CONVERGENCE' }

  // --- step using ψ (stage 1) or φ (stage 2) ---
  if state.stage == 1 and f ≤ state.fx and f > ftest:
    // modify fields by subtracting linear term
    fm  = f  - stp*state.gtest
    fxm = state.fx - state.stx*state.gtest
    fym = state.fy - state.sty*state.gtest
    gm  = g  - state.gtest
    gxm = state.gx - state.gtest
    gym = state.gy - state.gtest
    (stx, fxm, gxm, sty, fym, gym, stp, brackt) =
      dcstep(state.stx, fxm, gxm, state.sty, fym, gym, stp, fm, gm,
             state.brackt, state.stmin, state.stmax)
    // restore φ values
    state.fx = fxm + stx*state.gtest
    state.fy = fym + sty*state.gtest
    state.gx = gxm + state.gtest
    state.gy = gym + state.gtest
    state.stx = stx; state.sty = sty; state.brackt = brackt
  else:
    (state.stx, state.fx, state.gx, state.sty, state.fy, state.gy, stp, state.brackt) =
      dcstep(state.stx, state.fx, state.gx, state.sty, state.fy, state.gy, stp, f, g,
             state.brackt, state.stmin, state.stmax)

  // --- bisection safeguard (2/3 rule) ---
  if state.brackt:
    if |state.sty - state.stx| ≥ p66 * state.width1:
      stp = state.stx + p5*(state.sty - state.stx)
    state.width1 = state.width
    state.width  = |state.sty - state.stx|

  // --- update admissible [stmin, stmax] ---
  if state.brackt:
    state.stmin = min(state.stx, state.sty)
    state.stmax = max(state.stx, state.sty)
  else:
    state.stmin = stp + xtrapl*(stp - state.stx)
    state.stmax = stp + xtrapu*(stp - state.stx)

  // --- clamp and fallback ---
  stp = clamp(stp, stpmin, stpmax)
  if (state.brackt and (stp ≤ state.stmin or stp ≥ state.stmax))
     or (state.brackt and state.stmax - state.stmin ≤ xtol * state.stmax):
    stp = state.stx
  return { stp, task: 'FG' }
```

### 7.4 `dcstep` — safeguarded cubic/quadratic step

Four cases keyed on `fp`, `dp`, `fx`, `dx` (with `sgnd = dp · sign(dx)`):

**Hermite cubic minimiser between points $(a, f_a, f'_a)$ and $(b, f_b, f'_b)$:**
$$\theta = 3(f_a - f_b)/(b-a) + f'_a + f'_b,\quad s = \max(|\theta|,|f'_a|,|f'_b|),\quad \gamma = s\sqrt{(\theta/s)^2 - (f'_a/s)(f'_b/s)}.$$
Sign of $\gamma$ is chosen so the minimiser points into $(a,b)$. Minimiser $\text{stpc} = a + \dfrac{(\gamma - f'_a) + \theta}{2\gamma + f'_b - f'_a}(b-a)$.

**Case 1 — `fp > fx` (higher f ⇒ bracket located):** cubic on $(\text{stx}, \text{stp})$; quadratic using $f_x, f'_x, f_p$. Pick cubic if closer to `stx`, else cubic+quadratic average. Set `brackt = true`. `bound = true`.

**Case 2 — `fp ≤ fx` and `sgnd < 0` (derivatives flipped sign ⇒ bracketed):** cubic; secant from `dp, dx`. Pick **farther** from `stp`. `brackt = true`. `bound = false`.

**Case 3 — `fp ≤ fx`, same-sign derivs, `|dp| < |dx|`:** cubic (guard $\gamma$ with `max(0, ...)` because cubic might have no real minimum); secant. If bracketed pick **closer** to `stp`, else **farther**. `bound = true`. `brackt` unchanged.

**Case 4 — `fp ≤ fx`, same-sign, `|dp| ≥ |dx|` (slope not decreasing):** if bracketed, cubic through $(\text{sty}, f_y, d_y)$ and $(\text{stp}, f_p, d_p)$; else extrapolate to `stpmax` or `stpmin` depending on sign. `bound = false`.

**After cases — update interval endpoints:**
```
if fp > fx:                     // case 1
  sty, fy, dy = stp, fp, dp
else:
  if sgnd < 0:                   // case 2
    sty, fy, dy = stx, fx, dx    // swap first
  stx, fx, dx = stp, fp, dp      // new best
```

**Final clamp + safeguard:**
```
stp = clamp(stpf, stpmin, stpmax)
if brackt and bound:            // cases 1 and 3 only
  if sty > stx: stp = min(stx + 0.66*(sty-stx), stp)
  else:         stp = max(stx + 0.66*(sty-stx), stp)
```

### 7.5 Outer wrapper

```ts
function lineSearch(fg, x, d, f0, g0, alpha0, stpmax, maxLS, ftol, gtol, xtol):
  stp = alpha0
  dirDeriv0 = dot(g0, d)
  if dirDeriv0 >= 0: throw 'descent direction required'
  state = newLineSearchState()
  task = 'START'
  f = f0; dg = dirDeriv0
  xTrial = new Float64Array(n); gTrial = new Float64Array(n)
  for k in 0..maxLS:
    ({stp, task}) = dcsrch(stp, f, dg, state, ftol, gtol, xtol, 1e-20, stpmax, task)
    if task == 'FG':
      axpyInto(xTrial, x, stp, d)     // xTrial = x + stp*d
      f  = fg(xTrial, gTrial)
      dg = dot(gTrial, d)
      continue
    if task == 'CONVERGENCE' or task starts with 'WARNING':
      return { stp, f, g: gTrial, nfev: k+1, ok: task != 'WARNING_ROUNDING' }
    return { stp, f, g: gTrial, nfev: k+1, ok: false }    // ERROR
  return { stp, f, g: gTrial, nfev: maxLS, ok: false }
```

### 7.6 Initial step

- First outer iteration: `α0 = min(1, 1/||d||_∞)` (or `1/||g||` equivalently since `d = -g`).
- Subsequent iterations: `α0 = 1`.
- In both cases, clamp to `stpmax = maxFeasibleStep(x, d, l, u)` (see §8).

### 7.7 Alternatives

**Hager–Zhang** (approximate-Wolfe) is more robust when function values have noise near machine precision, but substantially more code and not what L-BFGS-B v3.0 ships. **Simple Armijo backtracking** is too weak — the BFGS update needs `sᵀy > 0`, which the curvature condition guarantees. **Stick with Moré–Thuente**; it is what SciPy, L-BFGS-B-C, and LBFGSpp default to and is required for byte-level agreement with SciPy output.

---

## 8. Bounds handling (`bounds.ts`)

```ts
export function classifyBounds(l: Vec, u: Vec): Uint8Array;
// returns nbd[i]: 0=free, 1=lower, 2=both, 3=upper

export function project(x: Vec, l: Vec, u: Vec, nbd: Uint8Array, out: Vec): void;
// xOut[i] = median(l[i], x[i], u[i]) (respecting ±Infinity via nbd)

export function projectedGradient(
  x: Vec, g: Vec, l: Vec, u: Vec, nbd: Uint8Array, out: Vec
): number;
// Computes p[i] = project(x - g, l, u)[i] - x[i] and returns ||p||_∞.
// Equivalent to the L-BFGS-B Fortran `projgr`.

export function maxFeasibleStep(
  x: Vec, d: Vec, l: Vec, u: Vec, nbd: Uint8Array
): number;
// stpmax = min over i of:
//   (u[i]-x[i])/d[i]  if d[i] > 0 and u[i] < +∞
//   (l[i]-x[i])/d[i]  if d[i] < 0 and l[i] > -∞
//   +∞ otherwise
// Default Infinity if no d[i] hits a bound.
```

The `nbd` encoding matches Nocedal's exactly and lets callers leave garbage in `l[i]` or `u[i]` when the corresponding side is infinite.

---

## 9. Numerical edge cases and robustness

**Infinite bounds.** Use `Number.POSITIVE_INFINITY`/`NEGATIVE_INFINITY`. `classifyBounds` uses `isFinite()` per side; breakpoint computation and projection skip the infinite side cleanly.

**Fixed variables (`l_i == u_i`).** Validate $l_i \le u_i$; equality is allowed. After the initial project, $x_i = l_i = u_i$. Breakpoint is 0 (variable is already active); it never joins the free set. The algorithm handles this automatically.

**Invalid bounds (`l_i > u_i`).** Throw `new Error('bound[i]: lower > upper')` before starting.

**Initial point infeasible.** Always `project(x0, l, u)` at step 3 of the driver. Many users pass `x0 = zeros` with positive-only bounds; silent projection is the scipy behaviour.

**Curvature-condition failure.** Implement as
```ts
const dr = stp * (dot(gNew, d) - dot(gOld, d));         // = sᵀy
const ddum = -stp * dot(gOld, d);                        // = -gOldᵀs > 0 by Wolfe
if (dr <= Number.EPSILON * ddum) { nskip++; /* do not update S,Y */ }
```
This is scale-invariant and matches the Fortran code exactly. **Do not** implement Powell damping — classical L-BFGS-B does not use it; it would diverge from SciPy.

**Ill-conditioning of M / T.** If `dpofa`-equivalent Cholesky on `T_k` fails (negative or tiny pivot), evict the oldest pair and retry. If all pairs have been evicted, restart with `col=0, theta=1` (pure steepest descent for one iteration).

**Very small steps.** `stpmin = 1e-20` is the MINPACK default. `stpmax = 1e20` default, overridden by `maxFeasibleStep` in bounded mode.

**NaN / Infinity in `f` or `g`.** If the user objective returns `NaN` at any trial point:
1. Inside the line search: treat as `f = +Infinity`, let `dcstep` back off — this works naturally because `f > ftest`.
2. At `x0`: throw.
3. Elsewhere: return status `Abnormal` with message `"Objective returned NaN"`.
Validate with `Number.isFinite(f) && gradientIsFinite(g)` after every `fg` call.

**Zero gradient at start.** `projGradInfNorm(x0, g0) == 0` trivially satisfies the test; return Converged with `iterations = 0`.

**Line-search failure recovery.** On `ok=false`: (1) reset `col = 0; theta = 1` (wipe memory); (2) recompute `d = -g` (steepest descent); (3) retry line search with `α0 = 1/||g||`. If that also fails, return `Abnormal` with the original trial `x, f, g`.

**Rounding at active bounds.** Post-step, call `project(x, l, u)` once before the convergence test — this absorbs up to ~$\varepsilon|x|$ boundary violations introduced by `x + stp*d` (a known scipy issue #10919).

---

## 10. Data structures and linear-algebra primitives (`linalg.ts`)

Required primitives, all taking `Float64Array` buffers and explicit length/stride parameters:

```ts
export function dot(a: Vec, b: Vec, n: number): number;               // aᵀb
export function nrm2(a: Vec, n: number): number;                       // ||a||₂
export function infNorm(a: Vec, n: number): number;                    // ||a||∞
export function scal(alpha: number, a: Vec, n: number): void;          // a ← α·a
export function axpy(alpha: number, x: Vec, y: Vec, n: number): void;  // y ← αx + y
export function copy(src: Vec, dst: Vec, n: number): void;             // dst ← src
export function gemv(
  trans: boolean, m: number, n: number, alpha: number,
  A: Vec, lda: number, x: Vec, beta: number, y: Vec
): void;                                                                // y ← α·op(A)·x + β·y
export function potrf(U: Vec, n: number): boolean;                     // upper Cholesky in place
export function trsv(
  upper: boolean, trans: boolean, n: number, A: Vec, lda: number, x: Vec
): void;                                                                // solve triangular
export function trsm(...): void;                                       // batched triangular solve
```

`potrf` implements the standard column-Cholesky in-place on the upper triangle; returns `false` on non-positive-definite to let callers fall back to Bunch–Kaufman. For the 2m×2m solves, `m` is small (≤ 20), so naive $O(m^3)$ unblocked algorithms are optimal — do not pull in a BLAS.

**Matrix layout.** All matrices are **column-major** `Float64Array` buffers. Column $j$ of an $M\times N$ matrix occupies positions `[j*lda, j*lda+M)`, with `lda ≥ M`. This matches the Fortran reference and makes it easy to cross-reference `lbfgsb.f`. For the $S$, $Y$ ring buffers specifically, each logical column $j$ of the ring maps to physical column `(head + j) mod m`.

**Preallocation.** All working buffers (`xc`, `c`, `p`, `Mp`, `gNew`, `d`, `wn`, `wt`, heap storage, line-search `xTrial`/`gTrial`) are allocated once inside `minimize()` and reused. No allocations occur inside hot loops.

---

## 11. Convergence criteria

Three termination tests; each with a distinct status code and message:

1. **Projected-gradient (primary).**
   $$\|\,P(x-g, l, u) - x\,\|_{\infty} \;\le\; \texttt{pgtol}.$$
   Status 0, message `"CONVERGENCE: NORM OF PROJECTED GRADIENT ≤ PGTOL"`.

2. **Relative function decrease (secondary).**
   $$\frac{|f_k - f_{k+1}|}{\max(|f_k|, |f_{k+1}|, 1)} \;\le\; \texttt{factr} \cdot \varepsilon_{\rm mach}.$$
   With default `factr = 1e7`, threshold ≈ `2.22e-9`. Status 0, message `"CONVERGENCE: REL_REDUCTION_OF_F ≤ FACTR*EPSMCH"`.

3. **Iteration / evaluation limit.** Status 1, message `"STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT"` or `"…FUNCTION EVALUATIONS…"`.

**Abnormal exits (status 2):**
- `"ABNORMAL_TERMINATION_IN_LNSRCH"` — line search failed twice consecutively.
- `"NaN_OR_INF_IN_OBJECTIVE"` — objective returned non-finite mid-optimisation.

Treat `converged = (status === 0)`.

---

## 12. Test plan

### 12.1 Unit tests per module

- **`linalg`.** Property tests for `dot`, `axpy`, `gemv` against hand-rolled references; Cholesky against random SPD matrices of size 2–40.
- **`bounds`.** `project` idempotent; `maxFeasibleStep` returns `Infinity` when no bound; returns the exact boundary distance on synthetic inputs.
- **`bfgsMat`.** After $k$ updates on a quadratic $f(x) = \tfrac12 x^T H x$, the compact product $B_k v$ matches the explicit $H v$ up to $10^{-10}$ for $k \ge n$.
- **`cauchy`.** On a 2D problem with known breakpoints, verify the final `xc` and `c = Wᵀ(xc−x)` to 14 digits against a symbolic computation.
- **`subspaceMin`.** Verify that on an unconstrained-equivalent problem (all bounds ±∞), the direction reduces to the L-BFGS two-loop result within `1e-12`.
- **`dcstep`.** Regression test against tabulated outputs from the Fortran `dcstep` on 100 random inputs per case.
- **`dcsrch`.** Verify Wolfe-condition satisfaction on a suite of 1-D functions (convex quadratic, cubic, periodic, noisy quadratic) from Moré–Thuente §5.

### 12.2 Standard benchmark problems

| Problem | Dim | Bounds | Known optimum |
|---|---|---|---|
| Rosenbrock | 2, 10, 100 | ±∞ | $x^* = (1,\ldots,1)$, $f^* = 0$ |
| Beale | 2 | ±∞ | $(3, 0.5)$, $f^* = 0$ |
| Booth | 2 | ±∞ | $(1, 3)$, $f^* = 0$ |
| Himmelblau | 2 | ±∞ | four minima, $f^* = 0$ |
| Extended Rosenbrock | 10, 100 | ±∞ | all-ones, $f^* = 0$ |
| Extended Powell | 16 | ±∞ | zero, $f^* = 0$ |
| Quadratic $\tfrac12 x^T A x - b^T x$ | 50 | ±∞ | $A^{-1} b$, f closed form |
| Bounded Rosenbrock | 2 | $[-1.5, 0.5]^2$ | boundary active at solution |
| Fixed variables | 5 | $l_i = u_i$ for $i \in \{2,4\}$ | reduces to 3-D problem |
| Half-bounded | 4 | $l = 0$, $u = +\infty$ | tests `nbd = 1` path |

### 12.3 SciPy parity tests

Run `scipy.optimize.minimize(method='L-BFGS-B')` with identical options on each benchmark, with seed-fixed random starts. Assert:
- Final $f$ agrees to `1e-8` absolute or $|f|\cdot 10^{-10}$ relative, whichever larger.
- Final $x$ agrees in each coordinate to `1e-6`.
- Iteration count within ±10% (exact parity not guaranteed due to FP non-associativity).
- On Rosenbrock-2D from $(-1.2, 1)$ with defaults, expect ≤ 25 iterations, ≤ 30 fevals.

### 12.4 Edge-case tests

- Start at optimum → `iterations = 0`, status 0.
- `n = 1` problem.
- `n = 0` problem → return immediately with empty result.
- All variables fixed (`l = u` everywhere) → `iterations = 0`, trivial projection.
- Unbounded, constant gradient (linear objective) → returns `IterationLimit` since projected gradient never vanishes.
- `fg` throws → propagates; solver is idempotent if caller retries.
- `fg` returns NaN at `x0` → throw with clear message.

### 12.5 Property-based tests

Using `fast-check`:
- **Feasibility invariant.** For any valid `bounds`, any random `x0`, after any number of iterations, `x` is always in `[l, u]` exactly.
- **Monotone descent.** For convex problems, `f_{k+1} ≤ f_k` holds at every accepted step (line search guarantees it via Armijo).
- **Scaling invariance.** Scaling `x` by a positive diagonal, and adjusting bounds accordingly, should reach a solution that rescales the same way (within tolerance) — tests for hidden dependence on absolute magnitudes.

### 12.6 Performance benchmarks

Report iter/sec on Rosenbrock-100D and Rosenbrock-1000D across Node 20+ LTS. Target: within 3× of SciPy wall-clock for $n \le 1000$. Primary optimisation targets for TS: (a) avoiding object allocations in the hot loop, (b) fusing `axpy` loops, (c) specialising the 2m×2m solve for m ≤ 20.

---

## 13. TypeScript-specific implementation notes

1. **Strict mode.** `"strict": true`, `"noUncheckedIndexedAccess": true`, `"noImplicitAny": true`, `"exactOptionalPropertyTypes": true` in `tsconfig.json`. No `any` in the public API.

2. **Functional organisation.** The solver is a pure function `minimize(...)` that owns a private `SolverState` object. No singleton state, no classes exposed to the user. `BFGSMat` is a namespace of pure functions operating on a plain-data state record — avoid `class` to keep V8 hidden classes stable and to simplify cloning for testing.

3. **Determinism.** No `Math.random` anywhere. All tie-breaks explicit. All internal buffers zero-initialised.

4. **No recursion in hot paths.** The two-loop recursion is a plain for-loop; `dcsrch` is iterative (state machine on `task`); Cauchy-point sweep is iterative.

5. **Heap implementation.** Provide an inline min-heap keyed by `Float64Array` of breakpoint values plus an `Int32Array` of indices. 40 lines of code; no external dep.

6. **`Float64Array` vs `number[]`.** Use `Float64Array` everywhere internal. Accept `ArrayLike<number>` on the boundary and `new Float64Array(...)` once. This gives: (i) contiguous memory; (ii) no NaN-tagged double issues; (iii) compatibility with WebGL/WebGPU adjacent code; (iv) V8 inlines arithmetic on typed arrays more aggressively. Cost: a one-time $O(n)$ copy at the boundary.

7. **ES modules.** Named exports only. An `index.ts` re-exports `minimize`, types, and `ExitStatus`. No default export to keep tree-shaking friendly and to avoid dual-package hazard.

8. **Isomorphism.** Solver itself uses no Node or browser globals. Tests should run in both environments (`vitest` supports both).

9. **Error reporting.** Programmer errors (invalid bounds, non-finite `x0`, `ftol ≥ gtol`) throw synchronously. Algorithmic non-convergence returns a result with `status !== 0`.

10. **Numerical helpers.** Implement the 10 primitives listed in §10 from scratch. The `linalg` module totals ~200 lines. Do **not** pull in `ndarray`, `gl-matrix`, or `scipy-like` packages — they add a large surface and their memory layouts differ. Keep a clean contract.

11. **Style.** Use `const` aggressively, prefer explicit narrow integer types (`Uint8Array` for `nbd`, `Int32Array` for heap indices). Annotate every loop bound as `const` before the loop. Use bracket-access on typed arrays (`a[i]`) rather than `.at()` for speed.

12. **Testing harness.** Prefer `vitest` for TS-native tests, `fast-check` for property tests, `benchmark` or `mitata` for micro-benchmarks. Cross-validate against SciPy by emitting JSON trajectories from Python and asserting TS matches.

---

## 14. References

- **Byrd, R. H.; Lu, P.; Nocedal, J.; Zhu, C.** (1995). *A Limited Memory Algorithm for Bound Constrained Optimization*. SIAM J. Sci. Comput. 16(5):1190–1208. Preprint: `users.iems.northwestern.edu/~nocedal/PDFfiles/limited.pdf`.
- **Zhu, C.; Byrd, R. H.; Lu, P.; Nocedal, J.** (1997). *Algorithm 778: L-BFGS-B*. ACM TOMS 23(4):550–560.
- **Morales, J. L.; Nocedal, J.** (2011). *Remark on "Algorithm 778: L-BFGS-B…"*. ACM TOMS 38(1):7:1–7:4. doi:10.1145/2049662.2049669.
- **Byrd, R. H.; Nocedal, J.; Schnabel, R. B.** (1994). *Representations of Quasi-Newton Matrices and their Use in Limited Memory Methods*. Math. Programming 63:129–156.
- **Moré, J. J.; Thuente, D. J.** (1994). *Line Search Algorithms with Guaranteed Sufficient Decrease*. ACM TOMS 20(3):286–307. doi:10.1145/192115.192132.
- **Nocedal, J.; Wright, S. J.** (2006). *Numerical Optimization*, 2nd ed., Springer. Chapter 7 (L-BFGS two-loop), Chapters 16–17 (gradient projection, line-search theory).
- **Reference Fortran 77 (v3.0, 2011).** `http://users.iems.northwestern.edu/~nocedal/lbfgsb.html`. BSD-3 licensed.

### Reference implementations consulted

- **SciPy** `_lbfgsb_py.py` and `__lbfgsb` C core (post-1.15 rewrite of the f2py Fortran wrapper): the canonical API shape, option defaults, and reference behaviour. BSD-3.
- **L-BFGS-B-C** by Stephen Becker (`github.com/stephenbeckr/L-BFGS-B-C`): f2c of v3.0 with integer `task` codes — closest line-by-line template for a high-level language port. BSD-3.
- **LBFGSpp** by Yixuan Qiu (`github.com/yixuan/LBFGSpp`): C++ header-only. `BFGSMat.h` is the clearest implementation of the compact representation. Note: its subspace minimisation uses the BOXCQP algorithm (Voglis & Lagaris 2004) rather than the classical direct primal method; do **not** follow it if scipy parity is required. MIT.
- **pure Python `lbfgsb`** on PyPI (`pypi.org/project/lbfgsb/`): pure-Python port of Algorithm 778 — useful analogue for a high-level language implementation. MIT.
- **Matlab reference** by bgranzow (`github.com/bgranzow/L-BFGS-B`): ~400-line readable structural reference.
- **jacobwilliams/lbfgsb** and **jonathanschilling/L-BFGS-B** — modern Fortran 2008 refactors with Doxygen comments; good for cross-checking the control flow of `setulb`/`mainlb`.

### Current JS/TS landscape

As of April 2026, **no mature JavaScript or TypeScript implementation of L-BFGS-B exists on npm.** The npm packages `fmin`, `optimization-js`, `ml-fmin`, and `bfgs-algorithm` provide unconstrained BFGS / L-BFGS only. No WASM port of Nocedal's Fortran v3.0 is published. This specification therefore fills a genuine gap and the implementation should aim to be the definitive TypeScript port.

---

## Appendix A. Workspace sizing (v3.0 reference)

For an implementation that mirrors the Fortran memory layout for easier parity-debugging:

- Real workspace `wa`: `2mn + 11m² + 5n + 8m` doubles.
- Integer workspace `iwa`: `3n` ints.
- Persistent state: `isave[44]`, `dsave[29]`, `lsave[4]`, `csave[60]` — in TS, collapse into a single `SolverState` record.

## Appendix B. Exit messages (scipy-compatible strings)

```
CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL
CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH
STOP: TOTAL_NO_OF_ITERATIONS_REACHED_LIMIT
STOP: TOTAL_NO_OF_F,G_EVALUATIONS_EXCEEDS_LIMIT
ABNORMAL_TERMINATION_IN_LNSRCH
ERROR: NO_OR_INFINITE_VALUE_IN_F_OR_G
ERROR: INVALID_BOUNDS (l > u)
```

Emit these verbatim so downstream tooling that parses the scipy message strings works unchanged.

## Appendix C. Parameter cross-reference (SciPy ↔ this API ↔ Fortran)

| Concept | SciPy (`minimize`) | This API | Fortran |
|---|---|---|---|
| memory | `maxcor` | `m` | `m` |
| function rel. tol. | `ftol` | `factr · εmach` | `factr · epsmch` |
| projected-grad tol. | `gtol` | `pgtol` | `pgtol` |
| max outer iter | `maxiter` | `maxIterations` | (driver) |
| max fevals | `maxfun` | `maxFunctionEvaluations` | (driver) |
| max line-search iter | `maxls` | `maxLineSearch` | `iback` (=20) |
| Armijo $c_1$ | (hard-coded) | `ftol` | `ftol` (1e-3) |
| curvature $c_2$ | (hard-coded) | `gtol` | `gtol` (0.9) |
| interval tol | (hard-coded) | `xtol` | `xtol` (0.1) |

End of specification.