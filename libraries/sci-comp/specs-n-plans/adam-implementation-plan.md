# Adam Optimizer — Implementation Plan

## Overview

Adam (Adaptive Moment Estimation) maintains per-parameter adaptive learning rates using exponentially decaying averages of the first moment (mean) and second moment (uncentered variance) of the gradient.

Reference: Kingma & Ba, 2014 — "Adam: A Method for Stochastic Optimization" ([arXiv:1412.6980](https://arxiv.org/abs/1412.6980)).

## Core Update Rule

```
gₜ      = ∇f(xₜ)                        // gradient (via finite differences)
mₜ      = β₁ · mₜ₋₁ + (1 − β₁) · gₜ    // first moment estimate
vₜ      = β₂ · vₜ₋₁ + (1 − β₂) · gₜ²   // second moment estimate
m̂ₜ      = mₜ / (1 − β₁ᵗ)               // bias-corrected first moment
v̂ₜ      = vₜ / (1 − β₂ᵗ)               // bias-corrected second moment
xₜ₊₁    = xₜ − lr · m̂ₜ / (√v̂ₜ + ε)     // parameter update
```

Bias correction is critical in early iterations: without it, `m` and `v` are biased toward zero because they are initialized at zero.

---

## Settings Interface

```typescript
export interface AdamSettings extends CommonSettings {
  /** Base learning rate α (default: 0.001). */
  learningRate?: number;
  /** Exponential decay rate for first moment β₁ (default: 0.9). */
  beta1?: number;
  /** Exponential decay rate for second moment β₂ (default: 0.999). */
  beta2?: number;
  /** Small constant for numerical stability (default: 1e-8). */
  epsilon?: number;
  /** Step size h for central finite-difference gradient (default: 1e-7). */
  finiteDiffStep?: number;
  /** Gradient norm convergence threshold (default: 1e-10). */
  gradTolerance?: number;
  /** How many iterations without improvement before stopping (default: 50). */
  noImprovementLimit?: number;
  /** Maximum gradient norm — clip before update (default: Infinity). */
  maxGradNorm?: number;
}
```

Note: no `lrSchedule`/`lrDecayRate`. Adam's adaptive per-parameter rates largely eliminate the need for manual LR scheduling. Can be added later if needed.

---

## File Structure

```
src/optimizers/adam.ts    — AdamSettings + Adam class
```

Single file. No changes to other files except:

- `index.ts` — add export and `registerOptimizer('adam', () => new Adam())`.

---

## Class Skeleton

```typescript
export class Adam extends Optimizer<AdamSettings> {
  constructor() { super('Adam'); }

  protected withDefaults(s: AdamSettings): AdamSettings { ... }
  protected runInternal(fn, x0, s): OptimizationResult { ... }
  protected async runInternalAsync(fn, x0, s): Promise<OptimizationResult> { ... }
}
```

Follows the same pattern as `GradientDescent`.

---

## Algorithm — Sync Path (`runInternal`)

### State

| Variable | Type | Description |
|---|---|---|
| `x` | `Float64Array(n)` | Current parameter vector |
| `m` | `Float64Array(n)` | First moment vector (initialized to 0) |
| `v` | `Float64Array(n)` | Second moment vector (initialized to 0) |
| `grad` | `Float64Array(n)` | Gradient buffer (reused each iteration) |
| `bestPoint` | `Float64Array(n)` | Best point found so far |
| `bestValue` | `number` | Best objective value |

### Per-Iteration Steps

```
1.  Evaluate f(x)
2.  Update bestPoint / bestValue if improved
3.  Record costHistory
4.  Check convergence (no-improvement streak)
5.  Fire onIteration callback
6.  Compute gradient via central finite differences       // reuse from GradientDescent
7.  Check gradient norm convergence (gradTolerance)
8.  Clip gradient if ‖∇f‖ > maxGradNorm
9.  Update moments:
      m[i] = β₁ · m[i] + (1 − β₁) · grad[i]
      v[i] = β₂ · v[i] + (1 − β₂) · grad[i]²
10. Bias correction:
      β₁ᵗ = β₁^(iteration+1)      // precompute as running product
      β₂ᵗ = β₂^(iteration+1)
      mHat = m[i] / (1 − β₁ᵗ)
      vHat = v[i] / (1 − β₂ᵗ)
11. Parameter update:
      x[i] -= lr · mHat / (√vHat + ε)
12. iteration++
```

### Bias Correction — Efficient Computation

Instead of `Math.pow(beta1, t)` every iteration, maintain running products:

```typescript
let beta1t = 1;  // β₁^t
let beta2t = 1;  // β₂^t

// each iteration:
beta1t *= beta1;
beta2t *= beta2;
const mCorr = 1 / (1 - beta1t);
const vCorr = 1 / (1 - beta2t);
```

This avoids `Math.pow` calls entirely.

---

## Algorithm — Async Path (`runInternalAsync`)

Identical logic, with `await fn(x)` and `await this.computeGradientAsync(...)`.

---

## Reuse from GradientDescent

Adam and GradientDescent share several concerns:

| Concern | Reuse strategy |
|---|---|
| `computeGradient` / `computeGradientAsync` | Duplicate into Adam (private methods). Extract to shared utility later if a third gradient method is added. |
| Convergence checks (no-improvement, grad norm) | Same inline logic pattern. Extract to base class later if warranted. |
| Gradient clipping | Same inline logic. |

Premature extraction adds complexity. Keep it duplicated for now (two copies); extract when the third gradient optimizer arrives.

---

## Convergence Criteria

Same two-gate approach as GradientDescent:

1. **No improvement**: `prevBest - bestValue > tolerance` not satisfied for `noImprovementLimit` consecutive iterations.
2. **Gradient norm**: `‖∇f‖ < gradTolerance`.

Either gate triggers `converged = true`.

---

## Default Values Summary

| Parameter | Default | Rationale |
|---|---|---|
| `learningRate` | `0.001` | Original paper recommendation |
| `beta1` | `0.9` | Original paper recommendation |
| `beta2` | `0.999` | Original paper recommendation |
| `epsilon` | `1e-8` | Original paper recommendation |
| `finiteDiffStep` | `1e-7` | Consistent with GradientDescent |
| `gradTolerance` | `1e-10` | Consistent with GradientDescent |
| `noImprovementLimit` | `50` | Consistent with GradientDescent |
| `maxGradNorm` | `Infinity` | No clipping by default |
| `maxIterations` | `10_000` | Consistent with GradientDescent |
| `tolerance` | `1e-8` | Consistent with GradientDescent |

---

## Registration

```typescript
// index.ts
export { Adam } from './optimizers/adam';
export type { AdamSettings } from './optimizers/adam';

import { Adam } from './optimizers/adam';
registerOptimizer('adam', () => new Adam());
```

---

## Testing Checklist

- [ ] Sphere function — should converge faster than GD with default settings
- [ ] Rosenbrock — verify it navigates the narrow valley
- [ ] Compare iteration count: Adam vs GD with momentum vs GD without momentum
- [ ] Verify bias correction: first few iterations should not have near-zero step sizes
- [ ] Async path: use `async (x) => fn(x)` wrapper to verify identical results
- [ ] Constraints via penalty: `minimize` with `boxConstraints` should respect bounds
- [ ] `maximize` / `maximizeAsync` — verify sign inversion works correctly
