# Bounded Optimization Benchmarks

6 optimizers, 7 problems. Compares **L-BFGS-B** (native box constraints via
`settings.bounds`, geometric handling through Cauchy point + subspace
minimisation) against other optimizers handling bounds **via the penalty layer**
(`boxConstraints()` → quadratic penalty).

Success requires **both** numerical accuracy and **exact feasibility**:

* ✅ converged, `feas_vio < 1e-6`, and `|Error| < 1e-3` or `Dist < 1e-3`
* ⚠️ converged but wrong answer (accuracy failure)
* ❌ did not converge **or** infeasible (any coord overshoots a bound by ≥ 1e-6)

Problems are sorted by `activeBoundsAtOptimum` ascending, so
"bounds inactive" sanity checks surface first — any regression there indicates
L-BFGS-B is mishandling a scenario that is effectively unconstrained.

> **Why feasibility matters.** Penalty-based methods trade bound violation
> against the original objective; near an active bound the optimum of `f + μ·P`
> typically sits *inside* the infeasible region by `O(1/√μ)`. The tables below
> show violations of `~1e-3` on problems with active bounds — small, but enough
> to disqualify the result when the caller needs strict feasibility (which, for
> box-constrained problems, is typically the whole point).

---

## Summary

| Optimizer | Mode | Success | Notes |
|-----------|------|--------:|-------|
| **L-BFGS-B** | native | **7/7** | Exact feasibility on every problem |
| L-BFGS | penalty | 3/7 | Fails active-bound problems (violates by ~1e-3) |
| PSO | penalty | 3/7 | Same failure mode |
| Adam | penalty | 2/7 | Plus occasional accuracy loss on sharp optima |
| GradientDescent | penalty | 2/7 | Same as Adam/LBFGS |
| Nelder-Mead | penalty | 2/7 | Lower accuracy on high-dim |

---

## Himmelblau (bounds inactive)

- **Dim:** 2 | **Type:** Bounds should not matter | **Active bounds at opt:** 0/2
- **Known min:** f(3, 2) = 0 | **x₀:** (0, 0) | **Bounds:** `[-4, 4]²`

|   | Method | Mode | Found Value | \|Error\| | Dist to Opt | Feas Vio | Conv | Iters | Fn Evals |
|---|--------|------|------------:|----------:|------------:|---------:|------|------:|---------:|
| ✅ | L-BFGS-B | native | 4.63e-12 | 4.63e-12 | 6.00e-7 | 0 | yes | 8 | 65 |
| ✅ | L-BFGS | penalty | 6.19e-12 | 6.19e-12 | 3.88e-7 | 0 | yes | 9 | 57 |
| ✅ | Adam | penalty | 3.76e-9 | 3.76e-9 | 1.65e-5 | 0 | yes | 232 | 1161 |
| ✅ | GradientDescent | penalty | 1.03e-12 | 1.03e-12 | 1.72e-7 | 0 | yes | 260 | 1301 |
| ✅ | PSO | penalty | 1.85e-13 | 1.85e-13 | 1.02e-7 | 0 | yes | 135 | 6800 |
| ✅ | Nelder-Mead | penalty | 1.43e-10 | 1.43e-10 | 1.90e-6 | 0 | yes | 86 | 168 |

All methods succeed — Himmelblau's minima lie inside the box. Sanity check
passes: L-BFGS-B does not degrade on effectively-unconstrained problems.

---

## Ackley (bounds inactive)

- **Dim:** 2 | **Type:** Multimodal, bounds inactive | **Active bounds at opt:** 0/2
- **Known min:** f(0, 0) = 0 | **x₀:** (1.5, −1.5) | **Bounds:** `[-2, 2]²`

|   | Method | Mode | Found Value | \|Error\| | Dist to Opt | Feas Vio | Conv | Iters | Fn Evals |
|---|--------|------|------------:|----------:|------------:|---------:|------|------:|---------:|
| ✅ | L-BFGS-B | native | 1.09e-10 | 1.09e-10 | 3.84e-11 | 0 | yes | 7 | 145 |
| ✅ | L-BFGS | penalty | 1.55e-10 | 1.55e-10 | 5.48e-11 | 0 | yes | 10 | 73 |
| ⚠️ | Adam | penalty | 3.5745 | 3.5745 | 1.3697 | 0 | yes | 117 | 586 |
| ⚠️ | GradientDescent | penalty | 3.5745 | 3.5745 | 1.3696 | 0 | yes | 198 | 991 |
| ✅ | PSO | penalty | 2.40e-10 | 2.40e-10 | 8.48e-11 | 0 | yes | 188 | 9450 |
| ⚠️ | Nelder-Mead | penalty | 5.3819 | 5.3819 | 2.1966 | 0 | yes | 42 | 84 |

Quasi-Newton methods (L-BFGS-B, L-BFGS) and PSO reach the global min; first-order
methods get stuck in local minima of Ackley's rugged surface.

---

## Bounded Rosenbrock

- **Dim:** 2 | **Type:** Upper bound active at solution | **Active bounds at opt:** 1/2
- **Known min:** f(0.5, 0.25) = 0.25 | **x₀:** (0, 0) | **Bounds:** `[-1.5, 0.5]²`

|   | Method | Mode | Found Value | \|Error\| | Dist to Opt | Feas Vio | Conv | Iters | Fn Evals |
|---|--------|------|------------:|----------:|------------:|---------:|------|------:|---------:|
| ✅ | L-BFGS-B | native | 0.2500 | 0 | 6.27e-13 | **0** | yes | 7 | 55 |
| ❌ | L-BFGS | penalty | 0.2498 | 2.50e-4 | 7.07e-4 | 5.00e-4 | yes | 19 | 113 |
| ❌ | Adam | penalty | 0.2500 | 7.92e-6 | 9.46e-4 | 7.93e-4 | yes | 125 | 626 |
| ❌ | GradientDescent | penalty | 0.2498 | 2.50e-4 | 7.07e-4 | 5.00e-4 | yes | 227 | 1136 |
| ❌ | PSO | penalty | 0.2498 | 2.50e-4 | 7.07e-4 | 4.99e-4 | yes | 151 | 7600 |
| ⚠️ | Nelder-Mead | penalty | 0.2516 | 0.0016 | 0.0038 | 0 | yes | 41 | 78 |

Every penalty method violates the upper bound of `x₀` by `~5·10⁻⁴`. L-BFGS-B
sits exactly on the bound.

---

## Bounded Beale

- **Dim:** 2 | **Type:** Upper bound just touches optimum | **Active bounds at opt:** 1/2
- **Known min:** f(3, 0.5) = 0 | **x₀:** (1, 1) | **Bounds:** `[0, 3]²`

|   | Method | Mode | Found Value | \|Error\| | Dist to Opt | Feas Vio | Conv | Iters | Fn Evals |
|---|--------|------|------------:|----------:|------------:|---------:|------|------:|---------:|
| ✅ | L-BFGS-B | native | 3.49e-13 | 3.49e-13 | 1.47e-6 | 0 | yes | 13 | 80 |
| ✅ | L-BFGS | penalty | 5.29e-11 | 5.29e-11 | 1.49e-6 | 3.17e-8 | yes | 19 | 115 |
| ✅ | Adam | penalty | 9.53e-10 | 9.53e-10 | 7.94e-5 | 0 | yes | 265 | 1326 |
| ✅ | GradientDescent | penalty | 1.17e-6 | 1.17e-6 | 0.0028 | 0 | yes | 1528 | 7641 |
| ✅ | PSO | penalty | 7.26e-13 | 7.26e-13 | 2.03e-6 | 0 | yes | 118 | 5950 |
| ✅ | Nelder-Mead | penalty | 1.81e-7 | 1.81e-7 | 3.62e-4 | 0 | yes | 48 | 92 |

All methods succeed because the upper bound only "touches" the optimum rather
than forcing it. L-BFGS and Beale agree because the gradient is small on the
boundary so the quadratic penalty barely activates.

---

## Fixed Variables Sphere

- **Dim:** 5 | **Type:** 2 fixed (l = u), 3 free | **Active bounds at opt:** 2/5
- **Known min:** f = 5 at (0, 2, 0, −1, 0) | **x₀:** (5, 2, 5, −1, 5)

|   | Method | Mode | Found Value | \|Error\| | Dist to Opt | Feas Vio | Conv | Iters | Fn Evals |
|---|--------|------|------------:|----------:|------------:|---------:|------|------:|---------:|
| ✅ | L-BFGS-B | native | 5.0000 | 6.04e-14 | 2.48e-7 | **0** | yes | 2 | 33 |
| ❌ | L-BFGS | penalty | 4.9950 | 0.0050 | 0.0022 | 0.0020 | yes | 11 | 138 |
| ❌ | Adam | penalty | 4.9950 | 0.0050 | 0.0022 | 0.0020 | yes | 245 | 2696 |
| ❌ | GradientDescent | penalty | 4.9950 | 0.0050 | 0.0022 | 0.0020 | yes | 434 | 4775 |
| ❌ | PSO | penalty | 4.9950 | 0.0050 | 0.0022 | 0.0020 | yes | 191 | 9600 |
| ❌ | Nelder-Mead | penalty | 4.9950 | 0.0050 | 0.0022 | 0.0020 | yes | 187 | 315 |

L-BFGS-B recognises `l = u` as equality constraints and fixes those coordinates
exactly. Penalty methods drift away from the equality by `~2·10⁻³`.

---

## Shifted Sphere, half-bounded

- **Dim:** 4 | **Type:** Lower-only bound; all coords pinned | **Active bounds at opt:** 4/4
- **Known min:** f(0, 0, 0, 0) = 4 | **x₀:** (2, 2, 2, 2) | **Bounds:** `l = 0, u = +∞`

|   | Method | Mode | Found Value | \|Error\| | Dist to Opt | Feas Vio | Conv | Iters | Fn Evals |
|---|--------|------|------------:|----------:|------------:|---------:|------|------:|---------:|
| ✅ | L-BFGS-B | native | 4.0000 | 0 | 0 | **0** | yes | 2 | 36 |
| ❌ | L-BFGS | penalty | 3.9960 | 0.0040 | 0.0020 | 9.99e-4 | yes | 9 | 111 |
| ❌ | Adam | penalty | 3.9960 | 0.0040 | 0.0020 | 9.99e-4 | yes | 177 | 1594 |
| ❌ | GradientDescent | penalty | 3.9960 | 0.0040 | 0.0020 | 9.80e-4 | yes | 158 | 1423 |
| ❌ | PSO | penalty | 3.9960 | 0.0040 | 0.0020 | 9.99e-4 | yes | 196 | 9850 |
| ❌ | Nelder-Mead | penalty | 4.3195 | 0.3195 | 0.1411 | 0.0038 | yes | 205 | 350 |

All four variables are pinned at the lower bound at the optimum — the cleanest
demonstration of native-vs-penalty bound handling. L-BFGS-B converges in 2
iterations with machine-precision feasibility; penalty methods trade 1e-3
feasibility for 1e-3 objective error.

---

## Bounded Sphere (all pinned)

- **Dim:** 5 | **Type:** All lower bounds active | **Active bounds at opt:** 5/5
- **Known min:** f = 1.25 at (0.5, …, 0.5) | **x₀:** (1, 1, 1, 1, 1) | **Bounds:** `[0.5, 10]⁵`

|   | Method | Mode | Found Value | \|Error\| | Dist to Opt | Feas Vio | Conv | Iters | Fn Evals |
|---|--------|------|------------:|----------:|------------:|---------:|------|------:|---------:|
| ✅ | L-BFGS-B | native | 1.2500 | 0 | 0 | **0** | yes | 1 | 33 |
| ❌ | L-BFGS | penalty | 1.2488 | 0.0012 | 0.0011 | 5.00e-4 | yes | 4 | 66 |
| ⚠️ | Adam | penalty | 1.2901 | 0.0401 | 0.0178 | 0 | yes | 55 | 606 |
| ❌ | GradientDescent | penalty | 1.2488 | 0.0012 | 0.0011 | 5.03e-4 | yes | 154 | 1695 |
| ❌ | PSO | penalty | 1.2488 | 0.0012 | 0.0011 | 5.00e-4 | yes | 210 | 10550 |
| ❌ | Nelder-Mead | penalty | 1.7046 | 0.4546 | 0.2720 | 3.91e-4 | yes | 825 | 1333 |

L-BFGS-B solves this in **one iteration** — the gradient points into the
infeasible region at every coordinate, so Cauchy immediately snaps each
variable to its lower bound and the search exits feasibility-optimal.

---

## Notes

* `maxFunctionEvaluations` is capped at 15000 (L-BFGS-B default). PSO's
  `maxIterations=10000, swarmSize=50` hits 10550 evals on one problem because
  PSO counts differently (50 evals per generation).
* Penalty parameter `μ = 1000` (library default). Tightening `μ` improves
  feasibility at the cost of conditioning — the structural advantage of native
  bounds does not disappear.
* L-BFGS-B uses `maxIterations=1000, historySize=10` — identical to
  L-BFGS — so the wallclock advantage it shows on pinned problems is algorithmic
  (fewer iterations needed), not hyperparameter-driven.

To regenerate this report: `npx tsx bounded-benchmarks.ts`.
