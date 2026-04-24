# Multi-Start Unconstrained Benchmarks

Complements [`unconstrained-benchmarks.md`](./unconstrained-benchmarks.md) (one x₀ per
problem) by running every optimizer from **three** starting points per problem,
specifically chosen to expose failure modes rather than measure average performance.

Each problem uses three x₀'s:

- **A (baseline)** — the x₀ used in the single-start benchmark
- **B (perturbed / adversarial)** — a small non-integer shift that breaks
  algorithm-specific luck (notably: zero-gradient contributions from `sin(kπxᵢ)`
  when xᵢ is integer- or half-integer aligned)
- **C (near-optimum)** — a probe close to the known global minimum, distinguishing
  *"no chance from this x₀"* from *"algorithm can't finish the job"*

Legend:

- ✅ found global (converged & `max(|error|, dist) < 1e-3`)
- ⚠️ converged but landed far from global (local minimum / inaccurate)
- ❌ did not converge

Each cell: `<status> <value> (<iterations> it, <fn-evals> fn)`.

---

## Sphere — unimodal, convex

Known min: f(0, 0) = 0. x₀A = (5, −3); x₀B = (−10, 10); x₀C = (100, 100).

| Method | x₀A (baseline) | x₀B (mirrored far) | x₀C (very far) |
|--------|----------------|--------------------|----------------|
| Nelder-Mead | ✅ 1.33e-10 (49 it, 97 fn) | ✅ 6.83e-10 (47 it, 87 fn) | ✅ 2.86e-10 (55 it, 103 fn) |
| PSO | ✅ 1.74e-12 (120 it, 6050 fn) | ✅ 1.16e-13 (128 it, 6450 fn) | ✅ 2.04e-13 (141 it, 7100 fn) |
| GradientDescent | ✅ 1.30e-8 (419 it, 2096 fn) | ✅ 1.30e-8 (452 it, 2261 fn) | ✅ 1.29e-8 (538 it, 2691 fn) |
| Adam | ✅ 1.37e-10 (243 it, 1216 fn) | ✅ 1.86e-14 (383 it, 1916 fn) | ✅ 2.76e-7 (3573 it, 17866 fn) |
| L-BFGS | ✅ 5.25e-16 (1 it, 11 fn) | ✅ 1.27e-15 (1 it, 11 fn) | ✅ 1.28e-11 (1 it, 11 fn) |

All optimizers are robust on the sphere; L-BFGS converges in one iteration from any
x₀ (isotropic quadratic).

---

## Rosenbrock — unimodal, narrow valley

Known min: f(1, 1) = 0. x₀A = (−1.2, 1.0); x₀B = (0, 0); x₀C = (2, −2).

| Method | x₀A (classic) | x₀B (origin) | x₀C (past optimum) |
|--------|---------------|--------------|---------------------|
| Nelder-Mead | ✅ 6.11e-11 (86 it, 162 fn) | ✅ 3.69e-10 (78 it, 146 fn) | ✅ 1.18e-9 (85 it, 160 fn) |
| PSO | ✅ 5.16e-9 (196 it, 9850 fn) | ✅ 1.08e-8 (180 it, 9050 fn) | ✅ 7.88e-11 (187 it, 9400 fn) |
| GradientDescent | ⚠️ 7.90e-7 (1414 it, 7071 fn) | ⚠️ 7.89e-7 (1390 it, 6951 fn) | ⚠️ 79.45 (11 it, 60 fn) |
| Adam | ✅ 1.85e-7 (1223 it, 6116 fn) | ⚠️ 1.65e-5 (308 it, 1541 fn) | ⚠️ 0.615 (78 it, 391 fn) |
| L-BFGS | ✅ 2.39e-12 (671 it, 3387 fn) | ✅ 3.73e-16 (26 it, 142 fn) | ✅ 5.38e-15 (28 it, 156 fn) |

Starting at the origin (B) is easier than the classic (A) for L-BFGS — only 26
iterations. GD is sensitive to x₀ in the narrow valley; past-optimum (C) sends it
upward into the banana.

---

## Beale — unimodal, non-convex

Known min: f(3, 0.5) = 0. x₀A = (0, 0); x₀B = (−2, −2); x₀C = (4, 0.4).

| Method | x₀A (baseline) | x₀B (opposite quadrant) | x₀C (near-optimum) |
|--------|----------------|-------------------------|---------------------|
| Nelder-Mead | ✅ 5.53e-10 (82 it, 161 fn) | ⚠️ 0.452 (2216 it, 4066 fn) | ✅ 2.22e-10 (36 it, 72 fn) |
| PSO | ✅ 1.64e-12 (121 it, 6100 fn) | ✅ 2.51e-14 (136 it, 6850 fn) | ✅ 4.70e-13 (118 it, 5950 fn) |
| GradientDescent | ⚠️ 1.17e-6 (1536 it, 7681 fn) | ⚠️ 3.77 (11 it, 60 fn) | ⚠️ 1.19e-6 (2392 it, 11961 fn) |
| Adam | ✅ 1.55e-8 (399 it, 1996 fn) | ⚠️ 7.60e-7 (2464 it, 12321 fn) | ✅ 3.00e-9 (177 it, 886 fn) |
| L-BFGS | ✅ 1.96e-11 (9 it, 54 fn) | ⚠️ 0.453 (224 it, 1233 fn) | ✅ 2.88e-15 (16 it, 91 fn) |

Opposite quadrant (B) is a different basin — all **local** methods (NM, GD, L-BFGS)
fail there; only PSO's sampling finds the global.

---

## Booth — unimodal, convex

Known min: f(1, 3) = 0. x₀A = (0, 0); x₀B = (−5, −5); x₀C = (10, 10).

| Method | x₀A (baseline) | x₀B (opposite quadrant) | x₀C (far) |
|--------|----------------|-------------------------|-----------|
| Nelder-Mead | ✅ 3.02e-10 (70 it, 139 fn) | ⚠️ 0.0389 (35 it, 69 fn) | ✅ 1.03e-9 (63 it, 123 fn) |
| PSO | ✅ 3.95e-12 (114 it, 5750 fn) | ✅ 1.12e-12 (119 it, 6000 fn) | ✅ 7.59e-14 (139 it, 7000 fn) |
| GradientDescent | ✅ 1.25e-8 (367 it, 1836 fn) | ✅ 1.25e-8 (367 it, 1836 fn) | ✅ 1.25e-8 (367 it, 1836 fn) |
| Adam | ✅ 1.40e-10 (297 it, 1486 fn) | ✅ 8.07e-8 (1020 it, 5101 fn) | ✅ 1.08e-7 (1149 it, 5746 fn) |
| L-BFGS | ✅ 9.89e-14 (8 it, 49 fn) | ✅ 6.00e-17 (6 it, 39 fn) | ✅ 3.19e-12 (5 it, 34 fn) |

Everyone finds the global (convex); only NM at B gets stuck with its noImprovementLimit
trigger before reaching full precision.

---

## Matyas — unimodal, convex

Known min: f(0, 0) = 0. x₀A = (5, −5); x₀B = (−8, 8); x₀C = (0.1, 0.1).

| Method | x₀A (baseline) | x₀B (mirrored) | x₀C (near-optimum) |
|--------|----------------|-----------------|---------------------|
| Nelder-Mead | ✅ 1.52e-10 (43 it, 82 fn) | ✅ 2.57e-10 (44 it, 84 fn) | ✅ 4.95e-10 (25 it, 47 fn) |
| PSO | ✅ 5.42e-13 (124 it, 6250 fn) | ✅ 1.08e-12 (124 it, 6250 fn) | ✅ 9.13e-13 (79 it, 4000 fn) |
| GradientDescent | ✅ 1.46e-7 (858 it, 4291 fn) | ✅ 1.46e-7 (900 it, 4501 fn) | ⚠️ 1.20e-5 (4379 it, 21896 fn) |
| Adam | ⚠️ 3.37e-6 (137 it, 686 fn) | ✅ 8.36e-10 (268 it, 1341 fn) | ✅ 2.50e-15 (51 it, 256 fn) |
| L-BFGS | ✅ 2.00e-16 (1 it, 10 fn) | ✅ 5.03e-16 (1 it, 10 fn) | ✅ 1.64e-19 (2 it, 15 fn) |

---

## Himmelblau — multimodal, 4 equivalent minima

Known global min: f(3, 2) = 0 (and 3 other equivalent minima).
x₀A = (0, 0); x₀B = (−4, −4); x₀C = (3, −3). "Distance to optimum" is measured to
the **nearest** of the 4 equivalent global minima, so ✅ is expected everywhere.

| Method | x₀A → basin (3,2) | x₀B → basin (−3.78,−3.28) | x₀C → basin (3.58,−1.85) |
|--------|-------------------|-----------------------------|------------------------------|
| Nelder-Mead | ✅ 1.43e-10 (86 it, 168 fn) | ✅ 2.67e-10 (41 it, 81 fn) | ✅ 1.09e-10 (46 it, 89 fn) |
| PSO | ✅ 5.50e-13 (130 it, 6550 fn) | ✅ 4.00e-15 (144 it, 7250 fn) | ✅ 1.13e-14 (144 it, 7250 fn) |
| GradientDescent | ✅ 1.03e-12 (260 it, 1301 fn) | ✅ 1.13e-6 (138 it, 691 fn) | ✅ 2.84e-10 (235 it, 1176 fn) |
| Adam | ✅ 3.76e-9 (232 it, 1161 fn) | ✅ 5.41e-11 (211 it, 1056 fn) | ✅ 1.57e-10 (226 it, 1131 fn) |
| L-BFGS | ✅ 6.19e-12 (9 it, 57 fn) | ✅ 1.05e-15 (7 it, 46 fn) | ✅ 1.30e-12 (8 it, 51 fn) |

---

## Three-Hump Camel — multimodal (3 local minima)

Known global min: f(0, 0) = 0. x₀A = (2, −1); x₀B = (−2, 1); x₀C = (0.5, −0.5).

| Method | x₀A (local basin) | x₀B (mirror, local basin) | x₀C (near global) |
|--------|-------------------|-----------------------------|--------------------|
| Nelder-Mead | ⚠️ 0.299 (33 it, 67 fn) | ⚠️ 0.299 (33 it, 67 fn) | ✅ 1.42e-10 (43 it, 84 fn) |
| PSO | ✅ 1.84e-13 (127 it, 6400 fn) | ✅ 3.39e-12 (118 it, 5950 fn) | ✅ 5.78e-14 (110 it, 5550 fn) |
| GradientDescent | ⚠️ 0.299 (289 it, 1446 fn) | ⚠️ 0.299 (289 it, 1446 fn) | ✅ 3.71e-8 (430 it, 2151 fn) |
| Adam | ⚠️ 0.299 (169 it, 846 fn) | ⚠️ 0.299 (169 it, 846 fn) | ✅ 3.26e-10 (167 it, 836 fn) |
| L-BFGS | ⚠️ 0.299 (9 it, 55 fn) | ⚠️ 0.299 (9 it, 55 fn) | ✅ 2.06e-11 (5 it, 31 fn) |

All gradient/simplex methods trap in the local basin near (1.75, −0.87) from A and B;
only PSO escapes.

---

## Rastrigin — highly multimodal

Known global min: f(0, 0) = 0. x₀A = (2.5, −3.5) *half-integer*; x₀B = (2.6, −3.4)
*+0.1 perturbation*; x₀C = (0.3, 0.4) *near global*.

| Method | x₀A (integer, "lucky") | x₀B (perturbed) | x₀C (near global) |
|--------|---------------------------|-----------------|---------------------|
| Nelder-Mead | ⚠️ 24.87 (42 it, 85 fn) | ⚠️ 17.91 (39 it, 79 fn) | ⚠️ 0.140 (22 it, 44 fn) |
| PSO | ✅ 4.26e-14 (141 it, 7100 fn) | ✅ 1.85e-11 (134 it, 6750 fn) | ✅ 5.95e-13 (145 it, 7300 fn) |
| GradientDescent | ⚠️ 12.93 (223 it, 1116 fn) | ⚠️ 17.91 (224 it, 1121 fn) | ✅ 9.52e-11 (242 it, 1211 fn) |
| Adam | ⚠️ 12.93 (243 it, 1216 fn) | ⚠️ 17.91 (228 it, 1141 fn) | ✅ 4.76e-11 (239 it, 1196 fn) |
| L-BFGS | ✅ 3.73e-14 (1 it, 11 fn) | ⚠️ 19.90 (8 it, 56 fn) | ⚠️ 4.97 (9 it, 61 fn) |

> **This is the key finding.** L-BFGS solves A in one iteration (global optimum,
> error 4e-14) because at half-integer xᵢ, `sin(2πxᵢ) = 0` and the cosine term drops
> out of the gradient — the landscape reduces to a quadratic and Armijo backtracking
> lands at the origin. **Shift x₀ by 0.1 (B), and the trick evaporates:** L-BFGS traps
> at the nearest well with value 19.90, worse than Adam/GD. From the near-global x₀
> (C) L-BFGS is also locked into a nearby shallow ripple (value 4.97), while gradient
> methods with small steps slide smoothly into the origin.

---

## Ackley — multimodal

Known global min: f(0, 0) = 0. x₀A = (2, −2); x₀B = (3.5, 3.5); x₀C = (0.5, 0.5).

| Method | x₀A (integer) | x₀B (mid-range non-integer) | x₀C (near global) |
|--------|---------------|------------------------------|--------------------|
| Nelder-Mead | ⚠️ 6.56 (30 it, 60 fn) | ⚠️ 10.12 (37 it, 74 fn) | ✅ 2.34e-9 (82 it, 162 fn) |
| PSO | ✅ 6.76e-11 (198 it, 9950 fn) | ✅ 8.22e-11 (208 it, 10450 fn) | ✅ 5.01e-11 (175 it, 8800 fn) |
| GradientDescent | ⚠️ 6.56 (127 it, 636 fn) | ⚠️ 9.00 (184 it, 921 fn) | ⚠️ 0.002 (124 it, 621 fn) |
| Adam | ⚠️ 6.56 (137 it, 686 fn) | ⚠️ 9.00 (147 it, 736 fn) | ⚠️ 0.0015 (88 it, 441 fn) |
| L-BFGS | ✅ 8.39e-12 (18 it, 107 fn) | ✅ 5.47e-13 (19 it, 111 fn) | ✅ 1.91e-10 (17 it, 100 fn) |

Ackley is L-BFGS-friendly *even off-integer*: the ripple amplitude decays with `exp(−0.2·r)`, so
a big Armijo step blows through the ripples into the smooth bowl near origin.

---

## Lévi N.13 — multimodal

Known global min: f(1, 1) = 0. x₀A = (−4, 5) *integer*; x₀B = (−3.9, 5.1) *+0.1*;
x₀C = (1.5, 1.5) *near global*.

| Method | x₀A (integer, "lucky") | x₀B (perturbed) | x₀C (near global) |
|--------|---------------------------|-----------------|---------------------|
| Nelder-Mead | ⚠️ 15.97 (60 it, 118 fn) | ⚠️ 20.22 (54 it, 113 fn) | ⚠️ 0.632 (34 it, 69 fn) |
| PSO | ✅ 9.78e-15 (129 it, 6500 fn) | ✅ 1.05e-13 (124 it, 6250 fn) | ✅ 1.31e-13 (133 it, 6700 fn) |
| GradientDescent | ⚠️ 0.110 (284 it, 1421 fn) | ⚠️ 0.632 (287 it, 1436 fn) | ⚠️ 0.110 (267 it, 1336 fn) |
| Adam | ⚠️ 16.97 (253 it, 1266 fn) | ⚠️ 16.97 (263 it, 1316 fn) | ⚠️ 0.110 (211 it, 1056 fn) |
| L-BFGS | ✅ 1.44e-14 (1 it, 11 fn) | ⚠️ 22.50 (9 it, 61 fn) | ⚠️ 0.110 (8 it, 47 fn) |

Same pattern as Rastrigin: at integer x₀, all `sin(kπxᵢ)` terms vanish, L-BFGS solves
in one iteration. +0.1 perturbation destroys the effect.

---

## Griewank — multimodal, shallow ripples

Known global min: f(0, 0) = 0. x₀A = (100, −200); x₀B = (5, 5); x₀C = (0.5, −0.5).

| Method | x₀A (far) | x₀B (moderate) | x₀C (near global) |
|--------|-----------|----------------|--------------------|
| Nelder-Mead | ⚠️ 3.74 (17 it, 33 fn) | ⚠️ 0.0074 (44 it, 84 fn) | ✅ 1.53e-10 (38 it, 73 fn) |
| PSO | ✅ 7.01e-13 (238 it, 11950 fn) | ✅ 3.01e-13 (126 it, 6350 fn) | ✅ 4.00e-15 (178 it, 8950 fn) |
| GradientDescent | ⚠️ 12.35 (978 it, 4891 fn) | ⚠️ 0.0074 (1257 it, 6286 fn) | ⚠️ 5.58e-7 (1123 it, 5616 fn) |
| Adam | ⚠️ 12.35 (228 it, 1141 fn) | ⚠️ 0.0074 (200 it, 1001 fn) | ✅ 3.05e-10 (170 it, 851 fn) |
| L-BFGS | ⚠️ 12.35 (7 it, 40 fn) | ⚠️ 0.0074 (8 it, 46 fn) | ✅ 1.62e-11 (4 it, 25 fn) |

---

## Styblinski-Tang — multimodal

Known global min: f(−2.9035, −2.9035) ≈ −78.3323.
x₀A = (0, 0); x₀B = (3, 3); x₀C = (−2, −2).

| Method | x₀A (baseline) | x₀B (opposite-sign basin) | x₀C (near global) |
|--------|----------------|-----------------------------|--------------------|
| Nelder-Mead | ⚠️ −64.20 (90 it, 177 fn) | ⚠️ −50.06 (36 it, 72 fn) | ✅ −78.33 (38 it, 71 fn) |
| PSO | ✅ −78.33 (124 it, 6250 fn) | ✅ −78.33 (147 it, 7400 fn) | ✅ −78.33 (127 it, 6400 fn) |
| GradientDescent | ✅ −78.33 (227 it, 1136 fn) | ⚠️ −50.06 (154 it, 771 fn) | ✅ −78.33 (199 it, 996 fn) |
| Adam | ⚠️ −78.33 (78 it, 391 fn) | ⚠️ −50.06 (80 it, 401 fn) | ✅ −78.33 (197 it, 986 fn) |
| L-BFGS | ✅ −78.33 (7 it, 44 fn) | ⚠️ −50.06 (6 it, 39 fn) | ✅ −78.33 (6 it, 39 fn) |

x₀B lies in the opposite basin (minimum at (+2.75, +2.75), value −50.06) — all local
methods are trapped; only PSO escapes.

---

## Easom — nearly flat, narrow peak

Known global min: f(π, π) = −1. x₀A = (1, 1); x₀B = (2.5, 2.5); x₀C = (3, 3).

| Method | x₀A (flat region) | x₀B (closer to peak) | x₀C (near peak) |
|--------|-------------------|-----------------------|-----------------|
| Nelder-Mead | ⚠️ −8.1e-5 (20 it, 40 fn) | ✅ −1.000 (34 it, 67 fn) | ⚠️ −0.9998 (7 it, 16 fn) |
| PSO | ✅ −1.000 (112 it, 5650 fn) | ✅ −1.000 (108 it, 5450 fn) | ✅ −1.000 (113 it, 5700 fn) |
| GradientDescent | ⚠️ −3.0e-5 (50 it, 251 fn) | ✅ −1.000 (245 it, 1226 fn) | ✅ −1.000 (181 it, 906 fn) |
| Adam | ⚠️ −8.1e-5 (58 it, 291 fn) | ✅ −1.000 (187 it, 936 fn) | ✅ −1.000 (128 it, 641 fn) |
| L-BFGS | ❌ −4.8e-5 (1000 it, 5005 fn) | ✅ −1.000 (4 it, 27 fn) | ✅ −1.000 (4 it, 26 fn) |

At x₀A, the gradient is ~1e-9 everywhere — **L-BFGS burns 1000 iterations and
doesn't converge** because the finite-difference gradient is below numerical noise.
Once you get within the peak's basin (B, C), L-BFGS is the fastest (4 iterations).
NM at C exits via noImprovement before reaching the sub-millidecimal precision.

---

## Goldstein-Price — multimodal

Known global min: f(0, −1) = 3. x₀A = (0, 0); x₀B = (1, 1); x₀C = (−0.5, −1).

| Method | x₀A (baseline) | x₀B (away from global) | x₀C (near global) |
|--------|----------------|-------------------------|--------------------|
| Nelder-Mead | ⚠️ 30.00 (75 it, 143 fn) | ⚠️ 840.11 (12 it, 24 fn) | ✅ 3.000 (49 it, 95 fn) |
| PSO | ✅ 3.000 (136 it, 6850 fn) | ✅ 3.000 (138 it, 6950 fn) | ✅ 3.000 (127 it, 6400 fn) |
| GradientDescent | ⚠️ 334.33 (5 it, 30 fn) | ⚠️ 1876.00 (3 it, 20 fn) | ⚠️ 279.13 (3 it, 20 fn) |
| Adam | ⚠️ 30.00 (255 it, 1276 fn) | ⚠️ 84.00 (257 it, 1286 fn) | ⚠️ 3.017 (61 it, 306 fn) |
| L-BFGS | ⚠️ 30.00 (11 it, 71 fn) | ⚠️ 84.00 (14 it, 89 fn) | ✅ 3.000 (12 it, 76 fn) |

GP has high-curvature plateaus that trap momentum-based methods early; GD's
`noImprovementLimit=50` triggers within a handful of iterations because the
quadratic-penalty gradient explodes. PSO's population sampling sidesteps this.

---

## McCormick — multimodal, unbounded domain

Known global min: f(−0.547, −1.547) ≈ −1.9133.
x₀A = (0, 0); x₀B = (−1, −2); x₀C = (2, 2).

| Method | x₀A (baseline) | x₀B (near global) | x₀C (opposite side) |
|--------|----------------|--------------------|----------------------|
| Nelder-Mead | ✅ −1.913 (79 it, 153 fn) | ✅ −1.913 (37 it, 70 fn) | ⚠️ 1.228 (33 it, 66 fn) |
| PSO | ❌ −38424 (10000 it, 500050 fn) | ❌ −27158 (10000 it, 500050 fn) | ❌ −75823 (10000 it, 500050 fn) |
| GradientDescent | ✅ −1.913 (417 it, 2086 fn) | ✅ −1.913 (419 it, 2096 fn) | ⚠️ 1.228 (322 it, 1611 fn) |
| Adam | ✅ −1.913 (193 it, 966 fn) | ✅ −1.913 (156 it, 781 fn) | ⚠️ 1.228 (190 it, 951 fn) |
| L-BFGS | ✅ −1.913 (7 it, 41 fn) | ✅ −1.913 (4 it, 25 fn) | ⚠️ 1.228 (5 it, 31 fn) |

**PSO diverges from every x₀** because McCormick is unbounded along certain directions
and the algorithm needs explicit box bounds. All local methods find the global from
A and B; x₀C sits on the wrong side of the `(x − y)²` trough and traps them.

---

## Summary

5 optimizers × 15 problems × 3 starting points = 225 runs.

| Optimizer | Global optima found | Success rate |
|-----------|---------------------|--------------|
| **PSO** | **42 / 45** | **93.3%** |
| L-BFGS | 31 / 45 | 68.9% |
| Nelder-Mead | 24 / 45 | 53.3% |
| Adam | 22 / 45 | 48.9% |
| GradientDescent | 19 / 45 | 42.2% |

PSO's sampling wins in aggregate — global by construction, at the cost of 40–700×
more function evaluations. L-BFGS is second overall but its results depend heavily
on x₀ — the "1-iteration wins" on Rastrigin / Lévi N.13 at integer x₀ flip to
failures under a 0.1 perturbation (as the B columns show). The takeaway is not that
one optimizer is universally best, but that **x₀ sensitivity must be considered**
when choosing a local method for a multimodal problem; for such problems either use
PSO directly or wrap a local optimizer in a multi-start strategy.

PSO's single blind spot is unbounded multimodal objectives (McCormick) — use box
constraints when the domain is not naturally bounded.

## Running

```bash
npx tsx src/optimization/single-objective/benchmarks/multistart-benchmarks.ts
```
