# Unconstrained Optimization Benchmarks

5 optimizers, 15 problems, default hyperparameters.

Legend:

* ✅ converged & max(\|Error\|, Dist) < 0.001
* ⚠️ converged but inaccurate
* ❌ did not converge

Problems specification: [link](unconstrained-benchmark-functions.md)

---

## Sphere

- **Dim:** 2 | **Type:** Unimodal, convex
- **Known min:** f(0, 0) = 0
- **x0:** (5, −3)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ✅ | Nelder-Mead | maxIterations=10000 | 1.33e-10 | (1.15e-5, 8.71e-7) | 1.33e-10 | 1.15e-5 | yes | 49 | 97 | 0.4 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 1.74e-12 | (−6.71e-7, −1.14e-6) | 1.74e-12 | 1.32e-6 | yes | 120 | 6050 | 2.9 |
| ✅ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 1.30e-8 | (9.78e-5, −5.87e-5) | 1.30e-8 | 1.14e-4 | yes | 419 | 2096 | 0.6 |
| ✅ | Adam | maxIterations=10000, learningRate=0.1 | 1.37e-10 | (9.31e-8, 1.17e-5) | 1.37e-10 | 1.17e-5 | yes | 243 | 1216 | 0.4 |
| ✅ | L-BFGS | maxIterations=1000, historySize=10 | 5.25e-16 | (−2.29e-8, −4.69e-10) | 5.25e-16 | 2.29e-8 | yes | 1 | 11 | 0.2 |

---

## Rosenbrock

- **Dim:** 2 | **Type:** Unimodal, narrow valley
- **Known min:** f(1, 1) = 0
- **x0:** (−1.2, 1.0)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ✅ | Nelder-Mead | maxIterations=10000 | 6.11e-11 | (1.0000, 1.0000) | 6.11e-11 | 1.74e-5 | yes | 86 | 162 | 0.2 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 5.16e-9 | (1.0001, 1.0001) | 5.16e-9 | 1.43e-4 | yes | 196 | 9850 | 3.2 |
| ⚠️ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 7.90e-7 | (0.9991, 0.9982) | 7.90e-7 | 0.00199 | yes | 1414 | 7071 | 0.6 |
| ✅ | Adam | maxIterations=10000, learningRate=0.1 | 1.85e-7 | (0.9996, 0.9991) | 1.85e-7 | 9.62e-4 | yes | 1223 | 6116 | 0.6 |
| ✅ | L-BFGS | maxIterations=1000, historySize=10 | 2.39e-12 | (1.0000, 1.0000) | 2.39e-12 | 2.23e-6 | yes | 671 | 3387 | 1.3 |

---

## Beale

- **Dim:** 2 | **Type:** Unimodal, non-convex
- **Known min:** f(3, 0.5) = 0
- **x0:** (0, 0)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ✅ | Nelder-Mead | maxIterations=10000 | 5.53e-10 | (2.9999, 0.5000) | 5.53e-10 | 6.00e-5 | yes | 82 | 161 | 0.2 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 1.64e-12 | (3.0000, 0.5000) | 1.64e-12 | 2.69e-6 | yes | 121 | 6100 | 1.1 |
| ⚠️ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 1.17e-6 | (2.9973, 0.4993) | 1.17e-6 | 0.00279 | yes | 1536 | 7681 | 0.2 |
| ✅ | Adam | maxIterations=10000, learningRate=0.1 | 1.55e-8 | (2.9997, 0.4999) | 1.55e-8 | 3.20e-4 | yes | 399 | 1996 | 0.1 |
| ✅ | L-BFGS | maxIterations=1000, historySize=10 | 1.96e-11 | (3.0000, 0.5000) | 1.96e-11 | 3.43e-6 | yes | 9 | 54 | 0.1 |

---

## Booth

- **Dim:** 2 | **Type:** Unimodal, convex
- **Known min:** f(1, 3) = 0
- **x0:** (0, 0)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ✅ | Nelder-Mead | maxIterations=10000 | 3.02e-10 | (1.0000, 3.0000) | 3.02e-10 | 1.29e-5 | yes | 70 | 139 | 0.3 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 3.95e-12 | (1.0000, 3.0000) | 3.95e-12 | 1.94e-6 | yes | 114 | 5750 | 0.9 |
| ✅ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 1.25e-8 | (1.0001, 2.9999) | 1.25e-8 | 1.12e-4 | yes | 367 | 1836 | 0.1 |
| ✅ | Adam | maxIterations=10000, learningRate=0.1 | 1.40e-10 | (1.0000, 3.0000) | 1.40e-10 | 3.97e-6 | yes | 297 | 1486 | 0.1 |
| ✅ | L-BFGS | maxIterations=1000, historySize=10 | 9.89e-14 | (1.0000, 3.0000) | 9.89e-14 | 2.27e-7 | yes | 8 | 49 | 0.0 |

---

## Matyas

- **Dim:** 2 | **Type:** Unimodal, convex
- **Known min:** f(0, 0) = 0
- **x0:** (5, −5)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ✅ | Nelder-Mead | maxIterations=10000 | 1.52e-10 | (6.10e-5, 5.05e-5) | 1.52e-10 | 7.92e-5 | yes | 43 | 82 | 0.1 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 5.42e-13 | (−3.41e-6, −3.75e-6) | 5.42e-13 | 5.07e-6 | yes | 124 | 6250 | 0.4 |
| ✅ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 1.46e-7 | (3.82e-4, −3.82e-4) | 1.46e-7 | 5.40e-4 | yes | 858 | 4291 | 0.1 |
| ⚠️ | Adam | maxIterations=10000, learningRate=0.1 | 3.37e-6 | (−0.0018, 0.0018) | 3.37e-6 | 0.00260 | yes | 137 | 686 | 0.3 |
| ✅ | L-BFGS | maxIterations=1000, historySize=10 | 2.00e-16 | (−5.14e-9, 2.29e-8) | 2.00e-16 | 2.35e-8 | yes | 1 | 10 | 0.0 |

---

## Himmelblau

- **Dim:** 2 | **Type:** Multimodal (4 minima)
- **Known min:** f(3, 2) = 0 (4 equivalent minima)
- **x0:** (0, 0)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ✅ | Nelder-Mead | maxIterations=10000 | 1.43e-10 | (3.0000, 2.0000) | 1.43e-10 | 1.90e-6 | yes | 86 | 168 | 0.1 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 5.50e-13 | (3.5844, −1.8481) | 5.50e-13 | 5.13e-7 | yes | 130 | 6550 | 0.5 |
| ✅ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 1.03e-12 | (3.0000, 2.0000) | 1.03e-12 | 1.72e-7 | yes | 260 | 1301 | 0.1 |
| ✅ | Adam | maxIterations=10000, learningRate=0.1 | 3.76e-9 | (3.0000, 2.0000) | 3.76e-9 | 1.65e-5 | yes | 232 | 1161 | 0.1 |
| ✅ | L-BFGS | maxIterations=1000, historySize=10 | 6.19e-12 | (3.0000, 2.0000) | 6.19e-12 | 3.88e-7 | yes | 9 | 57 | 0.0 |

---

## Three-Hump Camel

- **Dim:** 2 | **Type:** Multimodal (3 local minima)
- **Known min:** f(0, 0) = 0
- **x0:** (2, −1)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ⚠️ | Nelder-Mead | maxIterations=10000 | 0.298638 | (1.7476, −0.8738) | 0.298638 | 1.95382 | yes | 33 | 67 | 0.1 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 1.84e-13 | (1.51e-8, 4.20e-7) | 1.84e-13 | 4.21e-7 | yes | 127 | 6400 | 0.7 |
| ⚠️ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 0.298638 | (1.7476, −0.8739) | 0.298638 | 1.95389 | yes | 289 | 1446 | 0.1 |
| ⚠️ | Adam | maxIterations=10000, learningRate=0.1 | 0.298638 | (1.7476, −0.8738) | 0.298638 | 1.95382 | yes | 169 | 846 | 0.1 |
| ⚠️ | L-BFGS | maxIterations=1000, historySize=10 | 0.298638 | (1.7476, −0.8738) | 0.298638 | 1.95382 | yes | 9 | 55 | 0.0 |

---

## Rastrigin

- **Dim:** 2 | **Type:** Multimodal (highly)
- **Known min:** f(0, 0) = 0
- **x0:** (2.5, −3.5)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ⚠️ | Nelder-Mead | maxIterations=10000 | 24.873845 | (2.9849, −3.9798) | 24.873845 | 4.97474 | yes | 42 | 85 | 0.1 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 4.26e-14 | (−4.96e-9, −1.40e-8) | 4.26e-14 | 1.48e-8 | yes | 141 | 7100 | 1.1 |
| ⚠️ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 12.934432 | (1.9899, −2.9849) | 12.934432 | 3.58736 | yes | 223 | 1116 | 0.2 |
| ⚠️ | Adam | maxIterations=10000, learningRate=0.1 | 12.934432 | (1.9899, −2.9849) | 12.934432 | 3.58736 | yes | 243 | 1216 | 0.1 |
| ✅ | L-BFGS | maxIterations=1000, historySize=10 | 3.73e-14 | (6.31e-9, −1.24e-8) | 3.73e-14 | 1.39e-8 | yes | 1 | 11 | 0.0 |

---

## Ackley

- **Dim:** 2 | **Type:** Multimodal
- **Known min:** f(0, 0) = 0
- **x0:** (2, −2)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ⚠️ | Nelder-Mead | maxIterations=10000 | 6.559645 | (1.9745, −1.9744) | 6.559645 | 2.79229 | yes | 30 | 60 | 0.1 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 6.76e-11 | (−8.73e-12, 2.23e-11) | 6.76e-11 | 2.39e-11 | yes | 198 | 9950 | 1.0 |
| ⚠️ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 6.559645 | (1.9745, −1.9745) | 6.559645 | 2.79230 | yes | 127 | 636 | 0.1 |
| ⚠️ | Adam | maxIterations=10000, learningRate=0.1 | 6.559645 | (1.9745, −1.9745) | 6.559645 | 2.79230 | yes | 137 | 686 | 0.1 |
| ✅ | L-BFGS | maxIterations=1000, historySize=10 | 8.39e-12 | (2.10e-12, −2.10e-12) | 8.39e-12 | 2.97e-12 | yes | 18 | 107 | 0.0 |

---

## Levi N.13

- **Dim:** 2 | **Type:** Multimodal
- **Known min:** f(1, 1) = 0
- **x0:** (−4, 5)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ⚠️ | Nelder-Mead | maxIterations=10000 | 15.974616 | (1.0000, 4.9936) | 15.974616 | 3.99364 | yes | 60 | 118 | 0.1 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 9.78e-15 | (1.0000, 1.0000) | 9.78e-15 | 7.66e-8 | yes | 129 | 6500 | 0.7 |
| ⚠️ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 0.109874 | (0.6704, 1.0000) | 0.109874 | 0.32962 | yes | 284 | 1421 | 0.1 |
| ⚠️ | Adam | maxIterations=10000, learningRate=0.1 | 16.966524 | (0.0112, 4.9944) | 16.966524 | 4.11497 | yes | 253 | 1266 | 0.1 |
| ✅ | L-BFGS | maxIterations=1000, historySize=10 | 1.44e-14 | (1.0000, 1.0000) | 1.44e-14 | 1.69e-8 | yes | 1 | 11 | 0.0 |

---

## Griewank

- **Dim:** 2 | **Type:** Multimodal
- **Known min:** f(0, 0) = 0
- **x0:** (100, −200)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ⚠️ | Nelder-Mead | maxIterations=10000 | 3.740382 | (119.3555, −26.5625) | 3.740382 | 122.275 | yes | 17 | 33 | 0.1 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 7.01e-13 | (−1.07e-6, 7.27e-7) | 7.01e-13 | 1.29e-6 | yes | 238 | 11950 | 1.5 |
| ⚠️ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 12.352951 | (97.3404, −199.7305) | 12.352951 | 222.188 | yes | 978 | 4891 | 0.2 |
| ⚠️ | Adam | maxIterations=10000, learningRate=0.1 | 12.352950 | (97.3402, −199.7291) | 12.352950 | 222.186 | yes | 228 | 1141 | 0.1 |
| ⚠️ | L-BFGS | maxIterations=1000, historySize=10 | 12.352950 | (97.3402, −199.7291) | 12.352950 | 222.186 | yes | 7 | 40 | 0.0 |

---

## Styblinski-Tang

- **Dim:** 2 | **Type:** Multimodal
- **Known min:** f(−2.9035, −2.9035) = −78.33234
- **x0:** (0, 0)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ⚠️ | Nelder-Mead | maxIterations=10000 | −64.195612 | (2.7468, −2.9035) | 14.136728 | 5.65034 | yes | 90 | 177 | 0.1 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | −78.332331 | (−2.9035, −2.9035) | 8.59e-6 | 1.26e-7 | yes | 124 | 6250 | 0.6 |
| ✅ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | −78.332331 | (−2.9035, −2.9035) | 8.59e-6 | 4.05e-7 | yes | 227 | 1136 | 0.3 |
| ⚠️ | Adam | maxIterations=10000, learningRate=0.1 | −78.332150 | (−2.9058, −2.9058) | 1.90e-4 | 0.00324 | yes | 78 | 391 | 0.1 |
| ✅ | L-BFGS | maxIterations=1000, historySize=10 | −78.332331 | (−2.9035, −2.9035) | 8.59e-6 | 4.08e-8 | yes | 7 | 44 | 0.0 |

---

## Easom

- **Dim:** 2 | **Type:** Nearly flat, narrow peak
- **Known min:** f(pi, pi) = −1
- **x0:** (1, 1)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ⚠️ | Nelder-Mead | maxIterations=10000 | −8.11e-5 | (1.3054, 1.3048) | 0.999919 | 2.59716 | yes | 20 | 40 | 0.1 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | −1.000000 | (3.1416, 3.1416) | 3.75e-13 | 5.00e-7 | yes | 112 | 5650 | 0.7 |
| ⚠️ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | −3.03e-5 | (1.0000, 1.0000) | 0.999970 | 3.02862 | yes | 50 | 251 | 0.0 |
| ⚠️ | Adam | maxIterations=10000, learningRate=0.1 | −8.11e-5 | (1.3047, 1.3047) | 0.999919 | 2.59769 | yes | 58 | 291 | 0.0 |
| ❌ | L-BFGS | maxIterations=1000, historySize=10 | −4.81e-5 | (1.0940, 1.0940) | 0.999952 | 2.89571 | no | 1000 | 5005 | 0.8 |

---

## Goldstein-Price

- **Dim:** 2 | **Type:** Multimodal
- **Known min:** f(0, −1) = 3
- **x0:** (0, 0)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ⚠️ | Nelder-Mead | maxIterations=10000 | 30.000000 | (−0.6000, −0.4000) | 27.000000 | 0.84853 | yes | 75 | 143 | 0.1 |
| ✅ | PSO | maxIterations=10000, swarmSize=50, seed=42 | 3.000000 | (3.62e-8, −1.0000) | 3.85e-13 | 3.71e-8 | yes | 136 | 6850 | 0.6 |
| ⚠️ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | 334.333604 | (−0.7200, −0.7200) | 331.333604 | 0.77253 | yes | 5 | 30 | 0.0 |
| ⚠️ | Adam | maxIterations=10000, learningRate=0.1 | 30.000000 | (−0.6000, −0.4000) | 27.000000 | 0.84853 | yes | 255 | 1276 | 0.1 |
| ⚠️ | L-BFGS | maxIterations=1000, historySize=10 | 30.000000 | (−0.6000, −0.4000) | 27.000000 | 0.84853 | yes | 11 | 71 | 0.0 |

---

## McCormick

- **Dim:** 2 | **Type:** Multimodal
- **Known min:** f(−0.54719, −1.54719) = −1.9133
- **x0:** (0, 0)

|   | Method | Settings | Found Value | Found Point | \|Error\| | Dist to Opt | Conv | Iters | Fn Evals | Time (ms) |
|---|--------|----------|-------------|-------------|-----------|-------------|------|-------|----------|-----------|
| ✅ | Nelder-Mead | maxIterations=10000 | −1.913223 | (−0.5472, −1.5472) | 7.70e-5 | 9.79e-6 | yes | 79 | 153 | 0.1 |
| ❌ | PSO | maxIterations=10000, swarmSize=50, seed=42 | −38424.037 | (−38424.0, −38425.1) | 38422.123 | 54339.046 | no | 10000 | 500050 | 29.5 |
| ✅ | GradientDescent | maxIterations=10000, learningRate=0.001, momentum=0.9 | −1.913223 | (−0.5471, −1.5471) | 7.71e-5 | 1.62e-4 | yes | 417 | 2086 | 0.1 |
| ✅ | Adam | maxIterations=10000, learningRate=0.1 | −1.913223 | (−0.5472, −1.5472) | 7.70e-5 | 1.42e-5 | yes | 193 | 966 | 0.1 |
| ✅ | L-BFGS | maxIterations=1000, historySize=10 | −1.913223 | (−0.5472, −1.5472) | 7.70e-5 | 1.07e-5 | yes | 7 | 41 | 0.0 |

---

## Notes

- All optimizers use default hyperparameters (see Settings column)
- **Fn Evals** = total calls to objective function (via counter wrapper)
- **\|Error\|** = \|found_value - known_minimum\|
- **Dist to Opt** = \|\|found_point - known_point\|\|₂
- PSO uses seed=42 for reproducibility
- Himmelblau: Dist to Opt uses nearest of 4 known minima
- McCormick: PSO diverges without box constraints (unbounded domain)
- L-BFGS dominates smooth unimodal problems (Sphere, Matyas, Rosenbrock, Beale, Booth)
  with very few iterations; struggles on Easom (nearly-flat landscape starves the
  finite-difference gradient)
- **L-BFGS on Rastrigin / Lévi N.13 — caveat.** The global-minimum result (1 iteration,
  error ~10⁻¹⁴) is an artefact of the starting point, not a general property of the
  algorithm. At x₀ = (2.5, −3.5) for Rastrigin and x₀ = (−4, 5) for Lévi N.13, every
  `sin(2πxᵢ)` / `sin(3πxᵢ)` term vanishes, so the gradient reduces to its quadratic part
  pointing at the origin; Armijo backtracking halves a unit step to land exactly at the
  global optimum. Perturb x₀ by a non-integer offset (e.g. (2.6, −3.4)) and L-BFGS will
  fall into the nearest local well just like GD / Adam — it is a local optimizer. For
  multimodal problems use PSO or a multi-start wrapper.
