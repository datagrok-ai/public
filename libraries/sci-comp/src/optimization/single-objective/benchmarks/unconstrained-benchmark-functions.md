# Unconstrained Optimization Benchmark Functions

Source: [Test functions for optimization — Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization)

---

## 1. Sphere

| Property | Value |
|---|---|
| Formula | $f(\mathbf{x}) = \sum_{i=1}^{n} x_i^2$ |
| Dimension | $n$ (tested: 2) |
| Type | Unimodal, convex, differentiable |
| Search domain | $-\infty \le x_i \le \infty$ |
| Global minimum | $f(0, \ldots, 0) = 0$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 2. Rosenbrock

| Property | Value |
|---|---|
| Formula | $f(\mathbf{x}) = \sum_{i=1}^{n-1} \left[ 100 \left( x_{i+1} - x_i^2 \right)^2 + (1 - x_i)^2 \right]$ |
| Dimension | $n$ (tested: 2) |
| Type | Unimodal, non-convex, differentiable |
| Search domain | $-\infty \le x_i \le \infty$ |
| Global minimum | $f(1, \ldots, 1) = 0$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Rosenbrock_function) |

---

## 3. Beale

| Property | Value |
|---|---|
| Formula | $f(x, y) = (1.5 - x + xy)^2 + (2.25 - x + xy^2)^2 + (2.625 - x + xy^3)^2$ |
| Dimension | 2 |
| Type | Unimodal, non-convex, differentiable |
| Search domain | $-4.5 \le x, y \le 4.5$ |
| Global minimum | $f(3,\; 0.5) = 0$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 4. Booth

| Property | Value |
|---|---|
| Formula | $f(x, y) = (x + 2y - 7)^2 + (2x + y - 5)^2$ |
| Dimension | 2 |
| Type | Unimodal, convex, differentiable |
| Search domain | $-10 \le x, y \le 10$ |
| Global minimum | $f(1, 3) = 0$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 5. Matyas

| Property | Value |
|---|---|
| Formula | $f(x, y) = 0.26(x^2 + y^2) - 0.48xy$ |
| Dimension | 2 |
| Type | Unimodal, convex, differentiable |
| Search domain | $-10 \le x, y \le 10$ |
| Global minimum | $f(0, 0) = 0$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 6. Himmelblau

| Property | Value |
|---|---|
| Formula | $f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2$ |
| Dimension | 2 |
| Type | Multimodal (4 identical minima), differentiable |
| Search domain | $-5 \le x, y \le 5$ |
| Global minima | $f(3, 2) = 0$; $f(-2.805118,\; 3.131312) = 0$; $f(-3.779310,\; -3.283186) = 0$; $f(3.584428,\; -1.848126) = 0$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Himmelblau%27s_function) |

---

## 7. Three-Hump Camel

| Property | Value |
|---|---|
| Formula | $f(x, y) = 2x^2 - 1.05x^4 + \dfrac{x^6}{6} + xy + y^2$ |
| Dimension | 2 |
| Type | Multimodal (3 local minima), differentiable |
| Search domain | $-5 \le x, y \le 5$ |
| Global minimum | $f(0, 0) = 0$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 8. Rastrigin

| Property | Value |
|---|---|
| Formula | $f(\mathbf{x}) = 10n + \sum_{i=1}^{n} \left[ x_i^2 - 10 \cos(2\pi x_i) \right]$ |
| Dimension | $n$ (tested: 2) |
| Type | Multimodal (highly), differentiable |
| Search domain | $-5.12 \le x_i \le 5.12$ |
| Global minimum | $f(0, \ldots, 0) = 0$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Rastrigin_function) |

---

## 9. Ackley

| Property | Value |
|---|---|
| Formula | $f(x, y) = -20 \exp\!\left[ -0.2\sqrt{0.5(x^2 + y^2)} \right] - \exp\!\left[ 0.5(\cos 2\pi x + \cos 2\pi y) \right] + e + 20$ |
| Dimension | 2 |
| Type | Multimodal, differentiable |
| Search domain | $-5 \le x, y \le 5$ |
| Global minimum | $f(0, 0) = 0$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Ackley_function) |

---

## 10. Lévi N.13

| Property | Value |
|---|---|
| Formula | $f(x, y) = \sin^2(3\pi x) + (x - 1)^2 \left[1 + \sin^2(3\pi y)\right] + (y - 1)^2 \left[1 + \sin^2(2\pi y)\right]$ |
| Dimension | 2 |
| Type | Multimodal, differentiable |
| Search domain | $-10 \le x, y \le 10$ |
| Global minimum | $f(1, 1) = 0$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 11. Griewank

| Property | Value |
|---|---|
| Formula | $f(\mathbf{x}) = 1 + \dfrac{1}{4000} \sum_{i=1}^{n} x_i^2 - \prod_{i=1}^{n} \cos\!\left(\dfrac{x_i}{\sqrt{i}}\right)$ |
| Dimension | $n$ (tested: 2) |
| Type | Multimodal, differentiable |
| Search domain | $-600 \le x_i \le 600$ |
| Global minimum | $f(0, \ldots, 0) = 0$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Griewank_function) |

---

## 12. Styblinski-Tang

| Property | Value |
|---|---|
| Formula | $f(\mathbf{x}) = \dfrac{1}{2} \sum_{i=1}^{n} \left( x_i^4 - 16x_i^2 + 5x_i \right)$ |
| Dimension | $n$ (tested: 2) |
| Type | Multimodal, differentiable |
| Search domain | $-5 \le x_i \le 5$ |
| Global minimum | $f(-2.903534, \ldots, -2.903534) \approx -39.16617 \cdot n$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 13. Easom

| Property | Value |
|---|---|
| Formula | $f(x, y) = -\cos(x) \cdot \cos(y) \cdot \exp\!\left(-\left((x - \pi)^2 + (y - \pi)^2\right)\right)$ |
| Dimension | 2 |
| Type | Unimodal (nearly flat everywhere except near $(\pi, \pi)$), differentiable |
| Search domain | $-100 \le x, y \le 100$ |
| Global minimum | $f(\pi, \pi) = -1$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 14. Goldstein-Price

| Property | Value |
|---|---|
| Formula | $f(x, y) = \left[1 + (x + y + 1)^2 (19 - 14x + 3x^2 - 14y + 6xy + 3y^2)\right] \cdot \left[30 + (2x - 3y)^2 (18 - 32x + 12x^2 + 48y - 36xy + 27y^2)\right]$ |
| Dimension | 2 |
| Type | Multimodal, differentiable |
| Search domain | $-2 \le x, y \le 2$ |
| Global minimum | $f(0, -1) = 3$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 15. McCormick

| Property | Value |
|---|---|
| Formula | $f(x, y) = \sin(x + y) + (x - y)^2 - 1.5x + 2.5y + 1$ |
| Dimension | 2 |
| Type | Multimodal, differentiable |
| Search domain | $-1.5 \le x \le 4$, $-3 \le y \le 4$ |
| Global minimum | $f(-0.54719,\; -1.54719) \approx -1.9133$ |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |
