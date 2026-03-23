# Unconstrained Optimization Benchmark Functions

Source: [Test functions for optimization — Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization)

---

## 1. Sphere

| Property | Value |
|---|---|
| Formula | f(x) = Σ xᵢ² |
| Dimension | n (tested: 2) |
| Type | Unimodal, convex, differentiable |
| Search domain | −∞ ≤ xᵢ ≤ ∞ |
| Global minimum | f(0, …, 0) = 0 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 2. Rosenbrock

| Property | Value |
|---|---|
| Formula | f(x) = Σᵢ₌₁ⁿ⁻¹ [100·(xᵢ₊₁ − xᵢ²)² + (1 − xᵢ)²] |
| Dimension | n (tested: 2) |
| Type | Unimodal, non-convex, differentiable |
| Search domain | −∞ ≤ xᵢ ≤ ∞ |
| Global minimum | f(1, …, 1) = 0 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Rosenbrock_function) |

---

## 3. Beale

| Property | Value |
|---|---|
| Formula | f(x,y) = (1.5 − x + xy)² + (2.25 − x + xy²)² + (2.625 − x + xy³)² |
| Dimension | 2 |
| Type | Unimodal, non-convex, differentiable |
| Search domain | −4.5 ≤ x, y ≤ 4.5 |
| Global minimum | f(3, 0.5) = 0 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 4. Booth

| Property | Value |
|---|---|
| Formula | f(x,y) = (x + 2y − 7)² + (2x + y − 5)² |
| Dimension | 2 |
| Type | Unimodal, convex, differentiable |
| Search domain | −10 ≤ x, y ≤ 10 |
| Global minimum | f(1, 3) = 0 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 5. Matyas

| Property | Value |
|---|---|
| Formula | f(x,y) = 0.26·(x² + y²) − 0.48·x·y |
| Dimension | 2 |
| Type | Unimodal, convex, differentiable |
| Search domain | −10 ≤ x, y ≤ 10 |
| Global minimum | f(0, 0) = 0 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 6. Himmelblau

| Property | Value |
|---|---|
| Formula | f(x,y) = (x² + y − 11)² + (x + y² − 7)² |
| Dimension | 2 |
| Type | Multimodal (4 identical minima), differentiable |
| Search domain | −5 ≤ x, y ≤ 5 |
| Global minima | f(3, 2) = 0; f(−2.805118, 3.131312) = 0; f(−3.779310, −3.283186) = 0; f(3.584428, −1.848126) = 0 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Himmelblau%27s_function) |

---

## 7. Three-Hump Camel

| Property | Value |
|---|---|
| Formula | f(x,y) = 2x² − 1.05x⁴ + x⁶/6 + xy + y² |
| Dimension | 2 |
| Type | Multimodal (3 local minima), differentiable |
| Search domain | −5 ≤ x, y ≤ 5 |
| Global minimum | f(0, 0) = 0 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 8. Rastrigin

| Property | Value |
|---|---|
| Formula | f(x) = 10n + Σ [xᵢ² − 10·cos(2πxᵢ)] |
| Dimension | n (tested: 2) |
| Type | Multimodal (highly), differentiable |
| Search domain | −5.12 ≤ xᵢ ≤ 5.12 |
| Global minimum | f(0, …, 0) = 0 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Rastrigin_function) |

---

## 9. Ackley

| Property | Value |
|---|---|
| Formula | f(x,y) = −20·exp[−0.2·√(0.5·(x²+y²))] − exp[0.5·(cos(2πx)+cos(2πy))] + e + 20 |
| Dimension | 2 |
| Type | Multimodal, differentiable |
| Search domain | −5 ≤ x, y ≤ 5 |
| Global minimum | f(0, 0) = 0 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Ackley_function) |

---

## 10. Lévi N.13

| Property | Value |
|---|---|
| Formula | f(x,y) = sin²(3πx) + (x−1)²·(1+sin²(3πy)) + (y−1)²·(1+sin²(2πy)) |
| Dimension | 2 |
| Type | Multimodal, differentiable |
| Search domain | −10 ≤ x, y ≤ 10 |
| Global minimum | f(1, 1) = 0 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 11. Griewank

| Property | Value |
|---|---|
| Formula | f(x) = 1 + (1/4000)·Σ xᵢ² − Π cos(xᵢ/√i) |
| Dimension | n (tested: 2) |
| Type | Multimodal, differentiable |
| Search domain | −600 ≤ xᵢ ≤ 600 |
| Global minimum | f(0, …, 0) = 0 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Griewank_function) |

---

## 12. Styblinski-Tang

| Property | Value |
|---|---|
| Formula | f(x) = (1/2)·Σ (xᵢ⁴ − 16xᵢ² + 5xᵢ) |
| Dimension | n (tested: 2) |
| Type | Multimodal, differentiable |
| Search domain | −5 ≤ xᵢ ≤ 5 |
| Global minimum | f(−2.903534, …, −2.903534) ≈ −39.16617·n |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 13. Easom

| Property | Value |
|---|---|
| Formula | f(x,y) = −cos(x)·cos(y)·exp(−((x−π)²+(y−π)²)) |
| Dimension | 2 |
| Type | Unimodal (nearly flat everywhere except near (π,π)), differentiable |
| Search domain | −100 ≤ x, y ≤ 100 |
| Global minimum | f(π, π) = −1 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 14. Goldstein-Price

| Property | Value |
|---|---|
| Formula | f(x,y) = [1+(x+y+1)²·(19−14x+3x²−14y+6xy+3y²)] · [30+(2x−3y)²·(18−32x+12x²+48y−36xy+27y²)] |
| Dimension | 2 |
| Type | Multimodal, differentiable |
| Search domain | −2 ≤ x, y ≤ 2 |
| Global minimum | f(0, −1) = 3 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |

---

## 15. McCormick

| Property | Value |
|---|---|
| Formula | f(x,y) = sin(x+y) + (x−y)² − 1.5x + 2.5y + 1 |
| Dimension | 2 |
| Type | Multimodal, differentiable |
| Search domain | −1.5 ≤ x ≤ 4, −3 ≤ y ≤ 4 |
| Global minimum | f(−0.54719, −1.54719) ≈ −1.9133 |
| Reference | [Wikipedia](https://en.wikipedia.org/wiki/Test_functions_for_optimization) |
