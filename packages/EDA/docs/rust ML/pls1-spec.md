# Спецификация: PLS1-регрессия (без дефляции X) на Rust + WASM

**Версия:** 1.0
**Цель:** реализовать с нуля (на голом `Vec<f64>`) регрессию методом частных наименьших квадратов **PLS1** (один отклик `y`) по алгоритму Indahl (Algorithm 2, DOI 10.1002/cem.2589) — в варианте **без дефляции `X`**: рекуррентность весов `w` плюс ортогонализация scores `t`. Rust-ядро универсально; WASM — один из бэкендов. За основу взята референс-реализация на Eigen (`pls.cpp` / `PLS.h`) и её TS-обвязка (`pls-ml.ts`, `pls-constants.ts`).

> **Общий контракт с `linreg-spec` и `pca-nipals-spec` (вариант B).** Данные подаются как **список колонок** (`&[&[f64]]`). Центрирование и стандартизация выполняются на стороне **JS/TS**: в ядро приходят уже подготовленные предикторы `X` и отклик `y`. Статистики (μ/σ предикторов и отклика) передаются **опционально** и используются только для раскрутки регрессионных коэффициентов и вычисления свободного члена в исходном пространстве. Та же система ошибок и та же WASM-граница (плоская конкатенация колонок + `n_rows`).

> **Отличие от референса.** В `pls.cpp` центрирование и нормировка делаются внутри функции (вариант A). Здесь ядро их **не** делает (вариант B); сам PLS1-цикл, рекуррентность `w` и формула коэффициентов сохранены без изменений.

---

## 1. Область применения и ограничения

В объёме (in scope):
- PLS1: один скалярный отклик `y`, `A` латентных компонент.
- Алгоритм Indahl без дефляции `X`: веса `W`, X-scores `T`, X-loadings `P`, Y-loadings `q`, Y-scores `U`.
- Регрессионные коэффициенты `b = W·(PᵀW)⁻¹·q`, свободный член (bias), предсказание.
- Накопленная объяснённая дисперсия отклика по компонентам.
- WASM-обёртка через `wasm-bindgen`; порядок выходов согласован с `WASM_OUTPUT_IDX` из референса.

Вне объёма (out of scope) в v1.0:
- **PLS2** (многомерный отклик) — отдельная версия.
- Центрирование/стандартизация внутри ядра (делается на JS, вариант B).
- Квадратичный режим (`isQuadratic`): добавление квадратов признаков — это инжиниринг признаков на стороне JS (добавить колонки `xⱼ²` в список перед вызовом), ядро его не выполняет.
- Кросс-валидация и автоподбор числа компонент.
- Обработка пропусков: вход обязан быть полным.

---

## 2. Математическая модель

### 2.1. Постановка

Дано: матрица предикторов `X` размера `n × m` (n наблюдений, m признаков-колонок) и отклик `y` длины `n`, **оба уже центрированы и стандартизованы по столбцам** (на стороне JS, вариант B). PLS1 строит `A` латентных компонент, одновременно (а) захватывающих дисперсию `X` и (б) максимизирующих ковариацию с `y`, и выражает регрессию:

```
ŷ = b₀ + Σⱼ bⱼ·xⱼ
```

где `b` — вектор регрессионных коэффициентов в исходном пространстве признаков, `b₀` — свободный член.

### 2.2. Конвенция нормировки (для согласования раскрутки)

Референс использует популяционную нормировку: `σⱼ = ‖X_centered[:,j]‖ / √n`, `σ_y = √(‖y_centered‖²/n)`. После стандартизации `‖y‖² = n`. Эта же конвенция должна применяться на стороне JS, чтобы раскрутка коэффициентов (2.4) была корректной. Если JS использует выборочную нормировку (`/(n−1)`), формулы раскрутки нужно согласовать — это инвариант контракта.

### 2.3. Накопленная объяснённая дисперсия

Так как `‖y‖² = n`, а scores `t` ортонормированы:

```
expl_var(a) = Σ_{k≤a} q(k)² / n        # накопленная доля дисперсии y, ∈ [0, 1]
```

(источник формулы — тот же paper 10.1002/cem.2589; в `pls-ml.ts` это `computeExplVars`).

### 2.4. Раскрутка коэффициентов и свободный член

Коэффициенты `b` вычисляются в стандартизованном пространстве, затем раскручиваются (как в `pls.cpp` и `getRegrCoeffs` из `pls-ml.ts`):

```
b_origⱼ = bⱼ · (σ_y / σ_xⱼ)
b₀       = μ_y − Σⱼ b_origⱼ · μ_xⱼ
```

Без заданных статистик `coefficients()` возвращает `b` в стандартизованном пространстве и `b₀ = 0`.

---

## 3. Алгоритм PLS1 (Indahl, Algorithm 2, без дефляции X)

`X`, `y` — уже стандартизованные. Матрицы: `W` (m×A), `P` (m×A), `T` (n×A), вектор `q` (A). Ядро `X` **не мутирует** (дефляции нет).

### 3.1. Инициализация (компонента 0)

```
w = Xᵀ y
normV[0] = ‖w‖;            если 0 → MethodError
w = w / normV[0];          W[:,0] = w
t = X w
normTau[0] = ‖t‖;          если 0 → MethodError
t = t / normTau[0];        T[:,0] = t
p = Xᵀ t;                  P[:,0] = p
q[0] = tᵀ y
```

### 3.2. Итерации (компоненты a = 1 .. A−1)

```
# рекуррентность весов вместо дефляции X:
w = normV[a-1] · ( w − p / normTau[a-1] )
normV[a] = ‖w‖;            если 0 → MethodError
w = w / normV[a];          W[:,a] = w
t = X w
# ортогонализация scores против всех предыдущих (Gram–Schmidt, T ортонормирована):
t = t − T[:, 0..a] · ( T[:, 0..a]ᵀ t )
normTau[a] = ‖t‖;          если 0 → MethodError
t = t / normTau[a];        T[:,a] = t
p = Xᵀ t;                  P[:,a] = p
q[a] = tᵀ y
```

### 3.3. Постобработка

```
U = y · qᵀ / ‖q‖²                       # Y-scores (n×A)
H = Pᵀ W                                  # A×A
если det(H) ≈ 0 → MethodError
b = W · H⁻¹ · q                           # коэффициенты в стандартизованном пространстве
# раскрутка (если заданы статистики, см. 2.4):
для j: b_origⱼ = bⱼ · σ_y / σ_xⱼ
b₀ = μ_y − Σⱼ b_origⱼ · μ_xⱼ
prediction = b₀ + X_orig · b_orig         # на исходных (нестандартизованных) данных
expl_var[a] = накопленная Σ q²/n          # см. 2.3
```

> `H⁻¹` — обращение малой матрицы `A×A` (A — число компонент, по умолчанию 3). Реализуется методом Гаусса–Жордана на `Vec<f64>`; проверка `det ≈ 0` (вырожденность) → `MethodError`.

### 3.4. Численная устойчивость

- Деления на `‖w‖`, `‖t‖`, `‖q‖²`, на `det(H)` защищены проверкой на близость к нулю (`< eps`) → `MethodError` (соответствует `METHOD_ERROR` референса).
- Ортогонализация `t` против `T[:,0..a]` компенсирует отсутствие дефляции; при больших `A` можно повторить шаг Gram–Schmidt дважды (опция `reorthogonalize`, по умолчанию `false`).
- NaN/Inf на входе → `Err(NonFiniteValue)`.

---

## 4. Архитектура Rust-ядра

Общий стиль с предыдущими крейтами; чистый `Vec<f64>`, `std` ради `Vec`/`sqrt`.

### 4.1. Структура крейта

```
crates/
  pls-core/             # чистое ядро, без wasm
    src/
      lib.rs
      matrix.rs         # операции по колонкам (dot/axpy/norm), малый Gauss–Jordan для H⁻¹
      stats.rs          # FeatureStats: μ/σ предикторов + μ/σ отклика, раскрутка b
      pls1.rs           # алгоритм Indahl (без дефляции): W, T, P, q, U
      model.rs          # Pls1: fit / predict / coefficients / scores / loadings
      metrics.rs        # explained variance (накопленная), reconstruction
      error.rs          # типы ошибок (как в linreg + MethodError)
    tests/
      integration.rs    # см. раздел 7
  pls-wasm/             # обёртка wasm-bindgen, зависит от pls-core
    src/lib.rs
```

### 4.2. Представление данных — список колонок

Как в общих контрактах: `x_columns: &[&[f64]]`, `x_columns[j]` — предиктор `j` длиной `n_rows`; `m = x_columns.len()`. Отклик `y: &[f64]` длины `n_rows`. Инварианты: непустой список, равные длины колонок, `y.len() == n_rows`, отсутствие NaN/Inf, `n_rows ≥ 2`. Валидация компонент (из `pls-constants.ts` / `pls.cpp`): `1 ≤ A ≤ m` и `A ≤ n_rows`. Нарушение → соответствующая ошибка.

Колоночное хранение естественно: `Xᵀy` и `Xᵀt` — это `m` dot-product'ов «колонка · вектор»; `Xw` — аккумулятор по колонкам `t += w[j]·X[:,j]`. `X` не мутируется (дефляции нет) — ещё удобнее для срезов.

### 4.3. Ключевые типы (эскиз API)

```rust
// error.rs
#[derive(Debug, Clone, PartialEq)]
pub enum PlsError {
    DimensionMismatch { expected: usize, got: usize },
    EmptyDataset,
    NotFitted,
    InvalidHyperparameter(String),  // в т.ч. components вне [1, min(m, n)]
    NonFiniteValue,
    MethodError,                    // деление на ноль / вырожденная H (как METHOD_ERROR)
}

pub struct Config {
    pub components: usize,     // A; default 3 (как COMPONENTS.DEFAULT), min 1
    pub reorthogonalize: bool, // повторный Gram–Schmidt для t, default false
}

impl Default for Config { /* components=3, reorthogonalize=false */ }

// stats.rs — вариант B: статистики только для раскрутки коэффициентов.
pub struct FeatureStats {
    x_means: Vec<f64>, x_stds: Vec<f64>,  // длина m
    y_mean: f64,       y_std: f64,
}
impl FeatureStats {
    pub fn new(x_means: Vec<f64>, x_stds: Vec<f64>, y_mean: f64, y_std: f64)
        -> Result<Self, PlsError>; // валидирует: длины, конечность, σ>0
}

pub struct Pls1 {
    // всё в стандартизованном пространстве:
    weights: Vec<f64>,         // W, m×A (column-major по компонентам)
    x_loadings: Vec<f64>,      // P, m×A
    t_scores: Vec<f64>,        // T, n×A (кэш обучающих)
    u_scores: Vec<f64>,        // U, n×A
    y_loadings: Vec<f64>,      // q, длина A
    coeffs_std: Vec<f64>,      // b в стандартизованном пространстве, длина m
    explained_variance: Vec<f64>, // накопленная, длина A
    n_features: usize, n_rows: usize, components: usize,
    stats: Option<FeatureStats>,
    config: Config,
}

impl Pls1 {
    pub fn new(config: Config) -> Result<Self, PlsError>;

    // Опционально: μ/σ предикторов и отклика — ТОЛЬКО для раскрутки b и bias.
    pub fn set_feature_stats(
        &mut self, x_means: Vec<f64>, x_stds: Vec<f64>, y_mean: f64, y_std: f64
    ) -> Result<(), PlsError>;

    // КОНТРАКТ: x_columns и y — УЖЕ СТАНДАРТИЗОВАННЫЕ (вариант B).
    pub fn fit(&mut self, x_columns: &[&[f64]], y: &[f64]) -> Result<(), PlsError>;

    // Предсказание на ИСХОДНЫХ (нестандартизованных) данных:
    //   ŷ = b₀ + Σ b_origⱼ·xⱼ. Требует заданных статистик (иначе работает
    //   в стандартизованном пространстве с b₀=0).
    pub fn predict(&self, x_columns: &[&[f64]]) -> Result<Vec<f64>, PlsError>;

    // (m,) коэффициенты + b₀. С μ/σ — в исходном пространстве; иначе — стандартизованные, b₀=0.
    pub fn coefficients(&self) -> (Vec<f64>, f64);

    // Геттеры (NotFitted до fit):
    pub fn t_scores(&self) -> Result<&[f64], PlsError>;   // n×A
    pub fn u_scores(&self) -> Result<&[f64], PlsError>;   // n×A
    pub fn x_loadings(&self) -> Result<&[f64], PlsError>; // m×A
    pub fn y_loadings(&self) -> Result<&[f64], PlsError>; // A
    pub fn explained_variance(&self) -> &[f64];           // накопленная, A
    pub fn components(&self) -> usize;
}
```

> Как и в предыдущих крейтах: `fit`/внутренний цикл работают над поданными (стандартизованными) данными; `set_feature_stats` влияет только на `coefficients()`/`predict()` (раскрутка). Согласованность масштаба между `fit` и `predict` — ответственность вызывающей стороны.

---

## 5. Алгоритм обучения (псевдокод)

```
fit(x_columns, y):                         # уже стандартизованные X и y
    m = x_columns.len();  n = x_columns[0].len()
    валидировать инварианты входа (4.2) → EmptyDataset / DimensionMismatch / InvalidHyperparameter
    проверить NaN/Inf → NonFiniteValue
    A = config.components
    # --- компонента 0 ---
    w = [ dot(X[:,j], y) for j in 0..m ]              # Xᵀy
    nV = ‖w‖; if nV<eps: MethodError; w /= nV; W[:,0]=w
    t = accumulate( w[j]·X[:,j] )                     # Xw
    nT = ‖t‖; if nT<eps: MethodError; t /= nT; T[:,0]=t
    p = [ dot(X[:,j], t) for j in 0..m ]              # Xᵀt
    P[:,0]=p;  q[0]=dot(t,y)
    normV[0]=nV; normTau[0]=nT
    # --- компоненты 1..A-1 ---
    для a в 1..A:
        w = normV[a-1] · ( w − p / normTau[a-1] )     # рекуррентность весов
        nV = ‖w‖; if nV<eps: MethodError; w /= nV; W[:,a]=w; normV[a]=nV
        t = accumulate( w[j]·X[:,j] )
        # Gram–Schmidt против T[:,0..a]:
        для k в 0..a: c = dot(T[:,k], t); t -= c·T[:,k]
        if reorthogonalize: повторить блок Gram–Schmidt
        nT = ‖t‖; if nT<eps: MethodError; t /= nT; T[:,a]=t; normTau[a]=nT
        p = [ dot(X[:,j], t) for j in 0..m ]; P[:,a]=p
        q[a]=dot(t,y)
    # --- постобработка ---
    U = outer(y, q) / dot(q,q)                         # n×A
    H = Pᵀ·W (A×A);  if |det(H)|<eps: MethodError
    b_std = W · H⁻¹ · q                                # Gauss–Jordan для H⁻¹
    expl = накопленная Σ q²/n
    сохранить W,P,T,U,q,b_std,expl
```

> Все горячие операции (`Xᵀy`, `Xw`, `Xᵀt`, Gram–Schmidt, дот-произведения) — проходы по колонкам/векторам с последовательным доступом. Единственная плотная линейная алгебра — обращение `H` размера `A×A` (мало).

---

## 6. WASM-бэкенд

Граница как в предыдущих крейтах: предикторы — плоская конкатенация колонок (`Float64Array`) + `n_rows`; отклик — отдельный `Float64Array`. Выходы возвращаются плоско; порядок согласован с `WASM_OUTPUT_IDX` из `pls-constants.ts` (prediction, regr_coeffs, T, U, X-loadings, Y-loadings).

```rust
// pls-wasm/src/lib.rs
use wasm_bindgen::prelude::*;
use pls_core::{Pls1, Config};

#[wasm_bindgen]
pub struct WasmPls { inner: Pls1 }

#[wasm_bindgen]
impl WasmPls {
    #[wasm_bindgen(constructor)]
    pub fn new(components: usize) -> Result<WasmPls, JsError> {
        let cfg = Config { components, ..Default::default() };
        Ok(Self { inner: Pls1::new(cfg).map_err(to_js)? })
    }

    // μ/σ предикторов (Float64Array длины m) + отклика — для раскрутки b/bias.
    pub fn set_feature_stats(&mut self, x_means: &[f64], x_stds: &[f64],
                             y_mean: f64, y_std: f64) -> Result<(), JsError> {
        self.inner.set_feature_stats(x_means.to_vec(), x_stds.to_vec(), y_mean, y_std)
            .map_err(to_js)
    }

    // flat_x: конкатенация стандартизованных колонок; y: стандартизованный отклик.
    pub fn fit(&mut self, flat_x: &[f64], n_rows: usize, y: &[f64]) -> Result<(), JsError> {
        let cols = to_columns(flat_x, n_rows).map_err(to_js)?;
        self.inner.fit(&cols, y).map_err(to_js)
    }

    pub fn predict(&self, flat_x: &[f64], n_rows: usize) -> Result<Vec<f64>, JsError> {
        let cols = to_columns(flat_x, n_rows).map_err(to_js)?;
        self.inner.predict(&cols).map_err(to_js)
    }

    // Плоские выгрузки (формы знает JS): b(m), T(n×A), U(n×A), P(m×A), q(A), expl_var(A).
    pub fn regression_coefficients(&self) -> Vec<f64> { self.inner.coefficients().0 }
    pub fn bias(&self) -> f64 { self.inner.coefficients().1 }
    pub fn t_scores(&self) -> Result<Vec<f64>, JsError> { self.inner.t_scores().map(<[f64]>::to_vec).map_err(to_js) }
    pub fn u_scores(&self) -> Result<Vec<f64>, JsError> { self.inner.u_scores().map(<[f64]>::to_vec).map_err(to_js) }
    pub fn x_loadings(&self) -> Result<Vec<f64>, JsError> { self.inner.x_loadings().map(<[f64]>::to_vec).map_err(to_js) }
    pub fn y_loadings(&self) -> Result<Vec<f64>, JsError> { self.inner.y_loadings().map(<[f64]>::to_vec).map_err(to_js) }
    pub fn explained_variance(&self) -> Vec<f64> { self.inner.explained_variance().to_vec() }
}

fn to_columns(flat: &[f64], n_rows: usize) -> Result<Vec<&[f64]>, PlsError> {
    if n_rows == 0 || flat.len() % n_rows != 0 { /* -> DimensionMismatch */ }
    Ok(flat.chunks_exact(n_rows).collect())
}
```

Сборка: `wasm-pack build crates/pls-wasm --target web`.

---

## 7. Валидационные тесты

Стиль как в предыдущих документах: (а) unit, (б) числовые, (в) property, (г) WASM. Эталоны — аналитические/из свойств PLS, плюс сверка с эталонным запуском референса (см. 7.5).

### 7.1. Unit-тесты компонентов

| # | Что проверяем | Ожидание |
|---|---|---|
| U1 | `matrix`: dot/axpy/norm | совпадение с ручным расчётом |
| U2 | `H⁻¹` (Gauss–Jordan) на известной малой матрице | `H·H⁻¹ ≈ I` (tol 1e-9) |
| U3 | `H⁻¹` на вырожденной `H` | `Err(MethodError)` |
| U4 | `FeatureStats::new` с σ>0 / с σ≤0 или NaN | OK / `Err(InvalidHyperparameter|NonFiniteValue)` |
| U5 | Валидация компонент: `A=0`, `A>m`, `A>n` | `Err(InvalidHyperparameter)` |
| U6 | NaN/Inf на входе | `Err(NonFiniteValue)` |
| U7 | Колонки разной длины / `y.len()≠n` | `Err(DimensionMismatch)` |
| U8 | Пустой список / `n<2` | `Err(EmptyDataset)` |
| U9 | Геттеры до `fit` | `Err(NotFitted)` |
| U10 | Деление на ноль (`‖w‖`/`‖t‖`/`‖q‖²`≈0) на вырожденных данных | `Err(MethodError)` |

### 7.2. Корректность PLS1 (числовые тесты)

| # | Сценарий | Ожидание |
|---|---|---|
| C1 | **Ортонормированность scores**: `TᵀT ≈ I` | внедиагональ < 1e-8, диагональ ≈ 1 |
| C2 | **Ортонормированность весов**: `WᵀW ≈ I` | то же |
| C3 | **Эквивалентность обычной регрессии при `A=m`**: PLS с числом компонент = ранг X | коэффициенты совпадают с OLS на стандартизованных данных (tol 1e-6) |
| C4 | **Точное восстановление линейной связи без шума**: `y = Xβ`, `A` достаточно | `predict` на исходных данных ≈ `y` (tol 1e-6) |
| C5 | **Раскрутка коэффициентов**: заданы μ/σ | `b_origⱼ = b_stdⱼ·σ_y/σ_xⱼ`, `b₀ = μ_y − Σ b_origⱼμ_xⱼ` (tol 1e-9) |
| C6 | **Накопленная объяснённая дисперсия** монотонна и ∈[0,1] | `expl[a] ≤ expl[a+1]`, `0 ≤ expl ≤ 1` |
| C7 | **Согласованность с референсом** (`pls.cpp`): один датасет | `T`,`U`,`P`,`q`,`b`,`prediction` совпадают с эталоном (tol 1e-5, f32-референс) |
| C8 | **Y-scores**: `U = y·qᵀ/‖q‖²`, проверка `UᵀU` структуры | соответствует определению (tol 1e-8) |
| C9 | **Reorthogonalize улучшает `TᵀT`** при большом `A` | отклонение от `I` при `true` ≤ чем при `false` |
| C10 | **Детерминизм**: два `fit` на одних данных | побайтово идентичные выходы |
| C11 | **predict без статистик** | работает в стандартизованном пространстве, `b₀=0` |

### 7.3. Property-based тесты (`proptest`)

| # | Свойство |
|---|---|
| P1 | На случайных стандартизованных данных `TᵀT ≈ I` и `WᵀW ≈ I` при `A ≤ rank(X)` |
| P2 | `expl_var` невозрастающего приращения, в пределах [0,1] |
| P3 | `predict` не даёт NaN/Inf на конечных входах после успешного `fit` |
| P4 | Увеличение `A` (до ранга) не ухудшает MSE предсказания на обучающих данных |

### 7.4. Тесты WASM-границы

| # | Что проверяем |
|---|---|
| W1 | Плоская конкатенация колонок + `n_rows` + `y` читаются; `fit` отрабатывает |
| W2 | Выгрузки (`b`,`T`,`U`,`P`,`q`,`expl`) — `Float64Array` правильных длин (m, n·A, n·A, m·A, A, A) |
| W3 | Ошибки ядра (включая `MethodError`) прокидываются как исключение JS |
| W4 | Числовой паритет нативного ядра и WASM (tol 1e-9, f64) |
| W5 | Порядок выгрузок соответствует `WASM_OUTPUT_IDX` (prediction, regr_coeffs, T, U, X-load, Y-load) |
| W6 | `set_feature_stats` влияет только на `coefficients()`/`predict`, не на `T`/`P`/`q` |

### 7.5. Эталонная сверка с референсом (отдельной задачей CI)

Скомпилировать `pls.cpp` (Eigen) на фиксированном датасете, сохранить выходы `partialLeastSquareExtended` (T, U, P, q, b, prediction) в JSON-фикстуру и сверять Rust-результаты с допуском (с поправкой на f32 в референсе vs f64 в ядре). Это артефакт теста, не зависимость кода.

---

## 8. Гиперпараметры по умолчанию (сводка)

| Параметр | Default | Обоснование / источник |
|---|---|---|
| `components` (A) | `3` | `COMPONENTS.DEFAULT` из `pls-constants.ts`; диапазон `[1, min(m, n)]` |
| `reorthogonalize` | `false` | скорость; включать при большом `A` для строгой ортогональности `T` |

Центрирование/стандартизация — вне ядра (вариант B). μ/σ предикторов и отклика задаются опционально через `set_feature_stats`, используются только для раскрутки коэффициентов и `b₀`.

---

## 9. Этапы разработки

1. **Каркас** — два крейта, `Config::default`, ошибки, заглушки API. Тесты U5–U9.
2. **Матрица + H⁻¹ + FeatureStats** — dot/axpy/norm по колонкам, Gauss–Jordan, валидация статистик. Тесты U1–U4.
3. **PLS1 компонента 0 + инициализация** — `w,t,p,q`, защиты деления. Тест U10.
4. **Итерации с рекуррентностью w и Gram–Schmidt t** — `W,T,P,q`. Тесты C1–C2.
5. **Постобработка** — `U`, `H`, `b = W·H⁻¹·q`, накопленная дисперсия. Тесты C6, C8.
6. **Раскрутка + predict** — `b_orig`, `b₀`, предсказание на исходных данных. Тесты C5, C11.
7. **Согласование с референсом** — сверка чисел с `pls.cpp`. Тесты C3–C4, C7.
8. **Устойчивость** — reorthogonalize, детерминизм. Тесты C9–C10.
9. **Property-тесты.** P1–P4.
10. **WASM-обёртка** + `wasm-pack`, порядок выгрузок. Тесты W1–W6.
11. **CI** — `cargo test`, `wasm-pack test`, эталонная фикстура (7.5).

---

## 10. Критерии готовности (Definition of Done)

- Все тесты 7.1–7.4 зелёные нативно и в WASM.
- Ортонормированность `T`/`W` в пределах допусков (C1–C2, P1).
- Согласование с референсом `pls.cpp` в пределах f32-допуска (C7).
- Нет паник на некорректном вводе — только типизированные `Err`.
- `wasm-pack build --target web` без warning'ов; размер `.wasm` зафиксирован в README.
- Публичный API задокументирован (`cargo doc`), у каждого `pub` — пример.

---

## Приложение A. Соответствие контракту трёх спецификаций

| Аспект | linreg (Elastic Net) | PCA (NIPALS) | PLS1 (Indahl) |
|---|---|---|---|
| Тип | регрессия с регуляризацией | снижение размерности | супервизорная регрессия/латентные факторы |
| Отклик `y` | да | нет | да (один, PLS1) |
| Вход | список колонок | список колонок | список колонок + `y` |
| Масштабирование | вне ядра (B) | вне ядра (B) + центрирование на JS | вне ядра (B): центр.+станд. X и y на JS |
| Роль статистик | раскрутка коэффициентов | раскрутка нагрузок (σ) | раскрутка `b` + bias (μ/σ для X и y) |
| Дефляция | — | дефляция X (NIPALS) | **нет дефляции X** (рекуррентность w + Gram–Schmidt t) |
| Плотная ЛА | нет | нет | обращение `A×A` (H⁻¹) |
| Ошибки | `LinRegError` | `PcaError`(+`NotConverged`) | `PlsError`(+`MethodError`) |
| WASM-граница | плоская конкатенация + `n_rows` | то же | то же + `y` отдельным массивом |

## Приложение B. Замечания по референс-коду

- `pls.cpp` выполняет центрирование и нормировку **внутри** (вариант A). В нашем ядре это исключено (вариант B) — JS подаёт уже подготовленные данные; PLS1-цикл и формула `b = W·H⁻¹·q` перенесены без изменений.
- Референс возвращает `prediction = D·b` на сырых `D`; свободный член там вычисляется в TS-слое (`getRegrCoeffs`: `bias = μ_y − Σ bⱼ·μ_xⱼ`). У нас `predict` сразу учитывает `b₀`, что эквивалентно.
- Квадратичный режим (`isQuadratic` в `pls-constants.ts`) — это добавление колонок `xⱼ²` к предикторам до вызова ядра; делается на JS, ядро не меняется.
- Тип данных: референс использует `float` (f32); ядро — `f64`. При сверке (C7) учитывать разницу точности.
