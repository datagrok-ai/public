# Спецификация: Softmax-классификатор (multinomial logistic regression) на Rust + WASM

**Версия:** 1.0
**Цель:** реализовать с нуля (на голом `Vec<f64>`) многоклассовую логистическую регрессию (softmax), обучаемую полнобатчевым градиентным спуском с L2-регуляризацией весов. Rust-ядро универсально; WASM — один из бэкендов. За основу взята референс-реализация на Eigen (`softmax.h`, `softmax-api.cpp`) и её TS-обвязка (`softmax-classifier.ts`).

> **Общий контракт с `linreg`, `pca-nipals`, `pls1` (вариант B).** Данные подаются как **список колонок** (`&[&[f64]]`). Центрирование и стандартизация выполняются на стороне **JS/TS**: в ядро приходят уже стандартизованные признаки. μ/σ передаются **опционально** и используются только для нормализации новых данных при `predict` в исходном пространстве. Целевые метки — целочисленные индексы классов `0..c-1`. Та же система ошибок и та же WASM-граница.

> **Отличие от референса.** В `softmax.h` центрирование/стандартизация делаются внутри (вариант A). Здесь ядро их **не** делает (вариант B); сам цикл softmax-GD, формулы градиентов и обновлений сохранены.

---

## 1. Область применения и ограничения

В объёме (in scope):
- Softmax-регрессия: `c` классов, `n` признаков, обучение GD.
- Параметры: веса `W` (c×n) + смещения `B` (c).
- L2-штраф **только на `W`** через мультипликативный распад (как в референсе); bias не регуляризуется.
- Инициализация Xavier, ранняя остановка по `|Δloss| < tol`, обработка NaN-loss.
- Взвешивание классов по частоте — воспроизведение поведения референса **под флагом** (см. 2.5).
- Предсказание (argmax) и предсказание вероятностей.
- WASM-обёртка; формат вывода `params` = `W|B` (c×(n+1)), как в `softmax.h`.

Вне объёма (out of scope) в v1.0:
- Центрирование/стандартизация внутри ядра (вариант B, делается на JS).
- Мини-батчевый/стохастический GD (референс — полнобатчевый).
- Adam/momentum и расписания learning rate.
- Мультиметка (multi-label); тут ровно один класс на образец.

---

## 2. Математическая модель

### 2.1. Обозначения

`X` — стандартизованные признаки, `n` признаков × `m` образцов. `W` — веса `c×n`, `B` — смещения длины `c`. Для образца `x` (длины `n`):

```
zₖ = Σⱼ Wₖⱼ·xⱼ + Bₖ            # логит класса k,  k = 0..c-1
```

### 2.2. Softmax-вероятности

```
P(y=k | x) = exp(zₖ) / Σ_l exp(z_l)
```

Для численной устойчивости (улучшение относительно референса — см. 2.6) вычитаем максимум: `zₖ ← zₖ − max_l z_l` перед экспонированием. Это не меняет вероятностей.

### 2.3. Функция потерь

Усреднённая кросс-энтропия (negative log-likelihood) по `m` образцам, без явного L2-члена в значении loss (как в референсе — штраф входит только в обновление весов):

```
L = −(1/m) Σᵢ log P(yᵢ = targetᵢ | xᵢ)
```

### 2.4. Градиенты и обновления (как в `softmax.h`)

`Y` — one-hot матрица `c×m`. `XT` — стандартизованные признаки `m×n`.

```
dZ = Z − Y                          # c×m, Z — матрица вероятностей
dB = mean_по_образцам(dZ)           # длина c
dW = (dZ · XT) / m                  # c×n

# обновления: L2 только на W (мультипликативный распад), bias без штрафа
W ← W·(1 − lr·penalty/m) − lr·dW
B ← B − lr·dB
```

### 2.5. Взвешивание классов (поведение референса, под флагом)

Референс домножает экспоненты на частоту класса до нормализации:

```
Zₖⱼ = exp(zₖⱼ) · classWeightₖ,   затем нормализация по столбцу
```

где `classWeightₖ` — число образцов класса `k` в обучающей выборке. Это **нестандартно** (обычный взвешенный softmax кладёт вес в функцию потерь, а не в саму вероятность) и смещает предсказания в пользу частых классов. Чтобы воспроизвести эталон точно — поведение доступно под флагом `class_weighting` (по умолчанию `true` для паритета с референсом; `false` даёт канонический softmax).

> При `class_weighting=true` для точного совпадения с референсом нормализация устойчивого softmax (2.2) применяется к взвешенным экспонентам; реализация должна сохранять тот же порядок операций, иначе числа разойдутся с фикстурой.

### 2.6. Устойчивость (отличие от референса)

Референс экспонирует `z` напрямую (`exp(Z(i,j))`), что переполняется при больших логитах. Ядро вычитает построчный максимум перед `exp` (2.2). Математически вероятности идентичны; в тестах паритета с f32-референсом (7.5) это даёт расхождение лишь на уровне арифметики с плавающей точкой, не структурное.

### 2.7. Предсказание

```
argmax_k zₖ          # эквивалентно argmax вероятностей; нормализация не нужна
```

При `class_weighting=true` predict, как в референсе, берёт argmax по `exp(zₖ)·classWeightₖ`. Predict работает на **исходных** данных: ядро нормализует вход по сохранённым μ/σ (см. 4.4).

---

## 3. Алгоритм обучения (псевдокод)

```
fit(x_columns, targets):                    # x_columns УЖЕ СТАНДАРТИЗОВАНЫ; targets ∈ 0..c-1
    n = x_columns.len();  m = x_columns[0].len()
    валидировать инварианты (4.2) → EmptyDataset / DimensionMismatch / InvalidHyperparameter
    проверить NaN/Inf, диапазон targets ∈ [0,c) → NonFiniteValue / InvalidHyperparameter
    one_hot Y (c×m); classWeight[k] = частота класса k
    W = Xavier(c, n, scale=sqrt(6/(c+n)), seed);  B = zeros(c)
    lossPrev = 0
    для iter в 0..iterCount:
        # Z = W·X + B  (логиты), затем устойчивый softmax по столбцу
        для каждого образца j:
            z = [ Σₖ Wₖ·x_j + Bₖ ]            # длина c
            z -= max(z)                         # устойчивость (2.6)
            e = exp(z) · (class_weighting ? classWeight : 1)
            Z[:,j] = e / sum(e)
        loss = −(1/m) Σⱼ log Z[target_j, j]
        если |loss − lossPrev| < tol или loss is NaN: break
        lossPrev = loss
        dZ = Z − Y
        dB = mean_столбцов(dZ)                  # длина c
        dW = (dZ · Xᵀ) / m                       # c×n
        W = W·(1 − lr·penalty/m) − lr·dW
        B = B − lr·dB
    params = [W | B]                             # c×(n+1), как в softmax.h
```

> Горячие операции — `W·X` (логиты) и `dZ·Xᵀ` (градиент весов). При колоночном хранении `X` логиты считаются проходом по столбцам-образцам, а `dW` — аккумуляцией; обе операции последовательны по памяти. Плотных обращений матриц нет.

---

## 4. Архитектура Rust-ядра

Общий стиль с предыдущими крейтами; чистый `Vec<f64>`.

### 4.1. Структура крейта

```
crates/
  softmax-core/
    src/
      lib.rs
      matrix.rs         # операции по колонкам (dot/axpy), gemm-подобные проходы W·X, dZ·Xᵀ
      stats.rs          # FeatureStats: μ/σ признаков (для predict в исходном пространстве)
      init.rs           # детерминированный Xavier (seeded PRNG, без внешних крейтов)
      softmax.rs        # устойчивый softmax, кросс-энтропия, градиенты
      model.rs          # Softmax: fit / predict / predict_proba / params
      error.rs          # типы ошибок (как в предыдущих + диапазон меток)
    tests/
      integration.rs
  softmax-wasm/
    src/lib.rs
```

### 4.2. Представление данных — список колонок

`x_columns: &[&[f64]]`, `x_columns[j]` — признак `j` длиной `m` (число образцов); `n = x_columns.len()`. Метки `targets: &[i32]` (или `&[u32]`) длины `m`, значения `0..c-1`. Инварианты (ср. `softmax-api.cpp`): равные длины колонок, `targets.len()==m`, `1 ≤ c`, `n ≥ 1`, метки в диапазоне, отсутствие NaN/Inf. Нарушение → соответствующая типизированная ошибка.

### 4.3. Ключевые типы (эскиз API)

```rust
// error.rs
#[derive(Debug, Clone, PartialEq)]
pub enum SoftmaxError {
    DimensionMismatch { expected: usize, got: usize },  // ~ INCORRECT_SIZES
    EmptyDataset,
    NotFitted,
    InvalidHyperparameter(String),  // ~ ICORRECT_HYPER_PARAMS; в т.ч. метка вне [0,c)
    NonFiniteValue,
}

pub struct Config {
    pub classes: usize,        // c; обязателен (≥1)
    pub learning_rate: f64,    // default 1.0  (DEFAULT_LEARNING_RATE)
    pub iterations: usize,     // default 100  (DEFAULT_ITER_COUNT)
    pub penalty: f64,          // L2 на W; default 0.1 (DEFAULT_PENALTY)
    pub tolerance: f64,        // ранняя остановка; default 0.001 (DEFAULT_TOLERANCE)
    pub class_weighting: bool, // воспроизведение референса; default true
    pub seed: u64,             // для Xavier-инициализации; default 42
}

impl Default for Config { /* значения выше */ }

// stats.rs — вариант B: μ/σ только для нормализации входа в predict.
pub struct FeatureStats { means: Vec<f64>, stds: Vec<f64> } // длина n
impl FeatureStats {
    pub fn new(means: Vec<f64>, stds: Vec<f64>) -> Result<Self, SoftmaxError>;
}

pub struct Softmax {
    weights: Vec<f64>,      // W, c×n (row-major по классам)
    biases: Vec<f64>,       // B, длина c
    class_weights: Vec<f64>,// частоты классов с обучения (для predict при class_weighting)
    classes: usize, n_features: usize,
    loss_history: Vec<f64>,
    stats: Option<FeatureStats>,
    config: Config,
}

impl Softmax {
    pub fn new(config: Config) -> Result<Self, SoftmaxError>;

    // Опционально: μ/σ признаков — для нормализации входа в predict (вариант B).
    pub fn set_feature_stats(&mut self, means: Vec<f64>, stds: Vec<f64>)
        -> Result<(), SoftmaxError>;

    // КОНТРАКТ: x_columns УЖЕ СТАНДАРТИЗОВАНЫ; targets ∈ 0..c-1.
    pub fn fit(&mut self, x_columns: &[&[f64]], targets: &[i32]) -> Result<(), SoftmaxError>;

    // Предсказание на ИСХОДНЫХ данных: вход нормализуется по μ/σ (если заданы),
    // иначе ожидаются уже стандартизованные признаки. Возврат — индексы классов.
    pub fn predict(&self, x_columns: &[&[f64]]) -> Result<Vec<usize>, SoftmaxError>;
    // Матрица вероятностей m×c (row-major по образцам).
    pub fn predict_proba(&self, x_columns: &[&[f64]]) -> Result<Vec<f64>, SoftmaxError>;

    pub fn loss_history(&self) -> &[f64];
    // params в формате референса: c×(n+1), последний столбец — bias.
    pub fn params(&self) -> Vec<f64>;
    pub fn weights(&self) -> &[f64];
    pub fn biases(&self) -> &[f64];
}
```

> Как и прежде: `fit` работает над поданными (стандартизованными) данными; `set_feature_stats` влияет только на `predict`/`predict_proba` (нормализация входа). Согласованность масштаба — ответственность вызывающей стороны.

### 4.4. Нормализация входа в predict

Если μ/σ заданы, `predict`/`predict_proba` нормализуют каждый признак `xⱼ ← (xⱼ − μⱼ)/σⱼ` (при `σⱼ ≤ 0` колонка → 0, как в `normalized()` референса), затем считают логиты. Без μ/σ вход считается уже стандартизованным.

### 4.5. Детерминированная Xavier-инициализация

Референс использует `Matrix::Random` (недетерминирован). Для воспроизводимости тестов ядро применяет seeded PRNG (например, splitmix64/xorshift на `u64`, без внешних крейтов), генерируя `W ∈ [−scale, scale]`, `scale = sqrt(6/(c+n))`. Это означает, что точное совпадение весов с референсом невозможно (у него случайный старт); сверка (7.5) идёт по сошедшимся метрикам/предсказаниям, а не по сырым `W`.

---

## 5. WASM-бэкенд

Граница как в предыдущих крейтах: признаки — плоская конкатенация стандартизованных колонок (`Float64Array`) + `n_rows`(=m); метки — `Int32Array`. Выход `params` — плоский `Float64Array` `c×(n+1)`, как `P` в `softmax.h`. Порядок/формат согласован с `_fitSoftmax` из `softmax-classifier.ts`.

```rust
// softmax-wasm/src/lib.rs
use wasm_bindgen::prelude::*;
use softmax_core::{Softmax, Config};

#[wasm_bindgen]
pub struct WasmSoftmax { inner: Softmax }

#[wasm_bindgen]
impl WasmSoftmax {
    #[wasm_bindgen(constructor)]
    pub fn new(classes: usize, learning_rate: f64, iterations: usize,
               penalty: f64, tolerance: f64) -> Result<WasmSoftmax, JsError> {
        let cfg = Config { classes, learning_rate, iterations, penalty, tolerance,
                           ..Default::default() };
        Ok(Self { inner: Softmax::new(cfg).map_err(to_js)? })
    }

    pub fn set_feature_stats(&mut self, means: &[f64], stds: &[f64]) -> Result<(), JsError> {
        self.inner.set_feature_stats(means.to_vec(), stds.to_vec()).map_err(to_js)
    }

    // flat_x: конкатенация стандартизованных колонок; n_rows=m; targets: Int32Array.
    pub fn fit(&mut self, flat_x: &[f64], n_rows: usize, targets: &[i32]) -> Result<(), JsError> {
        let cols = to_columns(flat_x, n_rows).map_err(to_js)?;
        self.inner.fit(&cols, targets).map_err(to_js)
    }

    pub fn predict(&self, flat_x: &[f64], n_rows: usize) -> Result<Vec<u32>, JsError> {
        let cols = to_columns(flat_x, n_rows).map_err(to_js)?;
        Ok(self.inner.predict(&cols).map_err(to_js)?.iter().map(|&v| v as u32).collect())
    }

    pub fn params(&self) -> Vec<f64> { self.inner.params() }      // c×(n+1) плоско
    pub fn loss_history(&self) -> Vec<f64> { self.inner.loss_history().to_vec() }
}

fn to_columns(flat: &[f64], n_rows: usize) -> Result<Vec<&[f64]>, SoftmaxError> {
    if n_rows == 0 || flat.len() % n_rows != 0 { /* -> DimensionMismatch */ }
    Ok(flat.chunks_exact(n_rows).collect())
}
```

Сборка: `wasm-pack build crates/softmax-wasm --target web`.

---

## 6. Валидационные тесты

Стиль как в предыдущих документах: (а) unit, (б) числовые, (в) property, (г) WASM.

### 6.1. Unit-тесты компонентов

| # | Что проверяем | Ожидание |
|---|---|---|
| U1 | `matrix`: dot/axpy и проходы `W·X`, `dZ·Xᵀ` | совпадение с ручным расчётом |
| U2 | Устойчивый softmax: вычитание max | сумма вероятностей = 1; нет переполнения при больших логитах |
| U3 | Softmax инвариантен к сдвигу `z += const` | вероятности неизменны (tol 1e-12) |
| U4 | One-hot и частоты классов | корректные индексы и счётчики |
| U5 | Валидация: `lr≤0`, `iter<1`, `penalty<0`, `tol≤0`, `classes<1` | `Err(InvalidHyperparameter)` (ср. `ICORRECT_HYPER_PARAMS`) |
| U6 | Размеры: `targets.len()≠m`, колонки разной длины | `Err(DimensionMismatch)` (ср. `INCORRECT_SIZES`) |
| U7 | Метка вне `[0,c)` | `Err(InvalidHyperparameter)` |
| U8 | NaN/Inf на входе | `Err(NonFiniteValue)` |
| U9 | Геттеры/`predict` до `fit` | `Err(NotFitted)` |
| U10 | Xavier детерминирован по seed | два прогона — идентичные стартовые `W` |

### 6.2. Корректность обучения (числовые тесты)

| # | Сценарий | Ожидание |
|---|---|---|
| C1 | **Линейно разделимые 2 класса**, `class_weighting=false` | точность на обучении = 1.0; loss монотонно убывает |
| C2 | **3 класса, разделимые** | argmax-точность = 1.0 после сходимости |
| C3 | **Монотонность loss** при адекватном `lr` | `loss_history` не возрастает (до выхода на плато) |
| C4 | **Ранняя остановка по tol** | число итераций в истории < `iterations` |
| C5 | **NaN-loss прерывает обучение** (искусственно большой `lr`) | цикл выходит, без паники |
| C6 | **L2 на W сжимает веса**: `penalty` мал vs велик | `‖W‖` при большом penalty строго меньше; bias не затронут |
| C7 | **Bias не регуляризуется** | при сильном penalty `B` не загоняется к нулю мультипликативно |
| C8 | **predict = argmax логитов** | совпадает с argmax `predict_proba` |
| C9 | **Нормализация входа в predict** (заданы μ/σ): сырые данные | результат = predict на предварительно стандартизованных данных (tol 1e-9) |
| C10 | **class_weighting воспроизводит референс**: вероятности с домножением на частоты | соответствует определению 2.5 (tol 1e-8) |
| C11 | **Детерминизм** при фиксированном seed | побайтово идентичные `params` |
| C12 | **Эквивалентность канонического softmax** (`class_weighting=false`) аналитике для 2 классов = логистической регрессии | коэффициенты согласуются (tol 1e-4) |

### 6.3. Property-based тесты (`proptest`)

| # | Свойство |
|---|---|
| P1 | `predict_proba` строки суммируются к 1 (tol 1e-9), все ∈ [0,1] |
| P2 | softmax инвариантен к константному сдвигу логитов |
| P3 | `predict`/`predict_proba` без NaN/Inf на конечных входах после `fit` |
| P4 | Увеличение `iterations` (до плато) не увеличивает обучающий loss |

### 6.4. Тесты WASM-границы

| # | Что проверяем |
|---|---|
| W1 | Плоская конкатенация колонок + `n_rows` + `Int32Array` меток читаются; `fit` отрабатывает |
| W2 | `params` — `Float64Array` длины `c·(n+1)`; формат `[W|B]` как в `softmax.h` |
| W3 | Ошибки ядра прокидываются как исключение JS (включая `INCORRECT_SIZES`/`ICORRECT_HYPER_PARAMS`-аналоги) |
| W4 | Числовой паритет нативного ядра и WASM (tol 1e-9, f64) при фиксированном seed |
| W5 | `set_feature_stats` влияет только на `predict`, не на обучение |

### 6.5. Эталонная сверка с референсом (отдельной задачей CI)

Так как Xavier-старт референса случаен, сверка идёт **по сошедшимся характеристикам**, не по сырым `W`: на фиксированном разделимом датасете сравнить итоговые предсказания/точность и (при заданном одинаковом старте, если референс пропатчить на фиксированный seed) — значения `params` с допуском (f32 vs f64). Фикстура — JSON-выгрузка из скомпилированного `softmax.h`.

---

## 7. Гиперпараметры по умолчанию (сводка)

| Параметр | Default | Источник |
|---|---|---|
| `learning_rate` | `1.0` | `DEFAULT_LEARNING_RATE` |
| `iterations` | `100` | `DEFAULT_ITER_COUNT` |
| `penalty` | `0.1` | `DEFAULT_PENALTY` |
| `tolerance` | `0.001` | `DEFAULT_TOLERANCE` |
| `class_weighting` | `true` | паритет с референсом (см. 2.5); `false` = канонический softmax |
| `seed` | `42` | детерминированный Xavier |

Центрирование/стандартизация — вне ядра (вариант B). μ/σ задаются опционально через `set_feature_stats`, используются только для нормализации входа в `predict`.

---

## 8. Этапы разработки

1. **Каркас** — два крейта, `Config::default`, ошибки, заглушки. Тесты U5–U9.
2. **Матрица + Xavier + FeatureStats** — проходы `W·X`/`dZ·Xᵀ`, seeded PRNG, валидация статистик. Тесты U1, U10.
3. **Устойчивый softmax + loss** — вычитание max, кросс-энтропия, one-hot, частоты. Тесты U2–U4.
4. **Цикл GD** — градиенты `dW/dB`, обновления с L2 на W, ранняя остановка/NaN. Тесты C1–C5.
5. **Регуляризация и bias** — проверка распада весов и неприкосновенности bias. Тесты C6–C7.
6. **Predict + нормализация входа** — argmax, `predict_proba`, μ/σ. Тесты C8–C9, C12.
7. **class_weighting** — воспроизведение референса под флагом. Тест C10.
8. **Детерминизм + property.** C11, P1–P4.
9. **WASM-обёртка** + `wasm-pack`, формат `params`. Тесты W1–W5.
10. **CI** — `cargo test`, `wasm-pack test`, эталонная фикстура (6.5).

---

## 9. Критерии готовности (Definition of Done)

- Все тесты 6.1–6.4 зелёные нативно и в WASM.
- Устойчивый softmax не переполняется на больших логитах (U2); вероятности корректны (P1).
- Поведение `class_weighting=true` воспроизводит логику референса (C10).
- Нет паник на некорректном вводе — только типизированные `Err`.
- `wasm-pack build --target web` без warning'ов; размер `.wasm` зафиксирован в README.
- Публичный API задокументирован (`cargo doc`), у каждого `pub` — пример.

---

## Приложение A. Соответствие контракту четырёх спецификаций

| Аспект | linreg | PCA (NIPALS) | PLS1 | Softmax |
|---|---|---|---|---|
| Задача | регрессия | снижение размерности | супервиз. регрессия | многоклассовая классификация |
| Отклик | непрерывный | нет | непрерывный (1) | категориальный (`c` классов) |
| Вход | список колонок | список колонок | колонки + `y` | колонки + `targets:i32` |
| Масштабирование | вне ядра (B) | вне ядра (B) | вне ядра (B) | вне ядра (B) |
| Роль статистик | раскрутка коэф. | раскрутка нагрузок | раскрутка `b`+bias | нормализация входа в predict |
| Оптимизация | GD + prox (L1) | NIPALS | Indahl (без дефляции) | полнобатчевый GD |
| Регуляризация | L1+L2 (Elastic Net) | — | — | L2 на `W` (распад), bias свободен |
| Ошибки | `LinRegError` | `PcaError` | `PlsError` | `SoftmaxError` |
| WASM-граница | плоские колонки + `n_rows` | то же | + `y` | + `targets:Int32Array` |

## Приложение B. Замечания по референс-коду

- **Центрирование/стандартизация** в `softmax.h` (`XT.rowwise() - avgX`, деление на `stdDevX`) исключены из ядра — вариант B, делается на JS. Цикл GD и формулы градиентов перенесены без изменений.
- **L2 только на W**: `W = W·(1 − lr·penalty/m) − lr·dW`; `B` обновляется без штрафа — сохранено точно.
- **Взвешивание классов** (`classWeights` в экспоненте) — нестандартный приём референса; воспроизводится под флагом `class_weighting` (default `true`). Канонический softmax — `false`.
- **Устойчивость**: референс делает `exp(Z)` напрямую (риск overflow); ядро вычитает max — вероятности те же, но без переполнения. Это сознательное улучшение, отмеченное для сверки.
- **Xavier-инициализация** в референсе случайна (`Matrix::Random`); ядро использует seeded PRNG для детерминизма — точное совпадение сырых `W` с референсом невозможно, сверка по сошедшимся метрикам.
- **Формат вывода** `params` = `[W | B]`, `c×(n+1)` — последний столбец bias, как `P` в `softmax.h` и как ожидает `softmax-classifier.ts`.
- **predict** референса берёт argmax по `exp(score)·classWeight` без нормализации — эквивалентно argmax логитов (при `class_weighting=false`); сохранено.
