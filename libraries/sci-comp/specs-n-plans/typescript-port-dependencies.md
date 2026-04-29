# TypeScript Port — Анализ внешних зависимостей

Анализ Python-зависимостей пяти production-файлов из [typescript-port-plan.md](typescript-port-plan.md), которые нужно портировать на TypeScript в библиотеку `@datagrok-libraries/sci-comp`.

**Решение:** в качестве математического бэкенда используется **jstat 1.9.6** (уже подключён как `dependencies` в [package.json](../package.json), используется в [stats-utils.ts](../src/time-series/feature-extraction/stats-utils.ts) и [moead.ts](../src/optimization/multi-objectives/moead/moead.ts)). jstat покрывает все спецфункции, распределения и линалгебру; недостающие тесты надстраиваются поверх.

---

## 1. Сводка зависимостей по файлам

| Файл | numpy | scipy.stats | stdlib |
|------|-------|-------------|--------|
| `statistics_fixed.py` | ✓ | `ttest_ind`, `mannwhitneyu`, `fisher_exact`, `norm.cdf`, `spearmanr`, `dunnett` | — |
| `trend_test_incidence_modified.py` | ✓ | `norm.cdf`, `norm.sf` | — |
| `williams_fixed.py` | ✓ | `t.ppf`, `t.cdf` | `dataclasses`, `typing` |
| `williams_tables.py` | — | — | `math.sqrt`, `bisect` |
| `ancova.py` | ✓ | `f.cdf`, `t.cdf` | — |

---

## 2. numpy — что используется

Тривиально портируется на `Float64Array` + утилиты.

- **Создание/форма**: `np.array`, `np.asarray`, `np.zeros`, `np.zeros_like`, `np.arange`, `np.concatenate`
- **Маски и NaN**: `np.isnan`, булевая индексация (`a[~np.isnan(a)]`)
- **Редукции**: `np.sum`, `np.mean`, `np.var(ddof=1)`, `np.all`, `np.sqrt`
- **Broadcasting**: внешнее вычитание `a[:, None] - b[None, :]` (используется в `trend_test` для JT-статистики)
- **Линейная алгебра** (только в `ancova.py`): `np.linalg.lstsq`, `np.linalg.inv`, `np.linalg.pinv`, матричное умножение `@`, транспонирование `.T` — закрывается через jstat (см. §3.3).

---

## 3. jstat как бэкенд для scipy.stats

### 3.1. Распределения — закрыты целиком

| scipy | jstat | Использование |
|-------|-------|---------------|
| `norm.cdf`, `norm.sf` | `jStat.normal.cdf(x, 0, 1)`, `1 - jStat.normal.cdf(...)` | trend tests, JT, CA |
| `t.cdf` | `jStat.studentt.cdf(t, df)` | Williams, ANCOVA |
| `t.ppf` | `jStat.studentt.inv(p, df)` | Williams (критические значения) |
| `f.cdf` | `jStat.centralF.cdf(f, df1, df2)` | ANCOVA (slope homogeneity) |
| `chi2.cdf` (резерв) | `jStat.chisquare.cdf(x, df)` | — |

### 3.2. Спецфункции — закрыты целиком

`erf`, `erfc`, `erfcinv`, `gammafn`, `gammaln`, `loggam`, `gammap`, `gammapinv`, `lowRegGamma`, `betafn`, `betaln`, `betacf`, **`ibeta`**, **`ibetainv`**, `factorial`, `factorialln`, `combination`, **`combinationln`**, `permutation`.

→ **~250 строк numerical-recipes кода писать не надо.**

### 3.3. Линейная алгебра — закрыта целиком

jstat предоставляет: `inv`, `det`, `lu`, **`cholesky`**, `gauss_jordan`, `triaUpSolve`/`triaLowSolve`, `householder`, `jacobi` (eigen).

ANCOVA-OLS реализуется через нормальные уравнения `β = (XᵀX)⁻¹ Xᵀy`:

```ts
const XtX = matMul(transpose(X), X);
const Xty = matVec(transpose(X), y);
const beta = jStat.cholesky(XtX).then(L => triaLowSolve(L, Xty)...);  // или просто jStat.inv(XtX) для p≤6
```

Для размерностей `p × p ≤ 6` (типично для ANCOVA) — обычное обращение через `jStat.inv` абсолютно надёжно.

### 3.4. Численные методы — для Dunnett

jstat предоставляет: `simpson`, `romberg`, `gauss_quadrature`, `richardson`. Используются в реализации Dunnett 1955 (см. §4).

---

## 4. Что нужно реализовать поверх jstat

| Тест | Алгоритм | Объём | Комментарий |
|------|----------|-------|-------------|
| `welchTTest` | Welch–Satterthwaite df + `jStat.studentt.cdf` | ~15 строк | `jStat.ttest` — только paired/one-sample, поэтому Welch строим сами |
| `rank()` (утилита) | Сортировка + средний ранг для тай-групп + `tieCorrection` | ~30 строк | Нужна для Mann-Whitney и Spearman |
| `mannWhitneyU` | Ранжирование объединённой выборки → U → нормальная аппроксимация с поправкой на тайбрейки → `jStat.normal.cdf` | ~30 строк | |
| `spearman` | `rank(x)`, `rank(y)` → Pearson → t-аппроксимация → `jStat.studentt.cdf` | ~10 строк | |
| `fisherExact2x2` | `jStat.hypgeom` + `jStat.combinationln`, two-sided через "minlike" | ~30 строк | |
| **`dunnett`** | Формула Dunnett (1955): одномерное интегрирование стандартного нормального ядра, использует `jStat.simpson` или `jStat.gauss_quadrature` | ~80 строк | Единственный тест без прямой jstat-поддержки |
| OLS-обёртка (для ANCOVA) | `(XᵀX)⁻¹Xᵀy` через `jStat.inv` | ~10 строк | |
| `bonferroni` | Тривиальная итерация | ~5 строк | |

**Итого собственного кода:** ~210 строк против ~780, если бы писали всё с нуля. **Снижение в 3.7×.**

---

## 5. Что меняется по сравнению с до-jstat анализом

❌ **Больше не нужно писать:**
- Lanczos-аппроксимацию `gamma`/`lgamma`.
- Регуляризованную неполную бету `I_x(a,b)` через непрерывные дроби.
- Инверсию `ibetainv` (была одной из самых сложных и тонких частей).
- `erf`, `erfc`.
- Все CDF/PPF для `normal`, `studentt`, `centralF`, `chisquare`.
- Cholesky/LU/обращение матриц.
- Численные квадратуры (нужны для Dunnett).

✅ **Что остаётся писать:**
- Тесты-агрегаторы (Welch, MWU, Spearman, Fisher, Dunnett, Bonferroni) — поверх готовых распределений и спецфункций.
- Утилиту `rank()` с tie-correction.
- Тонкую обёртку OLS для ANCOVA.

---

## 6. Подводные камни jstat

1. **Нет TypeScript-типов** — `@ts-ignore: no types` уже стоит в [stats-utils.ts:5](../src/time-series/feature-extraction/stats-utils.ts#L5). Решение: оборачиваем все обращения к jstat в типизированные функции в `stats-utils.ts` (или новом `src/distributions.ts`), наружу выставляем чистый TS API.
2. **Точность `studentt.inv` при малых df** — известные исторические шероховатости. Проверяем в тестах CSV-фикстурами против scipy.
3. **`jStat.ttest`** — не покрывает Welch (только paired/one-sample), Welch строим сами.
4. **Версия 1.9.6** зафиксирована в `package-lock.json` — sanity-check на стабильность `ibetainv` для малых α (влияет на критические значения).

---

## 7. Архитектурное решение

**jstat используется как внутренний бэкенд.** Снаружи `sci-comp` экспортирует только чистый TS API:

```
src/stats/
  distributions.ts   # типизированные обёртки: normalCdf, studentTCdf, fCdf, ibeta, ...
  rank.ts            # утилита ранжирования с tie-correction
  tests/
    welch-t.ts
    mann-whitney.ts
    spearman.ts
    fisher-exact.ts
    dunnett.ts
    jonckheere.ts
    cochran-armitage.ts
    williams.ts      # + williams-tables.ts (статические таблицы)
    ancova.ts
  multiple-comparison/
    bonferroni.ts
```

Преимущества:
- **Независимость от backend** — можно поменять jstat на собственную реализацию без слома публичного API.
- **Единообразие со существующим стилем** [stats-utils.ts](../src/time-series/feature-extraction/stats-utils.ts) (`tdistCdf`, `twoTailPvalue`).
- **Валидируемость** — тесты против CSV-фикстур, сгенерированных Python-скриптами из [validate_*.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/).
- **Никаких новых внешних зависимостей** сверх уже подключённого jstat.

---

## 8. Внутренние Python-зависимости (тривиально)

- `math.sqrt`, `bisect.bisect_right` — есть в JS из коробки.
- `dataclasses` → обычные TS-интерфейсы/классы.
- `typing.Optional` → `T | null` или `T | undefined`.
- `from williams_tables import lookup_1971, lookup_1972` — внутренний импорт между файлами.
