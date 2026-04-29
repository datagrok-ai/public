# TypeScript Port Plan — SENDEX Statistical Methods

План портирования статистических методов SENDEX с Python на TypeScript. Документ включает исходный отчёт о валидации и полную инвентаризацию production- и validation-файлов, которые нужно перенести.

---

## 1. Validation Report (исходный)

> Скопировано из [c:/Users/vmaka/Datagrok/SEND/send-data-browser/docs/validate/stats/report.md](c:/Users/vmaka/Datagrok/SEND/send-data-browser/docs/validate/stats/report.md)

# Statistical Methods Validation Report

This report documents the validation of statistical methods used in SENDEX.

Sources:

* [statistics.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics.py)

* [williams.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/williams.py)

* [ancova.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/ancova.py)

Validation:

* [validate_all.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_all.py)

* [validate_fisher_boschloo.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_fisher_boschloo.py)

* [validate_ancova.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_ancova.py)

* [validate_dunnett.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_dunnett.py)

* [validate_fixed_williams.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_fixed_williams.py)

* [validate_hedges_g.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_hedges_g.py)

* [validate_helpers.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_helpers.py)

* [validate_trend_test.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_trend_test.py)

* [validate_trend_test_incidence.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_trend_test_incidence.py)

* [validate_trend_test_incidence_modified.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_trend_test_incidence_modified.py)

Updated & new files files:

* [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py)

* [williams_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/williams_fixed.py)

* [williams_tables.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/williams_tables.py)

* [trend_test_incidence_modified.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/trend_test_incidence_modified.py)

## Findings

### ⚠️ Dunnett's test

Original method doesn't track statistics at all.

| Aspect | Original (`statistics.py`) | Fixed (`statistics_fixed.py`) |
|--------|---------------------------|-------------------------------|
| `dunnett_pairwise` statistic | Always `None` | Actual value from `result.statistic[j]` |
| Fallback (Welch) statistic | Not captured | Captured from `welch_t_test()` |

Validated against `scipy.stats.dunnett` (direct call) using the Dunnett (1955) blood count dataset. Both p-values and test statistics are compared within tolerance. See [validate_dunnett.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_dunnett.py).

### ✅ Jonckheere-Terpstra trend test

Checked against two independent sources:

* manual JT via `scipy.stats.mannwhitneyu`
* `jonckheere-test` library

### ✅ Hedges' g

Checked:

* Cohen's d vs NIST.

* J correction: approximation vs exact.

* End-to-end - verifies that the full pipeline works correctly.

### ❌ Williams

[This audit](c:/Users/vmaka/Datagrok/SEND/send-data-browser/docs/validate/stats/williams/williams_audit_report.md) found:

✅ Correct:

* PAVA (isotonic regression) algorithm
* t-bar statistic formula
* Step-down logic

❌ Two critical errors:

* All 88 critical values in WILLIAMS_TABLE are wrong (see also [here](c:/Users/vmaka/Datagrok/SEND/send-data-browser/docs/validate/stats/williams/williams_table_audit.md)).

* Monte Carlo fallback gives incorrect values during step-down.

**Practical impact**: the test could be too liberal at lower doses.

**References** (original papers):

* [paper1](https://drive.google.com/file/d/1JG96l135hAyZILspqOkE4K5wJpCUBSzV/view?usp=drive_link)

* [paper2](https://drive.google.com/file/d/1YkgadKQF48wYeAyaTPFE_mt2Hvn7uwF8/view?usp=drive_link)

**Solution**:

* Extracted the tables from the papers: [williams_tables.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/williams_tables.py) and [williams_critical_values.xlsx](c:/Users/vmaka/Datagrok/SEND/send-data-browser/docs/validate/stats/williams/williams_critical_values.xlsx)

* Implemented the fixed version - [williams_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/williams_fixed.py)

* Validated against orignial papers - [validate_fixed_williams.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_fixed_williams.py) ✅

### ✅ ANCOVA

The `ancova.py` implementation is **mathematically correct**: all 17 components (OLS, adjusted means, SE, pairwise comparisons, slope homogeneity, effect decomposition, Hedges' g) match published formulas and reproduce SAS output to within 10⁻⁵.

More [details](c:/Users/vmaka/Datagrok/SEND/send-data-browser/docs/validate/stats/ancova/ancova_verification_summary.md)

### ✅ Fisher's Exact Test & Boschloo's Exact Test

Both `fisher_exact` and `boschloo_exact` from SciPy are safe to use in the SENDEX pipeline for 2×2 incidence table analysis.

More [details](c:/Users/vmaka/Datagrok/SEND/send-data-browser/docs/validate/stats/fisher-boschloo/fisher-boschloo-report.md)

### ✅ Cochran-Armitage trend test for incidence

Implemented the modified version [trend_test_incidence_modified.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/trend_test_incidence_modified.py). The modified function extends the original with additional parameters while preserving backward compatibility.

Two versions of the implementation are verified: the original (baseline) version with fixed interface, and the modified version with extended functionality. Both are validated

More [details](c:/Users/vmaka/Datagrok/SEND/send-data-browser/docs/validate/stats/Cochran-Armitage-trend-test-for-incidence/verification_report.md).

### ✅ Rest

`welch_t_test`, `mann_whitney_u`, `spearman_correlation`,`severity_trend`,`welch_pairwise`, and `bonferroni_correct`are validated against published numerical examples from peer-reviewed publications and reference software (R, SPSS, StatsDirect), complemented by cross-validation against direct `scipy.stats` calls and hand-calculated values.

More [details](c:/Users/vmaka/Datagrok/SEND/send-data-browser/docs/validate/stats/rest/statistics_verification_report.md).

---

## 2. Production Source Files (для портирования)

Канонические файлы, которые нужно перенести на TypeScript. Старые `statistics.py` и `williams.py` **портировать не нужно** — они содержат баги, описанные выше.

| Файл | Строк | Содержимое |
|------|-------|-----------|
| [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) | 285 | `welch_t_test`, `mann_whitney_u`, `fisher_exact_2x2`, `trend_test` (Jonckheere-Terpstra), `trend_test_incidence` (базовый Cochran-Armitage), `compute_effect_size` (Hedges' g), `spearman_correlation`, `severity_trend`, `dunnett_pairwise`, `welch_pairwise`, `bonferroni_correct` |
| [trend_test_incidence_modified.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/trend_test_incidence_modified.py) | 281 | `trend_test_incidence` (расширенная сигнатура), `threshold_test`, внутренние `_modified_test`, `_p_from_z`, `_degenerate_result` |
| [williams_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/williams_fixed.py) | 404 | `pava_increasing`, `pava_decreasing`, `_get_critical_value`, классы `WilliamsResult` / `WilliamsTestOutput`, `williams_test`, `williams_from_dose_groups`, `williams_from_group_stats` |
| [williams_tables.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/williams_tables.py) | 357 | Извлечённые из оригинальных статей таблицы критических значений (используется `williams_fixed.py`) |
| [ancova.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/ancova.py) | 283 | `_fit_ols`, `_f_compare`, `run_ancova`, `ancova_from_dose_groups` |

**Итого:** 5 файлов, ~1610 строк.

---

## 3. Validation Files (тесты для портирования)

| Файл | Тестируемые функции | Источник эталонов |
|------|---------------------|-------------------|
| [validate_helpers.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_helpers.py) | — | Общие утилиты (`check`, `report`) для остальных тестов |
| [validate_all.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_all.py) | — | Запускает все валидации последовательно |
| [validate_dunnett.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_dunnett.py) | `dunnett_pairwise` | `scipy.stats.dunnett` + датасет Dunnett (1955) |
| [validate_fisher_boschloo.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_fisher_boschloo.py) | `fisher_exact_2x2`, Boschloo | scipy + опубликованные примеры |
| [validate_hedges_g.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_hedges_g.py) | `compute_effect_size` | NIST/SEMATECH Dataplot + Fisher's Iris dataset |
| [validate_trend_test.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_trend_test.py) | `trend_test` (Jonckheere-Terpstra) | manual JT через `mannwhitneyu` + `jonckheere-test` library |
| [validate_trend_test_incidence.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_trend_test_incidence.py) | `trend_test_incidence` (базовый) | scipy + опубликованные примеры |
| [validate_trend_test_incidence_modified.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_trend_test_incidence_modified.py) | модифицированный `trend_test_incidence` + `threshold_test` | scipy + опубликованные примеры |
| [validate_ancova.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_ancova.py) | `run_ancova`, `ancova_from_dose_groups` | SAS PROC GLM (точность 10⁻⁵) |
| [validate_fixed_williams.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_fixed_williams.py) | `williams_test`, `williams_from_dose_groups`, `williams_from_group_stats` | Оригинальные статьи Williams (1971, 1972) |
| [validate_rest_statistics.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_rest_statistics.py) | `welch_t_test`, `mann_whitney_u`, `spearman_correlation`, `severity_trend`, `welch_pairwise`, `bonferroni_correct` | R, SPSS, StatsDirect + scipy + ручные расчёты |

**Итого:** 11 файлов.

---

## 4. Сводное соответствие «метод → реализация → тест»

| # | Метод | Source | Validation |
|---|-------|--------|-----------|
| 1 | Welch's t-test | [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) | [validate_rest_statistics.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_rest_statistics.py) |
| 2 | Mann-Whitney U | [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) | [validate_rest_statistics.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_rest_statistics.py) |
| 3 | Spearman correlation | [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) | [validate_rest_statistics.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_rest_statistics.py) |
| 4 | Severity trend | [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) | [validate_rest_statistics.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_rest_statistics.py) |
| 5 | Welch pairwise | [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) | [validate_rest_statistics.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_rest_statistics.py) |
| 6 | Bonferroni correction | [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) | [validate_rest_statistics.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_rest_statistics.py) |
| 7 | Fisher exact 2×2 | [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) | [validate_fisher_boschloo.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_fisher_boschloo.py) |
| 8 | Boschloo exact | (scipy wrapper) | [validate_fisher_boschloo.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_fisher_boschloo.py) |
| 9 | Jonckheere-Terpstra | [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) (`trend_test`) | [validate_trend_test.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_trend_test.py) |
| 10 | Cochran-Armitage (base) | [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) (`trend_test_incidence`) | [validate_trend_test_incidence.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_trend_test_incidence.py) |
| 11 | Cochran-Armitage (modified) + threshold | [trend_test_incidence_modified.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/trend_test_incidence_modified.py) | [validate_trend_test_incidence_modified.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_trend_test_incidence_modified.py) |
| 12 | Hedges' g (Cohen's d + J) | [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) (`compute_effect_size`) | [validate_hedges_g.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_hedges_g.py) |
| 13 | Dunnett's test | [statistics_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/statistics_fixed.py) (`dunnett_pairwise`) | [validate_dunnett.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_dunnett.py) |
| 14 | Williams test | [williams_fixed.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/williams_fixed.py) + [williams_tables.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/williams_tables.py) | [validate_fixed_williams.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_fixed_williams.py) |
| 15 | ANCOVA | [ancova.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/ancova.py) | [validate_ancova.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/validate_ancova.py) |

**Итого:** 15 методов, 5 source-файлов, 11 validation-файлов.
