# Baselines

`baseline.json` holds the numbers of the pre-Crux EDA implementation (sci-comp-ml, XGBoost 3.3.0 and LIBSVM 3.37
wasm) that the Crux methods are compared against: the thresholds of PLAN-ru section 4 (GROK-20932), asserted by the
`Crux thresholds` test category (`src/tests/crux-thresholds-tests.ts`). Datasets, folds and the file layout are
defined in `src/tests/baseline.ts`. On a stand the file is at `System:AppData/Eda/baselines/baseline.json`.

The file was produced once with the pre-Crux wasm and cannot be regenerated after the old wasm is removed. `fitMs`
is this machine's: Intel Core Ultra 7 255HX, 32 GB RAM, Windows 11, Chrome (headless, `grok test`) against the local
Docker stand. Produced on 2026-10-22 with EDA `1.9.0`. The 2x fit-time limit is asserted only for cases with
`fitMs` >= 100 ms; faster ones are dominated by jitter and only logged.

## Layout

- `cases["<method>/<case>"]`: `{score, fitMs, scoreName}` - the mean test-fold score (`accuracy` for
  classification, `r2` for regression) and the mean fit time in ms over 5 folds. Fold of row `i` is `i % 5`; each
  fold is tested by a model trained on the other four. `linreg/cars-price` also has `trainScore`: the mean `r2` of
  each fold's model on its own training rows (24 rows, 15 features: out-of-fold `r2` is not comparable there, decided
  2026-10-22); captured on 2026-10-22 with the same old wasm, before its removal.
- `pca.FF` / `TF` / `FT` / `TT`: `computePCA` on iris (`Sepal.Length`, `Sepal.Width`, `Petal.Length`,
  `Petal.Width`), 4 components, `center` / `scale` flags (first letter `center`, second `scale`); 4 x 150 scores.
- `pls`: `getPlsAnalysis` on cars, target `price`, 15 numeric features, 3 components: `tScores`, `uScores` (3 x 30),
  `xLoadings` (3 x 15), `yLoadings` (3), `coef` (15), `vip` (15).

## Cases

| Key | Dataset | Target | Features |
|---|---|---|---|
| `xgboost/cars-price`, `linreg/cars-price`, `pls/cars-price` | `System:DemoFiles/cars.csv` | `price` | 15 numeric columns (all but `model`, `price`) |
| `xgboost/iris-species`, `svm/iris-species`, `softmax/iris-species` | `System:DemoFiles/iris.csv` | `Species` | 4 numeric columns (all but `col 1`) |
| `svm/iris-petal` | `System:DemoFiles/iris.csv` | `Petal.Length` | 3 numeric columns (all but `col 1`, `Petal.Length`) |
| `xgboost/wine-type`, `svm/wine-type`, `softmax/wine-type` | `System:DemoFiles/winequality.csv` | `type` | 11 numeric columns (all but `quality`) |
| `xgboost/wine-quality`, `linreg/wine-quality` | `System:DemoFiles/winequality.csv` | `quality` | the same 11 numeric columns |
| `softmax/demog-sex` | `System:DemoFiles/demog.csv`, rows with `AGE`, `HEIGHT`, `WEIGHT`, `SEX` set | `SEX` | `AGE`, `HEIGHT`, `WEIGHT` |
| `linreg/demog-weight` | `System:DemoFiles/demog.csv`, rows with `AGE`, `HEIGHT`, `WEIGHT` set | `WEIGHT` | `AGE`, `HEIGHT` |
| `xgboost/synth-cls-100k` | `classificationDataset(100000, 10, true, 42)` | `Labels` (4 classes) | `#0`..`#9` |
| `xgboost/synth-reg-100k` | 100000 rows, `x0`..`x9` from `mulberry32(42)`, `y = sum (j + 1) * xj` | `y` | `x0`..`x9` |
| `svm/synth-cls-10k` | `classificationDataset(10000, 10, true, 42)` | `Labels` (4 classes) | `#0`..`#9` |
| `linreg/ill-linreg` | 200 rows, Float64: `t = i + 1`, `t_eps = t + 1e-6 * noise`, `z`; `noise`, `z`, `noise2` from `mulberry32(7)` in [-1, 1]; `y = 1 + 2 t - z + 0.05 * noise2` | `y` | `t`, `t_eps`, `z` |
| `linreg/onehot-linreg` | `System:DemoFiles/demog.csv`, rows with `AGE`, `HEIGHT`, `WEIGHT`, `SEX`, `RACE` set | `AGE` | `HEIGHT`, `WEIGHT`, one 0/1 INT column per category of `SEX` and `RACE` (`SEX=M`, ...) |

Hyperparameters are the defaults of the pre-Crux `train*` functions in `src/package.ts`; the Crux runs (`FITS` in
`src/tests/baseline.ts`) use the same values without `rate` and the Linear Regression `tolerance`:

| Method | Hyperparameters |
|---|---|
| `xgboost` | iterations 20, rate 0.3, max depth 6, lambda 1, alpha 0 |
| `svm` | RBF kernel, penalty 1, gamma 0 (= 1 / features), degree 3, offset 0, epsilon 0.1 |
| `softmax` | rate 1.0, iterations 100, penalty 0.1, tolerance 0.001 |
| `linreg` | rate 0.1, iterations 1000, L1 0, L2 0, tolerance 1e-7 |
| `pls` | components 3 |

Regression predictions of the old XGBoost and SVM on an integer target (`price`, `quality`) are rounded to integers,
so their `r2` is computed on rounded values.
