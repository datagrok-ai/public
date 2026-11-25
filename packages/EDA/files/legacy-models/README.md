# Legacy models

Models of the EDA ML methods in the pre-Crux format (sci-comp-ml, XGBoost and LIBSVM wasm). They back the
`Legacy models` test category (`src/tests/legacy-models-tests.ts`). Since the methods moved to Crux (GROK-20932 stage
2) the old classes are gone; the files are applied by `cruxPredict` through `src/legacy/` (stage 2b) and
parity-tested against `reference.json`.

The files were produced once with the pre-Crux wasm and cannot be regenerated after the old wasm is removed. On a
stand they are at `System:AppData/Eda/legacy-models/`.

The first five files were produced on 2026-10-21 with EDA `1.9.0`, the other five on 2026-10-29 by the old
implementation resurrected from `f427b38627`; hyperparameters = the defaults of the old `train*` functions in
`src/package.ts`. The `big` target of the `*-bool` files is a BOOL column, `Petal.Length > 2.5` per row (50 false,
100 true).

| File | Method (`mlname`) | Dataset | Target | Features | Hyperparameters | Rows | Bytes |
|---|---|---|---|---|---|---|---|
| `linear-regression.bin` | Linear Regression | `System:DemoFiles/cars.csv` | `price` | 15 numeric columns (all but `model`, `price`) | rate 0.1, iterations 1000, L1 0, L2 0 | 30 | 1001 |
| `pls-regression.bin` | PLS Regression | `System:DemoFiles/cars.csv` | `price` | 15 numeric columns (all but `model`, `price`) | components 3 | 30 | 5248 |
| `softmax.bin` | Softmax | `System:DemoFiles/iris.csv` | `Species` | `Sepal.Length`, `Sepal.Width`, `Petal.Length`, `Petal.Width` | rate 1.0, iterations 100, penalty 0.1, tolerance 0.001 | 150 | 1986 |
| `xgboost.bin` | XGBoost | `System:DemoFiles/iris.csv` | `Species` | `Sepal.Length`, `Sepal.Width`, `Petal.Length`, `Petal.Width` | iterations 20, rate 0.3, max depth 6, lambda 1, alpha 0 | 150 | 58944 |
| `svm.bin` | SVM | `System:DemoFiles/iris.csv` | `Species` | `Sepal.Length`, `Sepal.Width`, `Petal.Length`, `Petal.Width` | RBF kernel, penalty 1, gamma 0, degree 3, offset 0, epsilon 0.1 | 150 | 7424 |
| `xgboost-regression.bin` | XGBoost | `System:DemoFiles/cars.csv` | `price` (INT) | 15 numeric columns (all but `model`, `price`) | iterations 20, rate 0.3, max depth 6, lambda 1, alpha 0 | 30 | 33152 |
| `svm-regression.bin` | SVM | `System:DemoFiles/cars.csv` | `price` (INT) | 15 numeric columns (all but `model`, `price`) | RBF kernel, penalty 1, gamma 0, degree 3, offset 0, epsilon 0.1 | 30 | 9728 |
| `xgboost-bool.bin` | XGBoost | `System:DemoFiles/iris.csv` | `big` (BOOL) | `Sepal.Length`, `Sepal.Width`, `Petal.Length`, `Petal.Width` | iterations 20, rate 0.3, max depth 6, lambda 1, alpha 0 | 150 | 18688 |
| `svm-bool.bin` | SVM | `System:DemoFiles/iris.csv` | `big` (BOOL) | `Sepal.Length`, `Sepal.Width`, `Petal.Length`, `Petal.Width` | RBF kernel, penalty 1, gamma 0, degree 3, offset 0, epsilon 0.1 | 150 | 4928 |
| `softmax-bool.bin` | Softmax | `System:DemoFiles/iris.csv` | `big` (BOOL) | `Sepal.Length`, `Sepal.Width`, `Petal.Length`, `Petal.Width` | rate 1.0, iterations 100, penalty 0.1, tolerance 0.001 | 150 | 1756 |

## `reference.json`

The predictions of the old implementation (the pre-Crux classes and wasm resurrected from `f427b38627`, run in the
browser on 2026-10-29) for every `.bin` file above, on its training features: `{created, eda, fixtures}`, where
`fixtures.<file>` = `{dataset, target, exclude, type, prediction}` (`exclude` = the columns left out of the features
besides the target and non-numeric ones; `type` = the DG type of the old prediction column: `double`, `int`, `string`
or `bool`; `prediction` = its values, nulls as `null`). The `Legacy models` category compares the current predictions
with it. 9990 bytes.
