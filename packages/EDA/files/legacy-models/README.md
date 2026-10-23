# Legacy models

One model per EDA ML method in the pre-Crux format (sci-comp-ml, XGBoost and LIBSVM wasm). They back the
`Legacy models` test category (`src/tests/legacy-models-tests.ts`). Since the methods moved to Crux (GROK-20932 stage
2) the old classes are gone, and `cruxPredict` rejects every file with the "retrain the model" error. The files are
also kept for a possible later legacy support.

The files were produced once with the pre-Crux wasm and cannot be regenerated after the old wasm is removed. On a
stand they are at `System:AppData/Eda/legacy-models/`.

Produced on 2026-10-21 with EDA `1.9.0`, hyperparameters = the defaults of the `train*` functions in `src/package.ts`.

| File | Method (`mlname`) | Dataset | Target | Features | Hyperparameters | Rows | Bytes |
|---|---|---|---|---|---|---|---|
| `linear-regression.bin` | Linear Regression | `System:DemoFiles/cars.csv` | `price` | 15 numeric columns (all but `model`, `price`) | rate 0.1, iterations 1000, L1 0, L2 0 | 30 | 1001 |
| `pls-regression.bin` | PLS Regression | `System:DemoFiles/cars.csv` | `price` | 15 numeric columns (all but `model`, `price`) | components 3 | 30 | 5248 |
| `softmax.bin` | Softmax | `System:DemoFiles/iris.csv` | `Species` | `Sepal.Length`, `Sepal.Width`, `Petal.Length`, `Petal.Width` | rate 1.0, iterations 100, penalty 0.1, tolerance 0.001 | 150 | 1986 |
| `xgboost.bin` | XGBoost | `System:DemoFiles/iris.csv` | `Species` | `Sepal.Length`, `Sepal.Width`, `Petal.Length`, `Petal.Width` | iterations 20, rate 0.3, max depth 6, lambda 1, alpha 0 | 150 | 58944 |
| `svm.bin` | SVM | `System:DemoFiles/iris.csv` | `Species` | `Sepal.Length`, `Sepal.Width`, `Petal.Length`, `Petal.Width` | RBF kernel, penalty 1, gamma 0, degree 3, offset 0, epsilon 0.1 | 150 | 7424 |
