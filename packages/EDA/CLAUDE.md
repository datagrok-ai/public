# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

EDA (Exploratory Data Analysis) is a Datagrok package providing statistical analysis and machine learning tools. It includes dimensionality reduction (PCA, UMAP, t-SNE, SPE), supervised learning (linear regression, softmax, PLS, XGBoost, SVM, random forest, extra trees, hist gradient boosting, CatBoost, decision tree, k-nearest neighbors, naive Bayes, LDA, neural network, RealMLP, TabM - all via Crux), group comparison (two-sample t-test, ANOVA, control comparisons), missing data imputation, Pareto optimization, and probabilistic multi-parameter optimization (pMPO).

## Build Commands

```bash
grok setup                 # Once per checkout, at the repository root (pnpm workspace, see packages/BUILD.MD)
grok build                 # This package and the libraries it depends on, in order, cached
grok build --typecheck     # Build and type-check
pnpm run lint              # ESLint check
grok publish <server>      # Build and deploy a debug version

npm test                   # Run tests against localhost
grok test                  # Run all tests
grok test --test "TestName"      # Run specific test
grok test --category "CategoryName"  # Run tests in category
grok test --gui            # Run with visible browser
```

## Architecture

### High-Level Structure

The package combines TypeScript, the Crux WebAssembly module, and web workers for performance-critical computations:

- **Crux** (`src/crux/`, module in `crux/eda/`): PCA, PLS and every ML method run on Crux (scikit-learn, XGBoost and CatBoost ports, pytabkit's RealMLP and TabM), in its own module worker.
- **Web workers** (`src/workers/`): t-SNE and UMAP (manifold learning); the Crux module runs in its own worker (`src/crux/`).
- **TypeScript API** (`src/`): High-level interfaces, UI components, and integration with the Datagrok platform.

### Core Modules

#### Dimensionality Reduction
- **PCA** (`eda-tools.ts` `computePCA`): `ml.pca` (`svd_solver: 'auto'`, `random_state: 0`) on inputs standardised in TS (`standardizedRows` in `utils.ts`: population sigma; constant columns excluded and zero components appended). `center` / `scale` never touch the inputs: they apply to the **score columns** (`centerScaleDataFrame`), as in the pre-Crux NIPALS path - so `center = false` is not an uncentred PCA. Crux's PCA flips signs so the largest |loading| is positive (sklearn `svd_flip`)
- **UMAP/t-SNE** (`workers/umap-worker.ts`, `workers/tsne-worker.ts`): Manifold learning in workers
- Uses `@datagrok-libraries/ml` for multi-column dimensionality reduction with custom distance metrics

#### Partial Least Squares (`pls/`)
- **PLS regression and analysis** (`pls-tools.ts`, `pls-ml.ts`): Multivariate analysis for high-dimensional data
- `getPlsAnalysis` is one `ml.pls_analysis` call (component-major scores / loadings, `coef`, `vip`). Crux (scikit-learn) scales x and y by the sample std and does not normalise the scores; EDA keeps its own convention (doi.org/10.1002/cem.2589, Algorithm 2): x scores `t / ||t||`, x and y loadings `p, q * sqrt(n / (n - 1)) * ||t||`, y scores `U = y q^T / ||q||^2` from the population-standardised y, computed in TS. `coef` and VIP need no conversion. The prediction is `x . coef` without the intercept (`addMvaResults` adds the bias). Explained variances come from `q^2 / n` and `p^2 / n`
- `ml.pls_analysis` errors on data that cannot support `n_components` (`only k of A PLS components can be extracted; reduce n_components`, `y is uncorrelated with x ...`, `x is constant`, `y is constant`)
- Supports both linear and quadratic PLS models
- **Project persistence — do not revert:** the multivariate analysis never adds columns directly. `performMVA` reserves the result names (`getMvaNames`), then calls the `multivariateAnalysisTransform` package function (`meta.role: transform`) through a FuncCall with `{processed: false}`, and only then creates the viewers. That call goes into the table creation script, so a data-synced project replays it and restores the columns the saved layout refers to; without it the restored viewers crash on missing columns. The transform takes **column names**, not `column_list`: a replayed creation script passes such values as strings. The model (features, coefficients, bias) is stored in the `mvaModel` table tag, and the model scatter plot carries `initializationFunction: 'mvaModelInitFunction'`, which rebuilds the prediction-formula tooltip on project/layout opening. Regression tests: `category 'Multivariate analysis: projects'`.

#### ML methods on Crux (`src/crux/methods.ts`, `regression.ts`, `pls/pls-ml.ts`)
- Each `train*` in `package.ts` is one call into a method module, which maps the hyperparameters to a `cruxFit` call; each `apply*` is `cruxPredict(df, model)`: `Prediction` first, then one FLOAT `proba_<class>` column per class for every classifier but SVM. The model blob is the crux-model-v1 file with the feature names (`eda.features`); nothing else is stored. `mlname`, roles and the `train*` inputs are the Forge contract (Forge builds the hyperparameter form from them)
- **Legacy models** (`src/legacy/`): `cruxPredict` routes a blob that is not crux-model-v1 (`isCruxModel`) to `legacyPredict(df, bytes)`, which applies the pre-Crux containers with stored feature names: PLS Regression (any version with `Feature names`, with or without the `Variable Importance` column), Linear Regression from 1.8.0, Softmax, XGBoost (container v2) and SVM from 1.9.0. Prediction only: one `Prediction` column, no `proba_*`. Linear Regression / PLS (old f32 dot product) and Softmax (f32 standardisation, `exp` of the logits, argmax) are the old TS loops; XGBoost converts the stored XGBoost 3.3 UBJSON to JSON text (`ubjson.ts`) and scores it with Crux's `ml.xgboost_predict`, decoding objective 0 / 1 / 2 as before; SVM parses the LIBSVM text and predicts in TS (`libsvm.ts`, `c_svc` / `nu_svc` voting, `epsilon_svr` / `nu_svr`; only `c_svc` / `epsilon_svr` are fixture-tested) on f32-standardised rows, on the main thread (O(rows x SVs), as the old wasm). XGBoost / SVM regression on INT / BIG_INT targets is rounded as before; BOOL targets give BOOL. Rows with missing features are predicted by Linear Regression / PLS / Softmax / SVM from the DG null sentinel values, as the old code did (Crux models give null there; legacy XGBoost gets NaN, its missing value). One cache for the last blob (`preparedLegacyModel`, keyed by the blob's `blobId`) holds the prepared model: the container plus the XGBoost JSON text or the parsed LIBSVM model. A malformed inner model (UBJSON or LIBSVM) and an unknown XGBoost objective also throw `LEGACY_MODEL_ERROR`. Features are matched by the stored names only; blobs without names (Linear Regression before 1.8.0, Softmax / XGBoost before 1.9.0), XGBoost v1 containers and anything unrecognised throw `LEGACY_MODEL_ERROR`. `visualizePLSRegression` rebuilds the viewers for a legacy PLS blob from its component count. `src/crux/model-common.ts` holds what both `crux-model.ts` and `legacy-models.ts` use (`LEGACY_MODEL_ERROR`, `PREDICTION`, `featureColumns`, `blobId`, `runCrux`); the only import between the two is `cruxPredict` -> `legacyPredict`. **Do not revert:** the predictors keep the old `Float32Array` steps and loop order, and the SVM keeps libsvm's summation order - the parity test against `files/legacy-models/reference.json` is exact except `xgboost-regression` / `svm-regression` (FLOAT values within 1e-5 x max(1, |old|); INT after rounding, so exact in practice)
- `preprocess: 'standard'` (Crux's standard scaler, stored in the model's pipeline) for Linear Regression, Softmax, SVM, k-Nearest Neighbors and Neural Network; none for the tree methods, Naive Bayes, LDA (scale-invariant) and PLS Regression (Crux's PLS always scales), RealMLP and TabM (they scale inside the model). NaN features are allowed only for the seven `mlnativenans` methods (below)
- **Linear Regression** (`regression.ts`): L1 = L2 = 0 -> `linear` (exact OLS; Crux drops singular values below 1e-6 of the largest, which keeps one-hot designs finite), L1 = 0 -> `ridge {alpha: L2}`, else `elasticnet {alpha: L1 + L2, l1_ratio: L1 / (L1 + L2), max_iter: iterations}` (scikit-learn's definitions). Target INT or FLOAT
- **Softmax** (`methods.ts`): `logistic {C, max_iter: iterations, tol: tolerance}` with `C = 1 / penalty` for 3+ classes and `2 / penalty` for two classes (the old gradient-descent penalty on the two-row weights); `penalty` must be > 0 (`validate`)
- **XGBoost** (`methods.ts`): `xgboost {n_estimators: iterations, learning_rate: eta, max_depth, reg_lambda, reg_alpha}`; classification for STRING / BOOL targets, regression for INT / FLOAT (objective by class count)
- **SVM** (`methods.ts`): `svc` / `svr` with `kernel` from `{Linear: linear, Polynomial: poly, RBF: rbf, Sigmoid: sigmoid}`, `C`, `degree`, `coef0`, `svr_epsilon` (regression only); `gamma = 0` is sent as `1 / n_features` (libsvm's default). `isInteractive` stays far below XGBoost (SMO is superlinear). SVM returns `Prediction` alone, no `proba_*` columns: crux's `svc` gives probabilities only with `probability: true`, which EDA never sends (CRUX-ISSUES 21)
- **PLS Regression** (`pls/pls-ml.ts`): `trainPlsRegression` first makes the pre-Crux input checks (`checkDimensionReducerInputs` on the features, `checkColumnType` on the target) and runs `ml.pls_analysis` directly with the same `n_components` on the rows with a target (`targetRows`; an all-missing target throws `... has only missing values.`), rethrowing any error as `PLS Regression: ...` (`ml.fit pls` itself returns garbage or a NaN error when the data cannot support the components); then `ml.fit pls {n_components}` on the original rows. `visualizePLSRegression` rebuilds the viewers (loadings, coefficients, VIP, explained variance, scores with Hotelling ellipses) from a fresh `getPlsAnalysis` on the model's feature columns (`modelFeatures`: the manifest's `eda.features`) of the table Forge passes, on the rows with a target (`targetRows`), with `n_components` from the manifest
- Regression predictions of Crux models are FLOAT for every target type (INT targets are no longer rounded; legacy XGBoost / SVM models keep the old rounding). Rows with a missing target (numeric or string) are skipped at training (`training.n_rows`): by Crux's `ml.fit`, and for PLS Regression by EDA too (`targetRows` in `pls-ml.ts`), since `ml.pls_analysis` rejects a missing y; missing features are rejected except by the `mlnativenans` methods (below)
- XGBoost and SVM keep the pre-Crux input checks Crux does not make (`validate` in `methods.ts`, old messages): XGBoost iterations >= 1, eta in (0, 1], max depth >= 1, lambda / alpha >= 0, at least 2 classes; SVM degree a positive integer. Class counts ignore missing values (`classCount`, also Softmax's `C`)
- **Tree, neighbor, Bayes, LDA and neural-network methods** (`methods.ts`, crux model in brackets; `train*` = `cruxFit(taskOf(target), ...)`): **Random Forest** (`random_forest`) / **Extra Trees** (`extra_trees`) `{n_estimators: trees, max_depth, min_samples_leaf, max_features}`, `maxFeatures auto` sent as `sqrt` (classification) / `all` (regression); **Hist Gradient Boosting** (`hist_gradient_boosting`) `{n_estimators: iterations, learning_rate, max_leaf_nodes, min_samples_leaf, l2_regularization, early_stopping: 'off'}` plus `max_depth` only when > 0 (`auto` early stopping would hold out 10% of the rows above 10k); **CatBoost** (`catboost`, CatBoost's own kwargs) `{iterations, learning_rate, depth, l2_leaf_reg}`; **Decision Tree** (`decision_tree`) `{max_depth, min_samples_leaf}`; **k-Nearest Neighbors** (`knn`, `preprocess: 'standard'`, `mlhasrows`) `{n_neighbors, weights}`; **Naive Bayes** (`gaussian_nb`, no hyperparameters, class targets); **LDA** (`lda`, class targets) `{shrinkage}` only when > 0; **Neural Network** (`mlp`, `preprocess: 'standard'`) `{hidden_layer_sizes, activation, max_iter, learning_rate, alpha, tol: 1e-4}`. Defaults bound the model size (`maxDepth 10` for the forests and the tree, CatBoost 100 iterations)
- **RealMLP** (`methods.ts`, crux `realmlp`, pytabkit 1.7.3; no preprocess: it median-centres, robust-scales and smooth-clips itself) `{variant, n_epochs, hidden_sizes, batch_size, val_fraction}` plus `lr` only when > 0 (0 keeps RealMLP's 0.04 classification / 0.2 regression); `hidden_sizes` is a number array (`CruxFitParams.estimator` values are number, string, boolean or `number[]`). Holds out `val_fraction` of the rows it gets to pick the best epoch (upstream `n_cv 1`, `n_refit 0`) and trains on the rest; `training.n_rows` still counts every row. Rejects NaN features. Defaults TD, 128 epochs (upstream 256: the same quality on the probed tables, half the time), `256,256,256`, batch 256, fraction 0.2
- **TabM** (`methods.ts`, crux `tabm`, TabM 0.0.3 + rtdl_num_embeddings 0.0.12; no preprocess: quantile transform inside) `{k, d_block, embeddings, max_epochs, patience, lr, dropout, val_fraction}`; `mlnativenans` (a NaN becomes 0 after the quantile transform). EDA defaults k 8, width 128, 32 epochs instead of the paper's 32 / 512 / until early stopping, for speed in the wasm worker (the upstream defaults take 69-301 s at 1,000 x 10); `val_fraction` of the rows is held out for early stopping. TabM ignores unknown estimator keys silently (CRUX-ISSUES 17), so the manifest test pins every key
- `seed` (decorator default 0) on Random Forest, Extra Trees, Hist Gradient Boosting, CatBoost, Decision Tree, Neural Network, RealMLP and TabM goes through `cruxFit`'s `seed`
- `meta.mlnativenans: 'true'` on the `train` of the methods whose crux model takes NaN features: XGBoost, Random Forest, Extra Trees, Hist Gradient Boosting, CatBoost, Decision Tree, TabM (Forge then keeps those rows by default; an explicit Skip rows / Impute still applies); every other method rejects a NaN training row (`... missing feature values ...`) and leaves NaN rows null at predict
- No separate Logistic Regression / Ridge / Lasso / Elastic Net: Softmax is crux `logistic`, Linear Regression picks `linear` / `ridge` / `elasticnet`
- Checks Crux does not make (`validate`): kNN `neighbors` <= the rows with a target of the frame it gets (Forge trains per fold, so the message names the fold's rows; Crux would fit and fail at predict); CatBoost iterations >= 1, depth 1-16, learning rate in (0, 1], L2 >= 0 (Crux accepts iterations 0 and depth 17); Neural Network hidden layers as comma-separated positive integers, L2 >= 0; `alpha` always sent (Crux's classifier default is 1.0, not 1e-4); RealMLP epochs >= 1, batch size >= 1, learning rate >= 0, hidden layers text, validation fraction 0.05-0.5 (Crux accepts 0 epochs and fraction 0: an untrained model); TabM ensemble size / width / epochs / patience >= 1, learning rate > 0, dropout in [0, 1), validation fraction 0.05-0.5
- Thresholds (`cells = rows x features`; from `Crux methods benchmark`): Random Forest applicable `cells <= 1e6`, Extra Trees `<= 2e6`, both interactive `<= 2e4`; Hist Gradient Boosting / CatBoost / Decision Tree / Neural Network applicable `<= 1e7`; interactive HGB `rows <= 20000 && cells <= 2e5`, CatBoost `rows <= 20000 && cells <= 1e5`, Decision Tree `cells <= 2e5`, Neural Network `rows <= 2000 && features <= 50`; kNN applicable `rows <= 50000 && cells <= 2.5e6`, interactive `rows <= 10000 && cells <= 1e5`; Naive Bayes / LDA applicable `cells <= 1e7`, interactive `rows <= 100000` (LDA also `features <= 100`); RealMLP applicable `rows <= 10000 && cells <= 1e5`, never interactive (0.7-1 s per fit on iris, Forge runs six per retrain; ~45 s per fit on winequality); TabM applicable `rows <= 20000 && cells <= 2e5`, interactive `rows <= 200 && cells <= 2000`
- **Do not revert:** prediction matches features by NAME (`eda.features`), never by position: the platform does not hand `apply` the columns in the training order. `model-serialization-tests.ts` asserts reversed / rotated / superset tables and a missing column

#### Missing Values Imputation (`missing-values-imputation/`)
- **K-Nearest Neighbors** imputation (`knn-imputer.ts`)
- Handles both numerical and categorical features

#### ANOVA (`anova/`)
- **One-way ANOVA** (`anova-tools.ts`, `anova-ui.ts`) — two methods: **Fisher** (classical, requires equal variances) and **Welch** (default, robust to unequal variances).
- API: `oneWayAnova(cats, vals, alpha, {method, toValidate})` returns a discriminated union by `method`.
- **Numerical guarantees — do not revert:**
  - `FactorizedData.setStats` uses Welford accumulators (`means` + `m2`) instead of the naive `Σx²` form. Catastrophic cancellation otherwise on offset data.
  - p-values use `fSurvival(f, df1, df2)` → `jStat.ibeta(...)`. Never use `1 - jStat.centralF.cdf(...)` — it collapses to 0 in the tail.
  - Canary tests in `anova-tests.ts` (`Welford guard`, `ibeta guard`) will fail if either is rolled back.
- **NIST tests** (`category 'ANOVA: NIST StRD'`): fixtures in `src/tests/anova-fixtures.ts`, regenerated via `tools/generate-anova-fixtures.py` (requires scipy ≥ 1.16 + network to itl.nist.gov). Tolerances for SmLs07/09 are intentionally loose — Float64 Welford loses precision on 11-constant-digit stiff data; tests document observed accuracy, not gold-standard.
- **UI test caveat:** `DG.Column.fromStrings(['1','2',...])` auto-promotes purely numeric labels to an INT column, which breaks setStats indexing. Tests prefix labels with `g` (`factorCol` helper).

#### Control comparisons (`control-comparisons/`)
- **Many-to-one** comparison of each group against a single control (`control-comparisons-tools.ts`, `control-comparisons-ui.ts`). Two methods: **Dunnett** (default; pooled variance, joint reference distribution) and **Holm-Welch** (pairwise Welch's t-tests + Holm step-down over the k−1 control-vs-each p-values).
- API: `controlComparisons(cats, vals, controlCode, uniqueCount, {method, alpha})` → `ControlComparisonsReport`. Reuses the ANOVA `factorize`; Holm-Welch reuses `twoSampleTTestFromStats` (the package's Welch kernel) — no second Welch implementation.
- **Numerical guarantees — do not revert:**
  - `dunnettMaxCdf` integrates over the scale `s = S/σ ~ chi(df)/√df`. The outer window **must** be the adaptive Wilson–Hilferty bracket, never a fixed `[1e-3, 6]` grid: the peak narrows as `~1/√(2·df)`, so a fixed grid under-resolves it at large df and the adjusted p-values become garbage (e.g. raw p ≈ 2e-9 but adj p ≈ 0.2 on demog WEIGHT-by-RACE). Regression: `category 'Control comparisons: large df'` (union-bound invariant `raw ≤ adj ≤ m·raw`).
  - Holm is applied to exactly the **k−1** control-vs-each p-values, not all k(k−1)/2 pairs. Guard: `category 'Control comparisons: Holm family size'`.
  - One Γ-exact `hedgesCorrection` (from `ttest-tools.ts`) for **both** methods — do not mix in a second effect-size formula.
- **Fixtures** (`src/tests/control-comparisons-fixtures.ts`, hand-maintained): `FIXTURES` are synthetic/published cases (scipy.stats.dunnett, scipy ttest_ind, statsmodels Holm). `DEMOG_FIXTURES` cross-validate against real data — the test reads `System:DemoFiles/demog.csv` (`DEMOG_FILE`) so nothing large ships in the package. Caveats baked into the demog tolerances: Datagrok stores FLOAT columns as **float32** (so demog WEIGHT/HEIGHT agree only to ~1e-5; INT AGE is exact), and a group column with nulls loads as a literal `''` category (so only no-null group columns RACE/SEX/DIS_POP are used).

#### Pareto Optimization (`pareto-optimization/`)
- **ParetoOptimizer** (`pareto-optimizer.ts`): Multi-objective optimization
- **ParetoFrontViewer** (`pareto-front-viewer.ts`): Custom JsViewer for visualizing Pareto fronts
- Computes Pareto-optimal solutions from multiple objectives (minimize/maximize)

#### Probabilistic Scoring (`probabilistic-scoring/`)
- **pMPO** (`prob-scoring.ts`): Probabilistic Multi-Parameter Optimization for drug discovery
- Based on https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/
- Features:
  - Training with descriptor statistics and correlation filtering
  - Sigmoid-corrected desirability functions
  - ROC curve analysis and confusion matrix
  - Model evaluation and auto-tuning via Nelder-Mead optimization (`nelder-mead.ts`)
  - Export to MPO desirability profiles (integration with `@datagrok-libraries/statistics`)
- Sample data in `files/` directory for testing and demos

#### Crux bridge (`src/crux/`)
- Runs Crux functions (`ml.fit`, `ml.predict`, `ml.pca`, `ml.pls_analysis`, the metric functions) from EDA
- **Module files**: `crux/` (`crux-module-worker.js`, a verbatim copy of the crux-core worker, `.js` because the server's MIME map has no `mjs`; `eda/crux_wasm.js` + `eda/crux_wasm_bg.wasm` = `module-eda` (the scikit-learn registry plus XGBoost, CatBoost, RealMLP and TabM; 54 functions, 16.4 MB), committed; `eda/ENGINE_REV` and `engine.toml` = the pin: crux-core commit, feature list, build settings; build recipe in `crux/README.md`). `rspack.config.js` copies `crux/` verbatim (`info.minimized`, so the minimizer leaves them alone; without `README.md` and `engine.toml`) into `dist/crux/`, which `grok publish` uploads
- **Loader** (`crux-loader.ts`): `CruxModule` runs one module in a module worker (the `id/op/args` protocol of crux-core `crux-modules.mjs`; ops `init`, `callFn`, `fnResource`, `fnDescribe`; `callFn(name, ipc, params, resource?)` transfers `ipc`, so the caller's buffer is detached afterwards). Created lazily on the first call; URLs resolve from the running bundle (`document.currentScript`), so it works in `package-test.js` too. A load failure throws `Crux module "eda" failed to load from <url> (see crux/README.md)` and kills the worker; `terminate()` kills it on demand (the next call starts a new one). `callFn` registers a resource (the model) only when its id differs from the one live in the worker, and posts the registration and the call in one synchronous segment; a failed registration drops the resource entry, so the next call registers it again instead of running the previously live model. Crux memoises every registered model and never frees it, so after `MAX_LIVE_MODELS = 32` registrations the next call with no call pending restarts the worker. `cruxMl` is the `eda` module instance (the `ml.*` registry)
- **Bridge** (`crux-bridge.ts`): features go as one Arrow column `x: List<Float64>` built from `getRawData()` (row-major `Float64Array` + offsets, no per-row arrays); only INT / FLOAT features, DG nulls -> `NaN` (Crux's missing value). Target `y`: numeric -> Float64 (nulls -> `NaN`), STRING -> Utf8 from the category codes of `getRawData()` (`columnTexts`; the `''` category, i.e. a missing value, -> an Arrow null, never a valid `""` class), BOOL -> Bool. Crux encodes string / boolean classes itself and returns `class` as text. `toArrowIpc(features, target?)` = `matrixIpc(rowMajor, rows, cols, target?)`, the entry for a prepared matrix (PCA standardises in TS); `vectorsIpc(columns)` packs named vectors (`utf8Vector`, or arrow's `makeVector` of a typed array). Outputs: `scalar(table, name)` is the first value of a numeric column (undefined when null or not finite), `numbers(table, name)` its values (nulls stay null, infinities stay); `flatList(table, name)` reads a one-row `List<Float64>` cell (`ml.pls_analysis`), `flatBatch(table, name)` the concatenated child of a `List<Float64>` column with one list per row (`ml.pca`); both copy with `Float64Array.from` only when the values are not a `Float64Array`. `listColumns(table, name, names)` turns a `List<Float64>` column with one list per row into one FLOAT column per name (`ml.predict`'s `probas`): `flatBatch` when no row is null, per row otherwise; a null row is null in every column, non-finite values become DG nulls
- **Model** (`crux-model.ts`): `cruxFit(task, model, df, predictColumn, {estimator, preprocess, seed})` trains on all columns of `df` (`featurizer: vector`, `cv: 0`; `seed` goes to both `random_state` and `estimator.random_state`) and returns a crux-model-v1 (safetensors) blob with the DG feature names added as `__metadata__['eda.features']` (`withEdaMetadata`). A failed fit rejects with Crux's message (`ml.fit: ...`). `cruxPredict(df, bytes)` picks the features by those names (`modelFeatures(df, manifest)`) (any order, extra columns ignored, a missing one throws `... is missing`; a blob without `eda.features` throws `... has no feature names`, there is no positional fallback) and returns a DataFrame with one `Prediction` column: FLOAT for regression; for classification BOOL / STRING by the JSON type of the manifest's `classes`, numeric classes as INT when all are integers and FLOAT otherwise; for a classifier whose model gives probabilities (`hasProbabilities`, decided by the manifest, not the rows: every classifier but `svc`, whose `ml.predict` returns a `probas` column of null lists), one FLOAT column `proba_<class>` per manifest class follows, in the manifest's order, null where `Prediction` is null (also when every row is null, and for 0 rows) (`listColumns` in the bridge); none for regression, and none when two classes differ only in case (`DG.DataFrame.fromColumns` throws `Column named ... already exists` for such names). Rows Crux cannot predict (e.g. NaN features for a model without NaN support) stay null (`false` in a BOOL column). `readCruxManifest(bytes)` throws `LEGACY_MODEL_ERROR` ("... retrain the model") for any blob that is not crux-model-v1 (`isCruxModel`); `cruxPredict` sends such blobs to `legacyPredict` (Legacy models above), which throws the same error for anything it does not recognise. `cruxTerminate()` stops the worker (`cruxMl.terminate()`). `runCrux(name, ipc, params, resource?)` is the one call helper (`callFn` + `tableFromIPC`; defined in `model-common.ts`, re-exported). Model ids for the resource cache are the SHA-256 of the blob (`blobId`)
- **Metrics** (`crux-metrics.ts`): `classificationReport(actual, predicted, positiveClass?, probabilities?, curves)`, exposed as the headless `Eda:classificationReport` (`object` output; Forge calls it for every classifier). Over the rows where both columns have a value (classes compared as text: STRING categories with `''` as missing, BOOL as `true` / `false`, numbers as `String(v)`): accuracy, balanced accuracy, MCC, Cohen's kappa, precision / recall / F1 of `positiveClass` when the rows hold exactly two classes and it is one of them, else macro; the confusion matrix and each class's rows, precision, recall and F1; with `probabilities` (`{classes, columns}`) each class's one-vs-rest AUC-ROC and average precision, the overall ones (the positive class's, else the mean over the classes with both positive and negative rows) and, when `curves` and two classes with a positive class, its ROC and precision-recall curves (one point per distinct probability; the ROC's first threshold is +Infinity, the PR's last is null). `ml.classification_metrics` takes only integer labels: EDA encodes the classes to Int32 codes (`average: binary` with `pos_label` for two classes, `macro` otherwise; `binary` is refused for three). The labels and the per-class precision / recall / F1 come from one `ml.confusion_matrix` (both columns passed as text, labels = their sorted union, true-major). One `ml.roc_auc` per probability class on the labels binarized to 0 / 1 (never `pos_label`: an unknown one gives a curve of null rates instead of an error, CRUX-ISSUES 23); the curves (`ml.roc_curve` with `drop_intermediate: false`, `ml.precision_recall_curve`) on the same binarized labels. Rows with a null or non-finite probability are left out in TS (a NaN score is a Crux error). Package function: an empty `positiveClass` / `classes` and a null `probabilities` mean none; `classes` without `probabilities`, or of another length, throws `classificationReport: classes and probabilities differ in length`

### Function Registration

Functions are registered with `@grok.decorators.func()` / `@grok.decorators.param()` decorators:

```typescript
@grok.decorators.func({
  'top-menu': 'ML | Analyze | PCA...',
  'description': 'Principal component analysis (PCA)',
  'helpUrl': '/help/explore/dim-reduction#pca',
})
static async PCA(
  @grok.decorators.param({'type': 'dataframe', 'options': {'caption': 'Table'}}) table: DG.DataFrame,
  // ...
) { }
```

Run `grok api` to generate `package.g.ts` and `package-api.ts` from these annotations.

### External Dependencies

Key library integrations:
- **@datagrok-libraries/ml**: Distance metrics, dimensionality reduction, MCL clustering
- **@datagrok-libraries/statistics**: Statistical functions, MPO profile editor
- **@datagrok-libraries/math**: DBSCAN, mathematical utilities
- **@datagrok-libraries/test**: Test framework

Provided by the platform host (webpack `externals`, not bundled):
- `datagrok-api` (`dg`/`grok`/`ui`), `openchemlib/full.js`, `rxjs`, `cash-dom`, `dayjs`, `wu`, `exceljs`

Bundled into the artifact (not external):
- `mathjs` (matrix ops), `jstat` (statistical distributions), `umap-js`, `@keckelt/tsne`, `apache-arrow` (Crux IPC), `fast-sha256` (Crux model ids). The `crux/eda` module ships next to the bundle in `dist/crux/`. See `CREDITS.md` for attribution.

## Testing

Tests are organized by feature in `src/tests/`:
- `dim-reduction-tests.ts`: PCA, UMAP, t-SNE, SPE
- `linear-methods-tests.ts`: PCA (incl. the center / scale flags), PLS (incl. the rank error and PLS Regression on a target with missing values), linear regression (incl. the model choice), Hotelling ellipse
- `classifiers-tests.ts`: softmax, XGBoost, SVM (incl. the hyperparameters in the manifest, XGBoost's `Missing feature values are kept`, class probabilities; SVM returns none)
- `crux-methods-tests.ts`: Random Forest, Extra Trees, Hist Gradient Boosting, CatBoost, Decision Tree, k-Nearest Neighbors, Naive Bayes, LDA, Neural Network, RealMLP and TabM, table-driven (`METHODS`, one category per method): iris round trips, manifest params, missing feature values (native vs rejected), applicability (incl. shapes at each cap), bool target (on the first `boolRows` demog rows when set: RealMLP 1,000), seed, class probabilities (classifiers, `expectProbabilities` of `utils.ts`), method-specific validation; plus `Crux methods benchmark` (`benchmark: true`, one test per method, `BENCH` lines on cars / iris / winequality (`CASES` / `loadCase` of `baseline.ts`) / demog with gaps, and under `--benchmark` only also on synthetic 100k x 10 (`CASES`) and 20k x 100, built on the first benchmark test): every model <= 20 MB always, every case the method calls interactive fits in <= 1 s under `--benchmark` (a warm-up fit absorbs the worker start)
- `tests/utils.ts` holds the shared helpers (`roundTrip`, `applyTwice`, `expectFeatureMatching`, `expectRejected`, `expectGapsKept`, `demogFrame`, `irisMeasurements`, `expectProbabilities`, `IRIS_SPECIES`)
- `model-serialization-tests.ts`: repeated apply of the saved model **and apply-time feature matching** (softmax, linear regression, PLS, SVM, XGBoost) on iris/cars/winequality/demog. `expectFeatureMatching` checks each model is invariant to reversed / rotated column order and to extra ("decoy") columns in the table, and throws when a feature is missing
- `mis-vals-imputation-tests.ts`: KNN imputation
- `anova-tests.ts`: One-way ANOVA
- `ttest-tests.ts`: Two-sample t-test
- `control-comparisons-tests.ts`: Control comparisons (Dunnett, Holm-Welch)
- `pmpo-tests.ts`: Probabilistic scoring
- `pareto-tests.ts`: Pareto optimization
- `legacy-models-tests.ts`: the pre-Crux fixtures in `files/legacy-models/` against the old implementation's predictions in `reference.json` (one test per fixture; exact except `xgboost-regression` / `svm-regression`, 1e-5 x max(1, |old|) for FLOAT values: `TOLERANCE`), a PLS blob without VIP, unknown blobs (incl. an XGBoost v1 header, a bare `Float32Array`, Softmax without names) rejected with the "retrain" error, the cached prepared model, feature matching by name, legacy PLS visualize, UBJSON to JSON (category `Legacy models`)
- `crux-bridge-tests.ts`: the Crux bridge (module load, column types, nulls / NaN, boolean, missing string and string targets, class probabilities, save and read back, seed, legacy blob, no feature names, failed model registration retried, live-model cap; benchmark: Arrow packing time) (category `Crux bridge`)
- `crux-metrics-tests.ts`: `classificationReport` (category `Crux metrics`: scores vs hand-computed values, the curves, three classes, boolean columns, rows without values, the package function, iris parity, a 100k benchmark)
- `crux-thresholds-tests.ts` + `baseline.ts`: the Crux methods (`FITS`: the `train*` defaults) against the pre-Crux numbers in `files/baselines/baseline.json` on the same 5 folds (category `Crux thresholds`): XGBoost / SVM metrics within 1% (integer regression targets compared on rounded predictions, as the old methods rounded), Softmax accuracy not lower by more than 0.02, Linear Regression / PLS Regression R² not lower (1e-5; `linreg/cars-price`, 24 rows x 15 features, compares the training-set R²), PLS analysis and PCA outputs within 1e-6 up to a per-component sign; the time tests (`benchmark: true`) assert XGBoost / SVM fit time <= 2x (only cases with a baseline `fitMs` >= 100 ms; smaller ones are logged as informational) only under `DG.Test.isInBenchmark` (`--benchmark`; otherwise they only log the times, wall-clock is machine-dependent). `benchmark: true` tests also run without `--benchmark`; `--benchmark` runs only them

Test framework from `@datagrok-libraries/test`:
```typescript
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Feature Name', () => {
  test('Test description', async () => {
    // test code
  }, {timeout: 30000});
});
```

## Files

- **`files/`**: Demo/test data for pMPO (drugs-props-train.csv, drugs-props-test.csv, scores)
- **`files/legacy-models/`**: pre-Crux models (`.bin`: every method, plus regression and bool targets of XGBoost / SVM / Softmax) and `reference.json`, the old implementation's predictions for them; produced with the old wasm and not regenerable after it was removed; see its `README.md`. Published to `System:AppData/Eda/legacy-models/`
- **`files/baselines/baseline.json`**: scores, fit times, PCA and PLS outputs of the pre-Crux methods on fixed datasets and folds (`src/tests/baseline.ts`), produced once with the old wasm; the reference for the Crux thresholds. See its `README.md`. Published to `System:AppData/Eda/baselines/`
- **`css/`**: Custom stylesheets (pmpo.css for pMPO UI)

## Key Workflows

### Adding a New ML Method

1. Add a `train*` wrapper over `cruxFit` (hyperparameter map, `preprocess`) in `src/crux/methods.ts`; `apply*` is `cruxPredict`
2. Add function registration with `@grok.decorators.func()` in `package.ts` (`meta.mlname`, `train` / `apply` / `isApplicable` / `isInteractive` roles)
3. Add a row to `METHODS` in `crux-methods-tests.ts` and to `FITS` in `baseline.ts`; `meta.mlnativenans: 'true'` if the crux model takes NaN features; if the method replaces an existing one, compare it against a baseline as `Crux thresholds` does
4. Run `grok api` to generate wrappers
5. Add to package menu structure in `package.json` meta.menu

When copying `col.getRawData()` in bulk, ALWAYS slice to `col.length` first (`raw.subarray(0, col.length)`) - the raw array may be longer (capacity mechanism). A new Crux module revision: `crux/README.md`, then `grok build --force`.

### pMPO Model Development

The probabilistic scoring module is complex and tightly integrated:
- Model training in `Pmpo` class (`prob-scoring.ts`)
- Statistical computations in `stat-tools.ts`
- UI utilities and serialization in `pmpo-utils.ts`
- Optimization via Nelder-Mead in `nelder-mead.ts`
- All constants and types in `pmpo-defs.ts`
- Integration with `@datagrok-libraries/statistics` for MPO profile export
