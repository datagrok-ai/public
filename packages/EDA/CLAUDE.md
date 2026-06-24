# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

EDA (Exploratory Data Analysis) is a Datagrok package providing statistical analysis and machine learning tools. It includes dimensionality reduction (PCA, UMAP, t-SNE, SPE), supervised learning (linear regression, softmax, PLS, gradient boosting via XGBoost), group comparison (two-sample t-test, ANOVA, control comparisons), missing data imputation, Pareto optimization, and probabilistic multi-parameter optimization (pMPO).

## Build Commands

```bash
npm install                # Install dependencies
npm run build              # Full build: grok api && grok check --soft && webpack
npm run build-eda          # Webpack only
npm run build-all          # Build js-api, libraries, then this package

npm run link-all           # Link local datagrok-api and libraries

npm run debug-eda          # Build and publish to default server
npm run debug-eda-local    # Build and publish to local server
npm run release-eda        # Build and publish as release

npm test                   # Run tests against localhost
grok test                  # Run all tests
grok test --test "TestName"      # Run specific test
grok test --category "CategoryName"  # Run tests in category
grok test --gui            # Run with visible browser
```

## Architecture

### High-Level Structure

The package combines TypeScript, WASM modules, and web workers for performance-critical computations:

- **WASM modules** (`wasm/`): PCA, PLS, softmax and linear regression run on the Rust + WebAssembly `sci-comp-ml` backend (`wasm/sci_comp_ml.js` glue + `sci_comp_ml_bg.wasm`), wrapped by `wasm/eda-api.ts` and loaded lazily via `src/wasm-loader.ts`. XGBoost remains a C++/Emscripten module.
- **Web workers** (`src/workers/`): t-SNE and UMAP (manifold learning), plus PCA and PLS fits (`eda-ml-worker.ts`) — kept off the UI thread. `softmax-worker.ts` is a TypeScript fallback used only when the WASM softmax fit throws (see `softmax-classifier.ts`); the primary softmax path runs on the Rust+WASM backend on the main thread.
- **TypeScript API** (`src/`): High-level interfaces, UI components, and integration with the Datagrok platform. Input standardisation happens here (the kernels never standardise — "variant B").

### Core Modules

#### Dimensionality Reduction
- **PCA** (`eda-tools.ts` → `wasm/eda-api.ts` → `workers/eda-ml-worker.ts`): NIPALS Principal Component Analysis on the Rust+WASM backend, fit in a web worker
- **UMAP/t-SNE** (`workers/umap-worker.ts`, `workers/tsne-worker.ts`): Manifold learning in workers
- Uses `@datagrok-libraries/ml` for multi-column dimensionality reduction with custom distance metrics

#### Partial Least Squares (`pls/`)
- **PLS regression and analysis** (`pls-tools.ts`, `pls-ml.ts`): Multivariate analysis for high-dimensional data
- **Rust+WASM backend** (`wasm/eda-api.ts` → `workers/eda-ml-worker.ts`): PLS1 fit off the UI thread
- Supports both linear and quadratic PLS models

#### XGBoost (`xgbooster.ts`)
- **XGBoost** integration via WASM (`wasm/XGBoostAPI.wasm`)
- Binary classification and regression
- Requires initialization via `initXgboost()` in package init

#### Other Supervised Machine Techniques
- **Softmax Classifier** (`softmax-classifier.ts`): Multinomial logistic regression
  - Training on the Rust+WASM backend (`_fitSoftmax` in `eda-api.ts`, main thread); prediction in TypeScript
- **Linear Regression** (`regression.ts`): Elastic Net (L1/L2) on the Rust+WASM backend; defaults to ordinary least squares and exposes gradient-descent hyperparameters

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

### WASM Integration

- **sci-comp-ml** (PCA, PLS, softmax, linear regression): a Rust + WASM crate built with `wasm-pack`, copied into `wasm/` (`sci_comp_ml.js` + `sci_comp_ml_bg.wasm`). Initialised lazily on first use via `src/wasm-loader.ts`; the binary is emitted to `dist/` (via the `.wasm` file-loader rule + `experiments.asyncWebAssembly` in `webpack.config.js`) and located relative to the running bundle (`document.currentScript`), **not** `_package.webRoot` — so it resolves in the separate `package-test.js` bundle too. Wrapped by `wasm/eda-api.ts`; PCA/PLS fits run in `src/workers/eda-ml-worker.ts`.
- **XGBoost**: a C++/Emscripten module (`wasm/XGBoostAPI.*`), initialised via `initXgboost()` in `PackageFunctions.init()`.

Do not modify generated `.js`/`.wasm` files directly. The Rust source lives in the separate `sci-comp-rust` repo; refresh the artefact by rebuilding there and re-copying into `wasm/`.

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
- `mathjs` (matrix ops), `jstat` (statistical distributions), `umap-js`, `@keckelt/tsne`, and the `sci_comp_ml` wasm. See `CREDITS.md` for attribution.

## Testing

Tests are organized by feature in `src/tests/`:
- `dim-reduction-tests.ts`: PCA, UMAP, t-SNE, SPE
- `linear-methods-tests.ts`: Linear regression, PLS
- `classifiers-tests.ts`: softmax, XGBoost
- `model-serialization-tests.ts`: model pack/unpack round-trip (softmax, linear regression, PLS) on iris/cars/winequality
- `mis-vals-imputation-tests.ts`: KNN imputation
- `anova-tests.ts`: One-way ANOVA
- `ttest-tests.ts`: Two-sample t-test
- `control-comparisons-tests.ts`: Control comparisons (Dunnett, Holm-Welch)
- `pmpo-tests.ts`: Probabilistic scoring
- `pareto-tests.ts`: Pareto optimization

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
- **`css/`**: Custom stylesheets (pmpo.css for pMPO UI)

## Key Workflows

### Adding a New ML Method

1. Implement the algorithm in TypeScript (or Rust in the `sci-comp-rust` repo for a WASM backend; XGBoost is the only remaining C++/Emscripten kernel)
2. Add function registration with `@grok.decorators.func()`
3. Create tests in appropriate test file
4. Run `grok api` to generate wrappers
5. Add to package menu structure in `package.json` meta.menu

### Working with WASM

- **sci-comp-ml** (PCA/PLS/softmax/linreg): Rust source lives in the separate `sci-comp-rust` repo; the built `sci_comp_ml.js` + `sci_comp_ml_bg.wasm` are committed in `wasm/`. Loaded via `src/wasm-loader.ts`; PCA/PLS fits run in `src/workers/eda-ml-worker.ts` (the worker `init`s the wasm from a URL passed by the main thread).
- **XGBoost**: C++/Emscripten (`wasm/XGBoostAPI.*`), initialised in `PackageFunctions.init()`.
- The `.wasm` file-loader rule + `experiments.asyncWebAssembly` in `webpack.config.js` emit `sci_comp_ml_bg.wasm` to `dist/` with its original name.

### pMPO Model Development

The probabilistic scoring module is complex and tightly integrated:
- Model training in `Pmpo` class (`prob-scoring.ts`)
- Statistical computations in `stat-tools.ts`
- UI utilities and serialization in `pmpo-utils.ts`
- Optimization via Nelder-Mead in `nelder-mead.ts`
- All constants and types in `pmpo-defs.ts`
- Integration with `@datagrok-libraries/statistics` for MPO profile export
