# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

EDA (Exploratory Data Analysis) is a Datagrok package providing statistical analysis and machine learning tools. It includes dimensionality reduction (PCA, UMAP, t-SNE, SPE), supervised learning (SVM, linear regression, softmax, PLS, gradient boosting via XGBoost), ANOVA, missing data imputation, Pareto optimization, and probabilistic multi-parameter optimization (pMPO).

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

- **WASM Modules** (`wasm/`): C++ implementations compiled to WebAssembly for PCA, PLS, SVM, softmax regression, and XGBoost. These provide high-performance numerical computation.
- **Web Workers** (`src/workers/`): JavaScript workers for t-SNE, UMAP, and softmax training to avoid blocking the UI thread.
- **TypeScript API** (`src/`): High-level interfaces, UI components, and integration with Datagrok platform.

### Core Modules

#### Dimensionality Reduction
- **PCA** (`eda-tools.ts`, `wasm/PCA/`): Principal Component Analysis via WASM
- **UMAP/t-SNE** (`workers/umap-worker.ts`, `workers/tsne-worker.ts`): Manifold learning in workers
- Uses `@datagrok-libraries/ml` for multi-column dimensionality reduction with custom distance metrics

#### Partial Least Squares (`pls/`)
- **PLS regression and analysis** (`pls-tools.ts`, `pls-ml.ts`): Multivariate analysis for high-dimensional data
- **WASM backend** (`wasm/PLS/`): Performance-critical PLS computations
- Supports both linear and quadratic PLS models

#### Support Vector Machines (`svm.ts`)
- **LS-SVM** implementation with multiple kernels (linear, RBF, polynomial, sigmoid)
- Training/prediction via WASM (`wasm/svm.h`, `wasm/svmApi.cpp`)
- Interactive model training with progress tracking

#### XGBoost (`xgbooster.ts`)
- **XGBoost** integration via WASM (`wasm/XGBoostAPI.wasm`)
- Binary classification and regression
- Requires initialization via `initXgboost()` in package init

#### Other Supervised Machine Techniques
- **Softmax Classifier** (`softmax-classifier.ts`): Multinomial logistic regression
  - Training in web worker for non-blocking UI
  - WASM backend for prediction
- **Linear Regression** (`regression.ts`): Ordinary least squares

#### Missing Values Imputation (`missing-values-imputation/`)
- **K-Nearest Neighbors** imputation (`knn-imputer.ts`)
- Handles both numerical and categorical features

#### ANOVA (`anova/`)
- **One-way ANOVA** (`anova-tools.ts`, `anova-ui.ts`)
- Statistical analysis of variance with visual reports

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

WASM modules are initialized asynchronously in `PackageFunctions.init()`:
- `_initEDAAPI()` from `wasm/EDAAPI.js`
- `initXgboost()` from `wasm/xgbooster.js`

WASM source files (`wasm/*.cpp`, `wasm/*.h`) are C++ implementations compiled to WebAssembly. Do not modify generated `.js` and `.wasm` files directly.

### Function Registration

Functions are registered via JSDoc-style metadata comments (not TypeScript decorators):

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

External packages (must be in webpack externals):
- `mathjs`: Matrix operations and numerical computation
- `jstat`: Statistical distributions
- `umap-js`: UMAP implementation
- `@keckelt/tsne`: t-SNE implementation

## Testing

Tests are organized by feature in `src/tests/`:
- `dim-reduction-tests.ts`: PCA, UMAP, t-SNE, SPE
- `linear-methods-tests.ts`: Linear regression, PLS
- `classifiers-tests.ts`: SVM, softmax, XGBoost
- `mis-vals-imputation-tests.ts`: KNN imputation
- `anova-tests.ts`: One-way ANOVA
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

## Scripts and Files

- **`scripts/`**: Python scripts for exporting TypeScript constants from Python implementations
- **`files/`**: Demo/test data for pMPO (drugs-props-train.csv, drugs-props-test.csv, scores)
- **`css/`**: Custom stylesheets (pmpo.css for pMPO UI)

## Key Workflows

### Adding a New ML Method

1. Implement the algorithm in TypeScript (or C++ for WASM)
2. Add function registration with `@grok.decorators.func()`
3. Create tests in appropriate test file
4. Run `grok api` to generate wrappers
5. Add to package menu structure in `package.json` meta.menu

### Working with WASM

- WASM sources in `wasm/*.cpp` and `wasm/*.h`
- Generated files: `wasm/*.js`, `wasm/*.wasm`
- Initialization required in `PackageFunctions.init()`
- Web workers have separate WASM entry points (e.g., `EDAForWebWorker.js`)

### pMPO Model Development

The probabilistic scoring module is complex and tightly integrated:
- Model training in `Pmpo` class (`prob-scoring.ts`)
- Statistical computations in `stat-tools.ts`
- UI utilities and serialization in `pmpo-utils.ts`
- Optimization via Nelder-Mead in `nelder-mead.ts`
- All constants and types in `pmpo-defs.ts`
- Integration with `@datagrok-libraries/statistics` for MPO profile export
