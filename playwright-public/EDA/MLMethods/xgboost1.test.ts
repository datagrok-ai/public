import { test, expect } from '@playwright/test';
import {
  openDemoCsv, openTrainModelView, resetShell, setPredictColumn, trainEdaModelViaApi,
} from '../helpers';

// Test Track scenario: EDA/MLMethods/xgboost1.md
// 1. Open iris.csv from Demo files.
// 2. Top Menu > ML > Models > Train Model...
// 3. Predict = Species, Features = all except Species and col 1, Model Engine = Eda: XGBoost.
// 4. In the model result, vary Rate/Lambda/Alpha sliders and Iterations/Max depth clickers.
//
// UI scope and fallbacks (per xgboost1-run.md):
//   * UI: file open, Train Model view, Predict column.
//   * API: training, since the Train Model UI has no Model Engine dropdown or Train button
//     on dev and the canvas-based Features picker cannot toggle individual columns.
//   * Hyperparameter sliders live in the PredictiveModel result view; the API call returns
//     the raw model blob and does not open that view. Hyperparameter assertions are out of
//     scope here — they need the result-view to be openable from automation.

test.describe.serial('EDA / MLMethods / XGBoost (classification)', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('Train XGBoost classification on iris.csv predicting Species', async ({ page }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, 'iris.csv');
    await openTrainModelView(page);
    await setPredictColumn(page, 'Species');

    const trained = await trainEdaModelViaApi(page, 'eda:trainXGBooster', 'Species');
    expect(trained.ok, trained.error ?? 'xgboost training returned null').toBe(true);
  });
});
