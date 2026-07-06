import { test, expect } from '../helpers';
import {
  openDemoCsv, openTrainModelView, resetShell, setPredictColumn, trainEdaModelViaApi,
} from '../helpers';

// Test Track scenario: EDA/MLMethods/xgboost2.md
// 1. Open cars.csv from Demo files.
// 2. Top Menu > ML > Models > Train Model...
// 3. Predict = price, Features = all except Price and Model, Model Engine = Eda: XGBoost.
// 4. In the model result, vary Rate/Lambda/Alpha sliders and Iterations/Max depth clickers.
//
// UI scope and fallbacks (per xgboost2-run.md): same constraints as xgboost1 — UI covers
// file open / Train Model view / Predict; training itself goes through the same registered
// function the UI would call, and the hyperparameter sliders live in a result view that
// the API path does not open.

test.describe.serial('EDA / MLMethods / XGBoost (regression)', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('Train XGBoost regression on cars.csv predicting price', async ({ page }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, 'cars.csv');
    await openTrainModelView(page);
    await setPredictColumn(page, 'price');

    const trained = await trainEdaModelViaApi(page, 'eda:trainXGBooster', 'price', { numericOnly: true });
    expect(trained.ok, trained.error ?? 'xgboost training returned null').toBe(true);
  });
});
