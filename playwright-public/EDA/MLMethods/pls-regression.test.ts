import { test, expect } from '@playwright/test';
import {
  openDemoCsv, openTrainModelView, resetShell, setPredictColumn, trainEdaModelViaApi,
} from '../helpers';

// Test Track scenario: EDA/MLMethods/pls-regression.md
// 1. Open cars.csv from Demo files.
// 2. Top Menu > ML > Models > Train Model...
// 3. Predict = price, Features = all except Price and Model, Model Engine = Eda: PLS Regression.
// 4. Verify viewers: Loadings scatter, Regression coefficients bar, Explained variances bar,
//    Scores scatter.
//
// UI scope and fallbacks (per pls-regression-run.md):
//   * UI: file open, Train Model view open, Predict column.
//   * API: training itself — the Train Model UI has no Model Engine dropdown for PLS on dev,
//     and the canvas-based Features picker cannot toggle individual columns. The function
//     `eda:trainPLSRegression` is the same path the UI would have invoked.
//   * Viewer presence: the API call returns a raw model blob without opening the
//     PredictiveModel result view, so the four expected viewers cannot be asserted here.
//     The test verifies training succeeds; viewer mounting is tracked as a separate
//     scenario gap (see run notes).

test.describe.serial('EDA / MLMethods / PLS Regression', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('Train PLS Regression on cars.csv with 3 components predicting price', async ({ page }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, 'cars.csv');
    await openTrainModelView(page);
    await setPredictColumn(page, 'price');

    const trained = await trainEdaModelViaApi(page, 'eda:trainPLSRegression', 'price', {
      numericOnly: true,
      extraParams: { components: 3 },
    });
    expect(trained.ok, trained.error ?? 'training returned null').toBe(true);
  });
});
