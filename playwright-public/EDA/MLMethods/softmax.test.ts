import { test, expect } from '@playwright/test';
import {
  openDemoCsv, openTrainModelView, resetShell, setPredictColumn, trainEdaModelViaApi,
} from '../helpers';

// Test Track scenario: EDA/MLMethods/softmax.md
// 1. Open iris.csv from Demo files.
// 2. Top Menu > ML > Models > Train Model...
// 3. Predict = Species, Features = all except Species and col 1, One-hot encode strings on,
//    Model Engine = Eda: Softmax.
// 4. Vary Hyperparameters in the model result.
//
// `eda:trainSoftmax` uses every column of the passed `df` as a feature and requires them
// all to be numeric (it throws "incorrect features type" otherwise). Unlike XGBoost it does
// not tolerate the string target column. So training is invoked with `numericOnly: true`,
// which feeds only the numeric feature columns (excluding the string `Species` target) while
// `Species` is still passed as the predict column. UI scope and the API-training fallback
// match the other MLMethods tests (canvas-based Features picker, no Train button on dev).

test.describe.serial('EDA / MLMethods / Softmax', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('Train Softmax on iris.csv predicting Species', async ({ page }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, 'iris.csv');
    await openTrainModelView(page);
    await setPredictColumn(page, 'Species');

    const trained = await trainEdaModelViaApi(page, 'eda:trainSoftmax', 'Species', {
      numericOnly: true,
      extraParams: { rate: 0.1, iterations: 100, penalty: 0.01, tolerance: 0.001 },
    });
    expect(trained.ok, trained.error ?? 'softmax training returned null').toBe(true);
  });
});
