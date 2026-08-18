import { test, expect } from "./helpers";
import {
  openDemoCsv,
  openTrainModelView,
  resetShell,
  setPredictColumn,
  trainEdaModelViaApi,
} from "./helpers";

test.describe.serial("EDA / MLMethods / Softmax", () => {
  test.afterEach(async ({ page }) => {
    await resetShell(page);
  });

  test("Train Softmax on iris.csv predicting Species", async ({ page }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, "iris.csv");
    await openTrainModelView(page);
    await setPredictColumn(page, "Species");

    const trained = await trainEdaModelViaApi(
      page,
      "eda:trainSoftmax",
      "Species",
      {
        numericOnly: true,
        extraParams: {
          rate: 0.1,
          iterations: 100,
          penalty: 0.01,
          tolerance: 0.001,
        },
      },
    );
    expect(trained.ok, trained.error ?? "softmax training returned null").toBe(
      true,
    );
  });
});
