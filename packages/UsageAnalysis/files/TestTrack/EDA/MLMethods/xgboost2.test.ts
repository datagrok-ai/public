import { test, expect } from "./helpers";
import {
  openDemoCsv,
  openTrainModelView,
  resetShell,
  setPredictColumn,
  trainEdaModelViaApi,
} from "./helpers";

test.describe.serial("EDA / MLMethods / XGBoost (regression)", () => {
  test.afterEach(async ({ page }) => {
    await resetShell(page);
  });

  test("Train XGBoost regression on cars.csv predicting price", async ({
    page,
  }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, "cars.csv");
    await openTrainModelView(page);
    await setPredictColumn(page, "price");

    const trained = await trainEdaModelViaApi(
      page,
      "eda:trainXGBooster",
      "price",
      { numericOnly: true },
    );
    expect(trained.ok, trained.error ?? "xgboost training returned null").toBe(
      true,
    );
  });
});
