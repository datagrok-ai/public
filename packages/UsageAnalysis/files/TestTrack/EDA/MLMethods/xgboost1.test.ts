import { test, expect } from "./helpers";
import {
  openDemoCsv,
  openTrainModelView,
  resetShell,
  setPredictColumn,
  trainEdaModelViaApi,
} from "./helpers";

test.describe.serial("EDA / MLMethods / XGBoost (classification)", () => {
  test.afterEach(async ({ page }) => {
    await resetShell(page);
  });

  test("Train XGBoost classification on iris.csv predicting Species", async ({
    page,
  }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, "iris.csv");
    await openTrainModelView(page);
    await setPredictColumn(page, "Species");

    const trained = await trainEdaModelViaApi(
      page,
      "eda:trainXGBooster",
      "Species",
    );
    expect(trained.ok, trained.error ?? "xgboost training returned null").toBe(
      true,
    );
  });
});
