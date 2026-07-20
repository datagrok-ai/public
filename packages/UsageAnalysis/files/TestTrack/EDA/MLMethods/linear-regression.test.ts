import { test, expect } from "./helpers";
import {
  openDemoCsv,
  openTrainModelView,
  resetShell,
  setPredictColumn,
  trainEdaModelViaApi,
} from "./helpers";

// Test Track scenario: EDA/MLMethods/linear-regression.md
// 1. Open cars.csv from Demo files.
// 2. Top Menu > ML > Models > Train Model...
// 3. Set Predict = price, Features = all except Price and Model, Model Engine = Eda: Linear Regression.
// 4. Train the model; model result is interactive.
//
// UI scope and fallbacks (per linear-regression-run.md retrospective):
//   * File open, Train Model view open, and Predict input run through the UI.
//   * The "Select Columns" sub-dialog renders its checkboxes on a canvas — individual rows
//     are not DOM-toggleable. So the Features field is set indirectly: the underlying
//     `eda:trainLinearRegression` function is invoked with the documented column set, which
//     is what the UI's Train action would have produced. This is the "API as last resort"
//     fallback called out in the task contract.
//   * Model Engine dropdown is not always rendered in the Train Model UI on dev — the
//     engine is selected implicitly by the function name on the API call.

test.describe.serial("EDA / MLMethods / Linear Regression", () => {
  test.afterEach(async ({ page }) => {
    await resetShell(page);
  });

  test("Train Linear Regression on cars.csv predicting price", async ({
    page,
  }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, "cars.csv");
    await openTrainModelView(page);
    await setPredictColumn(page, "price");

    const trained = await trainEdaModelViaApi(
      page,
      "eda:trainLinearRegression",
      "price",
    );
    expect(trained.ok, trained.error ?? "training returned null").toBe(true);
  });
});
