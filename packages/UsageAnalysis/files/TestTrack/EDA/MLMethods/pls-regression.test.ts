import { test, expect } from "./helpers";
import {
  openDemoCsv,
  openTrainModelView,
  resetShell,
  setPredictColumn,
  trainEdaModelViaApi,
} from "./helpers";

test.describe.serial("EDA / MLMethods / PLS Regression", () => {
  test.afterEach(async ({ page }) => {
    await resetShell(page);
  });

  test("Train PLS Regression on cars.csv with 3 components predicting price", async ({
    page,
  }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, "cars.csv");
    await openTrainModelView(page);
    await setPredictColumn(page, "price");

    const trained = await trainEdaModelViaApi(
      page,
      "eda:trainPLSRegression",
      "price",
      {
        numericOnly: true,
        extraParams: { components: 3 },
      },
    );
    expect(trained.ok, trained.error ?? "training returned null").toBe(true);
  });
});
