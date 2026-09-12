/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/multivariate-analysis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.analyze.multivariate-analysis]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, selectWhereOneOf, selectedRowCount, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, viewHoldsViewers} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, readingIs, someHighlight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Multivariate analysis", () => {
  const session = feature(test, "features/analyze/multivariate-analysis.feature", import.meta.url);
  test("Multivariate analysis", {tag: ["@journey", "@eda", "@realizes:ml.menu.analyze.multivariate-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await run.scenario("Running the analysis adds the scores and the prediction, and docks its charts", async () => {
      await session.step(19, "When user picks \"ML > Analyze > Multivariate Analysis...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Multivariate Analysis..."));
      await session.step(20, "Then \"Multivariate Analysis (PLS)\" dialog should be visible", () => shouldBe(page, el("\"Multivariate Analysis (PLS)\" dialog"), "visible"));
      await session.step(21, "And editor of Predict input in \"Multivariate Analysis (PLS)\" dialog should have text \"price\"", () => shouldHaveText(page, el("editor of Predict input in \"Multivariate Analysis (PLS)\" dialog"), "price"));
      await session.step(22, "And editor of Using input in \"Multivariate Analysis (PLS)\" dialog should contain text \"(15)\"", () => shouldContainText(page, el("editor of Using input in \"Multivariate Analysis (PLS)\" dialog"), "(15)"));
      await session.step(23, "And Components input in \"Multivariate Analysis (PLS)\" dialog should have value \"3\"", () => shouldHaveValue(page, el("Components input in \"Multivariate Analysis (PLS)\" dialog"), "3"));
      await session.step(24, "When user clicks on RUN button in \"Multivariate Analysis (PLS)\" dialog", () => clickOn(page, el("RUN button in \"Multivariate Analysis (PLS)\" dialog")));
      await session.step(25, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(26, "And \"Multivariate Analysis (PLS)\" dialog should be hidden", () => shouldBe(page, el("\"Multivariate Analysis (PLS)\" dialog"), "hidden"));
      await session.step(27, "And 7 new columns should have been added", () => newColumnsCount(page, 7));
      await session.step(28, "And the table should have a column \"price (predicted)\"", () => hasColumn(page, "price (predicted)"));
      await session.step(29, "And the table should have a column \"x.score.t3\"", () => hasColumn(page, "x.score.t3"));
      await session.step(30, "And \"price (predicted)\" column should have no missing values", () => columnComplete(page, "price (predicted)"));
      await session.step(31, "And table \"cars(Features Analysis)\" should have 15 rows", () => tableRows(page, "cars(Features Analysis)", 15));
      await session.step(32, "And table \"cars(Explained Variance)\" should have 3 rows", () => tableRows(page, "cars(Explained Variance)", 3));
      await session.step(33, "And the current view should hold at least 7 viewers", () => viewHoldsViewers(page, 7));
      await session.step(34, "And \"Regression Coefficients\" tab should be visible", () => shouldBe(page, el("\"Regression Coefficients\" tab"), "visible"));
      await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A selection in the grid reaches Observed vs. Predicted and Scores", async () => {
      await session.step(39, "Then title of second scatter plot viewer should have text \"Observed vs. Predicted\"", () => shouldHaveText(page, el("title of second scatter plot viewer"), "Observed vs. Predicted"));
      await session.step(40, "And title of third scatter plot viewer should have text \"Scores\"", () => shouldHaveText(page, el("title of third scatter plot viewer"), "Scores"));
      await session.step(41, "When user selects rows where \"model\" is one of \"porsche, jaguar, mercedes\"", () => selectWhereOneOf(page, "model", "porsche, jaguar, mercedes"));
      await session.step(42, "Then 3 rows should be selected", () => selectedRowCount(page, 3));
      await session.step(43, "And the \"rows selected\" reading of second scatter plot viewer should be 3", () => readingIs(page, "rows selected", el("second scatter plot viewer"), 3));
      await session.step(44, "And the \"rows selected\" reading of third scatter plot viewer should be 3", () => readingIs(page, "rows selected", el("third scatter plot viewer"), 3));
      await session.step(45, "And second scatter plot viewer should show a selection highlight", () => someHighlight(page, el("second scatter plot viewer")));
      await session.step(46, "And third scatter plot viewer should show a selection highlight", () => someHighlight(page, el("third scatter plot viewer")));
      await session.step(47, "When user clears the row selection", () => clearSelection(page));
      await session.step(48, "Then the \"rows selected\" reading of second scatter plot viewer should be 0", () => readingIs(page, "rows selected", el("second scatter plot viewer"), 0));
    });
    await run.scenario("A bar of Regression Coefficients selects its predictor in Loadings", async () => {
      await session.step(51, "Then title of first scatter plot viewer should have text \"Loadings\"", () => shouldHaveText(page, el("title of first scatter plot viewer"), "Loadings"));
      await session.step(52, "And the \"rows selected\" reading of first scatter plot viewer should be 0", () => readingIs(page, "rows selected", el("first scatter plot viewer"), 0));
      await session.step(53, "When user clicks on \"Regression Coefficients\" tab", () => clickOn(page, el("\"Regression Coefficients\" tab")));
      await session.step(54, "Then title of first bar chart viewer should have text \"Regression Coefficients\"", () => shouldHaveText(page, el("title of first bar chart viewer"), "Regression Coefficients"));
      await session.step(55, "And the \"bars\" reading of first bar chart viewer should be 15", () => readingIs(page, "bars", el("first bar chart viewer"), 15));
      await session.step(56, "When user clicks on the \"bar horsepower\" area of first bar chart viewer", () => clickArea(page, "bar horsepower", el("first bar chart viewer")));
      await session.step(57, "Then the \"rows selected\" reading of first scatter plot viewer should be 1", () => readingIs(page, "rows selected", el("first scatter plot viewer"), 1));
      await session.step(58, "And first scatter plot viewer should show a selection highlight", () => someHighlight(page, el("first scatter plot viewer")));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
