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
import {columnComplete, distinctValues, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, onlyOfAnySelected, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColor, areaNotColor, boundTable, clickArea, clickAreaHolding, noBalloons, noErrors, noHighlight, readingIs, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Multivariate analysis", () => {
  const session = feature(test, "features/analyze/multivariate-analysis.feature", import.meta.url);
  test("Multivariate analysis", {tag: ["@journey", "@eda", "@realizes:ml.menu.analyze.multivariate-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await run.scenario("Running the analysis adds the scores and the prediction, and docks its charts", async () => {
      await session.step(21, "When user picks \"ML > Analyze > Multivariate Analysis...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Multivariate Analysis..."));
      await session.step(22, "Then \"Multivariate Analysis (PLS)\" dialog should be visible", () => shouldBe(page, el("\"Multivariate Analysis (PLS)\" dialog"), "visible"));
      await session.step(23, "And editor of Predict input in \"Multivariate Analysis (PLS)\" dialog should have text \"price\"", () => shouldHaveText(page, el("editor of Predict input in \"Multivariate Analysis (PLS)\" dialog"), "price"));
      await session.step(24, "And editor of Using input in \"Multivariate Analysis (PLS)\" dialog should contain text \"(15)\"", () => shouldContainText(page, el("editor of Using input in \"Multivariate Analysis (PLS)\" dialog"), "(15)"));
      await session.step(25, "And Components input in \"Multivariate Analysis (PLS)\" dialog should have value \"3\"", () => shouldHaveValue(page, el("Components input in \"Multivariate Analysis (PLS)\" dialog"), "3"));
      await session.step(26, "When user clicks on RUN button in \"Multivariate Analysis (PLS)\" dialog", () => clickOn(page, el("RUN button in \"Multivariate Analysis (PLS)\" dialog")));
      await session.step(27, "Then \"Multivariate Analysis (PLS)\" dialog should be hidden", () => shouldBe(page, el("\"Multivariate Analysis (PLS)\" dialog"), "hidden"));
      await session.step(28, "And 7 new columns should have been added", () => newColumnsCount(page, 7));
      await session.step(29, "And the table should have a column \"price (predicted)\"", () => hasColumn(page, "price (predicted)"));
      await session.step(30, "And the table should have a column \"x.score.t3\"", () => hasColumn(page, "x.score.t3"));
      await session.step(31, "And \"price (predicted)\" column should have no missing values", () => columnComplete(page, "price (predicted)"));
      await session.step(32, "And \"price (predicted)\" column should have at least 2 distinct values", () => distinctValues(page, "price (predicted)", 2));
      await session.step(33, "And table \"cars\" should have 30 rows", () => tableRows(page, "cars", 30));
      await session.step(34, "And table \"cars(Features Analysis)\" should have 15 rows", () => tableRows(page, "cars(Features Analysis)", 15));
      await session.step(35, "And table \"cars(Explained Variance)\" should have 3 rows", () => tableRows(page, "cars(Explained Variance)", 3));
      await session.step(36, "And the open tableview should have 1 grid viewer", () => viewerCount(page, 1, "grid"));
      await session.step(37, "And the open tableview should have 3 scatter plot viewers", () => viewerCount(page, 3, "scatter plot"));
      await session.step(38, "And the open tableview should have 3 bar chart viewers", () => viewerCount(page, 3, "bar chart"));
      await session.step(39, "And \"Regression Coefficients\" tab should be visible", () => shouldBe(page, el("\"Regression Coefficients\" tab"), "visible"));
      await session.step(40, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A selection in the grid reaches Observed vs. Predicted and Scores", async () => {
      await session.step(44, "Then title of second scatter plot viewer should have text \"Observed vs. Predicted\"", () => shouldHaveText(page, el("title of second scatter plot viewer"), "Observed vs. Predicted"));
      await session.step(45, "And title of third scatter plot viewer should have text \"Scores\"", () => shouldHaveText(page, el("title of third scatter plot viewer"), "Scores"));
      await session.step(46, "And grid viewer should be bound to table \"cars\"", () => boundTable(page, el("grid viewer"), "cars"));
      await session.step(47, "And second scatter plot viewer should be bound to table \"cars\"", () => boundTable(page, el("second scatter plot viewer"), "cars"));
      await session.step(48, "And third scatter plot viewer should be bound to table \"cars\"", () => boundTable(page, el("third scatter plot viewer"), "cars"));
      await session.step(49, "When user clears the row selection", () => clearSelection(page));
      await session.step(50, "Then second scatter plot viewer should show no selection highlight", () => noHighlight(page, el("second scatter plot viewer")));
      await session.step(51, "And the \"rows selected\" reading of third scatter plot viewer should be 0", () => readingIs(page, "rows selected", el("third scatter plot viewer"), 0));
      await session.step(52, "And the \"marker of row 1\" area of third scatter plot viewer should not contain the color \"#FF8C00\"", () => areaNotColor(page, "marker of row 1", el("third scatter plot viewer"), "#FF8C00"));
      await session.step(53, "And the \"marker of row 2\" area of third scatter plot viewer should not contain the color \"#FF8C00\"", () => areaNotColor(page, "marker of row 2", el("third scatter plot viewer"), "#FF8C00"));
      await session.step(54, "And the \"marker of row 3\" area of third scatter plot viewer should not contain the color \"#FF8C00\"", () => areaNotColor(page, "marker of row 3", el("third scatter plot viewer"), "#FF8C00"));
      await session.step(55, "When user clicks on the \"row header 1\" area of grid viewer", () => clickArea(page, "row header 1", el("grid viewer")));
      await session.step(56, "And user clicks on the \"row header 3\" area of grid viewer holding Shift", () => clickAreaHolding(page, "row header 3", el("grid viewer"), "Shift"));
      await session.step(57, "Then only rows where \"model\" is one of \"alfaromeo, audi, bmw\" should be selected", () => onlyOfAnySelected(page, "model", "alfaromeo, audi, bmw"));
      await session.step(58, "And the \"rows selected\" reading of second scatter plot viewer should be 3", () => readingIs(page, "rows selected", el("second scatter plot viewer"), 3));
      await session.step(59, "And the \"rows selected\" reading of third scatter plot viewer should be 3", () => readingIs(page, "rows selected", el("third scatter plot viewer"), 3));
      await session.step(60, "And the \"marker of row 1\" area of second scatter plot viewer should contain the color \"#FF8C00\"", () => areaColor(page, "marker of row 1", el("second scatter plot viewer"), "#FF8C00"));
      await session.step(61, "And the \"marker of row 2\" area of second scatter plot viewer should contain the color \"#FF8C00\"", () => areaColor(page, "marker of row 2", el("second scatter plot viewer"), "#FF8C00"));
      await session.step(62, "And the \"marker of row 3\" area of second scatter plot viewer should contain the color \"#FF8C00\"", () => areaColor(page, "marker of row 3", el("second scatter plot viewer"), "#FF8C00"));
      await session.step(63, "And the \"marker of row 1\" area of third scatter plot viewer should contain the color \"#FF8C00\"", () => areaColor(page, "marker of row 1", el("third scatter plot viewer"), "#FF8C00"));
      await session.step(64, "And the \"marker of row 2\" area of third scatter plot viewer should contain the color \"#FF8C00\"", () => areaColor(page, "marker of row 2", el("third scatter plot viewer"), "#FF8C00"));
      await session.step(65, "And the \"marker of row 3\" area of third scatter plot viewer should contain the color \"#FF8C00\"", () => areaColor(page, "marker of row 3", el("third scatter plot viewer"), "#FF8C00"));
      await session.step(66, "When user clears the row selection", () => clearSelection(page));
      await session.step(67, "Then the \"rows selected\" reading of second scatter plot viewer should be 0", () => readingIs(page, "rows selected", el("second scatter plot viewer"), 0));
      await session.step(68, "And the \"rows selected\" reading of third scatter plot viewer should be 0", () => readingIs(page, "rows selected", el("third scatter plot viewer"), 0));
      await session.step(69, "And second scatter plot viewer should show no selection highlight", () => noHighlight(page, el("second scatter plot viewer")));
      await session.step(70, "And the \"marker of row 1\" area of third scatter plot viewer should not contain the color \"#FF8C00\"", () => areaNotColor(page, "marker of row 1", el("third scatter plot viewer"), "#FF8C00"));
      await session.step(71, "And the \"marker of row 2\" area of third scatter plot viewer should not contain the color \"#FF8C00\"", () => areaNotColor(page, "marker of row 2", el("third scatter plot viewer"), "#FF8C00"));
      await session.step(72, "And the \"marker of row 3\" area of third scatter plot viewer should not contain the color \"#FF8C00\"", () => areaNotColor(page, "marker of row 3", el("third scatter plot viewer"), "#FF8C00"));
      await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A bar of Regression Coefficients selects its predictor in Loadings", async () => {
      await session.step(77, "Then title of first scatter plot viewer should have text \"Loadings\"", () => shouldHaveText(page, el("title of first scatter plot viewer"), "Loadings"));
      await session.step(78, "And first scatter plot viewer should be bound to table \"cars(Features Analysis)\"", () => boundTable(page, el("first scatter plot viewer"), "cars(Features Analysis)"));
      await session.step(79, "And the \"rows selected\" reading of first scatter plot viewer should be 0", () => readingIs(page, "rows selected", el("first scatter plot viewer"), 0));
      await session.step(80, "When user clicks on \"Regression Coefficients\" tab", () => clickOn(page, el("\"Regression Coefficients\" tab")));
      await session.step(81, "Then title of first bar chart viewer should have text \"Regression Coefficients\"", () => shouldHaveText(page, el("title of first bar chart viewer"), "Regression Coefficients"));
      await session.step(82, "And first bar chart viewer should be bound to table \"cars(Features Analysis)\"", () => boundTable(page, el("first bar chart viewer"), "cars(Features Analysis)"));
      await session.step(83, "And the \"bars\" reading of first bar chart viewer should be 15", () => readingIs(page, "bars", el("first bar chart viewer"), 15));
      await session.step(84, "When user clicks on the \"bar horsepower\" area of first bar chart viewer", () => clickArea(page, "bar horsepower", el("first bar chart viewer")));
      await session.step(85, "Then the \"rows selected\" reading of first scatter plot viewer should be 1", () => readingIs(page, "rows selected", el("first scatter plot viewer"), 1));
      await session.step(86, "And the \"marker of row 11\" area of first scatter plot viewer should contain the color \"#FF8C00\"", () => areaColor(page, "marker of row 11", el("first scatter plot viewer"), "#FF8C00"));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
