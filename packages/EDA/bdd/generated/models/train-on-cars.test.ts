/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/models/train-on-cars.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.models.train-model, eda.model.linear-regression, eda.model.pls-regression, eda.model.xg-boost]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragSliderTo, enterInto, hoverOver, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValueBetween, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, dragSelectionOverArea, noBalloons, noErrors, readingAtLeast, readingIs, readingReads, someHighlight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Training a model to predict the price of a car", () => {
  const session = feature(test, "features/models/train-on-cars.feature", import.meta.url);
  test("Training a model to predict the price of a car", {tag: ["@journey", "@eda", "@realizes:ml.menu.models.train-model", "@realizes:eda.model.linear-regression", "@realizes:eda.model.pls-regression", "@realizes:eda.model.xg-boost"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await run.scenario("The view takes price as the target and fifteen columns as the features", async () => {
      await session.step(27, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(28, "Then the \"Predictive model\" view should be current", () => viewIsCurrent(page, "Predictive model"));
      await session.step(29, "When user selects \"price\" in Predict input", () => selectIn(page, "price", el("Predict input")));
      await session.step(30, "Then editor of Predict input should have text \"price\"", () => shouldHaveText(page, el("editor of Predict input"), "price"));
      await session.step(31, "When user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(32, "And user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(33, "Then \"Select columns...\" dialog should contain text \"17 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "17 checked"));
      await session.step(34, "When user types \"price\" into Search input in \"Select columns...\" dialog", () => typeInto(page, "price", el("Search input in \"Select columns...\" dialog")));
      await session.step(35, "Then the \"text of cell 17 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"price\"", () => readingReads(page, "text of cell 17 of __name", el("grid viewer in \"Select columns...\" dialog"), "price"));
      await session.step(36, "When user clicks on the \"cell 17 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 17 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(37, "And user types \"model\" into Search input in \"Select columns...\" dialog", () => typeInto(page, "model", el("Search input in \"Select columns...\" dialog")));
      await session.step(38, "Then the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"model\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "model"));
      await session.step(39, "When user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(40, "Then \"Select columns...\" dialog should contain text \"15 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "15 checked"));
      await session.step(41, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(42, "Then editor of Features input should contain text \"(15)\"", () => shouldContainText(page, el("editor of Features input"), "(15)"));
      await session.step(43, "And \"Model Engine\" input should be visible", () => shouldBe(page, el("\"Model Engine\" input"), "visible"));
    });
    await run.scenario("Linear Regression trains on the fifteen features", async () => {
      await session.step(46, "When user selects \"Eda: Linear Regression\" in \"Model Engine\" input", () => selectIn(page, "Eda: Linear Regression", el("\"Model Engine\" input")));
      await session.step(47, "Then \"Eda: Linear Regression\" heading should be visible", () => shouldBe(page, el("\"Eda: Linear Regression\" heading"), "visible"));
      await session.step(48, "And \"R squared\" table row should be visible", () => shouldBe(page, el("\"R squared\" table row"), "visible"));
      await session.step(49, "And \"Predicted price vs Actual\" label should be visible", () => shouldBe(page, el("\"Predicted price vs Actual\" label"), "visible"));
      await session.step(50, "And the \"axes\" reading of pc plot viewer should be 17", () => readingIs(page, "axes", el("pc plot viewer"), 17));
      await session.step(51, "And the \"rows selected\" reading of first scatter plot viewer should be 0", () => readingIs(page, "rows selected", el("first scatter plot viewer"), 0));
      await session.step(52, "When user drags a selection box over the \"view\" area of first scatter plot viewer", () => dragSelectionOverArea(page, "view", el("first scatter plot viewer")));
      await session.step(53, "Then the \"rows selected\" reading of first scatter plot viewer should be at least 1", () => readingAtLeast(page, "rows selected", el("first scatter plot viewer"), 1));
      await session.step(54, "And first scatter plot viewer should show a selection highlight", () => someHighlight(page, el("first scatter plot viewer")));
      await session.step(55, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("PLS Regression trains with its components and shows its own charts", async () => {
      await session.step(59, "When user selects \"Eda: PLS Regression\" in \"Model Engine\" input", () => selectIn(page, "Eda: PLS Regression", el("\"Model Engine\" input")));
      await session.step(60, "And user enters \"3\" into Components input", () => enterInto(page, "3", el("Components input")));
      await session.step(61, "Then \"Eda: PLS Regression\" heading should be visible", () => shouldBe(page, el("\"Eda: PLS Regression\" heading"), "visible"));
      await session.step(62, "And \"components\" table row should contain text \"3\"", () => shouldContainText(page, el("\"components\" table row"), "3"));
      await session.step(63, "And \"R squared\" table row should be visible", () => shouldBe(page, el("\"R squared\" table row"), "visible"));
      await session.step(64, "And the \"axes\" reading of pc plot viewer should be 17", () => readingIs(page, "axes", el("pc plot viewer"), 17));
      await session.step(65, "And the \"bars\" reading of first bar chart viewer should be 15", () => readingIs(page, "bars", el("first bar chart viewer"), 15));
      await session.step(66, "And the \"bars\" reading of third bar chart viewer should be 3", () => readingIs(page, "bars", el("third bar chart viewer"), 3));
      await session.step(67, "And the \"rows shown\" reading of third scatter plot viewer should be 15", () => readingIs(page, "rows shown", el("third scatter plot viewer"), 15));
      await session.step(68, "And the \"rows shown\" reading of fourth scatter plot viewer should be 30", () => readingIs(page, "rows shown", el("fourth scatter plot viewer"), 30));
      await session.step(69, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("XGBoost retrains as its clickers and sliders move", async () => {
      await session.step(73, "When user selects \"Eda: XGBoost\" in \"Model Engine\" input", () => selectIn(page, "Eda: XGBoost", el("\"Model Engine\" input")));
      await session.step(74, "And user enters \"20\" into Iterations input", () => enterInto(page, "20", el("Iterations input")));
      await session.step(75, "And user enters \"6\" into \"Max Depth\" input", () => enterInto(page, "6", el("\"Max Depth\" input")));
      await session.step(76, "Then \"Eda: XGBoost\" heading should be visible", () => shouldBe(page, el("\"Eda: XGBoost\" heading"), "visible"));
      await session.step(77, "And \"iterations\" table row should contain text \"20\"", () => shouldContainText(page, el("\"iterations\" table row"), "20"));
      await session.step(78, "And \"maxDepth\" table row should contain text \"6\"", () => shouldContainText(page, el("\"maxDepth\" table row"), "6"));
      await session.step(79, "When user hovers over Iterations input", () => hoverOver(page, el("Iterations input")));
      await session.step(80, "And user clicks on plus icon in Iterations input", () => clickOn(page, el("plus icon in Iterations input")));
      await session.step(81, "Then \"iterations\" table row should contain text \"21\"", () => shouldContainText(page, el("\"iterations\" table row"), "21"));
      await session.step(82, "When user hovers over \"Max Depth\" input", () => hoverOver(page, el("\"Max Depth\" input")));
      await session.step(83, "And user clicks on minus icon in \"Max Depth\" input", () => clickOn(page, el("minus icon in \"Max Depth\" input")));
      await session.step(84, "Then \"maxDepth\" table row should contain text \"5\"", () => shouldContainText(page, el("\"maxDepth\" table row"), "5"));
      await session.step(85, "When user enters \"0.3\" into Rate input", () => enterInto(page, "0.3", el("Rate input")));
      await session.step(86, "And user enters \"1\" into Lambda input", () => enterInto(page, "1", el("Lambda input")));
      await session.step(87, "And user enters \"0\" into Alpha input", () => enterInto(page, "0", el("Alpha input")));
      await session.step(88, "Then \"eta\" table row should contain text \"0.30\"", () => shouldContainText(page, el("\"eta\" table row"), "0.30"));
      await session.step(89, "When user drags the slider of Rate input to 0.5", () => dragSliderTo(page, el("Rate input"), 0.5));
      await session.step(90, "Then Rate input should have a value between 0.48 and 0.52", () => shouldHaveValueBetween(page, el("Rate input"), 0.48, 0.52));
      await session.step(91, "And \"eta\" table row should not contain text \"eta0.30\"", () => shouldNotContainText(page, el("\"eta\" table row"), "eta0.30"));
      await session.step(92, "When user drags the slider of Lambda input to 50", () => dragSliderTo(page, el("Lambda input"), 50));
      await session.step(93, "Then Lambda input should have a value between 48 and 52", () => shouldHaveValueBetween(page, el("Lambda input"), 48, 52));
      await session.step(94, "And \"lambda\" table row should not contain text \"lambda1\"", () => shouldNotContainText(page, el("\"lambda\" table row"), "lambda1"));
      await session.step(95, "When user drags the slider of Alpha input to 40", () => dragSliderTo(page, el("Alpha input"), 40));
      await session.step(96, "Then Alpha input should have a value between 38 and 42", () => shouldHaveValueBetween(page, el("Alpha input"), 38, 42));
      await session.step(97, "And \"alpha\" table row should not contain text \"alpha0\"", () => shouldNotContainText(page, el("\"alpha\" table row"), "alpha0"));
      await session.step(98, "And \"R squared\" table row should be visible", () => shouldBe(page, el("\"R squared\" table row"), "visible"));
      await session.step(99, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
