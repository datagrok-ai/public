/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/models/train-on-iris.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.models.train-model, eda.model.softmax, eda.model.xg-boost]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragSliderTo, enterInto, hoverOver, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValueBetween, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Training a model to classify iris species", () => {
  const session = feature(test, "features/models/train-on-iris.feature", import.meta.url);
  test("Training a model to classify iris species", {tag: ["@journey", "@eda", "@realizes:ml.menu.models.train-model", "@realizes:eda.model.softmax", "@realizes:eda.model.xg-boost"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens iris dataset", () => openDataset(page, ds("iris")));
    await run.scenario("The view takes Species as the target and the four measurements as the features", async () => {
      await session.step(20, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(21, "Then the \"Predictive model\" view should be current", () => viewIsCurrent(page, "Predictive model"));
      await session.step(22, "When user selects \"Species\" in Predict input", () => selectIn(page, "Species", el("Predict input")));
      await session.step(23, "Then editor of Predict input should have text \"Species\"", () => shouldHaveText(page, el("editor of Predict input"), "Species"));
      await session.step(24, "When user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(25, "And user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(26, "Then \"Select columns...\" dialog should contain text \"6 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "6 checked"));
      await session.step(27, "And the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"col 1\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "col 1"));
      await session.step(28, "And the \"text of cell 6 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"Species\"", () => readingReads(page, "text of cell 6 of __name", el("grid viewer in \"Select columns...\" dialog"), "Species"));
      await session.step(29, "When user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(30, "And user clicks on the \"cell 6 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 6 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(31, "Then \"Select columns...\" dialog should contain text \"4 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "4 checked"));
      await session.step(32, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(33, "Then editor of Features input should contain text \"(4)\"", () => shouldContainText(page, el("editor of Features input"), "(4)"));
      await session.step(34, "And \"One-hot encoding\" input should be hidden", () => shouldBe(page, el("\"One-hot encoding\" input"), "hidden"));
      await session.step(35, "And \"Model Engine\" input should be visible", () => shouldBe(page, el("\"Model Engine\" input"), "visible"));
    });
    await run.scenario("Softmax classifies the species and retrains as its hyperparameters move", async () => {
      await session.step(38, "When user selects \"Eda: Softmax\" in \"Model Engine\" input", () => selectIn(page, "Eda: Softmax", el("\"Model Engine\" input")));
      await session.step(39, "And user enters \"100\" into Iterations input", () => enterInto(page, "100", el("Iterations input")));
      await session.step(40, "And user enters \"2\" into Rate input", () => enterInto(page, "2", el("Rate input")));
      await session.step(41, "And user enters \"0.1\" into Penalty input", () => enterInto(page, "0.1", el("Penalty input")));
      await session.step(42, "Then \"Eda: Softmax\" heading should be visible", () => shouldBe(page, el("\"Eda: Softmax\" heading"), "visible"));
      await session.step(43, "And \"Accuracy\" table row should be visible", () => shouldBe(page, el("\"Accuracy\" table row"), "visible"));
      await session.step(44, "And \"iterations\" table row should contain text \"100\"", () => shouldContainText(page, el("\"iterations\" table row"), "100"));
      await session.step(45, "And \"rate\" table row should contain text \"rate2\"", () => shouldContainText(page, el("\"rate\" table row"), "rate2"));
      await session.step(46, "And \"penalty\" table row should contain text \"penalty0.10\"", () => shouldContainText(page, el("\"penalty\" table row"), "penalty0.10"));
      await session.step(47, "And \"Predicted Species vs Actual\" label should be visible", () => shouldBe(page, el("\"Predicted Species vs Actual\" label"), "visible"));
      await session.step(48, "And \"Confusions\" label should be visible", () => shouldBe(page, el("\"Confusions\" label"), "visible"));
      await session.step(49, "When user drags the slider of Rate input to 10", () => dragSliderTo(page, el("Rate input"), 10));
      await session.step(50, "Then Rate input should have a value between 9.5 and 10.5", () => shouldHaveValueBetween(page, el("Rate input"), 9.5, 10.5));
      await session.step(51, "And \"rate\" table row should not contain text \"rate2\"", () => shouldNotContainText(page, el("\"rate\" table row"), "rate2"));
      await session.step(52, "When user drags the slider of Penalty input to 0.5", () => dragSliderTo(page, el("Penalty input"), 0.5));
      await session.step(53, "Then Penalty input should have a value between 0.48 and 0.52", () => shouldHaveValueBetween(page, el("Penalty input"), 0.48, 0.52));
      await session.step(54, "And \"penalty\" table row should not contain text \"penalty0.10\"", () => shouldNotContainText(page, el("\"penalty\" table row"), "penalty0.10"));
      await session.step(55, "And \"Accuracy\" table row should be visible", () => shouldBe(page, el("\"Accuracy\" table row"), "visible"));
      await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("XGBoost classifies the species and retrains as its clickers and sliders move", async () => {
      await session.step(60, "When user selects \"Eda: XGBoost\" in \"Model Engine\" input", () => selectIn(page, "Eda: XGBoost", el("\"Model Engine\" input")));
      await session.step(61, "And user enters \"20\" into Iterations input", () => enterInto(page, "20", el("Iterations input")));
      await session.step(62, "And user enters \"6\" into \"Max Depth\" input", () => enterInto(page, "6", el("\"Max Depth\" input")));
      await session.step(63, "And user enters \"0.3\" into Rate input", () => enterInto(page, "0.3", el("Rate input")));
      await session.step(64, "And user enters \"1\" into Lambda input", () => enterInto(page, "1", el("Lambda input")));
      await session.step(65, "And user enters \"0\" into Alpha input", () => enterInto(page, "0", el("Alpha input")));
      await session.step(66, "Then \"Eda: XGBoost\" heading should be visible", () => shouldBe(page, el("\"Eda: XGBoost\" heading"), "visible"));
      await session.step(67, "And \"Accuracy\" table row should be visible", () => shouldBe(page, el("\"Accuracy\" table row"), "visible"));
      await session.step(68, "And \"iterations\" table row should contain text \"20\"", () => shouldContainText(page, el("\"iterations\" table row"), "20"));
      await session.step(69, "And \"eta\" table row should contain text \"eta0.30\"", () => shouldContainText(page, el("\"eta\" table row"), "eta0.30"));
      await session.step(70, "And \"lambda\" table row should contain text \"lambda1\"", () => shouldContainText(page, el("\"lambda\" table row"), "lambda1"));
      await session.step(71, "And \"alpha\" table row should contain text \"alpha0\"", () => shouldContainText(page, el("\"alpha\" table row"), "alpha0"));
      await session.step(72, "When user hovers over Iterations input", () => hoverOver(page, el("Iterations input")));
      await session.step(73, "And user clicks on plus icon in Iterations input", () => clickOn(page, el("plus icon in Iterations input")));
      await session.step(74, "Then \"iterations\" table row should contain text \"21\"", () => shouldContainText(page, el("\"iterations\" table row"), "21"));
      await session.step(75, "When user hovers over \"Max Depth\" input", () => hoverOver(page, el("\"Max Depth\" input")));
      await session.step(76, "And user clicks on minus icon in \"Max Depth\" input", () => clickOn(page, el("minus icon in \"Max Depth\" input")));
      await session.step(77, "Then \"maxDepth\" table row should contain text \"5\"", () => shouldContainText(page, el("\"maxDepth\" table row"), "5"));
      await session.step(78, "When user drags the slider of Rate input to 0.5", () => dragSliderTo(page, el("Rate input"), 0.5));
      await session.step(79, "Then Rate input should have a value between 0.48 and 0.52", () => shouldHaveValueBetween(page, el("Rate input"), 0.48, 0.52));
      await session.step(80, "And \"eta\" table row should not contain text \"eta0.30\"", () => shouldNotContainText(page, el("\"eta\" table row"), "eta0.30"));
      await session.step(81, "When user drags the slider of Lambda input to 50", () => dragSliderTo(page, el("Lambda input"), 50));
      await session.step(82, "Then Lambda input should have a value between 48 and 52", () => shouldHaveValueBetween(page, el("Lambda input"), 48, 52));
      await session.step(83, "And \"lambda\" table row should not contain text \"lambda1\"", () => shouldNotContainText(page, el("\"lambda\" table row"), "lambda1"));
      await session.step(84, "When user drags the slider of Alpha input to 40", () => dragSliderTo(page, el("Alpha input"), 40));
      await session.step(85, "Then Alpha input should have a value between 38 and 42", () => shouldHaveValueBetween(page, el("Alpha input"), 38, 42));
      await session.step(86, "And \"alpha\" table row should not contain text \"alpha0\"", () => shouldNotContainText(page, el("\"alpha\" table row"), "alpha0"));
      await session.step(87, "And \"Accuracy\" table row should be visible", () => shouldBe(page, el("\"Accuracy\" table row"), "visible"));
      await session.step(88, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
