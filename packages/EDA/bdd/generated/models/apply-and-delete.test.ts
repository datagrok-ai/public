/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/models/apply-and-delete.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.models.train-model, ml.menu.models.apply-model, eda.model.pls-regression, eda.model.linear-regression]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, pressKeyIn, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, modelsOnServer, noModelOnServer, openDataset, openTableOf, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, pickFromContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Saved models applied to new data and deleted from the gallery", () => {
  const session = feature(test, "features/models/apply-and-delete.feature", import.meta.url);
  test("Saved models applied to new data and deleted from the gallery", {tag: ["@journey", "@eda", "@realizes:ml.menu.models.train-model", "@realizes:ml.menu.models.apply-model", "@realizes:eda.model.pls-regression", "@realizes:eda.model.linear-regression"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And no predictive model named \"BDD-Iris-PLS-{run}\" is on the server", () => noModelOnServer(page, session.text("BDD-Iris-PLS-{run}")));
    await session.step(28, "And no predictive model named \"BDD-Iris-LR-{run}\" is on the server", () => noModelOnServer(page, session.text("BDD-Iris-LR-{run}")));
    await run.scenario("A PLS model predicts Petal.Width from three other columns, and is saved", async () => {
      await session.step(31, "Given user opens iris dataset", () => openDataset(page, ds("iris")));
      await session.step(32, "Then the table should have 6 columns", () => columnCount(page, 6));
      await session.step(33, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(34, "Then the \"Predictive model\" view should be current", () => viewIsCurrent(page, "Predictive model"));
      await session.step(35, "When user selects \"Petal.Width\" in Predict input", () => selectIn(page, "Petal.Width", el("Predict input")));
      await session.step(36, "And user clicks on editor of Features input", () => clickOn(page, el("editor of Features input")));
      await session.step(37, "And user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
      await session.step(38, "Then the \"text of cell 2 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"Sepal.Length\"", () => readingReads(page, "text of cell 2 of __name", el("grid viewer in \"Select columns...\" dialog"), "Sepal.Length"));
      await session.step(39, "And the \"text of cell 4 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"Petal.Length\"", () => readingReads(page, "text of cell 4 of __name", el("grid viewer in \"Select columns...\" dialog"), "Petal.Length"));
      await session.step(40, "When user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(41, "And user clicks on the \"cell 3 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 3 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(42, "And user clicks on the \"cell 4 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 4 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(43, "Then \"Select columns...\" dialog should contain text \"3 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "3 checked"));
      await session.step(44, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(45, "Then editor of Features input should contain text \"(3)\"", () => shouldContainText(page, el("editor of Features input"), "(3)"));
      await session.step(46, "When user selects \"Eda: PLS Regression\" in \"Model Engine\" input", () => selectIn(page, "Eda: PLS Regression", el("\"Model Engine\" input")));
      await session.step(47, "Then Components input should have value \"3\"", () => shouldHaveValue(page, el("Components input"), "3"));
      await session.step(48, "And model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(49, "And \"Eda: PLS Regression\" heading should be visible", () => shouldBe(page, el("\"Eda: PLS Regression\" heading"), "visible"));
      await session.step(50, "And \"R squared\" table row should be visible", () => shouldBe(page, el("\"R squared\" table row"), "visible"));
      await session.step(51, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(52, "And user enters \"BDD-Iris-PLS-{run}\" into Name input in dialog", () => enterInto(page, session.text("BDD-Iris-PLS-{run}"), el("Name input in dialog")));
      await session.step(53, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(54, "Then dialog should be absent", () => shouldBe(page, el("dialog"), "absent"));
      await session.step(55, "And 1 predictive model named \"BDD-Iris-PLS-{run}\" should be on the server", () => modelsOnServer(page, 1, session.text("BDD-Iris-PLS-{run}")));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A linear regression on the same inputs is saved too", async () => {
      await session.step(59, "When user selects \"Eda: Linear Regression\" in \"Model Engine\" input", () => selectIn(page, "Eda: Linear Regression", el("\"Model Engine\" input")));
      await session.step(60, "Then model preview should be ready", () => shouldBe(page, el("model preview"), "ready"));
      await session.step(61, "And \"Eda: Linear Regression\" heading should be visible", () => shouldBe(page, el("\"Eda: Linear Regression\" heading"), "visible"));
      await session.step(62, "And \"R squared\" table row should be visible", () => shouldBe(page, el("\"R squared\" table row"), "visible"));
      await session.step(63, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(64, "And user enters \"BDD-Iris-LR-{run}\" into Name input in dialog", () => enterInto(page, session.text("BDD-Iris-LR-{run}"), el("Name input in dialog")));
      await session.step(65, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(66, "Then dialog should be absent", () => shouldBe(page, el("dialog"), "absent"));
      await session.step(67, "And 1 predictive model named \"BDD-Iris-LR-{run}\" should be on the server", () => modelsOnServer(page, 1, session.text("BDD-Iris-LR-{run}")));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Apply Model offers the newest model first, the earlier one below it, and adds the prediction", async () => {
      await session.step(71, "Given user opens iris dataset", () => openDataset(page, ds("iris")));
      await session.step(72, "And the context panel is open", () => contextPanelOpen(page));
      await session.step(73, "When user picks \"ML > Models > Apply Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Apply Model..."));
      await session.step(74, "Then \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("\"Apply predictive model\" dialog"), "visible"));
      await session.step(75, "When user presses ArrowDown in Model input in \"Apply predictive model\" dialog", () => pressKeyIn(page, "ArrowDown", el("Model input in \"Apply predictive model\" dialog")));
      await session.step(76, "Then the context panel should show \"BDD-Iris-PLS-{run}\"", () => contextPanelShows(page, session.text("BDD-Iris-PLS-{run}")));
      await session.step(77, "When user presses ArrowUp in Model input in \"Apply predictive model\" dialog", () => pressKeyIn(page, "ArrowUp", el("Model input in \"Apply predictive model\" dialog")));
      await session.step(78, "Then the context panel should show \"BDD-Iris-LR-{run}\"", () => contextPanelShows(page, session.text("BDD-Iris-LR-{run}")));
      await session.step(79, "When user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(80, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(81, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(82, "And the table should have 7 columns", () => columnCount(page, 7));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The model saved before it is one ArrowDown away and adds its prediction beside it", async () => {
      await session.step(86, "When user picks \"ML > Models > Apply Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Apply Model..."));
      await session.step(87, "Then \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("\"Apply predictive model\" dialog"), "visible"));
      await session.step(88, "When user presses ArrowDown in Model input in \"Apply predictive model\" dialog", () => pressKeyIn(page, "ArrowDown", el("Model input in \"Apply predictive model\" dialog")));
      await session.step(89, "Then the context panel should show \"BDD-Iris-PLS-{run}\"", () => contextPanelShows(page, session.text("BDD-Iris-PLS-{run}")));
      await session.step(90, "When user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(91, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(92, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(93, "And the table should have 8 columns", () => columnCount(page, 8));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A model applies to another table that has its input columns", async () => {
      await session.step(97, "Given user opens a table \"new readings\" with:", () => openTableOf(page, "new readings", [["Sepal.Length","Sepal.Width","Petal.Length"],["5.1","3.5","1.4"],["6.4","3.2","4.5"],["6.3","3.3","6.0"]]), [["Sepal.Length","Sepal.Width","Petal.Length"],["5.1","3.5","1.4"],["6.4","3.2","4.5"],["6.3","3.3","6.0"]]);
      await session.step(102, "When user picks \"ML > Models > Apply Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Apply Model..."));
      await session.step(103, "Then \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("\"Apply predictive model\" dialog"), "visible"));
      await session.step(104, "When user presses ArrowDown in Model input in \"Apply predictive model\" dialog", () => pressKeyIn(page, "ArrowDown", el("Model input in \"Apply predictive model\" dialog")));
      await session.step(105, "And user presses ArrowUp in Model input in \"Apply predictive model\" dialog", () => pressKeyIn(page, "ArrowUp", el("Model input in \"Apply predictive model\" dialog")));
      await session.step(106, "Then the context panel should show \"BDD-Iris-LR-{run}\"", () => contextPanelShows(page, session.text("BDD-Iris-LR-{run}")));
      await session.step(107, "When user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(108, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(109, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(110, "And the table should have 4 columns", () => columnCount(page, 4));
      await session.step(111, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The gallery shows the model's details and performance, and deletes it", async () => {
      await session.step(114, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(115, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(116, "When user clicks on \"Predictive models\" tree node inside browse tree", () => clickOn(page, el("\"Predictive models\" tree node inside browse tree")));
      await session.step(117, "Then the \"Models\" view should be current", () => viewIsCurrent(page, "Models"));
      await session.step(118, "When user clicks on \"BDD-Iris-LR-{run}\" label in gallery", () => clickOn(page, el(session.text("\"BDD-Iris-LR-{run}\" label in gallery"))));
      await session.step(119, "Then the context panel should show \"BDD-Iris-LR-{run}\"", () => contextPanelShows(page, session.text("BDD-Iris-LR-{run}")));
      await session.step(120, "And \"Details\" pane in context panel should be visible", () => shouldBe(page, el("\"Details\" pane in context panel"), "visible"));
      await session.step(121, "And \"Performance\" pane in context panel should be visible", () => shouldBe(page, el("\"Performance\" pane in context panel"), "visible"));
      await session.step(122, "And \"Sharing\" pane in context panel should be visible", () => shouldBe(page, el("\"Sharing\" pane in context panel"), "visible"));
      await session.step(123, "When user picks \"Delete\" from the context menu of \"BDD-Iris-LR-{run}\" label in gallery", () => pickFromContextMenu(page, "Delete", el(session.text("\"BDD-Iris-LR-{run}\" label in gallery"))));
      await session.step(124, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(125, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(126, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(127, "And \"BDD-Iris-LR-{run}\" label in gallery should be absent", () => shouldBe(page, el(session.text("\"BDD-Iris-LR-{run}\" label in gallery")), "absent"));
      await session.step(128, "And 0 predictive models named \"BDD-Iris-LR-{run}\" should be on the server", () => modelsOnServer(page, 0, session.text("BDD-Iris-LR-{run}")));
      await session.step(129, "When user picks \"Delete\" from the context menu of \"BDD-Iris-PLS-{run}\" label in gallery", () => pickFromContextMenu(page, "Delete", el(session.text("\"BDD-Iris-PLS-{run}\" label in gallery"))));
      await session.step(130, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(131, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(132, "And 0 predictive models named \"BDD-Iris-PLS-{run}\" should be on the server", () => modelsOnServer(page, 0, session.text("BDD-Iris-PLS-{run}")));
      await session.step(133, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
