/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/app/app-and-menus.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {gridShowsColumns, gridShowsNoText, hoverMethodCell, methodAbove, methodCellReads, openViewsNamed, topMenuGroupExactly} from '../../bindings/catalog.js';
import {noBalloonAtAll, openDemoFile} from '../../bindings/train.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, hoverOver, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu, topMenuLists} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, closeCurrentView, dialogCloses, switchTableView, switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Forge app, its menus, and the old tool beside it", () => {
  const session = feature(test, "features/app/app-and-menus.feature", import.meta.url);
  test("The Forge app, its menus, and the old tool beside it", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And user opens the demo file \"iris.csv\" as \"forge-test-iris-{time}\"", () => openDemoFile(page, "iris.csv", session.text("forge-test-iris-{time}")));
    await session.step(34, "Then no balloon of any colour should have been shown", () => noBalloonAtAll(page));
    await session.step(35, "And no errors should have been logged", () => noErrors(page));
    await run.scenario("Forge opens from Apps and from the ML menu; the old tool keeps its menu and its Browse node (1.1, 1.2, 6.4, 4.13)", async () => {
      await session.step(38, "Then the top menu should list:", () => topMenuLists(page, [["ML > Forge > Models"],["ML > Forge > Train..."],["ML > Forge > Apply..."],["ML > Models > Train Model..."],["ML > Models > Apply Model..."]]), [["ML > Forge > Models"],["ML > Forge > Train..."],["ML > Forge > Apply..."],["ML > Models > Train Model..."],["ML > Models > Apply Model..."]]);
      await session.step(44, "And the \"ML > Forge\" top menu group should list exactly \"Models, Train..., Apply...\"", () => topMenuGroupExactly(page, "ML > Forge", "Models, Train..., Apply..."));
      await session.step(45, "And the \"ML > Models\" top menu group should list exactly \"Train Model..., Apply Model...\"", () => topMenuGroupExactly(page, "ML > Models", "Train Model..., Apply Model..."));
      await session.step(46, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(47, "When user expands \"Apps\" tree node inside browse tree", () => expand(page, el("\"Apps\" tree node inside browse tree")));
      await session.step(48, "And user clicks on \"Apps > Forge\" tree node inside browse tree", () => clickOn(page, el("\"Apps > Forge\" tree node inside browse tree")));
      await session.step(49, "Then the \"Forge\" view should be current", () => viewIsCurrent(page, "Forge"));
      await session.step(50, "And \"Methods\" heading in forge app should be visible", () => shouldBe(page, el("\"Methods\" heading in forge app"), "visible"));
      await session.step(51, "And \"Models\" heading in forge app should be visible", () => shouldBe(page, el("\"Models\" heading in forge app"), "visible"));
      await session.step(52, "And first grid viewer should be visible", () => shouldBe(page, el("first grid viewer"), "visible"));
      await session.step(53, "And second grid viewer should be visible", () => shouldBe(page, el("second grid viewer"), "visible"));
      await session.step(54, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
      await session.step(55, "And user clicks on \"Predictive models\" tree node inside browse tree", () => clickOn(page, el("\"Predictive models\" tree node inside browse tree")));
      await session.step(56, "Then the \"Models\" view should be current", () => viewIsCurrent(page, "Models"));
      await session.step(57, "When user switches to the \"forge-test-iris-{time}\" table view", () => switchTableView(page, session.text("forge-test-iris-{time}")));
      await session.step(58, "And user picks \"ML > Forge > Models\" from the top menu", () => pickFromTopMenu(page, "ML > Forge > Models"));
      await session.step(59, "Then the \"Forge\" view should be current", () => viewIsCurrent(page, "Forge"));
      await session.step(60, "And there should be 1 open view named \"Forge\"", () => openViewsNamed(page, 1, "Forge"));
      await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Methods lists the installed methods with their package, type, roles and settings; Engine appears nowhere (1.3)", async () => {
      await session.step(65, "Given user switches to the \"Forge\" view", () => switchView(page, "Forge"));
      await session.step(66, "Then first grid viewer should show the columns \"Method, Package, Method type, Roles, Hyperparameters\"", () => gridShowsColumns(page, el("first grid viewer"), "Method, Package, Method type, Roles, Hyperparameters"));
      await session.step(67, "And the method \"CatBoost\" should be listed above the method \"Decision Tree\" in first grid viewer", () => methodAbove(page, "CatBoost", "Decision Tree", el("first grid viewer")));
      await session.step(68, "And the \"Package\" cell of the method \"CatBoost\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "CatBoost", el("first grid viewer"), "Eda"));
      await session.step(69, "And the \"Package\" cell of the method \"Decision Tree\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "Decision Tree", el("first grid viewer"), "Eda"));
      await session.step(70, "And the \"Package\" cell of the method \"Extra Trees\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "Extra Trees", el("first grid viewer"), "Eda"));
      await session.step(71, "And the \"Package\" cell of the method \"Hist Gradient Boosting\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "Hist Gradient Boosting", el("first grid viewer"), "Eda"));
      await session.step(72, "And the \"Package\" cell of the method \"k-Nearest Neighbors\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "k-Nearest Neighbors", el("first grid viewer"), "Eda"));
      await session.step(73, "And the \"Package\" cell of the method \"LDA\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "LDA", el("first grid viewer"), "Eda"));
      await session.step(74, "And the \"Package\" cell of the method \"Linear Regression\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "Linear Regression", el("first grid viewer"), "Eda"));
      await session.step(75, "And the \"Package\" cell of the method \"Naive Bayes\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "Naive Bayes", el("first grid viewer"), "Eda"));
      await session.step(76, "And the \"Package\" cell of the method \"Neural Network\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "Neural Network", el("first grid viewer"), "Eda"));
      await session.step(77, "And the \"Package\" cell of the method \"PLS Regression\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "PLS Regression", el("first grid viewer"), "Eda"));
      await session.step(78, "And the \"Package\" cell of the method \"Random Forest\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "Random Forest", el("first grid viewer"), "Eda"));
      await session.step(79, "And the \"Package\" cell of the method \"RealMLP\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "RealMLP", el("first grid viewer"), "Eda"));
      await session.step(80, "And the \"Package\" cell of the method \"Softmax\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "Softmax", el("first grid viewer"), "Eda"));
      await session.step(81, "And the \"Package\" cell of the method \"SVM\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "SVM", el("first grid viewer"), "Eda"));
      await session.step(82, "And the \"Package\" cell of the method \"TabM\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "TabM", el("first grid viewer"), "Eda"));
      await session.step(83, "And the \"Package\" cell of the method \"XGBoost\" in first grid viewer should read \"Eda\"", () => methodCellReads(page, "Package", "XGBoost", el("first grid viewer"), "Eda"));
      await session.step(84, "And the \"Method type\" cell of the method \"CatBoost\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "CatBoost", el("first grid viewer"), "function"));
      await session.step(85, "And the \"Method type\" cell of the method \"Decision Tree\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "Decision Tree", el("first grid viewer"), "function"));
      await session.step(86, "And the \"Method type\" cell of the method \"Extra Trees\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "Extra Trees", el("first grid viewer"), "function"));
      await session.step(87, "And the \"Method type\" cell of the method \"Hist Gradient Boosting\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "Hist Gradient Boosting", el("first grid viewer"), "function"));
      await session.step(88, "And the \"Method type\" cell of the method \"k-Nearest Neighbors\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "k-Nearest Neighbors", el("first grid viewer"), "function"));
      await session.step(89, "And the \"Method type\" cell of the method \"LDA\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "LDA", el("first grid viewer"), "function"));
      await session.step(90, "And the \"Method type\" cell of the method \"Linear Regression\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "Linear Regression", el("first grid viewer"), "function"));
      await session.step(91, "And the \"Method type\" cell of the method \"Naive Bayes\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "Naive Bayes", el("first grid viewer"), "function"));
      await session.step(92, "And the \"Method type\" cell of the method \"Neural Network\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "Neural Network", el("first grid viewer"), "function"));
      await session.step(93, "And the \"Method type\" cell of the method \"PLS Regression\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "PLS Regression", el("first grid viewer"), "function"));
      await session.step(94, "And the \"Method type\" cell of the method \"Random Forest\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "Random Forest", el("first grid viewer"), "function"));
      await session.step(95, "And the \"Method type\" cell of the method \"RealMLP\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "RealMLP", el("first grid viewer"), "function"));
      await session.step(96, "And the \"Method type\" cell of the method \"Softmax\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "Softmax", el("first grid viewer"), "function"));
      await session.step(97, "And the \"Method type\" cell of the method \"SVM\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "SVM", el("first grid viewer"), "function"));
      await session.step(98, "And the \"Method type\" cell of the method \"TabM\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "TabM", el("first grid viewer"), "function"));
      await session.step(99, "And the \"Method type\" cell of the method \"XGBoost\" in first grid viewer should read \"function\"", () => methodCellReads(page, "Method type", "XGBoost", el("first grid viewer"), "function"));
      await session.step(100, "And the \"Roles\" cell of the method \"XGBoost\" in first grid viewer should read \"train, apply, isApplicable, isInteractive\"", () => methodCellReads(page, "Roles", "XGBoost", el("first grid viewer"), "train, apply, isApplicable, isInteractive"));
      await session.step(101, "And the \"Hyperparameters\" cell of the method \"XGBoost\" in first grid viewer should read \"iterations, eta, maxDepth, lambda, alpha\"", () => methodCellReads(page, "Hyperparameters", "XGBoost", el("first grid viewer"), "iterations, eta, maxDepth, lambda, alpha"));
      await session.step(102, "And the \"Roles\" cell of the method \"PLS Regression\" in first grid viewer should read \"train, apply, isApplicable, isInteractive, visualize\"", () => methodCellReads(page, "Roles", "PLS Regression", el("first grid viewer"), "train, apply, isApplicable, isInteractive, visualize"));
      await session.step(103, "And the \"Hyperparameters\" cell of the method \"PLS Regression\" in first grid viewer should read \"components\"", () => methodCellReads(page, "Hyperparameters", "PLS Regression", el("first grid viewer"), "components"));
      await session.step(104, "And forge app should not contain text \"Engine\"", () => shouldNotContainText(page, el("forge app"), "Engine"));
      await session.step(105, "And no drawn cell or header of first grid viewer should contain \"Engine\"", () => gridShowsNoText(page, el("first grid viewer"), "Engine"));
      await session.step(106, "And no drawn cell or header of second grid viewer should contain \"Engine\"", () => gridShowsNoText(page, el("second grid viewer"), "Engine"));
      await session.step(107, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The catalog shows its captions and icons; the Methods grid explains itself in tooltips (1.4)", async () => {
      await session.step(111, "Given user switches to the \"Forge\" view", () => switchView(page, "Forge"));
      await session.step(112, "Then second grid viewer should show the columns \"Name, Method, Task, Target, Data storage, Training rows, Tags, Created\"", () => gridShowsColumns(page, el("second grid viewer"), "Name, Method, Task, Target, Data storage, Training rows, Tags, Created"));
      await session.step(113, "And \"Refresh\" icon in forge app should be enabled", () => shouldBe(page, el("\"Refresh\" icon in forge app"), "enabled"));
      await session.step(114, "And \"Apply model\" icon in forge app should be disabled", () => shouldBe(page, el("\"Apply model\" icon in forge app"), "disabled"));
      await session.step(115, "And \"Delete model\" icon in forge app should be disabled", () => shouldBe(page, el("\"Delete model\" icon in forge app"), "disabled"));
      await session.step(116, "And \"Compare in a new view\" icon in forge app should be disabled", () => shouldBe(page, el("\"Compare in a new view\" icon in forge app"), "disabled"));
      await session.step(117, "And \"Applicable to\" input in forge app should be visible", () => shouldBe(page, el("\"Applicable to\" input in forge app"), "visible"));
      await session.step(118, "And \"Applicable to\" input in forge app should have value \"\"", () => shouldHaveValue(page, el("\"Applicable to\" input in forge app"), ""));
      await session.step(119, "When user hovers over \"Refresh\" icon in forge app", () => hoverOver(page, el("\"Refresh\" icon in forge app")));
      await session.step(120, "Then tooltip should contain text \"Refresh\"", () => shouldContainText(page, el("tooltip"), "Refresh"));
      await session.step(121, "When user hovers over the \"header Roles\" area of first grid viewer", () => hoverArea(page, "header Roles", el("first grid viewer")));
      await session.step(122, "Then tooltip should contain text \"What the method can do.\"", () => shouldContainText(page, el("tooltip"), "What the method can do."));
      await session.step(123, "When user hovers over the \"Method type\" cell of the method \"Decision Tree\" in first grid viewer", () => hoverMethodCell(page, "Method type", "Decision Tree", el("first grid viewer")));
      await session.step(124, "Then tooltip should contain text \"A package function.\"", () => shouldContainText(page, el("tooltip"), "A package function."));
      await session.step(125, "When user hovers over the \"Roles\" cell of the method \"XGBoost\" in first grid viewer", () => hoverMethodCell(page, "Roles", "XGBoost", el("first grid viewer")));
      await session.step(126, "Then tooltip should contain text \"train: Trains a model on a table and a target column, returns the model\"", () => shouldContainText(page, el("tooltip"), "train: Trains a model on a table and a target column, returns the model"));
      await session.step(127, "And tooltip should contain text \"apply: Applies a trained model to a table, returns the predictions\"", () => shouldContainText(page, el("tooltip"), "apply: Applies a trained model to a table, returns the predictions"));
      await session.step(128, "And tooltip should contain text \"isApplicable: Tells whether the method can learn from the given data\"", () => shouldContainText(page, el("tooltip"), "isApplicable: Tells whether the method can learn from the given data"));
      await session.step(129, "And tooltip should contain text \"isInteractive: Tells whether training is fast enough to rerun on every change\"", () => shouldContainText(page, el("tooltip"), "isInteractive: Tells whether training is fast enough to rerun on every change"));
      await session.step(130, "When user hovers over the \"Roles\" cell of the method \"PLS Regression\" in first grid viewer", () => hoverMethodCell(page, "Roles", "PLS Regression", el("first grid viewer")));
      await session.step(131, "Then tooltip should contain text \"isInteractive: Tells whether training is fast enough to rerun on every change\"", () => shouldContainText(page, el("tooltip"), "isInteractive: Tells whether training is fast enough to rerun on every change"));
      await session.step(132, "And tooltip should contain text \"visualize: Shows method-specific views of a trained model\"", () => shouldContainText(page, el("tooltip"), "visualize: Shows method-specific views of a trained model"));
      await session.step(133, "When user closes the current view", () => closeCurrentView(page));
      await session.step(134, "Then forge app should be absent", () => shouldBe(page, el("forge app"), "absent"));
      await session.step(135, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(136, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The old tool's Train view and Apply dialog still open beside Forge (1.6, 2.9)", async () => {
      await session.step(139, "Given user switches to the \"forge-test-iris-{time}\" table view", () => switchTableView(page, session.text("forge-test-iris-{time}")));
      await session.step(140, "When user picks \"ML > Models > Train Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Train Model..."));
      await session.step(141, "Then the \"Predictive model\" view should be current", () => viewIsCurrent(page, "Predictive model"));
      await session.step(142, "And model preview should be visible", () => shouldBe(page, el("model preview"), "visible"));
      await session.step(143, "When user closes the current view", () => closeCurrentView(page));
      await session.step(144, "And user switches to the \"forge-test-iris-{time}\" table view", () => switchTableView(page, session.text("forge-test-iris-{time}")));
      await session.step(145, "And user picks \"ML > Models > Apply Model...\" from the top menu", () => pickFromTopMenu(page, "ML > Models > Apply Model..."));
      await session.step(146, "Then there should be 1 visible \"Apply predictive model\" dialog", () => visibleCount(page, 1, el("\"Apply predictive model\" dialog")));
      await session.step(147, "And Inputs input in \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("Inputs input in \"Apply predictive model\" dialog"), "visible"));
      await session.step(148, "And \"Columns\" group in \"Apply predictive model\" dialog should be absent", () => shouldBe(page, el("\"Columns\" group in \"Apply predictive model\" dialog"), "absent"));
      await session.step(149, "When user clicks on CANCEL button in \"Apply predictive model\" dialog", () => clickOn(page, el("CANCEL button in \"Apply predictive model\" dialog")));
      await session.step(150, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(151, "When user picks \"ML > Forge > Models\" from the top menu", () => pickFromTopMenu(page, "ML > Forge > Models"));
      await session.step(152, "Then the \"Forge\" view should be current", () => viewIsCurrent(page, "Forge"));
      await session.step(153, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(154, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
