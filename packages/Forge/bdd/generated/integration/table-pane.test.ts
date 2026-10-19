/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/integration/table-pane.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {applicationsCount, newestApplicationShouldHave, noModelsNamed, noRunsOf, shownIsNewestRun} from '../../bindings/server.js';
import {modelLinesList, paneComesAfter, tableIsCurrentObject} from '../../bindings/table-pane.js';
import {countTaskBarEntries, newTrainingShown, noBalloonAtAll, openDemoFile, panelShowsModel, rememberShownTraining, taskBarEntryEnded, trainThroughApiWith} from '../../bindings/train.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, hoverOver, pressKey, rightClickOn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {closeCurrentView, contextPanelOpen, dialogCloses, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, menuLists, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The table pane Forge models: the models that fit a table, their tooltip, menu and run icon, and Train", () => {
  const session = feature(test, "features/integration/table-pane.feature", import.meta.url);
  test("The table pane Forge models: the models that fit a table, their tooltip, menu and run icon, and Train", {tag: ["@journey", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(33, "Given user is logged in", () => loggedIn(page));
    await session.step(34, "And no Forge models named \"forge-test-pane-iris-{time}, forge-test-pane-flower-{time}, forge-test-pane-quality-{time}, forge-test-pane-demog-{time}\" are on the server", () => noModelsNamed(page, session.text("forge-test-pane-iris-{time}, forge-test-pane-flower-{time}, forge-test-pane-quality-{time}, forge-test-pane-demog-{time}")));
    await session.step(35, "And no Forge training runs of the tables \"forge-test-iris-{time}, forge-test-wine-{time}, forge-test-demog-{time}\" are on the server", () => noRunsOf(page, session.text("forge-test-iris-{time}, forge-test-wine-{time}, forge-test-demog-{time}")));
    await session.step(36, "And user opens the demo file \"winequality.csv\" as \"forge-test-wine-{time}\"", () => openDemoFile(page, "winequality.csv", session.text("forge-test-wine-{time}")));
    await session.step(37, "And a Forge model \"forge-test-pane-quality-{time}\" is trained on the current table through the API with:", () => trainThroughApiWith(page, session.text("forge-test-pane-quality-{time}"), [["target","quality"],["method","Linear Regression"]]), [["target","quality"],["method","Linear Regression"]]);
    await session.step(40, "And user opens the demo file \"demog.csv\" as \"forge-test-demog-{time}\"", () => openDemoFile(page, "demog.csv", session.text("forge-test-demog-{time}")));
    await session.step(41, "And a Forge model \"forge-test-pane-demog-{time}\" is trained on the current table through the API with:", () => trainThroughApiWith(page, session.text("forge-test-pane-demog-{time}"), [["target","AGE"],["features","HEIGHT, WEIGHT"],["method","Linear Regression"]]), [["target","AGE"],["features","HEIGHT, WEIGHT"],["method","Linear Regression"]]);
    await session.step(45, "And user opens the demo file \"iris.csv\" as \"forge-test-iris-{time}\"", () => openDemoFile(page, "iris.csv", session.text("forge-test-iris-{time}")));
    await session.step(46, "And a Forge model \"forge-test-pane-iris-{time}\" is trained on the current table through the API with:", () => trainThroughApiWith(page, session.text("forge-test-pane-iris-{time}"), [["target","Species"],["method","XGBoost"]]), [["target","Species"],["method","XGBoost"]]);
    await session.step(49, "And a Forge model \"forge-test-pane-flower-{time}\" is trained on the current table through the API with:", () => trainThroughApiWith(page, session.text("forge-test-pane-flower-{time}"), [["target","Species"],["method","XGBoost"]]), [["target","Species"],["method","XGBoost"]]);
    await session.step(52, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(53, "Then no balloon of any colour should have been shown", () => noBalloonAtAll(page));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
    await run.scenario("A table's context panel lists the models that fit it, newest first, after the built-in panes (6.1, 4.13)", async () => {
      await session.step(57, "Given the table \"forge-test-iris-{time}\" is the current object", () => tableIsCurrentObject(page, session.text("forge-test-iris-{time}")));
      await session.step(58, "Then \"Models\" pane in context panel should be visible", () => shouldBe(page, el("\"Models\" pane in context panel"), "visible"));
      await session.step(59, "And the \"Forge models\" pane in context panel should come after the \"Script\" pane", () => paneComesAfter(page, "Forge models", "Script"));
      await session.step(60, "When user expands \"Forge models\" pane in context panel", () => expand(page, el("\"Forge models\" pane in context panel")));
      await session.step(61, "Then the model lines of \"Forge models\" pane in context panel should list exactly \"forge-test-pane-flower-{time}, forge-test-pane-iris-{time}\" among those ending in \"-{time}\"", () => modelLinesList(page, el("\"Forge models\" pane in context panel"), session.text("forge-test-pane-flower-{time}, forge-test-pane-iris-{time}"), session.text("-{time}")));
      await session.step(62, "And \"forge-test-pane-quality-{time}\" model line in \"Forge models\" pane in context panel should be absent", () => shouldBe(page, el(session.text("\"forge-test-pane-quality-{time}\" model line in \"Forge models\" pane in context panel")), "absent"));
      await session.step(63, "And \"forge-test-pane-demog-{time}\" model line in \"Forge models\" pane in context panel should be absent", () => shouldBe(page, el(session.text("\"forge-test-pane-demog-{time}\" model line in \"Forge models\" pane in context panel")), "absent"));
      await session.step(64, "And menu arrow of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel should be visible", () => shouldBe(page, el(session.text("menu arrow of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel")), "visible"));
      await session.step(65, "And \"Train\" button in \"Forge models\" pane in context panel should be visible", () => shouldBe(page, el("\"Train\" button in \"Forge models\" pane in context panel"), "visible"));
      await session.step(66, "When user hovers over run icon of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel", () => hoverOver(page, el(session.text("run icon of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel"))));
      await session.step(67, "Then tooltip should contain text \"Run prediction\"", () => shouldContainText(page, el("tooltip"), "Run prediction"));
      await session.step(68, "When user hovers over \"Train\" button in \"Forge models\" pane in context panel", () => hoverOver(page, el("\"Train\" button in \"Forge models\" pane in context panel")));
      await session.step(69, "Then tooltip should contain text \"Train a new model on this table.\"", () => shouldContainText(page, el("tooltip"), "Train a new model on this table."));
      await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A line's name shows the model's tooltip and the row's menu (6.1)", async () => {
      await session.step(74, "When user hovers over name of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel", () => hoverOver(page, el(session.text("name of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel"))));
      await session.step(75, "Then tooltip should contain text \"forge-test-pane-iris-{time}\"", () => shouldContainText(page, el("tooltip"), session.text("forge-test-pane-iris-{time}")));
      await session.step(76, "And tooltip should contain text \"Created by\"", () => shouldContainText(page, el("tooltip"), "Created by"));
      await session.step(77, "And value of \"Trained on\" detail in tooltip should have the text \"forge-test-iris-{time} (150 rows)\"", () => shouldHaveText(page, el("value of \"Trained on\" detail in tooltip"), session.text("forge-test-iris-{time} (150 rows)")));
      await session.step(78, "And value of \"Last run\" detail in tooltip should have the text \"Never\"", () => shouldHaveText(page, el("value of \"Last run\" detail in tooltip"), "Never"));
      await session.step(79, "And value of \"Tags\" detail in tooltip should have the text \"\"", () => shouldHaveText(page, el("value of \"Tags\" detail in tooltip"), ""));
      await session.step(80, "And value of \"Inputs\" detail in tooltip should have the text \"Sepal.Length, Sepal.Width, Petal.Length, Petal.Width\"", () => shouldHaveText(page, el("value of \"Inputs\" detail in tooltip"), "Sepal.Length, Sepal.Width, Petal.Length, Petal.Width"));
      await session.step(81, "And value of \"Outputs\" detail in tooltip should have the text \"Species\"", () => shouldHaveText(page, el("value of \"Outputs\" detail in tooltip"), "Species"));
      await session.step(82, "And value of \"Method\" detail in tooltip should have the text \"XGBoost\"", () => shouldHaveText(page, el("value of \"Method\" detail in tooltip"), "XGBoost"));
      await session.step(83, "And value of \"Applicable to\" detail in tooltip should have the text \"forge-test-iris-{time}\"", () => shouldHaveText(page, el("value of \"Applicable to\" detail in tooltip"), session.text("forge-test-iris-{time}")));
      await session.step(84, "When user right-clicks on name of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel", () => rightClickOn(page, el(session.text("name of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel"))));
      await session.step(85, "Then the open menu should list \"Open\"", () => menuLists(page, "Open"));
      await session.step(86, "And the open menu should list \"Share...\"", () => menuLists(page, "Share..."));
      await session.step(87, "And the open menu should list \"History\"", () => menuLists(page, "History"));
      await session.step(88, "And the open menu should list \"Copy link\"", () => menuLists(page, "Copy link"));
      await session.step(89, "And the open menu should list \"Apply...\"", () => menuLists(page, "Apply..."));
      await session.step(90, "And the open menu should list \"Download\"", () => menuLists(page, "Download"));
      await session.step(91, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(92, "Then context menu should be absent", () => shouldBe(page, el("context menu"), "absent"));
      await session.step(93, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The name opens the model in the context panel (6.1)", async () => {
      await session.step(97, "When user clicks on name of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel", () => clickOn(page, el(session.text("name of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel"))));
      await session.step(98, "Then the context panel should show the Forge model \"forge-test-pane-iris-{time}\"", () => panelShowsModel(page, session.text("forge-test-pane-iris-{time}")));
      await session.step(99, "And \"Details\" pane in context panel should be visible", () => shouldBe(page, el("\"Details\" pane in context panel"), "visible"));
      await session.step(100, "Given the table \"forge-test-iris-{time}\" is the current object", () => tableIsCurrentObject(page, session.text("forge-test-iris-{time}")));
      await session.step(101, "Then \"Forge models\" pane in context panel should be visible", () => shouldBe(page, el("\"Forge models\" pane in context panel"), "visible"));
      await session.step(102, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The run icon applies the model to the table (6.1)", async () => {
      await session.step(106, "Given user counts the task bar entries", () => countTaskBarEntries(page));
      await session.step(107, "When user expands \"Forge models\" pane in context panel", () => expand(page, el("\"Forge models\" pane in context panel")));
      await session.step(108, "And user clicks on run icon of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel", () => clickOn(page, el(session.text("run icon of \"forge-test-pane-iris-{time}\" model line in \"Forge models\" pane in context panel"))));
      await session.step(109, "Then there should be 1 visible \"Apply predictive model\" dialog", () => visibleCount(page, 1, el("\"Apply predictive model\" dialog")));
      await session.step(110, "And Table input in \"Apply predictive model\" dialog should have value \"forge-test-iris-{time}\"", () => shouldHaveValue(page, el("Table input in \"Apply predictive model\" dialog"), session.text("forge-test-iris-{time}")));
      await session.step(111, "And Model input in \"Apply predictive model\" dialog should have value \"forge-test-pane-iris-{time}\"", () => shouldHaveValue(page, el("Model input in \"Apply predictive model\" dialog"), session.text("forge-test-pane-iris-{time}")));
      await session.step(112, "And summary of \"Columns\" group should have the text \"4 of 4 matched\"", () => shouldHaveText(page, el("summary of \"Columns\" group"), "4 of 4 matched"));
      await session.step(113, "When user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(114, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(115, "And the task bar entry that showed \"Predicting Species\" should have ended", () => taskBarEntryEnded(page, "Predicting Species"));
      await session.step(116, "And the table should have a column \"Species (predicted)\"", () => hasColumn(page, "Species (predicted)"));
      await session.step(117, "And an info balloon containing 'Added the column \"Species (predicted)\" to forge-test-iris-{time}.' should have been shown", () => infoBalloonText(page, session.text("Added the column \"Species (predicted)\" to forge-test-iris-{time}.")));
      await session.step(118, "And the model \"forge-test-pane-iris-{time}\" should have 1 application on the server", () => applicationsCount(page, session.text("forge-test-pane-iris-{time}"), 1));
      await session.step(119, "And the newest application of the model \"forge-test-pane-iris-{time}\" should have:", () => newestApplicationShouldHave(page, session.text("forge-test-pane-iris-{time}"), [["source","ui"],["status","completed"],["column_name","Species (predicted)"]]), [["source","ui"],["status","completed"],["column_name","Species (predicted)"]]);
      await session.step(123, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(124, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Each table lists only the models that fit it; a table no model fits says so (6.1)", async () => {
      await session.step(127, "When user switches to the \"forge-test-wine-{time}\" table view", () => switchTableView(page, session.text("forge-test-wine-{time}")));
      await session.step(128, "Given the table \"forge-test-wine-{time}\" is the current object", () => tableIsCurrentObject(page, session.text("forge-test-wine-{time}")));
      await session.step(129, "When user expands \"Forge models\" pane in context panel", () => expand(page, el("\"Forge models\" pane in context panel")));
      await session.step(130, "Then the model lines of \"Forge models\" pane in context panel should list exactly \"forge-test-pane-quality-{time}\" among those ending in \"-{time}\"", () => modelLinesList(page, el("\"Forge models\" pane in context panel"), session.text("forge-test-pane-quality-{time}"), session.text("-{time}")));
      await session.step(131, "When user switches to the \"forge-test-demog-{time}\" table view", () => switchTableView(page, session.text("forge-test-demog-{time}")));
      await session.step(132, "Given the table \"forge-test-demog-{time}\" is the current object", () => tableIsCurrentObject(page, session.text("forge-test-demog-{time}")));
      await session.step(133, "When user expands \"Forge models\" pane in context panel", () => expand(page, el("\"Forge models\" pane in context panel")));
      await session.step(134, "Then the model lines of \"Forge models\" pane in context panel should list exactly \"forge-test-pane-demog-{time}\" among those ending in \"-{time}\"", () => modelLinesList(page, el("\"Forge models\" pane in context panel"), session.text("forge-test-pane-demog-{time}"), session.text("-{time}")));
      await session.step(135, "Given user opens the demo file \"cars.csv\" as \"forge-test-cars-{time}\"", () => openDemoFile(page, "cars.csv", session.text("forge-test-cars-{time}")));
      await session.step(136, "And the table \"forge-test-cars-{time}\" is the current object", () => tableIsCurrentObject(page, session.text("forge-test-cars-{time}")));
      await session.step(137, "When user expands \"Forge models\" pane in context panel", () => expand(page, el("\"Forge models\" pane in context panel")));
      await session.step(138, "Then the model lines of \"Forge models\" pane in context panel should list exactly \"forge-test-pane-demog-{time}\" among those ending in \"-{time}\"", () => modelLinesList(page, el("\"Forge models\" pane in context panel"), session.text("forge-test-pane-demog-{time}"), session.text("-{time}")));
      await session.step(139, "Given user opens the demo file \"country-languages.csv\" as \"forge-test-langs-{time}\"", () => openDemoFile(page, "country-languages.csv", session.text("forge-test-langs-{time}")));
      await session.step(140, "And the table \"forge-test-langs-{time}\" is the current object", () => tableIsCurrentObject(page, session.text("forge-test-langs-{time}")));
      await session.step(141, "When user expands \"Forge models\" pane in context panel", () => expand(page, el("\"Forge models\" pane in context panel")));
      await session.step(142, "Then \"Forge models\" pane in context panel should contain text \"No model fits this table.\"", () => shouldContainText(page, el("\"Forge models\" pane in context panel"), "No model fits this table."));
      await session.step(143, "And the model lines of \"Forge models\" pane in context panel should list exactly \"\" among those ending in \"-{time}\"", () => modelLinesList(page, el("\"Forge models\" pane in context panel"), "", session.text("-{time}")));
      await session.step(144, "And \"Train\" button in \"Forge models\" pane in context panel should be visible", () => shouldBe(page, el("\"Train\" button in \"Forge models\" pane in context panel"), "visible"));
      await session.step(145, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(146, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Train in the pane opens the Train view on the pane's table, not the current view's (6.1)", async () => {
      await session.step(149, "When user switches to the \"forge-test-iris-{time}\" table view", () => switchTableView(page, session.text("forge-test-iris-{time}")));
      await session.step(150, "Given the table \"forge-test-demog-{time}\" is the current object", () => tableIsCurrentObject(page, session.text("forge-test-demog-{time}")));
      await session.step(151, "And user remembers the shown training", () => rememberShownTraining(page));
      await session.step(152, "When user expands \"Forge models\" pane in context panel", () => expand(page, el("\"Forge models\" pane in context panel")));
      await session.step(153, "And user clicks on \"Train\" button in \"Forge models\" pane in context panel", () => clickOn(page, el("\"Train\" button in \"Forge models\" pane in context panel")));
      await session.step(154, "Then the \"Predictive model\" view should be current", () => viewIsCurrent(page, "Predictive model"));
      await session.step(155, "And Table input should have value \"forge-test-demog-{time}\"", () => shouldHaveValue(page, el("Table input"), session.text("forge-test-demog-{time}")));
      await session.step(156, "And a new training's result should be shown", () => newTrainingShown(page));
      await session.step(157, "And the shown training should be the newest training run of the table \"forge-test-demog-{time}\", with:", () => shownIsNewestRun(page, session.text("forge-test-demog-{time}"), [["status","completed"],["engine_name","XGBoost"],["target_name","SEVERITY"],["source","ui"]]), [["status","completed"],["engine_name","XGBoost"],["target_name","SEVERITY"],["source","ui"]]);
      await session.step(162, "When user closes the current view", () => closeCurrentView(page));
      await session.step(163, "Then no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(164, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
