/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/apply/apply-winequality.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {maeMatches, namedTableNoColumnMatching} from '../../bindings/apply.js';
import {scrollGridToColumn} from '../../bindings/catalog.js';
import {applicationsCount, newestApplicationShouldHave, noModelsNamed, noRunsOf, startsWithNewestApplication} from '../../bindings/server.js';
import {cancelTaskBarEntry, clickCatalogRow, countTaskBarEntries, noBalloonAtAll, openDemoFile, panelShowsModel, taskBarEntryCount, taskBarEntryEnded, taskBarEntryNotShowedText, taskBarEntryShowedText, trainThroughApiWith} from '../../bindings/train.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {contextPanelOpen, dialogCloses, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, readingReads, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Applying the quality regressor to winequality, and cancelling an application", () => {
  const session = feature(test, "features/apply/apply-winequality.feature", import.meta.url);
  test("Applying the quality regressor to winequality, and cancelling an application", {tag: ["@journey", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And no Forge models named \"forge-test-apply-quality-{time}\" are on the server", () => noModelsNamed(page, session.text("forge-test-apply-quality-{time}")));
    await session.step(25, "And no Forge training runs of the table \"forge-test-wine-{time}\" are on the server", () => noRunsOf(page, session.text("forge-test-wine-{time}")));
    await session.step(26, "And user opens the demo file \"winequality.csv\" as \"forge-test-wine-{time}\"", () => openDemoFile(page, "winequality.csv", session.text("forge-test-wine-{time}")));
    await session.step(27, "And a Forge model \"forge-test-apply-quality-{time}\" is trained on the current table through the API with:", () => trainThroughApiWith(page, session.text("forge-test-apply-quality-{time}"), [["target","quality"],["method","XGBoost"]]), [["target","quality"],["method","XGBoost"]]);
    await session.step(30, "Then no balloon of any colour should have been shown", () => noBalloonAtAll(page));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
    await run.scenario("The quality model is preset and matches every column; its predictions are whole numbers with the training's MAE (3.2, 3.5)", async () => {
      await session.step(34, "Given user counts the task bar entries", () => countTaskBarEntries(page));
      await session.step(35, "When user picks \"ML > Forge > Apply...\" from the top menu", () => pickFromTopMenu(page, "ML > Forge > Apply..."));
      await session.step(36, "Then there should be 1 visible \"Apply predictive model\" dialog", () => visibleCount(page, 1, el("\"Apply predictive model\" dialog")));
      await session.step(37, "And Table input in \"Apply predictive model\" dialog should have value \"forge-test-wine-{time}\"", () => shouldHaveValue(page, el("Table input in \"Apply predictive model\" dialog"), session.text("forge-test-wine-{time}")));
      await session.step(38, "And Model input in \"Apply predictive model\" dialog should have value \"forge-test-apply-quality-{time}\"", () => shouldHaveValue(page, el("Model input in \"Apply predictive model\" dialog"), session.text("forge-test-apply-quality-{time}")));
      await session.step(39, "And summary of \"Columns\" group should have the text \"11 of 11 matched\"", () => shouldHaveText(page, el("summary of \"Columns\" group"), "11 of 11 matched"));
      await session.step(40, "And summary of \"Columns\" group should be valid", () => shouldBe(page, el("summary of \"Columns\" group"), "valid"));
      await session.step(41, "And \"Columns\" group should be collapsed", () => shouldBe(page, el("\"Columns\" group"), "collapsed"));
      await session.step(42, "And \"Missing values\" input in \"Apply predictive model\" dialog should be hidden", () => shouldBe(page, el("\"Missing values\" input in \"Apply predictive model\" dialog"), "hidden"));
      await session.step(43, "And \"More options\" group should be collapsed", () => shouldBe(page, el("\"More options\" group"), "collapsed"));
      await session.step(44, "And OK button in \"Apply predictive model\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Apply predictive model\" dialog"), "enabled"));
      await session.step(45, "When user expands \"Columns\" group", () => expand(page, el("\"Columns\" group")));
      await session.step(46, "Then editor of \"fixed.acidity\" input in \"Apply predictive model\" dialog should have text \"fixed.acidity\"", () => shouldHaveText(page, el("editor of \"fixed.acidity\" input in \"Apply predictive model\" dialog"), "fixed.acidity"));
      await session.step(47, "And editor of \"alcohol\" input in \"Apply predictive model\" dialog should have text \"alcohol\"", () => shouldHaveText(page, el("editor of \"alcohol\" input in \"Apply predictive model\" dialog"), "alcohol"));
      await session.step(48, "When user expands \"More options\" group", () => expand(page, el("\"More options\" group")));
      await session.step(49, "Then \"Batch size\" input in \"Apply predictive model\" dialog should have value \"10000\"", () => shouldHaveValue(page, el("\"Batch size\" input in \"Apply predictive model\" dialog"), "10000"));
      await session.step(50, "When user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(51, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(52, "And the task bar entry that showed \"Predicting quality\" should have ended", () => taskBarEntryEnded(page, "Predicting quality"));
      await session.step(53, "And the task bar should have shown 1 entry \"Predicting quality\"", () => taskBarEntryCount(page, 1, "Predicting quality"));
      await session.step(54, "And a new column \"quality (predicted)\" should have been added", () => newColumnNamed(page, "quality (predicted)"));
      await session.step(55, "And \"quality (predicted)\" column should have no missing values", () => columnComplete(page, "quality (predicted)"));
      await session.step(56, "And an info balloon containing 'Added the column \"quality (predicted)\" to forge-test-wine-{time}.' should have been shown", () => infoBalloonText(page, session.text("Added the column \"quality (predicted)\" to forge-test-wine-{time}.")));
      await session.step(57, "And every value of \"quality (predicted)\" column should match \"^-?\\d+$\"", () => everyValueMatches(page, "quality (predicted)", "^-?\\d+$"));
      await session.step(58, "And the mean absolute difference of \"quality\" and \"quality (predicted)\" columns of the table \"forge-test-wine-{time}\" should equal the Train MAE of the Forge model \"forge-test-apply-quality-{time}\" within 0.001", () => maeMatches(page, "quality", "quality (predicted)", session.text("forge-test-wine-{time}"), session.text("forge-test-apply-quality-{time}"), 0.001));
      await session.step(59, "And the model \"forge-test-apply-quality-{time}\" should have 1 application on the server", () => applicationsCount(page, session.text("forge-test-apply-quality-{time}"), 1));
      await session.step(60, "And the newest application of the model \"forge-test-apply-quality-{time}\" should have:", () => newestApplicationShouldHave(page, session.text("forge-test-apply-quality-{time}"), [["table_name",session.text("forge-test-wine-{time}")],["row_count","6497"],["skipped_rows","0"],["column_name","quality (predicted)"],["source","ui"],["status","completed"],["duration_ms","number"]]), [["table_name",session.text("forge-test-wine-{time}")],["row_count","6497"],["skipped_rows","0"],["column_name","quality (predicted)"],["source","ui"],["status","completed"],["duration_ms","number"]]);
      await session.step(68, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cancelled application adds no column and is recorded; the model's Activity and Last run show it (3.6, 4.4)", async () => {
      await session.step(72, "Given user counts the task bar entries", () => countTaskBarEntries(page));
      await session.step(73, "When user picks \"ML > Forge > Apply...\" from the top menu", () => pickFromTopMenu(page, "ML > Forge > Apply..."));
      await session.step(74, "Then there should be 1 visible \"Apply predictive model\" dialog", () => visibleCount(page, 1, el("\"Apply predictive model\" dialog")));
      await session.step(75, "And Model input in \"Apply predictive model\" dialog should have value \"forge-test-apply-quality-{time}\"", () => shouldHaveValue(page, el("Model input in \"Apply predictive model\" dialog"), session.text("forge-test-apply-quality-{time}")));
      await session.step(76, "When user expands \"More options\" group", () => expand(page, el("\"More options\" group")));
      await session.step(77, "And user enters \"1\" into \"Batch size\" input in \"Apply predictive model\" dialog", () => enterInto(page, "1", el("\"Batch size\" input in \"Apply predictive model\" dialog")));
      await session.step(78, "And user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(79, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(80, "When user cancels the \"Predicting quality\" task bar entry", () => cancelTaskBarEntry(page, "Predicting quality"));
      await session.step(81, "Then a warning balloon containing \"Application was cancelled.\" should have been shown", () => warningBalloonText(page, "Application was cancelled."));
      await session.step(82, "And the task bar entry that showed \"Predicting quality\" should have ended", () => taskBarEntryEnded(page, "Predicting quality"));
      await session.step(83, "And the task bar should have shown 1 entry \"Predicting quality\"", () => taskBarEntryCount(page, 1, "Predicting quality"));
      await session.step(84, "And the task bar entry that showed \"Predicting quality\" should have shown \"Rows 1-1 of 6497\"", () => taskBarEntryShowedText(page, "Predicting quality", "Rows 1-1 of 6497"));
      await session.step(85, "And the task bar entry that showed \"Predicting quality\" should not have shown \"Rows 6497-6497 of 6497\"", () => taskBarEntryNotShowedText(page, "Predicting quality", "Rows 6497-6497 of 6497"));
      await session.step(86, "And the table \"forge-test-wine-{time}\" should have no column matching \"^quality \\(predicted \\d+\\)$\"", () => namedTableNoColumnMatching(page, session.text("forge-test-wine-{time}"), "^quality \\(predicted \\d+\\)$"));
      await session.step(87, "And the newest application of the model \"forge-test-apply-quality-{time}\" should have:", () => newestApplicationShouldHave(page, session.text("forge-test-apply-quality-{time}"), [["table_name",session.text("forge-test-wine-{time}")],["status","cancelled"],["column_name","empty"]]), [["table_name",session.text("forge-test-wine-{time}")],["status","cancelled"],["column_name","empty"]]);
      await session.step(91, "And the model \"forge-test-apply-quality-{time}\" should have 2 applications on the server", () => applicationsCount(page, session.text("forge-test-apply-quality-{time}"), 2));
      await session.step(92, "When user picks \"ML > Forge > Models\" from the top menu", () => pickFromTopMenu(page, "ML > Forge > Models"));
      await session.step(93, "Then the \"Forge\" view should be current", () => viewIsCurrent(page, "Forge"));
      await session.step(94, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(95, "When user clicks on the Forge model \"forge-test-apply-quality-{time}\" in second grid viewer", () => clickCatalogRow(page, session.text("forge-test-apply-quality-{time}"), el("second grid viewer")));
      await session.step(96, "Then the context panel should show the Forge model \"forge-test-apply-quality-{time}\"", () => panelShowsModel(page, session.text("forge-test-apply-quality-{time}")));
      await session.step(97, "And value of \"Applications\" detail in context panel should have the text \"2\"", () => shouldHaveText(page, el("value of \"Applications\" detail in context panel"), "2"));
      await session.step(98, "And value of \"Last run\" detail in context panel should contain text \"(cancelled)\"", () => shouldContainText(page, el("value of \"Last run\" detail in context panel"), "(cancelled)"));
      await session.step(99, "And value of \"Last run\" detail in context panel should start with the time of the newest application of the model \"forge-test-apply-quality-{time}\"", () => startsWithNewestApplication(page, el("value of \"Last run\" detail in context panel"), session.text("forge-test-apply-quality-{time}")));
      await session.step(100, "When user expands \"Activity\" pane in context panel", () => expand(page, el("\"Activity\" pane in context panel")));
      await session.step(101, "Then \"Activity\" pane in context panel should contain text \"2 applications\"", () => shouldContainText(page, el("\"Activity\" pane in context panel"), "2 applications"));
      await session.step(102, "When user scrolls grid viewer in \"Activity\" pane in context panel to the \"Prediction column\" column", () => scrollGridToColumn(page, el("grid viewer in \"Activity\" pane in context panel"), "Prediction column"));
      await session.step(103, "Then the \"text of cell 1 of Prediction column\" reading of grid viewer in \"Activity\" pane in context panel should be \"\"", () => readingReads(page, "text of cell 1 of Prediction column", el("grid viewer in \"Activity\" pane in context panel"), ""));
      await session.step(104, "And the \"text of cell 2 of Prediction column\" reading of grid viewer in \"Activity\" pane in context panel should be \"quality (predicted)\"", () => readingReads(page, "text of cell 2 of Prediction column", el("grid viewer in \"Activity\" pane in context panel"), "quality (predicted)"));
      await session.step(105, "When user scrolls grid viewer in \"Activity\" pane in context panel to the \"Status\" column", () => scrollGridToColumn(page, el("grid viewer in \"Activity\" pane in context panel"), "Status"));
      await session.step(106, "Then the \"text of cell 1 of Status\" reading of grid viewer in \"Activity\" pane in context panel should be \"cancelled\"", () => readingReads(page, "text of cell 1 of Status", el("grid viewer in \"Activity\" pane in context panel"), "cancelled"));
      await session.step(107, "And the \"text of cell 2 of Status\" reading of grid viewer in \"Activity\" pane in context panel should be \"completed\"", () => readingReads(page, "text of cell 2 of Status", el("grid viewer in \"Activity\" pane in context panel"), "completed"));
      await session.step(108, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(109, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
