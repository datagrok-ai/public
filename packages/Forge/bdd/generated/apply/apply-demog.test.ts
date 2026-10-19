/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/apply/apply-demog.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {emptyExactlyInGaps, imputedPredictions, missingValueCount, radioChosen, rowsWithGaps, valueMatches} from '../../bindings/apply.js';
import {applicationsCount, newestApplicationShouldHave, newestRunShouldHave, noModelsNamed, noRunsOf} from '../../bindings/server.js';
import {countTaskBarEntries, noBalloonAtAll, openDemoFile, taskBarEntryCount, taskBarEntryEnded, trainThroughApiWith} from '../../bindings/train.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, hoverOver, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {dialogCloses} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Applying on data with natural gaps: demog's HEIGHT and WEIGHT, skipped by default or imputed", () => {
  const session = feature(test, "features/apply/apply-demog.feature", import.meta.url);
  test("Applying on data with natural gaps: demog's HEIGHT and WEIGHT, skipped by default or imputed", {tag: ["@journey", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And no Forge models named \"forge-test-apply-demog-{time}\" are on the server", () => noModelsNamed(page, session.text("forge-test-apply-demog-{time}")));
    await session.step(29, "And no Forge training runs of the table \"forge-test-demog-{time}\" are on the server", () => noRunsOf(page, session.text("forge-test-demog-{time}")));
    await session.step(30, "And user opens the demo file \"demog.csv\" as \"forge-test-demog-{time}\"", () => openDemoFile(page, "demog.csv", session.text("forge-test-demog-{time}")));
    await session.step(31, "And a Forge model \"forge-test-apply-demog-{time}\" is trained on the current table through the API with:", () => trainThroughApiWith(page, session.text("forge-test-apply-demog-{time}"), [["target","AGE"],["features","HEIGHT, WEIGHT"],["method","Linear Regression"]]), [["target","AGE"],["features","HEIGHT, WEIGHT"],["method","Linear Regression"]]);
    await session.step(35, "Then no balloon of any colour should have been shown", () => noBalloonAtAll(page));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
    await run.scenario("Rows with a natural gap in HEIGHT or WEIGHT get no prediction by default, and are counted (3.15, steps 2-3)", async () => {
      await session.step(39, "Given user counts the task bar entries", () => countTaskBarEntries(page));
      await session.step(40, "Then the newest training run of the table \"forge-test-demog-{time}\" should have:", () => newestRunShouldHave(page, session.text("forge-test-demog-{time}"), [["source","api"],["row_count","5097"],["options.missingValues.mode","skip"],["options.missingValues.skippedRows","753"]]), [["source","api"],["row_count","5097"],["options.missingValues.mode","skip"],["options.missingValues.skippedRows","753"]]);
      await session.step(45, "And \"HEIGHT\" column should have 751 missing values", () => missingValueCount(page, "HEIGHT", 751));
      await session.step(46, "And \"WEIGHT\" column should have 1 missing value", () => missingValueCount(page, "WEIGHT", 1));
      await session.step(47, "And the rows with a missing value in any of \"HEIGHT, WEIGHT\" columns should be 752", () => rowsWithGaps(page, "HEIGHT, WEIGHT", 752));
      await session.step(48, "And the value of \"AGE\" column in row 1193 should be \"\"", () => valueInRow(page, "AGE", 1193, ""));
      await session.step(49, "And the value of \"HEIGHT\" column in row 965 should be \"\"", () => valueInRow(page, "HEIGHT", 965, ""));
      await session.step(50, "And the value of \"WEIGHT\" column in row 2222 should be \"\"", () => valueInRow(page, "WEIGHT", 2222, ""));
      await session.step(51, "When user picks \"ML > Forge > Apply...\" from the top menu", () => pickFromTopMenu(page, "ML > Forge > Apply..."));
      await session.step(52, "Then there should be 1 visible \"Apply predictive model\" dialog", () => visibleCount(page, 1, el("\"Apply predictive model\" dialog")));
      await session.step(53, "And Table input in \"Apply predictive model\" dialog should have value \"forge-test-demog-{time}\"", () => shouldHaveValue(page, el("Table input in \"Apply predictive model\" dialog"), session.text("forge-test-demog-{time}")));
      await session.step(54, "And Model input in \"Apply predictive model\" dialog should have value \"forge-test-apply-demog-{time}\"", () => shouldHaveValue(page, el("Model input in \"Apply predictive model\" dialog"), session.text("forge-test-apply-demog-{time}")));
      await session.step(55, "And summary of \"Columns\" group should have the text \"2 of 2 matched\"", () => shouldHaveText(page, el("summary of \"Columns\" group"), "2 of 2 matched"));
      await session.step(56, "And summary of \"Columns\" group should be valid", () => shouldBe(page, el("summary of \"Columns\" group"), "valid"));
      await session.step(57, "And \"Columns\" group should be collapsed", () => shouldBe(page, el("\"Columns\" group"), "collapsed"));
      await session.step(58, "And \"Missing values\" input in \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("\"Missing values\" input in \"Apply predictive model\" dialog"), "visible"));
      await session.step(59, "And \"Missing values\" input in \"Apply predictive model\" dialog should have \"Skip rows\" chosen", () => radioChosen(page, el("\"Missing values\" input in \"Apply predictive model\" dialog"), "Skip rows"));
      await session.step(60, "And \"Neighbors\" input in \"Apply predictive model\" dialog should be hidden", () => shouldBe(page, el("\"Neighbors\" input in \"Apply predictive model\" dialog"), "hidden"));
      await session.step(61, "And \"Distance\" input in \"Apply predictive model\" dialog should be hidden", () => shouldBe(page, el("\"Distance\" input in \"Apply predictive model\" dialog"), "hidden"));
      await session.step(62, "And OK button in \"Apply predictive model\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Apply predictive model\" dialog"), "enabled"));
      await session.step(63, "When user hovers over label of \"Missing values\" input in \"Apply predictive model\" dialog", () => hoverOver(page, el("label of \"Missing values\" input in \"Apply predictive model\" dialog")));
      await session.step(64, "Then tooltip should contain text \"Rows with missing values in: HEIGHT: 751, WEIGHT: 1.\"", () => shouldContainText(page, el("tooltip"), "Rows with missing values in: HEIGHT: 751, WEIGHT: 1."));
      await session.step(65, "When user expands \"Columns\" group", () => expand(page, el("\"Columns\" group")));
      await session.step(66, "Then editor of \"HEIGHT\" input in \"Apply predictive model\" dialog should have text \"HEIGHT\"", () => shouldHaveText(page, el("editor of \"HEIGHT\" input in \"Apply predictive model\" dialog"), "HEIGHT"));
      await session.step(67, "And editor of \"WEIGHT\" input in \"Apply predictive model\" dialog should have text \"WEIGHT\"", () => shouldHaveText(page, el("editor of \"WEIGHT\" input in \"Apply predictive model\" dialog"), "WEIGHT"));
      await session.step(68, "When user hovers over label of \"HEIGHT\" input in \"Apply predictive model\" dialog", () => hoverOver(page, el("label of \"HEIGHT\" input in \"Apply predictive model\" dialog")));
      await session.step(69, "Then tooltip should contain text \"'HEIGHT' has 751 missing values.\"", () => shouldContainText(page, el("tooltip"), "'HEIGHT' has 751 missing values."));
      await session.step(70, "When user hovers over label of \"WEIGHT\" input in \"Apply predictive model\" dialog", () => hoverOver(page, el("label of \"WEIGHT\" input in \"Apply predictive model\" dialog")));
      await session.step(71, "Then tooltip should contain text \"'WEIGHT' has 1 missing value.\"", () => shouldContainText(page, el("tooltip"), "'WEIGHT' has 1 missing value."));
      await session.step(72, "When user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(73, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(74, "And the task bar entry that showed \"Predicting AGE\" should have ended", () => taskBarEntryEnded(page, "Predicting AGE"));
      await session.step(75, "And the task bar should have shown 1 entry \"Predicting AGE\"", () => taskBarEntryCount(page, 1, "Predicting AGE"));
      await session.step(76, "And a new column \"AGE (predicted)\" should have been added", () => newColumnNamed(page, "AGE (predicted)"));
      await session.step(77, "And an info balloon containing 'Added the column \"AGE (predicted)\" to forge-test-demog-{time}; 752 rows skipped (missing values).' should have been shown", () => infoBalloonText(page, session.text("Added the column \"AGE (predicted)\" to forge-test-demog-{time}; 752 rows skipped (missing values).")));
      await session.step(78, "And the table should have the columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, AGE (predicted)\"", () => columnsExactly(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, AGE (predicted)"));
      await session.step(79, "And the \"AGE (predicted)\" column should be empty exactly in the rows where any of \"HEIGHT, WEIGHT\" columns is empty", () => emptyExactlyInGaps(page, "AGE (predicted)", "HEIGHT, WEIGHT"));
      await session.step(80, "And the value of \"AGE (predicted)\" column in row 965 should be \"\"", () => valueInRow(page, "AGE (predicted)", 965, ""));
      await session.step(81, "And the value of \"AGE (predicted)\" column in row 2222 should be \"\"", () => valueInRow(page, "AGE (predicted)", 2222, ""));
      await session.step(82, "And the value of \"AGE (predicted)\" column in row 1193 should match \"^-?\\d\"", () => valueMatches(page, "AGE (predicted)", 1193, "^-?\\d"));
      await session.step(83, "And \"AGE (predicted)\" column should have 752 missing values", () => missingValueCount(page, "AGE (predicted)", 752));
      await session.step(84, "And \"HEIGHT\" column should have 751 missing values", () => missingValueCount(page, "HEIGHT", 751));
      await session.step(85, "And \"WEIGHT\" column should have 1 missing value", () => missingValueCount(page, "WEIGHT", 1));
      await session.step(86, "And the model \"forge-test-apply-demog-{time}\" should have 1 application on the server", () => applicationsCount(page, session.text("forge-test-apply-demog-{time}"), 1));
      await session.step(87, "And the newest application of the model \"forge-test-apply-demog-{time}\" should have:", () => newestApplicationShouldHave(page, session.text("forge-test-apply-demog-{time}"), [["table_name",session.text("forge-test-demog-{time}")],["row_count","5850"],["skipped_rows","752"],["column_name","AGE (predicted)"],["source","ui"],["status","completed"],["duration_ms","number"]]), [["table_name",session.text("forge-test-demog-{time}")],["row_count","5850"],["skipped_rows","752"],["column_name","AGE (predicted)"],["source","ui"],["status","completed"],["duration_ms","number"]]);
      await session.step(95, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Impute fills the gaps on a copy and predicts every row; the table keeps its gaps (3.15, step 4)", async () => {
      await session.step(99, "Given user counts the task bar entries", () => countTaskBarEntries(page));
      await session.step(100, "When user picks \"ML > Forge > Apply...\" from the top menu", () => pickFromTopMenu(page, "ML > Forge > Apply..."));
      await session.step(101, "Then there should be 1 visible \"Apply predictive model\" dialog", () => visibleCount(page, 1, el("\"Apply predictive model\" dialog")));
      await session.step(102, "And Model input in \"Apply predictive model\" dialog should have value \"forge-test-apply-demog-{time}\"", () => shouldHaveValue(page, el("Model input in \"Apply predictive model\" dialog"), session.text("forge-test-apply-demog-{time}")));
      await session.step(103, "And \"Missing values\" input in \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("\"Missing values\" input in \"Apply predictive model\" dialog"), "visible"));
      await session.step(104, "And \"Missing values\" input in \"Apply predictive model\" dialog should have \"Skip rows\" chosen", () => radioChosen(page, el("\"Missing values\" input in \"Apply predictive model\" dialog"), "Skip rows"));
      await session.step(105, "When user selects \"Impute\" in \"Missing values\" input in \"Apply predictive model\" dialog", () => selectIn(page, "Impute", el("\"Missing values\" input in \"Apply predictive model\" dialog")));
      await session.step(106, "Then \"Missing values\" input in \"Apply predictive model\" dialog should have \"Impute\" chosen", () => radioChosen(page, el("\"Missing values\" input in \"Apply predictive model\" dialog"), "Impute"));
      await session.step(107, "And \"Neighbors\" input in \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("\"Neighbors\" input in \"Apply predictive model\" dialog"), "visible"));
      await session.step(108, "And \"Neighbors\" input in \"Apply predictive model\" dialog should have value \"4\"", () => shouldHaveValue(page, el("\"Neighbors\" input in \"Apply predictive model\" dialog"), "4"));
      await session.step(109, "And \"Distance\" input in \"Apply predictive model\" dialog should be visible", () => shouldBe(page, el("\"Distance\" input in \"Apply predictive model\" dialog"), "visible"));
      await session.step(110, "And \"Distance\" input in \"Apply predictive model\" dialog should have value \"Euclidean\"", () => shouldHaveValue(page, el("\"Distance\" input in \"Apply predictive model\" dialog"), "Euclidean"));
      await session.step(111, "And OK button in \"Apply predictive model\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Apply predictive model\" dialog"), "enabled"));
      await session.step(112, "When user clicks on OK button in \"Apply predictive model\" dialog", () => clickOn(page, el("OK button in \"Apply predictive model\" dialog")));
      await session.step(113, "Then the \"Apply predictive model\" dialog should close", () => dialogCloses(page, "Apply predictive model"));
      await session.step(114, "And the task bar entry that showed \"Predicting AGE\" should have ended", () => taskBarEntryEnded(page, "Predicting AGE"));
      await session.step(115, "And the task bar should have shown 1 entry \"Predicting AGE\"", () => taskBarEntryCount(page, 1, "Predicting AGE"));
      await session.step(116, "And a new column \"AGE (predicted 2)\" should have been added", () => newColumnNamed(page, "AGE (predicted 2)"));
      await session.step(117, "And an info balloon containing 'Added the column \"AGE (predicted 2)\" to forge-test-demog-{time}.' should have been shown", () => infoBalloonText(page, session.text("Added the column \"AGE (predicted 2)\" to forge-test-demog-{time}.")));
      await session.step(118, "And the table should have the columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, AGE (predicted), AGE (predicted 2)\"", () => columnsExactly(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, AGE (predicted), AGE (predicted 2)"));
      await session.step(119, "And \"AGE (predicted 2)\" column should have no missing values", () => columnComplete(page, "AGE (predicted 2)"));
      await session.step(120, "And the value of \"AGE (predicted 2)\" column in row 965 should match \"^-?\\d\"", () => valueMatches(page, "AGE (predicted 2)", 965, "^-?\\d"));
      await session.step(121, "And the value of \"AGE (predicted 2)\" column in row 2222 should match \"^-?\\d\"", () => valueMatches(page, "AGE (predicted 2)", 2222, "^-?\\d"));
      await session.step(122, "And \"AGE (predicted)\" column should have 752 missing values", () => missingValueCount(page, "AGE (predicted)", 752));
      await session.step(123, "And \"HEIGHT\" column should have 751 missing values", () => missingValueCount(page, "HEIGHT", 751));
      await session.step(124, "And \"WEIGHT\" column should have 1 missing value", () => missingValueCount(page, "WEIGHT", 1));
      await session.step(125, "And the value of \"HEIGHT\" column in row 965 should be \"\"", () => valueInRow(page, "HEIGHT", 965, ""));
      await session.step(126, "And the value of \"WEIGHT\" column in row 2222 should be \"\"", () => valueInRow(page, "WEIGHT", 2222, ""));
      await session.step(127, "And the \"AGE (predicted)\" column should be empty exactly in the rows where any of \"HEIGHT, WEIGHT\" columns is empty", () => emptyExactlyInGaps(page, "AGE (predicted)", "HEIGHT, WEIGHT"));
      await session.step(128, "And the model \"forge-test-apply-demog-{time}\" should have 2 applications on the server", () => applicationsCount(page, session.text("forge-test-apply-demog-{time}"), 2));
      await session.step(129, "And the newest application of the model \"forge-test-apply-demog-{time}\" should have:", () => newestApplicationShouldHave(page, session.text("forge-test-apply-demog-{time}"), [["table_name",session.text("forge-test-demog-{time}")],["row_count","5850"],["skipped_rows","0"],["column_name","AGE (predicted 2)"],["source","ui"],["status","completed"],["duration_ms","number"]]), [["table_name",session.text("forge-test-demog-{time}")],["row_count","5850"],["skipped_rows","0"],["column_name","AGE (predicted 2)"],["source","ui"],["status","completed"],["duration_ms","number"]]);
      await session.step(137, "And the \"AGE (predicted 2)\" column should hold the predictions of the Forge model \"forge-test-apply-demog-{time}\" on \"HEIGHT, WEIGHT\" imputed with 4 neighbors and \"Euclidean\" distance", () => imputedPredictions(page, "AGE (predicted 2)", session.text("forge-test-apply-demog-{time}"), "HEIGHT, WEIGHT", 4, "Euclidean"));
      await session.step(138, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(139, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
