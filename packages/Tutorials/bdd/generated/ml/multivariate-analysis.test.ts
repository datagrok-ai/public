/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/ml/multivariate-analysis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.multivariate-analysis]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, noHintShown, packageInstalled, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Multivariate Analysis tutorial", () => {
  const session = feature(test, "features/ml/multivariate-analysis.feature", import.meta.url);
  test("A learner completes the Multivariate Analysis tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.multivariate-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(24, "And the \"Eda\" package is installed", () => packageInstalled(page, "Eda"));
    await session.step(25, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(26, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(27, "And the \"Multivariate Analysis\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Multivariate Analysis"));
    await session.step(28, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(31, "When user starts the \"Multivariate Analysis\" tutorial", () => startTutorial(page, "Multivariate Analysis"));
    await session.step(32, "Then the tutorial progress should be 1 of 7", () => tutorialProgress(page, 1, 7));
    await session.step(33, "When user picks \"ML > Analyze > Multivariate Analysis...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Multivariate Analysis..."));
    await session.step(34, "Then the tutorial step \"Click on \\\"ML | Analyze | Multivariate Analysis...\\\"\" should be done", () => stepDone(page, "Click on \"ML | Analyze | Multivariate Analysis...\""));
    await session.step(35, "And \"Multivariate Analysis (PLS)\" dialog should be visible", () => shouldBe(page, el("\"Multivariate Analysis (PLS)\" dialog"), "visible"));
    await session.step(37, "When user selects \"price\" in Predict input in \"Multivariate Analysis (PLS)\" dialog", () => selectIn(page, "price", el("Predict input in \"Multivariate Analysis (PLS)\" dialog")));
    await session.step(38, "Then the tutorial step \"Set \\\"Predict\\\" to \\\"price\\\"\" should be done", () => stepDone(page, "Set \"Predict\" to \"price\""));
    await session.step(39, "When user clicks on editor of Using input in \"Multivariate Analysis (PLS)\" dialog", () => clickOn(page, el("editor of Using input in \"Multivariate Analysis (PLS)\" dialog")));
    await session.step(40, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(41, "When user clicks on \"All\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"All\" link in \"Select columns...\" dialog")));
    await session.step(42, "Then \"Select columns...\" dialog should contain text \"16 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "16 checked"));
    await session.step(44, "When user types \"price\" into \"Search\" input in \"Select columns...\" dialog", () => typeInto(page, "price", el("\"Search\" input in \"Select columns...\" dialog")));
    await session.step(45, "And user toggles the \"price\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "price", el("\"Select columns...\" dialog")));
    await session.step(46, "Then \"Select columns...\" dialog should contain text \"15 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "15 checked"));
    await session.step(47, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(48, "Then the tutorial step \"Select all columns, except \\\"price\\\", as \\\"Using\\\"\" should be done", () => stepDone(page, "Select all columns, except \"price\", as \"Using\""));
    await session.step(49, "And editor of Using input in \"Multivariate Analysis (PLS)\" dialog should contain text \"(15)\"", () => shouldContainText(page, el("editor of Using input in \"Multivariate Analysis (PLS)\" dialog"), "(15)"));
    await session.step(50, "When user enters \"3\" into Components input in \"Multivariate Analysis (PLS)\" dialog", () => enterInto(page, "3", el("Components input in \"Multivariate Analysis (PLS)\" dialog")));
    await session.step(51, "Then the tutorial step \"Set the number of components to \\\"3\\\"\" should be done", () => stepDone(page, "Set the number of components to \"3\""));
    await session.step(52, "When user selects \"model\" in Names input in \"Multivariate Analysis (PLS)\" dialog", () => selectIn(page, "model", el("Names input in \"Multivariate Analysis (PLS)\" dialog")));
    await session.step(53, "Then the tutorial step \"Set \\\"Names\\\" to \\\"model\\\"\" should be done", () => stepDone(page, "Set \"Names\" to \"model\""));
    await session.step(55, "When user clicks on RUN button in \"Multivariate Analysis (PLS)\" dialog", () => clickOn(page, el("RUN button in \"Multivariate Analysis (PLS)\" dialog")));
    await session.step(56, "Then the tutorial step \"Click \\\"RUN\\\" and wait for the analysis to complete\" should be done", () => stepDone(page, "Click \"RUN\" and wait for the analysis to complete"));
    await session.step(57, "And the table should have a column \"price (predicted)\"", () => hasColumn(page, "price (predicted)"));
    await session.step(58, "And the table should have a column \"x.score.t3\"", () => hasColumn(page, "x.score.t3"));
    await session.step(59, "And table \"cars(Features Analysis)\" should have 15 rows", () => tableRows(page, "cars(Features Analysis)", 15));
    await session.step(60, "And the open tableview should have 3 scatter plot viewers", () => viewerCount(page, 3, "scatter plot"));
    await session.step(61, "And the open tableview should have 3 bar chart viewers", () => viewerCount(page, 3, "bar chart"));
    await session.step(63, "Then hint popup should contain text \"Observed vs. Predicted\"", () => shouldContainText(page, el("hint popup"), "Observed vs. Predicted"));
    await session.step(64, "When user clicks on \"next\" button in hint popup", () => clickOn(page, el("\"next\" button in hint popup")));
    await session.step(65, "Then hint popup should contain text \"Scores\"", () => shouldContainText(page, el("hint popup"), "Scores"));
    await session.step(66, "When user clicks on \"next\" button in hint popup", () => clickOn(page, el("\"next\" button in hint popup")));
    await session.step(67, "Then hint popup should contain text \"Loadings\"", () => shouldContainText(page, el("hint popup"), "Loadings"));
    await session.step(68, "When user clicks on \"next\" button in hint popup", () => clickOn(page, el("\"next\" button in hint popup")));
    await session.step(69, "Then hint popup should contain text \"Variable Importance\"", () => shouldContainText(page, el("hint popup"), "Variable Importance"));
    await session.step(70, "When user clicks on \"next\" button in hint popup", () => clickOn(page, el("\"next\" button in hint popup")));
    await session.step(71, "Then hint popup should contain text \"Explained Variance\"", () => shouldContainText(page, el("hint popup"), "Explained Variance"));
    await session.step(72, "When user clicks on \"done\" button in hint popup", () => clickOn(page, el("\"done\" button in hint popup")));
    await session.step(73, "Then the tutorial step \"Explore each viewer\" should be done", () => stepDone(page, "Explore each viewer"));
    await session.step(75, "And the \"Multivariate Analysis\" tutorial should be completed", () => tutorialCompleted(page, "Multivariate Analysis"));
    await session.step(76, "And the tutorial should have listed 7 steps", () => tutorialStepsListed(page, 7));
    await session.step(77, "And the tutorial progress should be 7 of 7", () => tutorialProgress(page, 7, 7));
    await session.step(78, "And no hint should be shown", () => noHintShown(page));
    await session.step(79, "And no errors should have been logged", () => noErrors(page));
  });
});
