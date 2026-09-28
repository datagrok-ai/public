/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/ml/scripting.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.scripting]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, stepDoneTimes, stepListedTimes, stepNotDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, typeIntoEditorLine, typeIntoFirstEmptyLine} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, noHintShown, standRunsService, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Scripting tutorial", () => {
  const session = feature(test, "features/ml/scripting.feature", import.meta.url);
  test("A learner completes the Scripting tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.scripting"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(22, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(23, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(24, "And the \"Scripting\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Scripting"));
    await session.step(25, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(29, "Given the stand runs the \"Jupyter\" service", () => standRunsService(page, "Jupyter"));
    await session.step(30, "When user starts the \"Scripting\" tutorial", () => startTutorial(page, "Scripting"));
    await session.step(31, "Then the tutorial progress should be 1 of 10", () => tutorialProgress(page, 1, 10));
    await session.step(32, "Given the tutorial step \"In the Browse Panel, click Platform > Functions > Scripts > New > Python Script. This opens a script editor.\" should not be done yet", () => stepNotDone(page, "In the Browse Panel, click Platform > Functions > Scripts > New > Python Script. This opens a script editor."));
    await session.step(33, "When user clicks on Platform---Functions---Scripts tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Scripts tree node inside browse tree")));
    await session.step(34, "And user clicks on NEW button", () => clickOn(page, el("NEW button")));
    await session.step(35, "And user picks \"Python Script\" from the open menu", () => pickFromOpenMenu(page, "Python Script"));
    await session.step(36, "Then the tutorial step \"In the Browse Panel, click Platform > Functions > Scripts > New > Python Script. This opens a script editor.\" should be done", () => stepDone(page, "In the Browse Panel, click Platform > Functions > Scripts > New > Python Script. This opens a script editor."));
    await session.step(37, "And code editor should be visible", () => shouldBe(page, el("code editor"), "visible"));
    await session.step(38, "Given the tutorial step \"Open a sample table for the script\" should not be done yet", () => stepNotDone(page, "Open a sample table for the script"));
    await session.step(39, "When user clicks on asterisk icon", () => clickOn(page, el("asterisk icon")));
    await session.step(40, "Then the tutorial step \"Open a sample table for the script\" should be done", () => stepDone(page, "Open a sample table for the script"));
    await session.step(41, "Given the tutorial step \"Run the script\" should not be done yet", () => stepNotDone(page, "Run the script"));
    await session.step(42, "When user clicks on play icon", () => clickOn(page, el("play icon")));
    await session.step(43, "Then the tutorial step \"Run the script\" should be done", () => stepDone(page, "Run the script"));
    await session.step(44, "And \"Template\" dialog should be visible", () => shouldBe(page, el("\"Template\" dialog"), "visible"));
    await session.step(45, "When user selects \"cars\" in \"Table\" input in \"Template\" dialog", () => selectIn(page, "cars", el("\"Table\" input in \"Template\" dialog")));
    await session.step(46, "Then the tutorial step \"Set \\\"Table\\\" to cars\" should be done", () => stepDone(page, "Set \"Table\" to cars"));
    await session.step(47, "When user clicks on OK button in \"Template\" dialog", () => clickOn(page, el("OK button in \"Template\" dialog")));
    await session.step(48, "Then the tutorial step \"Click \\\"OK\\\"\" should be done", () => stepDone(page, "Click \"OK\""));
    await session.step(49, "Given the tutorial step \"Find the results in the console\" should not be done yet", () => stepNotDone(page, "Find the results in the console"));
    await session.step(50, "When user clicks on \"Console\" status bar toggle", () => clickOn(page, el("\"Console\" status bar toggle")));
    await session.step(51, "Then the tutorial step \"Find the results in the console\" should be done", () => stepDone(page, "Find the results in the console"));
    await session.step(53, "Given the tutorial step \"Add the second output value to the script\" should not be done yet", () => stepNotDone(page, "Add the second output value to the script"));
    await session.step(54, "When user types \"#output: dataframe clone\" into the first empty line of code editor", () => typeIntoFirstEmptyLine(page, "#output: dataframe clone", el("code editor")));
    await session.step(55, "And user types \"clone = table\" into the last line of code editor", () => typeIntoEditorLine(page, "clone = table", "last", el("code editor")));
    await session.step(56, "Then the tutorial step \"Add the second output value to the script\" should be done", () => stepDone(page, "Add the second output value to the script"));
    await session.step(57, "Given the tutorial step \"Run the script\" should be listed 2 times", () => stepListedTimes(page, "Run the script", 2));
    await session.step(58, "When user clicks on play icon", () => clickOn(page, el("play icon")));
    await session.step(59, "Then the tutorial step \"Run the script\" should be done 2 times", () => stepDoneTimes(page, "Run the script", 2));
    await session.step(60, "When user selects \"cars\" in \"Table\" input in \"Template\" dialog", () => selectIn(page, "cars", el("\"Table\" input in \"Template\" dialog")));
    await session.step(61, "Then the tutorial step \"Set \\\"Table\\\" to cars\" should be done 2 times", () => stepDoneTimes(page, "Set \"Table\" to cars", 2));
    await session.step(62, "When user clicks on OK button in \"Template\" dialog", () => clickOn(page, el("OK button in \"Template\" dialog")));
    await session.step(63, "Then the tutorial step \"Click \\\"OK\\\"\" should be done 2 times", () => stepDoneTimes(page, "Click \"OK\"", 2));
    await session.step(65, "And the \"Scripting\" tutorial should be completed", () => tutorialCompleted(page, "Scripting"));
    await session.step(66, "And the tutorial should have listed 10 steps", () => tutorialStepsListed(page, 10));
    await session.step(67, "And the tutorial progress should be 10 of 10", () => tutorialProgress(page, 10, 10));
    await session.step(68, "And no hint should be shown", () => noHintShown(page));
    await session.step(69, "And no errors should have been logged", () => noErrors(page));
  });
});
