/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/calculated-columns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.calculated-columns]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, pressKey, pressKeyIn, replaceCode, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {displayedInRow, everyValueBetween, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {contextPanelShows, elementHinted, noHintShown, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, dragAreaBy, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Calculated Columns tutorial", () => {
  const session = feature(test, "features/transform/calculated-columns.feature", import.meta.url);
  test("A learner completes the Calculated Columns tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.calculated-columns"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(22, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(23, "And the \"Calculated Columns\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Calculated Columns"));
    await session.step(24, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(27, "When user starts the \"Calculated Columns\" tutorial", () => startTutorial(page, "Calculated Columns"));
    await session.step(28, "Then the tutorial progress should be 1 of 13", () => tutorialProgress(page, 1, 13));
    await session.step(29, "And add-new-column icon should be hinted", () => elementHinted(page, el("add-new-column icon")));
    await session.step(30, "When user clicks on add-new-column icon", () => clickOn(page, el("add-new-column icon")));
    await session.step(31, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
    await session.step(32, "And the tutorial step \"Open the \\\"Add New Column\\\" dialog\" should be done", () => stepDone(page, "Open the \"Add New Column\" dialog"));
    await session.step(34, "When user enters \"Height, m\" into Name input in \"Add New Column\" dialog", () => enterInto(page, "Height, m", el("Name input in \"Add New Column\" dialog")));
    await session.step(35, "Then the tutorial step \"Name a column \\\"Height, m\\\"\" should be done", () => stepDone(page, "Name a column \"Height, m\""));
    await session.step(36, "When user replaces the code of code editor in \"Add New Column\" dialog with \"Div(170, 100)\"", () => replaceCode(page, el("code editor in \"Add New Column\" dialog"), "Div(170, 100)"));
    await session.step(37, "Then the tutorial step \"Enter the expression \\\"Div(170, 100)\\\"\" should be done", () => stepDone(page, "Enter the expression \"Div(170, 100)\""));
    await session.step(38, "When user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
    await session.step(39, "Then the tutorial step \"Click \\\"OK\\\"\" should be done", () => stepDone(page, "Click \"OK\""));
    await session.step(40, "And the table should have a column \"Height, m\"", () => hasColumn(page, "Height, m"));
    await session.step(41, "And every value of \"Height, m\" column should lie between 1.699 and 1.701", () => everyValueBetween(page, "Height, m", 1.699, 1.701));
    await session.step(44, "When user drags the \"x scroll handle\" area of grid by 2000 pixels to the right", () => dragAreaBy(page, "x scroll handle", el("grid"), 2000, "right"));
    await session.step(45, "And user clicks on the \"header Height, m\" area of grid", () => clickArea(page, "header Height, m", el("grid")));
    await session.step(46, "Then the tutorial step \"Click on the \\\"Height, m\\\" column header\" should be done", () => stepDone(page, "Click on the \"Height, m\" column header"));
    await session.step(47, "And the context panel should show \"Height, m\"", () => contextPanelShows(page, "Height, m"));
    await session.step(48, "And \"Edit in dialog\" button in context panel should be hinted", () => elementHinted(page, el("\"Edit in dialog\" button in context panel")));
    await session.step(49, "When user clicks on \"Edit in dialog\" button in context panel", () => clickOn(page, el("\"Edit in dialog\" button in context panel")));
    await session.step(50, "Then \"Edit Column Formula\" dialog should be visible", () => shouldBe(page, el("\"Edit Column Formula\" dialog"), "visible"));
    await session.step(51, "And the tutorial step \"Click the \\\"Edit in dialog\\\" button under the formula field in the context panel\" should be done", () => stepDone(page, "Click the \"Edit in dialog\" button under the formula field in the context panel"));
    await session.step(53, "When user replaces the code of code editor in \"Edit Column Formula\" dialog with \"Div(${HEIGHT}, 100)\"", () => replaceCode(page, el("code editor in \"Edit Column Formula\" dialog"), "Div(${HEIGHT}, 100)"));
    await session.step(54, "And user clicks on OK button in \"Edit Column Formula\" dialog", () => clickOn(page, el("OK button in \"Edit Column Formula\" dialog")));
    await session.step(55, "Then the tutorial step \"Edit the formula to use the \\\"HEIGHT\\\" column values and click \\\"OK\\\"\" should be done", () => stepDone(page, "Edit the formula to use the \"HEIGHT\" column values and click \"OK\""));
    await session.step(56, "And the \"Height, m\" cell of row 2 should be displayed as \"1.636\"", () => displayedInRow(page, "Height, m", 2, "1.636"));
    await session.step(58, "When user double-clicks on the \"cell 1 of HEIGHT\" area of grid", () => doubleClickArea(page, "cell 1 of HEIGHT", el("grid")));
    await session.step(59, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
    await session.step(60, "And user types \"170\" into cell editor", () => typeInto(page, "170", el("cell editor")));
    await session.step(61, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(62, "Then the tutorial step \"Change the \\\"HEIGHT\\\" value in the first row to \\\"170\\\"\" should be done", () => stepDone(page, "Change the \"HEIGHT\" value in the first row to \"170\""));
    await session.step(63, "And the \"HEIGHT\" cell of row 1 should be displayed as \"170.000\"", () => displayedInRow(page, "HEIGHT", 1, "170.000"));
    await session.step(65, "And the \"Height, m\" cell of row 1 should be displayed as \"1.605\"", () => displayedInRow(page, "Height, m", 1, "1.605"));
    await session.step(67, "When user clicks on add-new-column icon", () => clickOn(page, el("add-new-column icon")));
    await session.step(68, "Then the tutorial step \"Add a new column that calculates BMI\" should be done", () => stepDone(page, "Add a new column that calculates BMI"));
    await session.step(69, "When user enters \"BMI\" into Name input in \"Add New Column\" dialog", () => enterInto(page, "BMI", el("Name input in \"Add New Column\" dialog")));
    await session.step(70, "Then the tutorial step \"Name a column \\\"BMI\\\"\" should be done", () => stepDone(page, "Name a column \"BMI\""));
    await session.step(71, "When user replaces the code of code editor in \"Add New Column\" dialog with \"Div(${WEIGHT}, Pow(${Height, m}, 2))\"", () => replaceCode(page, el("code editor in \"Add New Column\" dialog"), "Div(${WEIGHT}, Pow(${Height, m}, 2))"));
    await session.step(72, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
    await session.step(73, "Then the tutorial step \"Enter the BMI formula and click \\\"OK\\\"\" should be done", () => stepDone(page, "Enter the BMI formula and click \"OK\""));
    await session.step(74, "And the table should have a column \"BMI\"", () => hasColumn(page, "BMI"));
    await session.step(76, "And the \"BMI\" cell of row 2 should be displayed as \"34.73\"", () => displayedInRow(page, "BMI", 2, "34.73"));
    await session.step(78, "When user clicks on the \"header Height, m\" area of grid", () => clickArea(page, "header Height, m", el("grid")));
    await session.step(79, "Then the context panel should show \"Height, m\"", () => contextPanelShows(page, "Height, m"));
    await session.step(81, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
    await session.step(82, "When user replaces the code of code editor in Formula pane with \"RoundFloat(Div(${HEIGHT}, 100), 2)\"", () => replaceCode(page, el("code editor in Formula pane"), "RoundFloat(Div(${HEIGHT}, 100), 2)"));
    await session.step(83, "And user clicks on \"Apply\" button in Formula pane", () => clickOn(page, el("\"Apply\" button in Formula pane")));
    await session.step(84, "Then the tutorial step \"Update the formula for \\\"Height, m\\\" to round the values to 2 decimal places\" should be done", () => stepDone(page, "Update the formula for \"Height, m\" to round the values to 2 decimal places"));
    await session.step(85, "And the \"Height, m\" cell of row 2 should be displayed as \"1.640\"", () => displayedInRow(page, "Height, m", 2, "1.640"));
    await session.step(87, "And the \"BMI\" cell of row 2 should be displayed as \"34.58\"", () => displayedInRow(page, "BMI", 2, "34.58"));
    await session.step(89, "And the \"Calculated Columns\" tutorial should be completed", () => tutorialCompleted(page, "Calculated Columns"));
    await session.step(90, "And the tutorial should have listed 12 steps", () => tutorialStepsListed(page, 12));
    await session.step(91, "And the tutorial progress should be 13 of 13", () => tutorialProgress(page, 13, 13));
    await session.step(92, "And no hint should be shown", () => noHintShown(page));
    await session.step(93, "And no errors should have been logged", () => noErrors(page));
  });
});
