/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/chem/r-groups-analysis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.r-groups-analysis]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, finishedUpdating, hoverOver, pressKey, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectedRowCount, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelShows, contextPanelShowsCurrentCell, noHintShown, packageInstalled, sketcherIs, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, dragBetweenAreasHolding, noErrors, propertyShouldBe, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {innerPropertyShouldBe, pickInAreaSelector, pickInnerViewer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The R-Groups Analysis tutorial", () => {
  const session = feature(test, "features/chem/r-groups-analysis.feature", import.meta.url);
  test("A learner completes the R-Groups Analysis tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.r-groups-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(32, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(33, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(34, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(35, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(36, "And the \"R-Groups Analysis\" tutorial is not completed yet", () => tutorialNotCompleted(page, "R-Groups Analysis"));
    await session.step(37, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(40, "When user starts the \"R-Groups Analysis\" tutorial", () => startTutorial(page, "R-Groups Analysis"));
    await session.step(41, "Then the tutorial progress should be 1 of 15", () => tutorialProgress(page, 1, 15));
    await session.step(42, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
    await session.step(43, "Then the tutorial step \"On the Top Menu, click Chem > Analyze > R-Groups Analysis...\" should be done", () => stepDone(page, "On the Top Menu, click Chem > Analyze > R-Groups Analysis..."));
    await session.step(44, "And \"R-Groups Analysis\" dialog should be visible", () => shouldBe(page, el("\"R-Groups Analysis\" dialog"), "visible"));
    await session.step(45, "When user clicks on \"MCS\" button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("\"MCS\" button in \"R-Groups Analysis\" dialog")));
    await session.step(46, "Then the tutorial step \"Click MCS\" should be done", () => stepDone(page, "Click MCS"));
    await session.step(49, "And \"R-Groups Analysis\" dialog should have finished updating", () => finishedUpdating(page, el("\"R-Groups Analysis\" dialog")));
    await session.step(50, "When user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
    await session.step(51, "Then the tutorial step \"Click OK\" should be done", () => stepDone(page, "Click OK"));
    await session.step(52, "And the tutorial step \"Wait for the analysis to complete\" should be done", () => stepDone(page, "Wait for the analysis to complete"));
    await session.step(53, "And trellis plot viewer should be visible", () => shouldBe(page, el("trellis plot viewer"), "visible"));
    await session.step(54, "And the table should have a column \"R1\"", () => hasColumn(page, "R1"));
    await session.step(55, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Pie chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Pie chart"));
    await session.step(57, "When user hovers over trellis plot viewer", () => hoverOver(page, el("trellis plot viewer")));
    await session.step(58, "And user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
    await session.step(59, "Then the tutorial step \"In the trellis plot, click the gear icon for the embedded viewer\" should be done", () => stepDone(page, "In the trellis plot, click the gear icon for the embedded viewer"));
    await session.step(60, "And the context panel should show \"Trellis plot\"", () => contextPanelShows(page, "Trellis plot"));
    await session.step(61, "When user clicks on \"Pie chart\" tab in context panel", () => clickOn(page, el("\"Pie chart\" tab in context panel")));
    await session.step(62, "Then the tutorial step \"Go to Pie chart tab\" should be done", () => stepDone(page, "Go to Pie chart tab"));
    await session.step(63, "When user selects \"LC/MS\" in \"Category\" property", () => selectIn(page, "LC/MS", el("\"Category\" property")));
    await session.step(64, "Then the tutorial step \"Under Pie chart tab > Data, set Category to LC/MS\" should be done", () => stepDone(page, "Under Pie chart tab > Data, set Category to LC/MS"));
    await session.step(65, "And \"categoryColumnName\" inner property of trellis plot viewer should be \"LC/MS\"", () => innerPropertyShouldBe(page, "categoryColumnName", el("trellis plot viewer"), "LC/MS"));
    await session.step(67, "When user clicks on the \"cell body 1,1\" area of trellis plot viewer", () => clickArea(page, "cell body 1,1", el("trellis plot viewer")));
    await session.step(68, "Then the tutorial step \"Click any segment on a pie chart\" should be done", () => stepDone(page, "Click any segment on a pie chart"));
    await session.step(69, "And some rows should be selected", () => someSelected(page));
    await session.step(70, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(71, "Then the tutorial step \"Press Escape\" should be done", () => stepDone(page, "Press Escape"));
    await session.step(73, "When user clicks on the \"cell 1 of R1\" area of grid", () => clickArea(page, "cell 1 of R1", el("grid")));
    await session.step(74, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(75, "When user drags from the \"row header 1\" area to the \"row header 7\" area of grid holding Shift", () => dragBetweenAreasHolding(page, "row header 1", "row header 7", el("grid"), "Shift"));
    await session.step(76, "Then the tutorial step \"In the grid, press Shift+Drag Mouse Down\" should be done", () => stepDone(page, "In the grid, press Shift+Drag Mouse Down"));
    await session.step(77, "And 7 rows should be selected", () => selectedRowCount(page, 7));
    await session.step(78, "When user expands \"Distributions\" pane in context panel", () => expand(page, el("\"Distributions\" pane in context panel")));
    await session.step(79, "Then the tutorial step \"On the Context Panel, expand the Distributions pane\" should be done", () => stepDone(page, "On the Context Panel, expand the Distributions pane"));
    await session.step(80, "When user hovers over \"Distributions\" pane in context panel", () => hoverOver(page, el("\"Distributions\" pane in context panel")));
    await session.step(81, "Then the tutorial step \"In the pane, hover over line charts to see distributions\" should be done", () => stepDone(page, "In the pane, hover over line charts to see distributions"));
    await session.step(82, "And tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(84, "When user picks \"Histogram\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Histogram", el("trellis plot viewer")));
    await session.step(85, "Then the tutorial step \"In the top-left corner of the trellis plot, select Histogram.\" should be done", () => stepDone(page, "In the top-left corner of the trellis plot, select Histogram."));
    await session.step(86, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Histogram\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Histogram"));
    await session.step(87, "When user clicks on the \"cell 8 of R1\" area of grid", () => clickArea(page, "cell 8 of R1", el("grid")));
    await session.step(88, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(89, "When user hovers over trellis plot viewer", () => hoverOver(page, el("trellis plot viewer")));
    await session.step(90, "And user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
    await session.step(91, "And user clicks on \"Histogram\" tab in context panel", () => clickOn(page, el("\"Histogram\" tab in context panel")));
    await session.step(92, "And user selects \"In-vivo Activity\" in \"Value\" property", () => selectIn(page, "In-vivo Activity", el("\"Value\" property")));
    await session.step(93, "Then the tutorial step \"Set Value to In-Vivo Activity\" should be done", () => stepDone(page, "Set Value to In-Vivo Activity"));
    await session.step(94, "And \"valueColumnName\" inner property of trellis plot viewer should be \"In-vivo Activity\"", () => innerPropertyShouldBe(page, "valueColumnName", el("trellis plot viewer"), "In-vivo Activity"));
    await session.step(96, "When user picks \"R4\" in the column selector at the \"x selector 1\" area of trellis plot viewer", () => pickInAreaSelector(page, "R4", "x selector 1", el("trellis plot viewer")));
    await session.step(97, "Then the tutorial step \"Set the value for the X axis to R4\" should be done", () => stepDone(page, "Set the value for the X axis to R4"));
    await session.step(98, "And \"X Column Names\" property of trellis plot viewer should be \"R4\"", () => propertyShouldBe(page, "X Column Names", el("trellis plot viewer"), "R4"));
    await session.step(100, "And the \"R-Groups Analysis\" tutorial should be completed", () => tutorialCompleted(page, "R-Groups Analysis"));
    await session.step(101, "And the tutorial should have listed 15 steps", () => tutorialStepsListed(page, 15));
    await session.step(102, "And the tutorial progress should be 15 of 15", () => tutorialProgress(page, 15, 15));
    await session.step(103, "And no hint should be shown", () => noHintShown(page));
    await session.step(104, "And no errors should have been logged", () => noErrors(page));
  });
});
