/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/data-aggregation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.data-aggregation]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, noneOfFiltered, noneOfSelected, noneSelected, selectedPassFilter, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, noHintShown, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addToRow, clearSavedParameters, pickFromHistory} from '@datagrok-libraries/bdd/bindings/tiers/viewers/pivot-table';
import {clickArea, clickAreaHolding, closeContextMenu, noErrors, pickFromAreaContextMenu, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Data Aggregation tutorial", () => {
  const session = feature(test, "features/transform/data-aggregation.feature", import.meta.url);
  test("A learner completes the Data Aggregation tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.data-aggregation"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(22, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(23, "And the \"Data Aggregation\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Data Aggregation"));
    await session.step(24, "And user clears the saved pivot table parameters", () => clearSavedParameters(page));
    await session.step(25, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(28, "When user starts the \"Data Aggregation\" tutorial", () => startTutorial(page, "Data Aggregation"));
    await session.step(29, "Then the tutorial progress should be 1 of 11", () => tutorialProgress(page, 1, 11));
    await session.step(30, "When user picks \"Data > Aggregate Rows...\" from the top menu", () => pickFromTopMenu(page, "Data > Aggregate Rows..."));
    await session.step(31, "Then the tutorial step \"Open Aggregation Editor\" should be done", () => stepDone(page, "Open Aggregation Editor"));
    await session.step(32, "And the open tableview should have 1 pivot table viewer", () => viewerCount(page, 1, "pivot table"));
    await session.step(34, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE), avg(HEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE), avg(HEIGHT)"));
    await session.step(36, "When user adds \"RACE\" to the \"group by\" row of pivot table viewer", () => addToRow(page, "RACE", "group by"));
    await session.step(37, "Then the tutorial step \"Group rows by column \\\"RACE\\\"\" should be done", () => stepDone(page, "Group rows by column \"RACE\""));
    await session.step(38, "When user adds \"SEX\" to the \"group by\" row of pivot table viewer", () => addToRow(page, "SEX", "group by"));
    await session.step(39, "Then the tutorial step \"Group rows by column \\\"SEX\\\"\" should be done", () => stepDone(page, "Group rows by column \"SEX\""));
    await session.step(40, "And the \"group by\" reading of pivot table viewer should be \"RACE, SEX\"", () => readingReads(page, "group by", el("pivot table viewer"), "RACE, SEX"));
    await session.step(41, "When user adds \"DIS_POP\" to the \"pivot\" row of pivot table viewer", () => addToRow(page, "DIS_POP", "pivot"));
    await session.step(42, "Then the tutorial step \"Pivot data by column \\\"DIS_POP\\\"\" should be done", () => stepDone(page, "Pivot data by column \"DIS_POP\""));
    await session.step(43, "And the \"pivot\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "pivot", el("pivot table viewer"), "DIS_POP"));
    await session.step(45, "When user picks \"Remove others\" from the context menu of the \"aggregate chip avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Remove others", "aggregate chip avg(AGE)", el("pivot table viewer")));
    await session.step(46, "Then the tutorial step \"Leave only the \\\"avg(AGE)\\\" aggregation\" should be done", () => stepDone(page, "Leave only the \"avg(AGE)\" aggregation"));
    await session.step(47, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
    await session.step(48, "When user picks \"Column > WEIGHT\" from the context menu of the \"aggregate chip avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Column > WEIGHT", "aggregate chip avg(AGE)", el("pivot table viewer")));
    await session.step(49, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(50, "Then the tutorial step \"Change a column to \\\"WEIGHT\\\"\" should be done", () => stepDone(page, "Change a column to \"WEIGHT\""));
    await session.step(51, "And the \"aggregate\" reading of pivot table viewer should be \"avg(WEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(WEIGHT)"));
    await session.step(52, "When user picks \"Aggregation > med\" from the context menu of the \"aggregate chip avg(WEIGHT)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Aggregation > med", "aggregate chip avg(WEIGHT)", el("pivot table viewer")));
    await session.step(53, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(54, "Then the tutorial step \"Change the aggregation function to \\\"med\\\"\" should be done", () => stepDone(page, "Change the aggregation function to \"med\""));
    await session.step(55, "And the \"aggregate\" reading of pivot table viewer should be \"med(WEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "med(WEIGHT)"));
    await session.step(58, "When user clicks on the \"grid row header 1\" area of pivot table viewer holding Shift", () => clickAreaHolding(page, "grid row header 1", el("pivot table viewer"), "Shift"));
    await session.step(59, "Then the tutorial step \"Select rows in the source table with values of the first aggregated row\" should be done", () => stepDone(page, "Select rows in the source table with values of the first aggregated row"));
    await session.step(60, "And 37 rows should be selected", () => selectedRowCount(page, 37));
    await session.step(61, "And every selected row should pass the filter", () => selectedPassFilter(page));
    await session.step(62, "And no rows where \"RACE\" is \"Caucasian\" should be selected", () => noneOfSelected(page, "RACE", "Caucasian"));
    await session.step(63, "And no rows where \"RACE\" is \"Other\" should be selected", () => noneOfSelected(page, "RACE", "Other"));
    await session.step(64, "And no rows where \"SEX\" is \"M\" should be selected", () => noneOfSelected(page, "SEX", "M"));
    await session.step(65, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(66, "Then the tutorial step \"Remove selection by pressing \\\"Esc\\\"\" should be done", () => stepDone(page, "Remove selection by pressing \"Esc\""));
    await session.step(67, "And no rows should be selected", () => noneSelected(page));
    await session.step(70, "When user clicks on the \"grid cell 8 of RACE\" area of pivot table viewer", () => clickArea(page, "grid cell 8 of RACE", el("pivot table viewer")));
    await session.step(71, "Then the tutorial step \"Click on the last row in the aggregated table to filter by it\" should be done", () => stepDone(page, "Click on the last row in the aggregated table to filter by it"));
    await session.step(72, "And 75 rows should pass the filter", () => filterPasses(page, 75));
    await session.step(73, "And no rows where \"RACE\" is \"Asian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Asian"));
    await session.step(74, "And no rows where \"SEX\" is \"F\" should pass the filter", () => noneOfFiltered(page, "SEX", "F"));
    await session.step(76, "When user picks \"Save parameters\" from the history menu of pivot table viewer", () => pickFromHistory(page, "Save parameters"));
    await session.step(77, "Then the tutorial step \"Save parameters\" should be done", () => stepDone(page, "Save parameters"));
    await session.step(78, "And the \"history entries\" reading of pivot table viewer should contain \"med(WEIGHT)\"", () => readingContains(page, "history entries", el("pivot table viewer"), "med(WEIGHT)"));
    await session.step(80, "And the \"Data Aggregation\" tutorial should be completed", () => tutorialCompleted(page, "Data Aggregation"));
    await session.step(81, "And the tutorial should have listed 11 steps", () => tutorialStepsListed(page, 11));
    await session.step(82, "And the tutorial progress should be 11 of 11", () => tutorialProgress(page, 11, 11));
    await session.step(83, "And no hint should be shown", () => noHintShown(page));
    await session.step(84, "And no errors should have been logged", () => noErrors(page));
  });
});
