/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/eda/filters.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.filters]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, stepDoneTimes, stepNotDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKeyIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterIsExactlyAnyOf, filterIsExactlyCategory, filterPasses, filterPassesAll, noneOfFiltered, onlyOfSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {elementHinted, noHintShown, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {enterCardBound, forgetSavedState, pickCardIndicatorMenu, pickPanelMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {areaPainted, clickArea, hoverArea, noErrors, readingNotAsRemembered, readingReads, rememberReading, repainted, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Filters tutorial", () => {
  const session = feature(test, "features/eda/filters.feature", import.meta.url);
  test("A learner completes the Filters tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(23, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(24, "And the \"recentViewerSettings\" user settings are put back at feature end", () => userSettingsPutBack(page, "recentViewerSettings"));
    await session.step(25, "And the \"Filters\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Filters"));
    await session.step(26, "And no saved filter state \"AGE: [40,60]\" is kept, now or when the feature ends", () => forgetSavedState(page, "AGE: [40,60]"));
    await session.step(27, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(30, "When user starts the \"Filters\" tutorial", () => startTutorial(page, "Filters"));
    await session.step(31, "Then the tutorial progress should be 1 of 17", () => tutorialProgress(page, 1, 17));
    await session.step(32, "When user clicks on scatter-plot icon in toolbox", () => clickOn(page, el("scatter-plot icon in toolbox")));
    await session.step(33, "Then the tutorial step \"Open scatter plot\" should be done", () => stepDone(page, "Open scatter plot"));
    await session.step(34, "When user clicks on histogram icon in toolbox", () => clickOn(page, el("histogram icon in toolbox")));
    await session.step(35, "Then the tutorial step \"Open histogram\" should be done", () => stepDone(page, "Open histogram"));
    await session.step(36, "Then filter icon in toolbox should be hinted", () => elementHinted(page, el("filter icon in toolbox")));
    await session.step(37, "When user clicks on filter icon in toolbox", () => clickOn(page, el("filter icon in toolbox")));
    await session.step(38, "Then the tutorial step \"Open filters\" should be done", () => stepDone(page, "Open filters"));
    await session.step(39, "And filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
    await session.step(41, "When user clicks on the \"category AS of DIS_POP\" area of filter panel", () => clickArea(page, "category AS of DIS_POP", el("filter panel")));
    await session.step(42, "Then the tutorial step \"Click on the \\\"AS\\\" label within the \\\"DIS_POP\\\" filter\" should be done", () => stepDone(page, "Click on the \"AS\" label within the \"DIS_POP\" filter"));
    await session.step(43, "And the filter should pass exactly the rows where \"DIS_POP\" is \"AS\"", () => filterIsExactlyCategory(page, "DIS_POP", "AS"));
    await session.step(44, "When user presses ArrowDown in \"DIS_POP\" filter card", () => pressKeyIn(page, "ArrowDown", el("\"DIS_POP\" filter card")));
    await session.step(45, "Then the \"selected categories of DIS_POP\" reading of filter panel should be \"Indigestion\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "Indigestion"));
    await session.step(46, "When user presses ArrowDown in \"DIS_POP\" filter card", () => pressKeyIn(page, "ArrowDown", el("\"DIS_POP\" filter card")));
    await session.step(47, "Then the \"selected categories of DIS_POP\" reading of filter panel should be \"PsA\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "PsA"));
    await session.step(48, "When user presses ArrowDown in \"DIS_POP\" filter card", () => pressKeyIn(page, "ArrowDown", el("\"DIS_POP\" filter card")));
    await session.step(49, "Then the \"selected categories of DIS_POP\" reading of filter panel should be \"Psoriasis\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "Psoriasis"));
    await session.step(50, "When user presses ArrowDown in \"DIS_POP\" filter card", () => pressKeyIn(page, "ArrowDown", el("\"DIS_POP\" filter card")));
    await session.step(51, "Then the \"selected categories of DIS_POP\" reading of filter panel should be \"RA\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "RA"));
    await session.step(52, "And the tutorial step \"Filter the dataset by the most common disease\" should be done", () => stepDone(page, "Filter the dataset by the most common disease"));
    await session.step(53, "And 2550 rows should pass the filter", () => filterPasses(page, 2550));
    await session.step(55, "When user picks \"Invert all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Invert all", "DIS_POP"));
    await session.step(56, "Then the tutorial step \"Invert the category set for the \\\"DIS_POP\\\" filter\" should be done", () => stepDone(page, "Invert the category set for the \"DIS_POP\" filter"));
    await session.step(57, "And the filter should pass exactly the rows where \"DIS_POP\" is one of \"AS, Indigestion, PsA, Psoriasis, UC\"", () => filterIsExactlyAnyOf(page, "DIS_POP", "AS, Indigestion, PsA, Psoriasis, UC"));
    await session.step(59, "When user clicks on the \"count AS of DIS_POP\" area of filter panel", () => clickArea(page, "count AS of DIS_POP", el("filter panel")));
    await session.step(60, "Then the tutorial step \"Click on a non-empty row count\" should be done", () => stepDone(page, "Click on a non-empty row count"));
    await session.step(61, "And only rows where \"DIS_POP\" is \"AS\" should be selected", () => onlyOfSelected(page, "DIS_POP", "AS"));
    await session.step(63, "When user clicks on the \"category F of SEX\" area of filter panel", () => clickArea(page, "category F of SEX", el("filter panel")));
    await session.step(64, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
    await session.step(65, "And user clicks on the \"checkbox Black of RACE\" area of filter panel", () => clickArea(page, "checkbox Black of RACE", el("filter panel")));
    await session.step(66, "Then the tutorial step \"Filter the dataset to only females of Asian or Black origin\" should be done", () => stepDone(page, "Filter the dataset to only females of Asian or Black origin"));
    await session.step(67, "And the \"selected categories of SEX\" reading of filter panel should be \"F\"", () => readingReads(page, "selected categories of SEX", el("filter panel"), "F"));
    await session.step(68, "And the \"selected categories of RACE\" reading of filter panel should be \"Asian, Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Asian, Black"));
    await session.step(71, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
    await session.step(72, "Then \"Reset filter\" icon in filter panel should be hinted", () => elementHinted(page, el("\"Reset filter\" icon in filter panel")));
    await session.step(73, "When user clicks on \"Reset filter\" icon in filter panel", () => clickOn(page, el("\"Reset filter\" icon in filter panel")));
    await session.step(74, "Then the tutorial step \"Reset the filter\" should be done", () => stepDone(page, "Reset the filter"));
    await session.step(75, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(77, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
    await session.step(78, "And user hovers over the \"bin 3\" area of histogram viewer in filter panel", () => hoverArea(page, "bin 3", el("histogram viewer in filter panel")));
    await session.step(79, "Then the tutorial step \"Hover over the histogram bins\" should be done", () => stepDone(page, "Hover over the histogram bins"));
    await session.step(81, "And scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
    await session.step(83, "And the tutorial step \"Select one of the histogram bins\" should not be done yet", () => stepNotDone(page, "Select one of the histogram bins"));
    await session.step(84, "When user remembers the \"rows selected\" reading of scatter plot viewer", () => rememberReading(page, "rows selected", el("scatter plot viewer")));
    await session.step(85, "And user clicks on the \"bin 3\" area of histogram viewer in filter panel", () => clickArea(page, "bin 3", el("histogram viewer in filter panel")));
    await session.step(86, "Then the tutorial step \"Select one of the histogram bins\" should be done", () => stepDone(page, "Select one of the histogram bins"));
    await session.step(87, "And the \"rows selected\" reading of scatter plot viewer should not be as remembered", () => readingNotAsRemembered(page, "rows selected", el("scatter plot viewer")));
    await session.step(88, "And the \"selected bin 3\" area of histogram viewer in filter panel should be painted", () => areaPainted(page, "selected bin 3", el("histogram viewer in filter panel")));
    await session.step(90, "When user clicks on the \"cell 4 of AGE\" area of grid", () => clickArea(page, "cell 4 of AGE", el("grid")));
    await session.step(91, "Then the tutorial step \"Change the current row in the spreadsheet\" should be done", () => stepDone(page, "Change the current row in the spreadsheet"));
    await session.step(92, "And row 4 should be current", () => currentRowIs(page, 4));
    await session.step(94, "When user picks \"Min / max\" from the indicator menu of the \"AGE\" filter card", () => pickCardIndicatorMenu(page, "Min / max", "AGE"));
    await session.step(95, "And user enters \"40\" into the min field of the \"AGE\" filter card", () => enterCardBound(page, "40", "min", "AGE"));
    await session.step(96, "And user enters \"60\" into the max field of the \"AGE\" filter card", () => enterCardBound(page, "60", "max", "AGE"));
    await session.step(97, "Then the tutorial step \"Find records for people aged 40 to 60\" should be done", () => stepDone(page, "Find records for people aged 40 to 60"));
    await session.step(99, "And 2990 rows should pass the filter", () => filterPasses(page, 2990));
    await session.step(100, "And no rows where \"AGE\" is \"39\" should pass the filter", () => noneOfFiltered(page, "AGE", "39"));
    await session.step(101, "And no rows where \"AGE\" is \"61\" should pass the filter", () => noneOfFiltered(page, "AGE", "61"));
    await session.step(103, "When user picks \"Filter to Column...\" from the filter panel menu", () => pickPanelMenu(page, "Filter to Column..."));
    await session.step(104, "Then the tutorial step \"Save the current filter as a column\" should be done", () => stepDone(page, "Save the current filter as a column"));
    await session.step(105, "When user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
    await session.step(106, "Then the table should have a column \"AGE: [40,60]\"", () => hasColumn(page, "AGE: [40,60]"));
    await session.step(108, "When user picks \"Save or Apply > Save...\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply > Save..."));
    await session.step(109, "Then the tutorial step \"Save the filter configuration as \\\"AGE: [40,60]\\\"\" should be done", () => stepDone(page, "Save the filter configuration as \"AGE: [40,60]\""));
    await session.step(110, "When user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
    await session.step(112, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
    await session.step(113, "And user clicks on \"Reset filter\" icon in filter panel", () => clickOn(page, el("\"Reset filter\" icon in filter panel")));
    await session.step(114, "Then the tutorial step \"Reset the filter\" should be done 2 times", () => stepDoneTimes(page, "Reset the filter", 2));
    await session.step(115, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(116, "When user picks \"Save or Apply > AGE: [40,60]\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply > AGE: [40,60]"));
    await session.step(117, "Then the tutorial step \"Restore the filter state\" should be done", () => stepDone(page, "Restore the filter state"));
    await session.step(118, "And 2990 rows should pass the filter", () => filterPasses(page, 2990));
    await session.step(119, "And no rows where \"AGE\" is \"61\" should pass the filter", () => noneOfFiltered(page, "AGE", "61"));
    await session.step(121, "And the \"Filters\" tutorial should be completed", () => tutorialCompleted(page, "Filters"));
    await session.step(122, "And the tutorial should have listed 17 steps", () => tutorialStepsListed(page, 17));
    await session.step(123, "And the tutorial progress should be 17 of 17", () => tutorialProgress(page, 17, 17));
    await session.step(124, "And no hint should be shown", () => noHintShown(page));
    await session.step(125, "And no errors should have been logged", () => noErrors(page));
  });
});
