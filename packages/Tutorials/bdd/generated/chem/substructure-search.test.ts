/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/chem/substructure-search.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [tutorials.substructure-search-filtering]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {startTutorial, stepDone, tutorialCompleted, tutorialNotCompleted, tutorialProgress, tutorialStepsListed, tutorialsOpen} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, shouldBe, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPanelHas, filterPasses, filterPassesAll, selectedPassFilter, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, noHintShown, sketcherIs, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {enterCardBound, openCardSettings, pickCardIndicatorMenu, pickSearchType} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, hoverArea, noErrors, pickFromAreaContextMenu, readingNotAsRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Substructure Search and Filtering tutorial", () => {
  const session = feature(test, "features/chem/substructure-search.feature", import.meta.url);
  test("A learner completes the Substructure Search and Filtering tutorial", {tag: ["@tutorials", "@serial", "@realizes:tutorials.substructure-search-filtering"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(23, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(24, "And the \"tutorials\" user settings are put back at feature end", () => userSettingsPutBack(page, "tutorials"));
    await session.step(25, "And the \"achievement-badges\" user settings are put back at feature end", () => userSettingsPutBack(page, "achievement-badges"));
    await session.step(26, "And the \"Substructure Search and Filtering\" tutorial is not completed yet", () => tutorialNotCompleted(page, "Substructure Search and Filtering"));
    await session.step(27, "And the Tutorials app is open", () => tutorialsOpen(page));
    await session.step(30, "When user starts the \"Substructure Search and Filtering\" tutorial", () => startTutorial(page, "Substructure Search and Filtering"));
    await session.step(31, "Then the tutorial progress should be 1 of 11", () => tutorialProgress(page, 1, 11));
    await session.step(32, "When user picks \"Chem > Search > Substructure Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Substructure Search..."));
    await session.step(33, "Then the tutorial step \"Click Chem > Search > Substructure Search…\" should be done", () => stepDone(page, "Click Chem > Search > Substructure Search…"));
    await session.step(34, "And sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(37, "When user types \"c1ccc2ccccc2c1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccc2ccccc2c1", el("molecule input of sketcher dialog")));
    await session.step(38, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(39, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(40, "Then the tutorial step \"Change substructure\" should be done", () => stepDone(page, "Change substructure"));
    await session.step(41, "And 186 rows should pass the filter", () => filterPasses(page, 186));
    await session.step(42, "And the \"structure of smiles\" reading of filter panel should be \"c1ccc2ccccc2c1\"", () => readingReads(page, "structure of smiles", el("filter panel"), "c1ccc2ccccc2c1"));
    await session.step(45, "When user hovers over the \"card smiles\" area of filter panel", () => hoverArea(page, "card smiles", el("filter panel")));
    await session.step(46, "And user clicks on \"Clear\" button in filter panel", () => clickOn(page, el("\"Clear\" button in filter panel")));
    await session.step(47, "Then the tutorial step \"In the filter panel, click CLEAR\" should be done", () => stepDone(page, "In the filter panel, click CLEAR"));
    await session.step(48, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(50, "When user picks \"Current Value > Use as filter\" from the context menu of the \"cell 3 of smiles\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Use as filter", "cell 3 of smiles", el("grid")));
    await session.step(51, "Then the tutorial step \"In the grid, right-click any molecule and select Current Value > Use as filter\" should be done", () => stepDone(page, "In the grid, right-click any molecule and select Current Value > Use as filter"));
    await session.step(53, "And 1 row should pass the filter", () => filterPasses(page, 1));
    await session.step(55, "When user clicks on the \"card smiles\" area of filter panel", () => clickArea(page, "card smiles", el("filter panel")));
    await session.step(56, "Then the tutorial step \"On the Filter Panel, click the molecule and modify it in the sketcher\" should be done", () => stepDone(page, "On the Filter Panel, click the molecule and modify it in the sketcher"));
    await session.step(57, "And sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(58, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(59, "Then the tutorial step \"Click OK\" should be done", () => stepDone(page, "Click OK"));
    await session.step(61, "When user remembers the \"rows shown\" reading of filter panel", () => rememberReading(page, "rows shown", el("filter panel")));
    await session.step(62, "And user opens the settings of the \"smiles\" filter card", () => openCardSettings(page, "smiles"));
    await session.step(63, "And user picks search type \"Not contains\" in the \"smiles\" filter card", () => pickSearchType(page, "Not contains", "smiles"));
    await session.step(64, "Then the tutorial step \"Exclude the specified substructure from the view\" should be done", () => stepDone(page, "Exclude the specified substructure from the view"));
    await session.step(65, "And the \"search type of smiles\" reading of filter panel should be \"Not contains\"", () => readingReads(page, "search type of smiles", el("filter panel"), "Not contains"));
    await session.step(66, "And the \"rows shown\" reading of filter panel should not be as remembered", () => readingNotAsRemembered(page, "rows shown", el("filter panel")));
    await session.step(68, "When user unchecks checkbox of \"smiles\" filter card", () => uncheck(page, el("checkbox of \"smiles\" filter card")));
    await session.step(69, "Then the tutorial step \"In the filter panel, turn off the filter by clearing the checkbox\" should be done", () => stepDone(page, "In the filter panel, turn off the filter by clearing the checkbox"));
    await session.step(70, "And the \"enabled of smiles\" reading of filter panel should be \"false\"", () => readingReads(page, "enabled of smiles", el("filter panel"), "false"));
    await session.step(71, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(73, "When user picks \"Select Columns...\" from the viewer menu of filter panel", () => pickFromViewerMenu(page, "Select Columns...", el("filter panel")));
    await session.step(74, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(75, "When user types \"NOCount\" into \"Search\" input in \"Select columns...\" dialog", () => typeInto(page, "NOCount", el("\"Search\" input in \"Select columns...\" dialog")));
    await session.step(76, "And user toggles the \"NOCount\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "NOCount", el("\"Select columns...\" dialog")));
    await session.step(77, "And user types \"NumRotatableBonds\" into \"Search\" input in \"Select columns...\" dialog", () => typeInto(page, "NumRotatableBonds", el("\"Search\" input in \"Select columns...\" dialog")));
    await session.step(78, "And user toggles the \"NumRotatableBonds\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "NumRotatableBonds", el("\"Select columns...\" dialog")));
    await session.step(79, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(80, "Then the tutorial step \"Select columns to be used as filters\" should be done", () => stepDone(page, "Select columns to be used as filters"));
    await session.step(81, "And the filter panel should have a filter on \"NOCount\" column", () => filterPanelHas(page, "NOCount"));
    await session.step(82, "And the filter panel should have a filter on \"NumRotatableBonds\" column", () => filterPanelHas(page, "NumRotatableBonds"));
    await session.step(84, "When user picks \"Min / max\" from the indicator menu of the \"NOCount\" filter card", () => pickCardIndicatorMenu(page, "Min / max", "NOCount"));
    await session.step(85, "And user enters \"2\" into the min field of the \"NOCount\" filter card", () => enterCardBound(page, "2", "min", "NOCount"));
    await session.step(86, "And user enters \"4\" into the max field of the \"NOCount\" filter card", () => enterCardBound(page, "4", "max", "NOCount"));
    await session.step(87, "And user presses Control+A in grid", () => pressKeyIn(page, "Control+A", el("grid")));
    await session.step(88, "Then the tutorial step \"Interact with filters by changing their values. After that select rows of your interest\" should be done", () => stepDone(page, "Interact with filters by changing their values. After that select rows of your interest"));
    await session.step(90, "And 2499 rows should pass the filter", () => filterPasses(page, 2499));
    await session.step(91, "And 2499 rows should be selected", () => selectedRowCount(page, 2499));
    await session.step(92, "And every selected row should pass the filter", () => selectedPassFilter(page));
    await session.step(94, "And the \"Substructure Search and Filtering\" tutorial should be completed", () => tutorialCompleted(page, "Substructure Search and Filtering"));
    await session.step(95, "And the tutorial should have listed 10 steps", () => tutorialStepsListed(page, 10));
    await session.step(96, "And the tutorial progress should be 11 of 11", () => tutorialProgress(page, 11, 11));
    await session.step(97, "And no hint should be shown", () => noHintShown(page));
    await session.step(98, "And no errors should have been logged", () => noErrors(page));
  });
});
