/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/filters/card-state.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [filters.cp.chem-and-bio-filters, public#1473, public#1225, public#1216]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKeyIn, shouldBe, shouldHaveText, typeInto, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsRowMoleculeOf} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {openCardSettings, pickSearchType, setCutoff} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {noErrors, pickFromAreaContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu, scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A substructure card from a cell, after a reset and on a cloned view", () => {
  const session = feature(test, "features/filters/card-state.feature", import.meta.url);
  test("A substructure card from a cell, after a reset and on a cloned view", {tag: ["@journey", "@realizes:filters.cp.chem-and-bio-filters", "@realizes:public#1473", "@realizes:public#1225", "@realizes:public#1216"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await run.scenario("Use as filter adds one card holding the cell's molecule", async () => {
      await session.step(19, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(20, "And user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(21, "And user clicks on close of \"Structure\" filter card", () => clickOn(page, el("close of \"Structure\" filter card")));
      await session.step(22, "Then \"Structure\" filter card should be absent", () => shouldBe(page, el("\"Structure\" filter card"), "absent"));
      await session.step(23, "When user picks \"Current Value > Use as filter\" from the context menu of the \"cell 1 of Structure\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Use as filter", "cell 1 of Structure", el("grid")));
      await session.step(24, "Then there should be 1 visible \"Structure\" filter card", () => visibleCount(page, 1, el("\"Structure\" filter card")));
      await session.step(25, "And the \"type of Structure\" reading of filter panel should be \"Chem:substructureFilter\"", () => readingReads(page, "type of Structure", el("filter panel"), "Chem:substructureFilter"));
      await session.step(26, "And the \"structure of Structure\" reading of filter panel should be the molecule in row 1 of \"Structure\" column", () => readingIsRowMoleculeOf(page, "structure of Structure", el("filter panel"), 1, "Structure"));
      await session.step(27, "And 1 row should pass the filter", () => filterPasses(page, 1));
      await session.step(28, "And the \"filtering of Structure\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of Structure", el("filter panel"), "true"));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Use as filter on a Core cell fills the Core card and leaves the Structure card's molecule (public#1473, public#1225, public#1216)", async () => {
      await session.step(33, "When user scrolls the grid to the \"Core\" column", () => scrollGridTo(page, "Core"));
      await session.step(34, "And user picks \"Current Value > Use as filter\" from the context menu of the \"cell 1 of Core\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Use as filter", "cell 1 of Core", el("grid")));
      await session.step(35, "Then the \"filtering of Core\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of Core", el("filter panel"), "true"));
      await session.step(36, "And the \"structure of Core\" reading of filter panel should be the molecule in row 1 of \"Core\" column", () => readingIsRowMoleculeOf(page, "structure of Core", el("filter panel"), 1, "Core"));
      await session.step(37, "And the \"structure of Structure\" reading of filter panel should be the molecule in row 1 of \"Structure\" column", () => readingIsRowMoleculeOf(page, "structure of Structure", el("filter panel"), 1, "Structure"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The reset icon clears the card's structure and search type and keeps the card", async () => {
      await session.step(41, "When user picks \"Remove All\" from the viewer menu of filter panel", () => pickFromViewerMenu(page, "Remove All", el("filter panel")));
      await session.step(42, "And user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(43, "And user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(44, "Then there should be 1 visible \"Structure\" filter card", () => visibleCount(page, 1, el("\"Structure\" filter card")));
      await session.step(45, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(46, "When user clicks on \"Sketch\" text in \"Structure\" filter card", () => clickOn(page, el("\"Sketch\" text in \"Structure\" filter card")));
      await session.step(47, "And user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(48, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(49, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(50, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(51, "When user opens the settings of the \"Structure\" filter card", () => openCardSettings(page, "Structure"));
      await session.step(52, "And user picks search type \"Similar\" in the \"Structure\" filter card", () => pickSearchType(page, "Similar", "Structure"));
      await session.step(53, "Then the \"search type of Structure\" reading of filter panel should be \"Similar\"", () => readingReads(page, "search type of Structure", el("filter panel"), "Similar"));
      await session.step(54, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(55, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(56, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(57, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(58, "And \"Structure\" filter card should be visible", () => shouldBe(page, el("\"Structure\" filter card"), "visible"));
      await session.step(59, "And the \"type of Structure\" reading of filter panel should be \"Chem:substructureFilter\"", () => readingReads(page, "type of Structure", el("filter panel"), "Chem:substructureFilter"));
      await session.step(60, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(61, "And the \"search type of Structure\" reading of filter panel should be \"Contains\"", () => readingReads(page, "search type of Structure", el("filter panel"), "Contains"));
      await session.step(62, "And the \"filtering of Structure\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of Structure", el("filter panel"), "false"));
      await session.step(63, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(64, "And \"Sketch\" text in \"Structure\" filter card should be visible", () => shouldBe(page, el("\"Sketch\" text in \"Structure\" filter card"), "visible"));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Similar card at a non-default cutoff comes back the same in a cloned view", async () => {
      await session.step(68, "When user clicks on \"Sketch\" text in \"Structure\" filter card", () => clickOn(page, el("\"Sketch\" text in \"Structure\" filter card")));
      await session.step(69, "And user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(70, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(71, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(72, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(73, "When user opens the settings of the \"Structure\" filter card", () => openCardSettings(page, "Structure"));
      await session.step(74, "And user picks search type \"Similar\" in the \"Structure\" filter card", () => pickSearchType(page, "Similar", "Structure"));
      await session.step(75, "And user sets the similarity cutoff of the \"Structure\" filter card to 0.6", () => setCutoff(page, "Structure", 0.6));
      await session.step(76, "Then the \"search type of Structure\" reading of filter panel should be \"Similar\"", () => readingReads(page, "search type of Structure", el("filter panel"), "Similar"));
      await session.step(77, "And the \"similarity cutoff of Structure\" reading of filter panel should be 0.6", () => readingIs(page, "similarity cutoff of Structure", el("filter panel"), 0.6));
      await session.step(78, "And the \"fingerprint of Structure\" reading of filter panel should be \"Morgan\"", () => readingReads(page, "fingerprint of Structure", el("filter panel"), "Morgan"));
      await session.step(79, "When user picks \"View > Layout > Clone View\" from the top menu", () => pickFromTopMenu(page, "View > Layout > Clone View"));
      await session.step(80, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(81, "And the \"search type of Structure\" reading of filter panel should be \"Similar\"", () => readingReads(page, "search type of Structure", el("filter panel"), "Similar"));
      await session.step(82, "And the \"similarity cutoff of Structure\" reading of filter panel should be 0.6", () => readingIs(page, "similarity cutoff of Structure", el("filter panel"), 0.6));
      await session.step(83, "And the \"fingerprint of Structure\" reading of filter panel should be \"Morgan\"", () => readingReads(page, "fingerprint of Structure", el("filter panel"), "Morgan"));
      await session.step(84, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
