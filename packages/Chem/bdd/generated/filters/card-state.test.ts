/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/filters/card-state.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [filters.cp.chem-and-bio-filters]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openCardSettings, pickSearchType, setCutoff} from '../../bindings/filter-card.js';
import {readingIsRowMolecule} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKeyIn, shouldBe, shouldHaveText, typeInto, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterPassesAll, filterPassesFewer} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromAreaContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A substructure card from a cell, after a reset and on a cloned view", () => {
  const session = feature(test, "features/filters/card-state.feature", import.meta.url);
  test("A substructure card from a cell, after a reset and on a cloned view", {tag: ["@journey", "@realizes:filters.cp.chem-and-bio-filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await run.scenario("Use as filter adds one card holding the cell's molecule", async () => {
      await session.step(18, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(19, "And user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(20, "And user clicks on close of \"Structure\" filter card", () => clickOn(page, el("close of \"Structure\" filter card")));
      await session.step(21, "Then \"Structure\" filter card should be absent", () => shouldBe(page, el("\"Structure\" filter card"), "absent"));
      await session.step(22, "When user picks \"Current Value > Use as filter\" from the context menu of the \"cell 1 of Structure\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Use as filter", "cell 1 of Structure", el("grid")));
      await session.step(23, "Then there should be 1 visible \"Structure\" filter card", () => visibleCount(page, 1, el("\"Structure\" filter card")));
      await session.step(24, "And the \"type of Structure\" reading of filter panel should be \"Chem:substructureFilter\"", () => readingReads(page, "type of Structure", el("filter panel"), "Chem:substructureFilter"));
      await session.step(25, "And the \"structure of Structure\" reading of filter panel should be the molecule of row 1 of \"Structure\" column", () => readingIsRowMolecule(page, "structure of Structure", 1, "Structure"));
      await session.step(26, "And fewer than 100 rows should pass the filter", () => filterPassesFewer(page, 100));
      await session.step(27, "And the \"filtering of Structure\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of Structure", el("filter panel"), "true"));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The reset icon clears the card's structure and search type and keeps the card", async () => {
      await session.step(31, "When user picks \"Remove All\" from the viewer menu of filter panel", () => pickFromViewerMenu(page, "Remove All", el("filter panel")));
      await session.step(32, "And user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(33, "And user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(34, "Then there should be 1 visible \"Structure\" filter card", () => visibleCount(page, 1, el("\"Structure\" filter card")));
      await session.step(35, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(36, "When user clicks on \"Sketch\" text in \"Structure\" filter card", () => clickOn(page, el("\"Sketch\" text in \"Structure\" filter card")));
      await session.step(37, "And user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(38, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(39, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(40, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(41, "When user opens the settings of the \"Structure\" filter card", () => openCardSettings(page, "Structure"));
      await session.step(42, "And user picks search type \"Similar\" in the \"Structure\" filter card", () => pickSearchType(page, "Similar", "Structure"));
      await session.step(43, "Then the \"search type of Structure\" reading of filter panel should be \"Similar\"", () => readingReads(page, "search type of Structure", el("filter panel"), "Similar"));
      await session.step(44, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(45, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(46, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(47, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(48, "And \"Structure\" filter card should be visible", () => shouldBe(page, el("\"Structure\" filter card"), "visible"));
      await session.step(49, "And the \"type of Structure\" reading of filter panel should be \"Chem:substructureFilter\"", () => readingReads(page, "type of Structure", el("filter panel"), "Chem:substructureFilter"));
      await session.step(50, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(51, "And the \"search type of Structure\" reading of filter panel should be \"Contains\"", () => readingReads(page, "search type of Structure", el("filter panel"), "Contains"));
      await session.step(52, "And the \"filtering of Structure\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of Structure", el("filter panel"), "false"));
      await session.step(53, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(54, "And \"Sketch\" text in \"Structure\" filter card should be visible", () => shouldBe(page, el("\"Sketch\" text in \"Structure\" filter card"), "visible"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Similar card at a non-default cutoff comes back the same in a cloned view", async () => {
      await session.step(58, "When user clicks on \"Sketch\" text in \"Structure\" filter card", () => clickOn(page, el("\"Sketch\" text in \"Structure\" filter card")));
      await session.step(59, "And user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(60, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(61, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(62, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(63, "When user opens the settings of the \"Structure\" filter card", () => openCardSettings(page, "Structure"));
      await session.step(64, "And user picks search type \"Similar\" in the \"Structure\" filter card", () => pickSearchType(page, "Similar", "Structure"));
      await session.step(65, "And user sets the similarity cutoff of the \"Structure\" filter card to 0.6", () => setCutoff(page, "Structure", 0.6));
      await session.step(66, "Then the \"search type of Structure\" reading of filter panel should be \"Similar\"", () => readingReads(page, "search type of Structure", el("filter panel"), "Similar"));
      await session.step(67, "And the \"similarity cutoff of Structure\" reading of filter panel should be 0.6", () => readingIs(page, "similarity cutoff of Structure", el("filter panel"), 0.6));
      await session.step(68, "And the \"fingerprint of Structure\" reading of filter panel should be \"Morgan\"", () => readingReads(page, "fingerprint of Structure", el("filter panel"), "Morgan"));
      await session.step(69, "When user picks \"View > Layout > Clone View\" from the top menu", () => pickFromTopMenu(page, "View > Layout > Clone View"));
      await session.step(70, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(71, "And the \"search type of Structure\" reading of filter panel should be \"Similar\"", () => readingReads(page, "search type of Structure", el("filter panel"), "Similar"));
      await session.step(72, "And the \"similarity cutoff of Structure\" reading of filter panel should be 0.6", () => readingIs(page, "similarity cutoff of Structure", el("filter panel"), 0.6));
      await session.step(73, "And the \"fingerprint of Structure\" reading of filter panel should be \"Morgan\"", () => readingReads(page, "fingerprint of Structure", el("filter panel"), "Morgan"));
      await session.step(74, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
