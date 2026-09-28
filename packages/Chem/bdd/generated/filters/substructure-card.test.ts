/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/filters/substructure-card.feature
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
import {offersSearchTypes, openCardSettings, pickSearchType} from '../../bindings/filter-card.js';
import {filterPassesMatching} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKeyIn, shouldBe, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {panelHasNoCardOfType} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The substructure filter card and its search types", () => {
  const session = feature(test, "features/filters/substructure-card.feature", import.meta.url);
  test("The substructure filter card and its search types", {tag: ["@journey", "@realizes:filters.cp.chem-and-bio-filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await run.scenario("The panel opens with substructure cards, and their close icons remove them", async () => {
      await session.step(20, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(21, "Then \"Structure\" filter card should be visible", () => shouldBe(page, el("\"Structure\" filter card"), "visible"));
      await session.step(22, "And the \"type of Structure\" reading of filter panel should be \"Chem:substructureFilter\"", () => readingReads(page, "type of Structure", el("filter panel"), "Chem:substructureFilter"));
      await session.step(23, "And the \"type of Core\" reading of filter panel should be \"Chem:substructureFilter\"", () => readingReads(page, "type of Core", el("filter panel"), "Chem:substructureFilter"));
      await session.step(24, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(25, "And user clicks on close of \"Structure\" filter card", () => clickOn(page, el("close of \"Structure\" filter card")));
      await session.step(26, "And user hovers over \"Core\" filter card", () => hoverOver(page, el("\"Core\" filter card")));
      await session.step(27, "And user clicks on close of \"Core\" filter card", () => clickOn(page, el("close of \"Core\" filter card")));
      await session.step(28, "And user hovers over \"R1\" filter card", () => hoverOver(page, el("\"R1\" filter card")));
      await session.step(29, "And user clicks on close of \"R1\" filter card", () => clickOn(page, el("close of \"R1\" filter card")));
      await session.step(30, "And user hovers over \"R2\" filter card", () => hoverOver(page, el("\"R2\" filter card")));
      await session.step(31, "And user clicks on close of \"R2\" filter card", () => clickOn(page, el("close of \"R2\" filter card")));
      await session.step(32, "And user hovers over \"R3\" filter card", () => hoverOver(page, el("\"R3\" filter card")));
      await session.step(33, "And user clicks on close of \"R3\" filter card", () => clickOn(page, el("close of \"R3\" filter card")));
      await session.step(34, "And user hovers over \"R100\" filter card", () => hoverOver(page, el("\"R100\" filter card")));
      await session.step(35, "And user clicks on close of \"R100\" filter card", () => clickOn(page, el("close of \"R100\" filter card")));
      await session.step(36, "And user hovers over \"R101\" filter card", () => hoverOver(page, el("\"R101\" filter card")));
      await session.step(37, "And user clicks on close of \"R101\" filter card", () => clickOn(page, el("close of \"R101\" filter card")));
      await session.step(38, "Then the filter panel should have no \"Chem:substructureFilter\" filter card", () => panelHasNoCardOfType(page, "Chem:substructureFilter"));
      await session.step(39, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Add Filter | Substructure Filter... puts a card on every checked molecule column", async () => {
      await session.step(43, "When user picks \"Add Filter | Substructure Filter...\" from the viewer menu of filter panel", () => pickFromViewerMenu(page, "Add Filter | Substructure Filter...", el("filter panel")));
      await session.step(44, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(45, "When user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(46, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(47, "Then \"Structure\" filter card should be visible", () => shouldBe(page, el("\"Structure\" filter card"), "visible"));
      await session.step(48, "And the \"type of Structure\" reading of filter panel should be \"Chem:substructureFilter\"", () => readingReads(page, "type of Structure", el("filter panel"), "Chem:substructureFilter"));
      await session.step(49, "And the \"type of R101\" reading of filter panel should be \"Chem:substructureFilter\"", () => readingReads(page, "type of R101", el("filter panel"), "Chem:substructureFilter"));
      await session.step(50, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(51, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pyridine on the Structure card keeps the rows that contain it", async () => {
      await session.step(55, "When user clicks on \"Sketch\" text in \"Structure\" filter card", () => clickOn(page, el("\"Sketch\" text in \"Structure\" filter card")));
      await session.step(56, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(57, "When user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(58, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(59, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(60, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(61, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(62, "And the filter should pass exactly the molecules of \"Structure\" column containing \"c1ccncc1\"", () => filterPassesMatching(page, "Structure", "c1ccncc1"));
      await session.step(63, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The search types split and nest the rows", async () => {
      await session.step(67, "When user opens the settings of the \"Structure\" filter card", () => openCardSettings(page, "Structure"));
      await session.step(68, "Then the \"Structure\" filter card should offer search types \"Contains, Included in, Exact, Stereo agnostic, Similar, Not contains, Not included in\"", () => offersSearchTypes(page, "Structure", "Contains, Included in, Exact, Stereo agnostic, Similar, Not contains, Not included in"));
      await session.step(69, "When user picks search type \"Not contains\" in the \"Structure\" filter card", () => pickSearchType(page, "Not contains", "Structure"));
      await session.step(70, "Then 83 rows should pass the filter", () => filterPasses(page, 83));
      await session.step(71, "And the \"search type of Structure\" reading of filter panel should be \"Not contains\"", () => readingReads(page, "search type of Structure", el("filter panel"), "Not contains"));
      await session.step(72, "When user picks search type \"Included in\" in the \"Structure\" filter card", () => pickSearchType(page, "Included in", "Structure"));
      await session.step(73, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(74, "When user picks search type \"Not included in\" in the \"Structure\" filter card", () => pickSearchType(page, "Not included in", "Structure"));
      await session.step(75, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(76, "When user picks search type \"Exact\" in the \"Structure\" filter card", () => pickSearchType(page, "Exact", "Structure"));
      await session.step(77, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(78, "When user picks search type \"Contains\" in the \"Structure\" filter card", () => pickSearchType(page, "Contains", "Structure"));
      await session.step(79, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(80, "When user picks search type \"Similar\" in the \"Structure\" filter card", () => pickSearchType(page, "Similar", "Structure"));
      await session.step(81, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(82, "When user picks search type \"Not contains\" in the \"Structure\" filter card", () => pickSearchType(page, "Not contains", "Structure"));
      await session.step(83, "Then 83 rows should pass the filter", () => filterPasses(page, 83));
      await session.step(84, "When user picks search type \"Stereo agnostic\" in the \"Structure\" filter card", () => pickSearchType(page, "Stereo agnostic", "Structure"));
      await session.step(85, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(86, "When user picks search type \"Contains\" in the \"Structure\" filter card", () => pickSearchType(page, "Contains", "Structure"));
      await session.step(87, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(88, "And the \"search type of Structure\" reading of filter panel should be \"Contains\"", () => readingReads(page, "search type of Structure", el("filter panel"), "Contains"));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
