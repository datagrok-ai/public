/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/chem-state.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.int.project-save-reopen-chem-state]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {filterPassesMatching} from '../../bindings/molecules.js';
import {treeBuilt} from '../../bindings/scaffold-tree.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKeyIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterPassesFewer, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, closeAllViews, openDataset, openProject, saveAsProject, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingAtLeast, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A saved project brings back the Chem state it was saved with", () => {
  const session = feature(test, "features/projects/chem-state.feature", import.meta.url);
  test("A saved project brings back the Chem state it was saved with", {tag: ["@journey", "@realizes:chem.int.project-save-reopen-chem-state"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await run.scenario("The search and the tree narrow the table before it is saved", async () => {
      await session.step(16, "When user picks \"Chem > Search > Substructure Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Substructure Search..."));
      await session.step(17, "Then \"Substructure search\" dialog should be visible", () => shouldBe(page, el("\"Substructure search\" dialog"), "visible"));
      await session.step(18, "When user clicks on OK button in \"Substructure search\" dialog", () => clickOn(page, el("OK button in \"Substructure search\" dialog")));
      await session.step(19, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(20, "When user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(21, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(22, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(23, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(24, "When user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
      await session.step(25, "And user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(26, "And user clicks on \"Generate\" icon inside Scaffold Tree viewer", () => clickOn(page, el("\"Generate\" icon inside Scaffold Tree viewer")));
      await session.step(27, "Then Scaffold Tree viewer should have finished building its tree", () => treeBuilt(page, el("Scaffold Tree viewer")));
      await session.step(28, "When user clicks on the \"checkbox of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 1", el("Scaffold Tree viewer")));
      await session.step(29, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 1));
      await session.step(30, "And fewer than 17 rows should pass the filter", () => filterPassesFewer(page, 17));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The project comes back with the tree, the card and the rows", async () => {
      await session.step(34, "When user saves the current view as project \"chem-state-roundtrip\"", () => saveAsProject(page, "chem-state-roundtrip"));
      await session.step(35, "And user closes all views", () => closeAllViews(page));
      await session.step(36, "And user opens the \"chem-state-roundtrip\" project", () => openProject(page, "chem-state-roundtrip"));
      await session.step(37, "Then Scaffold Tree viewer should be visible", () => shouldBe(page, el("Scaffold Tree viewer"), "visible"));
      await session.step(38, "And the \"nodes\" reading of Scaffold Tree viewer should be at least 1", () => readingAtLeast(page, "nodes", el("Scaffold Tree viewer"), 1));
      await session.step(39, "And the \"checked nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 1));
      await session.step(40, "And \"Structure\" filter card should be visible", () => shouldBe(page, el("\"Structure\" filter card"), "visible"));
      await session.step(41, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(42, "And the table should have 100 rows", () => rowCount(page, 100));
      await session.step(43, "And fewer than 17 rows should pass the filter", () => filterPassesFewer(page, 17));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Unchecking the restored node leaves the substructure search alone", async () => {
      await session.step(47, "When user clicks on the \"checkbox of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 1", el("Scaffold Tree viewer")));
      await session.step(48, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 0));
      await session.step(49, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(50, "And the filter should pass exactly the molecules of \"Structure\" column containing \"c1ccncc1\"", () => filterPassesMatching(page, "Structure", "c1ccncc1"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
