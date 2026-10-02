/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scaffold-tree/functions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.scaffold-tree-add-filter]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, noErrors, propertyShouldBe, readingIs, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu, readingIncludes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scaffold Tree — added from the menu, and its filter card", () => {
  const session = feature(test, "features/scaffold-tree/functions.feature", import.meta.url);
  test("Scaffold Tree — added from the menu, and its filter card", {tag: ["@journey", "@realizes:chem.cp.scaffold-tree-add-filter"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
    await run.scenario("The menu adds the viewer in its empty state", async () => {
      await session.step(17, "When user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
      await session.step(18, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(19, "And Scaffold Tree viewer should be visible", () => shouldBe(page, el("Scaffold Tree viewer"), "visible"));
      await session.step(20, "And the open tableview should have 1 Scaffold Tree viewer", () => viewerCount(page, 1, "Scaffold Tree"));
      await session.step(21, "And Scaffold Tree viewer should be bound to table \"smiles-50\"", () => boundTable(page, el("Scaffold Tree viewer"), "smiles-50"));
      await session.step(22, "And the \"molecule column\" reading of Scaffold Tree viewer should be \"canonical_smiles\"", () => readingReads(page, "molecule column", el("Scaffold Tree viewer"), "canonical_smiles"));
      await session.step(23, "And the \"nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 0));
      await session.step(24, "And the \"message\" reading of Scaffold Tree viewer should include the text \"Scaffold Tree is empty\"", () => readingIncludes(page, "message", el("Scaffold Tree viewer"), "Scaffold Tree is empty"));
      await session.step(25, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer reports the column and the size it holds", async () => {
      await session.step(28, "Then \"size\" property of Scaffold Tree viewer should be \"large\"", () => propertyShouldBe(page, "size", el("Scaffold Tree viewer"), "large"));
      await session.step(29, "And \"molecule\" property of Scaffold Tree viewer should be \"canonical_smiles\"", () => propertyShouldBe(page, "molecule", el("Scaffold Tree viewer"), "canonical_smiles"));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Add Filter | Scaffold Tree Filter... offers the molecule column alone", async () => {
      await session.step(33, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(34, "And user picks \"Add Filter | Scaffold Tree Filter...\" from the viewer menu of filter panel", () => pickFromViewerMenu(page, "Add Filter | Scaffold Tree Filter...", el("filter panel")));
      await session.step(35, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(36, "And the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"canonical_smiles\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "canonical_smiles"));
      await session.step(37, "When user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(38, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(39, "Then there should be 2 visible \"canonical_smiles\" filter card", () => visibleCount(page, 2, el("\"canonical_smiles\" filter card")));
      await session.step(40, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
