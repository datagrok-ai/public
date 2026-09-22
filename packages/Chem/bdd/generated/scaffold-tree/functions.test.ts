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
import {treeBuilt} from '../../bindings/scaffold-tree.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, followingShouldBe, hoverOver, shouldBe, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPassesAll, filterPassesFewer, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, clickArea, noErrors, propertyShouldBe, readingAtLeast, readingIs, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu, readingIncludes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scaffold Tree — add, generate, filter, inspect", () => {
  const session = feature(test, "features/scaffold-tree/functions.feature", import.meta.url);
  test("Scaffold Tree — add, generate, filter, inspect", {tag: ["@journey", "@realizes:chem.cp.scaffold-tree-add-filter"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
    await run.scenario("The dataset opens with a molecule column", async () => {
      await session.step(18, "Then \"canonical_smiles\" column should have semantic type \"Molecule\"", () => columnSemType(page, "canonical_smiles", "Molecule"));
      await session.step(19, "And the table should have 50 rows", () => rowCount(page, 50));
      await session.step(20, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The menu adds the viewer in its empty state", async () => {
      await session.step(23, "When user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
      await session.step(24, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(25, "And Scaffold Tree viewer should be visible", () => shouldBe(page, el("Scaffold Tree viewer"), "visible"));
      await session.step(26, "And the open tableview should have 1 Scaffold Tree viewer", () => viewerCount(page, 1, "Scaffold Tree"));
      await session.step(27, "And Scaffold Tree viewer should be bound to table \"smiles-50\"", () => boundTable(page, el("Scaffold Tree viewer"), "smiles-50"));
      await session.step(28, "And the \"molecule column\" reading of Scaffold Tree viewer should be \"canonical_smiles\"", () => readingReads(page, "molecule column", el("Scaffold Tree viewer"), "canonical_smiles"));
      await session.step(29, "And the \"nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 0));
      await session.step(30, "And the \"message\" reading of Scaffold Tree viewer should include the text \"Scaffold Tree is empty\"", () => readingIncludes(page, "message", el("Scaffold Tree viewer"), "Scaffold Tree is empty"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The magic wand generates the scaffold hierarchy", async () => {
      await session.step(34, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(35, "Then \"Generate\" icon inside Scaffold Tree viewer should be enabled", () => shouldBe(page, el("\"Generate\" icon inside Scaffold Tree viewer"), "enabled"));
      await session.step(36, "When user clicks on \"Generate\" icon inside Scaffold Tree viewer", () => clickOn(page, el("\"Generate\" icon inside Scaffold Tree viewer")));
      await session.step(37, "Then Scaffold Tree viewer should have finished building its tree", () => treeBuilt(page, el("Scaffold Tree viewer")));
      await session.step(38, "Then the \"nodes\" reading of Scaffold Tree viewer should be at least 4", () => readingAtLeast(page, "nodes", el("Scaffold Tree viewer"), 4));
      await session.step(39, "And the \"root nodes\" reading of Scaffold Tree viewer should be at least 1", () => readingAtLeast(page, "root nodes", el("Scaffold Tree viewer"), 1));
      await session.step(40, "And the \"message\" reading of Scaffold Tree viewer should be \"\"", () => readingReads(page, "message", el("Scaffold Tree viewer"), ""));
      await session.step(41, "And the \"checked nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 0));
      await session.step(42, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on the first scaffold node filters the table", async () => {
      await session.step(46, "When user clicks on the \"checkbox of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 1", el("Scaffold Tree viewer")));
      await session.step(47, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 1));
      await session.step(48, "And the \"checked of node 1\" reading of Scaffold Tree viewer should be \"true\"", () => readingReads(page, "checked of node 1", el("Scaffold Tree viewer"), "true"));
      await session.step(49, "And fewer than 50 rows should pass the filter", () => filterPassesFewer(page, 50));
      await session.step(50, "And the \"hits of node 1\" reading of Scaffold Tree viewer should be at least 1", () => readingAtLeast(page, "hits of node 1", el("Scaffold Tree viewer"), 1));
      await session.step(51, "And the \"rows kept\" reading of Scaffold Tree viewer should be at least 1", () => readingAtLeast(page, "rows kept", el("Scaffold Tree viewer"), 1));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer's own toolbar offers its actions", async () => {
      await session.step(55, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(56, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Generate\" icon inside Scaffold Tree viewer"],["\"Sketch scaffolds manually\" icon inside Scaffold Tree viewer"],["\"Upload saved tree file\" icon inside Scaffold Tree viewer"],["\"Save this tree to disk\" icon inside Scaffold Tree viewer"],["\"Expand / collapse all\" icon inside Scaffold Tree viewer"],["\"Clear filter\" icon inside Scaffold Tree viewer"],["\"Drop all trees\" icon inside Scaffold Tree viewer"]]), [["\"Generate\" icon inside Scaffold Tree viewer"],["\"Sketch scaffolds manually\" icon inside Scaffold Tree viewer"],["\"Upload saved tree file\" icon inside Scaffold Tree viewer"],["\"Save this tree to disk\" icon inside Scaffold Tree viewer"],["\"Expand / collapse all\" icon inside Scaffold Tree viewer"],["\"Clear filter\" icon inside Scaffold Tree viewer"],["\"Drop all trees\" icon inside Scaffold Tree viewer"]]);
      await session.step(64, "When user clicks on \"Clear filter\" icon inside Scaffold Tree viewer", () => clickOn(page, el("\"Clear filter\" icon inside Scaffold Tree viewer")));
      await session.step(65, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 0));
      await session.step(66, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer reports the column and the size it was built with", async () => {
      await session.step(70, "Then \"size\" property of Scaffold Tree viewer should be \"large\"", () => propertyShouldBe(page, "size", el("Scaffold Tree viewer"), "large"));
      await session.step(71, "And \"molecule\" property of Scaffold Tree viewer should be \"canonical_smiles\"", () => propertyShouldBe(page, "molecule", el("Scaffold Tree viewer"), "canonical_smiles"));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Add Filter | Scaffold Tree Filter... offers the molecule column alone", async () => {
      await session.step(75, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(76, "And user picks \"Add Filter | Scaffold Tree Filter...\" from the viewer menu of filter panel", () => pickFromViewerMenu(page, "Add Filter | Scaffold Tree Filter...", el("filter panel")));
      await session.step(77, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(78, "And the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"canonical_smiles\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "canonical_smiles"));
      await session.step(79, "When user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(80, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(81, "Then there should be 2 visible \"canonical_smiles\" filter card", () => visibleCount(page, 2, el("\"canonical_smiles\" filter card")));
      await session.step(82, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
