/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scaffold-tree/scaffold-tree.feature
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
import {filterMatchesReading, filterPassesMatching, readingIsMolecule} from '../../bindings/molecules.js';
import {treeBuilt} from '../../bindings/scaffold-tree.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, followingShouldBe, hoverOver, pressKeyIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterPassesAll, filterPassesFewer} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, clickArea, hoverArea, noErrors, readingAtLeast, readingIs, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingIncludes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Scaffold Tree viewer — building, checking, editing and filtering", () => {
  const session = feature(test, "features/scaffold-tree/scaffold-tree.feature", import.meta.url);
  test("The Scaffold Tree viewer — building, checking, editing and filtering", {tag: ["@journey", "@realizes:chem.cp.scaffold-tree-add-filter"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await run.scenario("The menu adds the viewer in its empty state", async () => {
      await session.step(18, "When user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
      await session.step(19, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(20, "And Scaffold Tree viewer should be visible", () => shouldBe(page, el("Scaffold Tree viewer"), "visible"));
      await session.step(21, "And Scaffold Tree viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("Scaffold Tree viewer"), "spgi-100"));
      await session.step(22, "And the \"nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 0));
      await session.step(23, "And the \"message\" reading of Scaffold Tree viewer should include the text \"Scaffold Tree is empty\"", () => readingIncludes(page, "message", el("Scaffold Tree viewer"), "Scaffold Tree is empty"));
      await session.step(24, "And the \"generate blocked reason\" reading of Scaffold Tree viewer should be \"\"", () => readingReads(page, "generate blocked reason", el("Scaffold Tree viewer"), ""));
      await session.step(25, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The magic wand builds a tree from the molecule column", async () => {
      await session.step(29, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(30, "And user clicks on \"Generate\" icon inside Scaffold Tree viewer", () => clickOn(page, el("\"Generate\" icon inside Scaffold Tree viewer")));
      await session.step(31, "Then Scaffold Tree viewer should have finished building its tree", () => treeBuilt(page, el("Scaffold Tree viewer")));
      await session.step(32, "And the \"nodes\" reading of Scaffold Tree viewer should be at least 4", () => readingAtLeast(page, "nodes", el("Scaffold Tree viewer"), 4));
      await session.step(33, "And the \"root nodes\" reading of Scaffold Tree viewer should be at least 1", () => readingAtLeast(page, "root nodes", el("Scaffold Tree viewer"), 1));
      await session.step(34, "And the \"checked nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 0));
      await session.step(35, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Checking a node keeps the molecules that contain its scaffold", async () => {
      await session.step(39, "When user clicks on the \"checkbox of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 1", el("Scaffold Tree viewer")));
      await session.step(40, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 1));
      await session.step(41, "And the \"checked of node 1\" reading of Scaffold Tree viewer should be \"true\"", () => readingReads(page, "checked of node 1", el("Scaffold Tree viewer"), "true"));
      await session.step(42, "And fewer than 100 rows should pass the filter", () => filterPassesFewer(page, 100));
      await session.step(43, "And the \"hits of node 1\" reading of Scaffold Tree viewer should be at least 1", () => readingAtLeast(page, "hits of node 1", el("Scaffold Tree viewer"), 1));
      await session.step(44, "And the filter should pass exactly the molecules of \"Structure\" column containing the \"scaffold of node 1\" reading of Scaffold Tree viewer", () => filterMatchesReading(page, "Structure", "scaffold of node 1", el("Scaffold Tree viewer")));
      await session.step(45, "When user clicks on the \"checkbox of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 1", el("Scaffold Tree viewer")));
      await session.step(46, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 0));
      await session.step(47, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer's own toolbar offers its actions, and Clear filter unchecks the tree", async () => {
      await session.step(51, "When user clicks on the \"checkbox of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 1", el("Scaffold Tree viewer")));
      await session.step(52, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 1));
      await session.step(53, "And the filter should pass exactly the molecules of \"Structure\" column containing the \"scaffold of node 1\" reading of Scaffold Tree viewer", () => filterMatchesReading(page, "Structure", "scaffold of node 1", el("Scaffold Tree viewer")));
      await session.step(54, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(55, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Generate\" icon inside Scaffold Tree viewer"],["\"Sketch scaffolds manually\" icon inside Scaffold Tree viewer"],["\"Upload saved tree file\" icon inside Scaffold Tree viewer"],["\"Save this tree to disk\" icon inside Scaffold Tree viewer"],["\"Expand / collapse all\" icon inside Scaffold Tree viewer"],["\"Clear filter\" icon inside Scaffold Tree viewer"],["\"Drop all trees\" icon inside Scaffold Tree viewer"]]), [["\"Generate\" icon inside Scaffold Tree viewer"],["\"Sketch scaffolds manually\" icon inside Scaffold Tree viewer"],["\"Upload saved tree file\" icon inside Scaffold Tree viewer"],["\"Save this tree to disk\" icon inside Scaffold Tree viewer"],["\"Expand / collapse all\" icon inside Scaffold Tree viewer"],["\"Clear filter\" icon inside Scaffold Tree viewer"],["\"Drop all trees\" icon inside Scaffold Tree viewer"]]);
      await session.step(63, "When user clicks on \"Clear filter\" icon inside Scaffold Tree viewer", () => clickOn(page, el("\"Clear filter\" icon inside Scaffold Tree viewer")));
      await session.step(64, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 0));
      await session.step(65, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit scaffold replaces the structure and the filter follows", async () => {
      await session.step(69, "When user hovers over the \"node 1\" area of Scaffold Tree viewer", () => hoverArea(page, "node 1", el("Scaffold Tree viewer")));
      await session.step(70, "And user clicks on the \"edit icon of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "edit icon of node 1", el("Scaffold Tree viewer")));
      await session.step(71, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(72, "When user clears molecule input of sketcher dialog", () => clearField(page, el("molecule input of sketcher dialog")));
      await session.step(73, "And user types \"c1ccc2ncccc2c1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccc2ncccc2c1", el("molecule input of sketcher dialog")));
      await session.step(74, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(75, "And user clicks on \"Save\" button in sketcher dialog", () => clickOn(page, el("\"Save\" button in sketcher dialog")));
      await session.step(76, "Then the \"scaffold of node 1\" reading of Scaffold Tree viewer should be the molecule \"c1ccc2ncccc2c1\"", () => readingIsMolecule(page, "scaffold of node 1", el("Scaffold Tree viewer"), "c1ccc2ncccc2c1"));
      await session.step(77, "When user clicks on the \"checkbox of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 1", el("Scaffold Tree viewer")));
      await session.step(78, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 1));
      await session.step(79, "And the filter should pass exactly the molecules of \"Structure\" column containing \"c1ccc2ncccc2c1\"", () => filterPassesMatching(page, "Structure", "c1ccc2ncccc2c1"));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A scaffold sketched by hand keeps the molecules that contain it", async () => {
      await session.step(84, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(85, "And user clicks on \"Drop all trees\" icon inside Scaffold Tree viewer", () => clickOn(page, el("\"Drop all trees\" icon inside Scaffold Tree viewer")));
      await session.step(86, "And user clicks on \"Yes\" button in \"Delete Tree\" dialog", () => clickOn(page, el("\"Yes\" button in \"Delete Tree\" dialog")));
      await session.step(87, "Then the \"nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 0));
      await session.step(88, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(89, "And user clicks on \"Sketch scaffolds manually\" icon inside Scaffold Tree viewer", () => clickOn(page, el("\"Sketch scaffolds manually\" icon inside Scaffold Tree viewer")));
      await session.step(90, "And user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(91, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(92, "And user clicks on \"Add\" button in sketcher dialog", () => clickOn(page, el("\"Add\" button in sketcher dialog")));
      await session.step(93, "Then the \"nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 1));
      await session.step(94, "And the \"scaffold of node 1\" reading of Scaffold Tree viewer should be the molecule \"c1ccncc1\"", () => readingIsMolecule(page, "scaffold of node 1", el("Scaffold Tree viewer"), "c1ccncc1"));
      await session.step(95, "And the \"hits of node 1\" reading of Scaffold Tree viewer should be 17", () => readingIs(page, "hits of node 1", el("Scaffold Tree viewer"), 17));
      await session.step(96, "When user clicks on the \"checkbox of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 1", el("Scaffold Tree viewer")));
      await session.step(97, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(98, "And the filter should pass exactly the molecules of \"Structure\" column containing \"c1ccncc1\"", () => filterPassesMatching(page, "Structure", "c1ccncc1"));
    });
    await run.scenario("Cloning the view brings the same scaffold and the same filtered rows", async () => {
      await session.step(101, "When user picks \"View > Layout > Clone View\" from the top menu", () => pickFromTopMenu(page, "View > Layout > Clone View"));
      await session.step(102, "Then the open tableview should have 1 Scaffold Tree viewer", () => viewerCount(page, 1, "Scaffold Tree"));
      await session.step(103, "And the \"nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 1));
      await session.step(104, "And the \"scaffold of node 1\" reading of Scaffold Tree viewer should be the molecule \"c1ccncc1\"", () => readingIsMolecule(page, "scaffold of node 1", el("Scaffold Tree viewer"), "c1ccncc1"));
      await session.step(105, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
