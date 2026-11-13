/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scaffold-tree/scaffold-tree-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.scaffold-tree-add-filter, public#2315, public#2459, public#2154, public#2135, public#2964, public#3408, public#2322]
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
import {check, clearField, clickOn, hoverOver, pressKeyIn, shouldBe, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, hasColumn, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, filterPasses, filterPassesAll, resetFilter, selectFirstRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {addViewer, clickArea, hoverArea, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingAsRemembered, readingIs, readingReads, rememberReading, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scaffold Tree regressions from closed issues", () => {
  const session = feature(test, "features/scaffold-tree/scaffold-tree-regressions.feature", import.meta.url);
  test("Scaffold Tree regressions from closed issues", {tag: ["@journey", "@realizes:chem.cp.scaffold-tree-add-filter", "@realizes:public#2315", "@realizes:public#2459", "@realizes:public#2154", "@realizes:public#2135", "@realizes:public#2964", "@realizes:public#3408", "@realizes:public#2322"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await run.scenario("Order or Hide Columns lists the same columns after a tree with a checked scaffold is added (public#2315)", async () => {
      await session.step(23, "When user picks \"Order or Hide Columns...\" from the context menu of the \"header Id\" area of grid", () => pickFromAreaContextMenu(page, "Order or Hide Columns...", "header Id", el("grid")));
      await session.step(24, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
      await session.step(25, "And the \"rows shown\" reading of Grid viewer in Order or Hide Columns dialog should be 88", () => readingIs(page, "rows shown", el("Grid viewer in Order or Hide Columns dialog"), 88));
      await session.step(26, "When user remembers the \"rows shown\" reading of Grid viewer in Order or Hide Columns dialog", () => rememberReading(page, "rows shown", el("Grid viewer in Order or Hide Columns dialog")));
      await session.step(27, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
      await session.step(28, "And user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
      await session.step(29, "And user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(30, "And user clicks on \"Sketch scaffolds manually\" icon inside Scaffold Tree viewer", () => clickOn(page, el("\"Sketch scaffolds manually\" icon inside Scaffold Tree viewer")));
      await session.step(31, "And user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(32, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(33, "And user clicks on \"Add\" button in sketcher dialog", () => clickOn(page, el("\"Add\" button in sketcher dialog")));
      await session.step(34, "Then the \"hits of node 1\" reading of Scaffold Tree viewer should be 17", () => readingIs(page, "hits of node 1", el("Scaffold Tree viewer"), 17));
      await session.step(35, "When user clicks on the \"checkbox of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 1", el("Scaffold Tree viewer")));
      await session.step(36, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(37, "When user picks \"Order or Hide Columns...\" from the context menu of the \"header Id\" area of grid", () => pickFromAreaContextMenu(page, "Order or Hide Columns...", "header Id", el("grid")));
      await session.step(38, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
      await session.step(39, "And the \"rows shown\" reading of Grid viewer in Order or Hide Columns dialog should be as remembered", () => readingAsRemembered(page, "rows shown", el("Grid viewer in Order or Hide Columns dialog")));
      await session.step(40, "When user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
      await session.step(41, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The filter panel's reset unchecks the checked scaffold (public#2459, public#2154)", async () => {
      await session.step(45, "When user clicks on the \"checkbox of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 1", el("Scaffold Tree viewer")));
      await session.step(46, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 0));
      await session.step(47, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(48, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(49, "And user clicks on the \"checkbox of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 1", el("Scaffold Tree viewer")));
      await session.step(50, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 1));
      await session.step(51, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(52, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(53, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(54, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(55, "And the \"checked nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 0));
      await session.step(56, "And the \"checked of node 1\" reading of Scaffold Tree viewer should be \"false\"", () => readingReads(page, "checked of node 1", el("Scaffold Tree viewer"), "false"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Editing one scaffold keeps another one checked (public#2135)", async () => {
      await session.step(61, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(62, "And user clicks on \"Sketch scaffolds manually\" icon inside Scaffold Tree viewer", () => clickOn(page, el("\"Sketch scaffolds manually\" icon inside Scaffold Tree viewer")));
      await session.step(63, "And user types \"c1ccccc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccccc1", el("molecule input of sketcher dialog")));
      await session.step(64, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(65, "And user clicks on \"Add\" button in sketcher dialog", () => clickOn(page, el("\"Add\" button in sketcher dialog")));
      await session.step(66, "Then the \"nodes\" reading of Scaffold Tree viewer should be 2", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 2));
      await session.step(67, "And the \"scaffold of node 2\" reading of Scaffold Tree viewer should be the molecule \"c1ccccc1\"", () => readingIsMolecule(page, "scaffold of node 2", el("Scaffold Tree viewer"), "c1ccccc1"));
      await session.step(68, "And the \"hits of node 2\" reading of Scaffold Tree viewer should be 32", () => readingIs(page, "hits of node 2", el("Scaffold Tree viewer"), 32));
      await session.step(69, "When user clicks on the \"checkbox of node 2\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 2", el("Scaffold Tree viewer")));
      await session.step(70, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 1));
      await session.step(71, "And 32 rows should pass the filter", () => filterPasses(page, 32));
      await session.step(72, "When user hovers over the \"node 1\" area of Scaffold Tree viewer", () => hoverArea(page, "node 1", el("Scaffold Tree viewer")));
      await session.step(73, "And user clicks on the \"edit icon of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "edit icon of node 1", el("Scaffold Tree viewer")));
      await session.step(74, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(75, "When user clears molecule input of sketcher dialog", () => clearField(page, el("molecule input of sketcher dialog")));
      await session.step(76, "And user types \"c1ccc2ncccc2c1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccc2ncccc2c1", el("molecule input of sketcher dialog")));
      await session.step(77, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(78, "And user clicks on \"Save\" button in sketcher dialog", () => clickOn(page, el("\"Save\" button in sketcher dialog")));
      await session.step(79, "Then the \"scaffold of node 1\" reading of Scaffold Tree viewer should be the molecule \"c1ccc2ncccc2c1\"", () => readingIsMolecule(page, "scaffold of node 1", el("Scaffold Tree viewer"), "c1ccc2ncccc2c1"));
      await session.step(80, "And the \"checked of node 2\" reading of Scaffold Tree viewer should be \"true\"", () => readingReads(page, "checked of node 2", el("Scaffold Tree viewer"), "true"));
      await session.step(81, "And the \"checked nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 1));
      await session.step(82, "And 32 rows should pass the filter", () => filterPasses(page, 32));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Selection to Filter holds its filter while a scaffold tree is in the view (public#2964)", async () => {
      await session.step(87, "When user clicks on the \"checkbox of node 2\" area of Scaffold Tree viewer", () => clickArea(page, "checkbox of node 2", el("Scaffold Tree viewer")));
      await session.step(88, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 0));
      await session.step(89, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(90, "When user selects the first 1 rows", () => selectFirstRows(page, 1));
      await session.step(91, "And user picks \"Select > Selection to Filter\" from the top menu", () => pickFromTopMenu(page, "Select > Selection to Filter"));
      await session.step(92, "Then the \"checked nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "checked nodes", el("Scaffold Tree viewer"), 0));
      await session.step(93, "And the \"nodes\" reading of Scaffold Tree viewer should be 2", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 2));
      await session.step(94, "And 1 row should pass the filter", () => filterPasses(page, 1));
      await session.step(95, "When user clears the row selection", () => clearSelection(page));
      await session.step(96, "And user resets the filter", () => resetFilter(page));
      await session.step(97, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(98, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing a viewer with a colored scaffold drops the colors column and the highlight (public#3408)", async () => {
      await session.step(102, "When user hovers over the \"node 1\" area of Scaffold Tree viewer", () => hoverArea(page, "node 1", el("Scaffold Tree viewer")));
      await session.step(103, "And user clicks on the \"color icon of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "color icon of node 1", el("Scaffold Tree viewer")));
      await session.step(104, "Then the \"colored nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "colored nodes", el("Scaffold Tree viewer"), 1));
      await session.step(105, "And the table should have a column \"Structure colors\"", () => hasColumn(page, "Structure colors"));
      await session.step(106, "When user clicks on close icon of Scaffold Tree viewer", () => clickOn(page, el("close icon of Scaffold Tree viewer")));
      await session.step(107, "Then the open tableview should have 0 Scaffold Tree viewers", () => viewerCount(page, 0, "Scaffold Tree"));
      await session.step(108, "And the table should not have a column \"Structure colors\"", () => hasNoColumn(page, "Structure colors"));
      await session.step(109, "And \"Structure\" column should have tag \".chem-scaffold-tree-highlight\" equal to \"[]\"", () => columnTag(page, "Structure", ".chem-scaffold-tree-highlight", "[]"));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A property changed with a scaffold tree filter on the panel keeps the viewer's properties in the context panel (public#2322)", async () => {
      await session.step(114, "When user picks \"Add Filter | Scaffold Tree Filter...\" from the viewer menu of filter panel", () => pickFromViewerMenu(page, "Add Filter | Scaffold Tree Filter...", el("filter panel")));
      await session.step(115, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(116, "And the \"text of cell 1 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"Structure\"", () => readingReads(page, "text of cell 1 of __name", el("grid viewer in \"Select columns...\" dialog"), "Structure"));
      await session.step(117, "When user clicks on the \"cell 1 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 1 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(118, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(119, "And user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(120, "And user clicks on \"Sketch scaffolds manually\" icon inside filter panel", () => clickOn(page, el("\"Sketch scaffolds manually\" icon inside filter panel")));
      await session.step(121, "And user types \"c1ccccc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccccc1", el("molecule input of sketcher dialog")));
      await session.step(122, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(123, "And user clicks on \"Add\" button in sketcher dialog", () => clickOn(page, el("\"Add\" button in sketcher dialog")));
      await session.step(124, "Then sketcher dialog should be hidden", () => shouldBe(page, el("sketcher dialog"), "hidden"));
      await session.step(125, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(126, "When user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(127, "And user checks \"Show Filtered Out Points\" property in context panel", () => check(page, el("\"Show Filtered Out Points\" property in context panel")));
      await session.step(128, "Then \"showFilteredOutPoints\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "showFilteredOutPoints", el("scatter plot viewer"), "true"));
      await session.step(129, "And \"Show Filtered Out Points\" property in context panel should be visible", () => shouldBe(page, el("\"Show Filtered Out Points\" property in context panel"), "visible"));
      await session.step(130, "When user unchecks \"Show Filtered Out Points\" property in context panel", () => uncheck(page, el("\"Show Filtered Out Points\" property in context panel")));
      await session.step(131, "Then \"showFilteredOutPoints\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "showFilteredOutPoints", el("scatter plot viewer"), "false"));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
