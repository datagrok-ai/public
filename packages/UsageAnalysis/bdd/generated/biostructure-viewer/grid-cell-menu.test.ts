/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/biostructure-viewer/grid-cell-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [biostructureviewer.cell.molecule3d, biostructureviewer.viewer.biostructure]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, clipboardContains, doubleClickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {closeContextMenu, infoBalloonText, menuDoesNotList, menuLists, noBalloons, noErrors, pickFromAreaContextMenu, readingReads, rightClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {rightClickRowSpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid context menu on Molecule3D cells (GROK-14552)", () => {
  const session = feature(test, "features/biostructure-viewer/grid-cell-menu.feature", import.meta.url);
  test("Grid context menu on Molecule3D cells (GROK-14552)", {tag: ["@journey", "@viewers", "@realizes:biostructureviewer.cell.molecule3d", "@realizes:biostructureviewer.viewer.biostructure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(18, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(19, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(20, "When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree")));
    await session.step(21, "Then the \"pdb_data\" table view should open with 6 rows", () => tableViewOpened(page, "pdb_data", 6));
    await session.step(22, "And \"pdb\" column should have semantic type \"Molecule3D\"", () => columnSemType(page, "pdb", "Molecule3D"));
    await session.step(23, "And \"pdb_id\" column should have semantic type \"PDB_ID\"", () => columnSemType(page, "pdb_id", "PDB_ID"));
    await run.scenario("A Molecule3D cell's menu offers Copy, Download and Show > Biostructure / NGL", async () => {
      await session.step(26, "When user right-clicks on the \"cell 1 of pdb\" area of grid", () => rightClickArea(page, "cell 1 of pdb", el("grid")));
      await session.step(27, "Then the open menu should list \"Copy\"", () => menuLists(page, "Copy"));
      await session.step(28, "And the open menu should list \"Download\"", () => menuLists(page, "Download"));
      await session.step(29, "And the open menu should list \"Show > Biostructure\"", () => menuLists(page, "Show > Biostructure"));
      await session.step(30, "And the open menu should list \"Show > NGL\"", () => menuLists(page, "Show > NGL"));
      await session.step(31, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(32, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Copy puts the cell's PDB text on the clipboard", async () => {
      await session.step(35, "When user picks \"Copy\" from the context menu of the \"cell 1 of pdb\" area of grid", () => pickFromAreaContextMenu(page, "Copy", "cell 1 of pdb", el("grid")));
      await session.step(36, "Then an info balloon containing \"Value copied to clipboard\" should have been shown", () => infoBalloonText(page, "Value copied to clipboard"));
      await session.step(37, "And the clipboard should contain text \"HIV-1 PROTEASE INHIBITORS WIIH LOW NANOMOLAR POTENCY\"", () => clipboardContains(page, "HIV-1 PROTEASE INHIBITORS WIIH LOW NANOMOLAR POTENCY"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show > Biostructure docks a Biostructure viewer with the cell's structure", async () => {
      await session.step(41, "Then Biostructure viewer should be absent", () => shouldBe(page, el("Biostructure viewer"), "absent"));
      await session.step(42, "When user picks \"Show > Biostructure\" from the context menu of the \"cell 2 of pdb\" area of grid", () => pickFromAreaContextMenu(page, "Show > Biostructure", "cell 2 of pdb", el("grid")));
      await session.step(43, "Then Biostructure viewer should be visible", () => shouldBe(page, el("Biostructure viewer"), "visible"));
      await session.step(44, "And the \"structure loaded\" reading of Biostructure viewer should be \"true\"", () => readingReads(page, "structure loaded", el("Biostructure viewer"), "true"));
      await session.step(45, "And \"Reset Camera\" button in Biostructure viewer should be visible", () => shouldBe(page, el("\"Reset Camera\" button in Biostructure viewer"), "visible"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(48, "When user clicks on close icon of Biostructure viewer", () => clickOn(page, el("close icon of Biostructure viewer")));
      await session.step(49, "Then Biostructure viewer should be absent", () => shouldBe(page, el("Biostructure viewer"), "absent"));
    });
    await run.scenario("Show > NGL docks an NGL viewer with the cell's structure", async () => {
      await session.step(52, "Then NGL viewer should be absent", () => shouldBe(page, el("NGL viewer"), "absent"));
      await session.step(53, "When user picks \"Show > NGL\" from the context menu of the \"cell 2 of pdb\" area of grid", () => pickFromAreaContextMenu(page, "Show > NGL", "cell 2 of pdb", el("grid")));
      await session.step(54, "Then NGL viewer should be visible", () => shouldBe(page, el("NGL viewer"), "visible"));
      await session.step(55, "And the \"structure loaded\" reading of NGL viewer should be \"true\"", () => readingReads(page, "structure loaded", el("NGL viewer"), "true"));
      await session.step(56, "And \"Open...\" link in NGL viewer should be absent", () => shouldBe(page, el("\"Open...\" link in NGL viewer"), "absent"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
      await session.step(58, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(59, "When user clicks on close icon of NGL viewer", () => clickOn(page, el("close icon of NGL viewer")));
      await session.step(60, "Then NGL viewer should be absent", () => shouldBe(page, el("NGL viewer"), "absent"));
    });
    await run.scenario("A PDB_ID cell's menu has none of the Molecule3D items", async () => {
      await session.step(63, "When user right-clicks on the \"cell 1 of pdb_id\" area of grid", () => rightClickArea(page, "cell 1 of pdb_id", el("grid")));
      await session.step(64, "Then the open menu should list \"Properties...\"", () => menuLists(page, "Properties..."));
      await session.step(65, "And the open menu should not list \"Show\"", () => menuDoesNotList(page, "Show"));
      await session.step(66, "And the open menu should not list \"Copy\"", () => menuDoesNotList(page, "Copy"));
      await session.step(67, "And the open menu should not list \"Download\"", () => menuDoesNotList(page, "Download"));
      await session.step(68, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(69, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Right-clicking the row's empty space past the last column raises nothing (GROK-14552)", async () => {
      await session.step(72, "When user right-clicks on empty space of row 2 of grid", () => rightClickRowSpace(page, 2, el("grid")));
      await session.step(73, "Then the open menu should list \"Properties...\"", () => menuLists(page, "Properties..."));
      await session.step(74, "And the open menu should not list \"Show\"", () => menuDoesNotList(page, "Show"));
      await session.step(75, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(76, "Then no errors should have been logged", () => noErrors(page));
      await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
