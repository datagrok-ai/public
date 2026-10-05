/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/biostructure-viewer/grid-cell-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [biostructureviewer.cell.molecule3d, biostructureviewer.viewer.biostructure]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clipboardContains, doubleClickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {closeContextMenu, infoBalloonText, menuDoesNotList, menuLists, noBalloons, noErrors, pickFromAreaContextMenu, rightClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid context menu on Molecule3D cells (GROK-14552)", () => {
  const session = feature(test, "features/biostructure-viewer/grid-cell-menu.feature", import.meta.url);
  test("Grid context menu on Molecule3D cells (GROK-14552)", {tag: ["@journey", "@viewers", "@realizes:biostructureviewer.cell.molecule3d", "@realizes:biostructureviewer.viewer.biostructure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(24, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(25, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(26, "When user double-clicks on Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---pdb_data.csv tree node inside browse tree")));
    await session.step(27, "Then the \"pdb_data\" table view should open with 6 rows", () => tableViewOpened(page, "pdb_data", 6));
    await session.step(28, "And \"pdb\" column should have semantic type \"Molecule3D\"", () => columnSemType(page, "pdb", "Molecule3D"));
    await session.step(29, "And \"pdb_id\" column should have semantic type \"PDB_ID\"", () => columnSemType(page, "pdb_id", "PDB_ID"));
    await run.scenario("A Molecule3D cell's menu offers Copy, Download and Show > Biostructure / NGL", async () => {
      await session.step(32, "When user right-clicks on the \"cell 1 of pdb\" area of grid", () => rightClickArea(page, "cell 1 of pdb", el("grid")));
      await session.step(33, "Then the open menu should list \"Copy\"", () => menuLists(page, "Copy"));
      await session.step(34, "And the open menu should list \"Download\"", () => menuLists(page, "Download"));
      await session.step(35, "And the open menu should list \"Show > Biostructure\"", () => menuLists(page, "Show > Biostructure"));
      await session.step(36, "And the open menu should list \"Show > NGL\"", () => menuLists(page, "Show > NGL"));
      await session.step(37, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(38, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Copy puts the cell's PDB text on the clipboard", async () => {
      await session.step(41, "When user picks \"Copy\" from the context menu of the \"cell 1 of pdb\" area of grid", () => pickFromAreaContextMenu(page, "Copy", "cell 1 of pdb", el("grid")));
      await session.step(42, "Then an info balloon containing \"Value copied to clipboard\" should have been shown", () => infoBalloonText(page, "Value copied to clipboard"));
      await session.step(43, "And the clipboard should contain text \"HIV-1 PROTEASE INHIBITORS WIIH LOW NANOMOLAR POTENCY\"", () => clipboardContains(page, "HIV-1 PROTEASE INHIBITORS WIIH LOW NANOMOLAR POTENCY"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show > Biostructure docks a Biostructure viewer with the cell's structure", async () => {
      await session.step(47, "Then \"Reset Camera\" button inside open tableview should be absent", () => shouldBe(page, el("\"Reset Camera\" button inside open tableview"), "absent"));
      await session.step(48, "When user picks \"Show > Biostructure\" from the context menu of the \"cell 2 of pdb\" area of grid", () => pickFromAreaContextMenu(page, "Show > Biostructure", "cell 2 of pdb", el("grid")));
      await session.step(49, "Then \"Reset Camera\" button inside open tableview should be visible", () => shouldBe(page, el("\"Reset Camera\" button inside open tableview"), "visible"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
      await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Show > NGL on a Molecule3D cell raises no error", async () => {
      await session.step(54, "When user picks \"Show > NGL\" from the context menu of the \"cell 2 of pdb\" area of grid", () => pickFromAreaContextMenu(page, "Show > NGL", "cell 2 of pdb", el("grid")));
      await session.step(55, "Then no errors should have been logged", () => noErrors(page));
      await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A PDB_ID cell's menu has none of the Molecule3D items", async () => {
      await session.step(59, "When user right-clicks on the \"cell 1 of pdb_id\" area of grid", () => rightClickArea(page, "cell 1 of pdb_id", el("grid")));
      await session.step(60, "Then the open menu should list \"Properties...\"", () => menuLists(page, "Properties..."));
      await session.step(61, "And the open menu should not list \"Show\"", () => menuDoesNotList(page, "Show"));
      await session.step(62, "And the open menu should not list \"Copy\"", () => menuDoesNotList(page, "Copy"));
      await session.step(63, "And the open menu should not list \"Download\"", () => menuDoesNotList(page, "Download"));
      await session.step(64, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(65, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
