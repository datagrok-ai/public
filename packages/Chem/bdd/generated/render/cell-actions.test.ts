/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/render/cell-actions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.molecule-cell-actions]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clipboardMolecule, sortedBySimilarity} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clipboardContains, clipboardDiffers, clipboardImage, downloadContains, fileDownloaded, hoverOver, rememberClipboard, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, infoBalloonText, menuLists, noBalloons, noErrors, pickFromAreaContextMenu, rightClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Copy, export and sort from a molecule cell's context menu", () => {
  const session = feature(test, "features/render/cell-actions.feature", import.meta.url);
  test("Copy, export and sort from a molecule cell's context menu", {tag: ["@journey", "@realizes:chem.cp.molecule-cell-actions"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(11, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("The context menu offers the copy, export and sort actions", async () => {
      await session.step(14, "When user right-clicks on the \"cell 1 of canonical_smiles\" area of grid", () => rightClickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(15, "And user hovers over \"Current Value\" menu item", () => hoverOver(page, el("\"Current Value\" menu item")));
      await session.step(16, "Then the open menu should list \"Copy as SMILES\"", () => menuLists(page, "Copy as SMILES"));
      await session.step(17, "And the open menu should list \"Copy as MOLFILE V2000\"", () => menuLists(page, "Copy as MOLFILE V2000"));
      await session.step(18, "And the open menu should list \"Copy as MOLFILE V3000\"", () => menuLists(page, "Copy as MOLFILE V3000"));
      await session.step(19, "And the open menu should list \"Copy as SMARTS\"", () => menuLists(page, "Copy as SMARTS"));
      await session.step(20, "And the open menu should list \"Copy as Image\"", () => menuLists(page, "Copy as Image"));
      await session.step(21, "And the open menu should list \"Export as SVG\"", () => menuLists(page, "Export as SVG"));
      await session.step(22, "And the open menu should list \"Sort by similarity\"", () => menuLists(page, "Sort by similarity"));
      await session.step(23, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(24, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Each text copy holds the cell's molecule in its own notation", async () => {
      await session.step(27, "When user picks \"Current Value > Copy as SMILES\" from the context menu of the \"cell 1 of canonical_smiles\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Copy as SMILES", "cell 1 of canonical_smiles", el("grid")));
      await session.step(28, "Then an info balloon containing \"copied\" should have been shown", () => infoBalloonText(page, "copied"));
      await session.step(29, "And the clipboard should hold the molecule of row 1 of \"canonical_smiles\" column", () => clipboardMolecule(page, 1, "canonical_smiles"));
      await session.step(30, "When user remembers the clipboard text", () => rememberClipboard(page));
      await session.step(31, "And user picks \"Current Value > Copy as MOLFILE V2000\" from the context menu of the \"cell 1 of canonical_smiles\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Copy as MOLFILE V2000", "cell 1 of canonical_smiles", el("grid")));
      await session.step(32, "Then the clipboard should contain text \"V2000\"", () => clipboardContains(page, "V2000"));
      await session.step(33, "And the clipboard should contain text \"M  END\"", () => clipboardContains(page, "M  END"));
      await session.step(34, "And the clipboard should hold the molecule of row 1 of \"canonical_smiles\" column", () => clipboardMolecule(page, 1, "canonical_smiles"));
      await session.step(35, "And the clipboard text should differ from every remembered one", () => clipboardDiffers(page));
      await session.step(36, "When user remembers the clipboard text", () => rememberClipboard(page));
      await session.step(37, "And user picks \"Current Value > Copy as MOLFILE V3000\" from the context menu of the \"cell 1 of canonical_smiles\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Copy as MOLFILE V3000", "cell 1 of canonical_smiles", el("grid")));
      await session.step(38, "Then the clipboard should contain text \"V3000\"", () => clipboardContains(page, "V3000"));
      await session.step(39, "And the clipboard should contain text \"M  END\"", () => clipboardContains(page, "M  END"));
      await session.step(40, "And the clipboard should hold the molecule of row 1 of \"canonical_smiles\" column", () => clipboardMolecule(page, 1, "canonical_smiles"));
      await session.step(41, "And the clipboard text should differ from every remembered one", () => clipboardDiffers(page));
      await session.step(42, "When user remembers the clipboard text", () => rememberClipboard(page));
      await session.step(43, "And user picks \"Current Value > Copy as SMARTS\" from the context menu of the \"cell 1 of canonical_smiles\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Copy as SMARTS", "cell 1 of canonical_smiles", el("grid")));
      await session.step(44, "Then the clipboard should contain text \"#\"", () => clipboardContains(page, "#"));
      await session.step(45, "And the clipboard text should differ from every remembered one", () => clipboardDiffers(page));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Copy as Image puts a PNG on the clipboard", async () => {
      await session.step(49, "When user picks \"Current Value > Copy as Image\" from the context menu of the \"cell 1 of canonical_smiles\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Copy as Image", "cell 1 of canonical_smiles", el("grid")));
      await session.step(50, "Then an info balloon containing \"Image copied\" should have been shown", () => infoBalloonText(page, "Image copied"));
      await session.step(51, "And the clipboard should hold a PNG image of at least 500 bytes", () => clipboardImage(page, 500));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Export as SVG downloads a drawing of the molecule", async () => {
      await session.step(55, "Given user watches downloads", () => watchDownloads(page));
      await session.step(56, "When user picks \"Current Value > Export as SVG\" from the context menu of the \"cell 1 of canonical_smiles\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Export as SVG", "cell 1 of canonical_smiles", el("grid")));
      await session.step(57, "Then a file \"molecule.svg\" should have been downloaded", () => fileDownloaded(page, "molecule.svg"));
      await session.step(58, "And the downloaded file \"molecule.svg\" should contain text \"<svg\"", () => downloadContains(page, "molecule.svg", "<svg"));
      await session.step(59, "And the downloaded file \"molecule.svg\" should contain text \"<path\"", () => downloadContains(page, "molecule.svg", "<path"));
      await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sort by similarity puts the molecule first and the similar ones after it", async () => {
      await session.step(64, "When user picks \"Current Value > Sort by similarity\" from the context menu of the \"cell 3 of canonical_smiles\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Sort by similarity", "cell 3 of canonical_smiles", el("grid")));
      await session.step(65, "Then the first 5 rows of grid should be in falling similarity to row 3 of \"canonical_smiles\" column", () => sortedBySimilarity(page, 5, 3, "canonical_smiles"));
      await session.step(66, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
