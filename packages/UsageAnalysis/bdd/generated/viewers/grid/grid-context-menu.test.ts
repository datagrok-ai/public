/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-context-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/nx.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, clipboardHas, pressKey, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {currentRowValue, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, closeContextMenu, hasArea, hasNoArea, menuLists, noErrors, pickFromAreaContextMenu, readingIs, rightClickArea, showsRows, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The context menu acts on the right-clicked cell", () => {
  const session = feature(test, "features/viewers/grid/grid-context-menu.feature", import.meta.url);
  test("A right click below the current row makes the clicked row current and keeps the scroll", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(17, "When user clicks on the \"cell 5 of AGE\" area of grid", () => clickArea(page, "cell 5 of AGE", el("grid")));
    await session.step(18, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(19, "When user scrolls the mouse wheel down 60 times over the \"cell 5 of AGE\" area of grid", () => wheelOverAreaTimes(page, "down", 60, "cell 5 of AGE", el("grid")));
    await session.step(20, "Then grid should have a \"cell 1000 of AGE\" area", () => hasArea(page, el("grid"), "cell 1000 of AGE"));
    await session.step(21, "And grid should not have a \"cell 5 of AGE\" area", () => hasNoArea(page, el("grid"), "cell 5 of AGE"));
    await session.step(22, "When user right-clicks on the \"cell 1000 of AGE\" area of grid", () => rightClickArea(page, "cell 1000 of AGE", el("grid")));
    await session.step(23, "Then the open menu should list \"Current Column\"", () => menuLists(page, "Current Column"));
    await session.step(24, "And row 1000 should be current", () => currentRowIs(page, 1000));
    await session.step(25, "And the \"current row\" reading of grid should be 1000", () => readingIs(page, "current row", el("grid"), 1000));
    await session.step(26, "And grid should have a \"cell 1000 of AGE\" area", () => hasArea(page, el("grid"), "cell 1000 of AGE"));
    await session.step(27, "And grid should not have a \"cell 5 of AGE\" area", () => hasNoArea(page, el("grid"), "cell 5 of AGE"));
    await session.step(28, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(29, "Then row 1000 should be current", () => currentRowIs(page, 1000));
    await session.step(30, "And grid should have a \"cell 1000 of AGE\" area", () => hasArea(page, el("grid"), "cell 1000 of AGE"));
    await session.step(31, "And grid should not have a \"cell 5 of AGE\" area", () => hasNoArea(page, el("grid"), "cell 5 of AGE"));
    await session.step(32, "And no errors should have been logged", () => noErrors(page));
  });
  test("Copy as SMILES copies the right-clicked molecule, not the one that was current", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "Given user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(36, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(37, "When user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")));
    await session.step(38, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(39, "And \"canonical_smiles\" of the current row should be \"FC(F)(F)c1ccc(OC2CCNCC2)cc1\"", () => currentRowValue(page, "canonical_smiles", "FC(F)(F)c1ccc(OC2CCNCC2)cc1"));
    await session.step(40, "When user scrolls the mouse wheel down 200 times over the \"cell 5 of canonical_smiles\" area of grid", () => wheelOverAreaTimes(page, "down", 200, "cell 5 of canonical_smiles", el("grid")));
    await session.step(41, "Then the \"row order\" reading of grid should contain \"1000\"", () => readingContains(page, "row order", el("grid"), "1000"));
    await session.step(42, "When user picks \"Current Value > Copy as SMILES\" from the context menu of the \"cell 1000 of canonical_smiles\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Copy as SMILES", "cell 1000 of canonical_smiles", el("grid")));
    await session.step(43, "Then the clipboard should have the text \"CC(C)CCOC(=O)c1ccccc1N\"", () => clipboardHas(page, "CC(C)CCOC(=O)c1ccccc1N"));
    await session.step(44, "And row 1000 should be current", () => currentRowIs(page, 1000));
    await session.step(45, "And \"canonical_smiles\" of the current row should be \"CC(C)CCOC(=O)c1ccccc1N\"", () => currentRowValue(page, "canonical_smiles", "CC(C)CCOC(=O)c1ccccc1N"));
    await session.step(46, "And the \"row order\" reading of grid should contain \"1000\"", () => readingContains(page, "row order", el("grid"), "1000"));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
  test("Edit Helm... opens the right-clicked peptide, not the one that was current", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(50, "Given user opens helm-peptides dataset", () => openDataset(page, ds("helm-peptides")));
    await session.step(51, "Then the table should have 540 rows", () => rowCount(page, 540));
    await session.step(52, "When user clicks on the \"cell 5 of HELM\" area of grid", () => clickArea(page, "cell 5 of HELM", el("grid")));
    await session.step(53, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(54, "When user scrolls the mouse wheel down 200 times over the \"cell 5 of HELM\" area of grid", () => wheelOverAreaTimes(page, "down", 200, "cell 5 of HELM", el("grid")));
    await session.step(55, "Then the \"row order\" reading of grid should contain \"540\"", () => readingContains(page, "row order", el("grid"), "540"));
    await session.step(56, "When user picks \"Current Value > Edit Helm...\" from the context menu of the \"cell 540 of HELM\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Edit Helm...", "cell 540 of HELM", el("grid")));
    await session.step(57, "Then HELM notation tab should be visible", () => shouldBe(page, el("HELM notation tab"), "visible"));
    await session.step(58, "When user clicks on HELM notation tab", () => clickOn(page, el("HELM notation tab")));
    await session.step(59, "Then HELM notation should contain the text \"PEPTIDE1{meI.hHis.Hcy.Q.T.W.Q.Phe_4NH2.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.N.meK}$$$$\"", () => shouldContainText(page, el("HELM notation"), "PEPTIDE1{meI.hHis.Hcy.Q.T.W.Q.Phe_4NH2.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.N.meK}$$$$"));
    await session.step(60, "And row 540 should be current", () => currentRowIs(page, 540));
    await session.step(61, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(62, "Then HELM notation tab should be hidden", () => shouldBe(page, el("HELM notation tab"), "hidden"));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
  });
});
