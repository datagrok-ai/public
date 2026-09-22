/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/editor/open.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [helm.cell-editor.molecule, helm.action.edit-helm]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {helmInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, followingShouldBe, shouldBe, shouldHaveText, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits, currentRowIs, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, noBalloons, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Opening the HELM editor from a cell", () => {
  const session = feature(test, "features/editor/open.feature", import.meta.url);
  test("Opening the HELM editor from a cell", {tag: ["@journey", "@realizes:helm.cell-editor.molecule", "@realizes:helm.action.edit-helm"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the Helm package is initialized", () => helmInitialized(page));
    await session.step(15, "And user opens helm-showcase dataset", () => openDataset(page, ds("helm-showcase")));
    await session.step(16, "Then \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
    await run.scenario("A double-click opens the editor on the cell's sequence", async () => {
      await session.step(19, "When user double-clicks on the \"cell 1 of HELM\" area of grid", () => doubleClickArea(page, "cell 1 of HELM", el("grid")));
      await session.step(20, "Then HELM editor should be visible", () => shouldBe(page, el("HELM editor"), "visible"));
      await session.step(21, "And there should be 2 visible drawn monomers", () => visibleCount(page, 2, el("drawn monomers")));
      await session.step(22, "And notation pane should have text \"PEPTIDE1{A.C}$$$$V2.0\"", () => shouldHaveText(page, el("notation pane"), "PEPTIDE1{A.C}$$$$V2.0"));
      await session.step(23, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["editor canvas"],["palette search"],["Favorites palette tab"],["Peptides palette tab"],["RNA palette tab"],["Sequence tab"],["HELM tab"],["Properties tab"],["undo button"],["redo button"],["clean layout button"],["OK button in HELM editor"],["CANCEL button in HELM editor"]]), [["editor canvas"],["palette search"],["Favorites palette tab"],["Peptides palette tab"],["RNA palette tab"],["Sequence tab"],["HELM tab"],["Properties tab"],["undo button"],["redo button"],["clean layout button"],["OK button in HELM editor"],["CANCEL button in HELM editor"]]);
      await session.step(37, "And \"Properties\" tab in HELM editor should be visible", () => shouldBe(page, el("\"Properties\" tab in HELM editor"), "visible"));
      await session.step(38, "And \"Structure View\" tab in HELM editor should be absent", () => shouldBe(page, el("\"Structure View\" tab in HELM editor"), "absent"));
      await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cancel closes the editor and leaves the cell unchanged", async () => {
      await session.step(43, "When user clicks on CANCEL button in HELM editor", () => clickOn(page, el("CANCEL button in HELM editor")));
      await session.step(44, "Then HELM editor should be absent", () => shouldBe(page, el("HELM editor"), "absent"));
      await session.step(45, "And the value of \"HELM\" column in row 1 should be \"PEPTIDE1{A.C}$$$$\"", () => valueInRow(page, "HELM", 1, "PEPTIDE1{A.C}$$$$"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit Helm... on the current cell opens the same editor on its sequence", async () => {
      await session.step(49, "When user clicks on the \"cell 2 of HELM\" area of grid", () => clickArea(page, "cell 2 of HELM", el("grid")));
      await session.step(50, "Then row 2 should be current", () => currentRowIs(page, 2));
      await session.step(51, "When user picks \"Current Value > Edit Helm...\" from the context menu of the \"cell 2 of HELM\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Edit Helm...", "cell 2 of HELM", el("grid")));
      await session.step(52, "Then HELM editor should be visible", () => shouldBe(page, el("HELM editor"), "visible"));
      await session.step(53, "And there should be 10 visible drawn monomers", () => visibleCount(page, 10, el("drawn monomers")));
      await session.step(54, "And notation pane should have text \"PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$V2.0\"", () => shouldHaveText(page, el("notation pane"), "PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$V2.0"));
      await session.step(55, "When user clicks on CANCEL button in HELM editor", () => clickOn(page, el("CANCEL button in HELM editor")));
      await session.step(56, "Then HELM editor should be absent", () => shouldBe(page, el("HELM editor"), "absent"));
      await session.step(57, "And the value of \"HELM\" column in row 2 should be \"PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$\"", () => valueInRow(page, "HELM", 2, "PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$"));
      await session.step(58, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit Helm... on another cell opens that cell, not the current one", async () => {
      await session.step(64, "When user clicks on the \"cell 1 of HELM\" area of grid", () => clickArea(page, "cell 1 of HELM", el("grid")));
      await session.step(65, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(66, "When user picks \"Current Value > Edit Helm...\" from the context menu of the \"cell 2 of HELM\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Edit Helm...", "cell 2 of HELM", el("grid")));
      await session.step(67, "Then HELM editor should be visible", () => shouldBe(page, el("HELM editor"), "visible"));
      await session.step(68, "And notation pane should have text \"PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$V2.0\"", () => shouldHaveText(page, el("notation pane"), "PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$V2.0"));
      await session.step(69, "When user clicks on OK button in HELM editor", () => clickOn(page, el("OK button in HELM editor")));
      await session.step(70, "Then HELM editor should be absent", () => shouldBe(page, el("HELM editor"), "absent"));
      await session.step(71, "And the value of \"HELM\" column in row 1 should be \"PEPTIDE1{A.C}$$$$\"", () => valueInRow(page, "HELM", 1, "PEPTIDE1{A.C}$$$$"));
      await session.step(72, "And the value of \"HELM\" column in row 2 should be \"PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$V2.0\"", () => valueInRow(page, "HELM", 2, "PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$V2.0"));
      await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
