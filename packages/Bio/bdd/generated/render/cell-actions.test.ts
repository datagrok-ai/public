/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/render/cell-actions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.actions.copy-as, bio.panel.composition-analysis, bio.panel.monomer]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, clipboardHas, expand, shouldBe, shouldContainText, shouldHaveRows} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset, openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, pickFromAreaContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sequence cell actions and panels", () => {
  const session = feature(test, "features/render/cell-actions.feature", import.meta.url);
  test("Sequence cell actions and panels", {tag: ["@journey", "@realizes:bio.actions.copy-as", "@realizes:bio.panel.composition-analysis", "@realizes:bio.panel.monomer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens filter_FASTA dataset keeping the first 9 rows", () => openDatasetRows(page, ds("filter_FASTA"), 9));
    await session.step(15, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(16, "Then \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
    await session.step(17, "And the \"cell type of fasta\" reading of grid should be \"sequence\"", () => readingReads(page, "cell type of fasta", el("grid"), "sequence"));
    await run.scenario("Copy as puts the cell on the clipboard in the chosen notation", async () => {
      await session.step(20, "When user picks \"Copy > helm\" from the context menu of the \"cell 1 of fasta\" area of grid", () => pickFromAreaContextMenu(page, "Copy > helm", "cell 1 of fasta", el("grid")));
      await session.step(21, "Then the clipboard should have the text \"PEPTIDE1{M.D.Y.K.E.T.L.L.M.P.K.T.D.F.P.M.R.G.G.L.P.N.K.E.P.Q.I.Q.E.K.W}$$$$\"", () => clipboardHas(page, "PEPTIDE1{M.D.Y.K.E.T.L.L.M.P.K.T.D.F.P.M.R.G.G.L.P.N.K.E.P.Q.I.Q.E.K.W}$$$$"));
      await session.step(22, "When user picks \"Copy > separator\" from the context menu of the \"cell 1 of fasta\" area of grid", () => pickFromAreaContextMenu(page, "Copy > separator", "cell 1 of fasta", el("grid")));
      await session.step(23, "Then the clipboard should have the text \"M.D.Y.K.E.T.L.L.M.P.K.T.D.F.P.M.R.G.G.L.P.N.K.E.P.Q.I.Q.E.K.W\"", () => clipboardHas(page, "M.D.Y.K.E.T.L.L.M.P.K.T.D.F.P.M.R.G.G.L.P.N.K.E.P.Q.I.Q.E.K.W"));
      await session.step(24, "When user picks \"Copy > fasta\" from the context menu of the \"cell 1 of fasta\" area of grid", () => pickFromAreaContextMenu(page, "Copy > fasta", "cell 1 of fasta", el("grid")));
      await session.step(25, "Then the clipboard should have the text \"MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW\"", () => clipboardHas(page, "MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW"));
      await session.step(26, "When user picks \"Copy > biln\" from the context menu of the \"cell 1 of fasta\" area of grid", () => pickFromAreaContextMenu(page, "Copy > biln", "cell 1 of fasta", el("grid")));
      await session.step(27, "Then the clipboard should have the text \"M-D-Y-K-E-T-L-L-M-P-K-T-D-F-P-M-R-G-G-L-P-N-K-E-P-Q-I-Q-E-K-W\"", () => clipboardHas(page, "M-D-Y-K-E-T-L-L-M-P-K-T-D-F-P-M-R-G-G-L-P-N-K-E-P-Q-I-Q-E-K-W"));
      await session.step(28, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The current cell shows its composition on the context panel", async () => {
      await session.step(31, "When user clicks on the \"cell 2 of fasta\" area of grid", () => clickArea(page, "cell 2 of fasta", el("grid")));
      await session.step(32, "Then the \"current row\" reading of grid should be 2", () => readingIs(page, "current row", el("grid"), 2));
      await session.step(33, "And the \"current column\" reading of grid should be \"fasta\"", () => readingReads(page, "current column", el("grid"), "fasta"));
      await session.step(34, "And \"Composition analysis\" section should be visible", () => shouldBe(page, el("\"Composition analysis\" section"), "visible"));
      await session.step(35, "When user expands \"Composition analysis\" section", () => expand(page, el("\"Composition analysis\" section")));
      await session.step(36, "Then \"Composition analysis\" section should have 14 rows", () => shouldHaveRows(page, el("\"Composition analysis\" section"), 14));
      await session.step(37, "And \"Composition analysis\" section should contain text \"%\"", () => shouldContainText(page, el("\"Composition analysis\" section"), "%"));
      await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A monomer cell shows the monomer on the context panel", async () => {
      await session.step(41, "Given user opens filter_MSA dataset", () => openDataset(page, ds("filter_MSA")));
      await session.step(42, "When user picks \"Bio > Transform > Split to Monomers...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Split to Monomers..."));
      await session.step(43, "And user clicks on OK button in \"Split to Monomers\" dialog", () => clickOn(page, el("OK button in \"Split to Monomers\" dialog")));
      await session.step(44, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(45, "And the value of \"1\" column in row 1 should be \"meI\"", () => valueInRow(page, "1", 1, "meI"));
      await session.step(46, "And the \"cell type of 1\" reading of grid should be \"Monomer\"", () => readingReads(page, "cell type of 1", el("grid"), "Monomer"));
      await session.step(47, "When user clicks on the \"cell 1 of 1\" area of grid", () => clickArea(page, "cell 1 of 1", el("grid")));
      await session.step(48, "Then the \"current column\" reading of grid should be \"1\"", () => readingReads(page, "current column", el("grid"), "1"));
      await session.step(49, "And \"Monomer\" section should be visible", () => shouldBe(page, el("\"Monomer\" section"), "visible"));
      await session.step(50, "When user expands \"Monomer\" section", () => expand(page, el("\"Monomer\" section")));
      await session.step(51, "Then \"Monomer\" section should contain text \"meI\"", () => shouldContainText(page, el("\"Monomer\" section"), "meI"));
      await session.step(52, "And \"Monomer\" section should contain text \"N-Methyl-Isoleucine\"", () => shouldContainText(page, el("\"Monomer\" section"), "N-Methyl-Isoleucine"));
      await session.step(53, "And \"Monomer\" section should contain text \"Backbone\"", () => shouldContainText(page, el("\"Monomer\" section"), "Backbone"));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
