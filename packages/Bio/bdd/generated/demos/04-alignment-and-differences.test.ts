/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/04-alignment-and-differences.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridHeaderHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Alignment and differences", () => {
  const session = feature(test, "features/demos/04-alignment-and-differences.feature", import.meta.url);
  test("Align, then compare every sequence with the one you click", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(15, "And user opens FASTA_sample dataset", () => openDataset(page, ds("FASTA_sample")));
    await session.step(16, "Then \"Sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "Sequence", "Macromolecule"));
    await session.step(17, "Given the grid shows only the columns \"Entry, Sequence, Activity\"", () => gridShowsColumns(page, "Entry, Sequence, Activity"));
    await session.step(18, "And the \"Sequence\" column of the grid is 330 px wide", () => gridColumnWidth(page, "Sequence", 330));
    await session.step(19, "When user picks \"Bio > Analyze > MSA...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > MSA..."));
    await session.step(20, "Then MSA dialog should be visible", () => shouldBe(page, el("MSA dialog"), "visible"));
    await session.step(21, "When user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
    await session.step(22, "Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
    await session.step(23, "And a new column \"msa(Sequence)\" should have been added", () => newColumnNamed(page, "msa(Sequence)"));
    await session.step(24, "And \"msa(Sequence)\" column should have tag \"aligned\" equal to \"SEQ.MSA\"", () => columnTag(page, "msa(Sequence)", "aligned", "SEQ.MSA"));
    await session.step(25, "Given the grid shows only the columns \"Entry, msa(Sequence), Activity\"", () => gridShowsColumns(page, "Entry, msa(Sequence), Activity"));
    await session.step(26, "And the \"msa(Sequence)\" column of the grid is 820 px wide", () => gridColumnWidth(page, "msa(Sequence)", 820));
    await session.step(27, "And the monomers of the \"msa(Sequence)\" column are drawn at 16 px", () => sequenceFontSize(page, "msa(Sequence)", 16));
    await session.step(28, "And the column headers of the grid are 160 px tall", () => gridHeaderHeight(page, 160));
    await session.step(29, "When the video says \"Every position: its conservation and its composition\"", () => videoSays(page, "Every position: its conservation and its composition"));
    await session.step(30, "And user hovers over the \"position 20 of msa(Sequence) header\" area of grid", () => hoverArea(page, "position 20 of msa(Sequence) header", el("grid")));
    await session.step(31, "And user clicks on the \"cell 7 of msa(Sequence)\" area of grid", () => clickArea(page, "cell 7 of msa(Sequence)", el("grid")));
    await session.step(32, "And user hovers over the \"cell 7 of Entry\" area of grid", () => hoverArea(page, "cell 7 of Entry", el("grid")));
    await session.step(33, "And the video says \"Shared monomers fade: what differs from the clicked row stands out\"", () => videoSays(page, "Shared monomers fade: what differs from the clicked row stands out"));
    await session.step(34, "And user clicks on the \"cell 13 of msa(Sequence)\" area of grid", () => clickArea(page, "cell 13 of msa(Sequence)", el("grid")));
    await session.step(35, "And user hovers over the \"cell 13 of Entry\" area of grid", () => hoverArea(page, "cell 13 of Entry", el("grid")));
    await session.step(36, "Then no errors should have been logged", () => noErrors(page));
  });
});
