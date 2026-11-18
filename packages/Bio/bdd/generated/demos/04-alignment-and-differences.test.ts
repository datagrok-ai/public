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
import {clickOn, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag, everyValueSameLength} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridHeaderHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("4 · A peptide series · Alignment", () => {
  const session = feature(test, "features/demos/04-alignment-and-differences.feature", import.meta.url);
  test("Align a peptide family, then compare it with one peptide", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "And user opens FASTA_PT dataset", () => openDataset(page, ds("FASTA_PT")));
    await session.step(23, "Then \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
    await session.step(24, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(25, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(26, "And the grid shows only the columns \"id, sequence\"", () => gridShowsColumns(page, "id, sequence"));
    await session.step(27, "And the monomers of the \"sequence\" column are drawn at 16 px", () => sequenceFontSize(page, "sequence", 16));
    await session.step(28, "And the \"sequence\" column of the grid is 620 px wide", () => gridColumnWidth(page, "sequence", 620));
    await session.step(29, "When user picks \"Bio > Analyze > MSA...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > MSA..."));
    await session.step(30, "Then MSA dialog should be visible", () => shouldBe(page, el("MSA dialog"), "visible"));
    await session.step(31, "When user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")));
    await session.step(32, "Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
    await session.step(33, "And a new column \"msa(sequence)\" should have been added", () => newColumnNamed(page, "msa(sequence)"));
    await session.step(34, "And \"msa(sequence)\" column should have tag \"aligned\" equal to \"SEQ.MSA\"", () => columnTag(page, "msa(sequence)", "aligned", "SEQ.MSA"));
    await session.step(35, "And every value of \"msa(sequence)\" column should have the same length", () => everyValueSameLength(page, "msa(sequence)"));
    await session.step(36, "Given the grid shows only the columns \"id, msa(sequence)\"", () => gridShowsColumns(page, "id, msa(sequence)"));
    await session.step(37, "And the monomers of the \"msa(sequence)\" column are drawn at 16 px", () => sequenceFontSize(page, "msa(sequence)", 16));
    await session.step(38, "And the \"msa(sequence)\" column of the grid is 820 px wide", () => gridColumnWidth(page, "msa(sequence)", 820));
    await session.step(39, "And the column headers of the grid are 160 px tall", () => gridHeaderHeight(page, 160));
    await session.step(40, "When the video says \"Aligned: the shared motifs fall into columns\"", () => videoSays(page, "Aligned: the shared motifs fall into columns"));
    await session.step(41, "And user hovers over the \"position 8 of msa(sequence) header\" area of grid", () => hoverArea(page, "position 8 of msa(sequence) header", el("grid")));
    await session.step(42, "Then tooltip should contain text \"Position: 8\"", () => shouldContainText(page, el("tooltip"), "Position: 8"));
    await session.step(43, "When user hovers over the \"position 34 of msa(sequence) header\" area of grid", () => hoverArea(page, "position 34 of msa(sequence) header", el("grid")));
    await session.step(44, "Then tooltip should contain text \"Position: 34\"", () => shouldContainText(page, el("tooltip"), "Position: 34"));
    await session.step(45, "When the video says \"What differs from the clicked peptide stays in color\"", () => videoSays(page, "What differs from the clicked peptide stays in color"));
    await session.step(46, "And user clicks on the \"cell 7 of msa(sequence)\" area of grid", () => clickArea(page, "cell 7 of msa(sequence)", el("grid")));
    await session.step(47, "And user hovers over the \"cell 7 of id\" area of grid", () => hoverArea(page, "cell 7 of id", el("grid")));
    await session.step(48, "And user clicks on the \"cell 13 of msa(sequence)\" area of grid", () => clickArea(page, "cell 13 of msa(sequence)", el("grid")));
    await session.step(49, "And user hovers over the \"cell 13 of id\" area of grid", () => hoverArea(page, "cell 13 of id", el("grid")));
    await session.step(50, "Then no errors should have been logged", () => noErrors(page));
  });
});
