/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/09-sequences/02-align-sequences.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bioReady, columnStaysWide} from '../../../bindings/howto-sequences.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag, everyValueMatches, everyValueSameLength} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, hoverArea, noBalloons, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridHeaderHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Sequences", () => {
  const session = feature(test, "features/howto/09-sequences/02-align-sequences.feature", import.meta.url);
  test("Align sequences and split them into positions", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the Bio package is initialized", () => bioReady(page));
    await session.step(21, "And user opens FASTA_PT dataset", () => openDataset(page, ds("FASTA_PT")));
    await session.step(22, "Then \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
    await session.step(23, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(24, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(25, "And the grid shows only the columns \"id, sequence\"", () => gridShowsColumns(page, "id, sequence"));
    await session.step(26, "And the \"sequence\" column of the grid stays 520 px wide", () => columnStaysWide(page, "sequence", 520));
    await session.step(27, "And the column headers of the grid are 160 px tall", () => gridHeaderHeight(page, 160));
    await session.step(28, "When the video says \"Related peptides, not lined up: a shared motif like EAM.K sits at different places\"", () => videoSays(page, "Related peptides, not lined up: a shared motif like EAM.K sits at different places"));
    await session.step(30, "And user picks \"Bio > Analyze > MSA...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > MSA..."), undefined, "Bio > Analyze > MSA aligns them");
    await session.step(31, "Then MSA dialog should be visible", () => shouldBe(page, el("MSA dialog"), "visible"));
    await session.step(33, "When user clicks on OK button in MSA dialog", () => clickOn(page, el("OK button in MSA dialog")), undefined, "It runs in your browser");
    await session.step(34, "Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
    await session.step(35, "And a new column \"msa(sequence)\" should have been added", () => newColumnNamed(page, "msa(sequence)"));
    await session.step(36, "And \"msa(sequence)\" column should have tag \"aligned\" equal to \"SEQ.MSA\"", () => columnTag(page, "msa(sequence)", "aligned", "SEQ.MSA"));
    await session.step(37, "And every value of \"msa(sequence)\" column should have the same length", () => everyValueSameLength(page, "msa(sequence)"));
    await session.step(38, "When the video says \"Aligned: the shared motifs fall into columns\" at the bottom", () => videoSaysAt(page, "Aligned: the shared motifs fall into columns", "bottom"));
    await session.step(39, "And the video says \"The header: how conserved each position is, and its letters\" at the bottom", () => videoSaysAt(page, "The header: how conserved each position is, and its letters", "bottom"));
    await session.step(40, "And user hovers over the \"position 8 of msa(sequence) header\" area of grid", () => hoverArea(page, "position 8 of msa(sequence) header", el("grid")));
    await session.step(41, "Then tooltip should contain text \"Position: 8\"", () => shouldContainText(page, el("tooltip"), "Position: 8"));
    await session.step(42, "When user hovers over the \"position 3 of msa(sequence) header\" area of grid", () => hoverArea(page, "position 3 of msa(sequence) header", el("grid")));
    await session.step(43, "Then tooltip should contain text \"Position: 3\"", () => shouldContainText(page, el("tooltip"), "Position: 3"));
    await session.step(44, "When the video says \"Click a peptide: what differs from it stays in color\" at the bottom", () => videoSaysAt(page, "Click a peptide: what differs from it stays in color", "bottom"));
    await session.step(45, "And user clicks on the \"cell 7 of msa(sequence)\" area of grid", () => clickArea(page, "cell 7 of msa(sequence)", el("grid")));
    await session.step(46, "And user hovers over the \"cell 7 of id\" area of grid", () => hoverArea(page, "cell 7 of id", el("grid")));
    await session.step(47, "And user clicks on the \"cell 13 of msa(sequence)\" area of grid", () => clickArea(page, "cell 13 of msa(sequence)", el("grid")));
    await session.step(48, "And user hovers over the \"cell 13 of id\" area of grid", () => hoverArea(page, "cell 13 of id", el("grid")));
    await session.step(50, "And user picks \"Bio > Transform > Split to Monomers...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Split to Monomers..."), undefined, "Bio > Transform > Split to Monomers");
    await session.step(51, "Then \"Split to Monomers\" dialog should be visible", () => shouldBe(page, el("\"Split to Monomers\" dialog"), "visible"));
    await session.step(53, "When user selects \"msa(sequence)\" in Sequence input in \"Split to Monomers\" dialog", () => selectIn(page, "msa(sequence)", el("Sequence input in \"Split to Monomers\" dialog")), undefined, "Pick the aligned column");
    await session.step(54, "And user clicks on OK button in \"Split to Monomers\" dialog", () => clickOn(page, el("OK button in \"Split to Monomers\" dialog")));
    await session.step(55, "Then a new column \"1\" should have been added", () => newColumnNamed(page, "1"));
    await session.step(56, "And the \"cell type of 10\" reading of grid should be \"Monomer\"", () => readingReads(page, "cell type of 10", el("grid"), "Monomer"));
    await session.step(57, "And every value of \"8\" column should match \"^M$\"", () => everyValueMatches(page, "8", "^M$"));
    await session.step(58, "When the video says \"Every position now has a column of its own\" at the bottom", () => videoSaysAt(page, "Every position now has a column of its own", "bottom"));
    await session.step(60, "And user hovers over the \"cell 3 of 8\" area of grid", () => hoverArea(page, "cell 3 of 8", el("grid")), undefined, "Position 8 holds M in every peptide");
    await session.step(61, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(63, "When user double-clicks on the \"header 14\" area of grid", () => doubleClickArea(page, "header 14", el("grid")), undefined, "Double-click a position's header to sort by it");
    await session.step(64, "Then the \"sort column\" reading of grid should be \"14\"", () => readingReads(page, "sort column", el("grid"), "14"));
    await session.step(65, "When the video says \"Peptides with the same residue at position 14 now sit together\" at the bottom", () => videoSaysAt(page, "Peptides with the same residue at position 14 now sit together", "bottom"));
    await session.step(66, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(67, "And no errors should have been logged", () => noErrors(page));
  });
});
