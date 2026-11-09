/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/07-nucleotides-to-atoms.feature
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
import {noCurrentRow} from '../../bindings/demo-atomic.js';
import {atomicLevelOffCamera, bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {contextPanelClosed, openDatasetRows, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Nucleotides to atoms", () => {
  const session = feature(test, "features/demos/07-nucleotides-to-atoms.feature", import.meta.url);
  test("DNA and RNA, nucleotide by nucleotide", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(17, "And user opens FASTA_DNA dataset keeping the first 12 rows", () => openDatasetRows(page, ds("FASTA_DNA"), 12));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the monomers of the \"sequence\" column are drawn at 14 px", () => sequenceFontSize(page, "sequence", 14));
    await session.step(21, "And the \"sequence\" column of the grid is 390 px wide", () => gridColumnWidth(page, "sequence", 390));
    await session.step(22, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
    await session.step(23, "And user unchecks \"Non-linear\" checkbox in \"To Atomic Level\" dialog", () => uncheck(page, el("\"Non-linear\" checkbox in \"To Atomic Level\" dialog")));
    await session.step(24, "And user checks \"Highlight monomers\" checkbox in \"To Atomic Level\" dialog", () => check(page, el("\"Highlight monomers\" checkbox in \"To Atomic Level\" dialog")));
    await session.step(25, "And user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
    await session.step(26, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(27, "And a new column \"molfile(sequence)\" should have been added", () => newColumnNamed(page, "molfile(sequence)"));
    await session.step(28, "And \"molfile(sequence)\" column should have no missing values", () => columnComplete(page, "molfile(sequence)"));
    await session.step(29, "Given the grid shows only the columns \"molfile(sequence), sequence\"", () => gridShowsColumns(page, "molfile(sequence), sequence"));
    await session.step(30, "And the \"sequence\" column of the grid is 390 px wide", () => gridColumnWidth(page, "sequence", 390));
    await session.step(31, "And the \"molfile(sequence)\" column of the grid is 1110 px wide", () => gridColumnWidth(page, "molfile(sequence)", 1110));
    await session.step(32, "And the rows of the grid are 150 px tall", () => gridRowHeight(page, 150));
    await session.step(33, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(34, "When the video says \"DNA: deoxyribose, base and phosphate per letter\"", () => videoSays(page, "DNA: deoxyribose, base and phosphate per letter"));
    await session.step(35, "And user hovers over monomer 4 of the \"cell 2 of sequence\" area of grid", () => hoverMonomer(page, 4, "cell 2 of sequence", el("grid")));
    await session.step(36, "And user hovers over monomer 13 of the \"cell 2 of sequence\" area of grid", () => hoverMonomer(page, 13, "cell 2 of sequence", el("grid")));
    await session.step(37, "And user hovers over monomer 25 of the \"cell 2 of sequence\" area of grid", () => hoverMonomer(page, 25, "cell 2 of sequence", el("grid")));
    await session.step(38, "And the video says \"RNA: ribose instead of deoxyribose, uracil for thymine\"", () => videoSays(page, "RNA: ribose instead of deoxyribose, uracil for thymine"));
    await session.step(39, "Given user opens FASTA_RNA dataset keeping the first 12 rows", () => openDatasetRows(page, ds("FASTA_RNA"), 12));
    await session.step(40, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(41, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(42, "And the monomers of the \"sequence\" column are drawn at 14 px", () => sequenceFontSize(page, "sequence", 14));
    await session.step(43, "And the \"sequence\" column is converted to atomic level drawn straight with its monomers highlighted", () => atomicLevelOffCamera(page, "sequence"));
    await session.step(44, "Then \"molfile(sequence)\" column should have no missing values", () => columnComplete(page, "molfile(sequence)"));
    await session.step(45, "Given the grid shows only the columns \"molfile(sequence), sequence\"", () => gridShowsColumns(page, "molfile(sequence), sequence"));
    await session.step(46, "And the \"sequence\" column of the grid is 390 px wide", () => gridColumnWidth(page, "sequence", 390));
    await session.step(47, "And the \"molfile(sequence)\" column of the grid is 1110 px wide", () => gridColumnWidth(page, "molfile(sequence)", 1110));
    await session.step(48, "And the rows of the grid are 150 px tall", () => gridRowHeight(page, 150));
    await session.step(49, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(50, "When user hovers over monomer 7 of the \"cell 2 of sequence\" area of grid", () => hoverMonomer(page, 7, "cell 2 of sequence", el("grid")));
    await session.step(51, "And user hovers over monomer 16 of the \"cell 2 of sequence\" area of grid", () => hoverMonomer(page, 16, "cell 2 of sequence", el("grid")));
    await session.step(52, "Then no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
