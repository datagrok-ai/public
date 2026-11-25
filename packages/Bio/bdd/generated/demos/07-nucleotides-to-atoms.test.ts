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
import {regionExtracted} from '../../bindings/demo-chem.js';
import {atomicLevelOffCamera, bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDatasetRows, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("2 · Every sequence is a molecule · Nucleotides to atoms", () => {
  const session = feature(test, "features/demos/07-nucleotides-to-atoms.feature", import.meta.url);
  test("RNA and DNA, nucleotide by nucleotide", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(18, "And user opens FASTA_RNA dataset keeping the first 8 rows", () => openDatasetRows(page, ds("FASTA_RNA"), 8));
    await session.step(19, "And the region 1 to 15 of the \"sequence\" column is extracted as \"strand\"", () => regionExtracted(page, 1, 15, "sequence", "strand"));
    await session.step(20, "And user removes \"sequence\" column", () => removeColumn(page, "sequence"));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the monomers of the \"strand\" column are drawn at 16 px", () => sequenceFontSize(page, "strand", 16));
    await session.step(25, "And the \"strand\" column of the grid is 260 px wide", () => gridColumnWidth(page, "strand", 260));
    await session.step(26, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
    await session.step(27, "And user unchecks \"Non-linear\" checkbox in \"To Atomic Level\" dialog", () => uncheck(page, el("\"Non-linear\" checkbox in \"To Atomic Level\" dialog")));
    await session.step(28, "And user checks \"Highlight monomers\" checkbox in \"To Atomic Level\" dialog", () => check(page, el("\"Highlight monomers\" checkbox in \"To Atomic Level\" dialog")));
    await session.step(29, "And user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
    await session.step(30, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(31, "And a new column \"molfile(strand)\" should have been added", () => newColumnNamed(page, "molfile(strand)"));
    await session.step(32, "And \"molfile(strand)\" column should have no missing values", () => columnComplete(page, "molfile(strand)"));
    await session.step(33, "Given the grid shows only the columns \"molfile(strand), strand\"", () => gridShowsColumns(page, "molfile(strand), strand"));
    await session.step(34, "And the \"strand\" column of the grid is 230 px wide", () => gridColumnWidth(page, "strand", 230));
    await session.step(35, "And the \"molfile(strand)\" column of the grid is 1220 px wide", () => gridColumnWidth(page, "molfile(strand)", 1220));
    await session.step(36, "And the rows of the grid are 200 px tall", () => gridRowHeight(page, 200));
    await session.step(37, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(38, "When the video says \"RNA: ribose, base and phosphate per letter\"", () => videoSays(page, "RNA: ribose, base and phosphate per letter"));
    await session.step(39, "And user hovers over monomer 3 of the \"cell 2 of strand\" area of grid", () => hoverMonomer(page, 3, "cell 2 of strand", el("grid")));
    await session.step(40, "And user hovers over monomer 7 of the \"cell 2 of strand\" area of grid", () => hoverMonomer(page, 7, "cell 2 of strand", el("grid")));
    await session.step(41, "Given user opens FASTA_DNA dataset keeping the first 8 rows", () => openDatasetRows(page, ds("FASTA_DNA"), 8));
    await session.step(42, "And the region 1 to 15 of the \"sequence\" column is extracted as \"strand\"", () => regionExtracted(page, 1, 15, "sequence", "strand"));
    await session.step(43, "And user removes \"sequence\" column", () => removeColumn(page, "sequence"));
    await session.step(44, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(45, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(46, "And the monomers of the \"strand\" column are drawn at 16 px", () => sequenceFontSize(page, "strand", 16));
    await session.step(47, "And the \"strand\" column is converted to atomic level drawn straight with its monomers highlighted", () => atomicLevelOffCamera(page, "strand"));
    await session.step(48, "Then \"molfile(strand)\" column should have no missing values", () => columnComplete(page, "molfile(strand)"));
    await session.step(49, "Given the grid shows only the columns \"molfile(strand), strand\"", () => gridShowsColumns(page, "molfile(strand), strand"));
    await session.step(50, "And the \"strand\" column of the grid is 230 px wide", () => gridColumnWidth(page, "strand", 230));
    await session.step(51, "And the \"molfile(strand)\" column of the grid is 1220 px wide", () => gridColumnWidth(page, "molfile(strand)", 1220));
    await session.step(52, "And the rows of the grid are 200 px tall", () => gridRowHeight(page, 200));
    await session.step(53, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(55, "Then \"strand\" column should have tag \"alphabet\" equal to \"DNA\"", () => columnTag(page, "strand", "alphabet", "DNA"), undefined, "DNA: deoxyribose instead of ribose, thymine for uracil");
  });
});
