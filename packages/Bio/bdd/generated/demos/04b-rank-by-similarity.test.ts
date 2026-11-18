/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/04b-rank-by-similarity.feature
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
import {clickOn, enterInto, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, maxInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, hoverArea, noErrors, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridHeaderHeight, gridRowHeight, gridShowsColumns, scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("4 · A peptide series · Similarity to a lead", () => {
  const session = feature(test, "features/demos/04b-rank-by-similarity.feature", import.meta.url);
  test("Rank the peptides by their similarity to a lead", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(17, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(18, "Then \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
    await session.step(19, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the grid shows only the columns \"sequence_id, sequence, activity\"", () => gridShowsColumns(page, "sequence_id, sequence, activity"));
    await session.step(22, "And user sets \"defaultCellFont\" property of grid to 'normal normal 18px \"Roboto\"'", () => setProperty(page, "defaultCellFont", el("grid"), "normal normal 18px \"Roboto\""));
    await session.step(23, "And user sets \"colHeaderFont\" property of grid to 'normal bold 16px \"Roboto\"'", () => setProperty(page, "colHeaderFont", el("grid"), "normal bold 16px \"Roboto\""));
    await session.step(24, "And the column headers of the grid are 40 px tall", () => gridHeaderHeight(page, 40));
    await session.step(25, "And the rows of the grid are 40 px tall", () => gridRowHeight(page, 40));
    await session.step(26, "And the monomers of the \"sequence\" column are drawn at 28 px", () => sequenceFontSize(page, "sequence", 28));
    await session.step(27, "And the \"sequence_id\" column of the grid is 180 px wide", () => gridColumnWidth(page, "sequence_id", 180));
    await session.step(28, "And the \"sequence\" column of the grid is 640 px wide", () => gridColumnWidth(page, "sequence", 640));
    await session.step(29, "And the \"activity\" column of the grid is 140 px wide", () => gridColumnWidth(page, "activity", 140));
    await session.step(30, "When user picks \"Bio > Calculate > Similarity...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Similarity..."));
    await session.step(31, "Then Similarity dialog should be visible", () => shouldBe(page, el("Similarity dialog"), "visible"));
    await session.step(32, "When user enters \"DLKMMHQNTPGNMAQD\" into Reference input in Similarity dialog", () => enterInto(page, "DLKMMHQNTPGNMAQD", el("Reference input in Similarity dialog")));
    await session.step(33, "And user clicks on OK button in Similarity dialog", () => clickOn(page, el("OK button in Similarity dialog")));
    await session.step(34, "Then a new column \"Similarity\" should have been added", () => newColumnNamed(page, "Similarity"));
    await session.step(35, "And \"Similarity\" column should have its maximum in row 7", () => maxInRow(page, "Similarity", 7));
    await session.step(36, "Given the \"Similarity\" column of the grid is 160 px wide", () => gridColumnWidth(page, "Similarity", 160));
    await session.step(37, "When user double-clicks on the \"header Similarity\" area of grid", () => doubleClickArea(page, "header Similarity", el("grid")));
    await session.step(38, "Given user scrolls the grid to the \"sequence_id\" column", () => scrollGridTo(page, "sequence_id"));
    await session.step(39, "When user hovers over the \"cell 7 of sequence_id\" area of grid", () => hoverArea(page, "cell 7 of sequence_id", el("grid")));
    await session.step(40, "And the video says \"The lead first, then its closest analogs\"", () => videoSays(page, "The lead first, then its closest analogs"));
    await session.step(41, "And user clicks on the \"cell 7 of sequence\" area of grid", () => clickArea(page, "cell 7 of sequence", el("grid")));
    await session.step(42, "And user hovers over the \"cell 8 of activity\" area of grid", () => hoverArea(page, "cell 8 of activity", el("grid")));
    await session.step(43, "And the video says \"One substitution apart: activity 4.93 against 1.19\"", () => videoSays(page, "One substitution apart: activity 4.93 against 1.19"));
    await session.step(44, "Then no errors should have been logged", () => noErrors(page));
  });
});
