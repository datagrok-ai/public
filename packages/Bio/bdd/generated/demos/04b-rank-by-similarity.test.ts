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
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns, scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Similarity to a lead", () => {
  const session = feature(test, "features/demos/04b-rank-by-similarity.feature", import.meta.url);
  test("Rank the peptides by their similarity to a lead", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(15, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(16, "Then \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
    await session.step(17, "Given the grid shows only the columns \"sequence_id, sequence, activity\"", () => gridShowsColumns(page, "sequence_id, sequence, activity"));
    await session.step(18, "And the \"sequence\" column of the grid is 300 px wide", () => gridColumnWidth(page, "sequence", 300));
    await session.step(19, "And the monomers of the \"sequence\" column are drawn at 16 px", () => sequenceFontSize(page, "sequence", 16));
    await session.step(20, "When user picks \"Bio > Calculate > Similarity...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Similarity..."));
    await session.step(21, "Then Similarity dialog should be visible", () => shouldBe(page, el("Similarity dialog"), "visible"));
    await session.step(22, "When user enters \"DLKMMHQNTPGNMAQD\" into Reference input in Similarity dialog", () => enterInto(page, "DLKMMHQNTPGNMAQD", el("Reference input in Similarity dialog")));
    await session.step(23, "And user clicks on OK button in Similarity dialog", () => clickOn(page, el("OK button in Similarity dialog")));
    await session.step(24, "Then a new column \"Similarity\" should have been added", () => newColumnNamed(page, "Similarity"));
    await session.step(25, "And \"Similarity\" column should have its maximum in row 7", () => maxInRow(page, "Similarity", 7));
    await session.step(26, "When user double-clicks on the \"header Similarity\" area of grid", () => doubleClickArea(page, "header Similarity", el("grid")));
    await session.step(27, "Given user scrolls the grid to the \"sequence_id\" column", () => scrollGridTo(page, "sequence_id"));
    await session.step(28, "When the video says \"The lead first, then its closest analogs\"", () => videoSays(page, "The lead first, then its closest analogs"));
    await session.step(29, "And user clicks on the \"cell 7 of sequence\" area of grid", () => clickArea(page, "cell 7 of sequence", el("grid")));
    await session.step(30, "And user hovers over the \"cell 7 of sequence_id\" area of grid", () => hoverArea(page, "cell 7 of sequence_id", el("grid")));
    await session.step(31, "And the video says \"One substitution apart: activity 4.93 against 1.19\"", () => videoSays(page, "One substitution apart: activity 4.93 against 1.19"));
    await session.step(32, "Then no errors should have been logged", () => noErrors(page));
  });
});
