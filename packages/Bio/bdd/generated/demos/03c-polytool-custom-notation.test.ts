/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/03c-polytool-custom-notation.feature
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
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {dialogCloses, openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Custom notation with PolyTool", () => {
  const session = feature(test, "features/demos/03c-polytool-custom-notation.feature", import.meta.url);
  test("Custom notation to HELM and molecules", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(14, "And user opens cyclized dataset", () => openDataset(page, ds("cyclized")));
    await session.step(15, "Then \"seqs\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "seqs", "Macromolecule"));
    await session.step(16, "Given the \"seqs\" column of the grid is 520 px wide", () => gridColumnWidth(page, "seqs", 520));
    await session.step(17, "And the monomers of the \"seqs\" column are drawn at 15 px", () => sequenceFontSize(page, "seqs", 15));
    await session.step(18, "And the rows of the grid are 150 px tall", () => gridRowHeight(page, 150));
    await session.step(19, "When user picks \"Bio > PolyTool > Convert...\" from the top menu", () => pickFromTopMenu(page, "Bio > PolyTool > Convert..."));
    await session.step(20, "Then \"PolyTool Conversion\" dialog should be visible", () => shouldBe(page, el("\"PolyTool Conversion\" dialog"), "visible"));
    await session.step(21, "When user clicks on OK button in \"PolyTool Conversion\" dialog", () => clickOn(page, el("OK button in \"PolyTool Conversion\" dialog")));
    await session.step(22, "Then the \"PolyTool Conversion\" dialog should close", () => dialogCloses(page, "PolyTool Conversion"));
    await session.step(23, "And a new column \"transformed(seqs)\" should have been added", () => newColumnNamed(page, "transformed(seqs)"));
    await session.step(24, "And a new column \"molfile(seqs)\" should have been added", () => newColumnNamed(page, "molfile(seqs)"));
    await session.step(25, "Given the \"seqs\" column of the grid is 420 px wide", () => gridColumnWidth(page, "seqs", 420));
    await session.step(26, "And the \"transformed(seqs)\" column of the grid is 300 px wide", () => gridColumnWidth(page, "transformed(seqs)", 300));
    await session.step(27, "And the \"molfile(seqs)\" column of the grid is 340 px wide", () => gridColumnWidth(page, "molfile(seqs)", 340));
    await session.step(28, "When the video says \"Ring closures C(1)…C(1) become the HELM ring and the molecule's macrocycle\"", () => videoSays(page, "Ring closures C(1)…C(1) become the HELM ring and the molecule's macrocycle"));
    await session.step(29, "And user hovers over monomer 3 of the \"cell 1 of transformed(seqs)\" area of grid", () => hoverMonomer(page, 3, "cell 1 of transformed(seqs)", el("grid")));
    await session.step(30, "And user hovers over the \"cell 4 of n\" area of grid", () => hoverArea(page, "cell 4 of n", el("grid")));
    await session.step(31, "Then no errors should have been logged", () => noErrors(page));
  });
});
