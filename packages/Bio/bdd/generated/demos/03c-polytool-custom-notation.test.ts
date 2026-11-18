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
import {noCurrentRow} from '../../bindings/demo-atomic.js';
import {monomerNameLength} from '../../bindings/demo-notation.js';
import {gridScrolledToRow} from '../../bindings/demo-search.js';
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, contextPanelClosed, dialogCloses, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pointerAway} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("1 · Any notation · Custom notation with PolyTool", () => {
  const session = feature(test, "features/demos/03c-polytool-custom-notation.feature", import.meta.url);
  test("Custom notation to HELM and molecules", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And simple mode is off", () => simpleModeOff(page));
    await session.step(22, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(23, "And user opens cyclized dataset", () => openDataset(page, ds("cyclized")));
    await session.step(24, "Then \"seqs\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "seqs", "Macromolecule"));
    await session.step(25, "Given the browse panel is closed", () => browsePanelClosed(page));
    await session.step(26, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(27, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(28, "And the \"seqs\" column of the grid is 520 px wide", () => gridColumnWidth(page, "seqs", 520));
    await session.step(29, "And the monomers of the \"seqs\" column are drawn at 15 px", () => sequenceFontSize(page, "seqs", 15));
    await session.step(30, "And the monomer names of the \"seqs\" column are cut at 8 characters", () => monomerNameLength(page, "seqs", 8));
    await session.step(31, "And the rows of the grid are 150 px tall", () => gridRowHeight(page, 150));
    await session.step(32, "When user picks \"Bio > PolyTool > Convert...\" from the top menu", () => pickFromTopMenu(page, "Bio > PolyTool > Convert..."));
    await session.step(33, "Then \"PolyTool Conversion\" dialog should be visible", () => shouldBe(page, el("\"PolyTool Conversion\" dialog"), "visible"));
    await session.step(34, "When user clicks on OK button in \"PolyTool Conversion\" dialog", () => clickOn(page, el("OK button in \"PolyTool Conversion\" dialog")));
    await session.step(35, "Then the \"PolyTool Conversion\" dialog should close", () => dialogCloses(page, "PolyTool Conversion"));
    await session.step(36, "And a new column \"transformed(seqs)\" should have been added", () => newColumnNamed(page, "transformed(seqs)"));
    await session.step(37, "And a new column \"molfile(seqs)\" should have been added", () => newColumnNamed(page, "molfile(seqs)"));
    await session.step(38, "Given the grid shows only the columns \"seqs, transformed(seqs), molfile(seqs)\"", () => gridShowsColumns(page, "seqs, transformed(seqs), molfile(seqs)"));
    await session.step(39, "And the \"seqs\" column of the grid is 440 px wide", () => gridColumnWidth(page, "seqs", 440));
    await session.step(40, "And the \"transformed(seqs)\" column of the grid is 320 px wide", () => gridColumnWidth(page, "transformed(seqs)", 320));
    await session.step(41, "And the \"molfile(seqs)\" column of the grid is 720 px wide", () => gridColumnWidth(page, "molfile(seqs)", 720));
    await session.step(42, "And the rows of the grid are 300 px tall", () => gridRowHeight(page, 300));
    await session.step(43, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(44, "When the video says \"C(1)…C(1) closures become a HELM ring and macrocycle\"", () => videoSays(page, "C(1)…C(1) closures become a HELM ring and macrocycle"));
    await session.step(45, "And user hovers over monomer 3 of the \"cell 1 of seqs\" area of grid", () => hoverMonomer(page, 3, "cell 1 of seqs", el("grid")));
    await session.step(46, "And user hovers over monomer 1 of the \"cell 2 of seqs\" area of grid", () => hoverMonomer(page, 1, "cell 2 of seqs", el("grid")));
    await session.step(47, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(48, "And the video says \"Click chemistry: azide plus alkyne closes a triazole\"", () => videoSays(page, "Click chemistry: azide plus alkyne closes a triazole"));
    await session.step(49, "Given the rows of the grid are 370 px tall", () => gridRowHeight(page, 370));
    await session.step(50, "And the grid is scrolled down to row 7", () => gridScrolledToRow(page, 7));
    await session.step(51, "Then no errors should have been logged", () => noErrors(page));
  });
});
