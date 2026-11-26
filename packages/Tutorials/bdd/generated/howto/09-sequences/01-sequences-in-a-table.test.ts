/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/09-sequences/01-sequences-in-a-table.feature
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
import {bioReady} from '../../../bindings/howto-sequences.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, expand, isExpanded, selectIn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, contextPanelClosed, contextPanelShows, contextPanelShowsCurrentCell, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, dragAreaBy, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Sequences", () => {
  const session = feature(test, "features/howto/09-sequences/01-sequences-in-a-table.feature", import.meta.url);
  test("Open, inspect and convert sequences", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And the Bio package is initialized", () => bioReady(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(22, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(23, "And Files---App-Data---Bio tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Bio tree node inside browse tree")));
    await session.step(24, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(26, "When user expands Files---App-Data---Bio---samples tree node inside browse tree", () => expand(page, el("Files---App-Data---Bio---samples tree node inside browse tree")), undefined, "Sample sequences: Files > App Data > Bio > samples");
    await session.step(28, "And user double-clicks on Files---App-Data---Bio---samples---FASTA.fasta tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---FASTA.fasta tree node inside browse tree")), undefined, "A FASTA file opens as a table");
    await session.step(29, "Then the \"FASTA\" view should be current", () => viewIsCurrent(page, "FASTA"));
    await session.step(30, "And \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
    await session.step(31, "And \"sequence\" column should have units \"fasta\"", () => columnUnits(page, "sequence", "fasta"));
    await session.step(32, "When the video says \"Recognized as sequences: every residue in its color\"", () => videoSays(page, "Recognized as sequences: every residue in its color"));
    await session.step(34, "And user drags the \"column resizer sequence\" area of grid by 200 pixels to the right", () => dragAreaBy(page, "column resizer sequence", el("grid"), 200, "right"), undefined, "Drag the column's edge to see more of each sequence");
    await session.step(35, "And user hovers over monomer 6 of the \"cell 2 of sequence\" area of grid", () => hoverMonomer(page, 6, "cell 2 of sequence", el("grid")));
    await session.step(37, "And user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "The context panel: what Datagrok knows about the column");
    await session.step(38, "And user clicks on the \"header sequence\" area of grid", () => clickArea(page, "header sequence", el("grid")));
    await session.step(39, "Then the context panel should show \"sequence\"", () => contextPanelShows(page, "sequence"));
    await session.step(41, "When user clicks on \"Details\" accordion header in context panel", () => clickOn(page, el("\"Details\" accordion header in context panel")), undefined, "Details: FASTA notation, an alphabet of peptides");
    await session.step(43, "And user clicks on the \"cell 3 of sequence\" area of grid", () => clickArea(page, "cell 3 of sequence", el("grid")), undefined, "Click a sequence for panes about that one");
    await session.step(44, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(46, "When user clicks on \"Composition analysis\" accordion header in context panel", () => clickOn(page, el("\"Composition analysis\" accordion header in context panel")), undefined, "Composition: how often each residue occurs");
    await session.step(48, "And user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "Close the panel to make room");
    await session.step(50, "And user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."), undefined, "Bio > Transform > Convert Sequence Notation");
    await session.step(51, "Then \"Convert Sequence Notation\" dialog should be visible", () => shouldBe(page, el("\"Convert Sequence Notation\" dialog"), "visible"));
    await session.step(53, "When user selects \"separator\" in \"Convert to\" input in \"Convert Sequence Notation\" dialog", () => selectIn(page, "separator", el("\"Convert to\" input in \"Convert Sequence Notation\" dialog")), undefined, "Pick a notation: here, separator text");
    await session.step(54, "And user clicks on OK button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Sequence Notation\" dialog")));
    await session.step(55, "Then a new column \"separator(sequence)\" should have been added", () => newColumnNamed(page, "separator(sequence)"));
    await session.step(56, "And \"separator(sequence)\" column should have units \"separator\"", () => columnUnits(page, "separator(sequence)", "separator"));
    await session.step(57, "When the video says \"A new column: the same residues, with dashes between them\"", () => videoSays(page, "A new column: the same residues, with dashes between them"));
    await session.step(58, "And user hovers over monomer 6 of the \"cell 2 of separator(sequence)\" area of grid", () => hoverMonomer(page, 6, "cell 2 of separator(sequence)", el("grid")));
    await session.step(59, "Then no errors should have been logged", () => noErrors(page));
  });
});
