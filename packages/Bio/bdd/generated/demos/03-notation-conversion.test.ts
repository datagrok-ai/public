/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/03-notation-conversion.feature
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
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden, tooltipHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("1 · Any notation · Notation conversion", () => {
  const session = feature(test, "features/demos/03-notation-conversion.feature", import.meta.url);
  test("HELM to BILN and back to plain text", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(21, "And user opens HELM_cyclic dataset", () => openDataset(page, ds("HELM_cyclic")));
    await session.step(22, "Then \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
    await session.step(23, "And \"sequence\" column should have units \"helm\"", () => columnUnits(page, "sequence", "helm"));
    await session.step(24, "Given the browse panel is closed", () => browsePanelClosed(page));
    await session.step(25, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(26, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(27, "And the grid shows only the columns \"sequence_id, sequence\"", () => gridShowsColumns(page, "sequence_id, sequence"));
    await session.step(28, "And the \"sequence\" column of the grid is 300 px wide", () => gridColumnWidth(page, "sequence", 300));
    await session.step(29, "And the rows of the grid are 130 px tall", () => gridRowHeight(page, 130));
    await session.step(30, "When user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."));
    await session.step(31, "Then \"Convert Sequence Notation\" dialog should be visible", () => shouldBe(page, el("\"Convert Sequence Notation\" dialog"), "visible"));
    await session.step(32, "When user selects \"biln\" in \"Convert to\" input in \"Convert Sequence Notation\" dialog", () => selectIn(page, "biln", el("\"Convert to\" input in \"Convert Sequence Notation\" dialog")));
    await session.step(33, "And user clicks on OK button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Sequence Notation\" dialog")));
    await session.step(34, "Then a new column \"biln(sequence)\" should have been added", () => newColumnNamed(page, "biln(sequence)"));
    await session.step(35, "And \"biln(sequence)\" column should have units \"biln\"", () => columnUnits(page, "biln(sequence)", "biln"));
    await session.step(36, "Given the \"biln(sequence)\" column of the grid is 1080 px wide", () => gridColumnWidth(page, "biln(sequence)", 1080));
    await session.step(37, "And the monomers of the \"biln(sequence)\" column are drawn at 14 px", () => sequenceFontSize(page, "biln(sequence)", 14));
    await session.step(38, "And the monomer names of the \"biln(sequence)\" column are cut at 12 characters", () => monomerNameLength(page, "biln(sequence)", 12));
    await session.step(39, "When the video says \"BILN writes the head-to-tail bond into the line\"", () => videoSays(page, "BILN writes the head-to-tail bond into the line"));
    await session.step(40, "And user hovers over monomer 1 of the \"cell 1 of biln(sequence)\" area of grid", () => hoverMonomer(page, 1, "cell 1 of biln(sequence)", el("grid")));
    await session.step(41, "And the video says \"Separator: plain text, without the ring bond\"", () => videoSays(page, "Separator: plain text, without the ring bond"));
    await session.step(42, "Given user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."));
    await session.step(43, "And user selects \"separator\" in \"Convert to\" input in \"Convert Sequence Notation\" dialog", () => selectIn(page, "separator", el("\"Convert to\" input in \"Convert Sequence Notation\" dialog")));
    await session.step(44, "And user clicks on OK button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Sequence Notation\" dialog")));
    await session.step(45, "Then a new column \"separator(sequence)\" should have been added", () => newColumnNamed(page, "separator(sequence)"));
    await session.step(46, "And \"separator(sequence)\" column should have units \"separator\"", () => columnUnits(page, "separator(sequence)", "separator"));
    await session.step(47, "Given the tooltip is hidden", () => tooltipHidden(page));
    await session.step(48, "And the grid shows only the columns \"sequence, biln(sequence), separator(sequence)\"", () => gridShowsColumns(page, "sequence, biln(sequence), separator(sequence)"));
    await session.step(49, "And user sets \"showRowHeader\" property of grid to \"false\"", () => setProperty(page, "showRowHeader", el("grid"), "false"));
    await session.step(50, "And the \"sequence\" column of the grid is 210 px wide", () => gridColumnWidth(page, "sequence", 210));
    await session.step(51, "And the \"biln(sequence)\" column of the grid is 654 px wide", () => gridColumnWidth(page, "biln(sequence)", 654));
    await session.step(52, "And the \"separator(sequence)\" column of the grid is 654 px wide", () => gridColumnWidth(page, "separator(sequence)", 654));
    await session.step(53, "And the monomers of the \"biln(sequence)\" column are drawn at 12 px", () => sequenceFontSize(page, "biln(sequence)", 12));
    await session.step(54, "And the monomers of the \"separator(sequence)\" column are drawn at 12 px", () => sequenceFontSize(page, "separator(sequence)", 12));
    await session.step(55, "And the monomer names of the \"biln(sequence)\" column are cut at 7 characters", () => monomerNameLength(page, "biln(sequence)", 7));
    await session.step(56, "And the monomer names of the \"separator(sequence)\" column are cut at 7 characters", () => monomerNameLength(page, "separator(sequence)", 7));
    await session.step(57, "And the rows of the grid are 150 px tall", () => gridRowHeight(page, 150));
    await session.step(58, "And no row of the table is current", () => noCurrentRow(page));
    await session.step(59, "Then no errors should have been logged", () => noErrors(page));
  });
});
