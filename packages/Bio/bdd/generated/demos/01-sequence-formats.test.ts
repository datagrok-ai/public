/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/01-sequence-formats.feature
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
import {browsePanelNarrower, monomerNameLength} from '../../bindings/demo-notation.js';
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {doubleClickOn, isExpanded, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnUnits, setGridColumnFormat} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {browsePanelOpen, contextPanelClosed, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pointerAway} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("1 · Any notation · Sequence formats", () => {
  const session = feature(test, "features/demos/01-sequence-formats.feature", import.meta.url);
  test("Every notation opens as a Macromolecule column", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(21, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(22, "And Files---App-Data---Bio tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Bio tree node inside browse tree")));
    await session.step(23, "And Files---App-Data---Bio---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Bio---samples tree node inside browse tree")));
    await session.step(24, "And the browse panel is 50 px narrower", () => browsePanelNarrower(page, 50));
    await session.step(25, "When user double-clicks on Files---App-Data---Bio---samples---FASTA.fasta tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---FASTA.fasta tree node inside browse tree")));
    await session.step(26, "Then the \"FASTA\" view should be current", () => viewIsCurrent(page, "FASTA"));
    await session.step(27, "And \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
    await session.step(28, "And \"sequence\" column should have units \"fasta\"", () => columnUnits(page, "sequence", "fasta"));
    await session.step(29, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(30, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(31, "And the \"sequence\" column of the grid is 480 px wide", () => gridColumnWidth(page, "sequence", 480));
    await session.step(32, "And the monomers of the \"sequence\" column are drawn at 14 px", () => sequenceFontSize(page, "sequence", 14));
    await session.step(33, "When the video says \"FASTA: one letter per residue\"", () => videoSays(page, "FASTA: one letter per residue"));
    await session.step(34, "And user hovers over monomer 12 of the \"cell 2 of sequence\" area of grid", () => hoverMonomer(page, 12, "cell 2 of sequence", el("grid")));
    await session.step(35, "And user double-clicks on Files---App-Data---Bio---samples---SEPARATOR_PT.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---SEPARATOR_PT.csv tree node inside browse tree")));
    await session.step(36, "Then the \"SEPARATOR_PT\" view should be current", () => viewIsCurrent(page, "SEPARATOR_PT"));
    await session.step(37, "And \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
    await session.step(38, "And \"sequence\" column should have units \"separator\"", () => columnUnits(page, "sequence", "separator"));
    await session.step(39, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(40, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(41, "And the \"sequence\" column of the grid is 700 px wide", () => gridColumnWidth(page, "sequence", 700));
    await session.step(42, "And the monomers of the \"sequence\" column are drawn at 14 px", () => sequenceFontSize(page, "sequence", 14));
    await session.step(43, "When the video says \"Separator: monomers between dashes\"", () => videoSays(page, "Separator: monomers between dashes"));
    await session.step(44, "And user double-clicks on Files---App-Data---Bio---samples---HELM.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---HELM.csv tree node inside browse tree")));
    await session.step(45, "Then the \"HELM\" view should be current", () => viewIsCurrent(page, "HELM"));
    await session.step(46, "And \"HELM\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM", "Macromolecule"));
    await session.step(47, "And \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
    await session.step(48, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(49, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(50, "And the \"HELM\" column of the grid is 620 px wide", () => gridColumnWidth(page, "HELM", 620));
    await session.step(51, "And user sets the format of \"Activity\" column in the grid to \"0.00\"", () => setGridColumnFormat(page, "Activity", "0.00"));
    await session.step(52, "When the video says \"HELM: non-natural monomers, colored like their natural analogs\"", () => videoSays(page, "HELM: non-natural monomers, colored like their natural analogs"));
    await session.step(53, "And user hovers over monomer 9 of the \"cell 3 of HELM\" area of grid", () => hoverMonomer(page, 9, "cell 3 of HELM", el("grid")));
    await session.step(54, "And user double-clicks on Files---App-Data---Bio---samples---BILN_complex.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---BILN_complex.csv tree node inside browse tree")));
    await session.step(55, "Then the \"BILN_complex\" view should be current", () => viewIsCurrent(page, "BILN_complex"));
    await session.step(56, "And \"biln\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "biln", "Macromolecule"));
    await session.step(57, "And \"biln\" column should have units \"biln\"", () => columnUnits(page, "biln", "biln"));
    await session.step(58, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(59, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(60, "And the grid shows only the columns \"name, category, biln\"", () => gridShowsColumns(page, "name, category, biln"));
    await session.step(61, "And the \"biln\" column of the grid is 620 px wide", () => gridColumnWidth(page, "biln", 620));
    await session.step(62, "And the monomers of the \"biln\" column are drawn at 14 px", () => sequenceFontSize(page, "biln", 14));
    await session.step(63, "And the monomer names of the \"biln\" column are cut at 8 characters", () => monomerNameLength(page, "biln", 8));
    await session.step(64, "And the rows of the grid are 56 px tall", () => gridRowHeight(page, 56));
    await session.step(65, "When the video says \"BILN: rings and cross-links written in the line\"", () => videoSays(page, "BILN: rings and cross-links written in the line"));
    await session.step(66, "And user hovers over monomer 2 of the \"cell 2 of biln\" area of grid", () => hoverMonomer(page, 2, "cell 2 of biln", el("grid")));
    await session.step(67, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(68, "Then no errors should have been logged", () => noErrors(page));
  });
});
