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
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {doubleClickOn, isExpanded, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {browsePanelOpen, contextPanelClosed, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Sequence formats", () => {
  const session = feature(test, "features/demos/01-sequence-formats.feature", import.meta.url);
  test("Every notation opens as a Macromolecule column", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(14, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(15, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(16, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(17, "And Files---App-Data---Bio tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Bio tree node inside browse tree")));
    await session.step(18, "And Files---App-Data---Bio---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Bio---samples tree node inside browse tree")));
    await session.step(19, "When user double-clicks on Files---App-Data---Bio---samples---FASTA.fasta tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---FASTA.fasta tree node inside browse tree")));
    await session.step(20, "Then the \"FASTA\" view should be current", () => viewIsCurrent(page, "FASTA"));
    await session.step(21, "And \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
    await session.step(22, "And \"sequence\" column should have units \"fasta\"", () => columnUnits(page, "sequence", "fasta"));
    await session.step(23, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(24, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(25, "And the \"sequence\" column of the grid is 420 px wide", () => gridColumnWidth(page, "sequence", 420));
    await session.step(26, "When the video says \"FASTA: one letter per residue\"", () => videoSays(page, "FASTA: one letter per residue"));
    await session.step(27, "And user hovers over monomer 12 of the \"cell 2 of sequence\" area of grid", () => hoverMonomer(page, 12, "cell 2 of sequence", el("grid")));
    await session.step(28, "And user double-clicks on Files---App-Data---Bio---samples---SEPARATOR_PT.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---SEPARATOR_PT.csv tree node inside browse tree")));
    await session.step(29, "Then the \"SEPARATOR_PT\" view should be current", () => viewIsCurrent(page, "SEPARATOR_PT"));
    await session.step(30, "And \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
    await session.step(31, "And \"sequence\" column should have units \"separator\"", () => columnUnits(page, "sequence", "separator"));
    await session.step(32, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(33, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(34, "And the \"sequence\" column of the grid is 620 px wide", () => gridColumnWidth(page, "sequence", 620));
    await session.step(35, "When the video says \"Separator: monomers between dashes\"", () => videoSays(page, "Separator: monomers between dashes"));
    await session.step(36, "And user double-clicks on Files---App-Data---Bio---samples---HELM.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---HELM.csv tree node inside browse tree")));
    await session.step(37, "Then the \"HELM\" view should be current", () => viewIsCurrent(page, "HELM"));
    await session.step(38, "And \"HELM\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM", "Macromolecule"));
    await session.step(39, "And \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
    await session.step(40, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(41, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(42, "And the \"HELM\" column of the grid is 620 px wide", () => gridColumnWidth(page, "HELM", 620));
    await session.step(43, "When the video says \"HELM: non-natural monomers, colored like their natural analogs\"", () => videoSays(page, "HELM: non-natural monomers, colored like their natural analogs"));
    await session.step(44, "And user hovers over monomer 9 of the \"cell 3 of HELM\" area of grid", () => hoverMonomer(page, 9, "cell 3 of HELM", el("grid")));
    await session.step(45, "And user double-clicks on Files---App-Data---Bio---samples---BILN.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---BILN.csv tree node inside browse tree")));
    await session.step(46, "Then the \"BILN\" view should be current", () => viewIsCurrent(page, "BILN"));
    await session.step(47, "And \"biln\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "biln", "Macromolecule"));
    await session.step(48, "And \"biln\" column should have units \"biln\"", () => columnUnits(page, "biln", "biln"));
    await session.step(49, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(50, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(51, "And the \"biln\" column of the grid is 620 px wide", () => gridColumnWidth(page, "biln", 620));
    await session.step(52, "When the video says \"BILN: chains and ring bonds in one line\"", () => videoSays(page, "BILN: chains and ring bonds in one line"));
    await session.step(53, "And user double-clicks on Files---App-Data---Bio---samples---MSA.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Bio---samples---MSA.csv tree node inside browse tree")));
    await session.step(54, "Then the \"MSA\" view should be current", () => viewIsCurrent(page, "MSA"));
    await session.step(55, "And \"MSA\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "MSA", "Macromolecule"));
    await session.step(56, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(57, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(58, "And the \"MSA\" column of the grid is 620 px wide", () => gridColumnWidth(page, "MSA", 620));
    await session.step(59, "When the video says \"Aligned: monomers line up position by position\"", () => videoSays(page, "Aligned: monomers line up position by position"));
    await session.step(60, "And user hovers over monomer 10 of the \"cell 4 of MSA\" area of grid", () => hoverMonomer(page, 10, "cell 4 of MSA", el("grid")));
    await session.step(61, "Then no errors should have been logged", () => noErrors(page));
  });
});
