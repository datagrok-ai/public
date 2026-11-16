/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/06d-monomers-as-molecules.feature
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
import {clickOn, expand, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelClosed, collapseInnerPane, contextPanelClosed, contextPanelOpen, contextPanelShowsCurrentCell, expandInnerPane, openDatasetRows, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("2 · Every sequence is a molecule · Monomers as molecules", () => {
  const session = feature(test, "features/demos/06d-monomers-as-molecules.feature", import.meta.url);
  test("Monomers are molecules too", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(17, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(18, "And user opens HELM_cyclic dataset keeping the first 8 rows", () => openDatasetRows(page, ds("HELM_cyclic"), 8));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the grid shows only the columns \"sequence_id, sequence\"", () => gridShowsColumns(page, "sequence_id, sequence"));
    await session.step(23, "And the \"sequence\" column of the grid is 260 px wide", () => gridColumnWidth(page, "sequence", 260));
    await session.step(24, "And the rows of the grid are 110 px tall", () => gridRowHeight(page, 110));
    await session.step(25, "When user picks \"Bio > Transform > Split to Monomers...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Split to Monomers..."));
    await session.step(26, "Then \"Split to Monomers\" dialog should be visible", () => shouldBe(page, el("\"Split to Monomers\" dialog"), "visible"));
    await session.step(27, "When user clicks on OK button in \"Split to Monomers\" dialog", () => clickOn(page, el("OK button in \"Split to Monomers\" dialog")));
    await session.step(28, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(29, "And 16 new columns should have been added", () => newColumnsCount(page, 16));
    await session.step(30, "And \"7\" column should have semantic type \"Monomer\"", () => columnSemType(page, "7", "Monomer"));
    await session.step(31, "And the value of \"7\" column in row 1 should be \"Phe_34diCl\"", () => valueInRow(page, "7", 1, "Phe_34diCl"));
    await session.step(32, "Given the grid shows only the columns \"sequence_id, 1, 2, 3, 4, 5, 6, 7, 8\"", () => gridShowsColumns(page, "sequence_id, 1, 2, 3, 4, 5, 6, 7, 8"));
    await session.step(33, "And the \"1\" column of the grid is 96 px wide", () => gridColumnWidth(page, "1", 96));
    await session.step(34, "And the \"2\" column of the grid is 96 px wide", () => gridColumnWidth(page, "2", 96));
    await session.step(35, "And the \"3\" column of the grid is 96 px wide", () => gridColumnWidth(page, "3", 96));
    await session.step(36, "And the \"4\" column of the grid is 96 px wide", () => gridColumnWidth(page, "4", 96));
    await session.step(37, "And the \"5\" column of the grid is 96 px wide", () => gridColumnWidth(page, "5", 96));
    await session.step(38, "And the \"6\" column of the grid is 96 px wide", () => gridColumnWidth(page, "6", 96));
    await session.step(39, "And the \"7\" column of the grid is 96 px wide", () => gridColumnWidth(page, "7", 96));
    await session.step(40, "And the \"8\" column of the grid is 96 px wide", () => gridColumnWidth(page, "8", 96));
    await session.step(41, "And the rows of the grid are 36 px tall", () => gridRowHeight(page, 36));
    await session.step(42, "When the video says \"One column per position, a monomer per cell\"", () => videoSays(page, "One column per position, a monomer per cell"));
    await session.step(43, "And user clicks on the \"cell 1 of 7\" area of grid", () => clickArea(page, "cell 1 of 7", el("grid")));
    await session.step(44, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(45, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(46, "When user expands \"Monomer\" section", () => expand(page, el("\"Monomer\" section")));
    await session.step(47, "And the video says \"A monomer is a molecule: every chemistry panel applies\"", () => videoSays(page, "A monomer is a molecule: every chemistry panel applies"));
    await session.step(48, "And user collapses the \"Details\" pane inside the \"Monomer\" section", () => collapseInnerPane(page, "Details", "Monomer"));
    await session.step(49, "And user expands the \"Chemistry\" pane inside the \"Molecule\" section", () => expandInnerPane(page, "Chemistry", "Molecule"));
    await session.step(50, "And user expands the \"Properties\" pane inside the \"Chemistry\" section", () => expandInnerPane(page, "Properties", "Chemistry"));
    await session.step(51, "Then no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
