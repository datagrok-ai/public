/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/08-sirna.feature
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
import {clickPaneMolecule, expandInnerPane, hoverDuplexLinkage, hoverDuplexMonomer, moveLens, paneShowsDrawing} from '../../bindings/demo-atomic.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {expand, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {call} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {balloonsClosed, contextPanelClosed, contextPanelOpen, contextPanelShowsCurrentCell, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · siRNA", () => {
  const session = feature(test, "features/demos/08-sirna.feature", import.meta.url);
  test("siRNA duplexes, every modification in view", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(19, "And user calls \"Bio:demoBioSiRNA\" function", () => call(page, "Bio:demoBioSiRNA"));
    await session.step(20, "And the notification balloons are closed", () => balloonsClosed(page));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the \"oligo_helm\" column of the grid is 640 px wide", () => gridColumnWidth(page, "oligo_helm", 640));
    await session.step(24, "When user hovers over monomer 1 of the duplex in the \"cell 1 of oligo_helm\" area of grid", () => hoverDuplexMonomer(page, 1, "cell 1 of oligo_helm", el("grid")));
    await session.step(25, "And user hovers over monomer 2 of the duplex in the \"cell 1 of oligo_helm\" area of grid", () => hoverDuplexMonomer(page, 2, "cell 1 of oligo_helm", el("grid")));
    await session.step(26, "And user hovers over the linkage after monomer 2 of the duplex in the \"cell 1 of oligo_helm\" area of grid", () => hoverDuplexLinkage(page, 2, "cell 1 of oligo_helm", el("grid")));
    await session.step(27, "And user hovers over monomer 20 of the duplex in the \"cell 1 of oligo_helm\" area of grid", () => hoverDuplexMonomer(page, 20, "cell 1 of oligo_helm", el("grid")));
    await session.step(28, "And the video says \"Every modification counted, each strand drawn atom by atom\"", () => videoSays(page, "Every modification counted, each strand drawn atom by atom"));
    await session.step(29, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(30, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(31, "When user expands \"Oligo-Nucleotide\" section", () => expand(page, el("\"Oligo-Nucleotide\" section")));
    await session.step(32, "And user expands \"Oligo Structures\" section", () => expand(page, el("\"Oligo Structures\" section")));
    await session.step(33, "And user expands the \"Sense\" pane inside the \"Oligo Structures\" section", () => expandInnerPane(page, "Sense", "Oligo Structures"));
    await session.step(34, "Then the \"Sense\" section should show a drawn molecule", () => paneShowsDrawing(page, "Sense"));
    await session.step(35, "When user clicks on the molecule in the \"Sense\" section", () => clickPaneMolecule(page, "Sense"));
    await session.step(36, "And user hovers over the molecule in the \"Molecule\" dialog at 30% across and 45% down", () => moveLens(page, "Molecule", 30, 45));
    await session.step(37, "And user hovers over the molecule in the \"Molecule\" dialog at 73% across and 38% down", () => moveLens(page, "Molecule", 73, 38));
    await session.step(38, "And user hovers over the molecule in the \"Molecule\" dialog at 69% across and 59% down", () => moveLens(page, "Molecule", 69, 59));
    await session.step(39, "Then no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
