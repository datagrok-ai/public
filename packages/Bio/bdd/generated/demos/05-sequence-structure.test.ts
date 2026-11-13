/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/05-sequence-structure.feature
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
import {clickPaneMolecule, moveLens, paneShowsDrawing, turn3d} from '../../bindings/demo-atomic.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {close, expand, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelOpen, contextPanelShowsCurrentCell, openDatasetRows, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Sequence structure", () => {
  const session = feature(test, "features/demos/05-sequence-structure.feature", import.meta.url);
  test("A sequence's structure in 2D and 3D", {tag: ["@guide", "@demo", "@viewport:1440x810"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(17, "And user opens HELM_cyclic dataset keeping the first 8 rows", () => openDatasetRows(page, ds("HELM_cyclic"), 8));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(20, "And the grid shows only the columns \"sequence_id, sequence, activity\"", () => gridShowsColumns(page, "sequence_id, sequence, activity"));
    await session.step(21, "And the \"sequence\" column of the grid is 400 px wide", () => gridColumnWidth(page, "sequence", 400));
    await session.step(22, "And the rows of the grid are 190 px tall", () => gridRowHeight(page, 190));
    await session.step(23, "When user clicks on the \"cell 2 of sequence\" area of grid", () => clickArea(page, "cell 2 of sequence", el("grid")));
    await session.step(24, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(25, "When user expands \"Molecular Structure\" section", () => expand(page, el("\"Molecular Structure\" section")));
    await session.step(26, "Then the \"Molecular Structure\" section should show a drawn molecule", () => paneShowsDrawing(page, "Molecular Structure"));
    await session.step(27, "When the video says \"Built atom by atom from the monomer library\"", () => videoSays(page, "Built atom by atom from the monomer library"));
    await session.step(28, "And user clicks on the molecule in the \"Molecular Structure\" section", () => clickPaneMolecule(page, "Molecular Structure"));
    await session.step(29, "And user hovers over the molecule in the \"Molecule\" dialog at 32% across and 30% down", () => moveLens(page, "Molecule", 32, 30));
    await session.step(30, "And user hovers over the molecule in the \"Molecule\" dialog at 62% across and 28% down", () => moveLens(page, "Molecule", 62, 28));
    await session.step(31, "And user hovers over the molecule in the \"Molecule\" dialog at 58% across and 75% down", () => moveLens(page, "Molecule", 58, 75));
    await session.step(32, "And user closes \"Molecule\" dialog", () => close(page, el("\"Molecule\" dialog")));
    await session.step(33, "And user expands \"Molecular 3D Structure\" section", () => expand(page, el("\"Molecular 3D Structure\" section")));
    await session.step(34, "Then the \"Molecular 3D Structure\" section should show a drawn molecule", () => paneShowsDrawing(page, "Molecular 3D Structure"));
    await session.step(35, "When the video says \"3D coordinates computed in the browser\"", () => videoSays(page, "3D coordinates computed in the browser"));
    await session.step(36, "And user turns the 3D structure in the \"Molecular 3D Structure\" section", () => turn3d(page, "Molecular 3D Structure"));
    await session.step(37, "Then no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
