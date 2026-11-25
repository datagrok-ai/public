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
import {contextPanelWider} from '../../bindings/demo-chem.js';
import {monomerNameLength} from '../../bindings/demo-notation.js';
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {close, collapse, expand, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelOpen, contextPanelShowsCurrentCell, expandInnerPane, openDatasetRows, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("2 · Every sequence is a molecule · Sequence structure", () => {
  const session = feature(test, "features/demos/05-sequence-structure.feature", import.meta.url);
  test("A peptide drug in 2D and 3D", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(21, "And user opens BILN_complex dataset keeping the first 16 rows", () => openDatasetRows(page, ds("BILN_complex"), 16));
    await session.step(22, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(23, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(24, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(25, "And the context panel is 260 px wider", () => contextPanelWider(page, 260));
    await session.step(26, "And the grid shows only the columns \"name, biln\"", () => gridShowsColumns(page, "name, biln"));
    await session.step(27, "And the \"biln\" column of the grid is 420 px wide", () => gridColumnWidth(page, "biln", 420));
    await session.step(28, "And the monomers of the \"biln\" column are drawn at 16 px", () => sequenceFontSize(page, "biln", 16));
    await session.step(29, "And the monomer names of the \"biln\" column are cut at 8 characters", () => monomerNameLength(page, "biln", 8));
    await session.step(30, "And the rows of the grid are 64 px tall", () => gridRowHeight(page, 64));
    await session.step(31, "When user clicks on the \"cell 2 of biln\" area of grid", () => clickArea(page, "cell 2 of biln", el("grid")));
    await session.step(32, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(33, "When user expands \"Molecular Structure\" section", () => expand(page, el("\"Molecular Structure\" section")));
    await session.step(34, "Then the \"Molecular Structure\" section should show a drawn molecule", () => paneShowsDrawing(page, "Molecular Structure"));
    await session.step(35, "When the video says \"Octreotide, built atom by atom from its monomers\"", () => videoSays(page, "Octreotide, built atom by atom from its monomers"));
    await session.step(36, "And user clicks on the molecule in the \"Molecular Structure\" section", () => clickPaneMolecule(page, "Molecular Structure"));
    await session.step(37, "And user hovers over the molecule in the \"Molecule\" dialog at 36% across and 45% down", () => moveLens(page, "Molecule", 36, 45));
    await session.step(38, "And user hovers over the molecule in the \"Molecule\" dialog at 62% across and 22% down", () => moveLens(page, "Molecule", 62, 22));
    await session.step(39, "And user closes \"Molecule\" dialog", () => close(page, el("\"Molecule\" dialog")));
    await session.step(40, "And user collapses \"Molecular Structure\" section", () => collapse(page, el("\"Molecular Structure\" section")));
    await session.step(41, "And user expands \"Molecular 3D Structure\" section", () => expand(page, el("\"Molecular 3D Structure\" section")));
    await session.step(42, "Then the \"Molecular 3D Structure\" section should show a drawn molecule", () => paneShowsDrawing(page, "Molecular 3D Structure"));
    await session.step(43, "When the video says \"3D coordinates computed in the browser\"", () => videoSays(page, "3D coordinates computed in the browser"));
    await session.step(44, "And user turns the 3D structure in the \"Molecular 3D Structure\" section", () => turn3d(page, "Molecular 3D Structure"));
    await session.step(45, "And user collapses \"Molecular 3D Structure\" section", () => collapse(page, el("\"Molecular 3D Structure\" section")));
    await session.step(46, "And user expands \"Molecular Structure\" section", () => expand(page, el("\"Molecular Structure\" section")));
    await session.step(47, "Then the \"Molecular Structure\" section should show a drawn molecule", () => paneShowsDrawing(page, "Molecular Structure"));
    await session.step(48, "When the video says \"Every chemistry panel works on the peptide\"", () => videoSays(page, "Every chemistry panel works on the peptide"));
    await session.step(49, "And user expands the \"Explore\" pane inside the \"Molecular Structure\" section", () => expandInnerPane(page, "Explore", "Molecular Structure"));
    await session.step(50, "And user expands the \"Chemistry\" pane inside the \"Explore\" section", () => expandInnerPane(page, "Chemistry", "Explore"));
    await session.step(51, "And user expands the \"Properties\" pane inside the \"Chemistry\" section", () => expandInnerPane(page, "Properties", "Chemistry"));
    await session.step(52, "Then no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
