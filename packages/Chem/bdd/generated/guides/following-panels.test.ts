/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/following-panels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {engineGate} from '../../bindings/engine.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {expand, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, contextPanelOpen, contextPanelShowsCurrentCell, dragPaneOut, openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noLoader} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dockBesideViewer, dockToViewEdge, dragRightBorder, dragTopBorder} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Keep a molecule's panels in view while clicking through molecules", () => {
  const session = feature(test, "features/guides/following-panels.feature", import.meta.url);
  test("Dock Properties, Gasteiger and 3D under the grid", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(18, "And user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
    await session.step(19, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(21, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")), undefined, "Click a molecule");
    await session.step(23, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page), undefined, "Its panels are in the Context Panel");
    await session.step(25, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")), undefined, "Open the Chemistry group");
    await session.step(27, "And user drags Properties accordion header in context panel out of the accordion", () => dragPaneOut(page, el("Properties accordion header in context panel")), undefined, "Drag Properties out of the Context Panel: it becomes a viewer");
    await session.step(29, "Then \"Properties\" floating panel should contain the text \"C24H23NO3\"", () => shouldContainText(page, el("\"Properties\" floating panel"), "C24H23NO3"), undefined, "It shows the properties of this molecule");
    await session.step(31, "When user docks info panel viewer in \"Properties\" floating panel to the bottom edge of the view", () => dockToViewEdge(page, el("info panel viewer in \"Properties\" floating panel"), "bottom"), undefined, "Dock it under the grid");
    await session.step(33, "And user drags \"Gasteiger Partial Charges\" accordion header in context panel out of the accordion", () => dragPaneOut(page, el("\"Gasteiger Partial Charges\" accordion header in context panel")), undefined, "Drag out Gasteiger Partial Charges");
    await session.step(35, "And user docks info panel viewer in \"Gasteiger Partial Charges\" floating panel to the right side of info panel viewer in \"Properties\" dock panel", () => dockBesideViewer(page, el("info panel viewer in \"Gasteiger Partial Charges\" floating panel"), "right", el("info panel viewer in \"Properties\" dock panel")), undefined, "Dock it next to Properties");
    await session.step(37, "And user expands Structure accordion header in context panel", () => expand(page, el("Structure accordion header in context panel")), undefined, "Open the Structure group");
    await session.step(39, "And user drags \"3D Structure\" accordion header in context panel out of the accordion", () => dragPaneOut(page, el("\"3D Structure\" accordion header in context panel")), undefined, "Drag out 3D Structure");
    await session.step(41, "And user docks info panel viewer in \"3D Structure\" floating panel to the right side of info panel viewer in \"Gasteiger Partial Charges\" dock panel", () => dockBesideViewer(page, el("info panel viewer in \"3D Structure\" floating panel"), "right", el("info panel viewer in \"Gasteiger Partial Charges\" dock panel")), undefined, "Dock it next to Gasteiger Partial Charges");
    await session.step(43, "And user drags the top border of info panel viewer in \"Properties\" dock panel by 210 pixels up", () => dragTopBorder(page, el("info panel viewer in \"Properties\" dock panel"), 210), undefined, "Make the docked panels taller");
    await session.step(45, "And user drags the right border of info panel viewer in \"Properties\" dock panel by 132 pixels left", () => dragRightBorder(page, el("info panel viewer in \"Properties\" dock panel"), 132), undefined, "Share the width between them");
    await session.step(47, "And user drags the right border of info panel viewer in \"Gasteiger Partial Charges\" dock panel by 200 pixels left", () => dragRightBorder(page, el("info panel viewer in \"Gasteiger Partial Charges\" dock panel"), 200), undefined, "And with 3D Structure");
    await session.step(49, "And user clicks on the \"cell 3 of canonical_smiles\" area of grid", () => clickArea(page, "cell 3 of canonical_smiles", el("grid")), undefined, "Click another molecule");
    await session.step(51, "Then \"Properties\" dock panel should contain the text \"C22H23NO3\"", () => shouldContainText(page, el("\"Properties\" dock panel"), "C22H23NO3"), undefined, "Properties recalculates for it");
    await session.step(53, "And no loading indicator should be visible", () => noLoader(page), undefined, "So do the charges and the 3D structure");
    await session.step(55, "When user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")), undefined, "Click another one");
    await session.step(57, "Then \"Properties\" dock panel should contain the text \"C12H14F3NO\"", () => shouldContainText(page, el("\"Properties\" dock panel"), "C12H14F3NO"), undefined, "Its formula in Properties");
    await session.step(59, "And no loading indicator should be visible", () => noLoader(page), undefined, "The panels stay with the molecule you click");
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
});
