/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/apps/flow-crux.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, openApp, packageInstalled, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Flow's Sketcher Input node, drawn in Crux", () => {
  const session = feature(test, "features/apps/flow-crux.feature", import.meta.url);
  test("Each stroke in the node's Crux is one edit of the flow, and opening the node's sketcher again is none", {tag: ["@crux", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(15, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "Given user opens the \"Flow\" app", () => openApp(page, "Flow"));
    await session.step(21, "When user clicks on flow blank canvas card", () => clickOn(page, el("flow blank canvas card")));
    await session.step(22, "And user clicks on \"Inputs\" pane header in toolbox", () => clickOn(page, el("\"Inputs\" pane header in toolbox")));
    await session.step(23, "And user double-clicks on flow sketcher input item", () => doubleClickOn(page, el("flow sketcher input item")));
    await session.step(24, "Then the \"parameter edits\" reading of flow editor widget should be 0", () => readingIs(page, "parameter edits", el("flow editor widget"), 0));
    await session.step(25, "When user clicks on flow sketcher preview", () => clickOn(page, el("flow sketcher preview")));
    await session.step(26, "Then the \"ready\" reading of crux sketcher widget should be \"true\"", () => readingReads(page, "ready", el("crux sketcher widget"), "true"));
    await session.step(27, "And the \"parameter edits\" reading of flow editor widget should be 0", () => readingIs(page, "parameter edits", el("flow editor widget"), 0));
    await session.step(28, "When user clicks on crux benzene tool", () => clickOn(page, el("crux benzene tool")));
    await session.step(29, "And user clicks on crux canvas", () => clickOn(page, el("crux canvas")));
    await session.step(30, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"c1ccccc1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "c1ccccc1"));
    await session.step(31, "And the \"parameter edits\" reading of flow editor widget should be 1", () => readingIs(page, "parameter edits", el("flow editor widget"), 1));
    await session.step(32, "When user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
    await session.step(33, "And user clicks on the \"atom 0\" area of crux sketcher widget", () => clickArea(page, "atom 0", el("crux sketcher widget")));
    await session.step(34, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"Cc1ccccc1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "Cc1ccccc1"));
    await session.step(35, "And the \"parameter edits\" reading of flow editor widget should be 2", () => readingIs(page, "parameter edits", el("flow editor widget"), 2));
    await session.step(36, "When user clicks on flow sketcher Done button", () => clickOn(page, el("flow sketcher Done button")));
    await session.step(37, "And user clicks on flow sketcher preview", () => clickOn(page, el("flow sketcher preview")));
    await session.step(38, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"Cc1ccccc1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "Cc1ccccc1"));
    await session.step(39, "When user clicks on crux nitrogen tool", () => clickOn(page, el("crux nitrogen tool")));
    await session.step(40, "And user clicks on the \"atom 1\" area of crux sketcher widget", () => clickArea(page, "atom 1", el("crux sketcher widget")));
    await session.step(41, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"Cc1ccccn1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "Cc1ccccn1"));
    await session.step(42, "And the \"parameter edits\" reading of flow editor widget should be 3", () => readingIs(page, "parameter edits", el("flow editor widget"), 3));
  });
});
