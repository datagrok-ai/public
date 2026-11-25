/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/reactions-edited-crux.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '../../bindings/reaction-export.js';
import '../../bindings/reactions-app-text.js';
import '../../bindings/sketcher-app.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {engineGate} from '../../bindings/engine.js';
import {appSides, openReactionsApp} from '../../bindings/reaction-editor.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelClosed, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Chemistry · Reactions edited in Crux Sketch", () => {
  const session = feature(test, "features/demos/reactions-edited-crux.feature", import.meta.url);
  test("The Reactions app draws what you type, and writes what you draw", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(17, "And user opens Chem's Reactions app", () => openReactionsApp(page));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "When the video says \"Type a reaction: Crux Sketch draws it\"", () => videoSays(page, "Type a reaction: Crux Sketch draws it"));
    await session.step(21, "And user enters \"Benzylamine acetylated\" into \"Title\" input", () => enterInto(page, "Benzylamine acetylated", el("\"Title\" input")));
    await session.step(22, "And user types \"CC(=O)Cl.NCc1ccccc1>>CC(=O)NCc1ccccc1\" into reactions app text", () => typeInto(page, "CC(=O)Cl.NCc1ccccc1>>CC(=O)NCc1ccccc1", el("reactions app text")));
    await session.step(23, "Then the \"atoms\" reading of crux sketcher widget should be 23", () => readingIs(page, "atoms", el("crux sketcher widget"), 23));
    await session.step(24, "When the video says \"Draw on it: the text and the picture follow\"", () => videoSays(page, "Draw on it: the text and the picture follow"));
    await session.step(25, "And user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
    await session.step(27, "And user clicks on the \"atom 4\" area of crux sketcher widget", () => clickArea(page, "atom 4", el("crux sketcher widget")));
    await session.step(28, "And user clicks on the \"atom 15\" area of crux sketcher widget", () => clickArea(page, "atom 15", el("crux sketcher widget")));
    await session.step(29, "Then the Reactions app's text should have the reactants \"CC(=O)Cl.CNCc1ccccc1\" and the products \"CC(=O)N(C)Cc1ccccc1\"", () => appSides(page, "CC(=O)Cl.CNCc1ccccc1", "CC(=O)N(C)Cc1ccccc1"));
    await session.step(30, "When user hovers over reactions app drawing", () => hoverOver(page, el("reactions app drawing")));
    await session.step(31, "Then no errors should have been logged", () => noErrors(page));
  });
});
