/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/drawing.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/datasets.js';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {cruxKey, cruxOpenEmpty} from '../../../bindings/crux.js';
import {sketcherHolds} from '../../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, dragAreaBy} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Draw a molecule in Crux with bonds, a chain and hotkeys", () => {
  const session = feature(test, "features/guides/crux/drawing.feature", import.meta.url);
  test("Draw bonds and a chain, then change atoms and bonds with keys", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the Crux sketcher is open", () => cruxOpenEmpty(page));
    await session.step(17, "When user clicks on Crux single bond tool", () => clickOn(page, el("Crux single bond tool")), undefined, "Pick the bond tool");
    await session.step(19, "And user clicks on Crux canvas", () => clickOn(page, el("Crux canvas")), undefined, "Click the empty canvas: a bond appears");
    await session.step(21, "And user clicks on the \"atom 1\" area of Crux sketcher widget", () => clickArea(page, "atom 1", el("Crux sketcher widget")), undefined, "Click an end of the bond to add another");
    await session.step(23, "And user clicks on Crux chain tool", () => clickOn(page, el("Crux chain tool")), undefined, "Pick the chain tool");
    await session.step(25, "And user drags the \"atom 2\" area of Crux sketcher widget by 175 pixels to the right", () => dragAreaBy(page, "atom 2", el("Crux sketcher widget"), 175, "right"), undefined, "Drag from the end of the chain to draw more carbons");
    await session.step(27, "Then the sketcher in sketcher dialog should hold the molecule \"CCCCCCCC\"", () => sketcherHolds(page, el("sketcher dialog"), "CCCCCCCC"), undefined, "An octane chain");
    await session.step(29, "When user presses the \"o\" key over the \"atom 7\" area of Crux sketcher widget", () => cruxKey(page, "o", "atom 7", el("Crux sketcher widget")), undefined, "Point at the last carbon and press o: it becomes an oxygen");
    await session.step(31, "And user presses the \"2\" key over the \"bond 0\" area of Crux sketcher widget", () => cruxKey(page, "2", "bond 0", el("Crux sketcher widget")), undefined, "Point at the first bond and press 2: it becomes a double bond");
    await session.step(33, "And user presses the \"3\" key over the \"atom 3\" area of Crux sketcher widget", () => cruxKey(page, "3", "atom 3", el("Crux sketcher widget")), undefined, "Point at a carbon and press 3: a phenyl ring sprouts from it");
    await session.step(35, "Then the sketcher in sketcher dialog should hold the molecule \"C=CCC(CCCO)c1ccccc1\"", () => sketcherHolds(page, el("sketcher dialog"), "C=CCC(CCCO)c1ccccc1"), undefined, "4-Phenylhept-6-en-1-ol, drawn with three keys");
  });
});
