/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/select-move-rotate.feature
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
import {cruxOpenOn, newCoordinates, rememberCoordinates} from '../../../bindings/crux.js';
import {sketcherHolds} from '../../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragBy} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea, dragAreaBy} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Select, move, turn, flip and clean up in Crux", () => {
  const session = feature(test, "features/guides/crux/select-move-rotate.feature", import.meta.url);
  test("Select the naphthalene, move it, turn it, flip it and clean up", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And the Crux sketcher is open on \"c1ccc2ccccc2c1.OCC\"", () => cruxOpenOn(page, "c1ccc2ccccc2c1.OCC"));
    await session.step(15, "And user remembers the coordinates of the Crux sketcher's molblock", () => rememberCoordinates(page));
    await session.step(17, "When user clicks on Crux select tool", () => clickOn(page, el("Crux select tool")), undefined, "Pick the selection tool");
    await session.step(19, "And user double-clicks on the \"atom 0\" area of Crux sketcher widget", () => doubleClickArea(page, "atom 0", el("Crux sketcher widget")), undefined, "Double-click the naphthalene: all of it is selected");
    await session.step(21, "And user drags the \"atom 0\" area of Crux sketcher widget by 120 pixels to the left", () => dragAreaBy(page, "atom 0", el("Crux sketcher widget"), 120, "left"), undefined, "Drag it to the left");
    await session.step(23, "And user drags Crux rotate handle by 90 and 0 pixels", () => dragBy(page, el("Crux rotate handle"), 90, 0), undefined, "Drag the round handle to turn it");
    await session.step(25, "And user clicks on Crux flip horizontal button", () => clickOn(page, el("Crux flip horizontal button")), undefined, "Flip it");
    await session.step(27, "And user clicks on Crux clean up button", () => clickOn(page, el("Crux clean up button")), undefined, "Clean Up tidies the drawing");
    await session.step(29, "Then the sketcher in sketcher dialog should hold the molecule \"c1ccc2ccccc2c1.CCO\"", () => sketcherHolds(page, el("sketcher dialog"), "c1ccc2ccccc2c1.CCO"), undefined, "Only the drawing moved: still naphthalene and ethanol");
    await session.step(30, "And the Crux sketcher's molblock should have new coordinates", () => newCoordinates(page));
  });
});
