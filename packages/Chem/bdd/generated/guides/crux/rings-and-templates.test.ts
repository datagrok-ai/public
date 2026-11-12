/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/rings-and-templates.feature
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
import {cruxOpenEmpty, cruxSpot} from '../../../bindings/crux.js';
import {sketcherHolds} from '../../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Rings and templates in Crux", () => {
  const session = feature(test, "features/guides/crux/rings-and-templates.feature", import.meta.url);
  test("Place a ring, fuse another, add a spiro ring and a template from the Structure Library", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And the Crux sketcher is open", () => cruxOpenEmpty(page));
    await session.step(16, "When user clicks on Crux benzene tool", () => clickOn(page, el("Crux benzene tool")), undefined, "Pick benzene from the ring bar");
    await session.step(18, "And user clicks on Crux canvas", () => clickOn(page, el("Crux canvas")), undefined, "Click the canvas to place it");
    await session.step(20, "And user clicks on Crux cyclohexane tool", () => clickOn(page, el("Crux cyclohexane tool")), undefined, "Pick cyclohexane");
    await session.step(22, "And user clicks on the \"bond 0\" area of Crux sketcher widget", () => clickArea(page, "bond 0", el("Crux sketcher widget")), undefined, "Click a bond of the benzene: the new ring is fused onto it");
    await session.step(24, "Then the sketcher in sketcher dialog should hold the molecule \"c1ccc2c(c1)CCCC2\"", () => sketcherHolds(page, el("sketcher dialog"), "c1ccc2c(c1)CCCC2"), undefined, "Tetralin: a benzene fused with a cyclohexane");
    await session.step(26, "When user clicks on Crux cyclopropane tool", () => clickOn(page, el("Crux cyclopropane tool")), undefined, "Pick cyclopropane");
    await session.step(28, "And user clicks on the \"atom 7\" area of Crux sketcher widget", () => clickArea(page, "atom 7", el("Crux sketcher widget")), undefined, "Click a CH2 of the new ring: a spiro ring grows there");
    await session.step(30, "Then the sketcher in sketcher dialog should hold the molecule \"c1ccc2c(c1)CCC1(CC1)C2\"", () => sketcherHolds(page, el("sketcher dialog"), "c1ccc2c(c1)CCC1(CC1)C2"), undefined, "A spiro compound: the two rings share one carbon");
    await session.step(32, "When user clicks on Crux structure library button", () => clickOn(page, el("Crux structure library button")), undefined, "Open the Structure Library");
    await session.step(34, "And user types \"azulene\" into Crux structure library search", () => typeInto(page, "azulene", el("Crux structure library search")), undefined, "Search for azulene: its group opens with the match");
    await session.step(36, "And user clicks on \"Azulene\" button inside Crux structure library", () => clickOn(page, el("\"Azulene\" button inside Crux structure library")), undefined, "Pick azulene");
    await session.step(38, "And user clicks on Crux canvas 80% across and 25% down", () => cruxSpot(page, 80, 25), undefined, "Click an empty spot of the canvas to place it");
    await session.step(40, "Then the sketcher in sketcher dialog should hold the molecule \"c1ccc2c(c1)CCC1(CC1)C2.c1ccc2cccc-2cc1\"", () => sketcherHolds(page, el("sketcher dialog"), "c1ccc2c(c1)CCC1(CC1)C2.c1ccc2cccc-2cc1"), undefined, "Azulene sits beside the spiro compound");
  });
});
