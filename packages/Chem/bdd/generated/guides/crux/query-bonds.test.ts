/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/query-bonds.feature
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
import {cruxMenu, cruxOpenOn, cruxSmarts} from '../../../bindings/crux.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Query bonds in Crux", () => {
  const session = feature(test, "features/guides/crux/query-bonds.feature", import.meta.url);
  test("Make one bond any bond and another a chain bond", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And the Crux sketcher is open on \"NCCc1ccccc1\"", () => cruxOpenOn(page, "NCCc1ccccc1"));
    await session.step(15, "When user clicks on Crux settings button", () => clickOn(page, el("Crux settings button")), undefined, "Open the gear menu");
    await session.step(17, "And user clicks on Crux query mode item", () => clickOn(page, el("Crux query mode item")), undefined, "Switch on Query mode");
    await session.step(19, "And user clicks on Crux query bond tool", () => clickOn(page, el("Crux query bond tool")), undefined, "Pick the query bond tool: Any bond");
    await session.step(21, "And user clicks on the \"bond 2\" area of Crux sketcher widget", () => clickArea(page, "bond 2", el("Crux sketcher widget")), undefined, "Click the bond to the ring: it may now be any bond");
    await session.step(23, "Then the Crux sketcher should hold the query \"[#7]-[#6]-[#6]~[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1\"", () => cruxSmarts(page, "[#7]-[#6]-[#6]~[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"), undefined, "The query takes any bond there (~)");
    await session.step(25, "When user opens the Crux context menu on the \"bond 1\" area", () => cruxMenu(page, "bond 1"), undefined, "Right-click the C–C bond of the side chain");
    await session.step(27, "And user clicks on Crux topology item", () => clickOn(page, el("Crux topology item")), undefined, "Open Topology");
    await session.step(29, "And user clicks on Crux chain topology item", () => clickOn(page, el("Crux chain topology item")), undefined, "Choose Chain: the bond must not be in a ring");
    await session.step(31, "Then the Crux sketcher should hold the query \"[#7]-[#6]-&!@[#6]~[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1\"", () => cruxSmarts(page, "[#7]-[#6]-&!@[#6]~[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"), undefined, "The query asks for a chain bond there (!@)");
  });
});
