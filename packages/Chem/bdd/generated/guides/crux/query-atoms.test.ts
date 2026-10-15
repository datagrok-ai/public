/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/query-atoms.feature
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
import {cruxOpenOn, cruxSmarts} from '../../../bindings/crux.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Query atoms in Crux", () => {
  const session = feature(test, "features/guides/crux/query-atoms.feature", import.meta.url);
  test("Switch to query mode, make a ring atom \"N or O\" and the acid's OH any heteroatom", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And the Crux sketcher is open on \"OC(=O)c1ccccc1\"", () => cruxOpenOn(page, "OC(=O)c1ccccc1"));
    await session.step(16, "When user clicks on Crux settings button", () => clickOn(page, el("Crux settings button")), undefined, "Open the gear menu");
    await session.step(18, "And user clicks on Crux query mode item", () => clickOn(page, el("Crux query mode item")), undefined, "Switch on Query mode");
    await session.step(20, "And user clicks on Crux periodic table button", () => clickOn(page, el("Crux periodic table button")), undefined, "Open the periodic table");
    await session.step(22, "And user clicks on Crux periodic table list button", () => clickOn(page, el("Crux periodic table list button")), undefined, "Choose List");
    await session.step(24, "And user clicks on \"Nitrogen\" button inside Crux periodic table", () => clickOn(page, el("\"Nitrogen\" button inside Crux periodic table")), undefined, "Choose nitrogen");
    await session.step(26, "And user clicks on \"Oxygen\" button inside Crux periodic table", () => clickOn(page, el("\"Oxygen\" button inside Crux periodic table")), undefined, "… and oxygen");
    await session.step(28, "And user clicks on Crux periodic table Add button", () => clickOn(page, el("Crux periodic table Add button")), undefined, "Add the list");
    await session.step(30, "And user clicks on the \"atom 4\" area of Crux sketcher widget", () => clickArea(page, "atom 4", el("Crux sketcher widget")), undefined, "Click a ring atom: it may now be N or O");
    await session.step(32, "Then the Crux sketcher should hold the query \"[#8]-[#6](=[#8])-[#6]1:[#7,#8]:[#6]:[#6]:[#6]:[#6]:1\"", () => cruxSmarts(page, "[#8]-[#6](=[#8])-[#6]1:[#7,#8]:[#6]:[#6]:[#6]:[#6]:1"), undefined, "The query asks for N or O at that ring position");
    await session.step(34, "When user clicks on Crux periodic table button", () => clickOn(page, el("Crux periodic table button")), undefined, "Open the periodic table again");
    await session.step(36, "And user clicks on Crux periodic table Q button", () => clickOn(page, el("Crux periodic table Q button")), undefined, "Choose Q: any atom but carbon and hydrogen");
    await session.step(38, "And user clicks on Crux periodic table Add button", () => clickOn(page, el("Crux periodic table Add button")), undefined, "Add it");
    await session.step(40, "And user clicks on the \"atom 0\" area of Crux sketcher widget", () => clickArea(page, "atom 0", el("Crux sketcher widget")), undefined, "Click the acid's OH oxygen: it becomes Q");
    await session.step(42, "Then the Crux sketcher should hold the query \"[!#6&!#1]-[#6](=[#8])-[#6]1:[#7,#8]:[#6]:[#6]:[#6]:[#6]:1\"", () => cruxSmarts(page, "[!#6&!#1]-[#6](=[#8])-[#6]1:[#7,#8]:[#6]:[#6]:[#6]:[#6]:1"), undefined, "The query now takes any heteroatom there");
  });
});
