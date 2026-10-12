/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/substructure-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/crux.js';
import '../../../bindings/datasets.js';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter a table by a query drawn in Crux", () => {
  const session = feature(test, "features/guides/crux/substructure-filter.feature", import.meta.url);
  test("Filter by a benzene ring, then widen the query with an atom list", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")), undefined, "Open a table of 100 molecules");
    await session.step(17, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")), undefined, "Open the filter panel");
    await session.step(19, "And user clicks on \"Sketch\" text in \"Structure\" filter card", () => clickOn(page, el("\"Sketch\" text in \"Structure\" filter card")), undefined, "Click Sketch on the Structure filter: Crux opens in query mode");
    await session.step(21, "And user clicks on Crux benzene tool", () => clickOn(page, el("Crux benzene tool")), undefined, "Pick benzene");
    await session.step(23, "And user clicks on Crux canvas", () => clickOn(page, el("Crux canvas")), undefined, "Draw it: the table filters as you draw");
    await session.step(25, "Then 32 rows should pass the filter", () => filterPasses(page, 32), undefined, "32 molecules hold a benzene ring");
    await session.step(27, "When user clicks on Crux periodic table button", () => clickOn(page, el("Crux periodic table button")), undefined, "Open the periodic table");
    await session.step(29, "And user clicks on Crux periodic table list button", () => clickOn(page, el("Crux periodic table list button")), undefined, "Choose List");
    await session.step(31, "And user clicks on \"Carbon\" button inside Crux periodic table", () => clickOn(page, el("\"Carbon\" button inside Crux periodic table")), undefined, "Choose carbon …");
    await session.step(33, "And user clicks on \"Nitrogen\" button inside Crux periodic table", () => clickOn(page, el("\"Nitrogen\" button inside Crux periodic table")), undefined, "… and nitrogen");
    await session.step(35, "And user clicks on Crux periodic table Add button", () => clickOn(page, el("Crux periodic table Add button")), undefined, "Add the list");
    await session.step(37, "And user clicks on the \"atom 0\" area of Crux sketcher widget", () => clickArea(page, "atom 0", el("Crux sketcher widget")), undefined, "Click a ring atom: it may now be C or N");
    await session.step(39, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")), undefined, "OK");
    await session.step(41, "Then 44 rows should pass the filter", () => filterPasses(page, 44), undefined, "44 molecules hold a benzene or a pyridine ring");
  });
});
