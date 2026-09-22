/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/export-filtered-rows.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, downloadThrough, downloadedContains, downloadedNotContains} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, pickFromOpenMenu, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Export only the filtered rows", () => {
  const session = feature(test, "features/guides/export-filtered-rows.feature", import.meta.url);
  test("Extract the filtered rows into their own table, then export that table", {tag: ["@guide", "@help:transform"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(13, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
    await session.step(14, "And user clicks on the \"category RA of DIS_POP\" area of filter panel", () => clickArea(page, "category RA of DIS_POP", el("filter panel")));
    await session.step(15, "Then 2550 rows should pass the filter", () => filterPasses(page, 2550));
    await session.step(16, "When user clicks on \"Filtered: 2,550\" text in status bar", () => clickOn(page, el("\"Filtered: 2,550\" text in status bar")));
    await session.step(17, "And user picks \"Extract Rows\" from the open menu", () => pickFromOpenMenu(page, "Extract Rows"));
    await session.step(18, "Then grid should show 2550 rows", () => showsRows(page, el("grid"), 2550));
    await session.step(19, "When user clicks on Export icon in toolbar", () => clickOn(page, el("Export icon in toolbar")));
    await session.step(20, "And user downloads a file through \"As CSV\" text in toolbar", () => downloadThrough(page, el("\"As CSV\" text in toolbar")));
    await session.step(21, "Then the downloaded file should contain \"RA\"", () => downloadedContains(page, "RA"));
    await session.step(22, "And the downloaded file should not contain \"Psoriasis\"", () => downloadedNotContains(page, "Psoriasis"));
  });
});
