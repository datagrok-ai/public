/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/link-tables-selection-to-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectedRowCount, tableFilterCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dragAreaToArea, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter one table by the rows selected in another", () => {
  const session = feature(test, "features/guides/link-tables-selection-to-filter.feature", import.meta.url);
  test("Link two tables so that selecting rows in one filters the other", {tag: ["@guide", "@help:transform"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await session.step(13, "And user opens spgi-linked1 dataset", () => openDataset(page, ds("spgi-linked1")));
    await session.step(14, "When user clicks on spgi-100 tab", () => clickOn(page, el("spgi-100 tab")));
    await session.step(15, "And user picks \"Data > Link Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Link Tables..."));
    await session.step(16, "And user clicks on \"New Link\" text in \"Link Tables\" dialog", () => clickOn(page, el("\"New Link\" text in \"Link Tables\" dialog")));
    await session.step(17, "And user selects \"selection to filter\" in Link Type input in \"Link Tables\" dialog", () => selectIn(page, "selection to filter", el("Link Type input in \"Link Tables\" dialog")));
    await session.step(18, "And user clicks on LINK button in \"Link Tables\" dialog", () => clickOn(page, el("LINK button in \"Link Tables\" dialog")));
    await session.step(19, "Then \"spgi-100 -> SPGI-linked1\" text in \"Link Tables\" dialog should be visible", () => shouldBe(page, el("\"spgi-100 -> SPGI-linked1\" text in \"Link Tables\" dialog"), "visible"));
    await session.step(20, "When user clicks on CLOSE button in \"Link Tables\" dialog", () => clickOn(page, el("CLOSE button in \"Link Tables\" dialog")));
    await session.step(21, "And user drags the \"row header 1\" area of grid to the \"row header 5\" area", () => dragAreaToArea(page, "row header 1", el("grid"), "row header 5"));
    await session.step(22, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
    await session.step(23, "And 9 rows of table \"SPGI-linked1\" should pass the filter", () => tableFilterCount(page, 9, "SPGI-linked1"));
    await session.step(24, "When user clicks on SPGI-linked1 tab", () => clickOn(page, el("SPGI-linked1 tab")));
    await session.step(25, "Then grid should show 9 rows", () => showsRows(page, el("grid"), 9));
  });
});
