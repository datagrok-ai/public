/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/apps/usage-analysis-tabs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openAppOnView} from '../../bindings/demo-oligo.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Apps > Usage Analysis: a tab picked while the app starts stays, with its own sub-views", () => {
  const session = feature(test, "features/apps/usage-analysis-tabs.feature", import.meta.url);
  test("Functions, picked before the toolbox is built, stays current and gets its Execution time sub-view", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "Given user opens the \"Usage Analysis\" app on its \"Overview\" view", () => openAppOnView(page, "Usage Analysis", "Overview"));
    await session.step(14, "When user clicks on Functions tab", () => clickOn(page, el("Functions tab")));
    await session.step(15, "Then the \"Functions\" view should be current", () => viewIsCurrent(page, "Functions"));
    await session.step(16, "And Apply button should be visible", () => shouldBe(page, el("Apply button"), "visible"));
    await session.step(17, "And \"Execution time\" button should be visible", () => shouldBe(page, el("\"Execution time\" button"), "visible"));
    await session.step(18, "And the \"Functions\" view should be current", () => viewIsCurrent(page, "Functions"));
    await session.step(19, "And no errors should have been logged", () => noErrors(page));
  });
});
