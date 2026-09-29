/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/open-demo-dataset.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/nx.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Open a demo dataset", () => {
  const session = feature(test, "features/guides/open-demo-dataset.feature", import.meta.url);
  test("Open demog.csv from the Demo files", {tag: ["@guide", "@help:access/files"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And simple mode is off", () => simpleModeOff(page));
    await session.step(10, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(11, "When user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(12, "And user expands Files---Demo tree node inside browse tree", () => expand(page, el("Files---Demo tree node inside browse tree")));
    await session.step(13, "And user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(14, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(15, "And grid should show 5850 rows", () => showsRows(page, el("grid"), 5850));
  });
});
