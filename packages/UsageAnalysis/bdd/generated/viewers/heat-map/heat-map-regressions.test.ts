/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map, core#931]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, loadLayout, noErrors, painted, readingIs, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {stackViewer, switchTab, viewersTabbed} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map regressions of closed issues", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map-regressions.feature", import.meta.url);
  test("A layout brings the stacked heat map back on its own tab (core#931)", {tag: ["@viewers", "@realizes:viewers.heat-map", "@realizes:core#931"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(13, "Then the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(14, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(18, "Given user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(19, "When user stacks heat map viewer onto histogram viewer as a tab", () => stackViewer(page, el("heat map viewer"), el("histogram viewer")));
    await session.step(20, "And user switches the tabbed panel to the \"Histogram\" tab", () => switchTab(page, "Histogram"));
    await session.step(21, "And user switches the tabbed panel to the \"Heat map\" tab", () => switchTab(page, "Heat map"));
    await session.step(22, "Then heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(23, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(24, "And user closes all views", () => closeAllViews(page));
    await session.step(25, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "And user loads the saved layout", () => loadLayout(page));
    await session.step(27, "Then the viewers \"Histogram, Heat map\" should be the tabs of one panel", () => viewersTabbed(page, "Histogram, Heat map"));
    await session.step(28, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(29, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(30, "And histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
});
