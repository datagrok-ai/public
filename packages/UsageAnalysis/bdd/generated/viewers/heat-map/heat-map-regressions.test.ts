/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map, core#928, core#931]
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
import {addViewer, hasArea, hasNoArea, loadLayout, noErrors, painted, readingIs, saveLayoutToServer, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {notShowsText, showsText, stackViewer, switchTab, viewersTabbed} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map regressions of closed issues", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map-regressions.feature", import.meta.url);
  test("Max Heatmap Columns 0 says why the map is empty, and 100 takes the hint away (core#928)", {tag: ["@viewers", "@realizes:viewers.heat-map", "@realizes:core#928"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(15, "Then the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(16, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(20, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"0\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "0"));
    await session.step(21, "Then heat map viewer should show the text \"No columns to show: Max Heatmap Columns is 0\"", () => showsText(page, el("heat map viewer"), "No columns to show: Max Heatmap Columns is 0"));
    await session.step(22, "And the \"columns shown\" reading of heat map viewer should be 0", () => readingIs(page, "columns shown", el("heat map viewer"), 0));
    await session.step(23, "And heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
    await session.step(24, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"100\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "100"));
    await session.step(25, "Then the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(26, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(27, "And heat map viewer should not show the text \"No columns to show\"", () => notShowsText(page, el("heat map viewer"), "No columns to show"));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout brings the stacked heat map back on its own tab (core#931)", {tag: ["@viewers", "@realizes:viewers.heat-map", "@realizes:core#931"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(15, "Then the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(16, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(32, "Given user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(33, "When user stacks heat map viewer onto histogram viewer as a tab", () => stackViewer(page, el("heat map viewer"), el("histogram viewer")));
    await session.step(34, "And user switches the tabbed panel to the \"Histogram\" tab", () => switchTab(page, "Histogram"));
    await session.step(35, "And user switches the tabbed panel to the \"Heat map\" tab", () => switchTab(page, "Heat map"));
    await session.step(36, "Then heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(37, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(38, "And user closes all views", () => closeAllViews(page));
    await session.step(39, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(40, "And user loads the saved layout", () => loadLayout(page));
    await session.step(41, "Then the viewers \"Histogram, Heat map\" should be the tabs of one panel", () => viewersTabbed(page, "Histogram, Heat map"));
    await session.step(42, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(43, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(44, "And histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
});
