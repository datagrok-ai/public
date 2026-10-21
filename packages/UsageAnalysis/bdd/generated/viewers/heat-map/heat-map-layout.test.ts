/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map]
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
import {closeAllViews, openDataset, openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, dragAreaBy, hasArea, loadLayout, noErrors, painted, readingAsRemembered, readingBetween, readingIs, readingReads, rememberReading, reportsNoError, saveLayoutToServer, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {stackViewer, switchTab, viewersTabbed} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map layouts that used to load broken", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map-layout.feature", import.meta.url);
  test("A zoomed heat map stacked as a tab is not collapsed after the layout loads (public#2927)", {tag: ["@viewers", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(21, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(22, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(23, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(26, "Given user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(27, "When user stacks heat map viewer onto histogram viewer as a tab", () => stackViewer(page, el("heat map viewer"), el("histogram viewer")));
    await session.step(28, "And user switches the tabbed panel to the \"Heat map\" tab", () => switchTab(page, "Heat map"));
    await session.step(29, "And user drags the \"x scroll max handle\" area of heat map viewer by 200 pixels to the left", () => dragAreaBy(page, "x scroll max handle", el("heat map viewer"), 200, "left"));
    await session.step(30, "Then the \"x scroll span\" reading of heat map viewer should be between 0.4 and 0.9", () => readingBetween(page, "x scroll span", el("heat map viewer"), 0.4, 0.9));
    await session.step(31, "And the \"columns shown\" reading of heat map viewer should be between 4 and 10", () => readingBetween(page, "columns shown", el("heat map viewer"), 4, 10));
    await session.step(32, "When user remembers the \"columns shown\" reading of heat map viewer", () => rememberReading(page, "columns shown", el("heat map viewer")));
    await session.step(33, "And user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(34, "And user closes all views", () => closeAllViews(page));
    await session.step(35, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(36, "And user loads the saved layout", () => loadLayout(page));
    await session.step(37, "Then the viewers \"Histogram, Heat map\" should be the tabs of one panel", () => viewersTabbed(page, "Histogram, Heat map"));
    await session.step(38, "When user switches the tabbed panel to the \"Heat map\" tab", () => switchTab(page, "Heat map"));
    await session.step(39, "Then the \"columns shown\" reading of heat map viewer should be as remembered", () => readingAsRemembered(page, "columns shown", el("heat map viewer")));
    await session.step(40, "And the \"x scroll span\" reading of heat map viewer should be between 0.4 and 0.9", () => readingBetween(page, "x scroll span", el("heat map viewer"), 0.4, 0.9));
    await session.step(41, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(42, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(43, "And heat map viewer should report no error", () => reportsNoError(page, el("heat map viewer")));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout saved after Is Heatmap off and on loads as the whole heat map (public#2901)", {tag: ["@viewers", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(21, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(22, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(23, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(47, "When user sets \"isHeatmap\" property of heat map viewer to \"false\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "false"));
    await session.step(48, "Then the \"row height\" reading of heat map viewer should be between 20 and 40", () => readingBetween(page, "row height", el("heat map viewer"), 20, 40));
    await session.step(49, "When user sets \"isHeatmap\" property of heat map viewer to \"true\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "true"));
    await session.step(50, "And user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(51, "And user closes all views", () => closeAllViews(page));
    await session.step(52, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(53, "And user loads the saved layout", () => loadLayout(page));
    await session.step(54, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(55, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(56, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
    await session.step(57, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
    await session.step(58, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(59, "And heat map viewer should report no error", () => reportsNoError(page, el("heat map viewer")));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout saved in grid mode with Max Heatmap Columns 200 turns back into the heat map (public#2901)", {tag: ["@viewers", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(21, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(22, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(23, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(63, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"200\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "200"));
    await session.step(64, "And user sets \"isHeatmap\" property of heat map viewer to \"false\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "false"));
    await session.step(65, "And user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(66, "And user closes all views", () => closeAllViews(page));
    await session.step(67, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(68, "And user loads the saved layout", () => loadLayout(page));
    await session.step(69, "Then the \"is heatmap\" reading of heat map viewer should be \"false\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "false"));
    await session.step(70, "And the \"row height\" reading of heat map viewer should be between 20 and 40", () => readingBetween(page, "row height", el("heat map viewer"), 20, 40));
    await session.step(71, "And heat map viewer should have a \"cell 1 of AGE\" area", () => hasArea(page, el("heat map viewer"), "cell 1 of AGE"));
    await session.step(72, "When user sets \"isHeatmap\" property of heat map viewer to \"true\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "true"));
    await session.step(73, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(74, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(75, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
    await session.step(76, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
    await session.step(77, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(78, "And heat map viewer should report no error", () => reportsNoError(page, el("heat map viewer")));
    await session.step(79, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout saved on 1000 rows loads on a 100-row table with all 100 rows shown (GROK-16239)", {tag: ["@viewers", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(21, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(22, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(23, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(82, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(83, "And user closes all views", () => closeAllViews(page));
    await session.step(84, "And user opens demog-1000 dataset keeping the first 100 rows", () => openDatasetRows(page, ds("demog-1000"), 100));
    await session.step(85, "And user loads the saved layout", () => loadLayout(page));
    await session.step(86, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(87, "And the \"rows\" reading of heat map viewer should be 100", () => readingIs(page, "rows", el("heat map viewer"), 100));
    await session.step(88, "And the \"rows shown\" reading of heat map viewer should be 100", () => readingIs(page, "rows shown", el("heat map viewer"), 100));
    await session.step(89, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(90, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
    await session.step(91, "And heat map viewer should have a \"cell 1 of AGE\" area", () => hasArea(page, el("heat map viewer"), "cell 1 of AGE"));
    await session.step(92, "And heat map viewer should have a \"cell 100 of AGE\" area", () => hasArea(page, el("heat map viewer"), "cell 100 of AGE"));
    await session.step(93, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(94, "And heat map viewer should report no error", () => reportsNoError(page, el("heat map viewer")));
    await session.step(95, "And no errors should have been logged", () => noErrors(page));
  });
});
