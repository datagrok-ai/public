/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map, core#931, core#1202]
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
import {clickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, dragAreaBy, hasArea, hasNoArea, loadLayout, noErrors, painted, readingAsRemembered, readingBetween, readingIs, readingReads, rememberReading, saveLayoutToServer, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {stackViewer, switchTab, toggleInColumnList, viewersTabbed} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map layouts that used to load broken", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map-layout.feature", import.meta.url);
  test("A zoomed heat map stacked as a tab is not collapsed after the layout loads, on the tab it was saved on (public#2927, core#931)", {tag: ["@viewers", "@realizes:viewers.heat-map", "@realizes:core#931"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(27, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(28, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(29, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(30, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(34, "Given user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(35, "When user stacks heat map viewer onto histogram viewer as a tab", () => stackViewer(page, el("heat map viewer"), el("histogram viewer")));
    await session.step(36, "And user switches the tabbed panel to the \"Heat map\" tab", () => switchTab(page, "Heat map"));
    await session.step(37, "And user drags the \"x scroll max handle\" area of heat map viewer by 200 pixels to the left", () => dragAreaBy(page, "x scroll max handle", el("heat map viewer"), 200, "left"));
    await session.step(38, "Then the \"x scroll span\" reading of heat map viewer should be between 0.4 and 0.9", () => readingBetween(page, "x scroll span", el("heat map viewer"), 0.4, 0.9));
    await session.step(39, "And the \"columns shown\" reading of heat map viewer should be between 4 and 10", () => readingBetween(page, "columns shown", el("heat map viewer"), 4, 10));
    await session.step(40, "When user remembers the \"columns shown\" reading of heat map viewer", () => rememberReading(page, "columns shown", el("heat map viewer")));
    await session.step(41, "And user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(42, "And user closes all views", () => closeAllViews(page));
    await session.step(43, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(44, "And user loads the saved layout", () => loadLayout(page));
    await session.step(45, "Then the viewers \"Histogram, Heat map\" should be the tabs of one panel", () => viewersTabbed(page, "Histogram, Heat map"));
    await session.step(46, "And histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
    await session.step(47, "And the \"columns shown\" reading of heat map viewer should be as remembered", () => readingAsRemembered(page, "columns shown", el("heat map viewer")));
    await session.step(48, "And the \"x scroll span\" reading of heat map viewer should be between 0.4 and 0.9", () => readingBetween(page, "x scroll span", el("heat map viewer"), 0.4, 0.9));
    await session.step(49, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(50, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout saved after Is Heatmap off and on loads as the whole heat map (public#2901)", {tag: ["@viewers", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(27, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(28, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(29, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(30, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(54, "When user sets \"isHeatmap\" property of heat map viewer to \"false\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "false"));
    await session.step(55, "Then the \"row height\" reading of heat map viewer should be between 20 and 40", () => readingBetween(page, "row height", el("heat map viewer"), 20, 40));
    await session.step(56, "When user sets \"isHeatmap\" property of heat map viewer to \"true\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "true"));
    await session.step(57, "And user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(58, "And user closes all views", () => closeAllViews(page));
    await session.step(59, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(60, "And user loads the saved layout", () => loadLayout(page));
    await session.step(61, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(62, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(63, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
    await session.step(64, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
    await session.step(65, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout saved in grid mode with Max Heatmap Columns 200 turns back into the heat map (public#2901)", {tag: ["@viewers", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(27, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(28, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(29, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(30, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(69, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"200\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "200"));
    await session.step(70, "And user sets \"isHeatmap\" property of heat map viewer to \"false\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "false"));
    await session.step(71, "And user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(72, "And user closes all views", () => closeAllViews(page));
    await session.step(73, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(74, "And user loads the saved layout", () => loadLayout(page));
    await session.step(75, "Then the \"is heatmap\" reading of heat map viewer should be \"false\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "false"));
    await session.step(76, "And the \"row height\" reading of heat map viewer should be between 20 and 40", () => readingBetween(page, "row height", el("heat map viewer"), 20, 40));
    await session.step(77, "And heat map viewer should have a \"cell 1 of AGE\" area", () => hasArea(page, el("heat map viewer"), "cell 1 of AGE"));
    await session.step(78, "When user sets \"isHeatmap\" property of heat map viewer to \"true\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "true"));
    await session.step(79, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(80, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(81, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
    await session.step(82, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
    await session.step(83, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(84, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout saved on 1000 rows loads on a 100-row table with all 100 rows shown (GROK-16239)", {tag: ["@viewers", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(27, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(28, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(29, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(30, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(87, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(88, "And user closes all views", () => closeAllViews(page));
    await session.step(89, "And user opens demog-1000 dataset keeping the first 100 rows", () => openDatasetRows(page, ds("demog-1000"), 100));
    await session.step(90, "And user loads the saved layout", () => loadLayout(page));
    await session.step(91, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(92, "And the \"rows\" reading of heat map viewer should be 100", () => readingIs(page, "rows", el("heat map viewer"), 100));
    await session.step(93, "And the \"rows shown\" reading of heat map viewer should be 100", () => readingIs(page, "rows shown", el("heat map viewer"), 100));
    await session.step(94, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(95, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
    await session.step(96, "And heat map viewer should have a \"cell 1 of AGE\" area", () => hasArea(page, el("heat map viewer"), "cell 1 of AGE"));
    await session.step(97, "And heat map viewer should have a \"cell 100 of AGE\" area", () => hasArea(page, el("heat map viewer"), "cell 100 of AGE"));
    await session.step(98, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(99, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column hidden under Max Heatmap Columns stays hidden when the cap is raised (core#1202)", {tag: ["@viewers", "@realizes:viewers.heat-map", "@known-failure:115", "@realizes:core#1202"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(27, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(28, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(29, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(30, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await knownFailure(async () => {
      await session.step(103, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"5\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "5"));
      await session.step(104, "Then the \"columns shown\" reading of heat map viewer should be 5", () => readingIs(page, "columns shown", el("heat map viewer"), 5));
      await session.step(105, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
      await session.step(106, "When user clicks on settings icon of heat map viewer", () => clickOn(page, el("settings icon of heat map viewer")));
      await session.step(107, "Given \"Columns\" category in context panel is expanded", () => isExpanded(page, el("\"Columns\" category in context panel")));
      await session.step(108, "When user clicks on first button in \"Columns\" property in context panel", () => clickOn(page, el("first button in \"Columns\" property in context panel")));
      await session.step(109, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
      await session.step(110, "When user toggles the \"AGE\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "AGE", el("Order or Hide Columns dialog")));
      await session.step(111, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
      await session.step(112, "Then the \"columns shown\" reading of heat map viewer should be 4", () => readingIs(page, "columns shown", el("heat map viewer"), 4));
      await session.step(113, "And heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
      await session.step(114, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"6\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "6"));
      await session.step(115, "Then heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
    }, 115);
  });
});
