/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-legend-and-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {categoricalColorIs, colorCategorical, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickLegendItem, legendItemColor, legendItemsDiffer, legendLists, legendSide, loadLayout, noErrors, painted, readingIs, readingReads, repainted, reportsNoError, saveLayoutToServer, setProperty, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart legend, category colours and what survives a round-trip", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-legend-and-persistence.feature", import.meta.url);
  test("Line chart legend, category colours and what survives a round-trip", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(18, "And user colors \"Stereo Category\" column categorically:", () => colorCategorical(page, "Stereo Category", [["R_ONE","#FF0000"],["S_ABS","#00FF00"],["S_ACHIR","#0000FF"],["S_PART","#FFFF00"],["S_UNKN","#FF00FF"]]));
    await session.step(24, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","Chemical Space X"],["yColumnNames","Chemical Space Y"],["splitColumnNames","Stereo Category"]]));
    await session.step(28, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(29, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
    await session.step(30, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
    await session.step(31, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
    await session.step(32, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("The legend comes with the split and goes when the split does", async () => {
      await session.step(35, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(36, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(37, "Then the legend of line chart viewer should be hidden", () => shouldBe(page, el("the legend of line chart viewer"), "hidden"));
      await session.step(38, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(39, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(40, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Position moves the legend to the side it names", async () => {
      await session.step(44, "When user sets \"legendPosition\" property of line chart viewer to \"Left\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Left"));
      await session.step(45, "Then the legend of line chart viewer should be on the left", () => legendSide(page, el("line chart viewer"), "left"));
      await session.step(46, "When user sets \"legendPosition\" property of line chart viewer to \"Top\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Top"));
      await session.step(47, "Then the legend of line chart viewer should be on the top", () => legendSide(page, el("line chart viewer"), "top"));
      await session.step(48, "When user sets \"legendPosition\" property of line chart viewer to \"Bottom\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Bottom"));
      await session.step(49, "Then the legend of line chart viewer should be on the bottom", () => legendSide(page, el("line chart viewer"), "bottom"));
      await session.step(50, "When user sets \"legendPosition\" property of line chart viewer to \"Right\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Right"));
      await session.step(51, "Then the legend of line chart viewer should be on the right", () => legendSide(page, el("line chart viewer"), "right"));
      await session.step(52, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(53, "When user sets \"legendPosition\" property of line chart viewer to \"Auto\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Auto"));
      await session.step(54, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Visibility Never hides the legend and Auto brings it back", async () => {
      await session.step(57, "When user sets \"legendVisibility\" property of line chart viewer to \"Never\"", () => setProperty(page, "legendVisibility", el("line chart viewer"), "Never"));
      await session.step(58, "Then the legend of line chart viewer should be hidden", () => shouldBe(page, el("the legend of line chart viewer"), "hidden"));
      await session.step(59, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(60, "When user sets \"legendVisibility\" property of line chart viewer to \"Always\"", () => setProperty(page, "legendVisibility", el("line chart viewer"), "Always"));
      await session.step(61, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(62, "When user sets \"legendVisibility\" property of line chart viewer to \"Auto\"", () => setProperty(page, "legendVisibility", el("line chart viewer"), "Auto"));
      await session.step(63, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking a legend category filters the chart down to it and back", async () => {
      await session.step(67, "When user clicks on \"S_ABS\" item in the legend of line chart viewer", () => clickLegendItem(page, "S_ABS", el("line chart viewer")));
      await session.step(68, "Then the \"rows shown\" reading of line chart viewer should be 2", () => readingIs(page, "rows shown", el("line chart viewer"), 2));
      await session.step(69, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(70, "And 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(71, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
      await session.step(72, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(73, "When user clicks on \"S_ABS\" item in the legend of line chart viewer", () => clickLegendItem(page, "S_ABS", el("line chart viewer")));
      await session.step(74, "Then the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(75, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The legend items keep the colours the column was given", async () => {
      await session.step(79, "Then the \"R_ONE\" item in the legend of line chart viewer should be colored \"#FF0000\"", () => legendItemColor(page, "R_ONE", el("line chart viewer"), "#FF0000"));
      await session.step(80, "And the \"S_ABS\" item in the legend of line chart viewer should be colored \"#00FF00\"", () => legendItemColor(page, "S_ABS", el("line chart viewer"), "#00FF00"));
      await session.step(81, "And the \"R_ONE\" and \"S_ABS\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "R_ONE", "S_ABS", el("line chart viewer")));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The chart configuration and the category colour come back from a saved layout (GROK-17278)", async () => {
      await session.step(85, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(86, "And user colors \"Stereo Category\" column categorically:", () => colorCategorical(page, "Stereo Category", [["R_ONE","#00AAFF"]]));
      await session.step(88, "Then the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#00AAFF\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#00AAFF"));
      await session.step(89, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(90, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(91, "When user loads the saved layout", () => loadLayout(page));
      await session.step(92, "Then the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(93, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(94, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
      await session.step(95, "And the \"x column\" reading of line chart viewer should be \"Chemical Space X\"", () => readingReads(page, "x column", el("line chart viewer"), "Chemical Space X"));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The chart and its colours come back from a saved project (GROK-19825)", async () => {
      await session.step(99, "When user saves the current view as project \"zz-linechart-legend-colors\"", () => saveAsProject(page, "zz-linechart-legend-colors"));
      await session.step(100, "And user closes all views", () => closeAllViews(page));
      await session.step(101, "And user opens the \"zz-linechart-legend-colors\" project", () => openProject(page, "zz-linechart-legend-colors"));
      await session.step(102, "Then line chart viewer should be added to the open tableview", () => viewerAdded(page, "line chart"));
      await session.step(103, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(104, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(105, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
      await session.step(106, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
