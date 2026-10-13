/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-legend-and-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
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
import {clickOn, followingShouldBe, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {categoricalColorIs, colorCategorical, filterPasses, filterTo, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, clickLegendItem, hasArea, hasNoArea, legendItemColor, legendItemsDiffer, legendLists, legendSide, loadLayout, noErrors, painted, pickColorSwatch, propertyShouldBe, readingIs, readingReads, repainted, reportsNoError, resizeTo, restoreSize, saveLayoutToServer, setProperties, setProperty, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart legend, category colours and what survives a round-trip", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-legend-and-persistence.feature", import.meta.url);
  test("Line chart legend, category colours and what survives a round-trip", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 15, page);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(42, "And user colors \"Stereo Category\" column categorically:", () => colorCategorical(page, "Stereo Category", [["R_ONE","#FF0000"],["S_ABS","#00FF00"],["S_ACHIR","#0000FF"],["S_PART","#FFFF00"],["S_UNKN","#FF00FF"]]), [["R_ONE","#FF0000"],["S_ABS","#00FF00"],["S_ACHIR","#0000FF"],["S_PART","#FFFF00"],["S_UNKN","#FF00FF"]]);
    await session.step(48, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","Chemical Space X"],["yColumnNames","Chemical Space Y"],["splitColumnNames","Stereo Category"]]), [["xColumnName","Chemical Space X"],["yColumnNames","Chemical Space Y"],["splitColumnNames","Stereo Category"]]);
    await session.step(52, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(53, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
    await session.step(54, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
    await session.step(55, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
    await session.step(56, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("The legend comes with the split and goes when the split does", async () => {
      await session.step(59, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(60, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(61, "Then the legend of line chart viewer should be hidden", () => shouldBe(page, el("the legend of line chart viewer"), "hidden"));
      await session.step(62, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(63, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(64, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Position moves the legend to the side it names", async () => {
      await session.step(68, "When user sets \"legendPosition\" property of line chart viewer to \"Left\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Left"));
      await session.step(69, "Then the legend of line chart viewer should be on the left", () => legendSide(page, el("line chart viewer"), "left"));
      await session.step(70, "When user sets \"legendPosition\" property of line chart viewer to \"Top\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Top"));
      await session.step(71, "Then the legend of line chart viewer should be on the top", () => legendSide(page, el("line chart viewer"), "top"));
      await session.step(72, "When user sets \"legendPosition\" property of line chart viewer to \"Bottom\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Bottom"));
      await session.step(73, "Then the legend of line chart viewer should be on the bottom", () => legendSide(page, el("line chart viewer"), "bottom"));
      await session.step(74, "When user sets \"legendPosition\" property of line chart viewer to \"Right\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Right"));
      await session.step(75, "Then the legend of line chart viewer should be on the right", () => legendSide(page, el("line chart viewer"), "right"));
      await session.step(76, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(77, "When user sets \"legendPosition\" property of line chart viewer to \"Auto\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Auto"));
      await session.step(78, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Visibility Never hides the legend and Auto brings it back", async () => {
      await session.step(81, "When user sets \"legendVisibility\" property of line chart viewer to \"Never\"", () => setProperty(page, "legendVisibility", el("line chart viewer"), "Never"));
      await session.step(82, "Then the legend of line chart viewer should be hidden", () => shouldBe(page, el("the legend of line chart viewer"), "hidden"));
      await session.step(83, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(84, "When user sets \"legendVisibility\" property of line chart viewer to \"Always\"", () => setProperty(page, "legendVisibility", el("line chart viewer"), "Always"));
      await session.step(85, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(86, "When user sets \"legendVisibility\" property of line chart viewer to \"Auto\"", () => setProperty(page, "legendVisibility", el("line chart viewer"), "Auto"));
      await session.step(87, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking a legend category filters the chart down to it and back", async () => {
      await session.step(91, "When user clicks on \"S_ABS\" item in the legend of line chart viewer", () => clickLegendItem(page, "S_ABS", el("line chart viewer")));
      await session.step(92, "Then the \"rows shown\" reading of line chart viewer should be 2", () => readingIs(page, "rows shown", el("line chart viewer"), 2));
      await session.step(93, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(94, "And 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(95, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
      await session.step(96, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(97, "When user clicks on \"S_ABS\" item in the legend of line chart viewer", () => clickLegendItem(page, "S_ABS", el("line chart viewer")));
      await session.step(98, "Then the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(99, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The legend items keep the colours the column was given", async () => {
      await session.step(103, "Then the \"R_ONE\" item in the legend of line chart viewer should be colored \"#FF0000\"", () => legendItemColor(page, "R_ONE", el("line chart viewer"), "#FF0000"));
      await session.step(104, "And the \"S_ABS\" item in the legend of line chart viewer should be colored \"#00FF00\"", () => legendItemColor(page, "S_ABS", el("line chart viewer"), "#00FF00"));
      await session.step(105, "And the \"R_ONE\" and \"S_ABS\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "R_ONE", "S_ABS", el("line chart viewer")));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The chart configuration and the category colour come back from a saved layout (GROK-17278)", async () => {
      await session.step(109, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(110, "And user colors \"Stereo Category\" column categorically:", () => colorCategorical(page, "Stereo Category", [["R_ONE","#00AAFF"]]), [["R_ONE","#00AAFF"]]);
      await session.step(112, "Then the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#00AAFF\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#00AAFF"));
      await session.step(113, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(114, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(115, "When user loads the saved layout", () => loadLayout(page));
      await session.step(116, "Then the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(117, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(118, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
      await session.step(119, "And the \"x column\" reading of line chart viewer should be \"Chemical Space X\"", () => readingReads(page, "x column", el("line chart viewer"), "Chemical Space X"));
      await session.step(120, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The chart and its colours come back from a saved project (GROK-19825)", async () => {
      await session.step(123, "When user saves the current view as project \"zz-linechart-legend-colors\"", () => saveAsProject(page, "zz-linechart-legend-colors"));
      await session.step(124, "And user closes all views", () => closeAllViews(page));
      await session.step(125, "And user opens the \"zz-linechart-legend-colors\" project", () => openProject(page, "zz-linechart-legend-colors"));
      await session.step(126, "Then line chart viewer should be added to the open tableview", () => viewerAdded(page, "line chart"));
      await session.step(127, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(128, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(129, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
      await session.step(130, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(131, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The default palette gives neighbouring split categories different colors", async () => {
      await session.step(134, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(135, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","STARTED"],["yColumnNames","WEIGHT"],["splitColumnNames","DIS_POP"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["xColumnName","STARTED"],["yColumnNames","WEIGHT"],["splitColumnNames","DIS_POP"],["Legend Visibility","Always"],["Legend Position","Right"]]);
      await session.step(141, "Then the legend of line chart viewer should list 6 items", () => legendLists(page, el("line chart viewer"), 6));
      await session.step(142, "And the \"AS\" and \"Indigestion\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "AS", "Indigestion", el("line chart viewer")));
      await session.step(143, "And the \"Indigestion\" and \"PsA\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "Indigestion", "PsA", el("line chart viewer")));
      await session.step(144, "And the \"PsA\" and \"Psoriasis\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "PsA", "Psoriasis", el("line chart viewer")));
      await session.step(145, "And the \"Psoriasis\" and \"RA\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "Psoriasis", "RA", el("line chart viewer")));
      await session.step(146, "And the \"RA\" and \"UC\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "RA", "UC", el("line chart viewer")));
      await session.step(147, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Multi Axis gives every Y column its own block of legend items", async () => {
      await session.step(150, "When user sets \"yColumnNames\" property of line chart viewer to \"WEIGHT, HEIGHT\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "WEIGHT, HEIGHT"));
      await session.step(151, "Then the legend of line chart viewer should list 6 items", () => legendLists(page, el("line chart viewer"), 6));
      await session.step(152, "When user sets \"multiAxis\" property of line chart viewer to \"true\"", () => setProperty(page, "multiAxis", el("line chart viewer"), "true"));
      await session.step(153, "Then the legend of line chart viewer should list 12 items", () => legendLists(page, el("line chart viewer"), 12));
      await session.step(156, "And the \"lines\" reading of line chart viewer should be 11", () => readingIs(page, "lines", el("line chart viewer"), 11));
      await session.step(157, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"HEIGHT / AS\" legend item in legend of line chart viewer"],["\"HEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"HEIGHT / PsA\" legend item in legend of line chart viewer"],["\"HEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"HEIGHT / RA\" legend item in legend of line chart viewer"],["\"HEIGHT / UC\" legend item in legend of line chart viewer"]]), [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"HEIGHT / AS\" legend item in legend of line chart viewer"],["\"HEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"HEIGHT / PsA\" legend item in legend of line chart viewer"],["\"HEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"HEIGHT / RA\" legend item in legend of line chart viewer"],["\"HEIGHT / UC\" legend item in legend of line chart viewer"]]);
      await session.step(170, "And the \"WEIGHT / RA\" and \"HEIGHT / RA\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "WEIGHT / RA", "HEIGHT / RA", el("line chart viewer")));
      await session.step(171, "And the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(172, "And line chart viewer should have a \"y2 axis\" area", () => hasArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(173, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Three Y columns under Multi Axis share one scale", async () => {
      await session.step(176, "When user sets \"yColumnNames\" property of line chart viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(177, "Then the \"y columns\" reading of line chart viewer should be \"AGE, HEIGHT, WEIGHT\"", () => readingReads(page, "y columns", el("line chart viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(178, "And the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(179, "And the \"y axes\" reading of line chart viewer should be 0", () => readingIs(page, "y axes", el("line chart viewer"), 0));
      await session.step(180, "And line chart viewer should not have a \"y2 axis\" area", () => hasNoArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(181, "And the legend of line chart viewer should list 18 items", () => legendLists(page, el("line chart viewer"), 18));
      await session.step(182, "And \"AGE / RA\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"AGE / RA\" legend item in legend of line chart viewer"), "visible"));
      await session.step(183, "And \"HEIGHT / RA\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"HEIGHT / RA\" legend item in legend of line chart viewer"), "visible"));
      await session.step(184, "And \"WEIGHT / RA\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"WEIGHT / RA\" legend item in legend of line chart viewer"), "visible"));
      await session.step(185, "When user sets \"yColumnNames\" property of line chart viewer to \"WEIGHT, HEIGHT\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "WEIGHT, HEIGHT"));
      await session.step(186, "Then the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(187, "And line chart viewer should have a \"y2 axis\" area", () => hasArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(188, "And the legend of line chart viewer should list 12 items", () => legendLists(page, el("line chart viewer"), 12));
      await session.step(189, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Replacing a Y column replaces its block and leaves the other one", async () => {
      await session.step(192, "Then the \"WEIGHT / RA\" item in the legend of line chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "WEIGHT / RA", el("line chart viewer"), "#9467BD"));
      await session.step(193, "And the \"WEIGHT / AS\" item in the legend of line chart viewer should be colored \"#1F77B4\"", () => legendItemColor(page, "WEIGHT / AS", el("line chart viewer"), "#1F77B4"));
      await session.step(194, "When user sets \"yColumnNames\" property of line chart viewer to \"WEIGHT, AGE\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "WEIGHT, AGE"));
      await session.step(195, "Then the \"y columns\" reading of line chart viewer should be \"WEIGHT, AGE\"", () => readingReads(page, "y columns", el("line chart viewer"), "WEIGHT, AGE"));
      await session.step(196, "And the legend of line chart viewer should list 12 items", () => legendLists(page, el("line chart viewer"), 12));
      await session.step(197, "And the \"lines\" reading of line chart viewer should be 12", () => readingIs(page, "lines", el("line chart viewer"), 12));
      await session.step(198, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]), [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]);
      await session.step(211, "And \"HEIGHT / RA\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("\"HEIGHT / RA\" legend item in legend of line chart viewer"), "absent"));
      await session.step(212, "And the \"WEIGHT / RA\" item in the legend of line chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "WEIGHT / RA", el("line chart viewer"), "#9467BD"));
      await session.step(213, "And the \"WEIGHT / AS\" item in the legend of line chart viewer should be colored \"#1F77B4\"", () => legendItemColor(page, "WEIGHT / AS", el("line chart viewer"), "#1F77B4"));
      await session.step(214, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Split, Multi Axis and the Y blocks come back from a saved layout", async () => {
      await session.step(217, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(218, "And user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["multiAxis","false"],["splitColumnNames",""],["yColumnNames","WEIGHT"]]), [["multiAxis","false"],["splitColumnNames",""],["yColumnNames","WEIGHT"]]);
      await session.step(222, "Then \"multiAxis\" property of line chart viewer should be \"false\"", () => propertyShouldBe(page, "multiAxis", el("line chart viewer"), "false"));
      await session.step(223, "And the \"y columns\" reading of line chart viewer should be \"WEIGHT\"", () => readingReads(page, "y columns", el("line chart viewer"), "WEIGHT"));
      await session.step(224, "And the legend of line chart viewer should be hidden", () => shouldBe(page, el("the legend of line chart viewer"), "hidden"));
      await session.step(225, "When user loads the saved layout", () => loadLayout(page));
      await session.step(226, "Then \"multiAxis\" property of line chart viewer should be \"true\"", () => propertyShouldBe(page, "multiAxis", el("line chart viewer"), "true"));
      await session.step(227, "And the \"y columns\" reading of line chart viewer should be \"WEIGHT, AGE\"", () => readingReads(page, "y columns", el("line chart viewer"), "WEIGHT, AGE"));
      await session.step(228, "And the \"lines\" reading of line chart viewer should be 12", () => readingIs(page, "lines", el("line chart viewer"), 12));
      await session.step(229, "And the legend of line chart viewer should list 12 items", () => legendLists(page, el("line chart viewer"), 12));
      await session.step(230, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]), [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]);
      await session.step(243, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Split, Multi Axis and the Y blocks come back from a saved project", async () => {
      await session.step(246, "When user saves the current view as project \"bdd-line-chart-legend-blocks\"", () => saveAsProject(page, "bdd-line-chart-legend-blocks"));
      await session.step(247, "And user closes all views", () => closeAllViews(page));
      await session.step(248, "And user opens the \"bdd-line-chart-legend-blocks\" project", () => openProject(page, "bdd-line-chart-legend-blocks"));
      await session.step(249, "Then line chart viewer should be added to the open tableview", () => viewerAdded(page, "line chart"));
      await session.step(250, "And \"multiAxis\" property of line chart viewer should be \"true\"", () => propertyShouldBe(page, "multiAxis", el("line chart viewer"), "true"));
      await session.step(251, "And the \"y columns\" reading of line chart viewer should be \"WEIGHT, AGE\"", () => readingReads(page, "y columns", el("line chart viewer"), "WEIGHT, AGE"));
      await session.step(252, "And the \"lines\" reading of line chart viewer should be 12", () => readingIs(page, "lines", el("line chart viewer"), 12));
      await session.step(253, "And the legend of line chart viewer should list 12 items", () => legendLists(page, el("line chart viewer"), 12));
      await session.step(254, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]), [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]);
      await session.step(267, "And \"HEIGHT / RA\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("\"HEIGHT / RA\" legend item in legend of line chart viewer"), "absent"));
      await session.step(268, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An explicit marker column keeps its legend items above the Auto markers threshold", async () => {
      await session.step(271, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["multiAxis","false"],["yColumnNames","AGE"],["splitColumnNames","RACE"],["markersColumnName","SEX"]]), [["multiAxis","false"],["yColumnNames","AGE"],["splitColumnNames","RACE"],["markersColumnName","SEX"]]);
      await session.step(276, "Then the legend of line chart viewer should list 6 items", () => legendLists(page, el("line chart viewer"), 6));
      await session.step(277, "And \"F\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"F\" legend item in legend of line chart viewer"), "visible"));
      await session.step(278, "And \"M\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"M\" legend item in legend of line chart viewer"), "visible"));
      await session.step(279, "When user filters rows where \"RACE\" is \"Other\"", () => filterTo(page, "RACE", "Other"));
      await session.step(280, "Then the legend of line chart viewer should list 3 items", () => legendLists(page, el("line chart viewer"), 3));
      await session.step(281, "And \"F\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"F\" legend item in legend of line chart viewer"), "visible"));
      await session.step(282, "And \"M\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"M\" legend item in legend of line chart viewer"), "visible"));
      await session.step(283, "When user resets the filter", () => resetFilter(page));
      await session.step(284, "And user sets \"markersColumnName\" property of line chart viewer to \"RACE\"", () => setProperty(page, "markersColumnName", el("line chart viewer"), "RACE"));
      await session.step(285, "Then the legend of line chart viewer should list 4 items", () => legendLists(page, el("line chart viewer"), 4));
      await session.step(286, "And marker of \"Black\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("marker of \"Black\" legend item in legend of line chart viewer"), "visible"));
      await session.step(287, "And \"M\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("\"M\" legend item in legend of line chart viewer"), "absent"));
      await session.step(288, "When user sets \"markersColumnName\" property of line chart viewer to \"\"", () => setProperty(page, "markersColumnName", el("line chart viewer"), ""));
      await session.step(289, "Then the legend of line chart viewer should list 4 items", () => legendLists(page, el("line chart viewer"), 4));
      await session.step(290, "And marker of \"Black\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("marker of \"Black\" legend item in legend of line chart viewer"), "absent"));
      await session.step(291, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Under a filter a two-column split recolors its legend as it recolors the lines, and a picked color stays", async () => {
      await session.step(294, "When user sets \"splitColumnNames\" property of line chart viewer to \"RACE, SEVERITY\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "RACE, SEVERITY"));
      await session.step(295, "Then the legend of line chart viewer should list 17 items", () => legendLists(page, el("line chart viewer"), 17));
      await session.step(296, "And the \"Asian / High\" item in the legend of line chart viewer should be colored \"#1F77B4\"", () => legendItemColor(page, "Asian / High", el("line chart viewer"), "#1F77B4"));
      await session.step(297, "When user filters rows where \"RACE\" is \"Other\"", () => filterTo(page, "RACE", "Other"));
      await session.step(298, "Then the legend of line chart viewer should list 4 items", () => legendLists(page, el("line chart viewer"), 4));
      await session.step(299, "And the \"Other / High\" item in the legend of line chart viewer should be colored \"#1F77B4\"", () => legendItemColor(page, "Other / High", el("line chart viewer"), "#1F77B4"));
      await session.step(300, "And the \"Other / Low\" item in the legend of line chart viewer should be colored \"#FF7F0E\"", () => legendItemColor(page, "Other / Low", el("line chart viewer"), "#FF7F0E"));
      await session.step(301, "And the \"Other / Medium\" item in the legend of line chart viewer should be colored \"#2CA02C\"", () => legendItemColor(page, "Other / Medium", el("line chart viewer"), "#2CA02C"));
      await session.step(302, "And the \"Other / None\" item in the legend of line chart viewer should be colored \"#D62728\"", () => legendItemColor(page, "Other / None", el("line chart viewer"), "#D62728"));
      await session.step(303, "And the \"view\" area of line chart viewer should contain the color \"#1F77B4\"", () => areaColor(page, "view", el("line chart viewer"), "#1F77B4"));
      await session.step(304, "And the \"view\" area of line chart viewer should contain the color \"#D62728\"", () => areaColor(page, "view", el("line chart viewer"), "#D62728"));
      await session.step(305, "When user hovers over \"Other / None\" legend item in legend of line chart viewer", () => hoverOver(page, el("\"Other / None\" legend item in legend of line chart viewer")));
      await session.step(306, "And user clicks on color picker icon", () => clickOn(page, el("color picker icon")));
      await session.step(307, "Then \"Other / None\" dialog should be visible", () => shouldBe(page, el("\"Other / None\" dialog"), "visible"));
      await session.step(308, "When user picks the color \"#BCBD22\" in the color picker dialog", () => pickColorSwatch(page, "#BCBD22"));
      await session.step(309, "And user clicks on OK button in \"Other / None\" dialog", () => clickOn(page, el("OK button in \"Other / None\" dialog")));
      await session.step(310, "Then the \"Other / None\" item in the legend of line chart viewer should be colored \"#BCBD22\"", () => legendItemColor(page, "Other / None", el("line chart viewer"), "#BCBD22"));
      await session.step(311, "And the \"view\" area of line chart viewer should contain the color \"#BCBD22\"", () => areaColor(page, "view", el("line chart viewer"), "#BCBD22"));
      await session.step(312, "When user resets the filter", () => resetFilter(page));
      await session.step(313, "Then the legend of line chart viewer should list 17 items", () => legendLists(page, el("line chart viewer"), 17));
      await session.step(314, "And the \"Other / None\" item in the legend of line chart viewer should be colored \"#BCBD22\"", () => legendItemColor(page, "Other / None", el("line chart viewer"), "#BCBD22"));
      await session.step(315, "When user resizes line chart viewer to 700 by 420", () => resizeTo(page, el("line chart viewer"), 700, 420));
      await session.step(316, "Then the \"Other / None\" item in the legend of line chart viewer should be colored \"#BCBD22\"", () => legendItemColor(page, "Other / None", el("line chart viewer"), "#BCBD22"));
      await session.step(317, "When user restores the size of line chart viewer", () => restoreSize(page, el("line chart viewer")));
      await session.step(318, "And user sets \"splitColumnNames\" property of line chart viewer to \"DIS_POP\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "DIS_POP"));
      await session.step(319, "Then the legend of line chart viewer should list 6 items", () => legendLists(page, el("line chart viewer"), 6));
      await session.step(320, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
