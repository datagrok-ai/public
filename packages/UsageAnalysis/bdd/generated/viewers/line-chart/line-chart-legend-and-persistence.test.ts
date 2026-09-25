/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-legend-and-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/nx.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {followingShouldBe, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {categoricalColorIs, colorCategorical, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickLegendItem, hasArea, hasNoArea, legendItemColor, legendItemsDiffer, legendLists, legendSide, loadLayout, noErrors, painted, propertyShouldBe, readingIs, readingReads, repainted, reportsNoError, saveLayoutToServer, setProperties, setProperty, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart legend, category colours and what survives a round-trip", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-legend-and-persistence.feature", import.meta.url);
  test("Line chart legend, category colours and what survives a round-trip", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 13, page);
    await session.step(33, "Given user is logged in", () => loggedIn(page));
    await session.step(34, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(35, "And user colors \"Stereo Category\" column categorically:", () => colorCategorical(page, "Stereo Category", [["R_ONE","#FF0000"],["S_ABS","#00FF00"],["S_ACHIR","#0000FF"],["S_PART","#FFFF00"],["S_UNKN","#FF00FF"]]), [["R_ONE","#FF0000"],["S_ABS","#00FF00"],["S_ACHIR","#0000FF"],["S_PART","#FFFF00"],["S_UNKN","#FF00FF"]]);
    await session.step(41, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","Chemical Space X"],["yColumnNames","Chemical Space Y"],["splitColumnNames","Stereo Category"]]), [["xColumnName","Chemical Space X"],["yColumnNames","Chemical Space Y"],["splitColumnNames","Stereo Category"]]);
    await session.step(45, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(46, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
    await session.step(47, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
    await session.step(48, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
    await session.step(49, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("The legend comes with the split and goes when the split does", async () => {
      await session.step(52, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(53, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(54, "Then the legend of line chart viewer should be hidden", () => shouldBe(page, el("the legend of line chart viewer"), "hidden"));
      await session.step(55, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(56, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(57, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Position moves the legend to the side it names", async () => {
      await session.step(61, "When user sets \"legendPosition\" property of line chart viewer to \"Left\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Left"));
      await session.step(62, "Then the legend of line chart viewer should be on the left", () => legendSide(page, el("line chart viewer"), "left"));
      await session.step(63, "When user sets \"legendPosition\" property of line chart viewer to \"Top\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Top"));
      await session.step(64, "Then the legend of line chart viewer should be on the top", () => legendSide(page, el("line chart viewer"), "top"));
      await session.step(65, "When user sets \"legendPosition\" property of line chart viewer to \"Bottom\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Bottom"));
      await session.step(66, "Then the legend of line chart viewer should be on the bottom", () => legendSide(page, el("line chart viewer"), "bottom"));
      await session.step(67, "When user sets \"legendPosition\" property of line chart viewer to \"Right\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Right"));
      await session.step(68, "Then the legend of line chart viewer should be on the right", () => legendSide(page, el("line chart viewer"), "right"));
      await session.step(69, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(70, "When user sets \"legendPosition\" property of line chart viewer to \"Auto\"", () => setProperty(page, "legendPosition", el("line chart viewer"), "Auto"));
      await session.step(71, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Visibility Never hides the legend and Auto brings it back", async () => {
      await session.step(74, "When user sets \"legendVisibility\" property of line chart viewer to \"Never\"", () => setProperty(page, "legendVisibility", el("line chart viewer"), "Never"));
      await session.step(75, "Then the legend of line chart viewer should be hidden", () => shouldBe(page, el("the legend of line chart viewer"), "hidden"));
      await session.step(76, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(77, "When user sets \"legendVisibility\" property of line chart viewer to \"Always\"", () => setProperty(page, "legendVisibility", el("line chart viewer"), "Always"));
      await session.step(78, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(79, "When user sets \"legendVisibility\" property of line chart viewer to \"Auto\"", () => setProperty(page, "legendVisibility", el("line chart viewer"), "Auto"));
      await session.step(80, "Then the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking a legend category filters the chart down to it and back", async () => {
      await session.step(84, "When user clicks on \"S_ABS\" item in the legend of line chart viewer", () => clickLegendItem(page, "S_ABS", el("line chart viewer")));
      await session.step(85, "Then the \"rows shown\" reading of line chart viewer should be 2", () => readingIs(page, "rows shown", el("line chart viewer"), 2));
      await session.step(86, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(87, "And 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(88, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
      await session.step(89, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(90, "When user clicks on \"S_ABS\" item in the legend of line chart viewer", () => clickLegendItem(page, "S_ABS", el("line chart viewer")));
      await session.step(91, "Then the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(92, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The legend items keep the colours the column was given", async () => {
      await session.step(96, "Then the \"R_ONE\" item in the legend of line chart viewer should be colored \"#FF0000\"", () => legendItemColor(page, "R_ONE", el("line chart viewer"), "#FF0000"));
      await session.step(97, "And the \"S_ABS\" item in the legend of line chart viewer should be colored \"#00FF00\"", () => legendItemColor(page, "S_ABS", el("line chart viewer"), "#00FF00"));
      await session.step(98, "And the \"R_ONE\" and \"S_ABS\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "R_ONE", "S_ABS", el("line chart viewer")));
      await session.step(99, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The chart configuration and the category colour come back from a saved layout (GROK-17278)", async () => {
      await session.step(102, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(103, "And user colors \"Stereo Category\" column categorically:", () => colorCategorical(page, "Stereo Category", [["R_ONE","#00AAFF"]]), [["R_ONE","#00AAFF"]]);
      await session.step(105, "Then the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#00AAFF\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#00AAFF"));
      await session.step(106, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(107, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(108, "When user loads the saved layout", () => loadLayout(page));
      await session.step(109, "Then the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(110, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(111, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
      await session.step(112, "And the \"x column\" reading of line chart viewer should be \"Chemical Space X\"", () => readingReads(page, "x column", el("line chart viewer"), "Chemical Space X"));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The chart and its colours come back from a saved project (GROK-19825)", async () => {
      await session.step(116, "When user saves the current view as project \"zz-linechart-legend-colors\"", () => saveAsProject(page, "zz-linechart-legend-colors"));
      await session.step(117, "And user closes all views", () => closeAllViews(page));
      await session.step(118, "And user opens the \"zz-linechart-legend-colors\" project", () => openProject(page, "zz-linechart-legend-colors"));
      await session.step(119, "Then line chart viewer should be added to the open tableview", () => viewerAdded(page, "line chart"));
      await session.step(120, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(121, "And the legend of line chart viewer should list 5 items", () => legendLists(page, el("line chart viewer"), 5));
      await session.step(122, "And the categorical color of \"R_ONE\" in \"Stereo Category\" column should be \"#FF0000\"", () => categoricalColorIs(page, "R_ONE", "Stereo Category", "#FF0000"));
      await session.step(123, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(124, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The default palette gives neighbouring split categories different colors", async () => {
      await session.step(127, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(128, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","STARTED"],["yColumnNames","WEIGHT"],["splitColumnNames","DIS_POP"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["xColumnName","STARTED"],["yColumnNames","WEIGHT"],["splitColumnNames","DIS_POP"],["Legend Visibility","Always"],["Legend Position","Right"]]);
      await session.step(134, "Then the legend of line chart viewer should list 6 items", () => legendLists(page, el("line chart viewer"), 6));
      await session.step(135, "And the \"AS\" and \"Indigestion\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "AS", "Indigestion", el("line chart viewer")));
      await session.step(136, "And the \"Indigestion\" and \"PsA\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "Indigestion", "PsA", el("line chart viewer")));
      await session.step(137, "And the \"PsA\" and \"Psoriasis\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "PsA", "Psoriasis", el("line chart viewer")));
      await session.step(138, "And the \"Psoriasis\" and \"RA\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "Psoriasis", "RA", el("line chart viewer")));
      await session.step(139, "And the \"RA\" and \"UC\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "RA", "UC", el("line chart viewer")));
      await session.step(140, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Multi Axis gives every Y column its own block of legend items", async () => {
      await session.step(143, "When user sets \"yColumnNames\" property of line chart viewer to \"WEIGHT, HEIGHT\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "WEIGHT, HEIGHT"));
      await session.step(144, "Then the legend of line chart viewer should list 6 items", () => legendLists(page, el("line chart viewer"), 6));
      await session.step(145, "When user sets \"multiAxis\" property of line chart viewer to \"true\"", () => setProperty(page, "multiAxis", el("line chart viewer"), "true"));
      await session.step(146, "Then the legend of line chart viewer should list 12 items", () => legendLists(page, el("line chart viewer"), 12));
      await session.step(149, "And the \"lines\" reading of line chart viewer should be 11", () => readingIs(page, "lines", el("line chart viewer"), 11));
      await session.step(150, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"HEIGHT / AS\" legend item in legend of line chart viewer"],["\"HEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"HEIGHT / PsA\" legend item in legend of line chart viewer"],["\"HEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"HEIGHT / RA\" legend item in legend of line chart viewer"],["\"HEIGHT / UC\" legend item in legend of line chart viewer"]]), [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"HEIGHT / AS\" legend item in legend of line chart viewer"],["\"HEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"HEIGHT / PsA\" legend item in legend of line chart viewer"],["\"HEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"HEIGHT / RA\" legend item in legend of line chart viewer"],["\"HEIGHT / UC\" legend item in legend of line chart viewer"]]);
      await session.step(163, "And the \"WEIGHT / RA\" and \"HEIGHT / RA\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "WEIGHT / RA", "HEIGHT / RA", el("line chart viewer")));
      await session.step(164, "And the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(165, "And line chart viewer should have a \"y2 axis\" area", () => hasArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(166, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Three Y columns under Multi Axis share one scale", async () => {
      await session.step(169, "When user sets \"yColumnNames\" property of line chart viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(170, "Then the \"y columns\" reading of line chart viewer should be \"AGE, HEIGHT, WEIGHT\"", () => readingReads(page, "y columns", el("line chart viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(171, "And the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(172, "And the \"y axes\" reading of line chart viewer should be 0", () => readingIs(page, "y axes", el("line chart viewer"), 0));
      await session.step(173, "And line chart viewer should not have a \"y2 axis\" area", () => hasNoArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(174, "And the legend of line chart viewer should list 18 items", () => legendLists(page, el("line chart viewer"), 18));
      await session.step(175, "And \"AGE / RA\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"AGE / RA\" legend item in legend of line chart viewer"), "visible"));
      await session.step(176, "And \"HEIGHT / RA\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"HEIGHT / RA\" legend item in legend of line chart viewer"), "visible"));
      await session.step(177, "And \"WEIGHT / RA\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"WEIGHT / RA\" legend item in legend of line chart viewer"), "visible"));
      await session.step(178, "When user sets \"yColumnNames\" property of line chart viewer to \"WEIGHT, HEIGHT\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "WEIGHT, HEIGHT"));
      await session.step(179, "Then the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
      await session.step(180, "And line chart viewer should have a \"y2 axis\" area", () => hasArea(page, el("line chart viewer"), "y2 axis"));
      await session.step(181, "And the legend of line chart viewer should list 12 items", () => legendLists(page, el("line chart viewer"), 12));
      await session.step(182, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Replacing a Y column replaces its block and leaves the other one", async () => {
      await session.step(185, "Then the \"WEIGHT / RA\" item in the legend of line chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "WEIGHT / RA", el("line chart viewer"), "#9467BD"));
      await session.step(186, "And the \"WEIGHT / AS\" item in the legend of line chart viewer should be colored \"#1F77B4\"", () => legendItemColor(page, "WEIGHT / AS", el("line chart viewer"), "#1F77B4"));
      await session.step(187, "When user sets \"yColumnNames\" property of line chart viewer to \"WEIGHT, AGE\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "WEIGHT, AGE"));
      await session.step(188, "Then the \"y columns\" reading of line chart viewer should be \"WEIGHT, AGE\"", () => readingReads(page, "y columns", el("line chart viewer"), "WEIGHT, AGE"));
      await session.step(189, "And the legend of line chart viewer should list 12 items", () => legendLists(page, el("line chart viewer"), 12));
      await session.step(190, "And the \"lines\" reading of line chart viewer should be 12", () => readingIs(page, "lines", el("line chart viewer"), 12));
      await session.step(191, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]), [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]);
      await session.step(204, "And \"HEIGHT / RA\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("\"HEIGHT / RA\" legend item in legend of line chart viewer"), "absent"));
      await session.step(205, "And the \"WEIGHT / RA\" item in the legend of line chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "WEIGHT / RA", el("line chart viewer"), "#9467BD"));
      await session.step(206, "And the \"WEIGHT / AS\" item in the legend of line chart viewer should be colored \"#1F77B4\"", () => legendItemColor(page, "WEIGHT / AS", el("line chart viewer"), "#1F77B4"));
      await session.step(207, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Split, Multi Axis and the Y blocks come back from a saved layout", async () => {
      await session.step(210, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(211, "And user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["multiAxis","false"],["splitColumnNames",""],["yColumnNames","WEIGHT"]]), [["multiAxis","false"],["splitColumnNames",""],["yColumnNames","WEIGHT"]]);
      await session.step(215, "Then \"multiAxis\" property of line chart viewer should be \"false\"", () => propertyShouldBe(page, "multiAxis", el("line chart viewer"), "false"));
      await session.step(216, "And the \"y columns\" reading of line chart viewer should be \"WEIGHT\"", () => readingReads(page, "y columns", el("line chart viewer"), "WEIGHT"));
      await session.step(217, "And the legend of line chart viewer should be hidden", () => shouldBe(page, el("the legend of line chart viewer"), "hidden"));
      await session.step(218, "When user loads the saved layout", () => loadLayout(page));
      await session.step(219, "Then \"multiAxis\" property of line chart viewer should be \"true\"", () => propertyShouldBe(page, "multiAxis", el("line chart viewer"), "true"));
      await session.step(220, "And the \"y columns\" reading of line chart viewer should be \"WEIGHT, AGE\"", () => readingReads(page, "y columns", el("line chart viewer"), "WEIGHT, AGE"));
      await session.step(221, "And the \"lines\" reading of line chart viewer should be 12", () => readingIs(page, "lines", el("line chart viewer"), 12));
      await session.step(222, "And the legend of line chart viewer should list 12 items", () => legendLists(page, el("line chart viewer"), 12));
      await session.step(223, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]), [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]);
      await session.step(236, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Split, Multi Axis and the Y blocks come back from a saved project", async () => {
      await session.step(239, "When user saves the current view as project \"bdd-line-chart-legend-blocks\"", () => saveAsProject(page, "bdd-line-chart-legend-blocks"));
      await session.step(240, "And user closes all views", () => closeAllViews(page));
      await session.step(241, "And user opens the \"bdd-line-chart-legend-blocks\" project", () => openProject(page, "bdd-line-chart-legend-blocks"));
      await session.step(242, "Then line chart viewer should be added to the open tableview", () => viewerAdded(page, "line chart"));
      await session.step(243, "And \"multiAxis\" property of line chart viewer should be \"true\"", () => propertyShouldBe(page, "multiAxis", el("line chart viewer"), "true"));
      await session.step(244, "And the \"y columns\" reading of line chart viewer should be \"WEIGHT, AGE\"", () => readingReads(page, "y columns", el("line chart viewer"), "WEIGHT, AGE"));
      await session.step(245, "And the \"lines\" reading of line chart viewer should be 12", () => readingIs(page, "lines", el("line chart viewer"), 12));
      await session.step(246, "And the legend of line chart viewer should list 12 items", () => legendLists(page, el("line chart viewer"), 12));
      await session.step(247, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]), [["\"WEIGHT / AS\" legend item in legend of line chart viewer"],["\"WEIGHT / Indigestion\" legend item in legend of line chart viewer"],["\"WEIGHT / PsA\" legend item in legend of line chart viewer"],["\"WEIGHT / Psoriasis\" legend item in legend of line chart viewer"],["\"WEIGHT / RA\" legend item in legend of line chart viewer"],["\"WEIGHT / UC\" legend item in legend of line chart viewer"],["\"AGE / AS\" legend item in legend of line chart viewer"],["\"AGE / Indigestion\" legend item in legend of line chart viewer"],["\"AGE / PsA\" legend item in legend of line chart viewer"],["\"AGE / Psoriasis\" legend item in legend of line chart viewer"],["\"AGE / RA\" legend item in legend of line chart viewer"],["\"AGE / UC\" legend item in legend of line chart viewer"]]);
      await session.step(260, "And \"HEIGHT / RA\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("\"HEIGHT / RA\" legend item in legend of line chart viewer"), "absent"));
      await session.step(261, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
