/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-styling.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart, core#825, public#3463, public#1651, public#3210, public#3709]
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
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, clearSelection, removeColumn, selectWhereIs, selectWhereOneOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaShorter, areaTaller, closeContextMenu, hoverArea, legendLists, lessInk, moreInk, noErrors, notRepainted, painted, pickFromAreaContextMenu, pointerAway, propertyShouldBe, readingFinite, readingIs, readingReads, repainted, reportsNoError, setProperties, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart axis types, label orientation, line styling and a chart type per chart", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-styling.feature", import.meta.url);
  test("Line chart axis types, label orientation, line styling and a chart type per chart", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart", "@realizes:core#825", "@realizes:public#3463", "@realizes:public#1651", "@realizes:public#3210", "@realizes:public#3709"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"]]);
    await session.step(26, "Then line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(27, "And the 'y axis max of \"WEIGHT\"' reading of line chart viewer should be a finite number", () => readingFinite(page, "y axis max of \"WEIGHT\"", el("line chart viewer")));
    await run.scenario("A logarithmic Y axis redraws the chart and a linear one draws it back", async () => {
      await session.step(30, "When user sets \"yAxisType\" property of line chart viewer to \"logarithmic\"", () => setProperty(page, "yAxisType", el("line chart viewer"), "logarithmic"));
      await session.step(31, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(32, "And the 'y axis max of \"WEIGHT\"' reading of line chart viewer should be a finite number", () => readingFinite(page, "y axis max of \"WEIGHT\"", el("line chart viewer")));
      await session.step(33, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(34, "When user sets \"yAxisType\" property of line chart viewer to \"linear\"", () => setProperty(page, "yAxisType", el("line chart viewer"), "linear"));
      await session.step(35, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(36, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Vertical X axis labels take a taller axis strip, Auto gives the room back", async () => {
      await session.step(40, "When user sets \"xAxisLabelOrientation\" property of line chart viewer to \"Vert\"", () => setProperty(page, "xAxisLabelOrientation", el("line chart viewer"), "Vert"));
      await session.step(41, "Then the \"x axis\" area of line chart viewer should be taller than before", () => areaTaller(page, "x axis", el("line chart viewer")));
      await session.step(42, "When user sets \"xAxisLabelOrientation\" property of line chart viewer to \"Auto\"", () => setProperty(page, "xAxisLabelOrientation", el("line chart viewer"), "Auto"));
      await session.step(43, "Then the \"x axis\" area of line chart viewer should be shorter than before", () => areaShorter(page, "x axis", el("line chart viewer")));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Line width, transparency and colouring type restyle the line", async () => {
      await session.step(47, "When user sets \"lineWidth\" property of line chart viewer to \"3\"", () => setProperty(page, "lineWidth", el("line chart viewer"), "3"));
      await session.step(48, "Then line chart viewer should have more ink than before", () => moreInk(page, el("line chart viewer")));
      await session.step(49, "When user sets \"lineTransparency\" property of line chart viewer to \"0.5\"", () => setProperty(page, "lineTransparency", el("line chart viewer"), "0.5"));
      await session.step(50, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(51, "When user sets \"lineColoringType\" property of line chart viewer to \"Custom\"", () => setProperty(page, "lineColoringType", el("line chart viewer"), "Custom"));
      await session.step(52, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(53, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(54, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["lineWidth","1"],["lineTransparency","0"],["lineColoringType","Auto"]]), [["lineWidth","1"],["lineTransparency","0"],["lineColoringType","Auto"]]);
      await session.step(58, "Then line chart viewer should have less ink than before", () => lessInk(page, el("line chart viewer")));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("In a multi-axis chart one chart's type changes through its own menu and the others keep theirs", async () => {
      await session.step(62, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","AGE, HEIGHT, WEIGHT"],["multiAxis","true"]]), [["yColumnNames","AGE, HEIGHT, WEIGHT"],["multiAxis","true"]]);
      await session.step(65, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(66, "And \"chartTypes\" property of line chart viewer should be \"Line Chart, Line Chart, Line Chart\"", () => propertyShouldBe(page, "chartTypes", el("line chart viewer"), "Line Chart, Line Chart, Line Chart"));
      await session.step(67, "When user picks \"HEIGHT > Chart type > Area Chart\" from the context menu of the \"plot\" area of line chart viewer", () => pickFromAreaContextMenu(page, "HEIGHT > Chart type > Area Chart", "plot", el("line chart viewer")));
      await session.step(68, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(69, "Then \"chartTypes\" property of line chart viewer should be \"Line Chart, Area Chart, Line Chart\"", () => propertyShouldBe(page, "chartTypes", el("line chart viewer"), "Line Chart, Area Chart, Line Chart"));
      await session.step(70, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(71, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["multiAxis","false"],["yColumnNames","WEIGHT"]]), [["multiAxis","false"],["yColumnNames","WEIGHT"]]);
      await session.step(74, "Then the \"y columns\" reading of line chart viewer should be \"WEIGHT\"", () => readingReads(page, "y columns", el("line chart viewer"), "WEIGHT"));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The split legend lists no Y / split pair that has no data (core#825)", async () => {
      await session.step(79, "When user adds a calculated column \"H_F\" with formula \"if(${SEX} == \\\"F\\\", ${HEIGHT}, null)\"", () => addCalculated(page, "H_F", "if(${SEX} == \"F\", ${HEIGHT}, null)"));
      await session.step(80, "And user adds a calculated column \"W_M\" with formula \"if(${SEX} == \\\"M\\\", ${WEIGHT}, null)\"", () => addCalculated(page, "W_M", "if(${SEX} == \"M\", ${WEIGHT}, null)"));
      await session.step(81, "And user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","H_F, W_M"],["splitColumnNames","SEX"],["multiAxis","true"],["yGlobalScale","true"]]), [["yColumnNames","H_F, W_M"],["splitColumnNames","SEX"],["multiAxis","true"],["yGlobalScale","true"]]);
      await session.step(86, "Then the legend of line chart viewer should list 2 items", () => legendLists(page, el("line chart viewer"), 2));
      await session.step(87, "And \"H_F / F\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"H_F / F\" legend item in legend of line chart viewer"), "visible"));
      await session.step(88, "And \"W_M / M\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"W_M / M\" legend item in legend of line chart viewer"), "visible"));
      await session.step(89, "And \"H_F / M\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("\"H_F / M\" legend item in legend of line chart viewer"), "absent"));
      await session.step(90, "And \"W_M / F\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("\"W_M / F\" legend item in legend of line chart viewer"), "absent"));
      await session.step(91, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","WEIGHT"],["splitColumnNames",""],["multiAxis","false"],["yGlobalScale","false"]]), [["yColumnNames","WEIGHT"],["splitColumnNames",""],["multiAxis","false"],["yGlobalScale","false"]]);
      await session.step(96, "And user removes \"H_F\" column", () => removeColumn(page, "H_F"));
      await session.step(97, "And user removes \"W_M\" column", () => removeColumn(page, "W_M"));
      await session.step(98, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Under Row Source = Selected the split legend lists the categories of the selected rows (public#1651)", async () => {
      await session.step(102, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["splitColumnNames","RACE"],["rowSource","Selected"]]), [["splitColumnNames","RACE"],["rowSource","Selected"]]);
      await session.step(105, "And user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(106, "Then the legend of line chart viewer should list 1 item", () => legendLists(page, el("line chart viewer"), 1));
      await session.step(107, "And \"Asian\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"Asian\" legend item in legend of line chart viewer"), "visible"));
      await session.step(108, "When user selects rows where \"RACE\" is one of \"Asian, Black\"", () => selectWhereOneOf(page, "RACE", "Asian, Black"));
      await session.step(109, "Then the legend of line chart viewer should list 2 items", () => legendLists(page, el("line chart viewer"), 2));
      await session.step(110, "And \"Black\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"Black\" legend item in legend of line chart viewer"), "visible"));
      await session.step(111, "And \"Caucasian\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("\"Caucasian\" legend item in legend of line chart viewer"), "absent"));
      await session.step(112, "When user clears the row selection", () => clearSelection(page));
      await session.step(113, "And user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["splitColumnNames",""],["rowSource","Filtered"]]), [["splitColumnNames",""],["rowSource","Filtered"]]);
      await session.step(116, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A split with a category per row keeps the chart working (public#3210)", async () => {
      await session.step(120, "When user sets \"splitColumnNames\" property of line chart viewer to \"USUBJID\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "USUBJID"));
      await session.step(121, "Then the \"categories\" reading of line chart viewer should be 1000", () => readingIs(page, "categories", el("line chart viewer"), 1000));
      await session.step(122, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(123, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(124, "Then the \"categories\" reading of line chart viewer should be 1", () => readingIs(page, "categories", el("line chart viewer"), 1));
      await session.step(125, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(126, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A row group hovered in a bar chart repaints the line chart only while Show Mouse Over Row Group is on (public#3709)", async () => {
      await session.step(130, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
      await session.step(132, "And user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(133, "And user hovers over the \"bar Asian\" area of bar chart viewer", () => hoverArea(page, "bar Asian", el("bar chart viewer")));
      await session.step(134, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(135, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(136, "And user sets \"showMouseOverRowGroup\" property of line chart viewer to \"false\"", () => setProperty(page, "showMouseOverRowGroup", el("line chart viewer"), "false"));
      await session.step(137, "And user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(138, "And user hovers over the \"bar Asian\" area of bar chart viewer", () => hoverArea(page, "bar Asian", el("bar chart viewer")));
      await session.step(139, "Then bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
      await session.step(140, "And line chart viewer should not have repainted", () => notRepainted(page, el("line chart viewer")));
      await session.step(141, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(142, "And user sets \"showMouseOverRowGroup\" property of line chart viewer to \"true\"", () => setProperty(page, "showMouseOverRowGroup", el("line chart viewer"), "true"));
      await session.step(143, "And user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(144, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
