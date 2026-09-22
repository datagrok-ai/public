/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-styling.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaShorter, areaTaller, lessInk, moreInk, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingFinite, readingIs, readingReads, repainted, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart axis types, label orientation, line styling and a chart type per chart", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-styling.feature", import.meta.url);
  test("Line chart axis types, label orientation, line styling and a chart type per chart", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"]]);
    await session.step(20, "Then line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(21, "And the 'y axis max of \"WEIGHT\"' reading of line chart viewer should be a finite number", () => readingFinite(page, "y axis max of \"WEIGHT\"", el("line chart viewer")));
    await run.scenario("A logarithmic Y axis redraws the chart and a linear one draws it back", async () => {
      await session.step(24, "When user sets \"yAxisType\" property of line chart viewer to \"logarithmic\"", () => setProperty(page, "yAxisType", el("line chart viewer"), "logarithmic"));
      await session.step(25, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(26, "And the 'y axis max of \"WEIGHT\"' reading of line chart viewer should be a finite number", () => readingFinite(page, "y axis max of \"WEIGHT\"", el("line chart viewer")));
      await session.step(27, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(28, "When user sets \"yAxisType\" property of line chart viewer to \"linear\"", () => setProperty(page, "yAxisType", el("line chart viewer"), "linear"));
      await session.step(29, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(30, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Vertical X axis labels take a taller axis strip, Auto gives the room back", async () => {
      await session.step(34, "When user sets \"xAxisLabelOrientation\" property of line chart viewer to \"Vert\"", () => setProperty(page, "xAxisLabelOrientation", el("line chart viewer"), "Vert"));
      await session.step(35, "Then the \"x axis\" area of line chart viewer should be taller than before", () => areaTaller(page, "x axis", el("line chart viewer")));
      await session.step(36, "When user sets \"xAxisLabelOrientation\" property of line chart viewer to \"Auto\"", () => setProperty(page, "xAxisLabelOrientation", el("line chart viewer"), "Auto"));
      await session.step(37, "Then the \"x axis\" area of line chart viewer should be shorter than before", () => areaShorter(page, "x axis", el("line chart viewer")));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Line width, transparency and colouring type restyle the line", async () => {
      await session.step(41, "When user sets \"lineWidth\" property of line chart viewer to \"3\"", () => setProperty(page, "lineWidth", el("line chart viewer"), "3"));
      await session.step(42, "Then line chart viewer should have more ink than before", () => moreInk(page, el("line chart viewer")));
      await session.step(43, "When user sets \"lineTransparency\" property of line chart viewer to \"0.5\"", () => setProperty(page, "lineTransparency", el("line chart viewer"), "0.5"));
      await session.step(44, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(45, "When user sets \"lineColoringType\" property of line chart viewer to \"Custom\"", () => setProperty(page, "lineColoringType", el("line chart viewer"), "Custom"));
      await session.step(46, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(47, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(48, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["lineWidth","1"],["lineTransparency","0"],["lineColoringType","Auto"]]), [["lineWidth","1"],["lineTransparency","0"],["lineColoringType","Auto"]]);
      await session.step(52, "Then line chart viewer should have less ink than before", () => lessInk(page, el("line chart viewer")));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("In a multi-axis chart one chart's type changes through its own menu and the others keep theirs", async () => {
      await session.step(56, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","AGE, HEIGHT, WEIGHT"],["multiAxis","true"]]), [["yColumnNames","AGE, HEIGHT, WEIGHT"],["multiAxis","true"]]);
      await session.step(59, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(60, "And \"chartTypes\" property of line chart viewer should be \"Line Chart, Line Chart, Line Chart\"", () => propertyShouldBe(page, "chartTypes", el("line chart viewer"), "Line Chart, Line Chart, Line Chart"));
      await session.step(61, "When user picks \"HEIGHT > Chart type > Area Chart\" from the context menu of the \"plot\" area of line chart viewer", () => pickFromAreaContextMenu(page, "HEIGHT > Chart type > Area Chart", "plot", el("line chart viewer")));
      await session.step(62, "Then \"chartTypes\" property of line chart viewer should be \"Line Chart, Area Chart, Line Chart\"", () => propertyShouldBe(page, "chartTypes", el("line chart viewer"), "Line Chart, Area Chart, Line Chart"));
      await session.step(63, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(64, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["multiAxis","false"],["yColumnNames","WEIGHT"]]), [["multiAxis","false"],["yColumnNames","WEIGHT"]]);
      await session.step(67, "Then the \"y columns\" reading of line chart viewer should be \"WEIGHT\"", () => readingReads(page, "y columns", el("line chart viewer"), "WEIGHT"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
