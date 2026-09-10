/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaNarrower, areaPainted, areaRepainted, areaShorter, closeContextMenu, hasArea, hasNoArea, menuLists, moreInk, noErrors, openContextMenu, painted, pickFromContextMenu, propertyShouldBe, propertyShouldContain, readingIs, readingReads, repainted, reportsError, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {descriptionAbove, descriptionBelow} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart chrome, chart types and the empty chart", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart.feature", import.meta.url);
  test("Line chart chrome, chart types and the empty chart", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(19, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]));
    await session.step(22, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(23, "And the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
    await session.step(24, "And the \"x column\" reading of line chart viewer should be \"CAST Idea ID\"", () => readingReads(page, "x column", el("line chart viewer"), "CAST Idea ID"));
    await session.step(25, "And the \"y columns\" reading of line chart viewer should be \"Chemical Space X\"", () => readingReads(page, "y columns", el("line chart viewer"), "Chemical Space X"));
    await session.step(26, "And the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
    await session.step(27, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
    await session.step(28, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(29, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(30, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
    await run.scenario("With no Y column the chart says so and draws no chart box", async () => {
      await session.step(33, "Then line chart viewer should have a \"chart 1\" area", () => hasArea(page, el("line chart viewer"), "chart 1"));
      await session.step(34, "When user sets \"yColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), ""));
      await session.step(35, "Then line chart viewer should report the error \"No Y columns selected\"", () => reportsError(page, el("line chart viewer"), "No Y columns selected"));
      await session.step(36, "And line chart viewer should not have a \"chart 1\" area", () => hasNoArea(page, el("line chart viewer"), "chart 1"));
      await session.step(37, "And the \"charts\" reading of line chart viewer should be 0", () => readingIs(page, "charts", el("line chart viewer"), 0));
      await session.step(38, "And the \"lines\" reading of line chart viewer should be 0", () => readingIs(page, "lines", el("line chart viewer"), 0));
      await session.step(39, "And the \"markers drawn\" reading of line chart viewer should be 0", () => readingIs(page, "markers drawn", el("line chart viewer"), 0));
      await session.step(40, "And the \"y axes\" reading of line chart viewer should be 0", () => readingIs(page, "y axes", el("line chart viewer"), 0));
      await session.step(41, "When user sets \"yColumnNames\" property of line chart viewer to \"Chemical Space X\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "Chemical Space X"));
      await session.step(42, "Then line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(43, "And line chart viewer should have a \"chart 1\" area", () => hasArea(page, el("line chart viewer"), "chart 1"));
      await session.step(44, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Chart Type redraws the same five series four ways", async () => {
      await session.step(48, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(49, "Then the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(50, "When user picks \"Chart Type > Area Chart\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Chart Type > Area Chart", el("line chart viewer")));
      await session.step(51, "Then \"chartTypes\" property of line chart viewer should contain \"Area Chart\"", () => propertyShouldContain(page, "chartTypes", el("line chart viewer"), "Area Chart"));
      await session.step(52, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(53, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(54, "When user picks \"Chart Type > Stacked Area Chart\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Chart Type > Stacked Area Chart", el("line chart viewer")));
      await session.step(55, "Then \"chartTypes\" property of line chart viewer should contain \"Stacked Area Chart\"", () => propertyShouldContain(page, "chartTypes", el("line chart viewer"), "Stacked Area Chart"));
      await session.step(56, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(57, "When user picks \"Chart Type > Stacked Bar Chart\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Chart Type > Stacked Bar Chart", el("line chart viewer")));
      await session.step(58, "Then \"chartTypes\" property of line chart viewer should contain \"Stacked Bar Chart\"", () => propertyShouldContain(page, "chartTypes", el("line chart viewer"), "Stacked Bar Chart"));
      await session.step(59, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(60, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(61, "When user picks \"Chart Type > Line Chart\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Chart Type > Line Chart", el("line chart viewer")));
      await session.step(62, "Then \"chartTypes\" property of line chart viewer should contain \"Line Chart\"", () => propertyShouldContain(page, "chartTypes", el("line chart viewer"), "Line Chart"));
      await session.step(63, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(64, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(65, "Then the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Interpolation and line width repaint the same hundred points", async () => {
      await session.step(69, "When user sets \"interpolation\" property of line chart viewer to \"Spline\"", () => setProperty(page, "interpolation", el("line chart viewer"), "Spline"));
      await session.step(70, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(71, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(72, "When user sets \"splineTension\" property of line chart viewer to \"1\"", () => setProperty(page, "splineTension", el("line chart viewer"), "1"));
      await session.step(73, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(74, "When user sets \"lineWidth\" property of line chart viewer to \"5\"", () => setProperty(page, "lineWidth", el("line chart viewer"), "5"));
      await session.step(75, "Then line chart viewer should have more ink than before", () => moreInk(page, el("line chart viewer")));
      await session.step(76, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["lineWidth","1"],["interpolation","None"]]));
      await session.step(79, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(80, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The left histogram takes a box of its own out of the plot", async () => {
      await session.step(84, "Then line chart viewer should not have a \"left panel\" area", () => hasNoArea(page, el("line chart viewer"), "left panel"));
      await session.step(85, "And the \"left panel\" reading of line chart viewer should be \"None\"", () => readingReads(page, "left panel", el("line chart viewer"), "None"));
      await session.step(86, "When user sets \"leftPanel\" property of line chart viewer to \"Histogram\"", () => setProperty(page, "leftPanel", el("line chart viewer"), "Histogram"));
      await session.step(87, "Then line chart viewer should have a \"left panel\" area", () => hasArea(page, el("line chart viewer"), "left panel"));
      await session.step(88, "And the \"left panel\" reading of line chart viewer should be \"Histogram\"", () => readingReads(page, "left panel", el("line chart viewer"), "Histogram"));
      await session.step(89, "And the \"left panel\" area of line chart viewer should be painted", () => areaPainted(page, "left panel", el("line chart viewer")));
      await session.step(90, "And the \"plot\" area of line chart viewer should be narrower than before", () => areaNarrower(page, "plot", el("line chart viewer")));
      await session.step(91, "When user sets \"leftPanel\" property of line chart viewer to \"None\"", () => setProperty(page, "leftPanel", el("line chart viewer"), "None"));
      await session.step(92, "Then line chart viewer should not have a \"left panel\" area", () => hasNoArea(page, el("line chart viewer"), "left panel"));
      await session.step(93, "And the \"left panel\" reading of line chart viewer should be \"None\"", () => readingReads(page, "left panel", el("line chart viewer"), "None"));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The overview strip appears under the plot and Overview None takes it away", async () => {
      await session.step(97, "Then line chart viewer should not have an \"overview\" area", () => hasNoArea(page, el("line chart viewer"), "overview"));
      await session.step(98, "When user picks \"Overview > Line Chart\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Overview > Line Chart", el("line chart viewer")));
      await session.step(99, "Then \"overviewType\" property of line chart viewer should be \"Line Chart\"", () => propertyShouldBe(page, "overviewType", el("line chart viewer"), "Line Chart"));
      await session.step(100, "And line chart viewer should have an \"overview\" area", () => hasArea(page, el("line chart viewer"), "overview"));
      await session.step(101, "And the \"overview\" area of line chart viewer should be painted", () => areaPainted(page, "overview", el("line chart viewer")));
      await session.step(102, "And the \"plot\" area of line chart viewer should be shorter than before", () => areaShorter(page, "plot", el("line chart viewer")));
      await session.step(103, "When user picks \"Overview > Area Chart\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Overview > Area Chart", el("line chart viewer")));
      await session.step(104, "Then \"overviewType\" property of line chart viewer should be \"Area Chart\"", () => propertyShouldBe(page, "overviewType", el("line chart viewer"), "Area Chart"));
      await session.step(105, "And the \"overview\" area of line chart viewer should have repainted", () => areaRepainted(page, "overview", el("line chart viewer")));
      await session.step(106, "When user picks \"Overview > None\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Overview > None", el("line chart viewer")));
      await session.step(107, "Then \"overviewType\" property of line chart viewer should be \"None\"", () => propertyShouldBe(page, "overviewType", el("line chart viewer"), "None"));
      await session.step(108, "And line chart viewer should not have an \"overview\" area", () => hasNoArea(page, el("line chart viewer"), "overview"));
      await session.step(109, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show X Axis and Show Y Axis take the axis boxes away", async () => {
      await session.step(112, "Then line chart viewer should have an \"x axis\" area", () => hasArea(page, el("line chart viewer"), "x axis"));
      await session.step(113, "And line chart viewer should have a \"y axis\" area", () => hasArea(page, el("line chart viewer"), "y axis"));
      await session.step(114, "And the \"y axes\" reading of line chart viewer should be 1", () => readingIs(page, "y axes", el("line chart viewer"), 1));
      await session.step(115, "When user sets \"showXAxis\" property of line chart viewer to \"false\"", () => setProperty(page, "showXAxis", el("line chart viewer"), "false"));
      await session.step(116, "Then line chart viewer should not have an \"x axis\" area", () => hasNoArea(page, el("line chart viewer"), "x axis"));
      await session.step(117, "And line chart viewer should have a \"y axis\" area", () => hasArea(page, el("line chart viewer"), "y axis"));
      await session.step(118, "When user sets \"showYAxis\" property of line chart viewer to \"false\"", () => setProperty(page, "showYAxis", el("line chart viewer"), "false"));
      await session.step(119, "Then line chart viewer should not have a \"y axis\" area", () => hasNoArea(page, el("line chart viewer"), "y axis"));
      await session.step(120, "And the \"y axes\" reading of line chart viewer should be 0", () => readingIs(page, "y axes", el("line chart viewer"), 0));
      await session.step(121, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["showXAxis","true"],["showYAxis","true"]]));
      await session.step(124, "Then line chart viewer should have an \"x axis\" area", () => hasArea(page, el("line chart viewer"), "x axis"));
      await session.step(125, "And the \"y axes\" reading of line chart viewer should be 1", () => readingIs(page, "y axes", el("line chart viewer"), 1));
      await session.step(126, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The four selector flags are read as the auto-layout resolved them", async () => {
      await session.step(129, "Then the \"x selector shown\" reading of line chart viewer should be \"true\"", () => readingReads(page, "x selector shown", el("line chart viewer"), "true"));
      await session.step(130, "And the \"y selectors shown\" reading of line chart viewer should be \"true\"", () => readingReads(page, "y selectors shown", el("line chart viewer"), "true"));
      await session.step(131, "And the \"split selector shown\" reading of line chart viewer should be \"true\"", () => readingReads(page, "split selector shown", el("line chart viewer"), "true"));
      await session.step(132, "And the \"aggr selector shown\" reading of line chart viewer should be \"false\"", () => readingReads(page, "aggr selector shown", el("line chart viewer"), "false"));
      await session.step(133, "When user sets \"showXSelector\" property of line chart viewer to \"false\"", () => setProperty(page, "showXSelector", el("line chart viewer"), "false"));
      await session.step(134, "Then the \"x selector shown\" reading of line chart viewer should be \"false\"", () => readingReads(page, "x selector shown", el("line chart viewer"), "false"));
      await session.step(135, "And the \"y selectors shown\" reading of line chart viewer should be \"true\"", () => readingReads(page, "y selectors shown", el("line chart viewer"), "true"));
      await session.step(136, "When user sets \"showYSelectors\" property of line chart viewer to \"false\"", () => setProperty(page, "showYSelectors", el("line chart viewer"), "false"));
      await session.step(137, "Then the \"y selectors shown\" reading of line chart viewer should be \"false\"", () => readingReads(page, "y selectors shown", el("line chart viewer"), "false"));
      await session.step(138, "When user sets \"showSplitSelector\" property of line chart viewer to \"false\"", () => setProperty(page, "showSplitSelector", el("line chart viewer"), "false"));
      await session.step(139, "Then the \"split selector shown\" reading of line chart viewer should be \"false\"", () => readingReads(page, "split selector shown", el("line chart viewer"), "false"));
      await session.step(140, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["showXSelector","true"],["showYSelectors","true"],["showSplitSelector","true"]]));
      await session.step(144, "Then the \"x selector shown\" reading of line chart viewer should be \"true\"", () => readingReads(page, "x selector shown", el("line chart viewer"), "true"));
      await session.step(145, "And the \"split selector shown\" reading of line chart viewer should be \"true\"", () => readingReads(page, "split selector shown", el("line chart viewer"), "true"));
      await session.step(146, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Asking for the aggregation selector is not enough to get it", async () => {
      await session.step(149, "Then \"showAggrTypeSelector\" property of line chart viewer should be \"true\"", () => propertyShouldBe(page, "showAggrTypeSelector", el("line chart viewer"), "true"));
      await session.step(150, "And the \"aggr selector shown\" reading of line chart viewer should be \"false\"", () => readingReads(page, "aggr selector shown", el("line chart viewer"), "false"));
      await session.step(151, "And the \"aggregated\" reading of line chart viewer should be \"false\"", () => readingReads(page, "aggregated", el("line chart viewer"), "false"));
      await session.step(152, "When user sets \"showAggrTypeSelector\" property of line chart viewer to \"false\"", () => setProperty(page, "showAggrTypeSelector", el("line chart viewer"), "false"));
      await session.step(153, "Then the \"aggr selector shown\" reading of line chart viewer should be \"false\"", () => readingReads(page, "aggr selector shown", el("line chart viewer"), "false"));
      await session.step(154, "When user sets \"showAggrTypeSelector\" property of line chart viewer to \"true\"", () => setProperty(page, "showAggrTypeSelector", el("line chart viewer"), "true"));
      await session.step(155, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The description sits above the plot and Description Position moves it below", async () => {
      await session.step(158, "When user sets \"description\" property of line chart viewer to \"Chemical space over the idea ids\"", () => setProperty(page, "description", el("line chart viewer"), "Chemical space over the idea ids"));
      await session.step(159, "Then the description of line chart viewer should be above its content", () => descriptionAbove(page, el("line chart viewer")));
      await session.step(160, "When user sets \"descriptionPosition\" property of line chart viewer to \"Bottom\"", () => setProperty(page, "descriptionPosition", el("line chart viewer"), "Bottom"));
      await session.step(161, "Then the description of line chart viewer should be below its content", () => descriptionBelow(page, el("line chart viewer")));
      await session.step(162, "When user sets \"descriptionVisibilityMode\" property of line chart viewer to \"Never\"", () => setProperty(page, "descriptionVisibilityMode", el("line chart viewer"), "Never"));
      await session.step(163, "Then the description of line chart viewer should be hidden", () => shouldBe(page, el("the description of line chart viewer"), "hidden"));
      await session.step(164, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["descriptionVisibilityMode","Always"],["descriptionPosition","Top"],["description",""]]));
      await session.step(168, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The chart-area context menu offers the groups the chart is configured through", async () => {
      await session.step(171, "When user opens the context menu of line chart viewer", () => openContextMenu(page, el("line chart viewer")));
      await session.step(172, "Then the open menu should list \"Reset View\"", () => menuLists(page, "Reset View"));
      await session.step(173, "And the open menu should list \"Tools\"", () => menuLists(page, "Tools"));
      await session.step(174, "And the open menu should list \"Data\"", () => menuLists(page, "Data"));
      await session.step(175, "And the open menu should list \"Markers\"", () => menuLists(page, "Markers"));
      await session.step(176, "And the open menu should list \"Chart Type\"", () => menuLists(page, "Chart Type"));
      await session.step(177, "And the open menu should list \"Overview\"", () => menuLists(page, "Overview"));
      await session.step(178, "And the open menu should list \"Selection\"", () => menuLists(page, "Selection"));
      await session.step(179, "And the open menu should list \"Controls\"", () => menuLists(page, "Controls"));
      await session.step(180, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(181, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
