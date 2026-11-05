/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-axes-and-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart, public#2485]
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
import {clickOn, hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addRangeFilter, filterBetween, filterPasses, filterToAnyOf, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaLessInk, areaMoreInk, hoverArea, noBalloons, noErrors, painted, propertyShouldBe, readingAtLeast, readingBetween, readingIs, readingLower, readingReads, repainted, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart axes, date bucketing and following the filter", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-axes-and-filter.feature", import.meta.url);
  test("Line chart axes, date bucketing and following the filter", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart", "@realizes:public#2485"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(23, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
    await session.step(26, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(27, "And \"axesFollowFilter\" property of line chart viewer should be \"true\"", () => propertyShouldBe(page, "axesFollowFilter", el("line chart viewer"), "true"));
    await session.step(28, "And the \"x axis min\" reading of line chart viewer should be between 634780 and 634782", () => readingBetween(page, "x axis min", el("line chart viewer"), 634780, 634782));
    await session.step(29, "And the \"x axis max\" reading of line chart viewer should be between 634886 and 634888", () => readingBetween(page, "x axis max", el("line chart viewer"), 634886, 634888));
    await session.step(30, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
    await session.step(31, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(32, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("Axes Follow Filter pulls the X axis onto the filtered rows", async () => {
      await session.step(35, "When user adds a range filter on \"CAST Idea ID\" from 634800 to 634850", () => addRangeFilter(page, "CAST Idea ID", 634800, 634850));
      await session.step(36, "Then 49 rows should pass the filter", () => filterPasses(page, 49));
      await session.step(37, "And the \"rows shown\" reading of line chart viewer should be 49", () => readingIs(page, "rows shown", el("line chart viewer"), 49));
      await session.step(38, "And the \"markers drawn\" reading of line chart viewer should be 49", () => readingIs(page, "markers drawn", el("line chart viewer"), 49));
      await session.step(39, "And the \"x axis span\" reading of line chart viewer should be 52", () => readingIs(page, "x axis span", el("line chart viewer"), 52));
      await session.step(40, "And the \"x axis min\" reading of line chart viewer should be between 634798 and 634800", () => readingBetween(page, "x axis min", el("line chart viewer"), 634798, 634800));
      await session.step(41, "And the \"x axis max\" reading of line chart viewer should be between 634850 and 634852", () => readingBetween(page, "x axis max", el("line chart viewer"), 634850, 634852));
      await session.step(42, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(43, "When user hovers over \"CAST Idea ID\" filter card", () => hoverOver(page, el("\"CAST Idea ID\" filter card")));
      await session.step(44, "And user clicks on close of \"CAST Idea ID\" filter card", () => clickOn(page, el("close of \"CAST Idea ID\" filter card")));
      await session.step(45, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(46, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Axes Follow Filter off the same filter leaves the axis where it was", async () => {
      await session.step(50, "When user sets \"axesFollowFilter\" property of line chart viewer to \"false\"", () => setProperty(page, "axesFollowFilter", el("line chart viewer"), "false"));
      await session.step(51, "And user adds a range filter on \"CAST Idea ID\" from 634800 to 634850", () => addRangeFilter(page, "CAST Idea ID", 634800, 634850));
      await session.step(52, "Then 49 rows should pass the filter", () => filterPasses(page, 49));
      await session.step(53, "And the \"rows shown\" reading of line chart viewer should be 49", () => readingIs(page, "rows shown", el("line chart viewer"), 49));
      await session.step(54, "And the \"markers drawn\" reading of line chart viewer should be 49", () => readingIs(page, "markers drawn", el("line chart viewer"), 49));
      await session.step(55, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(56, "When user sets \"axesFollowFilter\" property of line chart viewer to \"true\"", () => setProperty(page, "axesFollowFilter", el("line chart viewer"), "true"));
      await session.step(57, "Then the \"x axis span\" reading of line chart viewer should be 52", () => readingIs(page, "x axis span", el("line chart viewer"), 52));
      await session.step(58, "When user hovers over \"CAST Idea ID\" filter card", () => hoverOver(page, el("\"CAST Idea ID\" filter card")));
      await session.step(59, "And user clicks on close of \"CAST Idea ID\" filter card", () => clickOn(page, el("close of \"CAST Idea ID\" filter card")));
      await session.step(60, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(61, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic X axis over a column that starts below zero drops the non-positive points", async () => {
      await session.step(65, "When user sets \"xColumnName\" property of line chart viewer to \"Chemical Space X\"", () => setProperty(page, "xColumnName", el("line chart viewer"), "Chemical Space X"));
      await session.step(66, "Then the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(67, "And the \"x axis min\" reading of line chart viewer should be between -4.7 and -4.6", () => readingBetween(page, "x axis min", el("line chart viewer"), -4.7, -4.6));
      await session.step(68, "When user sets \"xAxisType\" property of line chart viewer to \"logarithmic\"", () => setProperty(page, "xAxisType", el("line chart viewer"), "logarithmic"));
      await session.step(69, "Then the \"markers drawn\" reading of line chart viewer should be 59", () => readingIs(page, "markers drawn", el("line chart viewer"), 59));
      await session.step(70, "And the \"x axis min\" reading of line chart viewer should be between 0 and 1", () => readingBetween(page, "x axis min", el("line chart viewer"), 0, 1));
      await session.step(71, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(72, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(73, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xAxisType","linear"],["xColumnName","CAST Idea ID"]]), [["xAxisType","linear"],["xColumnName","CAST Idea ID"]]);
      await session.step(76, "Then the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hovering an empty logarithmic chart raises nothing (github-2574)", async () => {
      await session.step(80, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xColumnName","Chemical Space X"],["xAxisType","logarithmic"]]), [["xColumnName","Chemical Space X"],["xAxisType","logarithmic"]]);
      await session.step(83, "And user filters rows where \"TPSA\" is between 0 and 1", () => filterBetween(page, "TPSA", 0, 1));
      await session.step(84, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(85, "And the \"rows shown\" reading of line chart viewer should be 0", () => readingIs(page, "rows shown", el("line chart viewer"), 0));
      await session.step(86, "And the \"lines\" reading of line chart viewer should be 0", () => readingIs(page, "lines", el("line chart viewer"), 0));
      await session.step(87, "And the \"markers drawn\" reading of line chart viewer should be 0", () => readingIs(page, "markers drawn", el("line chart viewer"), 0));
      await session.step(88, "When user hovers over the \"plot\" area of line chart viewer", () => hoverArea(page, "plot", el("line chart viewer")));
      await session.step(89, "Then no errors should have been logged", () => noErrors(page));
      await session.step(90, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(91, "When user resets the filter", () => resetFilter(page));
      await session.step(92, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(93, "And the \"markers drawn\" reading of line chart viewer should be 59", () => readingIs(page, "markers drawn", el("line chart viewer"), 59));
      await session.step(94, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xAxisType","linear"],["xColumnName","CAST Idea ID"]]), [["xAxisType","linear"],["xColumnName","CAST Idea ID"]]);
      await session.step(97, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A date column buckets into as many X positions as the mapping asks for", async () => {
      await session.step(100, "When user sets \"xColumnName\" property of line chart viewer to \"Competition assay Date\"", () => setProperty(page, "xColumnName", el("line chart viewer"), "Competition assay Date"));
      await session.step(101, "Then the \"aggregated\" reading of line chart viewer should be \"true\"", () => readingReads(page, "aggregated", el("line chart viewer"), "true"));
      await session.step(102, "And the \"x categories\" reading of line chart viewer should be 50", () => readingIs(page, "x categories", el("line chart viewer"), 50));
      await session.step(103, "And the \"x column\" reading of line chart viewer should be \"Competition assay Date\"", () => readingReads(page, "x column", el("line chart viewer"), "Competition assay Date"));
      await session.step(104, "When user sets \"xMap\" property of line chart viewer to \"year\"", () => setProperty(page, "xMap", el("line chart viewer"), "year"));
      await session.step(105, "Then the \"x column\" reading of line chart viewer should be \"Competition assay Date year\"", () => readingReads(page, "x column", el("line chart viewer"), "Competition assay Date year"));
      await session.step(106, "And the \"x categories\" reading of line chart viewer should be 3", () => readingIs(page, "x categories", el("line chart viewer"), 3));
      await session.step(107, "And the \"markers drawn\" reading of line chart viewer should be 3", () => readingIs(page, "markers drawn", el("line chart viewer"), 3));
      await session.step(108, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(109, "When user sets \"xMap\" property of line chart viewer to \"month\"", () => setProperty(page, "xMap", el("line chart viewer"), "month"));
      await session.step(110, "Then the \"x categories\" reading of line chart viewer should be 12", () => readingIs(page, "x categories", el("line chart viewer"), 12));
      await session.step(111, "When user sets \"xMap\" property of line chart viewer to \"year quarter\"", () => setProperty(page, "xMap", el("line chart viewer"), "year quarter"));
      await session.step(112, "Then the \"x categories\" reading of line chart viewer should be 11", () => readingIs(page, "x categories", el("line chart viewer"), 11));
      await session.step(113, "And the \"aggregation\" reading of line chart viewer should be \"avg\"", () => readingReads(page, "aggregation", el("line chart viewer"), "avg"));
      await session.step(114, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(115, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xMap",""],["xColumnName","CAST Idea ID"]]), [["xMap",""],["xColumnName","CAST Idea ID"]]);
      await session.step(118, "Then the \"aggregated\" reading of line chart viewer should be \"false\"", () => readingReads(page, "aggregated", el("line chart viewer"), "false"));
      await session.step(119, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Filtering while the X axis is year-quarter buckets keeps the chart drawing (GROK-18375)", async () => {
      await session.step(122, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xColumnName","Competition assay Date"],["xMap","year quarter"]]), [["xColumnName","Competition assay Date"],["xMap","year quarter"]]);
      await session.step(125, "Then the \"x categories\" reading of line chart viewer should be 11", () => readingIs(page, "x categories", el("line chart viewer"), 11));
      await session.step(126, "When user filters rows where \"Series\" is one of \"Aminopiperidines, Pyrrolidines\"", () => filterToAnyOf(page, "Series", "Aminopiperidines, Pyrrolidines"));
      await session.step(127, "Then 25 rows should pass the filter", () => filterPasses(page, 25));
      await session.step(128, "And the \"rows shown\" reading of line chart viewer should be 25", () => readingIs(page, "rows shown", el("line chart viewer"), 25));
      await session.step(129, "And the \"x categories\" reading of line chart viewer should be lower than before", () => readingLower(page, "x categories", el("line chart viewer")));
      await session.step(130, "And the \"markers drawn\" reading of line chart viewer should be at least 1", () => readingAtLeast(page, "markers drawn", el("line chart viewer"), 1));
      await session.step(131, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(132, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(133, "When user resets the filter", () => resetFilter(page));
      await session.step(134, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(135, "And the \"x categories\" reading of line chart viewer should be 11", () => readingIs(page, "x categories", el("line chart viewer"), 11));
      await session.step(136, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xMap",""],["xColumnName","CAST Idea ID"]]), [["xMap",""],["xColumnName","CAST Idea ID"]]);
      await session.step(139, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Narrowing the X column to its middle half keeps 29 rows (GROK-20185)", async () => {
      await session.step(142, "When user sets \"xColumnName\" property of line chart viewer to \"Chemical Space X\"", () => setProperty(page, "xColumnName", el("line chart viewer"), "Chemical Space X"));
      await session.step(143, "Then the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(144, "When user adds a range filter on \"Chemical Space X\" from 0.5731 to 10.2274", () => addRangeFilter(page, "Chemical Space X", 0.5731, 10.2274));
      await session.step(145, "Then 29 rows should pass the filter", () => filterPasses(page, 29));
      await session.step(146, "And the \"rows shown\" reading of line chart viewer should be 29", () => readingIs(page, "rows shown", el("line chart viewer"), 29));
      await session.step(147, "And the \"markers drawn\" reading of line chart viewer should be 29", () => readingIs(page, "markers drawn", el("line chart viewer"), 29));
      await session.step(148, "And the \"x axis span\" reading of line chart viewer should be lower than before", () => readingLower(page, "x axis span", el("line chart viewer")));
      await session.step(149, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(150, "When user hovers over \"Chemical Space X\" filter card", () => hoverOver(page, el("\"Chemical Space X\" filter card")));
      await session.step(151, "And user clicks on close of \"Chemical Space X\" filter card", () => clickOn(page, el("close of \"Chemical Space X\" filter card")));
      await session.step(152, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(153, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(154, "When user sets \"xColumnName\" property of line chart viewer to \"CAST Idea ID\"", () => setProperty(page, "xColumnName", el("line chart viewer"), "CAST Idea ID"));
      await session.step(155, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Inverting the axis and dropping the grid lines redraw without moving the range", async () => {
      await session.step(158, "When user sets \"invertXAxis\" property of line chart viewer to \"true\"", () => setProperty(page, "invertXAxis", el("line chart viewer"), "true"));
      await session.step(159, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(160, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(161, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(162, "When user sets \"showVerticalGridLines\" property of line chart viewer to \"false\"", () => setProperty(page, "showVerticalGridLines", el("line chart viewer"), "false"));
      await session.step(163, "Then the \"plot\" area of line chart viewer should have less ink than before", () => areaLessInk(page, "plot", el("line chart viewer")));
      await session.step(164, "When user sets \"showHorizontalGridLines\" property of line chart viewer to \"false\"", () => setProperty(page, "showHorizontalGridLines", el("line chart viewer"), "false"));
      await session.step(165, "Then the \"plot\" area of line chart viewer should have less ink than before", () => areaLessInk(page, "plot", el("line chart viewer")));
      await session.step(166, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["invertXAxis","false"],["showVerticalGridLines","true"],["showHorizontalGridLines","true"]]), [["invertXAxis","false"],["showVerticalGridLines","true"],["showHorizontalGridLines","true"]]);
      await session.step(170, "Then the \"plot\" area of line chart viewer should have more ink than before", () => areaMoreInk(page, "plot", el("line chart viewer")));
      await session.step(171, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("MinMax tickmarks leave the axis with fewer labels than Auto", async () => {
      await session.step(174, "When user sets \"xAxisTickmarksMode\" property of line chart viewer to \"MinMax\"", () => setProperty(page, "xAxisTickmarksMode", el("line chart viewer"), "MinMax"));
      await session.step(175, "Then the \"x axis\" area of line chart viewer should have less ink than before", () => areaLessInk(page, "x axis", el("line chart viewer")));
      await session.step(176, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(177, "When user sets \"yAxisTickmarksMode\" property of line chart viewer to \"MinMax\"", () => setProperty(page, "yAxisTickmarksMode", el("line chart viewer"), "MinMax"));
      await session.step(178, "Then the \"y axis\" area of line chart viewer should have less ink than before", () => areaLessInk(page, "y axis", el("line chart viewer")));
      await session.step(179, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xAxisTickmarksMode","Auto"],["yAxisTickmarksMode","Auto"]]), [["xAxisTickmarksMode","Auto"],["yAxisTickmarksMode","Auto"]]);
      await session.step(182, "Then the \"x axis\" area of line chart viewer should have more ink than before", () => areaMoreInk(page, "x axis", el("line chart viewer")));
      await session.step(183, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic Y axis fits the data instead of leaving an empty band above it (public#2485)", async () => {
      await session.step(187, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xColumnName","Num Heavy Atoms"],["yColumnNames","Whole blood assay 1"],["yAxisType","logarithmic"]]), [["xColumnName","Num Heavy Atoms"],["yColumnNames","Whole blood assay 1"],["yAxisType","logarithmic"]]);
      await session.step(191, "Then the \"y axis max of \\\"Whole blood assay 1\\\"\" reading of line chart viewer should be between 31.5 and 37.8", () => readingBetween(page, "y axis max of \"Whole blood assay 1\"", el("line chart viewer"), 31.5, 37.8));
      await session.step(192, "And the \"y axis min of \\\"Whole blood assay 1\\\"\" reading of line chart viewer should be between 0.04 and 0.14", () => readingBetween(page, "y axis min of \"Whole blood assay 1\"", el("line chart viewer"), 0.04, 0.14));
      await session.step(193, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(194, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yAxisType","linear"],["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["yAxisType","linear"],["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
      await session.step(198, "Then the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(199, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
