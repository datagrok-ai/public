/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-axes-and-filter.feature
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
import {clickOn, hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addRangeFilter, filterBetween, filterPasses, filterToAnyOf, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaLessInk, areaMoreInk, hoverArea, noBalloons, noErrors, painted, propertyShouldBe, readingAtLeast, readingBetween, readingIs, readingLower, readingReads, repainted, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart axes, date bucketing and following the filter", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-axes-and-filter.feature", import.meta.url);
  test("Line chart axes, date bucketing and following the filter", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(21, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]));
    await session.step(24, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(25, "And \"axesFollowFilter\" property of line chart viewer should be \"true\"", () => propertyShouldBe(page, "axesFollowFilter", el("line chart viewer"), "true"));
    await session.step(26, "And the \"x axis min\" reading of line chart viewer should be between 634780 and 634782", () => readingBetween(page, "x axis min", el("line chart viewer"), 634780, 634782));
    await session.step(27, "And the \"x axis max\" reading of line chart viewer should be between 634886 and 634888", () => readingBetween(page, "x axis max", el("line chart viewer"), 634886, 634888));
    await session.step(28, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
    await session.step(29, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(30, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("Axes Follow Filter pulls the X axis onto the filtered rows", async () => {
      await session.step(33, "When user adds a range filter on \"CAST Idea ID\" from 634800 to 634850", () => addRangeFilter(page, "CAST Idea ID", 634800, 634850));
      await session.step(34, "Then 49 rows should pass the filter", () => filterPasses(page, 49));
      await session.step(35, "And the \"rows shown\" reading of line chart viewer should be 49", () => readingIs(page, "rows shown", el("line chart viewer"), 49));
      await session.step(36, "And the \"markers drawn\" reading of line chart viewer should be 49", () => readingIs(page, "markers drawn", el("line chart viewer"), 49));
      await session.step(37, "And the \"x axis span\" reading of line chart viewer should be 52", () => readingIs(page, "x axis span", el("line chart viewer"), 52));
      await session.step(38, "And the \"x axis min\" reading of line chart viewer should be between 634798 and 634800", () => readingBetween(page, "x axis min", el("line chart viewer"), 634798, 634800));
      await session.step(39, "And the \"x axis max\" reading of line chart viewer should be between 634850 and 634852", () => readingBetween(page, "x axis max", el("line chart viewer"), 634850, 634852));
      await session.step(40, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(41, "When user hovers over \"CAST Idea ID\" filter card", () => hoverOver(page, el("\"CAST Idea ID\" filter card")));
      await session.step(42, "And user clicks on close of \"CAST Idea ID\" filter card", () => clickOn(page, el("close of \"CAST Idea ID\" filter card")));
      await session.step(43, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(44, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Axes Follow Filter off the same filter leaves the axis where it was", async () => {
      await session.step(48, "When user sets \"axesFollowFilter\" property of line chart viewer to \"false\"", () => setProperty(page, "axesFollowFilter", el("line chart viewer"), "false"));
      await session.step(49, "And user adds a range filter on \"CAST Idea ID\" from 634800 to 634850", () => addRangeFilter(page, "CAST Idea ID", 634800, 634850));
      await session.step(50, "Then 49 rows should pass the filter", () => filterPasses(page, 49));
      await session.step(51, "And the \"rows shown\" reading of line chart viewer should be 49", () => readingIs(page, "rows shown", el("line chart viewer"), 49));
      await session.step(52, "And the \"markers drawn\" reading of line chart viewer should be 49", () => readingIs(page, "markers drawn", el("line chart viewer"), 49));
      await session.step(53, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(54, "When user sets \"axesFollowFilter\" property of line chart viewer to \"true\"", () => setProperty(page, "axesFollowFilter", el("line chart viewer"), "true"));
      await session.step(55, "Then the \"x axis span\" reading of line chart viewer should be 52", () => readingIs(page, "x axis span", el("line chart viewer"), 52));
      await session.step(56, "When user hovers over \"CAST Idea ID\" filter card", () => hoverOver(page, el("\"CAST Idea ID\" filter card")));
      await session.step(57, "And user clicks on close of \"CAST Idea ID\" filter card", () => clickOn(page, el("close of \"CAST Idea ID\" filter card")));
      await session.step(58, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(59, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic X axis over a column that starts below zero drops the non-positive points", async () => {
      await session.step(63, "When user sets \"xColumnName\" property of line chart viewer to \"Chemical Space X\"", () => setProperty(page, "xColumnName", el("line chart viewer"), "Chemical Space X"));
      await session.step(64, "Then the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(65, "And the \"x axis min\" reading of line chart viewer should be between -4.7 and -4.6", () => readingBetween(page, "x axis min", el("line chart viewer"), -4.7, -4.6));
      await session.step(66, "When user sets \"xAxisType\" property of line chart viewer to \"logarithmic\"", () => setProperty(page, "xAxisType", el("line chart viewer"), "logarithmic"));
      await session.step(67, "Then the \"markers drawn\" reading of line chart viewer should be 59", () => readingIs(page, "markers drawn", el("line chart viewer"), 59));
      await session.step(68, "And the \"x axis min\" reading of line chart viewer should be between 0 and 1", () => readingBetween(page, "x axis min", el("line chart viewer"), 0, 1));
      await session.step(69, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(70, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(71, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xAxisType","linear"],["xColumnName","CAST Idea ID"]]));
      await session.step(74, "Then the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hovering an empty logarithmic chart raises nothing (github-2574)", async () => {
      await session.step(78, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xColumnName","Chemical Space X"],["xAxisType","logarithmic"]]));
      await session.step(81, "And user filters rows where \"TPSA\" is between 0 and 1", () => filterBetween(page, "TPSA", 0, 1));
      await session.step(82, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(83, "And the \"rows shown\" reading of line chart viewer should be 0", () => readingIs(page, "rows shown", el("line chart viewer"), 0));
      await session.step(84, "And the \"lines\" reading of line chart viewer should be 0", () => readingIs(page, "lines", el("line chart viewer"), 0));
      await session.step(85, "And the \"markers drawn\" reading of line chart viewer should be 0", () => readingIs(page, "markers drawn", el("line chart viewer"), 0));
      await session.step(86, "When user hovers over the \"plot\" area of line chart viewer", () => hoverArea(page, "plot", el("line chart viewer")));
      await session.step(87, "Then no errors should have been logged", () => noErrors(page));
      await session.step(88, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(89, "When user resets the filter", () => resetFilter(page));
      await session.step(90, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(91, "And the \"markers drawn\" reading of line chart viewer should be 59", () => readingIs(page, "markers drawn", el("line chart viewer"), 59));
      await session.step(92, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xAxisType","linear"],["xColumnName","CAST Idea ID"]]));
      await session.step(95, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A date column buckets into as many X positions as the mapping asks for", async () => {
      await session.step(98, "When user sets \"xColumnName\" property of line chart viewer to \"Competition assay Date\"", () => setProperty(page, "xColumnName", el("line chart viewer"), "Competition assay Date"));
      await session.step(99, "Then the \"aggregated\" reading of line chart viewer should be \"true\"", () => readingReads(page, "aggregated", el("line chart viewer"), "true"));
      await session.step(100, "And the \"x categories\" reading of line chart viewer should be 50", () => readingIs(page, "x categories", el("line chart viewer"), 50));
      await session.step(101, "And the \"x column\" reading of line chart viewer should be \"Competition assay Date\"", () => readingReads(page, "x column", el("line chart viewer"), "Competition assay Date"));
      await session.step(102, "When user sets \"xMap\" property of line chart viewer to \"year\"", () => setProperty(page, "xMap", el("line chart viewer"), "year"));
      await session.step(103, "Then the \"x column\" reading of line chart viewer should be \"Competition assay Date year\"", () => readingReads(page, "x column", el("line chart viewer"), "Competition assay Date year"));
      await session.step(104, "And the \"x categories\" reading of line chart viewer should be 3", () => readingIs(page, "x categories", el("line chart viewer"), 3));
      await session.step(105, "And the \"markers drawn\" reading of line chart viewer should be 3", () => readingIs(page, "markers drawn", el("line chart viewer"), 3));
      await session.step(106, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(107, "When user sets \"xMap\" property of line chart viewer to \"month\"", () => setProperty(page, "xMap", el("line chart viewer"), "month"));
      await session.step(108, "Then the \"x categories\" reading of line chart viewer should be 12", () => readingIs(page, "x categories", el("line chart viewer"), 12));
      await session.step(109, "When user sets \"xMap\" property of line chart viewer to \"year quarter\"", () => setProperty(page, "xMap", el("line chart viewer"), "year quarter"));
      await session.step(110, "Then the \"x categories\" reading of line chart viewer should be 11", () => readingIs(page, "x categories", el("line chart viewer"), 11));
      await session.step(111, "And the \"aggregation\" reading of line chart viewer should be \"avg\"", () => readingReads(page, "aggregation", el("line chart viewer"), "avg"));
      await session.step(112, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(113, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xMap",""],["xColumnName","CAST Idea ID"]]));
      await session.step(116, "Then the \"aggregated\" reading of line chart viewer should be \"false\"", () => readingReads(page, "aggregated", el("line chart viewer"), "false"));
      await session.step(117, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Filtering while the X axis is year-quarter buckets keeps the chart drawing (GROK-18375)", async () => {
      await session.step(120, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xColumnName","Competition assay Date"],["xMap","year quarter"]]));
      await session.step(123, "Then the \"x categories\" reading of line chart viewer should be 11", () => readingIs(page, "x categories", el("line chart viewer"), 11));
      await session.step(124, "When user filters rows where \"Series\" is one of \"Aminopiperidines, Pyrrolidines\"", () => filterToAnyOf(page, "Series", "Aminopiperidines, Pyrrolidines"));
      await session.step(125, "Then 25 rows should pass the filter", () => filterPasses(page, 25));
      await session.step(126, "And the \"rows shown\" reading of line chart viewer should be 25", () => readingIs(page, "rows shown", el("line chart viewer"), 25));
      await session.step(127, "And the \"x categories\" reading of line chart viewer should be lower than before", () => readingLower(page, "x categories", el("line chart viewer")));
      await session.step(128, "And the \"markers drawn\" reading of line chart viewer should be at least 1", () => readingAtLeast(page, "markers drawn", el("line chart viewer"), 1));
      await session.step(129, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(130, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(131, "When user resets the filter", () => resetFilter(page));
      await session.step(132, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(133, "And the \"x categories\" reading of line chart viewer should be 11", () => readingIs(page, "x categories", el("line chart viewer"), 11));
      await session.step(134, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xMap",""],["xColumnName","CAST Idea ID"]]));
      await session.step(137, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Narrowing the X column to its middle half keeps 29 rows (GROK-20185)", async () => {
      await session.step(140, "When user sets \"xColumnName\" property of line chart viewer to \"Chemical Space X\"", () => setProperty(page, "xColumnName", el("line chart viewer"), "Chemical Space X"));
      await session.step(141, "Then the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(142, "When user adds a range filter on \"Chemical Space X\" from 0.5731 to 10.2274", () => addRangeFilter(page, "Chemical Space X", 0.5731, 10.2274));
      await session.step(143, "Then 29 rows should pass the filter", () => filterPasses(page, 29));
      await session.step(144, "And the \"rows shown\" reading of line chart viewer should be 29", () => readingIs(page, "rows shown", el("line chart viewer"), 29));
      await session.step(145, "And the \"markers drawn\" reading of line chart viewer should be 29", () => readingIs(page, "markers drawn", el("line chart viewer"), 29));
      await session.step(146, "And the \"x axis span\" reading of line chart viewer should be lower than before", () => readingLower(page, "x axis span", el("line chart viewer")));
      await session.step(147, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(148, "When user hovers over \"Chemical Space X\" filter card", () => hoverOver(page, el("\"Chemical Space X\" filter card")));
      await session.step(149, "And user clicks on close of \"Chemical Space X\" filter card", () => clickOn(page, el("close of \"Chemical Space X\" filter card")));
      await session.step(150, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(151, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(152, "When user sets \"xColumnName\" property of line chart viewer to \"CAST Idea ID\"", () => setProperty(page, "xColumnName", el("line chart viewer"), "CAST Idea ID"));
      await session.step(153, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Inverting the axis and dropping the grid lines redraw without moving the range", async () => {
      await session.step(156, "When user sets \"invertXAxis\" property of line chart viewer to \"true\"", () => setProperty(page, "invertXAxis", el("line chart viewer"), "true"));
      await session.step(157, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(158, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(159, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(160, "When user sets \"showVerticalGridLines\" property of line chart viewer to \"false\"", () => setProperty(page, "showVerticalGridLines", el("line chart viewer"), "false"));
      await session.step(161, "Then the \"plot\" area of line chart viewer should have less ink than before", () => areaLessInk(page, "plot", el("line chart viewer")));
      await session.step(162, "When user sets \"showHorizontalGridLines\" property of line chart viewer to \"false\"", () => setProperty(page, "showHorizontalGridLines", el("line chart viewer"), "false"));
      await session.step(163, "Then the \"plot\" area of line chart viewer should have less ink than before", () => areaLessInk(page, "plot", el("line chart viewer")));
      await session.step(164, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["invertXAxis","false"],["showVerticalGridLines","true"],["showHorizontalGridLines","true"]]));
      await session.step(168, "Then the \"plot\" area of line chart viewer should have more ink than before", () => areaMoreInk(page, "plot", el("line chart viewer")));
      await session.step(169, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("MinMax tickmarks leave the axis with fewer labels than Auto", async () => {
      await session.step(172, "When user sets \"xAxisTickmarksMode\" property of line chart viewer to \"MinMax\"", () => setProperty(page, "xAxisTickmarksMode", el("line chart viewer"), "MinMax"));
      await session.step(173, "Then the \"x axis\" area of line chart viewer should have less ink than before", () => areaLessInk(page, "x axis", el("line chart viewer")));
      await session.step(174, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(175, "When user sets \"yAxisTickmarksMode\" property of line chart viewer to \"MinMax\"", () => setProperty(page, "yAxisTickmarksMode", el("line chart viewer"), "MinMax"));
      await session.step(176, "Then the \"y axis\" area of line chart viewer should have less ink than before", () => areaLessInk(page, "y axis", el("line chart viewer")));
      await session.step(177, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xAxisTickmarksMode","Auto"],["yAxisTickmarksMode","Auto"]]));
      await session.step(180, "Then the \"x axis\" area of line chart viewer should have more ink than before", () => areaMoreInk(page, "x axis", el("line chart viewer")));
      await session.step(181, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
