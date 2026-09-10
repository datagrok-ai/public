/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-aggregation-and-markers.feature
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
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hasNoArea, lessInk, moreInk, noErrors, propertyShouldBe, readingHigher, readingIs, readingLower, readingReads, repainted, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart aggregation, whiskers and markers", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-aggregation-and-markers.feature", import.meta.url);
  test("Line chart aggregation, whiskers and markers", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(18, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","Competition assay Date"],["xMap","year quarter"],["yColumnNames","Chemical Space X"]]));
    await session.step(22, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(23, "And the \"aggregated\" reading of line chart viewer should be \"true\"", () => readingReads(page, "aggregated", el("line chart viewer"), "true"));
    await session.step(24, "And the \"x categories\" reading of line chart viewer should be 11", () => readingIs(page, "x categories", el("line chart viewer"), 11));
    await session.step(25, "And the \"aggregation\" reading of line chart viewer should be \"avg\"", () => readingReads(page, "aggregation", el("line chart viewer"), "avg"));
    await session.step(26, "And the \"markers drawn\" reading of line chart viewer should be 11", () => readingIs(page, "markers drawn", el("line chart viewer"), 11));
    await session.step(27, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("The aggregation function is named, and changing it moves the points", async () => {
      await session.step(30, "When user sets \"aggrType\" property of line chart viewer to \"max\"", () => setProperty(page, "aggrType", el("line chart viewer"), "max"));
      await session.step(31, "Then the \"aggregation\" reading of line chart viewer should be \"max\"", () => readingReads(page, "aggregation", el("line chart viewer"), "max"));
      await session.step(32, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be higher than before", () => readingHigher(page, "y axis max of \"Chemical Space X\"", el("line chart viewer")));
      await session.step(33, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(34, "When user sets \"aggrType\" property of line chart viewer to \"min\"", () => setProperty(page, "aggrType", el("line chart viewer"), "min"));
      await session.step(35, "Then the \"aggregation\" reading of line chart viewer should be \"min\"", () => readingReads(page, "aggregation", el("line chart viewer"), "min"));
      await session.step(36, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be lower than before", () => readingLower(page, "y axis max of \"Chemical Space X\"", el("line chart viewer")));
      await session.step(37, "When user sets \"aggrType\" property of line chart viewer to \"avg\"", () => setProperty(page, "aggrType", el("line chart viewer"), "avg"));
      await session.step(38, "Then the \"aggregation\" reading of line chart viewer should be \"avg\"", () => readingReads(page, "aggregation", el("line chart viewer"), "avg"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Whiskers draw the spread the aggregation collapsed", async () => {
      await session.step(42, "Then line chart viewer should not have a \"whiskers\" area", () => hasNoArea(page, el("line chart viewer"), "whiskers"));
      await session.step(43, "And the \"whiskers\" reading of line chart viewer should be \"None\"", () => readingReads(page, "whiskers", el("line chart viewer"), "None"));
      await session.step(44, "When user sets \"whiskersType\" property of line chart viewer to \"Avg | ±StError\"", () => setProperty(page, "whiskersType", el("line chart viewer"), "Avg | ±StError"));
      await session.step(45, "Then the \"whiskers\" reading of line chart viewer should be \"Avg | ±StError\"", () => readingReads(page, "whiskers", el("line chart viewer"), "Avg | ±StError"));
      await session.step(46, "And line chart viewer should have a \"whiskers\" area", () => hasArea(page, el("line chart viewer"), "whiskers"));
      await session.step(47, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be higher than before", () => readingHigher(page, "y axis max of \"Chemical Space X\"", el("line chart viewer")));
      await session.step(48, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(49, "When user sets \"whiskersType\" property of line chart viewer to \"None\"", () => setProperty(page, "whiskersType", el("line chart viewer"), "None"));
      await session.step(50, "Then line chart viewer should not have a \"whiskers\" area", () => hasNoArea(page, el("line chart viewer"), "whiskers"));
      await session.step(51, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be lower than before", () => readingLower(page, "y axis max of \"Chemical Space X\"", el("line chart viewer")));
      await session.step(52, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Whiskers on a chart that is not aggregated are not drawn whatever the property says", async () => {
      await session.step(56, "When user sets \"whiskersType\" property of line chart viewer to \"Avg | ±StError\"", () => setProperty(page, "whiskersType", el("line chart viewer"), "Avg | ±StError"));
      await session.step(57, "Then line chart viewer should have a \"whiskers\" area", () => hasArea(page, el("line chart viewer"), "whiskers"));
      await session.step(58, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xMap",""],["xColumnName","CAST Idea ID"]]));
      await session.step(61, "Then the \"aggregated\" reading of line chart viewer should be \"false\"", () => readingReads(page, "aggregated", el("line chart viewer"), "false"));
      await session.step(62, "And \"whiskersType\" property of line chart viewer should be \"Avg | ±StError\"", () => propertyShouldBe(page, "whiskersType", el("line chart viewer"), "Avg | ±StError"));
      await session.step(63, "And the \"whiskers\" reading of line chart viewer should be \"None\"", () => readingReads(page, "whiskers", el("line chart viewer"), "None"));
      await session.step(64, "And line chart viewer should not have a \"whiskers\" area", () => hasNoArea(page, el("line chart viewer"), "whiskers"));
      await session.step(65, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["whiskersType",""],["xColumnName","Competition assay Date"],["xMap","year quarter"]]));
      await session.step(69, "Then the \"aggregated\" reading of line chart viewer should be \"true\"", () => readingReads(page, "aggregated", el("line chart viewer"), "true"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Markers are drawn unless the chart is told never to draw them", async () => {
      await session.step(73, "Then the \"markers drawn\" reading of line chart viewer should be 11", () => readingIs(page, "markers drawn", el("line chart viewer"), 11));
      await session.step(74, "When user sets \"showMarkers\" property of line chart viewer to \"Never\"", () => setProperty(page, "showMarkers", el("line chart viewer"), "Never"));
      await session.step(75, "Then the \"markers drawn\" reading of line chart viewer should be 0", () => readingIs(page, "markers drawn", el("line chart viewer"), 0));
      await session.step(76, "And line chart viewer should have less ink than before", () => lessInk(page, el("line chart viewer")));
      await session.step(77, "When user sets \"showMarkers\" property of line chart viewer to \"Always\"", () => setProperty(page, "showMarkers", el("line chart viewer"), "Always"));
      await session.step(78, "Then the \"markers drawn\" reading of line chart viewer should be 11", () => readingIs(page, "markers drawn", el("line chart viewer"), 11));
      await session.step(79, "And line chart viewer should have more ink than before", () => moreInk(page, el("line chart viewer")));
      await session.step(80, "When user sets \"showMarkers\" property of line chart viewer to \"Auto\"", () => setProperty(page, "showMarkers", el("line chart viewer"), "Auto"));
      await session.step(81, "Then the \"markers drawn\" reading of line chart viewer should be 11", () => readingIs(page, "markers drawn", el("line chart viewer"), 11));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A size column is reported and a marker type change redraws", async () => {
      await session.step(85, "Then the \"marker size column\" reading of line chart viewer should be \"\"", () => readingReads(page, "marker size column", el("line chart viewer"), ""));
      await session.step(86, "When user sets \"markerType\" property of line chart viewer to \"square\"", () => setProperty(page, "markerType", el("line chart viewer"), "square"));
      await session.step(87, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(88, "When user sets \"markersSizeColumnName\" property of line chart viewer to \"Chemical Space Y\"", () => setProperty(page, "markersSizeColumnName", el("line chart viewer"), "Chemical Space Y"));
      await session.step(89, "Then the \"marker size column\" reading of line chart viewer should be \"Chemical Space Y\"", () => readingReads(page, "marker size column", el("line chart viewer"), "Chemical Space Y"));
      await session.step(90, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(91, "And the \"markers drawn\" reading of line chart viewer should be 11", () => readingIs(page, "markers drawn", el("line chart viewer"), 11));
      await session.step(92, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["markersSizeColumnName",""],["markerType","circle"]]));
      await session.step(95, "Then the \"marker size column\" reading of line chart viewer should be \"\"", () => readingReads(page, "marker size column", el("line chart viewer"), ""));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Splitting an aggregated chart aggregates within every category", async () => {
      await session.step(99, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(100, "Then the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(101, "And the \"categories\" reading of line chart viewer should be 5", () => readingIs(page, "categories", el("line chart viewer"), 5));
      await session.step(102, "And the \"markers drawn\" reading of line chart viewer should be 32", () => readingIs(page, "markers drawn", el("line chart viewer"), 32));
      await session.step(103, "And the \"x categories\" reading of line chart viewer should be 11", () => readingIs(page, "x categories", el("line chart viewer"), 11));
      await session.step(104, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(105, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(106, "Then the \"markers drawn\" reading of line chart viewer should be 11", () => readingIs(page, "markers drawn", el("line chart viewer"), 11));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The aggregation selector shows up only for a multi-axis aggregated split", async () => {
      await session.step(110, "Then the \"aggr selector shown\" reading of line chart viewer should be \"false\"", () => readingReads(page, "aggr selector shown", el("line chart viewer"), "false"));
      await session.step(111, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(112, "Then the \"aggr selector shown\" reading of line chart viewer should be \"false\"", () => readingReads(page, "aggr selector shown", el("line chart viewer"), "false"));
      await session.step(113, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","Chemical Space X, TPSA"],["multiAxis","true"]]));
      await session.step(116, "Then the \"aggr selector shown\" reading of line chart viewer should be \"true\"", () => readingReads(page, "aggr selector shown", el("line chart viewer"), "true"));
      await session.step(117, "And the \"aggregation\" reading of line chart viewer should be \"avg, avg\"", () => readingReads(page, "aggregation", el("line chart viewer"), "avg, avg"));
      await session.step(118, "And the \"lines\" reading of line chart viewer should be 10", () => readingIs(page, "lines", el("line chart viewer"), 10));
      await session.step(119, "When user sets \"multiAxis\" property of line chart viewer to \"false\"", () => setProperty(page, "multiAxis", el("line chart viewer"), "false"));
      await session.step(120, "Then the \"aggr selector shown\" reading of line chart viewer should be \"false\"", () => readingReads(page, "aggr selector shown", el("line chart viewer"), "false"));
      await session.step(121, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","Chemical Space X"],["splitColumnNames",""]]));
      await session.step(124, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
