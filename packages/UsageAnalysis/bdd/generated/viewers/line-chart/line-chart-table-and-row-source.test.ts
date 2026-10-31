/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-table-and-row-source.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart, public#1650, public#1628]
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
import {clearSelection, filterBetween, filterPasses, resetFilter, selectWhereIs, selectWhereOneOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, boundTable, noErrors, painted, propertyShouldBe, readingAtLeast, readingDoesNotRead, readingHigher, readingIs, readingLower, readingReads, repainted, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart table binding, Row Source and the viewer's own filter", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-table-and-row-source.feature", import.meta.url);
  test("Line chart table binding, Row Source and the viewer's own filter", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart", "@realizes:public#1650", "@realizes:public#1628"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(27, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
    await session.step(30, "Then line chart viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("line chart viewer"), "spgi-100"));
    await session.step(31, "And the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
    await session.step(32, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(33, "And \"rowSource\" property of line chart viewer should be \"Filtered\"", () => propertyShouldBe(page, "rowSource", el("line chart viewer"), "Filtered"));
    await session.step(34, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("Row Source Selected draws the selected rows and nothing else", async () => {
      await session.step(37, "When user selects rows where \"Stereo Category\" is \"R_ONE\"", () => selectWhereIs(page, "Stereo Category", "R_ONE"));
      await session.step(38, "And user sets \"rowSource\" property of line chart viewer to \"Selected\"", () => setProperty(page, "rowSource", el("line chart viewer"), "Selected"));
      await session.step(39, "Then the \"rows shown\" reading of line chart viewer should be 36", () => readingIs(page, "rows shown", el("line chart viewer"), 36));
      await session.step(40, "And the \"markers drawn\" reading of line chart viewer should be 36", () => readingIs(page, "markers drawn", el("line chart viewer"), 36));
      await session.step(41, "And the \"rows selected\" reading of line chart viewer should be 36", () => readingIs(page, "rows selected", el("line chart viewer"), 36));
      await session.step(42, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(43, "When user sets \"rowSource\" property of line chart viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("line chart viewer"), "Filtered"));
      await session.step(44, "Then the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(45, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(46, "When user clears the row selection", () => clearSelection(page));
      await session.step(47, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source Selected with nothing selected leaves the chart empty, split or not, then draws exactly the selected rows (public#1650)", async () => {
      await session.step(51, "Given user clears the row selection", () => clearSelection(page));
      await session.step(52, "When user sets \"rowSource\" property of line chart viewer to \"Selected\"", () => setProperty(page, "rowSource", el("line chart viewer"), "Selected"));
      await session.step(53, "Then the \"rows shown\" reading of line chart viewer should be 0", () => readingIs(page, "rows shown", el("line chart viewer"), 0));
      await session.step(54, "And the \"markers drawn\" reading of line chart viewer should be 0", () => readingIs(page, "markers drawn", el("line chart viewer"), 0));
      await session.step(55, "And the \"lines\" reading of line chart viewer should be 0", () => readingIs(page, "lines", el("line chart viewer"), 0));
      await session.step(56, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(57, "When user sets \"rowSource\" property of line chart viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("line chart viewer"), "Filtered"));
      await session.step(58, "And user sets \"splitColumnNames\" property of line chart viewer to \"Series\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Series"));
      await session.step(59, "Then the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(60, "When user sets \"rowSource\" property of line chart viewer to \"Selected\"", () => setProperty(page, "rowSource", el("line chart viewer"), "Selected"));
      await session.step(61, "Then the \"rows shown\" reading of line chart viewer should be 0", () => readingIs(page, "rows shown", el("line chart viewer"), 0));
      await session.step(62, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
      await session.step(64, "When user selects rows where \"Series\" is \"Pyrrolidines\"", () => selectWhereIs(page, "Series", "Pyrrolidines"));
      await session.step(65, "Then the \"rows shown\" reading of line chart viewer should be 21", () => readingIs(page, "rows shown", el("line chart viewer"), 21));
      await session.step(66, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
      await session.step(67, "When user selects rows where \"Series\" is one of \"Pyrrolidines, Aminopiperidines\"", () => selectWhereOneOf(page, "Series", "Pyrrolidines, Aminopiperidines"));
      await session.step(68, "Then the \"rows shown\" reading of line chart viewer should be 25", () => readingIs(page, "rows shown", el("line chart viewer"), 25));
      await session.step(69, "And the \"lines\" reading of line chart viewer should be 2", () => readingIs(page, "lines", el("line chart viewer"), 2));
      await session.step(70, "When user clears the row selection", () => clearSelection(page));
      await session.step(71, "And user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["splitColumnNames",""],["rowSource","Filtered"]]), [["splitColumnNames",""],["rowSource","Filtered"]]);
      await session.step(74, "Then the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source Filtered follows the table filter where All ignores it", async () => {
      await session.step(78, "When user filters rows where \"Average Mass\" is between 200 and 400", () => filterBetween(page, "Average Mass", 200, 400));
      await session.step(79, "Then 51 rows should pass the filter", () => filterPasses(page, 51));
      await session.step(80, "And the \"rows shown\" reading of line chart viewer should be 51", () => readingIs(page, "rows shown", el("line chart viewer"), 51));
      await session.step(81, "And the \"markers drawn\" reading of line chart viewer should be 51", () => readingIs(page, "markers drawn", el("line chart viewer"), 51));
      await session.step(82, "When user sets \"rowSource\" property of line chart viewer to \"All\"", () => setProperty(page, "rowSource", el("line chart viewer"), "All"));
      await session.step(83, "Then the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(84, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(85, "When user sets \"rowSource\" property of line chart viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("line chart viewer"), "Filtered"));
      await session.step(86, "Then the \"rows shown\" reading of line chart viewer should be 51", () => readingIs(page, "rows shown", el("line chart viewer"), 51));
      await session.step(87, "When user resets the filter", () => resetFilter(page));
      await session.step(88, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(89, "And the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The chart's own filter expression narrows it and leaves the table alone", async () => {
      await session.step(93, "When user sets \"filter\" property of line chart viewer to \"${CAST Idea ID} < 634834\"", () => setProperty(page, "filter", el("line chart viewer"), "${CAST Idea ID} < 634834"));
      await session.step(94, "Then the \"rows shown\" reading of line chart viewer should be 49", () => readingIs(page, "rows shown", el("line chart viewer"), 49));
      await session.step(95, "And the \"markers drawn\" reading of line chart viewer should be 49", () => readingIs(page, "markers drawn", el("line chart viewer"), 49));
      await session.step(96, "And 100 rows should pass the filter", () => filterPasses(page, 100));
      await session.step(97, "And the \"x axis span\" reading of line chart viewer should be lower than before", () => readingLower(page, "x axis span", el("line chart viewer")));
      await session.step(98, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(99, "When user sets \"filter\" property of line chart viewer to \"\"", () => setProperty(page, "filter", el("line chart viewer"), ""));
      await session.step(100, "Then the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(101, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The chart's own filter and the table filter narrow it together", async () => {
      await session.step(105, "When user sets \"filter\" property of line chart viewer to \"${CAST Idea ID} < 634834\"", () => setProperty(page, "filter", el("line chart viewer"), "${CAST Idea ID} < 634834"));
      await session.step(106, "And user filters rows where \"Average Mass\" is between 200 and 400", () => filterBetween(page, "Average Mass", 200, 400));
      await session.step(107, "Then 51 rows should pass the filter", () => filterPasses(page, 51));
      await session.step(108, "And the \"rows shown\" reading of line chart viewer should be lower than before", () => readingLower(page, "rows shown", el("line chart viewer")));
      await session.step(109, "And the \"rows shown\" reading of line chart viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("line chart viewer"), 1));
      await session.step(110, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(111, "When user sets \"filter\" property of line chart viewer to \"\"", () => setProperty(page, "filter", el("line chart viewer"), ""));
      await session.step(112, "And user resets the filter", () => resetFilter(page));
      await session.step(113, "Then the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bound to the other table the chart draws that table's rows and re-picks its columns", async () => {
      await session.step(117, "When user sets \"table\" property of line chart viewer to \"demog-1000\"", () => setProperty(page, "table", el("line chart viewer"), "demog-1000"));
      await session.step(118, "Then line chart viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("line chart viewer"), "demog-1000"));
      await session.step(119, "And the \"rows shown\" reading of line chart viewer should be 1000", () => readingIs(page, "rows shown", el("line chart viewer"), 1000));
      await session.step(120, "And the \"rows shown\" reading of line chart viewer should be higher than before", () => readingHigher(page, "rows shown", el("line chart viewer")));
      await session.step(121, "And the \"x column\" reading of line chart viewer should not be \"CAST Idea ID\"", () => readingDoesNotRead(page, "x column", el("line chart viewer"), "CAST Idea ID"));
      await session.step(122, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(123, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(124, "When user sets \"table\" property of line chart viewer to \"spgi-100\"", () => setProperty(page, "table", el("line chart viewer"), "spgi-100"));
      await session.step(125, "Then line chart viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("line chart viewer"), "spgi-100"));
      await session.step(126, "And the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
      await session.step(127, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
      await session.step(130, "Then the \"x column\" reading of line chart viewer should be \"CAST Idea ID\"", () => readingReads(page, "x column", el("line chart viewer"), "CAST Idea ID"));
      await session.step(131, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
