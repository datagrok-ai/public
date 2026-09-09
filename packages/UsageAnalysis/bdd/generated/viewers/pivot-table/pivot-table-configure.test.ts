/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-configure.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pivot-table]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addToRow, dragAreaOntoWidget, pivotedAggregationMatches} from '../../../bindings/pivot-table.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, closeContextMenu, hasArea, hasNoArea, menuLists, noErrors, pickFromAreaContextMenu, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — configuring the tag rows", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-configure.feature", import.meta.url);
  test("Pivot table — configuring the tag rows", {tag: ["@journey", "@viewers", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(19, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
    await session.step(20, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
    await session.step(21, "And the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
    await session.step(22, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await run.scenario("The + picker of the Group by row adds a second key column", async () => {
      await session.step(25, "When user adds \"SEX\" to the \"group by\" row of pivot table viewer", () => addToRow(page, "SEX", "group by"));
      await session.step(26, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP, SEX\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP, SEX"));
      await session.step(27, "And the \"key columns\" reading of pivot table viewer should be \"DIS_POP, SEX\"", () => readingReads(page, "key columns", el("pivot table viewer"), "DIS_POP, SEX"));
      await session.step(28, "And pivot table viewer should have a \"group by chip SEX\" area", () => hasArea(page, el("pivot table viewer"), "group by chip SEX"));
      await session.step(29, "And the \"aggregated rows\" reading of pivot table viewer should be 12", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 12));
      await session.step(30, "And the \"text of grid cell 1 of SEX\" reading of pivot table viewer should be \"F\"", () => readingReads(page, "text of grid cell 1 of SEX", el("pivot table viewer"), "F"));
      await session.step(31, "When user clicks on the \"remove group by chip SEX\" area of pivot table viewer", () => clickArea(page, "remove group by chip SEX", el("pivot table viewer")));
      await session.step(32, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
      await session.step(33, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column header dragged from the main grid onto the Group by row becomes a key", async () => {
      await session.step(37, "When user drags the \"header RACE\" area of grid onto the \"group by row\" area of pivot table viewer", () => dragAreaOntoWidget(page, "header RACE", el("grid"), "group by row", el("pivot table viewer")));
      await session.step(38, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP, RACE\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP, RACE"));
      await session.step(39, "And pivot table viewer should have a \"group by chip RACE\" area", () => hasArea(page, el("pivot table viewer"), "group by chip RACE"));
      await session.step(40, "And the \"aggregated rows\" reading of pivot table viewer should be 24", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 24));
      await session.step(41, "When user clicks on the \"remove group by chip RACE\" area of pivot table viewer", () => clickArea(page, "remove group by chip RACE", el("pivot table viewer")));
      await session.step(42, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
      await session.step(43, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column header dragged onto the Aggregate row becomes an aggregation", async () => {
      await session.step(47, "When user drags the \"header HEIGHT\" area of grid onto the \"aggregate row\" area of pivot table viewer", () => dragAreaOntoWidget(page, "header HEIGHT", el("grid"), "aggregate row", el("pivot table viewer")));
      await session.step(48, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE), avg(HEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE), avg(HEIGHT)"));
      await session.step(49, "And the \"aggregated columns\" reading of pivot table viewer should be 11", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 11));
      await session.step(50, "When user clicks on the \"remove aggregate chip avg(HEIGHT)\" area of pivot table viewer", () => clickArea(page, "remove aggregate chip avg(HEIGHT)", el("pivot table viewer")));
      await session.step(51, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(52, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The chip menu takes two aggregation picks without closing", async () => {
      await session.step(56, "When user picks \"Aggregation > sum\" from the context menu of the \"aggregate chip avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Aggregation > sum", "aggregate chip avg(AGE)", el("pivot table viewer")));
      await session.step(57, "Then the \"aggregate\" reading of pivot table viewer should be \"sum(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "sum(AGE)"));
      await session.step(58, "And the \"aggregations\" reading of pivot table viewer should be \"sum(AGE)\"", () => readingReads(page, "aggregations", el("pivot table viewer"), "sum(AGE)"));
      await session.step(59, "And the aggregated values of pivot table viewer should match \"sum(AGE)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "sum(AGE)", "DIS_POP", "SEVERITY"));
      await session.step(60, "And the open menu should list \"Aggregation > med\"", () => menuLists(page, "Aggregation > med"));
      await session.step(61, "When user clicks on \"med\" menu item in context menu", () => clickOn(page, el("\"med\" menu item in context menu")));
      await session.step(62, "Then the \"aggregate\" reading of pivot table viewer should be \"med(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "med(AGE)"));
      await session.step(63, "And the aggregated values of pivot table viewer should match \"med(AGE)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "med(AGE)", "DIS_POP", "SEVERITY"));
      await session.step(64, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(65, "And user picks \"Aggregation > avg\" from the context menu of the \"aggregate chip med(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Aggregation > avg", "aggregate chip med(AGE)", el("pivot table viewer")));
      await session.step(66, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(67, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(68, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Repointing the chip at another column rebuilds the offered aggregations", async () => {
      await session.step(71, "When user picks \"Column > HEIGHT\" from the context menu of the \"aggregate chip avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Column > HEIGHT", "aggregate chip avg(AGE)", el("pivot table viewer")));
      await session.step(72, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(HEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(HEIGHT)"));
      await session.step(73, "And the open menu should list \"Aggregation > geomean\"", () => menuLists(page, "Aggregation > geomean"));
      await session.step(74, "And the open menu should list \"Aggregation > stdev\"", () => menuLists(page, "Aggregation > stdev"));
      await session.step(75, "And \"Aggregation > avg\" menu item in context menu should be selected", () => shouldBe(page, el("\"Aggregation > avg\" menu item in context menu"), "selected"));
      await session.step(76, "And \"Aggregation > sum\" menu item in context menu should not be selected", () => shouldNotBe(page, el("\"Aggregation > sum\" menu item in context menu"), "selected"));
      await session.step(77, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(78, "Then the aggregated values of pivot table viewer should match \"avg(HEIGHT)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "avg(HEIGHT)", "DIS_POP", "SEVERITY"));
      await session.step(79, "When user picks \"Column > AGE\" from the context menu of the \"aggregate chip avg(HEIGHT)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Column > AGE", "aggregate chip avg(HEIGHT)", el("pivot table viewer")));
      await session.step(80, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(81, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Remove others leaves the chip it was opened on", async () => {
      await session.step(85, "When user adds \"HEIGHT\" to the \"aggregate\" row of pivot table viewer", () => addToRow(page, "HEIGHT", "aggregate"));
      await session.step(86, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE), avg(HEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE), avg(HEIGHT)"));
      await session.step(87, "When user picks \"Remove others\" from the context menu of the \"aggregate chip avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Remove others", "aggregate chip avg(AGE)", el("pivot table viewer")));
      await session.step(88, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(89, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(90, "And pivot table viewer should not have a \"aggregate chip avg(HEIGHT)\" area", () => hasNoArea(page, el("pivot table viewer"), "aggregate chip avg(HEIGHT)"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the last aggregate hides the Pivot row and clears the pivot columns", async () => {
      await session.step(94, "When user clicks on the \"remove aggregate chip avg(AGE)\" area of pivot table viewer", () => clickArea(page, "remove aggregate chip avg(AGE)", el("pivot table viewer")));
      await session.step(95, "Then the \"aggregate\" reading of pivot table viewer should be \"\"", () => readingReads(page, "aggregate", el("pivot table viewer"), ""));
      await session.step(96, "And the \"pivot\" reading of pivot table viewer should be \"\"", () => readingReads(page, "pivot", el("pivot table viewer"), ""));
      await session.step(97, "And the \"pivot row shown\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "pivot row shown", el("pivot table viewer"), "false"));
      await session.step(98, "And pivot table viewer should not have a \"pivot row\" area", () => hasNoArea(page, el("pivot table viewer"), "pivot row"));
      await session.step(99, "And the \"aggregated columns\" reading of pivot table viewer should be 1", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 1));
      await session.step(100, "And the \"text of grid cell 5 of DIS_POP\" reading of pivot table viewer should be \"RA\"", () => readingReads(page, "text of grid cell 5 of DIS_POP", el("pivot table viewer"), "RA"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Re-adding an aggregate brings the Pivot row back, empty", async () => {
      await session.step(104, "When user adds \"AGE\" to the \"aggregate\" row of pivot table viewer", () => addToRow(page, "AGE", "aggregate"));
      await session.step(105, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(106, "And the \"pivot row shown\" reading of pivot table viewer should be \"true\"", () => readingReads(page, "pivot row shown", el("pivot table viewer"), "true"));
      await session.step(107, "And pivot table viewer should have a \"pivot row\" area", () => hasArea(page, el("pivot table viewer"), "pivot row"));
      await session.step(108, "And the \"pivot\" reading of pivot table viewer should be \"\"", () => readingReads(page, "pivot", el("pivot table viewer"), ""));
      await session.step(109, "And the \"aggregated columns\" reading of pivot table viewer should be 2", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 2));
      await session.step(110, "When user adds \"SEVERITY\" to the \"pivot\" row of pivot table viewer", () => addToRow(page, "SEVERITY", "pivot"));
      await session.step(111, "Then the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
      await session.step(112, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Refresh reseeds the first two numerical columns and clears the rest", async () => {
      await session.step(116, "When user clicks on the \"refresh\" area of pivot table viewer", () => clickArea(page, "refresh", el("pivot table viewer")));
      await session.step(117, "Then the \"group by\" reading of pivot table viewer should be \"\"", () => readingReads(page, "group by", el("pivot table viewer"), ""));
      await session.step(118, "And the \"pivot\" reading of pivot table viewer should be \"\"", () => readingReads(page, "pivot", el("pivot table viewer"), ""));
      await session.step(119, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE), avg(HEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE), avg(HEIGHT)"));
      await session.step(120, "And the \"aggregated rows\" reading of pivot table viewer should be 1", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 1));
      await session.step(121, "And the \"aggregated columns\" reading of pivot table viewer should be 2", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 2));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The column lists written as properties rebuild the chips", async () => {
      await session.step(125, "When user sets \"Group By Column Names\" property of pivot table viewer to \"DIS_POP\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "DIS_POP"));
      await session.step(126, "And user sets \"Aggregate Column Names\" property of pivot table viewer to \"WEIGHT\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "WEIGHT"));
      await session.step(127, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(128, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(129, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
      await session.step(130, "And the \"aggregate\" reading of pivot table viewer should be \"avg(WEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(WEIGHT)"));
      await session.step(131, "And the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
      await session.step(132, "And pivot table viewer should have a \"aggregate chip avg(WEIGHT)\" area", () => hasArea(page, el("pivot table viewer"), "aggregate chip avg(WEIGHT)"));
      await session.step(133, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(134, "And the aggregated values of pivot table viewer should match \"avg(WEIGHT)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "avg(WEIGHT)", "DIS_POP", "SEVERITY"));
      await session.step(135, "When user sets \"Aggregate Column Names\" property of pivot table viewer to \"AGE\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "AGE"));
      await session.step(136, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(137, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(138, "And the \"text of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"52.30\"", () => readingReads(page, "text of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "52.30"));
      await session.step(139, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
