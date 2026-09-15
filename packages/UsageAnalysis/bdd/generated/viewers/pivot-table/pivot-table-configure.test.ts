/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-configure.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pivot-table]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addToRow, aggregationMatches, pivotedAggregationMatches} from '../../../bindings/pivot-table.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, focusOn, pressKey, shouldBe, shouldNotBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, closeContextMenu, hasArea, hasNoArea, menuLists, noErrors, pickFromAreaContextMenu, pointerAway, readingIs, readingReads, setProperties} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {checkedInColumnList, columnListStartsWith, dragAreaOntoWidget, toggleInColumnList, uncheckedInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — configuring the tag rows", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-configure.feature", import.meta.url);
  test("Pivot table — configuring the tag rows", {tag: ["@journey", "@viewers", "@realizes:viewers.pivot-table", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(25, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(26, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
    await session.step(27, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
    await session.step(28, "And the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
    await session.step(29, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await run.scenario("The + picker of the Group by row adds a second key column", async () => {
      await session.step(32, "When user clicks on the \"add group by\" area of pivot table viewer", () => clickArea(page, "add group by", el("pivot table viewer")));
      await session.step(33, "And user moves the pointer away from pivot table viewer", () => pointerAway(page, el("pivot table viewer")));
      await session.step(34, "Then column picker popup should be visible", () => shouldBe(page, el("column picker popup"), "visible"));
      await session.step(35, "When user clicks on status bar", () => clickOn(page, el("status bar")));
      await session.step(36, "Then column picker popup should be absent", () => shouldBe(page, el("column picker popup"), "absent"));
      await session.step(37, "And the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
      await session.step(39, "When user adds \"SEX\" to the \"group by\" row of pivot table viewer", () => addToRow(page, "SEX", "group by"));
      await session.step(40, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP, SEX\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP, SEX"));
      await session.step(41, "And the \"key columns\" reading of pivot table viewer should be \"DIS_POP, SEX\"", () => readingReads(page, "key columns", el("pivot table viewer"), "DIS_POP, SEX"));
      await session.step(42, "And pivot table viewer should have a \"group by chip SEX\" area", () => hasArea(page, el("pivot table viewer"), "group by chip SEX"));
      await session.step(43, "And the \"aggregated rows\" reading of pivot table viewer should be 12", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 12));
      await session.step(44, "And the \"text of grid cell 1 of SEX\" reading of pivot table viewer should be \"F\"", () => readingReads(page, "text of grid cell 1 of SEX", el("pivot table viewer"), "F"));
      await session.step(45, "When user clicks on the \"remove group by chip SEX\" area of pivot table viewer", () => clickArea(page, "remove group by chip SEX", el("pivot table viewer")));
      await session.step(46, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
      await session.step(47, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column header dragged from the main grid onto the Group by row becomes a key", async () => {
      await session.step(51, "When user drags the \"header RACE\" area of grid onto the \"group by row\" area of pivot table viewer", () => dragAreaOntoWidget(page, "header RACE", el("grid"), "group by row", el("pivot table viewer")));
      await session.step(52, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP, RACE\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP, RACE"));
      await session.step(53, "And pivot table viewer should have a \"group by chip RACE\" area", () => hasArea(page, el("pivot table viewer"), "group by chip RACE"));
      await session.step(54, "And the \"aggregated rows\" reading of pivot table viewer should be 24", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 24));
      await session.step(55, "When user clicks on the \"remove group by chip RACE\" area of pivot table viewer", () => clickArea(page, "remove group by chip RACE", el("pivot table viewer")));
      await session.step(56, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
      await session.step(57, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column header dragged onto the Aggregate row becomes an aggregation", async () => {
      await session.step(61, "When user drags the \"header HEIGHT\" area of grid onto the \"aggregate row\" area of pivot table viewer", () => dragAreaOntoWidget(page, "header HEIGHT", el("grid"), "aggregate row", el("pivot table viewer")));
      await session.step(62, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE), avg(HEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE), avg(HEIGHT)"));
      await session.step(63, "And the \"aggregated columns\" reading of pivot table viewer should be 11", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 11));
      await session.step(64, "When user clicks on the \"remove aggregate chip avg(HEIGHT)\" area of pivot table viewer", () => clickArea(page, "remove aggregate chip avg(HEIGHT)", el("pivot table viewer")));
      await session.step(65, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(66, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The chip menu takes two aggregation picks without closing", async () => {
      await session.step(70, "When user picks \"Aggregation > sum\" from the context menu of the \"aggregate chip avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Aggregation > sum", "aggregate chip avg(AGE)", el("pivot table viewer")));
      await session.step(71, "Then the \"aggregate\" reading of pivot table viewer should be \"sum(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "sum(AGE)"));
      await session.step(72, "And the \"aggregations\" reading of pivot table viewer should be \"sum(AGE)\"", () => readingReads(page, "aggregations", el("pivot table viewer"), "sum(AGE)"));
      await session.step(73, "And the aggregated values of pivot table viewer should match \"sum(AGE)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "sum(AGE)", "DIS_POP", "SEVERITY"));
      await session.step(74, "And the open menu should list \"Aggregation > med\"", () => menuLists(page, "Aggregation > med"));
      await session.step(75, "When user clicks on \"med\" menu item in context menu", () => clickOn(page, el("\"med\" menu item in context menu")));
      await session.step(76, "Then the \"aggregate\" reading of pivot table viewer should be \"med(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "med(AGE)"));
      await session.step(77, "And the aggregated values of pivot table viewer should match \"med(AGE)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "med(AGE)", "DIS_POP", "SEVERITY"));
      await session.step(78, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(79, "And user picks \"Aggregation > avg\" from the context menu of the \"aggregate chip med(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Aggregation > avg", "aggregate chip med(AGE)", el("pivot table viewer")));
      await session.step(80, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(81, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(82, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Repointing the chip at another column rebuilds the offered aggregations", async () => {
      await session.step(85, "When user picks \"Column > HEIGHT\" from the context menu of the \"aggregate chip avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Column > HEIGHT", "aggregate chip avg(AGE)", el("pivot table viewer")));
      await session.step(86, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(HEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(HEIGHT)"));
      await session.step(87, "And the open menu should list \"Aggregation > geomean\"", () => menuLists(page, "Aggregation > geomean"));
      await session.step(88, "And the open menu should list \"Aggregation > stdev\"", () => menuLists(page, "Aggregation > stdev"));
      await session.step(89, "And \"Aggregation > avg\" menu item in context menu should be selected", () => shouldBe(page, el("\"Aggregation > avg\" menu item in context menu"), "selected"));
      await session.step(90, "And \"Aggregation > sum\" menu item in context menu should not be selected", () => shouldNotBe(page, el("\"Aggregation > sum\" menu item in context menu"), "selected"));
      await session.step(91, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(92, "Then the aggregated values of pivot table viewer should match \"avg(HEIGHT)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "avg(HEIGHT)", "DIS_POP", "SEVERITY"));
      await session.step(93, "When user picks \"Column > AGE\" from the context menu of the \"aggregate chip avg(HEIGHT)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Column > AGE", "aggregate chip avg(HEIGHT)", el("pivot table viewer")));
      await session.step(94, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(95, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Remove others leaves the chip it was opened on", async () => {
      await session.step(99, "When user adds \"HEIGHT\" to the \"aggregate\" row of pivot table viewer", () => addToRow(page, "HEIGHT", "aggregate"));
      await session.step(100, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE), avg(HEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE), avg(HEIGHT)"));
      await session.step(101, "When user picks \"Remove others\" from the context menu of the \"aggregate chip avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Remove others", "aggregate chip avg(AGE)", el("pivot table viewer")));
      await session.step(102, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(103, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(104, "And pivot table viewer should not have a \"aggregate chip avg(HEIGHT)\" area", () => hasNoArea(page, el("pivot table viewer"), "aggregate chip avg(HEIGHT)"));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the last aggregate hides the Pivot row and clears the pivot columns", async () => {
      await session.step(108, "When user clicks on the \"remove aggregate chip avg(AGE)\" area of pivot table viewer", () => clickArea(page, "remove aggregate chip avg(AGE)", el("pivot table viewer")));
      await session.step(109, "Then the \"aggregate\" reading of pivot table viewer should be \"\"", () => readingReads(page, "aggregate", el("pivot table viewer"), ""));
      await session.step(110, "And the \"pivot\" reading of pivot table viewer should be \"\"", () => readingReads(page, "pivot", el("pivot table viewer"), ""));
      await session.step(111, "And the \"pivot row shown\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "pivot row shown", el("pivot table viewer"), "false"));
      await session.step(112, "And pivot table viewer should not have a \"pivot row\" area", () => hasNoArea(page, el("pivot table viewer"), "pivot row"));
      await session.step(113, "And the \"aggregated columns\" reading of pivot table viewer should be 1", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 1));
      await session.step(114, "And the \"text of grid cell 5 of DIS_POP\" reading of pivot table viewer should be \"RA\"", () => readingReads(page, "text of grid cell 5 of DIS_POP", el("pivot table viewer"), "RA"));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Re-adding an aggregate brings the Pivot row back, empty", async () => {
      await session.step(118, "When user adds \"AGE\" to the \"aggregate\" row of pivot table viewer", () => addToRow(page, "AGE", "aggregate"));
      await session.step(119, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(120, "And the \"pivot row shown\" reading of pivot table viewer should be \"true\"", () => readingReads(page, "pivot row shown", el("pivot table viewer"), "true"));
      await session.step(121, "And pivot table viewer should have a \"pivot row\" area", () => hasArea(page, el("pivot table viewer"), "pivot row"));
      await session.step(122, "And the \"pivot\" reading of pivot table viewer should be \"\"", () => readingReads(page, "pivot", el("pivot table viewer"), ""));
      await session.step(123, "And the \"aggregated columns\" reading of pivot table viewer should be 2", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 2));
      await session.step(124, "When user adds \"SEVERITY\" to the \"pivot\" row of pivot table viewer", () => addToRow(page, "SEVERITY", "pivot"));
      await session.step(125, "Then the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
      await session.step(126, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Refresh reseeds the first two numerical columns and clears the rest", async () => {
      await session.step(130, "When user clicks on the \"refresh\" area of pivot table viewer", () => clickArea(page, "refresh", el("pivot table viewer")));
      await session.step(131, "Then the \"group by\" reading of pivot table viewer should be \"\"", () => readingReads(page, "group by", el("pivot table viewer"), ""));
      await session.step(132, "And the \"pivot\" reading of pivot table viewer should be \"\"", () => readingReads(page, "pivot", el("pivot table viewer"), ""));
      await session.step(133, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE), avg(HEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE), avg(HEIGHT)"));
      await session.step(134, "And the \"aggregated rows\" reading of pivot table viewer should be 1", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 1));
      await session.step(135, "And the \"aggregated columns\" reading of pivot table viewer should be 2", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 2));
      await session.step(136, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The column lists edited in the context panel rebuild the chips", async () => {
      await session.step(139, "When user sets properties of pivot table viewer:", () => setProperties(page, el("pivot table viewer"), [["Group By Column Names","DIS_POP"],["Aggregate Column Names","AGE"],["Aggregate Agg Types","avg"],["Pivot Column Names","SEVERITY"]]));
      await session.step(144, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(145, "And the \"text of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"52.30\"", () => readingReads(page, "text of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "52.30"));
      await session.step(146, "When user clicks on settings icon of pivot table viewer", () => clickOn(page, el("settings icon of pivot table viewer")));
      await session.step(147, "And user clicks on \"...\" button in \"Aggregate\" property", () => clickOn(page, el("\"...\" button in \"Aggregate\" property")));
      await session.step(148, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(149, "When user types \"AGE\" into \"Search\" input in \"Select columns...\" dialog", () => typeInto(page, "AGE", el("\"Search\" input in \"Select columns...\" dialog")));
      await session.step(150, "Then the column list of \"Select columns...\" dialog should start with \"AGE\"", () => columnListStartsWith(page, el("\"Select columns...\" dialog"), "AGE"));
      await session.step(151, "And the \"AGE\" column should be checked in the column list of \"Select columns...\" dialog", () => checkedInColumnList(page, "AGE", el("\"Select columns...\" dialog")));
      await session.step(152, "When user toggles the \"AGE\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "AGE", el("\"Select columns...\" dialog")));
      await session.step(153, "And user focuses on \"Search\" input in \"Select columns...\" dialog", () => focusOn(page, el("\"Search\" input in \"Select columns...\" dialog")));
      await session.step(154, "And user types \"WEIGHT\" into \"Search\" input in \"Select columns...\" dialog", () => typeInto(page, "WEIGHT", el("\"Search\" input in \"Select columns...\" dialog")));
      await session.step(155, "Then the column list of \"Select columns...\" dialog should start with \"WEIGHT\"", () => columnListStartsWith(page, el("\"Select columns...\" dialog"), "WEIGHT"));
      await session.step(156, "And the \"WEIGHT\" column should not be checked in the column list of \"Select columns...\" dialog", () => uncheckedInColumnList(page, "WEIGHT", el("\"Select columns...\" dialog")));
      await session.step(157, "When user toggles the \"WEIGHT\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "WEIGHT", el("\"Select columns...\" dialog")));
      await session.step(158, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(159, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
      await session.step(160, "And the \"aggregate\" reading of pivot table viewer should be \"avg(WEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(WEIGHT)"));
      await session.step(161, "And pivot table viewer should have a \"aggregate chip avg(WEIGHT)\" area", () => hasArea(page, el("pivot table viewer"), "aggregate chip avg(WEIGHT)"));
      await session.step(162, "And the \"pivot\" reading of pivot table viewer should be \"\"", () => readingReads(page, "pivot", el("pivot table viewer"), ""));
      await session.step(163, "And the aggregated values of pivot table viewer should match \"avg(WEIGHT)\" grouped by \"DIS_POP\"", () => aggregationMatches(page, el("pivot table viewer"), "avg(WEIGHT)", "DIS_POP"));
      await session.step(164, "When user clicks on \"...\" button in \"Group By\" property", () => clickOn(page, el("\"...\" button in \"Group By\" property")));
      await session.step(165, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(166, "When user types \"DIS_POP\" into \"Search\" input in \"Select columns...\" dialog", () => typeInto(page, "DIS_POP", el("\"Search\" input in \"Select columns...\" dialog")));
      await session.step(167, "Then the column list of \"Select columns...\" dialog should start with \"DIS_POP\"", () => columnListStartsWith(page, el("\"Select columns...\" dialog"), "DIS_POP"));
      await session.step(168, "When user toggles the \"DIS_POP\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "DIS_POP", el("\"Select columns...\" dialog")));
      await session.step(169, "And user focuses on \"Search\" input in \"Select columns...\" dialog", () => focusOn(page, el("\"Search\" input in \"Select columns...\" dialog")));
      await session.step(170, "And user types \"SEX\" into \"Search\" input in \"Select columns...\" dialog", () => typeInto(page, "SEX", el("\"Search\" input in \"Select columns...\" dialog")));
      await session.step(171, "Then the column list of \"Select columns...\" dialog should start with \"SEX\"", () => columnListStartsWith(page, el("\"Select columns...\" dialog"), "SEX"));
      await session.step(172, "When user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
      await session.step(173, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(174, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
      await session.step(175, "And the \"group by\" reading of pivot table viewer should be \"SEX\"", () => readingReads(page, "group by", el("pivot table viewer"), "SEX"));
      await session.step(176, "And pivot table viewer should have a \"group by chip SEX\" area", () => hasArea(page, el("pivot table viewer"), "group by chip SEX"));
      await session.step(177, "And the \"aggregated rows\" reading of pivot table viewer should be 2", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 2));
      await session.step(178, "And the aggregated values of pivot table viewer should match \"avg(WEIGHT)\" grouped by \"SEX\"", () => aggregationMatches(page, el("pivot table viewer"), "avg(WEIGHT)", "SEX"));
      await session.step(179, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Escape closes the + picker of the Group by row without taking a column (GROK-20900)", async () => {
      await session.step(183, "Given the \"group by\" reading of pivot table viewer should be \"SEX\"", () => readingReads(page, "group by", el("pivot table viewer"), "SEX"));
      await session.step(184, "When user clicks on the \"add group by\" area of pivot table viewer", () => clickArea(page, "add group by", el("pivot table viewer")));
      await session.step(185, "And user moves the pointer away from pivot table viewer", () => pointerAway(page, el("pivot table viewer")));
      await session.step(186, "Then column picker popup should be visible", () => shouldBe(page, el("column picker popup"), "visible"));
      await session.step(187, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(188, "Then column picker popup should be absent", () => shouldBe(page, el("column picker popup"), "absent"));
      await session.step(189, "And the \"group by\" reading of pivot table viewer should be \"SEX\"", () => readingReads(page, "group by", el("pivot table viewer"), "SEX"));
      await session.step(190, "And no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true});
    run.finish();
  });
});
