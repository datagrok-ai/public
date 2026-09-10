/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/histogram/histogram-selection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.histogram]
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
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, filterPasses, filterPassesAll, filterTo, noneSelected, onlyOfAnySelected, onlyOfSelected, resetFilter, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaRepainted, clickArea, clickAreaHolding, dragSelectionBetweenAreas, eventFired, hasArea, hasNoArea, hoverArea, listenFor, moreHighlight, noErrors, noHighlight, notRepainted, pointerAway, repainted, repaintedBy, resizeTo, setProperty, showsRows, someHighlight, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Histogram bin selection, row markers and mouse-over", () => {
  const session = feature(test, "features/viewers/histogram/histogram-selection.feature", import.meta.url);
  test("Histogram bin selection, row markers and mouse-over", {tag: ["@journey", "@viewers", "@realizes:viewers.histogram"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"]]));
    await session.step(15, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(16, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(17, "And histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
    await session.step(18, "And histogram viewer should show no selection highlight", () => noHighlight(page, el("histogram viewer")));
    await session.step(19, "And histogram viewer should not have a \"selected bin 8\" area", () => hasNoArea(page, el("histogram viewer"), "selected bin 8"));
    await run.scenario("A click on a bin selects the rows of that bin", async () => {
      await session.step(22, "Given user listens for \"d4-histogram-select-bins\" event on histogram viewer", () => listenFor(page, "d4-histogram-select-bins", el("histogram viewer")));
      await session.step(23, "When user clicks on the \"bin 8\" area of histogram viewer", () => clickArea(page, "bin 8", el("histogram viewer")));
      await session.step(24, "Then \"d4-histogram-select-bins\" event should have fired on histogram viewer", () => eventFired(page, "d4-histogram-select-bins", el("histogram viewer")));
      await session.step(25, "And only rows where \"AGE\" is one of \"43, 44, 45, 46\" should be selected", () => onlyOfAnySelected(page, "AGE", "43, 44, 45, 46"));
      await session.step(26, "And 99 rows should be selected", () => selectedRowCount(page, 99));
      await session.step(27, "And histogram viewer should show a selection highlight", () => someHighlight(page, el("histogram viewer")));
      await session.step(28, "And histogram viewer should have a \"selected bin 8\" area", () => hasArea(page, el("histogram viewer"), "selected bin 8"));
      await session.step(29, "And the \"selected bin 8\" area of histogram viewer should contain the color \"#FF8C00\"", () => areaColor(page, "selected bin 8", el("histogram viewer"), "#FF8C00"));
      await session.step(30, "And histogram viewer should not have a \"selected bin 1\" area", () => hasNoArea(page, el("histogram viewer"), "selected bin 1"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control adds a bin to the selection", async () => {
      await session.step(34, "When user clicks on the \"bin 9\" area of histogram viewer holding Control", () => clickAreaHolding(page, "bin 9", el("histogram viewer"), "Control"));
      await session.step(35, "Then only rows where \"AGE\" is one of \"43, 44, 45, 46, 47, 48, 49\" should be selected", () => onlyOfAnySelected(page, "AGE", "43, 44, 45, 46, 47, 48, 49"));
      await session.step(36, "And 180 rows should be selected", () => selectedRowCount(page, 180));
      await session.step(37, "And histogram viewer should show more selection highlight than before", () => moreHighlight(page, el("histogram viewer")));
      await session.step(38, "And histogram viewer should have a \"selected bin 9\" area", () => hasArea(page, el("histogram viewer"), "selected bin 9"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clearing the selection takes the overlay with it", async () => {
      await session.step(42, "When user clears the row selection", () => clearSelection(page));
      await session.step(43, "Then no rows should be selected", () => noneSelected(page));
      await session.step(44, "And histogram viewer should show no selection highlight", () => noHighlight(page, el("histogram viewer")));
      await session.step(45, "And histogram viewer should not have a \"selected bin 8\" area", () => hasNoArea(page, el("histogram viewer"), "selected bin 8"));
      await session.step(46, "And histogram viewer should not have a \"selected bin 9\" area", () => hasNoArea(page, el("histogram viewer"), "selected bin 9"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Shift drag selects every bin it crosses", async () => {
      await session.step(50, "When user drags a selection box from the \"bin 6\" area to the \"bin 7\" area of histogram viewer", () => dragSelectionBetweenAreas(page, "bin 6", "bin 7", el("histogram viewer")));
      await session.step(51, "Then only rows where \"AGE\" is one of \"36, 37, 38, 39, 40, 41, 42\" should be selected", () => onlyOfAnySelected(page, "AGE", "36, 37, 38, 39, 40, 41, 42"));
      await session.step(52, "And histogram viewer should show a selection highlight", () => someHighlight(page, el("histogram viewer")));
      await session.step(53, "And histogram viewer should have a \"selected bin 6\" area", () => hasArea(page, el("histogram viewer"), "selected bin 6"));
      await session.step(54, "And histogram viewer should have a \"selected bin 7\" area", () => hasArea(page, el("histogram viewer"), "selected bin 7"));
      await session.step(55, "When user clears the row selection", () => clearSelection(page));
      await session.step(56, "Then no rows should be selected", () => noneSelected(page));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A selected row source shows the selection instead of overlaying it", async () => {
      await session.step(60, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(61, "And user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(62, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(63, "And histogram viewer should show 553 rows", () => showsRows(page, el("histogram viewer"), 553));
      await session.step(64, "And histogram viewer should show a selection highlight", () => someHighlight(page, el("histogram viewer")));
      await session.step(65, "When user sets \"Row Source\" property of histogram viewer to \"Selected\"", () => setProperty(page, "Row Source", el("histogram viewer"), "Selected"));
      await session.step(66, "Then histogram viewer should show 15 rows", () => showsRows(page, el("histogram viewer"), 15));
      await session.step(67, "And histogram viewer should show no selection highlight", () => noHighlight(page, el("histogram viewer")));
      await session.step(68, "And only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(69, "When user sets \"Row Source\" property of histogram viewer to \"All\"", () => setProperty(page, "Row Source", el("histogram viewer"), "All"));
      await session.step(70, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
      await session.step(71, "And histogram viewer should show a selection highlight", () => someHighlight(page, el("histogram viewer")));
      await session.step(72, "When user sets \"Row Source\" property of histogram viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("histogram viewer"), "Filtered"));
      await session.step(73, "And user resets the filter", () => resetFilter(page));
      await session.step(74, "And user clears the row selection", () => clearSelection(page));
      await session.step(75, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(76, "And no rows should be selected", () => noneSelected(page));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The current row is a dot on the baseline", async () => {
      await session.step(80, "When user makes row 8 current", () => makeRowCurrent(page, 8));
      await session.step(81, "Then row 8 should be current", () => currentRowIs(page, 8));
      await session.step(82, "And histogram viewer should have a \"current row marker\" area", () => hasArea(page, el("histogram viewer"), "current row marker"));
      await session.step(83, "And the \"current row marker\" area of histogram viewer should contain the color \"#38B738\"", () => areaColor(page, "current row marker", el("histogram viewer"), "#38B738"));
      await session.step(84, "When user sets \"Show Current Row\" property of histogram viewer to \"false\"", () => setProperty(page, "Show Current Row", el("histogram viewer"), "false"));
      await session.step(85, "Then histogram viewer should not have a \"current row marker\" area", () => hasNoArea(page, el("histogram viewer"), "current row marker"));
      await session.step(86, "And histogram viewer should have repainted", () => repainted(page, el("histogram viewer")));
      await session.step(87, "When user sets \"Show Current Row\" property of histogram viewer to \"true\"", () => setProperty(page, "Show Current Row", el("histogram viewer"), "true"));
      await session.step(88, "Then histogram viewer should have a \"current row marker\" area", () => hasArea(page, el("histogram viewer"), "current row marker"));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The mouse-over row marker follows the grid", async () => {
      await session.step(92, "Given histogram viewer should not have a \"mouse over row marker\" area", () => hasNoArea(page, el("histogram viewer"), "mouse over row marker"));
      await session.step(93, "When user hovers over the \"cell 8 of AGE\" area of grid", () => hoverArea(page, "cell 8 of AGE", el("grid")));
      await session.step(94, "Then histogram viewer should have a \"mouse over row marker\" area", () => hasArea(page, el("histogram viewer"), "mouse over row marker"));
      await session.step(95, "And the \"mouse over row marker\" area of histogram viewer should contain the color \"#AAAAAA\"", () => areaColor(page, "mouse over row marker", el("histogram viewer"), "#AAAAAA"));
      await session.step(96, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(97, "Then histogram viewer should not have a \"mouse over row marker\" area", () => hasNoArea(page, el("histogram viewer"), "mouse over row marker"));
      await session.step(98, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hovering a bin shows its tooltip and dims the other bins", async () => {
      await session.step(101, "Given user listens for \"d4-histogram-mouse-over-bins\" event on histogram viewer", () => listenFor(page, "d4-histogram-mouse-over-bins", el("histogram viewer")));
      await session.step(102, "When user hovers over the \"bin 8\" area of histogram viewer", () => hoverArea(page, "bin 8", el("histogram viewer")));
      await session.step(103, "Then \"d4-histogram-mouse-over-bins\" event should have fired on histogram viewer", () => eventFired(page, "d4-histogram-mouse-over-bins", el("histogram viewer")));
      await session.step(104, "And the \"bin 1\" area of histogram viewer should have repainted", () => areaRepainted(page, "bin 1", el("histogram viewer")));
      await session.step(105, "And histogram viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("histogram viewer"), 500));
      await session.step(106, "And tooltip should contain text \"total: 99\"", () => shouldContainText(page, el("tooltip"), "total: 99"));
      await session.step(107, "And tooltip should contain text \"filtered: 99\"", () => shouldContainText(page, el("tooltip"), "filtered: 99"));
      await session.step(108, "When user moves the pointer away from histogram viewer", () => pointerAway(page, el("histogram viewer")));
      await session.step(109, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A mouse-over row group in another viewer repaints the bins", async () => {
      await session.step(112, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]));
      await session.step(114, "And user takes a snapshot of histogram viewer", () => takeSnapshot(page, el("histogram viewer")));
      await session.step(115, "And user hovers over the \"bar Asian\" area of bar chart viewer", () => hoverArea(page, "bar Asian", el("bar chart viewer")));
      await session.step(116, "Then histogram viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("histogram viewer"), 500));
      await session.step(117, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
      await session.step(118, "And user moves the pointer away from histogram viewer", () => pointerAway(page, el("histogram viewer")));
      await session.step(119, "And user sets \"Show Mouse Over Row Group\" property of histogram viewer to \"false\"", () => setProperty(page, "Show Mouse Over Row Group", el("histogram viewer"), "false"));
      await session.step(120, "And user hovers over the \"bar Black\" area of bar chart viewer", () => hoverArea(page, "bar Black", el("bar chart viewer")));
      await session.step(121, "Then histogram viewer should not have repainted", () => notRepainted(page, el("histogram viewer")));
      await session.step(122, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
      await session.step(123, "And user sets \"Show Mouse Over Row Group\" property of histogram viewer to \"true\"", () => setProperty(page, "Show Mouse Over Row Group", el("histogram viewer"), "true"));
      await session.step(124, "And user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(125, "Then bar chart viewer should be absent", () => shouldBe(page, el("bar chart viewer"), "absent"));
      await session.step(126, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
