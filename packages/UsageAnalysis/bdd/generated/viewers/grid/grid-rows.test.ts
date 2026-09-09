/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-rows.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {allRowsSelected, clearSelection, columnsSelected, filterBetween, filterPasses, filterPassesAll, filterTo, noColumnsSelected, noneSelected, onlyOfSelected, resetFilter, rowsRangeSelected, selectFirstRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, doubleClickArea, dragSelectionBetweenAreas, noErrors, pickFromAreaContextMenu, readingAsRemembered, readingDiffers, readingHigher, readingIs, readingReads, rememberReading, repainted, setProperty, showsRows, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid row selection, navigation and row source", () => {
  const session = feature(test, "features/viewers/grid/grid-rows.feature", import.meta.url);
  test("Grid row selection, navigation and row source", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(14, "And no rows should be selected", () => noneSelected(page));
    await run.scenario("A click sets the current row and selects nothing", async () => {
      await session.step(17, "When user clicks on the \"cell 6 of USUBJID\" area of grid", () => clickArea(page, "cell 6 of USUBJID", el("grid")));
      await session.step(18, "Then row 6 should be current", () => currentRowIs(page, 6));
      await session.step(19, "And the \"current row\" reading of grid should be 6", () => readingIs(page, "current row", el("grid"), 6));
      await session.step(20, "And the \"current column\" reading of grid should be \"USUBJID\"", () => readingReads(page, "current column", el("grid"), "USUBJID"));
      await session.step(21, "And no rows should be selected", () => noneSelected(page));
      await session.step(22, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shift-dragging the row-number strip selects a range", async () => {
      await session.step(25, "When user drags a selection box from the \"row header 6\" area to the \"row header 11\" area of grid", () => dragSelectionBetweenAreas(page, "row header 6", "row header 11", el("grid")));
      await session.step(26, "Then 6 rows should be selected", () => selectedRowCount(page, 6));
      await session.step(27, "And rows 6 to 11 should be selected", () => rowsRangeSelected(page, 6, 11));
      await session.step(28, "And no columns should be selected", () => noColumnsSelected(page));
      await session.step(29, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(30, "When user clears the row selection", () => clearSelection(page));
      await session.step(31, "Then no rows should be selected", () => noneSelected(page));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control+A selects every row and every column and Escape clears both", async () => {
      await session.step(35, "When user clicks on the \"cell 1 of AGE\" area of grid", () => clickArea(page, "cell 1 of AGE", el("grid")));
      await session.step(36, "And user presses Control+A", () => pressKey(page, "Control+A"));
      await session.step(37, "Then all rows should be selected", () => allRowsSelected(page));
      await session.step(38, "And columns \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\" should be selected", () => columnsSelected(page, "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(39, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(40, "Then no rows should be selected", () => noneSelected(page));
      await session.step(41, "And no columns should be selected", () => noColumnsSelected(page));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Six Control+clicks on headers select six columns without an error", async () => {
      await session.step(45, "When user clicks on the \"header AGE\" area of grid holding Control", () => clickAreaHolding(page, "header AGE", el("grid"), "Control"));
      await session.step(46, "And user clicks on the \"header HEIGHT\" area of grid holding Control", () => clickAreaHolding(page, "header HEIGHT", el("grid"), "Control"));
      await session.step(47, "And user clicks on the \"header WEIGHT\" area of grid holding Control", () => clickAreaHolding(page, "header WEIGHT", el("grid"), "Control"));
      await session.step(48, "And user clicks on the \"header SEX\" area of grid holding Control", () => clickAreaHolding(page, "header SEX", el("grid"), "Control"));
      await session.step(49, "And user clicks on the \"header RACE\" area of grid holding Control", () => clickAreaHolding(page, "header RACE", el("grid"), "Control"));
      await session.step(50, "And user clicks on the \"header DIS_POP\" area of grid holding Control", () => clickAreaHolding(page, "header DIS_POP", el("grid"), "Control"));
      await session.step(51, "Then columns \"AGE, HEIGHT, WEIGHT, SEX, RACE, DIS_POP\" should be selected", () => columnsSelected(page, "AGE, HEIGHT, WEIGHT, SEX, RACE, DIS_POP"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
      await session.step(53, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(54, "Then no columns should be selected", () => noColumnsSelected(page));
    });
    await run.scenario("Arrows and Control+Home, Control+End and PageDown move the current row", async () => {
      await session.step(57, "When user clicks on the \"cell 5 of AGE\" area of grid", () => clickArea(page, "cell 5 of AGE", el("grid")));
      await session.step(58, "And user presses ArrowDown", () => pressKey(page, "ArrowDown"));
      await session.step(59, "And user presses ArrowDown", () => pressKey(page, "ArrowDown"));
      await session.step(60, "And user presses ArrowRight", () => pressKey(page, "ArrowRight"));
      await session.step(61, "And user presses ArrowRight", () => pressKey(page, "ArrowRight"));
      await session.step(62, "Then row 7 should be current", () => currentRowIs(page, 7));
      await session.step(63, "And the \"current column\" reading of grid should be \"RACE\"", () => readingReads(page, "current column", el("grid"), "RACE"));
      await session.step(64, "When user presses Control+Home", () => pressKey(page, "Control+Home"));
      await session.step(65, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(66, "When user presses Control+End", () => pressKey(page, "Control+End"));
      await session.step(67, "Then row 1000 should be current", () => currentRowIs(page, 1000));
      await session.step(68, "When user presses Control+Home", () => pressKey(page, "Control+Home"));
      await session.step(69, "And user takes a snapshot of grid", () => takeSnapshot(page, el("grid")));
      await session.step(70, "And user presses PageDown", () => pressKey(page, "PageDown"));
      await session.step(71, "Then the \"current row\" reading of grid should be higher than before", () => readingHigher(page, "current row", el("grid")));
      await session.step(72, "When user takes a snapshot of grid", () => takeSnapshot(page, el("grid")));
      await session.step(73, "And user presses PageDown", () => pressKey(page, "PageDown"));
      await session.step(74, "Then the \"current row\" reading of grid should be higher than before", () => readingHigher(page, "current row", el("grid")));
      await session.step(75, "When user presses Control+Home", () => pressKey(page, "Control+Home"));
      await session.step(76, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Space selects the current row and Shift+Enter every row of the same value", async () => {
      await session.step(80, "When user clicks on the \"cell 4 of AGE\" area of grid", () => clickArea(page, "cell 4 of AGE", el("grid")));
      await session.step(81, "And user presses Space", () => pressKey(page, "Space"));
      await session.step(82, "Then 1 row should be selected", () => selectedRowCount(page, 1));
      await session.step(83, "And only rows where \"USUBJID\" is \"X0273T21000400002\" should be selected", () => onlyOfSelected(page, "USUBJID", "X0273T21000400002"));
      await session.step(84, "When user presses Space", () => pressKey(page, "Space"));
      await session.step(85, "Then no rows should be selected", () => noneSelected(page));
      await session.step(86, "When user clicks on the \"cell 4 of SEX\" area of grid", () => clickArea(page, "cell 4 of SEX", el("grid")));
      await session.step(87, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
      await session.step(88, "Then only rows where \"SEX\" is \"M\" should be selected", () => onlyOfSelected(page, "SEX", "M"));
      await session.step(89, "And 447 rows should be selected", () => selectedRowCount(page, 447));
      await session.step(90, "When user clears the row selection", () => clearSelection(page));
      await session.step(91, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Allow Row Selection off keeps the caret moving and the selection empty", async () => {
      await session.step(94, "When user sets \"Allow Row Selection\" property of grid to \"false\"", () => setProperty(page, "Allow Row Selection", el("grid"), "false"));
      await session.step(95, "And user clicks on the \"cell 5 of AGE\" area of grid", () => clickArea(page, "cell 5 of AGE", el("grid")));
      await session.step(96, "Then row 5 should be current", () => currentRowIs(page, 5));
      await session.step(97, "When user presses ArrowDown", () => pressKey(page, "ArrowDown"));
      await session.step(98, "Then row 6 should be current", () => currentRowIs(page, 6));
      await session.step(99, "When user presses Space", () => pressKey(page, "Space"));
      await session.step(100, "Then no rows should be selected", () => noneSelected(page));
      await session.step(101, "When user sets \"Allow Row Selection\" property of grid to \"true\"", () => setProperty(page, "Allow Row Selection", el("grid"), "true"));
      await session.step(102, "And user presses Space", () => pressKey(page, "Space"));
      await session.step(103, "Then 1 row should be selected", () => selectedRowCount(page, 1));
      await session.step(104, "When user clears the row selection", () => clearSelection(page));
      await session.step(105, "Then no rows should be selected", () => noneSelected(page));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source decides how many rows the grid shows", async () => {
      await session.step(109, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
      await session.step(110, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(111, "When user selects the first 5 rows", () => selectFirstRows(page, 5));
      await session.step(112, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
      await session.step(113, "When user sets \"Row Source\" property of grid to \"Filtered\"", () => setProperty(page, "Row Source", el("grid"), "Filtered"));
      await session.step(114, "Then grid should show 447 rows", () => showsRows(page, el("grid"), 447));
      await session.step(115, "When user sets \"Row Source\" property of grid to \"Selected\"", () => setProperty(page, "Row Source", el("grid"), "Selected"));
      await session.step(116, "Then grid should show 5 rows", () => showsRows(page, el("grid"), 5));
      await session.step(117, "When user sets \"Row Source\" property of grid to \"All\"", () => setProperty(page, "Row Source", el("grid"), "All"));
      await session.step(118, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
      await session.step(119, "And 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(120, "When user resets the filter", () => resetFilter(page));
      await session.step(121, "And user clears the row selection", () => clearSelection(page));
      await session.step(122, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(123, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter and a sort share one order", async () => {
      await session.step(126, "When user filters rows where \"AGE\" is between 31 and 999", () => filterBetween(page, "AGE", 31, 999));
      await session.step(127, "Then 844 rows should pass the filter", () => filterPasses(page, 844));
      await session.step(128, "When user sets \"Row Source\" property of grid to \"Filtered\"", () => setProperty(page, "Row Source", el("grid"), "Filtered"));
      await session.step(129, "Then grid should show 844 rows", () => showsRows(page, el("grid"), 844));
      await session.step(130, "When user remembers the \"row order\" reading of grid", () => rememberReading(page, "row order", el("grid")));
      await session.step(131, "And user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
      await session.step(132, "And user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
      await session.step(133, "Then the \"sort direction\" reading of grid should be \"ascending\"", () => readingReads(page, "sort direction", el("grid"), "ascending"));
      await session.step(134, "And grid should show 844 rows", () => showsRows(page, el("grid"), 844));
      await session.step(135, "And the \"row order\" reading of grid should differ from before", () => readingDiffers(page, "row order", el("grid")));
      await session.step(136, "When user picks \"Sort > Reset\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Sort > Reset", "header AGE", el("grid")));
      await session.step(137, "Then the \"row order\" reading of grid should be as remembered", () => readingAsRemembered(page, "row order", el("grid")));
      await session.step(138, "When user sets \"Row Source\" property of grid to \"All\"", () => setProperty(page, "Row Source", el("grid"), "All"));
      await session.step(139, "And user resets the filter", () => resetFilter(page));
      await session.step(140, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(141, "And grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
      await session.step(142, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
