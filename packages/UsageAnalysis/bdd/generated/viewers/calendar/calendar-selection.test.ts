/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/calendar/calendar-selection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.calendar]
--- */
import {test} from '@playwright/test';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/nx.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, clearSelection, filterPasses, noneSelected, resetFilter, selectedRowCount, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, clickAreaHolding, noErrors, painted, propertyShouldBe, readingHigher, readingIs, readingLower, readingReads, repainted, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Calendar selection, On Click and Show Filtered Only", () => {
  const session = feature(test, "features/viewers/calendar/calendar-selection.feature", import.meta.url);
  test("Calendar selection, On Click and Show Filtered Only", {tag: ["@journey", "@viewers", "@realizes:viewers.calendar"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a calendar viewer", () => addViewer(page, "calendar"));
    await session.step(19, "Then the \"date column\" reading of calendar viewer should be \"STARTED\"", () => readingReads(page, "date column", el("calendar viewer"), "STARTED"));
    await session.step(20, "And the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
    await session.step(21, "And \"On Click\" property of calendar viewer should be \"Select\"", () => propertyShouldBe(page, "On Click", el("calendar viewer"), "Select"));
    await session.step(22, "And calendar viewer should be painted", () => painted(page, el("calendar viewer")));
    await run.scenario("Clicking a day selects exactly the rows dated that day", async () => {
      await session.step(25, "Given user clears the row selection", () => clearSelection(page));
      await session.step(26, "When user clicks on the \"day 1989-12-21\" area of calendar viewer", () => clickArea(page, "day 1989-12-21", el("calendar viewer")));
      await session.step(27, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
      await session.step(28, "When user clears the row selection", () => clearSelection(page));
      await session.step(29, "And user clicks on the \"busiest day\" area of calendar viewer", () => clickArea(page, "busiest day", el("calendar viewer")));
      await session.step(30, "Then some rows should be selected", () => someSelected(page));
      await session.step(31, "When user clears the row selection", () => clearSelection(page));
      await session.step(32, "Then no rows should be selected", () => noneSelected(page));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking a weekday header selects every row that falls on it, Sunday included", async () => {
      await session.step(36, "Given user clears the row selection", () => clearSelection(page));
      await session.step(37, "When user clicks on the \"weekday Sunday\" area of calendar viewer", () => clickArea(page, "weekday Sunday", el("calendar viewer")));
      await session.step(38, "Then 153 rows should be selected", () => selectedRowCount(page, 153));
      await session.step(39, "When user clears the row selection", () => clearSelection(page));
      await session.step(40, "And user clicks on the \"weekday Saturday\" area of calendar viewer", () => clickArea(page, "weekday Saturday", el("calendar viewer")));
      await session.step(41, "Then 117 rows should be selected", () => selectedRowCount(page, 117));
      await session.step(42, "When user clears the row selection", () => clearSelection(page));
      await session.step(43, "And user clicks on the \"weekday Thursday\" area of calendar viewer", () => clickArea(page, "weekday Thursday", el("calendar viewer")));
      await session.step(44, "Then 154 rows should be selected", () => selectedRowCount(page, 154));
      await session.step(45, "When user clears the row selection", () => clearSelection(page));
      await session.step(46, "Then no rows should be selected", () => noneSelected(page));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shift extends the weekday selection and Ctrl toggles it back off", async () => {
      await session.step(50, "Given user clears the row selection", () => clearSelection(page));
      await session.step(51, "When user clicks on the \"weekday Monday\" area of calendar viewer", () => clickArea(page, "weekday Monday", el("calendar viewer")));
      await session.step(52, "Then 151 rows should be selected", () => selectedRowCount(page, 151));
      await session.step(53, "When user clicks on the \"weekday Tuesday\" area of calendar viewer holding Shift", () => clickAreaHolding(page, "weekday Tuesday", el("calendar viewer"), "Shift"));
      await session.step(54, "Then 287 rows should be selected", () => selectedRowCount(page, 287));
      await session.step(55, "When user clicks on the \"weekday Tuesday\" area of calendar viewer holding Control", () => clickAreaHolding(page, "weekday Tuesday", el("calendar viewer"), "Control"));
      await session.step(56, "Then 151 rows should be selected", () => selectedRowCount(page, 151));
      await session.step(57, "When user clears the row selection", () => clearSelection(page));
      await session.step(58, "Then no rows should be selected", () => noneSelected(page));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking a month band selects that month, and On Click = Filter filters it instead", async () => {
      await session.step(62, "Given user clears the row selection", () => clearSelection(page));
      await session.step(63, "When user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(64, "Then 32 rows should be selected", () => selectedRowCount(page, 32));
      await session.step(65, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(66, "When user clears the row selection", () => clearSelection(page));
      await session.step(67, "And user sets \"onClick\" property of calendar viewer to \"Filter\"", () => setProperty(page, "onClick", el("calendar viewer"), "Filter"));
      await session.step(68, "And user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(69, "Then 32 rows should pass the filter", () => filterPasses(page, 32));
      await session.step(70, "And no rows should be selected", () => noneSelected(page));
      await session.step(71, "And calendar viewer should have repainted", () => repainted(page, el("calendar viewer")));
      await session.step(72, "When user resets the filter", () => resetFilter(page));
      await session.step(73, "And user sets \"onClick\" property of calendar viewer to \"Select\"", () => setProperty(page, "onClick", el("calendar viewer"), "Select"));
      await session.step(74, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(75, "And \"On Click\" property of calendar viewer should be \"Select\"", () => propertyShouldBe(page, "On Click", el("calendar viewer"), "Select"));
      await session.step(76, "And the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Filtered Only decides whether the discs count the filtered rows or all of them", async () => {
      await session.step(80, "Then the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
      await session.step(81, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(82, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(83, "And the \"rows shown\" reading of calendar viewer should be 447", () => readingIs(page, "rows shown", el("calendar viewer"), 447));
      await session.step(84, "And the \"rows of weekday Sunday\" reading of calendar viewer should be lower than before", () => readingLower(page, "rows of weekday Sunday", el("calendar viewer")));
      await session.step(85, "And the \"days drawn\" reading of calendar viewer should be lower than before", () => readingLower(page, "days drawn", el("calendar viewer")));
      await session.step(86, "And calendar viewer should have repainted", () => repainted(page, el("calendar viewer")));
      await session.step(87, "When user sets \"showFilteredOnly\" property of calendar viewer to \"false\"", () => setProperty(page, "showFilteredOnly", el("calendar viewer"), "false"));
      await session.step(88, "Then the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
      await session.step(89, "And the \"rows of weekday Sunday\" reading of calendar viewer should be 153", () => readingIs(page, "rows of weekday Sunday", el("calendar viewer"), 153));
      await session.step(90, "And the \"days drawn\" reading of calendar viewer should be higher than before", () => readingHigher(page, "days drawn", el("calendar viewer")));
      await session.step(91, "And 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(92, "And calendar viewer should have repainted", () => repainted(page, el("calendar viewer")));
      await session.step(93, "When user sets \"showFilteredOnly\" property of calendar viewer to \"true\"", () => setProperty(page, "showFilteredOnly", el("calendar viewer"), "true"));
      await session.step(94, "Then the \"rows shown\" reading of calendar viewer should be 447", () => readingIs(page, "rows shown", el("calendar viewer"), 447));
      await session.step(95, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(96, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(97, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(98, "And the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
      await session.step(99, "And the \"rows of weekday Sunday\" reading of calendar viewer should be 153", () => readingIs(page, "rows of weekday Sunday", el("calendar viewer"), 153));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
