/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/calendar/calendar.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.calendar]
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
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, filterPasses, selectFirstRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaColor, areaNotColor, areaTaller, hasArea, hasNoArea, hoverArea, noErrors, oneTooltip, openToolbox, painted, pointerAway, propertyShouldBe, readingAsRemembered, readingAtLeast, readingHigher, readingIs, readingReads, rememberReading, repainted, setProperties, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Calendar counting, chrome and tooltips", () => {
  const session = feature(test, "features/viewers/calendar/calendar.feature", import.meta.url);
  test("Calendar counting, chrome and tooltips", {tag: ["@journey", "@viewers", "@realizes:viewers.calendar"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(28, "And user adds a calendar viewer", () => addViewer(page, "calendar"));
    await session.step(29, "Then the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
    await session.step(30, "And calendar viewer should be painted", () => painted(page, el("calendar viewer")));
    await run.scenario("The date column is picked up by itself and every row is counted on a weekday", async () => {
      await session.step(33, "Then the \"date column\" reading of calendar viewer should be \"STARTED\"", () => readingReads(page, "date column", el("calendar viewer"), "STARTED"));
      await session.step(34, "And \"Date\" property of calendar viewer should be \"STARTED\"", () => propertyShouldBe(page, "Date", el("calendar viewer"), "STARTED"));
      await session.step(35, "And the \"days drawn\" reading of calendar viewer should be at least 1", () => readingAtLeast(page, "days drawn", el("calendar viewer"), 1));
      await session.step(36, "And the \"weeks\" reading of calendar viewer should be at least 1", () => readingAtLeast(page, "weeks", el("calendar viewer"), 1));
      await session.step(37, "And the \"rows of weekday Sunday\" reading of calendar viewer should be 153", () => readingIs(page, "rows of weekday Sunday", el("calendar viewer"), 153));
      await session.step(38, "And the \"rows of weekday Monday\" reading of calendar viewer should be 151", () => readingIs(page, "rows of weekday Monday", el("calendar viewer"), 151));
      await session.step(39, "And the \"rows of weekday Tuesday\" reading of calendar viewer should be 136", () => readingIs(page, "rows of weekday Tuesday", el("calendar viewer"), 136));
      await session.step(40, "And the \"rows of weekday Wednesday\" reading of calendar viewer should be 148", () => readingIs(page, "rows of weekday Wednesday", el("calendar viewer"), 148));
      await session.step(41, "And the \"rows of weekday Thursday\" reading of calendar viewer should be 154", () => readingIs(page, "rows of weekday Thursday", el("calendar viewer"), 154));
      await session.step(42, "And the \"rows of weekday Friday\" reading of calendar viewer should be 141", () => readingIs(page, "rows of weekday Friday", el("calendar viewer"), 141));
      await session.step(43, "And the \"rows of weekday Saturday\" reading of calendar viewer should be 117", () => readingIs(page, "rows of weekday Saturday", el("calendar viewer"), 117));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell is reported only when the frame drew it and something fell on it", async () => {
      await session.step(47, "Then calendar viewer should have a \"day 1989-12-21\" area", () => hasArea(page, el("calendar viewer"), "day 1989-12-21"));
      await session.step(48, "And the \"rows of day 1989-12-21\" reading of calendar viewer should be 5", () => readingIs(page, "rows of day 1989-12-21", el("calendar viewer"), 5));
      await session.step(49, "And calendar viewer should not have a \"day 1989-12-23\" area", () => hasNoArea(page, el("calendar viewer"), "day 1989-12-23"));
      await session.step(50, "And calendar viewer should have a \"month 1989-12\" area", () => hasArea(page, el("calendar viewer"), "month 1989-12"));
      await session.step(51, "And the \"rows of month 1989-12\" reading of calendar viewer should be 43", () => readingIs(page, "rows of month 1989-12", el("calendar viewer"), 43));
      await session.step(52, "And the \"rows of month 1990-01\" reading of calendar viewer should be 32", () => readingIs(page, "rows of month 1990-01", el("calendar viewer"), 32));
      await session.step(53, "And calendar viewer should have a \"busiest day\" area", () => hasArea(page, el("calendar viewer"), "busiest day"));
      await session.step(54, "And the \"rows in busiest day\" reading of calendar viewer should be at least 1", () => readingAtLeast(page, "rows in busiest day", el("calendar viewer"), 1));
      await session.step(55, "And calendar viewer should have a \"days\" area", () => hasArea(page, el("calendar viewer"), "days"));
      await session.step(56, "And calendar viewer should have a \"months\" area", () => hasArea(page, el("calendar viewer"), "months"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The busiest drawn day's tooltip names it and counts its rows", async () => {
      await session.step(60, "When user hovers over the \"busiest day\" area of calendar viewer", () => hoverArea(page, "busiest day", el("calendar viewer")));
      await session.step(61, "Then tooltip should contain text \"rows\"", () => shouldContainText(page, el("tooltip"), "rows"));
      await session.step(62, "And tooltip should contain text \"Click to select\"", () => shouldContainText(page, el("tooltip"), "Click to select"));
      await session.step(63, "And exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(64, "When user moves the pointer away from calendar viewer", () => pointerAway(page, el("calendar viewer")));
      await session.step(65, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Header takes the year caption away and gives the days its space", async () => {
      await session.step(68, "Then calendar viewer should have a \"header\" area", () => hasArea(page, el("calendar viewer"), "header"));
      await session.step(69, "When user remembers the \"weeks\" reading of calendar viewer", () => rememberReading(page, "weeks", el("calendar viewer")));
      await session.step(70, "And user sets \"showHeader\" property of calendar viewer to \"false\"", () => setProperty(page, "showHeader", el("calendar viewer"), "false"));
      await session.step(71, "Then calendar viewer should not have a \"header\" area", () => hasNoArea(page, el("calendar viewer"), "header"));
      await session.step(72, "And the \"weeks\" reading of calendar viewer should be higher than before", () => readingHigher(page, "weeks", el("calendar viewer")));
      await session.step(73, "And the \"days\" area of calendar viewer should be taller than before", () => areaTaller(page, "days", el("calendar viewer")));
      await session.step(74, "And the \"rows of weekday Sunday\" reading of calendar viewer should be 153", () => readingIs(page, "rows of weekday Sunday", el("calendar viewer"), 153));
      await session.step(75, "And calendar viewer should have repainted", () => repainted(page, el("calendar viewer")));
      await session.step(76, "When user sets \"showHeader\" property of calendar viewer to \"true\"", () => setProperty(page, "showHeader", el("calendar viewer"), "true"));
      await session.step(77, "Then calendar viewer should have a \"header\" area", () => hasArea(page, el("calendar viewer"), "header"));
      await session.step(78, "And the \"weeks\" reading of calendar viewer should be as remembered", () => readingAsRemembered(page, "weeks", el("calendar viewer")));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Red Weekends is the only red on the grid of days", async () => {
      await session.step(82, "Then the \"days\" area of calendar viewer should contain the color \"#FF0000\"", () => areaColor(page, "days", el("calendar viewer"), "#FF0000"));
      await session.step(83, "When user sets \"redWeekends\" property of calendar viewer to \"false\"", () => setProperty(page, "redWeekends", el("calendar viewer"), "false"));
      await session.step(84, "Then the \"days\" area of calendar viewer should not contain the color \"#FF0000\"", () => areaNotColor(page, "days", el("calendar viewer"), "#FF0000"));
      await session.step(85, "And calendar viewer should have repainted", () => repainted(page, el("calendar viewer")));
      await session.step(86, "When user sets \"redWeekends\" property of calendar viewer to \"true\"", () => setProperty(page, "redWeekends", el("calendar viewer"), "true"));
      await session.step(87, "Then the \"days\" area of calendar viewer should contain the color \"#FF0000\"", () => areaColor(page, "days", el("calendar viewer"), "#FF0000"));
      await session.step(88, "And calendar viewer should have repainted", () => repainted(page, el("calendar viewer")));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Odd and Even Month Color fill the day cells of alternate months", async () => {
      await session.step(92, "Then the \"days\" area of calendar viewer should not contain the color \"#00AA00\"", () => areaNotColor(page, "days", el("calendar viewer"), "#00AA00"));
      await session.step(93, "And the \"days\" area of calendar viewer should not contain the color \"#AA00AA\"", () => areaNotColor(page, "days", el("calendar viewer"), "#AA00AA"));
      await session.step(94, "When user sets \"oddMonthColor\" property of calendar viewer to \"#00AA00\"", () => setProperty(page, "oddMonthColor", el("calendar viewer"), "#00AA00"));
      await session.step(95, "Then the \"days\" area of calendar viewer should contain the color \"#00AA00\"", () => areaColor(page, "days", el("calendar viewer"), "#00AA00"));
      await session.step(96, "And the \"days\" area of calendar viewer should not contain the color \"#AA00AA\"", () => areaNotColor(page, "days", el("calendar viewer"), "#AA00AA"));
      await session.step(97, "When user sets \"evenMonthColor\" property of calendar viewer to \"#AA00AA\"", () => setProperty(page, "evenMonthColor", el("calendar viewer"), "#AA00AA"));
      await session.step(98, "Then the \"days\" area of calendar viewer should contain the color \"#AA00AA\"", () => areaColor(page, "days", el("calendar viewer"), "#AA00AA"));
      await session.step(99, "And the \"days\" area of calendar viewer should contain the color \"#00AA00\"", () => areaColor(page, "days", el("calendar viewer"), "#00AA00"));
      await session.step(100, "When user sets properties of calendar viewer:", () => setProperties(page, el("calendar viewer"), [["oddMonthColor","#FFFFFF"],["evenMonthColor","#F7F5F3"]]), [["oddMonthColor","#FFFFFF"],["evenMonthColor","#F7F5F3"]]);
      await session.step(103, "Then the \"days\" area of calendar viewer should not contain the color \"#00AA00\"", () => areaNotColor(page, "days", el("calendar viewer"), "#00AA00"));
      await session.step(104, "And the \"days\" area of calendar viewer should not contain the color \"#AA00AA\"", () => areaNotColor(page, "days", el("calendar viewer"), "#AA00AA"));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Back Color fills the header behind the year caption", async () => {
      await session.step(108, "Then the \"header\" area of calendar viewer should not contain the color \"#FFCC00\"", () => areaNotColor(page, "header", el("calendar viewer"), "#FFCC00"));
      await session.step(109, "When user sets \"backColor\" property of calendar viewer to \"#FFCC00\"", () => setProperty(page, "backColor", el("calendar viewer"), "#FFCC00"));
      await session.step(110, "Then the \"header\" area of calendar viewer should contain the color \"#FFCC00\"", () => areaColor(page, "header", el("calendar viewer"), "#FFCC00"));
      await session.step(111, "And the \"rows of weekday Sunday\" reading of calendar viewer should be 153", () => readingIs(page, "rows of weekday Sunday", el("calendar viewer"), 153));
      await session.step(112, "When user sets \"backColor\" property of calendar viewer to \"#FFFFFF\"", () => setProperty(page, "backColor", el("calendar viewer"), "#FFFFFF"));
      await session.step(113, "Then the \"header\" area of calendar viewer should not contain the color \"#FFCC00\"", () => areaNotColor(page, "header", el("calendar viewer"), "#FFCC00"));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source Selected counts the selected rows and nothing else", async () => {
      await session.step(117, "When user selects the first 20 rows", () => selectFirstRows(page, 20));
      await session.step(118, "Then 20 rows should be selected", () => selectedRowCount(page, 20));
      await session.step(119, "When user sets \"rowSource\" property of calendar viewer to \"Selected\"", () => setProperty(page, "rowSource", el("calendar viewer"), "Selected"));
      await session.step(120, "Then the \"rows shown\" reading of calendar viewer should be 20", () => readingIs(page, "rows shown", el("calendar viewer"), 20));
      await session.step(121, "And the \"rows of weekday Sunday\" reading of calendar viewer should be 2", () => readingIs(page, "rows of weekday Sunday", el("calendar viewer"), 2));
      await session.step(122, "And calendar viewer should have repainted", () => repainted(page, el("calendar viewer")));
      await session.step(123, "When user sets \"rowSource\" property of calendar viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("calendar viewer"), "Filtered"));
      await session.step(124, "And user clears the row selection", () => clearSelection(page));
      await session.step(125, "Then the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
      await session.step(126, "And the \"rows of weekday Sunday\" reading of calendar viewer should be 153", () => readingIs(page, "rows of weekday Sunday", el("calendar viewer"), 153));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer's own filter counts only the rows it keeps and leaves the table alone", async () => {
      await session.step(130, "When user sets \"filter\" property of calendar viewer to \"${AGE} > 40\"", () => setProperty(page, "filter", el("calendar viewer"), "${AGE} > 40"));
      await session.step(131, "Then the \"rows shown\" reading of calendar viewer should be 635", () => readingIs(page, "rows shown", el("calendar viewer"), 635));
      await session.step(132, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(133, "And the \"rows of weekday Sunday\" reading of calendar viewer should be 102", () => readingIs(page, "rows of weekday Sunday", el("calendar viewer"), 102));
      await session.step(134, "And calendar viewer should have repainted", () => repainted(page, el("calendar viewer")));
      await session.step(135, "When user sets \"filter\" property of calendar viewer to \"\"", () => setProperty(page, "filter", el("calendar viewer"), ""));
      await session.step(136, "Then the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
      await session.step(137, "And the \"rows of weekday Sunday\" reading of calendar viewer should be 153", () => readingIs(page, "rows of weekday Sunday", el("calendar viewer"), 153));
      await session.step(138, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Calendar icon on the Viewers pane of the toolbox adds the viewer", async () => {
      await session.step(141, "When user clicks on close icon of calendar viewer", () => clickOn(page, el("close icon of calendar viewer")));
      await session.step(142, "Then calendar viewer should be absent", () => shouldBe(page, el("calendar viewer"), "absent"));
      await session.step(143, "When user opens toolbox", () => openToolbox(page));
      await session.step(144, "And user clicks on calendar icon in toolbox", () => clickOn(page, el("calendar icon in toolbox")));
      await session.step(145, "Then calendar viewer should be visible", () => shouldBe(page, el("calendar viewer"), "visible"));
      await session.step(146, "And the open tableview should have 1 calendar viewer", () => viewerCount(page, 1, "calendar"));
      await session.step(147, "And the \"date column\" reading of calendar viewer should be \"STARTED\"", () => readingReads(page, "date column", el("calendar viewer"), "STARTED"));
      await session.step(148, "And the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
      await session.step(149, "And calendar viewer should be painted", () => painted(page, el("calendar viewer")));
      await session.step(150, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title bar closes the calendar", async () => {
      await session.step(153, "When user clicks on close icon of calendar viewer", () => clickOn(page, el("close icon of calendar viewer")));
      await session.step(154, "Then calendar viewer should be absent", () => shouldBe(page, el("calendar viewer"), "absent"));
      await session.step(155, "And the open tableview should have 0 calendar viewers", () => viewerCount(page, 0, "calendar"));
      await session.step(156, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
