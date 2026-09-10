/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/calendar/calendar.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.calendar]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaColor, areaNotColor, areaTaller, hasArea, hasNoArea, hoverArea, noErrors, oneTooltip, painted, pointerAway, propertyShouldBe, readingAsRemembered, readingAtLeast, readingHigher, readingIs, readingReads, rememberReading, repainted, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Calendar counting, chrome and tooltips", () => {
  const session = feature(test, "features/viewers/calendar/calendar.feature", import.meta.url);
  test("Calendar counting, chrome and tooltips", {tag: ["@journey", "@viewers", "@realizes:viewers.calendar"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "And user adds a calendar viewer", () => addViewer(page, "calendar"));
    await session.step(23, "Then the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
    await session.step(24, "And calendar viewer should be painted", () => painted(page, el("calendar viewer")));
    await run.scenario("The date column is picked up by itself and every row is counted on a weekday", async () => {
      await session.step(27, "Then the \"date column\" reading of calendar viewer should be \"STARTED\"", () => readingReads(page, "date column", el("calendar viewer"), "STARTED"));
      await session.step(28, "And \"Date\" property of calendar viewer should be \"STARTED\"", () => propertyShouldBe(page, "Date", el("calendar viewer"), "STARTED"));
      await session.step(29, "And the \"days drawn\" reading of calendar viewer should be at least 1", () => readingAtLeast(page, "days drawn", el("calendar viewer"), 1));
      await session.step(30, "And the \"weeks\" reading of calendar viewer should be at least 1", () => readingAtLeast(page, "weeks", el("calendar viewer"), 1));
      await session.step(31, "And the \"rows of weekday Sunday\" reading of calendar viewer should be 153", () => readingIs(page, "rows of weekday Sunday", el("calendar viewer"), 153));
      await session.step(32, "And the \"rows of weekday Monday\" reading of calendar viewer should be 151", () => readingIs(page, "rows of weekday Monday", el("calendar viewer"), 151));
      await session.step(33, "And the \"rows of weekday Tuesday\" reading of calendar viewer should be 136", () => readingIs(page, "rows of weekday Tuesday", el("calendar viewer"), 136));
      await session.step(34, "And the \"rows of weekday Wednesday\" reading of calendar viewer should be 148", () => readingIs(page, "rows of weekday Wednesday", el("calendar viewer"), 148));
      await session.step(35, "And the \"rows of weekday Thursday\" reading of calendar viewer should be 154", () => readingIs(page, "rows of weekday Thursday", el("calendar viewer"), 154));
      await session.step(36, "And the \"rows of weekday Friday\" reading of calendar viewer should be 141", () => readingIs(page, "rows of weekday Friday", el("calendar viewer"), 141));
      await session.step(37, "And the \"rows of weekday Saturday\" reading of calendar viewer should be 117", () => readingIs(page, "rows of weekday Saturday", el("calendar viewer"), 117));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell is reported only when the frame drew it and something fell on it", async () => {
      await session.step(41, "Then calendar viewer should have a \"day 1989-12-21\" area", () => hasArea(page, el("calendar viewer"), "day 1989-12-21"));
      await session.step(42, "And the \"rows of day 1989-12-21\" reading of calendar viewer should be 5", () => readingIs(page, "rows of day 1989-12-21", el("calendar viewer"), 5));
      await session.step(43, "And calendar viewer should not have a \"day 1989-12-23\" area", () => hasNoArea(page, el("calendar viewer"), "day 1989-12-23"));
      await session.step(44, "And calendar viewer should have a \"month 1989-12\" area", () => hasArea(page, el("calendar viewer"), "month 1989-12"));
      await session.step(45, "And the \"rows of month 1989-12\" reading of calendar viewer should be 43", () => readingIs(page, "rows of month 1989-12", el("calendar viewer"), 43));
      await session.step(46, "And the \"rows of month 1990-01\" reading of calendar viewer should be 32", () => readingIs(page, "rows of month 1990-01", el("calendar viewer"), 32));
      await session.step(47, "And calendar viewer should have a \"busiest day\" area", () => hasArea(page, el("calendar viewer"), "busiest day"));
      await session.step(48, "And the \"rows in busiest day\" reading of calendar viewer should be at least 1", () => readingAtLeast(page, "rows in busiest day", el("calendar viewer"), 1));
      await session.step(49, "And calendar viewer should have a \"days\" area", () => hasArea(page, el("calendar viewer"), "days"));
      await session.step(50, "And calendar viewer should have a \"months\" area", () => hasArea(page, el("calendar viewer"), "months"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The busiest drawn day's tooltip names it and counts its rows", async () => {
      await session.step(54, "When user hovers over the \"busiest day\" area of calendar viewer", () => hoverArea(page, "busiest day", el("calendar viewer")));
      await session.step(55, "Then tooltip should contain text \"rows\"", () => shouldContainText(page, el("tooltip"), "rows"));
      await session.step(56, "And tooltip should contain text \"Click to select\"", () => shouldContainText(page, el("tooltip"), "Click to select"));
      await session.step(57, "And exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(58, "When user moves the pointer away from calendar viewer", () => pointerAway(page, el("calendar viewer")));
      await session.step(59, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Header takes the year caption away and gives the days its space", async () => {
      await session.step(62, "Then calendar viewer should have a \"header\" area", () => hasArea(page, el("calendar viewer"), "header"));
      await session.step(63, "When user remembers the \"weeks\" reading of calendar viewer", () => rememberReading(page, "weeks", el("calendar viewer")));
      await session.step(64, "And user sets \"showHeader\" property of calendar viewer to \"false\"", () => setProperty(page, "showHeader", el("calendar viewer"), "false"));
      await session.step(65, "Then calendar viewer should not have a \"header\" area", () => hasNoArea(page, el("calendar viewer"), "header"));
      await session.step(66, "And the \"weeks\" reading of calendar viewer should be higher than before", () => readingHigher(page, "weeks", el("calendar viewer")));
      await session.step(67, "And the \"days\" area of calendar viewer should be taller than before", () => areaTaller(page, "days", el("calendar viewer")));
      await session.step(68, "And the \"rows of weekday Sunday\" reading of calendar viewer should be 153", () => readingIs(page, "rows of weekday Sunday", el("calendar viewer"), 153));
      await session.step(69, "And calendar viewer should have repainted", () => repainted(page, el("calendar viewer")));
      await session.step(70, "When user sets \"showHeader\" property of calendar viewer to \"true\"", () => setProperty(page, "showHeader", el("calendar viewer"), "true"));
      await session.step(71, "Then calendar viewer should have a \"header\" area", () => hasArea(page, el("calendar viewer"), "header"));
      await session.step(72, "And the \"weeks\" reading of calendar viewer should be as remembered", () => readingAsRemembered(page, "weeks", el("calendar viewer")));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Red Weekends is the only red on the grid of days", async () => {
      await session.step(76, "Then the \"days\" area of calendar viewer should contain the color \"#FF0000\"", () => areaColor(page, "days", el("calendar viewer"), "#FF0000"));
      await session.step(77, "When user sets \"redWeekends\" property of calendar viewer to \"false\"", () => setProperty(page, "redWeekends", el("calendar viewer"), "false"));
      await session.step(78, "Then the \"days\" area of calendar viewer should not contain the color \"#FF0000\"", () => areaNotColor(page, "days", el("calendar viewer"), "#FF0000"));
      await session.step(79, "And calendar viewer should have repainted", () => repainted(page, el("calendar viewer")));
      await session.step(80, "When user sets \"redWeekends\" property of calendar viewer to \"true\"", () => setProperty(page, "redWeekends", el("calendar viewer"), "true"));
      await session.step(81, "Then the \"days\" area of calendar viewer should contain the color \"#FF0000\"", () => areaColor(page, "days", el("calendar viewer"), "#FF0000"));
      await session.step(82, "And calendar viewer should have repainted", () => repainted(page, el("calendar viewer")));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title bar closes the calendar", async () => {
      await session.step(86, "When user clicks on close icon of calendar viewer", () => clickOn(page, el("close icon of calendar viewer")));
      await session.step(87, "Then calendar viewer should be absent", () => shouldBe(page, el("calendar viewer"), "absent"));
      await session.step(88, "And the open tableview should have 0 calendar viewers", () => viewerCount(page, 0, "calendar"));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
