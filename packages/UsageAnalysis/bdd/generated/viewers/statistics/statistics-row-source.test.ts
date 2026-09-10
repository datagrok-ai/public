/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/statistics/statistics-row-source.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.stats-viewer]
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
import {addCategoricalFilter, filterPasses, filterTo, openEmptyFilterPanel, resetFilter, selectNoRows, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, noErrors, propertyShouldBe, readingIs, readingReads, reportsNoError, resizeTo, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Which rows the statistics are computed over, and which columns get a row", () => {
  const session = feature(test, "features/viewers/statistics/statistics-row-source.feature", import.meta.url);
  test("Which rows the statistics are computed over, and which columns get a row", {tag: ["@journey", "@viewers", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(22, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(23, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(24, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(25, "And \"rowSource\" property of statistics viewer should be \"Filtered\"", () => propertyShouldBe(page, "rowSource", el("statistics viewer"), "Filtered"));
    await session.step(26, "And the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(27, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(28, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
    await run.scenario("With Row Source Filtered the table's filter moves every number", async () => {
      await session.step(31, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(32, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(33, "And the \"rows shown\" reading of statistics viewer should be 447", () => readingIs(page, "rows shown", el("statistics viewer"), 447));
      await session.step(34, "And the \"values of AGE\" reading of statistics viewer should be \"447\"", () => readingReads(page, "values of AGE", el("statistics viewer"), "447"));
      await session.step(35, "And the \"avg of AGE\" reading of statistics viewer should be \"44.55\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "44.55"));
      await session.step(36, "And the \"unique of SEX\" reading of statistics viewer should be \"1\"", () => readingReads(page, "unique of SEX", el("statistics viewer"), "1"));
      await session.step(37, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(38, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(39, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(40, "And the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
      await session.step(41, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
      await session.step(42, "And the \"unique of SEX\" reading of statistics viewer should be \"2\"", () => readingReads(page, "unique of SEX", el("statistics viewer"), "2"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Row Source Selected the selection is what is counted", async () => {
      await session.step(46, "When user sets \"rowSource\" property of statistics viewer to \"Selected\"", () => setProperty(page, "rowSource", el("statistics viewer"), "Selected"));
      await session.step(47, "And user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
      await session.step(48, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(49, "And the \"rows shown\" reading of statistics viewer should be 553", () => readingIs(page, "rows shown", el("statistics viewer"), 553));
      await session.step(50, "And the \"values of AGE\" reading of statistics viewer should be \"553\"", () => readingReads(page, "values of AGE", el("statistics viewer"), "553"));
      await session.step(51, "And the \"avg of AGE\" reading of statistics viewer should be \"46.59\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "46.59"));
      await session.step(52, "And the \"max of AGE\" reading of statistics viewer should be \"80.00\"", () => readingReads(page, "max of AGE", el("statistics viewer"), "80.00"));
      await session.step(53, "When user selects rows where \"SEX\" is \"M\"", () => selectWhereIs(page, "SEX", "M"));
      await session.step(54, "Then the \"rows shown\" reading of statistics viewer should be 447", () => readingIs(page, "rows shown", el("statistics viewer"), 447));
      await session.step(55, "And the \"avg of AGE\" reading of statistics viewer should be \"44.55\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "44.55"));
      await session.step(56, "And the \"max of AGE\" reading of statistics viewer should be \"89.00\"", () => readingReads(page, "max of AGE", el("statistics viewer"), "89.00"));
      await session.step(57, "When user selects no rows", () => selectNoRows(page));
      await session.step(58, "And user sets \"rowSource\" property of statistics viewer to \"All\"", () => setProperty(page, "rowSource", el("statistics viewer"), "All"));
      await session.step(59, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
      await session.step(60, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
      await session.step(61, "When user sets \"rowSource\" property of statistics viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("statistics viewer"), "Filtered"));
      await session.step(62, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Row Source All the table's filter leaves the statistics where they are", async () => {
      await session.step(65, "When user sets \"rowSource\" property of statistics viewer to \"All\"", () => setProperty(page, "rowSource", el("statistics viewer"), "All"));
      await session.step(66, "And user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
      await session.step(67, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(68, "And the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
      await session.step(69, "And the \"values of AGE\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "values of AGE", el("statistics viewer"), "1000"));
      await session.step(70, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
      await session.step(71, "When user sets \"rowSource\" property of statistics viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("statistics viewer"), "Filtered"));
      await session.step(72, "Then the \"rows shown\" reading of statistics viewer should be 447", () => readingIs(page, "rows shown", el("statistics viewer"), 447));
      await session.step(73, "And the \"avg of AGE\" reading of statistics viewer should be \"44.55\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "44.55"));
      await session.step(74, "When user resets the filter", () => resetFilter(page));
      await session.step(75, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
      await session.step(76, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Columns decides which of the table's columns gets a row", async () => {
      await session.step(80, "Then the \"columns\" reading of statistics viewer should be 11", () => readingIs(page, "columns", el("statistics viewer"), 11));
      await session.step(81, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
      await session.step(82, "And statistics viewer should have a \"row HEIGHT\" area", () => hasArea(page, el("statistics viewer"), "row HEIGHT"));
      await session.step(83, "When user sets \"columnNames\" property of statistics viewer to \"USUBJID, AGE, SEX, RACE, DIS_POP, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => setProperty(page, "columnNames", el("statistics viewer"), "USUBJID, AGE, SEX, RACE, DIS_POP, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(84, "Then the \"columns shown\" reading of statistics viewer should be 10", () => readingIs(page, "columns shown", el("statistics viewer"), 10));
      await session.step(85, "And the \"columns\" reading of statistics viewer should be 11", () => readingIs(page, "columns", el("statistics viewer"), 11));
      await session.step(86, "And statistics viewer should not have a \"row HEIGHT\" area", () => hasNoArea(page, el("statistics viewer"), "row HEIGHT"));
      await session.step(87, "And statistics viewer should have a \"row WEIGHT\" area", () => hasArea(page, el("statistics viewer"), "row WEIGHT"));
      await session.step(88, "When user sets \"columnNames\" property of statistics viewer to \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => setProperty(page, "columnNames", el("statistics viewer"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(89, "Then the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
      await session.step(90, "And statistics viewer should have a \"row HEIGHT\" area", () => hasArea(page, el("statistics viewer"), "row HEIGHT"));
      await session.step(91, "And the \"values of HEIGHT\" reading of statistics viewer should be \"872\"", () => readingReads(page, "values of HEIGHT", el("statistics viewer"), "872"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
