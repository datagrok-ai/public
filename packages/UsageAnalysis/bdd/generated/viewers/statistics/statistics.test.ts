/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/statistics/statistics.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.stats-viewer, entities.viewer.action.close-viewer]
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
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaPainted, hasArea, hasNoArea, noErrors, readingIs, readingReads, reportsNoError, resizeTo, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The statistics table — one row per column, one column per aggregation", () => {
  const session = feature(test, "features/viewers/statistics/statistics.feature", import.meta.url);
  test("The statistics table — one row per column, one column per aggregation", {tag: ["@journey", "@viewers", "@realizes:viewers.stats-viewer", "@realizes:entities.viewer.action.close-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(20, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(21, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(22, "And the \"stats\" reading of statistics viewer should be \"values, nulls, unique, min, max, avg, med, stdev\"", () => readingReads(page, "stats", el("statistics viewer"), "values, nulls, unique, min, max, avg, med, stdev"));
    await session.step(23, "And the \"columns\" reading of statistics viewer should be 11", () => readingIs(page, "columns", el("statistics viewer"), 11));
    await session.step(24, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(25, "And the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(26, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
    await run.scenario("Every column of the table gets a row, and the count statistics are filled in for all of them", async () => {
      await session.step(29, "Then the \"values of USUBJID\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "values of USUBJID", el("statistics viewer"), "1000"));
      await session.step(30, "And the \"values of AGE\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "values of AGE", el("statistics viewer"), "1000"));
      await session.step(31, "And the \"values of SEX\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "values of SEX", el("statistics viewer"), "1000"));
      await session.step(32, "And the \"unique of USUBJID\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "unique of USUBJID", el("statistics viewer"), "1000"));
      await session.step(33, "And the \"unique of SEX\" reading of statistics viewer should be \"2\"", () => readingReads(page, "unique of SEX", el("statistics viewer"), "2"));
      await session.step(34, "And the \"unique of RACE\" reading of statistics viewer should be \"4\"", () => readingReads(page, "unique of RACE", el("statistics viewer"), "4"));
      await session.step(35, "And the \"unique of STARTED\" reading of statistics viewer should be \"541\"", () => readingReads(page, "unique of STARTED", el("statistics viewer"), "541"));
      await session.step(36, "And the \"nulls of AGE\" reading of statistics viewer should be \"0\"", () => readingReads(page, "nulls of AGE", el("statistics viewer"), "0"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A blank is counted out of values and into nulls", async () => {
      await session.step(40, "Then the \"values of HEIGHT\" reading of statistics viewer should be \"872\"", () => readingReads(page, "values of HEIGHT", el("statistics viewer"), "872"));
      await session.step(41, "And the \"nulls of HEIGHT\" reading of statistics viewer should be \"128\"", () => readingReads(page, "nulls of HEIGHT", el("statistics viewer"), "128"));
      await session.step(42, "And the \"values of WEIGHT\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "values of WEIGHT", el("statistics viewer"), "1000"));
      await session.step(43, "And the \"nulls of WEIGHT\" reading of statistics viewer should be \"0\"", () => readingReads(page, "nulls of WEIGHT", el("statistics viewer"), "0"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The numerical aggregations are computed for a numerical column", async () => {
      await session.step(47, "Then the \"min of AGE\" reading of statistics viewer should be \"18.00\"", () => readingReads(page, "min of AGE", el("statistics viewer"), "18.00"));
      await session.step(48, "And the \"max of AGE\" reading of statistics viewer should be \"89.00\"", () => readingReads(page, "max of AGE", el("statistics viewer"), "89.00"));
      await session.step(49, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
      await session.step(50, "And the \"med of AGE\" reading of statistics viewer should be \"45.00\"", () => readingReads(page, "med of AGE", el("statistics viewer"), "45.00"));
      await session.step(51, "And the \"stdev of AGE\" reading of statistics viewer should be \"13.45\"", () => readingReads(page, "stdev of AGE", el("statistics viewer"), "13.45"));
      await session.step(52, "And the \"min of HEIGHT\" reading of statistics viewer should be \"137.32\"", () => readingReads(page, "min of HEIGHT", el("statistics viewer"), "137.32"));
      await session.step(53, "And the \"max of HEIGHT\" reading of statistics viewer should be \"198.86\"", () => readingReads(page, "max of HEIGHT", el("statistics viewer"), "198.86"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("For a categorical column the numerical aggregations are empty and the counts are not", async () => {
      await session.step(57, "Then the \"values of SEX\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "values of SEX", el("statistics viewer"), "1000"));
      await session.step(58, "And the \"nulls of SEX\" reading of statistics viewer should be \"0\"", () => readingReads(page, "nulls of SEX", el("statistics viewer"), "0"));
      await session.step(59, "And the \"unique of SEX\" reading of statistics viewer should be \"2\"", () => readingReads(page, "unique of SEX", el("statistics viewer"), "2"));
      await session.step(60, "And the \"min of SEX\" reading of statistics viewer should be \"\"", () => readingReads(page, "min of SEX", el("statistics viewer"), ""));
      await session.step(61, "And the \"max of SEX\" reading of statistics viewer should be \"\"", () => readingReads(page, "max of SEX", el("statistics viewer"), ""));
      await session.step(62, "And the \"avg of SEX\" reading of statistics viewer should be \"\"", () => readingReads(page, "avg of SEX", el("statistics viewer"), ""));
      await session.step(63, "And the \"stdev of SEX\" reading of statistics viewer should be \"\"", () => readingReads(page, "stdev of SEX", el("statistics viewer"), ""));
      await session.step(64, "And the \"avg of RACE\" reading of statistics viewer should be \"\"", () => readingReads(page, "avg of RACE", el("statistics viewer"), ""));
      await session.step(65, "And the \"avg of DIS_POP\" reading of statistics viewer should be \"\"", () => readingReads(page, "avg of DIS_POP", el("statistics viewer"), ""));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An aggregation has a header and a cell per column; the grid's own service columns have neither", async () => {
      await session.step(69, "Then statistics viewer should have a \"header values\" area", () => hasArea(page, el("statistics viewer"), "header values"));
      await session.step(70, "And statistics viewer should have a \"header avg\" area", () => hasArea(page, el("statistics viewer"), "header avg"));
      await session.step(71, "And statistics viewer should have a \"header stdev\" area", () => hasArea(page, el("statistics viewer"), "header stdev"));
      await session.step(72, "And statistics viewer should have a \"cell avg of AGE\" area", () => hasArea(page, el("statistics viewer"), "cell avg of AGE"));
      await session.step(73, "And statistics viewer should have a \"cell avg of STARTED\" area", () => hasArea(page, el("statistics viewer"), "cell avg of STARTED"));
      await session.step(74, "And statistics viewer should have a \"row AGE\" area", () => hasArea(page, el("statistics viewer"), "row AGE"));
      await session.step(75, "And statistics viewer should have a \"row SEVERITY\" area", () => hasArea(page, el("statistics viewer"), "row SEVERITY"));
      await session.step(76, "And statistics viewer should not have a \"header name\" area", () => hasNoArea(page, el("statistics viewer"), "header name"));
      await session.step(77, "And statistics viewer should not have a \"cell name of AGE\" area", () => hasNoArea(page, el("statistics viewer"), "cell name of AGE"));
      await session.step(78, "And the \"row AGE\" area of statistics viewer should be painted", () => areaPainted(page, "row AGE", el("statistics viewer")));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title bar closes the statistics viewer, and adding it again computes the same table", async () => {
      await session.step(82, "When user clicks on close icon of statistics viewer", () => clickOn(page, el("close icon of statistics viewer")));
      await session.step(83, "Then statistics viewer should be absent", () => shouldBe(page, el("statistics viewer"), "absent"));
      await session.step(84, "And the open tableview should have 0 statistics viewers", () => viewerCount(page, 0, "statistics"));
      await session.step(85, "When user adds a statistics viewer", () => addViewer(page, "statistics"));
      await session.step(86, "Then the open tableview should have 1 statistics viewer", () => viewerCount(page, 1, "statistics"));
      await session.step(87, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
      await session.step(88, "And the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
      await session.step(89, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
