/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/statistics/statistics-date-columns.feature
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, noErrors, pickFromAreaContextMenu, readingDoesNotRead, readingIs, readingReads, readingsDiffer, reportsNoError, resizeTo, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Statistics for a date column, and the aggregations added for one from the submenu", () => {
  const session = feature(test, "features/viewers/statistics/statistics-date-columns.feature", import.meta.url);
  test("Statistics for a date column, and the aggregations added for one from the submenu", {tag: ["@journey", "@viewers", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(18, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(19, "Then the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(20, "And statistics viewer should have a \"row STARTED\" area", () => hasArea(page, el("statistics viewer"), "row STARTED"));
    await session.step(21, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
    await run.scenario("The count statistics are filled in for a date column", async () => {
      await session.step(24, "Then the \"values of STARTED\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "values of STARTED", el("statistics viewer"), "1000"));
      await session.step(25, "And the \"nulls of STARTED\" reading of statistics viewer should be \"0\"", () => readingReads(page, "nulls of STARTED", el("statistics viewer"), "0"));
      await session.step(26, "And the \"unique of STARTED\" reading of statistics viewer should be \"541\"", () => readingReads(page, "unique of STARTED", el("statistics viewer"), "541"));
      await session.step(27, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A date column gets the value aggregations too, and they are not empty", async () => {
      await session.step(30, "Then the \"min of STARTED\" reading of statistics viewer should not be \"\"", () => readingDoesNotRead(page, "min of STARTED", el("statistics viewer"), ""));
      await session.step(31, "And the \"max of STARTED\" reading of statistics viewer should not be \"\"", () => readingDoesNotRead(page, "max of STARTED", el("statistics viewer"), ""));
      await session.step(32, "And the \"avg of STARTED\" reading of statistics viewer should not be \"\"", () => readingDoesNotRead(page, "avg of STARTED", el("statistics viewer"), ""));
      await session.step(33, "And the \"min of USUBJID\" reading of statistics viewer should be \"\"", () => readingReads(page, "min of USUBJID", el("statistics viewer"), ""));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("min and max are added for a date column from the Statistics submenu", async () => {
      await session.step(37, "When user sets \"stats\" property of statistics viewer to \"values, nulls, unique\"", () => setProperty(page, "stats", el("statistics viewer"), "values, nulls, unique"));
      await session.step(38, "Then the \"stats\" reading of statistics viewer should be \"values, nulls, unique\"", () => readingReads(page, "stats", el("statistics viewer"), "values, nulls, unique"));
      await session.step(39, "And statistics viewer should not have a \"header min\" area", () => hasNoArea(page, el("statistics viewer"), "header min"));
      await session.step(40, "And statistics viewer should not have a \"header max\" area", () => hasNoArea(page, el("statistics viewer"), "header max"));
      await session.step(41, "When user picks \"Statistics > min\" from the context menu of the \"row STARTED\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Statistics > min", "row STARTED", el("statistics viewer")));
      await session.step(42, "Then the \"stats\" reading of statistics viewer should contain \"min\"", () => readingContains(page, "stats", el("statistics viewer"), "min"));
      await session.step(43, "And statistics viewer should have a \"header min\" area", () => hasArea(page, el("statistics viewer"), "header min"));
      await session.step(44, "And the \"min of STARTED\" reading of statistics viewer should not be \"\"", () => readingDoesNotRead(page, "min of STARTED", el("statistics viewer"), ""));
      await session.step(45, "When user picks \"Statistics > max\" from the context menu of the \"row STARTED\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Statistics > max", "row STARTED", el("statistics viewer")));
      await session.step(46, "Then the \"stats\" reading of statistics viewer should contain \"max\"", () => readingContains(page, "stats", el("statistics viewer"), "max"));
      await session.step(47, "And statistics viewer should have a \"header max\" area", () => hasArea(page, el("statistics viewer"), "header max"));
      await session.step(48, "And the \"max of STARTED\" reading of statistics viewer should not be \"\"", () => readingDoesNotRead(page, "max of STARTED", el("statistics viewer"), ""));
      await session.step(49, "And the \"min of STARTED\" and \"max of STARTED\" readings of statistics viewer should differ", () => readingsDiffer(page, "min of STARTED", "max of STARTED", el("statistics viewer")));
      await session.step(50, "When user sets \"stats\" property of statistics viewer to \"values, nulls, unique, min, max, avg, med, stdev\"", () => setProperty(page, "stats", el("statistics viewer"), "values, nulls, unique, min, max, avg, med, stdev"));
      await session.step(51, "Then the \"stats\" reading of statistics viewer should be \"values, nulls, unique, min, max, avg, med, stdev\"", () => readingReads(page, "stats", el("statistics viewer"), "values, nulls, unique, min, max, avg, med, stdev"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Column visibility — a column dropped from Columns loses its row and gets it back", async () => {
      await session.step(55, "Then statistics viewer should have a \"row HEIGHT\" area", () => hasArea(page, el("statistics viewer"), "row HEIGHT"));
      await session.step(56, "And the \"values of HEIGHT\" reading of statistics viewer should be \"872\"", () => readingReads(page, "values of HEIGHT", el("statistics viewer"), "872"));
      await session.step(57, "When user sets \"columnNames\" property of statistics viewer to \"USUBJID, AGE, SEX, RACE, DIS_POP, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => setProperty(page, "columnNames", el("statistics viewer"), "USUBJID, AGE, SEX, RACE, DIS_POP, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(58, "Then the \"columns shown\" reading of statistics viewer should be 10", () => readingIs(page, "columns shown", el("statistics viewer"), 10));
      await session.step(59, "And statistics viewer should not have a \"row HEIGHT\" area", () => hasNoArea(page, el("statistics viewer"), "row HEIGHT"));
      await session.step(60, "And statistics viewer should have a \"row STARTED\" area", () => hasArea(page, el("statistics viewer"), "row STARTED"));
      await session.step(61, "When user sets \"columnNames\" property of statistics viewer to \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => setProperty(page, "columnNames", el("statistics viewer"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(62, "Then the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
      await session.step(63, "And statistics viewer should have a \"row HEIGHT\" area", () => hasArea(page, el("statistics viewer"), "row HEIGHT"));
      await session.step(64, "And the \"values of HEIGHT\" reading of statistics viewer should be \"872\"", () => readingReads(page, "values of HEIGHT", el("statistics viewer"), "872"));
      await session.step(65, "And the \"nulls of HEIGHT\" reading of statistics viewer should be \"128\"", () => readingReads(page, "nulls of HEIGHT", el("statistics viewer"), "128"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A date statistic reads as the date the cell shows", async () => {
      await session.step(76, "Then the \"min of STARTED\" reading of statistics viewer should be \"12/3/1989\"", () => readingReads(page, "min of STARTED", el("statistics viewer"), "12/3/1989"));
      await session.step(77, "And the \"max of STARTED\" reading of statistics viewer should be \"11/30/1991\"", () => readingReads(page, "max of STARTED", el("statistics viewer"), "11/30/1991"));
    });
    run.finish();
  });
});
