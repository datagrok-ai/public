/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-row-source.feature
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
import {aggregationMatches, filteredAggregationMatches} from '../../../bindings/pivot-table.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {addCategoricalFilter, filterPasses, filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, readingIs, readingReads, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — what Row Source aggregates", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-row-source.feature", import.meta.url);
  test("Pivot table — what Row Source aggregates", {tag: ["@journey", "@viewers", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(14, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), ""));
    await session.step(15, "And user sets \"Row Source\" property of pivot table viewer to \"All\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "All"));
    await session.step(16, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(17, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
    await run.scenario("At Row Source All the source filter changes nothing", async () => {
      await session.step(20, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(21, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(22, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(23, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(24, "And the \"text of grid cell 5 of avg(AGE)\" reading of pivot table viewer should be \"51.60\"", () => readingReads(page, "text of grid cell 5 of avg(AGE)", el("pivot table viewer"), "51.60"));
      await session.step(25, "And the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"DIS_POP\"", () => aggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "DIS_POP"));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("At Row Source Filtered the pivot re-aggregates over the filtered rows", async () => {
      await session.step(29, "When user sets \"Row Source\" property of pivot table viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "Filtered"));
      await session.step(30, "Then the \"rows shown\" reading of pivot table viewer should be 447", () => readingIs(page, "rows shown", el("pivot table viewer"), 447));
      await session.step(31, "And pivot table viewer should show 447 rows", () => showsRows(page, el("pivot table viewer"), 447));
      await session.step(32, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(33, "And the \"text of grid cell 5 of avg(AGE)\" reading of pivot table viewer should be \"52.10\"", () => readingReads(page, "text of grid cell 5 of avg(AGE)", el("pivot table viewer"), "52.10"));
      await session.step(34, "And the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"DIS_POP\" over the filtered rows", () => filteredAggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "DIS_POP"));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Lifting the filter takes the pivot back to the whole table", async () => {
      await session.step(38, "When user adds a categorical filter on \"SEX\" keeping \"F, M\"", () => addCategoricalFilter(page, "SEX", "F, M"));
      await session.step(39, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(40, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(41, "And the \"text of grid cell 5 of avg(AGE)\" reading of pivot table viewer should be \"51.60\"", () => readingReads(page, "text of grid cell 5 of avg(AGE)", el("pivot table viewer"), "51.60"));
      await session.step(42, "And the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"DIS_POP\"", () => aggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "DIS_POP"));
      await session.step(43, "When user sets \"Row Source\" property of pivot table viewer to \"All\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "All"));
      await session.step(44, "Then the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
