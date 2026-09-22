/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-row-source.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pivot-table]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {aggregationMatches, filteredAggregationMatches} from '../../../bindings/pivot-table.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {addCategoricalFilter, filterIsExactlyCategory, filterPasses, filterPassesAll, noneOfFiltered, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, noErrors, propertyShouldBe, readingIs, readingReads, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — what Row Source aggregates", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-row-source.feature", import.meta.url);
  test("Pivot table — what Row Source aggregates", {tag: ["@journey", "@viewers", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(25, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(26, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), ""));
    await session.step(27, "And user sets \"Row Source\" property of pivot table viewer to \"All\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "All"));
    await session.step(28, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(29, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
    await run.scenario("At Row Source All the source filter changes nothing", async () => {
      await session.step(32, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(33, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(34, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(35, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(36, "And the \"text of grid cell 5 of avg(AGE)\" reading of pivot table viewer should be \"51.60\"", () => readingReads(page, "text of grid cell 5 of avg(AGE)", el("pivot table viewer"), "51.60"));
      await session.step(37, "And the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"DIS_POP\"", () => aggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "DIS_POP"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("At Row Source Filtered the pivot re-aggregates over the filtered rows", async () => {
      await session.step(41, "When user sets \"Row Source\" property of pivot table viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "Filtered"));
      await session.step(42, "Then the \"rows shown\" reading of pivot table viewer should be 447", () => readingIs(page, "rows shown", el("pivot table viewer"), 447));
      await session.step(43, "And pivot table viewer should show 447 rows", () => showsRows(page, el("pivot table viewer"), 447));
      await session.step(44, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(45, "And the \"text of grid cell 5 of avg(AGE)\" reading of pivot table viewer should be \"52.10\"", () => readingReads(page, "text of grid cell 5 of avg(AGE)", el("pivot table viewer"), "52.10"));
      await session.step(46, "And the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"DIS_POP\" over the filtered rows", () => filteredAggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "DIS_POP"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Lifting the filter takes the pivot back to the whole table", async () => {
      await session.step(50, "When user adds a categorical filter on \"SEX\" keeping \"F, M\"", () => addCategoricalFilter(page, "SEX", "F, M"));
      await session.step(51, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(52, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(53, "And the \"text of grid cell 5 of avg(AGE)\" reading of pivot table viewer should be \"51.60\"", () => readingReads(page, "text of grid cell 5 of avg(AGE)", el("pivot table viewer"), "51.60"));
      await session.step(54, "And the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"DIS_POP\"", () => aggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "DIS_POP"));
      await session.step(55, "When user sets \"Row Source\" property of pivot table viewer to \"All\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "All"));
      await session.step(56, "Then the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("At Row Source All a key cell click takes the whole group, and a panel card stacks on it", async () => {
      await session.step(60, "Given \"Row Source\" property of pivot table viewer should be \"All\"", () => propertyShouldBe(page, "Row Source", el("pivot table viewer"), "All"));
      await session.step(61, "And \"Filtering Enabled\" property of pivot table viewer should be \"true\"", () => propertyShouldBe(page, "Filtering Enabled", el("pivot table viewer"), "true"));
      await session.step(62, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(63, "When user clicks on the \"grid cell 5 of DIS_POP\" area of pivot table viewer", () => clickArea(page, "grid cell 5 of DIS_POP", el("pivot table viewer")));
      await session.step(64, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(65, "And the filter should pass exactly the rows where \"DIS_POP\" is \"RA\"", () => filterIsExactlyCategory(page, "DIS_POP", "RA"));
      await session.step(66, "And the \"filter label\" reading of pivot table viewer should be \"DIS_POP in [RA]\"", () => readingReads(page, "filter label", el("pivot table viewer"), "DIS_POP in [RA]"));
      await session.step(67, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(68, "Then 104 rows should pass the filter", () => filterPasses(page, 104));
      await session.step(69, "And no rows where \"SEX\" is \"F\" should pass the filter", () => noneOfFiltered(page, "SEX", "F"));
      await session.step(70, "And no rows where \"DIS_POP\" is \"UC\" should pass the filter", () => noneOfFiltered(page, "DIS_POP", "UC"));
      await session.step(71, "And the \"filter label\" reading of pivot table viewer should be \"DIS_POP in [RA]\"", () => readingReads(page, "filter label", el("pivot table viewer"), "DIS_POP in [RA]"));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell under a pivot column takes the group crossed with that pivot value", async () => {
      await session.step(75, "When user adds a categorical filter on \"SEX\" keeping \"F, M\"", () => addCategoricalFilter(page, "SEX", "F, M"));
      await session.step(76, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(77, "And user sets \"Row Source\" property of pivot table viewer to \"All\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "All"));
      await session.step(78, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(79, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(80, "And the \"text of grid cell 5 of DIS_POP\" reading of pivot table viewer should be \"RA\"", () => readingReads(page, "text of grid cell 5 of DIS_POP", el("pivot table viewer"), "RA"));
      await session.step(81, "When user clicks on the \"grid cell 5 of None avg(AGE)\" area of pivot table viewer", () => clickArea(page, "grid cell 5 of None avg(AGE)", el("pivot table viewer")));
      await session.step(82, "Then 254 rows should pass the filter", () => filterPasses(page, 254));
      await session.step(83, "And no rows where \"DIS_POP\" is \"UC\" should pass the filter", () => noneOfFiltered(page, "DIS_POP", "UC"));
      await session.step(84, "And no rows where \"SEVERITY\" is \"Low\" should pass the filter", () => noneOfFiltered(page, "SEVERITY", "Low"));
      await session.step(85, "And the \"filter label\" reading of pivot table viewer should be \"DIS_POP in [RA], None avg(AGE) = 52.30\"", () => readingReads(page, "filter label", el("pivot table viewer"), "DIS_POP in [RA], None avg(AGE) = 52.30"));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("At Row Source Filtered the same click writes no filter at all", async () => {
      await session.step(89, "When user sets \"Row Source\" property of pivot table viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "Filtered"));
      await session.step(90, "And user resets the filter", () => resetFilter(page));
      await session.step(91, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(92, "When user clicks on the \"grid cell 5 of DIS_POP\" area of pivot table viewer", () => clickArea(page, "grid cell 5 of DIS_POP", el("pivot table viewer")));
      await session.step(93, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
