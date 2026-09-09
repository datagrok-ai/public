/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-links.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pivot-table]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {allOfFiltered, allOfSelected, filterPasses, filterPassesAll, noneOfFiltered, noneOfSelected, noneSelected, onlyOfSelected, resetFilter, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, clickAreaHolding, noErrors, propertyShouldBe, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — the filter and the selection a click sends back", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-links.feature", import.meta.url);
  test("Pivot table — the filter and the selection a click sends back", {tag: ["@journey", "@viewers", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(17, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), ""));
    await session.step(18, "And user sets \"Row Source\" property of pivot table viewer to \"All\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "All"));
    await session.step(19, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(20, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
    await session.step(21, "And the \"text of grid cell 5 of DIS_POP\" reading of pivot table viewer should be \"RA\"", () => readingReads(page, "text of grid cell 5 of DIS_POP", el("pivot table viewer"), "RA"));
    await run.scenario("Selecting a row of the aggregated grid selects the group's source rows", async () => {
      await session.step(24, "Given no rows should be selected", () => noneSelected(page));
      await session.step(25, "When user clicks on the \"grid row header 5\" area of pivot table viewer holding Control", () => clickAreaHolding(page, "grid row header 5", el("pivot table viewer"), "Control"));
      await session.step(26, "Then 434 rows should be selected", () => selectedRowCount(page, 434));
      await session.step(27, "And only rows where \"DIS_POP\" is \"RA\" should be selected", () => onlyOfSelected(page, "DIS_POP", "RA"));
      await session.step(28, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(29, "When user clicks on the \"grid row header 5\" area of pivot table viewer holding Control", () => clickAreaHolding(page, "grid row header 5", el("pivot table viewer"), "Control"));
      await session.step(30, "Then no rows should be selected", () => noneSelected(page));
      await session.step(31, "And the \"text of grid cell 6 of DIS_POP\" reading of pivot table viewer should be \"UC\"", () => readingReads(page, "text of grid cell 6 of DIS_POP", el("pivot table viewer"), "UC"));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second selected row adds its group to the selection", async () => {
      await session.step(35, "When user clicks on the \"grid row header 4\" area of pivot table viewer holding Control", () => clickAreaHolding(page, "grid row header 4", el("pivot table viewer"), "Control"));
      await session.step(36, "Then 204 rows should be selected", () => selectedRowCount(page, 204));
      await session.step(37, "And only rows where \"DIS_POP\" is \"Psoriasis\" should be selected", () => onlyOfSelected(page, "DIS_POP", "Psoriasis"));
      await session.step(38, "When user clicks on the \"grid row header 5\" area of pivot table viewer holding Control", () => clickAreaHolding(page, "grid row header 5", el("pivot table viewer"), "Control"));
      await session.step(39, "Then 638 rows should be selected", () => selectedRowCount(page, 638));
      await session.step(40, "And all rows where \"DIS_POP\" is \"RA\" should be selected", () => allOfSelected(page, "DIS_POP", "RA"));
      await session.step(41, "And all rows where \"DIS_POP\" is \"Psoriasis\" should be selected", () => allOfSelected(page, "DIS_POP", "Psoriasis"));
      await session.step(42, "And no rows where \"DIS_POP\" is \"UC\" should be selected", () => noneOfSelected(page, "DIS_POP", "UC"));
      await session.step(43, "When user clicks on the \"grid row header 4\" area of pivot table viewer holding Control", () => clickAreaHolding(page, "grid row header 4", el("pivot table viewer"), "Control"));
      await session.step(44, "And user clicks on the \"grid row header 5\" area of pivot table viewer holding Control", () => clickAreaHolding(page, "grid row header 5", el("pivot table viewer"), "Control"));
      await session.step(45, "Then no rows should be selected", () => noneSelected(page));
      await session.step(46, "And the \"text of grid cell 6 of DIS_POP\" reading of pivot table viewer should be \"UC\"", () => readingReads(page, "text of grid cell 6 of DIS_POP", el("pivot table viewer"), "UC"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell click filters the source table down to the clicked group", async () => {
      await session.step(50, "Then \"Filtering Enabled\" property of pivot table viewer should be \"true\"", () => propertyShouldBe(page, "Filtering Enabled", el("pivot table viewer"), "true"));
      await session.step(51, "And the \"filter label\" reading of pivot table viewer should be \"\"", () => readingReads(page, "filter label", el("pivot table viewer"), ""));
      await session.step(52, "When user clicks on the \"grid cell 5 of DIS_POP\" area of pivot table viewer", () => clickArea(page, "grid cell 5 of DIS_POP", el("pivot table viewer")));
      await session.step(53, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(54, "And all rows where \"DIS_POP\" is \"RA\" should pass the filter", () => allOfFiltered(page, "DIS_POP", "RA"));
      await session.step(55, "And no rows where \"DIS_POP\" is \"UC\" should pass the filter", () => noneOfFiltered(page, "DIS_POP", "UC"));
      await session.step(56, "And the \"filter label\" reading of pivot table viewer should be \"DIS_POP in [RA]\"", () => readingReads(page, "filter label", el("pivot table viewer"), "DIS_POP in [RA]"));
      await session.step(57, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(58, "When user sets \"Filtering Enabled\" property of pivot table viewer to \"false\"", () => setProperty(page, "Filtering Enabled", el("pivot table viewer"), "false"));
      await session.step(59, "And user resets the filter", () => resetFilter(page));
      await session.step(60, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(61, "And the \"text of grid cell 6 of DIS_POP\" reading of pivot table viewer should be \"UC\"", () => readingReads(page, "text of grid cell 6 of DIS_POP", el("pivot table viewer"), "UC"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The next click replaces the filter instead of narrowing it further", async () => {
      await session.step(65, "When user sets \"Filtering Enabled\" property of pivot table viewer to \"true\"", () => setProperty(page, "Filtering Enabled", el("pivot table viewer"), "true"));
      await session.step(66, "And user clicks on the \"grid cell 5 of DIS_POP\" area of pivot table viewer", () => clickArea(page, "grid cell 5 of DIS_POP", el("pivot table viewer")));
      await session.step(67, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(68, "When user clicks on the \"grid cell 4 of DIS_POP\" area of pivot table viewer", () => clickArea(page, "grid cell 4 of DIS_POP", el("pivot table viewer")));
      await session.step(69, "Then 204 rows should pass the filter", () => filterPasses(page, 204));
      await session.step(70, "And all rows where \"DIS_POP\" is \"Psoriasis\" should pass the filter", () => allOfFiltered(page, "DIS_POP", "Psoriasis"));
      await session.step(71, "And no rows where \"DIS_POP\" is \"RA\" should pass the filter", () => noneOfFiltered(page, "DIS_POP", "RA"));
      await session.step(72, "And the \"filter label\" reading of pivot table viewer should be \"DIS_POP in [Psoriasis]\"", () => readingReads(page, "filter label", el("pivot table viewer"), "DIS_POP in [Psoriasis]"));
      await session.step(73, "When user sets \"Filtering Enabled\" property of pivot table viewer to \"false\"", () => setProperty(page, "Filtering Enabled", el("pivot table viewer"), "false"));
      await session.step(74, "And user resets the filter", () => resetFilter(page));
      await session.step(75, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(76, "And the \"text of grid cell 6 of DIS_POP\" reading of pivot table viewer should be \"UC\"", () => readingReads(page, "text of grid cell 6 of DIS_POP", el("pivot table viewer"), "UC"));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filtering Enabled off a click leaves the filter alone", async () => {
      await session.step(80, "Then \"Filtering Enabled\" property of pivot table viewer should be \"false\"", () => propertyShouldBe(page, "Filtering Enabled", el("pivot table viewer"), "false"));
      await session.step(81, "When user clicks on the \"grid cell 5 of DIS_POP\" area of pivot table viewer", () => clickArea(page, "grid cell 5 of DIS_POP", el("pivot table viewer")));
      await session.step(82, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(83, "And the \"filter label\" reading of pivot table viewer should be \"\"", () => readingReads(page, "filter label", el("pivot table viewer"), ""));
      await session.step(84, "When user clicks on the \"grid cell 4 of DIS_POP\" area of pivot table viewer", () => clickArea(page, "grid cell 4 of DIS_POP", el("pivot table viewer")));
      await session.step(85, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(86, "And the \"text of grid cell 6 of DIS_POP\" reading of pivot table viewer should be \"UC\"", () => readingReads(page, "text of grid cell 6 of DIS_POP", el("pivot table viewer"), "UC"));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("At Row Source Filtered a click leaves the filter alone as well", async () => {
      await session.step(90, "When user sets \"Row Source\" property of pivot table viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "Filtered"));
      await session.step(91, "And user sets \"Filtering Enabled\" property of pivot table viewer to \"true\"", () => setProperty(page, "Filtering Enabled", el("pivot table viewer"), "true"));
      await session.step(92, "And user clicks on the \"grid cell 5 of DIS_POP\" area of pivot table viewer", () => clickArea(page, "grid cell 5 of DIS_POP", el("pivot table viewer")));
      await session.step(93, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(94, "And the \"filter label\" reading of pivot table viewer should be \"\"", () => readingReads(page, "filter label", el("pivot table viewer"), ""));
      await session.step(95, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(96, "When user clicks on the \"grid cell 2 of DIS_POP\" area of pivot table viewer", () => clickArea(page, "grid cell 2 of DIS_POP", el("pivot table viewer")));
      await session.step(97, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(98, "When user sets \"Row Source\" property of pivot table viewer to \"All\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "All"));
      await session.step(99, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
