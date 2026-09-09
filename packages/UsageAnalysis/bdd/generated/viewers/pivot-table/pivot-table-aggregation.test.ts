/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-aggregation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pivot-table]
--- */
import {test} from '@playwright/test';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {aggregationMatches, clearSavedParameters, pickFromHistory, pivotedAggregationMatches} from '../../../bindings/pivot-table.js';
import {readingContains, readingNotContains} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType, everyValueBetween, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {renameColumn, rowCount, tableColumnComplete, tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, openDataset, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, noErrors, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — the aggregation it publishes and the parameters it remembers", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-aggregation.feature", import.meta.url);
  test("Pivot table — the aggregation it publishes and the parameters it remembers", {tag: ["@journey", "@viewers", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user clears the saved pivot table parameters", () => clearSavedParameters(page));
    await session.step(18, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(19, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
    await session.step(20, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
    await session.step(21, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await run.scenario("ADD publishes the cross tab it drew", async () => {
      await session.step(24, "Then the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "DIS_POP", "SEVERITY"));
      await session.step(25, "When user clicks on the \"add to workspace\" area of pivot table viewer", () => clickArea(page, "add to workspace", el("pivot table viewer")));
      await session.step(26, "Then table \"demog-1000 aggregation\" should be open", () => tableOpen(page, "demog-1000 aggregation"));
      await session.step(27, "And table \"demog-1000 aggregation\" should have 6 rows", () => tableRows(page, "demog-1000 aggregation", 6));
      await session.step(28, "And table \"demog-1000 aggregation\" should have no missing values in \"DIS_POP\" column", () => tableColumnComplete(page, "demog-1000 aggregation", "DIS_POP"));
      await session.step(29, "And the \"demog-1000 aggregation\" view should be current", () => viewIsCurrent(page, "demog-1000 aggregation"));
      await session.step(30, "And the value of \"DIS_POP\" column in row 1 should be \"AS\"", () => valueInRow(page, "DIS_POP", 1, "AS"));
      await session.step(31, "And the value of \"DIS_POP\" column in row 5 should be \"RA\"", () => valueInRow(page, "DIS_POP", 5, "RA"));
      await session.step(32, "And table \"demog-1000 aggregation\" should have no missing values in \"None avg(AGE)\" column", () => tableColumnComplete(page, "demog-1000 aggregation", "None avg(AGE)"));
      await session.step(33, "And every value of \"None avg(AGE)\" column should lie between 36.91 and 52.31", () => everyValueBetween(page, "None avg(AGE)", 36.91, 52.31));
      await session.step(34, "When user closes the current view", () => closeCurrentView(page));
      await session.step(35, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(36, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The published key column keeps the type of the column it groups", async () => {
      await session.step(40, "When user sets \"Pivot Column Names\" property of pivot table viewer to \"\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), ""));
      await session.step(41, "Then the \"aggregated columns\" reading of pivot table viewer should be 2", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 2));
      await session.step(42, "When user clicks on the \"add to workspace\" area of pivot table viewer", () => clickArea(page, "add to workspace", el("pivot table viewer")));
      await session.step(43, "Then the \"demog-1000 aggregation\" view should be current", () => viewIsCurrent(page, "demog-1000 aggregation"));
      await session.step(44, "And \"DIS_POP\" column should have type \"string\"", () => columnType(page, "DIS_POP", "string"));
      await session.step(45, "And \"avg(AGE)\" column should have type \"double\"", () => columnType(page, "avg(AGE)", "double"));
      await session.step(46, "And the table should have 6 rows", () => rowCount(page, 6));
      await session.step(47, "When user closes the current view", () => closeCurrentView(page));
      await session.step(48, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(49, "Then the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"DIS_POP\"", () => aggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "DIS_POP"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Grouping by an identifier column makes one row per identifier", async () => {
      await session.step(53, "When user sets \"Group By Column Names\" property of pivot table viewer to \"USUBJID\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "USUBJID"));
      await session.step(54, "Then the \"aggregated rows\" reading of pivot table viewer should be 1000", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 1000));
      await session.step(55, "And the \"limited\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "limited", el("pivot table viewer"), "false"));
      await session.step(56, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(57, "And the \"aggregated columns\" reading of pivot table viewer should be 2", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 2));
      await session.step(58, "When user sets \"Group By Column Names\" property of pivot table viewer to \"DIS_POP\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "DIS_POP"));
      await session.step(59, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(60, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Save parameters, and the history menu puts the configuration back", async () => {
      await session.step(64, "When user sets \"Group By Column Names\" property of pivot table viewer to \"RACE\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "RACE"));
      await session.step(65, "And user sets \"Aggregate Column Names\" property of pivot table viewer to \"WEIGHT\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "WEIGHT"));
      await session.step(66, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(67, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), ""));
      await session.step(68, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(WEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(WEIGHT)"));
      await session.step(69, "And the \"aggregated rows\" reading of pivot table viewer should be 4", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 4));
      await session.step(70, "When user picks \"Save parameters\" from the history menu of pivot table viewer", () => pickFromHistory(page, "Save parameters"));
      await session.step(71, "Then the \"history entries\" reading of pivot table viewer should be \"key(RACE),avg(WEIGHT)\"", () => readingReads(page, "history entries", el("pivot table viewer"), "key(RACE),avg(WEIGHT)"));
      await session.step(72, "When user sets \"Group By Column Names\" property of pivot table viewer to \"DIS_POP\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "DIS_POP"));
      await session.step(73, "And user sets \"Aggregate Column Names\" property of pivot table viewer to \"AGE\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "AGE"));
      await session.step(74, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(75, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(76, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(77, "When user picks \"key(RACE),avg(WEIGHT)\" from the history menu of pivot table viewer", () => pickFromHistory(page, "key(RACE),avg(WEIGHT)"));
      await session.step(78, "Then the \"group by\" reading of pivot table viewer should be \"RACE\"", () => readingReads(page, "group by", el("pivot table viewer"), "RACE"));
      await session.step(79, "And the \"aggregate\" reading of pivot table viewer should be \"avg(WEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(WEIGHT)"));
      await session.step(80, "And the \"aggregated rows\" reading of pivot table viewer should be 4", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 4));
      await session.step(81, "And the \"text of grid cell 1 of RACE\" reading of pivot table viewer should be \"Asian\"", () => readingReads(page, "text of grid cell 1 of RACE", el("pivot table viewer"), "Asian"));
      await session.step(82, "And the \"text of grid cell 1 of avg(WEIGHT)\" reading of pivot table viewer should be \"70.52\"", () => readingReads(page, "text of grid cell 1 of avg(WEIGHT)", el("pivot table viewer"), "70.52"));
      await session.step(83, "And the aggregated values of pivot table viewer should match \"avg(WEIGHT)\" grouped by \"RACE\"", () => aggregationMatches(page, el("pivot table viewer"), "avg(WEIGHT)", "RACE"));
      await session.step(84, "When user clears the saved pivot table parameters", () => clearSavedParameters(page));
      await session.step(85, "And user sets \"Group By Column Names\" property of pivot table viewer to \"DIS_POP\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "DIS_POP"));
      await session.step(86, "And user sets \"Aggregate Column Names\" property of pivot table viewer to \"AGE\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "AGE"));
      await session.step(87, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(88, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(89, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved configuration is offered only while the table still has its columns", async () => {
      await session.step(93, "When user sets \"Group By Column Names\" property of pivot table viewer to \"RACE\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "RACE"));
      await session.step(94, "And user sets \"Aggregate Column Names\" property of pivot table viewer to \"WEIGHT\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "WEIGHT"));
      await session.step(95, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(96, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), ""));
      await session.step(97, "And user picks \"Save parameters\" from the history menu of pivot table viewer", () => pickFromHistory(page, "Save parameters"));
      await session.step(98, "Then the \"history entries\" reading of pivot table viewer should contain \"avg(WEIGHT)\"", () => readingContains(page, "history entries", el("pivot table viewer"), "avg(WEIGHT)"));
      await session.step(99, "When user sets \"Group By Column Names\" property of pivot table viewer to \"DIS_POP\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "DIS_POP"));
      await session.step(100, "And user sets \"Aggregate Column Names\" property of pivot table viewer to \"AGE\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "AGE"));
      await session.step(101, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(102, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(103, "And user renames \"WEIGHT\" column to \"MASS\"", () => renameColumn(page, "WEIGHT", "MASS"));
      await session.step(104, "And user clicks on close icon of pivot table viewer", () => clickOn(page, el("close icon of pivot table viewer")));
      await session.step(105, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
      await session.step(106, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(107, "And the \"history entries\" reading of pivot table viewer should not contain \"WEIGHT\"", () => readingNotContains(page, "history entries", el("pivot table viewer"), "WEIGHT"));
      await session.step(108, "When user renames \"MASS\" column to \"WEIGHT\"", () => renameColumn(page, "MASS", "WEIGHT"));
      await session.step(109, "And user clicks on close icon of pivot table viewer", () => clickOn(page, el("close icon of pivot table viewer")));
      await session.step(110, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
      await session.step(111, "Then the \"history entries\" reading of pivot table viewer should contain \"avg(WEIGHT)\"", () => readingContains(page, "history entries", el("pivot table viewer"), "avg(WEIGHT)"));
      await session.step(112, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(113, "When user clears the saved pivot table parameters", () => clearSavedParameters(page));
      await session.step(114, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
