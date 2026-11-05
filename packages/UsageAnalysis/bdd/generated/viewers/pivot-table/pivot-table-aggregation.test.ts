/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-aggregation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pivot-table, public#3470]
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
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnType, everyValueBetween, setColumnSemType, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {renameColumn, rowCount, setLinkKeyPair, setLinkTables, tableColumnComplete, tableFilterCount, tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, openDataset, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {aggregationMatches, clearSavedParameters, pickFromHistory, pivotedAggregationMatches} from '@datagrok-libraries/bdd/bindings/tiers/viewers/pivot-table';
import {addViewer, clickArea, noBalloons, noErrors, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains, readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — the aggregation it publishes and the parameters it remembers", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-aggregation.feature", import.meta.url);
  test("Pivot table — the aggregation it publishes and the parameters it remembers", {tag: ["@journey", "@viewers", "@realizes:viewers.pivot-table", "@realizes:public#3470"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user clears the saved pivot table parameters", () => clearSavedParameters(page));
    await session.step(20, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(21, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
    await session.step(22, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
    await session.step(23, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await run.scenario("ADD publishes the cross tab it drew", async () => {
      await session.step(26, "Then the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "DIS_POP", "SEVERITY"));
      await session.step(27, "When user clicks on the \"add to workspace\" area of pivot table viewer", () => clickArea(page, "add to workspace", el("pivot table viewer")));
      await session.step(28, "Then table \"demog-1000 aggregation\" should be open", () => tableOpen(page, "demog-1000 aggregation"));
      await session.step(29, "And table \"demog-1000 aggregation\" should have 6 rows", () => tableRows(page, "demog-1000 aggregation", 6));
      await session.step(30, "And table \"demog-1000 aggregation\" should have no missing values in \"DIS_POP\" column", () => tableColumnComplete(page, "demog-1000 aggregation", "DIS_POP"));
      await session.step(31, "And the \"demog-1000 aggregation\" view should be current", () => viewIsCurrent(page, "demog-1000 aggregation"));
      await session.step(32, "And the value of \"DIS_POP\" column in row 1 should be \"AS\"", () => valueInRow(page, "DIS_POP", 1, "AS"));
      await session.step(33, "And the value of \"DIS_POP\" column in row 5 should be \"RA\"", () => valueInRow(page, "DIS_POP", 5, "RA"));
      await session.step(34, "And table \"demog-1000 aggregation\" should have no missing values in \"None avg(AGE)\" column", () => tableColumnComplete(page, "demog-1000 aggregation", "None avg(AGE)"));
      await session.step(35, "And every value of \"None avg(AGE)\" column should lie between 36.91 and 52.31", () => everyValueBetween(page, "None avg(AGE)", 36.91, 52.31));
      await session.step(36, "When user closes the current view", () => closeCurrentView(page));
      await session.step(37, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(38, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The published key column keeps the type of the column it groups", async () => {
      await session.step(42, "When user sets \"Pivot Column Names\" property of pivot table viewer to \"\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), ""));
      await session.step(43, "Then the \"aggregated columns\" reading of pivot table viewer should be 2", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 2));
      await session.step(44, "When user sets the semantic type of \"DIS_POP\" column to \"Text\"", () => setColumnSemType(page, "DIS_POP", "Text"));
      await session.step(45, "Then \"DIS_POP\" column should have semantic type \"Text\"", () => columnSemType(page, "DIS_POP", "Text"));
      await session.step(46, "When user clicks on the \"add to workspace\" area of pivot table viewer", () => clickArea(page, "add to workspace", el("pivot table viewer")));
      await session.step(47, "Then the \"demog-1000 aggregation\" view should be current", () => viewIsCurrent(page, "demog-1000 aggregation"));
      await session.step(48, "And \"DIS_POP\" column should have type \"string\"", () => columnType(page, "DIS_POP", "string"));
      await session.step(49, "And \"DIS_POP\" column should have semantic type \"Text\"", () => columnSemType(page, "DIS_POP", "Text"));
      await session.step(50, "And \"avg(AGE)\" column should have type \"double\"", () => columnType(page, "avg(AGE)", "double"));
      await session.step(51, "And the table should have 6 rows", () => rowCount(page, 6));
      await session.step(52, "When user closes the current view", () => closeCurrentView(page));
      await session.step(53, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(54, "And user sets the semantic type of \"DIS_POP\" column to \"\"", () => setColumnSemType(page, "DIS_POP", ""));
      await session.step(55, "Then the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"DIS_POP\"", () => aggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "DIS_POP"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Grouping by an identifier column makes one row per identifier", async () => {
      await session.step(59, "When user sets \"Group By Column Names\" property of pivot table viewer to \"USUBJID\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "USUBJID"));
      await session.step(60, "Then the \"aggregated rows\" reading of pivot table viewer should be 1000", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 1000));
      await session.step(61, "And the \"limited\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "limited", el("pivot table viewer"), "false"));
      await session.step(62, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(63, "And the \"aggregated columns\" reading of pivot table viewer should be 2", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 2));
      await session.step(64, "When user sets \"Group By Column Names\" property of pivot table viewer to \"DIS_POP\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "DIS_POP"));
      await session.step(65, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(66, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Save parameters, and the history menu puts the configuration back", async () => {
      await session.step(70, "When user sets \"Group By Column Names\" property of pivot table viewer to \"RACE\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "RACE"));
      await session.step(71, "And user sets \"Aggregate Column Names\" property of pivot table viewer to \"WEIGHT\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "WEIGHT"));
      await session.step(72, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(73, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), ""));
      await session.step(74, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(WEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(WEIGHT)"));
      await session.step(75, "And the \"aggregated rows\" reading of pivot table viewer should be 4", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 4));
      await session.step(76, "When user picks \"Save parameters\" from the history menu of pivot table viewer", () => pickFromHistory(page, "Save parameters"));
      await session.step(77, "Then the \"history entries\" reading of pivot table viewer should be \"key(RACE),avg(WEIGHT)\"", () => readingReads(page, "history entries", el("pivot table viewer"), "key(RACE),avg(WEIGHT)"));
      await session.step(78, "When user sets \"Group By Column Names\" property of pivot table viewer to \"DIS_POP\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "DIS_POP"));
      await session.step(79, "And user sets \"Aggregate Column Names\" property of pivot table viewer to \"AGE\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "AGE"));
      await session.step(80, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(81, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(82, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(83, "When user picks \"key(RACE),avg(WEIGHT)\" from the history menu of pivot table viewer", () => pickFromHistory(page, "key(RACE),avg(WEIGHT)"));
      await session.step(84, "Then the \"group by\" reading of pivot table viewer should be \"RACE\"", () => readingReads(page, "group by", el("pivot table viewer"), "RACE"));
      await session.step(85, "And the \"aggregate\" reading of pivot table viewer should be \"avg(WEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(WEIGHT)"));
      await session.step(86, "And the \"aggregated rows\" reading of pivot table viewer should be 4", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 4));
      await session.step(87, "And the \"text of grid cell 1 of RACE\" reading of pivot table viewer should be \"Asian\"", () => readingReads(page, "text of grid cell 1 of RACE", el("pivot table viewer"), "Asian"));
      await session.step(88, "And the \"text of grid cell 1 of avg(WEIGHT)\" reading of pivot table viewer should be \"70.52\"", () => readingReads(page, "text of grid cell 1 of avg(WEIGHT)", el("pivot table viewer"), "70.52"));
      await session.step(89, "And the aggregated values of pivot table viewer should match \"avg(WEIGHT)\" grouped by \"RACE\"", () => aggregationMatches(page, el("pivot table viewer"), "avg(WEIGHT)", "RACE"));
      await session.step(90, "When user clears the saved pivot table parameters", () => clearSavedParameters(page));
      await session.step(91, "And user sets \"Group By Column Names\" property of pivot table viewer to \"DIS_POP\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "DIS_POP"));
      await session.step(92, "And user sets \"Aggregate Column Names\" property of pivot table viewer to \"AGE\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "AGE"));
      await session.step(93, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(94, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(95, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved configuration is offered only while the table still has its columns", async () => {
      await session.step(99, "When user sets \"Group By Column Names\" property of pivot table viewer to \"RACE\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "RACE"));
      await session.step(100, "And user sets \"Aggregate Column Names\" property of pivot table viewer to \"WEIGHT\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "WEIGHT"));
      await session.step(101, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(102, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), ""));
      await session.step(103, "And user picks \"Save parameters\" from the history menu of pivot table viewer", () => pickFromHistory(page, "Save parameters"));
      await session.step(104, "Then the \"history entries\" reading of pivot table viewer should contain \"avg(WEIGHT)\"", () => readingContains(page, "history entries", el("pivot table viewer"), "avg(WEIGHT)"));
      await session.step(105, "When user sets \"Group By Column Names\" property of pivot table viewer to \"DIS_POP\"", () => setProperty(page, "Group By Column Names", el("pivot table viewer"), "DIS_POP"));
      await session.step(106, "And user sets \"Aggregate Column Names\" property of pivot table viewer to \"AGE\"", () => setProperty(page, "Aggregate Column Names", el("pivot table viewer"), "AGE"));
      await session.step(107, "And user sets \"Aggregate Agg Types\" property of pivot table viewer to \"avg\"", () => setProperty(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(108, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(109, "And user renames \"WEIGHT\" column to \"MASS\"", () => renameColumn(page, "WEIGHT", "MASS"));
      await session.step(110, "And user clicks on close icon of pivot table viewer", () => clickOn(page, el("close icon of pivot table viewer")));
      await session.step(111, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
      await session.step(112, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(113, "And the \"history entries\" reading of pivot table viewer should not contain \"WEIGHT\"", () => readingNotContains(page, "history entries", el("pivot table viewer"), "WEIGHT"));
      await session.step(114, "When user renames \"MASS\" column to \"WEIGHT\"", () => renameColumn(page, "MASS", "WEIGHT"));
      await session.step(115, "And user clicks on close icon of pivot table viewer", () => clickOn(page, el("close icon of pivot table viewer")));
      await session.step(116, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
      await session.step(117, "Then the \"history entries\" reading of pivot table viewer should contain \"avg(WEIGHT)\"", () => readingContains(page, "history entries", el("pivot table viewer"), "avg(WEIGHT)"));
      await session.step(118, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(119, "When user clears the saved pivot table parameters", () => clearSavedParameters(page));
      await session.step(120, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Link Tables links the source with the table the pivot published (public#3470)", async () => {
      await session.step(124, "When user clicks on the \"add to workspace\" area of pivot table viewer", () => clickArea(page, "add to workspace", el("pivot table viewer")));
      await session.step(125, "Then table \"demog-1000 aggregation\" should be open", () => tableOpen(page, "demog-1000 aggregation"));
      await session.step(126, "When user picks \"Data > Link Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Link Tables..."));
      await session.step(127, "And user sets the tables of the Link Tables dialog to \"demog-1000\" and \"demog-1000 aggregation\"", () => setLinkTables(page, "demog-1000", "demog-1000 aggregation"));
      await session.step(128, "And user sets key columns 1 of the Link Tables dialog to \"DIS_POP\" and \"DIS_POP\"", () => setLinkKeyPair(page, 1, "DIS_POP", "DIS_POP"));
      await session.step(129, "And user clicks on LINK button in \"Link Tables\" dialog", () => clickOn(page, el("LINK button in \"Link Tables\" dialog")));
      await session.step(130, "Then UNLINK button in \"Link Tables\" dialog should be visible", () => shouldBe(page, el("UNLINK button in \"Link Tables\" dialog"), "visible"));
      await session.step(131, "And 1 rows of table \"demog-1000 aggregation\" should pass the filter", () => tableFilterCount(page, 1, "demog-1000 aggregation"));
      await session.step(132, "When user clicks on CLOSE button in \"Link Tables\" dialog", () => clickOn(page, el("CLOSE button in \"Link Tables\" dialog")));
      await session.step(133, "Then no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(134, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
