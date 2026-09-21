/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-inner-grid.feature
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
import {addViewerColumn} from '../../../bindings/pivot-table.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clickPlainCheckbox, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaColor, clickArea, dragAreaBy, hasArea, hasNoArea, noErrors, pickFromAreaContextMenu, readingDoesNotRead, readingHigher, readingIs, readingLower, readingReads, readingsDiffer, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — the inner grid that shows the aggregation", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-inner-grid.feature", import.meta.url);
  test("Pivot table — the inner grid that shows the aggregation", {tag: ["@journey", "@viewers", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(20, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(21, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
    await session.step(22, "And the \"text of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"52.30\"", () => readingReads(page, "text of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "52.30"));
    await run.scenario("Linear colour coding paints the value column from its own numbers", async () => {
      await session.step(25, "Then the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Off\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Off"));
      await session.step(26, "And the \"color of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"#ffffff\"", () => readingReads(page, "color of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "#ffffff"));
      await session.step(27, "When user picks \"Grid > Color Coding > Linear\" from the context menu of the \"grid header None avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Color Coding > Linear", "grid header None avg(AGE)", el("pivot table viewer")));
      await session.step(28, "Then the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Linear"));
      await session.step(29, "And the \"color of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should not be \"#ffffff\"", () => readingDoesNotRead(page, "color of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "#ffffff"));
      await session.step(30, "And the \"color of grid cell 2 of None avg(AGE)\" reading of pivot table viewer should be \"#0000ff\"", () => readingReads(page, "color of grid cell 2 of None avg(AGE)", el("pivot table viewer"), "#0000ff"));
      await session.step(31, "And the \"color of grid cell 5 of None avg(AGE)\" and \"color of grid cell 2 of None avg(AGE)\" readings of pivot table viewer should differ", () => readingsDiffer(page, "color of grid cell 5 of None avg(AGE)", "color of grid cell 2 of None avg(AGE)", el("pivot table viewer")));
      await session.step(32, "And the \"grid cell 2 of None avg(AGE)\" area of pivot table viewer should contain the color \"#0000FF\"", () => areaColor(page, "grid cell 2 of None avg(AGE)", el("pivot table viewer"), "#0000FF"));
      await session.step(33, "And the \"color coding of High avg(AGE)\" reading of pivot table viewer should be \"Off\"", () => readingReads(page, "color coding of High avg(AGE)", el("pivot table viewer"), "Off"));
      await session.step(34, "When user picks \"Grid > Color Coding > Off\" from the context menu of the \"grid header None avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Color Coding > Off", "grid header None avg(AGE)", el("pivot table viewer")));
      await session.step(35, "Then the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Off\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Off"));
      await session.step(36, "And the \"color of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"#ffffff\"", () => readingReads(page, "color of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "#ffffff"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Grid > Hide takes a column away and Order or Hide Columns brings it back", async () => {
      await session.step(40, "Then the \"column visible of None avg(AGE)\" reading of pivot table viewer should be \"true\"", () => readingReads(page, "column visible of None avg(AGE)", el("pivot table viewer"), "true"));
      await session.step(41, "And the \"columns shown\" reading of pivot table viewer should be 7", () => readingIs(page, "columns shown", el("pivot table viewer"), 7));
      await session.step(42, "When user picks \"Grid > Hide\" from the context menu of the \"grid header None avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Hide", "grid header None avg(AGE)", el("pivot table viewer")));
      await session.step(43, "Then the \"column visible of None avg(AGE)\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "column visible of None avg(AGE)", el("pivot table viewer"), "false"));
      await session.step(44, "And the \"columns shown\" reading of pivot table viewer should be 6", () => readingIs(page, "columns shown", el("pivot table viewer"), 6));
      await session.step(45, "And pivot table viewer should not have a \"grid header None avg(AGE)\" area", () => hasNoArea(page, el("pivot table viewer"), "grid header None avg(AGE)"));
      await session.step(46, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(47, "When user picks \"Grid > Order or Hide Columns...\" from the context menu of the \"grid header Low avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Order or Hide Columns...", "grid header Low avg(AGE)", el("pivot table viewer")));
      await session.step(48, "Then \"Order or Hide Columns\" dialog should be visible", () => shouldBe(page, el("\"Order or Hide Columns\" dialog"), "visible"));
      await session.step(49, "When user clicks the plain checkbox in the \"Order or Hide Columns\" dialog", () => clickPlainCheckbox(page, "Order or Hide Columns"));
      await session.step(50, "Then the \"column visible of None avg(AGE)\" reading of pivot table viewer should be \"true\"", () => readingReads(page, "column visible of None avg(AGE)", el("pivot table viewer"), "true"));
      await session.step(51, "And the \"columns shown\" reading of pivot table viewer should be 7", () => readingIs(page, "columns shown", el("pivot table viewer"), 7));
      await session.step(52, "When user clicks on CLOSE button in \"Order or Hide Columns\" dialog", () => clickOn(page, el("CLOSE button in \"Order or Hide Columns\" dialog")));
      await session.step(53, "Then \"Order or Hide Columns\" dialog should be absent", () => shouldBe(page, el("\"Order or Hide Columns\" dialog"), "absent"));
      await session.step(54, "And pivot table viewer should have a \"grid header None avg(AGE)\" area", () => hasArea(page, el("pivot table viewer"), "grid header None avg(AGE)"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer picker adds one in-cell viewer column per pivot category", async () => {
      await session.step(58, "Then the \"viewer columns\" reading of pivot table viewer should be 0", () => readingIs(page, "viewer columns", el("pivot table viewer"), 0));
      await session.step(59, "And the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
      await session.step(60, "When user adds a \"Scatter plot\" viewer column to pivot table viewer", () => addViewerColumn(page, "Scatter plot"));
      await session.step(61, "Then the \"viewer columns\" reading of pivot table viewer should be 5", () => readingIs(page, "viewer columns", el("pivot table viewer"), 5));
      await session.step(62, "And the \"viewer column names\" reading of pivot table viewer should be \"Critical, High, Low, Medium, None\"", () => readingReads(page, "viewer column names", el("pivot table viewer"), "Critical, High, Low, Medium, None"));
      await session.step(63, "And the \"aggregate\" reading of pivot table viewer should contain \"Scatter plot\"", () => readingContains(page, "aggregate", el("pivot table viewer"), "Scatter plot"));
      await session.step(64, "When user clicks on the \"remove aggregate chip Scatter plot\" area of pivot table viewer", () => clickArea(page, "remove aggregate chip Scatter plot", el("pivot table viewer")));
      await session.step(65, "Then the \"viewer columns\" reading of pivot table viewer should be 0", () => readingIs(page, "viewer columns", el("pivot table viewer"), 0));
      await session.step(66, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With two pivot columns the picker adds no viewer column at all", async () => {
      await session.step(70, "When user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY, SEX\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY, SEX"));
      await session.step(71, "Then the \"aggregated columns\" reading of pivot table viewer should be 10", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 10));
      await session.step(72, "When user adds a \"Scatter plot\" viewer column to pivot table viewer", () => addViewerColumn(page, "Scatter plot"));
      await session.step(73, "Then the \"aggregate\" reading of pivot table viewer should contain \"Scatter plot\"", () => readingContains(page, "aggregate", el("pivot table viewer"), "Scatter plot"));
      await session.step(74, "And the \"viewer columns\" reading of pivot table viewer should be 0", () => readingIs(page, "viewer columns", el("pivot table viewer"), 0));
      await session.step(75, "And the \"viewer column names\" reading of pivot table viewer should be \"\"", () => readingReads(page, "viewer column names", el("pivot table viewer"), ""));
      await session.step(76, "When user clicks on the \"remove aggregate chip Scatter plot\" area of pivot table viewer", () => clickArea(page, "remove aggregate chip Scatter plot", el("pivot table viewer")));
      await session.step(77, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(78, "Then the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(79, "And the \"viewer columns\" reading of pivot table viewer should be 0", () => readingIs(page, "viewer columns", el("pivot table viewer"), 0));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column resizer of the inner grid widens the column it belongs to", async () => {
      await session.step(83, "When user drags the \"grid column resizer None avg(AGE)\" area of pivot table viewer by 60 pixels to the right", () => dragAreaBy(page, "grid column resizer None avg(AGE)", el("pivot table viewer"), 60, "right"));
      await session.step(84, "Then the \"column width of None avg(AGE)\" reading of pivot table viewer should be higher than before", () => readingHigher(page, "column width of None avg(AGE)", el("pivot table viewer")));
      await session.step(85, "And the \"text of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"52.30\"", () => readingReads(page, "text of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "52.30"));
      await session.step(86, "When user drags the \"grid column resizer None avg(AGE)\" area of pivot table viewer by 60 pixels to the left", () => dragAreaBy(page, "grid column resizer None avg(AGE)", el("pivot table viewer"), 60, "left"));
      await session.step(87, "Then the \"column width of None avg(AGE)\" reading of pivot table viewer should be lower than before", () => readingLower(page, "column width of None avg(AGE)", el("pivot table viewer")));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
