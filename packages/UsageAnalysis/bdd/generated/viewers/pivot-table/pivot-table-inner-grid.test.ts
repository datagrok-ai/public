/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-inner-grid.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pivot-table, public#2606]
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
import {addCategoricalFilter, selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {clickPlainCheckbox, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerColumn} from '@datagrok-libraries/bdd/bindings/tiers/viewers/pivot-table';
import {addViewer, areaColor, clickArea, dragAreaBy, hasArea, hasNoArea, noErrors, pickFromAreaContextMenu, readingDoesNotRead, readingHigher, readingIs, readingLower, readingReads, readingsDiffer, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — the inner grid that shows the aggregation", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-inner-grid.feature", import.meta.url);
  test("Pivot table — the inner grid that shows the aggregation", {tag: ["@journey", "@viewers", "@realizes:viewers.pivot-table", "@realizes:public#2606"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(22, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(23, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
    await session.step(24, "And the \"text of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"52.30\"", () => readingReads(page, "text of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "52.30"));
    await run.scenario("Linear colour coding paints the value column from its own numbers", async () => {
      await session.step(27, "Then the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Off\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Off"));
      await session.step(28, "And the \"color of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"#ffffff\"", () => readingReads(page, "color of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "#ffffff"));
      await session.step(29, "When user picks \"Grid > Color Coding > Linear\" from the context menu of the \"grid header None avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Color Coding > Linear", "grid header None avg(AGE)", el("pivot table viewer")));
      await session.step(30, "Then the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Linear"));
      await session.step(31, "And the \"color of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should not be \"#ffffff\"", () => readingDoesNotRead(page, "color of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "#ffffff"));
      await session.step(32, "And the \"color of grid cell 2 of None avg(AGE)\" reading of pivot table viewer should be \"#0000ff\"", () => readingReads(page, "color of grid cell 2 of None avg(AGE)", el("pivot table viewer"), "#0000ff"));
      await session.step(33, "And the \"color of grid cell 5 of None avg(AGE)\" and \"color of grid cell 2 of None avg(AGE)\" readings of pivot table viewer should differ", () => readingsDiffer(page, "color of grid cell 5 of None avg(AGE)", "color of grid cell 2 of None avg(AGE)", el("pivot table viewer")));
      await session.step(34, "And the \"grid cell 2 of None avg(AGE)\" area of pivot table viewer should contain the color \"#0000FF\"", () => areaColor(page, "grid cell 2 of None avg(AGE)", el("pivot table viewer"), "#0000FF"));
      await session.step(35, "And the \"color coding of High avg(AGE)\" reading of pivot table viewer should be \"Off\"", () => readingReads(page, "color coding of High avg(AGE)", el("pivot table viewer"), "Off"));
      await session.step(36, "When user picks \"Grid > Color Coding > Off\" from the context menu of the \"grid header None avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Color Coding > Off", "grid header None avg(AGE)", el("pivot table viewer")));
      await session.step(37, "Then the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Off\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Off"));
      await session.step(38, "And the \"color of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"#ffffff\"", () => readingReads(page, "color of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "#ffffff"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Grid > Hide takes a column away and Order or Hide Columns brings it back", async () => {
      await session.step(42, "Then the \"column visible of None avg(AGE)\" reading of pivot table viewer should be \"true\"", () => readingReads(page, "column visible of None avg(AGE)", el("pivot table viewer"), "true"));
      await session.step(43, "And the \"columns shown\" reading of pivot table viewer should be 7", () => readingIs(page, "columns shown", el("pivot table viewer"), 7));
      await session.step(44, "When user picks \"Grid > Hide\" from the context menu of the \"grid header None avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Hide", "grid header None avg(AGE)", el("pivot table viewer")));
      await session.step(45, "Then the \"column visible of None avg(AGE)\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "column visible of None avg(AGE)", el("pivot table viewer"), "false"));
      await session.step(46, "And the \"columns shown\" reading of pivot table viewer should be 6", () => readingIs(page, "columns shown", el("pivot table viewer"), 6));
      await session.step(47, "And pivot table viewer should not have a \"grid header None avg(AGE)\" area", () => hasNoArea(page, el("pivot table viewer"), "grid header None avg(AGE)"));
      await session.step(48, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(49, "When user picks \"Grid > Order or Hide Columns...\" from the context menu of the \"grid header Low avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Order or Hide Columns...", "grid header Low avg(AGE)", el("pivot table viewer")));
      await session.step(50, "Then \"Order or Hide Columns\" dialog should be visible", () => shouldBe(page, el("\"Order or Hide Columns\" dialog"), "visible"));
      await session.step(51, "When user clicks the plain checkbox in the \"Order or Hide Columns\" dialog", () => clickPlainCheckbox(page, "Order or Hide Columns"));
      await session.step(52, "Then the \"column visible of None avg(AGE)\" reading of pivot table viewer should be \"true\"", () => readingReads(page, "column visible of None avg(AGE)", el("pivot table viewer"), "true"));
      await session.step(53, "And the \"columns shown\" reading of pivot table viewer should be 7", () => readingIs(page, "columns shown", el("pivot table viewer"), 7));
      await session.step(54, "When user clicks on CLOSE button in \"Order or Hide Columns\" dialog", () => clickOn(page, el("CLOSE button in \"Order or Hide Columns\" dialog")));
      await session.step(55, "Then \"Order or Hide Columns\" dialog should be absent", () => shouldBe(page, el("\"Order or Hide Columns\" dialog"), "absent"));
      await session.step(56, "And pivot table viewer should have a \"grid header None avg(AGE)\" area", () => hasArea(page, el("pivot table viewer"), "grid header None avg(AGE)"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer picker adds one in-cell viewer column per pivot category", async () => {
      await session.step(60, "Then the \"viewer columns\" reading of pivot table viewer should be 0", () => readingIs(page, "viewer columns", el("pivot table viewer"), 0));
      await session.step(61, "And the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
      await session.step(62, "When user adds a \"Scatter plot\" viewer column to pivot table viewer", () => addViewerColumn(page, "Scatter plot"));
      await session.step(63, "Then the \"viewer columns\" reading of pivot table viewer should be 5", () => readingIs(page, "viewer columns", el("pivot table viewer"), 5));
      await session.step(64, "And the \"viewer column names\" reading of pivot table viewer should be \"Critical, High, Low, Medium, None\"", () => readingReads(page, "viewer column names", el("pivot table viewer"), "Critical, High, Low, Medium, None"));
      await session.step(65, "And the \"aggregate\" reading of pivot table viewer should contain \"Scatter plot\"", () => readingContains(page, "aggregate", el("pivot table viewer"), "Scatter plot"));
      await session.step(66, "When user clicks on the \"remove aggregate chip Scatter plot\" area of pivot table viewer", () => clickArea(page, "remove aggregate chip Scatter plot", el("pivot table viewer")));
      await session.step(67, "Then the \"viewer columns\" reading of pivot table viewer should be 0", () => readingIs(page, "viewer columns", el("pivot table viewer"), 0));
      await session.step(68, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With two pivot columns the picker adds no viewer column at all", async () => {
      await session.step(72, "When user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY, SEX\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY, SEX"));
      await session.step(73, "Then the \"aggregated columns\" reading of pivot table viewer should be 10", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 10));
      await session.step(74, "When user adds a \"Scatter plot\" viewer column to pivot table viewer", () => addViewerColumn(page, "Scatter plot"));
      await session.step(75, "Then the \"aggregate\" reading of pivot table viewer should contain \"Scatter plot\"", () => readingContains(page, "aggregate", el("pivot table viewer"), "Scatter plot"));
      await session.step(76, "And the \"viewer columns\" reading of pivot table viewer should be 0", () => readingIs(page, "viewer columns", el("pivot table viewer"), 0));
      await session.step(77, "And the \"viewer column names\" reading of pivot table viewer should be \"\"", () => readingReads(page, "viewer column names", el("pivot table viewer"), ""));
      await session.step(78, "When user clicks on the \"remove aggregate chip Scatter plot\" area of pivot table viewer", () => clickArea(page, "remove aggregate chip Scatter plot", el("pivot table viewer")));
      await session.step(79, "And user sets \"Pivot Column Names\" property of pivot table viewer to \"SEVERITY\"", () => setProperty(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(80, "Then the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(81, "And the \"viewer columns\" reading of pivot table viewer should be 0", () => readingIs(page, "viewer columns", el("pivot table viewer"), 0));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column resizer of the inner grid widens the column it belongs to", async () => {
      await session.step(85, "When user drags the \"grid column resizer None avg(AGE)\" area of pivot table viewer by 60 pixels to the right", () => dragAreaBy(page, "grid column resizer None avg(AGE)", el("pivot table viewer"), 60, "right"));
      await session.step(86, "Then the \"column width of None avg(AGE)\" reading of pivot table viewer should be higher than before", () => readingHigher(page, "column width of None avg(AGE)", el("pivot table viewer")));
      await session.step(87, "And the \"text of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"52.30\"", () => readingReads(page, "text of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "52.30"));
      await session.step(88, "When user drags the \"grid column resizer None avg(AGE)\" area of pivot table viewer by 60 pixels to the left", () => dragAreaBy(page, "grid column resizer None avg(AGE)", el("pivot table viewer"), 60, "left"));
      await session.step(89, "Then the \"column width of None avg(AGE)\" reading of pivot table viewer should be lower than before", () => readingLower(page, "column width of None avg(AGE)", el("pivot table viewer")));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("At Row Source Filtered the inner grid keeps its colouring and its hidden column through a filter (public#2606)", async () => {
      await session.step(94, "When user sets \"Row Source\" property of pivot table viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "Filtered"));
      await session.step(95, "Then the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(96, "When user picks \"Grid > Color Coding > Linear\" from the context menu of the \"grid header None avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Color Coding > Linear", "grid header None avg(AGE)", el("pivot table viewer")));
      await session.step(97, "Then the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Linear"));
      await session.step(98, "When user picks \"Grid > Hide\" from the context menu of the \"grid header High avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Hide", "grid header High avg(AGE)", el("pivot table viewer")));
      await session.step(99, "Then the \"column visible of High avg(AGE)\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "column visible of High avg(AGE)", el("pivot table viewer"), "false"));
      await session.step(100, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(101, "Then the \"rows shown\" reading of pivot table viewer should be 447", () => readingIs(page, "rows shown", el("pivot table viewer"), 447));
      await session.step(102, "And the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Linear"));
      await session.step(103, "And the \"color of grid cell 2 of None avg(AGE)\" reading of pivot table viewer should not be \"#ffffff\"", () => readingDoesNotRead(page, "color of grid cell 2 of None avg(AGE)", el("pivot table viewer"), "#ffffff"));
      await session.step(104, "And the \"column visible of High avg(AGE)\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "column visible of High avg(AGE)", el("pivot table viewer"), "false"));
      await session.step(105, "When user adds a categorical filter on \"SEX\" keeping \"F, M\"", () => addCategoricalFilter(page, "SEX", "F, M"));
      await session.step(106, "And user picks \"Grid > Color Coding > Off\" from the context menu of the \"grid header None avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Color Coding > Off", "grid header None avg(AGE)", el("pivot table viewer")));
      await session.step(107, "And user picks \"Grid > Order or Hide Columns...\" from the context menu of the \"grid header Low avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Order or Hide Columns...", "grid header Low avg(AGE)", el("pivot table viewer")));
      await session.step(108, "And user clicks the plain checkbox in the \"Order or Hide Columns\" dialog", () => clickPlainCheckbox(page, "Order or Hide Columns"));
      await session.step(109, "And user clicks on CLOSE button in \"Order or Hide Columns\" dialog", () => clickOn(page, el("CLOSE button in \"Order or Hide Columns\" dialog")));
      await session.step(110, "Then the \"column visible of High avg(AGE)\" reading of pivot table viewer should be \"true\"", () => readingReads(page, "column visible of High avg(AGE)", el("pivot table viewer"), "true"));
      await session.step(111, "And the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Off\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Off"));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("At Row Source Selected the inner grid keeps its colouring and its hidden column through a selection (public#2606)", async () => {
      await session.step(116, "Given user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
      await session.step(117, "When user sets \"Row Source\" property of pivot table viewer to \"Selected\"", () => setProperty(page, "Row Source", el("pivot table viewer"), "Selected"));
      await session.step(118, "Then the \"rows shown\" reading of pivot table viewer should be 553", () => readingIs(page, "rows shown", el("pivot table viewer"), 553));
      await session.step(119, "When user picks \"Grid > Color Coding > Linear\" from the context menu of the \"grid header None avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Color Coding > Linear", "grid header None avg(AGE)", el("pivot table viewer")));
      await session.step(120, "Then the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Linear"));
      await session.step(121, "When user picks \"Grid > Hide\" from the context menu of the \"grid header High avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Hide", "grid header High avg(AGE)", el("pivot table viewer")));
      await session.step(122, "Then the \"column visible of High avg(AGE)\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "column visible of High avg(AGE)", el("pivot table viewer"), "false"));
      await session.step(123, "When user selects rows where \"SEX\" is \"M\"", () => selectWhereIs(page, "SEX", "M"));
      await session.step(124, "Then the \"rows shown\" reading of pivot table viewer should be 447", () => readingIs(page, "rows shown", el("pivot table viewer"), 447));
      await session.step(125, "And the \"color coding of None avg(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None avg(AGE)", el("pivot table viewer"), "Linear"));
      await session.step(126, "And the \"color of grid cell 2 of None avg(AGE)\" reading of pivot table viewer should not be \"#ffffff\"", () => readingDoesNotRead(page, "color of grid cell 2 of None avg(AGE)", el("pivot table viewer"), "#ffffff"));
      await session.step(127, "And the \"column visible of High avg(AGE)\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "column visible of High avg(AGE)", el("pivot table viewer"), "false"));
      await session.step(128, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
