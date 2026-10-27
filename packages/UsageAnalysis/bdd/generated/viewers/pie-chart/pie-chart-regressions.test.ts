/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart, core#982]
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
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, addCategoricalFilter, allRowsSelected, clearSelection, filterIsExactlyCategory, filterPasses, filterPassesAll, filterTo, noneOfFiltered, removeColumn, resetFilter, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openDatasetRowsAs, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, boundTable, clickArea, hasArea, hasNoArea, noErrors, propertiesShouldBe, readingBetween, readingIs, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart regressions from the pie chart audit", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-regressions.feature", import.meta.url);
  test("Under the pie's own Filter a wedge click leaves only that category in the table (core#982)", {tag: ["@viewers", "@realizes:viewers.pie-chart", "@realizes:core#982"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(23, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["On Click","Filter"],["Filter","${SEX} == \"F\""]]), [["Category","RACE"],["On Click","Filter"],["Filter","${SEX} == \"F\""]]);
    await session.step(27, "Then pie chart viewer should show 553 rows", () => showsRows(page, el("pie chart viewer"), 553));
    await session.step(28, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(29, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
    await session.step(30, "Then the filter should pass exactly the rows where \"RACE\" is \"Caucasian\"", () => filterIsExactlyCategory(page, "RACE", "Caucasian"));
    await session.step(31, "And no rows where \"RACE\" is \"Asian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Asian"));
    await session.step(32, "And no rows where \"RACE\" is \"Black\" should pass the filter", () => noneOfFiltered(page, "RACE", "Black"));
    await session.step(33, "When user clicks on the \"empty space\" area of pie chart viewer", () => clickArea(page, "empty space", el("pie chart viewer")));
    await session.step(34, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(35, "When user clicks on close icon of pie chart viewer", () => clickOn(page, el("close icon of pie chart viewer")));
    await session.step(36, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("Changing Category while the pie filters releases its filter (core#982)", {tag: ["@viewers", "@realizes:viewers.pie-chart", "@realizes:core#982"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(42, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Row Source","All"],["On Click","Filter"]]), [["Category","RACE"],["Row Source","All"],["On Click","Filter"]]);
    await session.step(46, "And user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
    await session.step(47, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
    await session.step(48, "When user sets \"Category\" property of pie chart viewer to \"SEX\"", () => setProperty(page, "Category", el("pie chart viewer"), "SEX"));
    await session.step(49, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(50, "And the \"slices\" reading of pie chart viewer should be 2", () => readingIs(page, "slices", el("pie chart viewer"), 2));
    await session.step(51, "When user clicks on close icon of pie chart viewer", () => clickOn(page, el("close icon of pie chart viewer")));
    await session.step(52, "Then no errors should have been logged", () => noErrors(page));
  });
  test("After On Click goes to Select and back the old categories do not come back on a filter change (core#982)", {tag: ["@viewers", "@realizes:viewers.pie-chart", "@realizes:core#982"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(57, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Row Source","All"],["On Click","Filter"]]), [["Category","RACE"],["Row Source","All"],["On Click","Filter"]]);
    await session.step(61, "And user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
    await session.step(62, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
    await session.step(63, "When user sets \"On Click\" property of pie chart viewer to \"Select\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Select"));
    await session.step(64, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(65, "When user sets \"On Click\" property of pie chart viewer to \"Filter\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Filter"));
    await session.step(66, "And user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
    await session.step(67, "Then the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
    await session.step(68, "When user adds a categorical filter on \"SEX\" keeping \"F, M\"", () => addCategoricalFilter(page, "SEX", "F, M"));
    await session.step(69, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(70, "When user clicks on close icon of pie chart viewer", () => clickOn(page, el("close icon of pie chart viewer")));
    await session.step(71, "Then no errors should have been logged", () => noErrors(page));
  });
  test("Row Source All honours Include Nulls off (core#982)", {tag: ["@viewers", "@realizes:viewers.pie-chart", "@realizes:core#982"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(76, "When user adds a calculated column \"RACE_GAPS\" with formula \"if(Mod(${AGE}, 10) == 0, null, ${RACE})\"", () => addCalculated(page, "RACE_GAPS", "if(Mod(${AGE}, 10) == 0, null, ${RACE})"));
    await session.step(77, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE_GAPS"],["Row Source","All"]]), [["Category","RACE_GAPS"],["Row Source","All"]]);
    await session.step(80, "Then pie chart viewer should have a \"slice (no value)\" area", () => hasArea(page, el("pie chart viewer"), "slice (no value)"));
    await session.step(81, "And the \"slices\" reading of pie chart viewer should be 5", () => readingIs(page, "slices", el("pie chart viewer"), 5));
    await session.step(82, "When user sets \"Include Nulls\" property of pie chart viewer to \"false\"", () => setProperty(page, "Include Nulls", el("pie chart viewer"), "false"));
    await session.step(83, "Then pie chart viewer should not have a \"slice (no value)\" area", () => hasNoArea(page, el("pie chart viewer"), "slice (no value)"));
    await session.step(84, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(85, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 793", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 793));
    await session.step(86, "When user clicks on close icon of pie chart viewer", () => clickOn(page, el("close icon of pie chart viewer")));
    await session.step(87, "And user removes \"RACE_GAPS\" column", () => removeColumn(page, "RACE_GAPS"));
    await session.step(88, "Then no errors should have been logged", () => noErrors(page));
  });
  test("The filtered-rows overlay follows the filter, not the selection (core#982)", {tag: ["@viewers", "@realizes:viewers.pie-chart", "@realizes:core#982"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(93, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Row Source","All"]]), [["Category","RACE"],["Row Source","All"]]);
    await session.step(96, "Then pie chart viewer should have a \"slice Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "slice Caucasian"));
    await session.step(97, "And pie chart viewer should not have a \"filtered Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "filtered Caucasian"));
    await session.step(98, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
    await session.step(99, "Then pie chart viewer should have a \"filtered Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "filtered Caucasian"));
    await session.step(100, "When user selects all rows", () => selectAllRows(page));
    await session.step(101, "Then all rows should be selected", () => allRowsSelected(page));
    await session.step(102, "And pie chart viewer should have a \"filtered Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "filtered Caucasian"));
    await session.step(103, "When user clears the row selection", () => clearSelection(page));
    await session.step(104, "And user resets the filter", () => resetFilter(page));
    await session.step(105, "Then pie chart viewer should not have a \"filtered Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "filtered Caucasian"));
    await session.step(106, "When user clicks on close icon of pie chart viewer", () => clickOn(page, el("close icon of pie chart viewer")));
    await session.step(107, "Then no errors should have been logged", () => noErrors(page));
  });
  test("Moved to another table with the same columns the pie keeps its Category and angle aggregation (core#982)", {tag: ["@viewers", "@realizes:viewers.pie-chart", "@realizes:core#982"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(112, "Given user opens demog-1000 dataset keeping the first 500 rows as \"demog-500\"", () => openDatasetRowsAs(page, ds("demog-1000"), 500, "demog-500"));
    await session.step(113, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
    await session.step(114, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","avg"]]), [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","avg"]]);
    await session.step(118, "Then the \"angle value of Caucasian\" reading of pie chart viewer should be between 79.487 and 79.488", () => readingBetween(page, "angle value of Caucasian", el("pie chart viewer"), 79.487, 79.488));
    await session.step(119, "When user sets \"Table\" property of pie chart viewer to \"demog-500\"", () => setProperty(page, "Table", el("pie chart viewer"), "demog-500"));
    await session.step(120, "Then pie chart viewer should be bound to table \"demog-500\"", () => boundTable(page, el("pie chart viewer"), "demog-500"));
    await session.step(121, "And pie chart viewer should show 500 rows", () => showsRows(page, el("pie chart viewer"), 500));
    await session.step(122, "And properties of pie chart viewer should be:", () => propertiesShouldBe(page, el("pie chart viewer"), [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","avg"]]), [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","avg"]]);
    await session.step(126, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(127, "And the \"angle value of Caucasian\" reading of pie chart viewer should be between 75.665 and 75.666", () => readingBetween(page, "angle value of Caucasian", el("pie chart viewer"), 75.665, 75.666));
    await session.step(128, "When user clicks on close icon of pie chart viewer", () => clickOn(page, el("close icon of pie chart viewer")));
    await session.step(129, "Then no errors should have been logged", () => noErrors(page));
  });
});
