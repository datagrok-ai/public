/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-interaction.feature
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
import {clickOn, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, addRangeFilter, allOfFiltered, clearSelection, filterIsExactly, filterIsExactlyCategory, filterPasses, filterPassesAll, noneOfFiltered, noneSelected, onlyOfAnySelected, onlyOfSelected, rowCount, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, clickArea, clickAreaHolding, hasArea, hasNoArea, moreHighlight, noErrors, propertyShouldBe, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart clicks — Select and Filter", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-interaction.feature", import.meta.url);
  test("Pie chart clicks — Select and Filter", {tag: ["@journey", "@viewers", "@realizes:viewers.pie-chart", "@realizes:core#982"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Row Source","All"]]), [["Category","RACE"],["Row Source","All"]]);
    await session.step(26, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(27, "And the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(28, "And no rows should be selected", () => noneSelected(page));
    await run.scenario("On Click Select selects exactly the wedge's category, and the filtered-rows overlay does not follow the selection (core#982)", async () => {
      await session.step(32, "Then \"On Click\" property of pie chart viewer should be \"Select\"", () => propertyShouldBe(page, "On Click", el("pie chart viewer"), "Select"));
      await session.step(33, "And pie chart viewer should not have a \"selected Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Caucasian"));
      await session.step(34, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(35, "Then only rows where \"RACE\" is \"Caucasian\" should be selected", () => onlyOfSelected(page, "RACE", "Caucasian"));
      await session.step(36, "And 896 rows should be selected", () => selectedRowCount(page, 896));
      await session.step(37, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(38, "And pie chart viewer should have a \"selected Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "selected Caucasian"));
      await session.step(39, "And pie chart viewer should not have a \"filtered Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "filtered Caucasian"));
      await session.step(40, "And the \"pie\" area of pie chart viewer should contain the color \"#FF8C00\"", () => areaColor(page, "pie", el("pie chart viewer"), "#FF8C00"));
      await session.step(41, "And pie chart viewer should show more selection highlight than before", () => moreHighlight(page, el("pie chart viewer")));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Another wedge switches the selection and Escape clears it", async () => {
      await session.step(45, "When user clicks on the \"slice Asian\" area of pie chart viewer", () => clickArea(page, "slice Asian", el("pie chart viewer")));
      await session.step(46, "Then only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(47, "And 15 rows should be selected", () => selectedRowCount(page, 15));
      await session.step(48, "And pie chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(49, "And pie chart viewer should not have a \"selected Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Caucasian"));
      await session.step(50, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(51, "Then no rows should be selected", () => noneSelected(page));
      await session.step(52, "And pie chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control adds a category to the selection", async () => {
      await session.step(56, "When user clicks on the \"slice Asian\" area of pie chart viewer", () => clickArea(page, "slice Asian", el("pie chart viewer")));
      await session.step(57, "And user clicks on the \"slice Black\" area of pie chart viewer holding Control", () => clickAreaHolding(page, "slice Black", el("pie chart viewer"), "Control"));
      await session.step(58, "Then only rows where \"RACE\" is one of \"Asian, Black\" should be selected", () => onlyOfAnySelected(page, "RACE", "Asian, Black"));
      await session.step(59, "And 42 rows should be selected", () => selectedRowCount(page, 42));
      await session.step(60, "And pie chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(61, "And pie chart viewer should have a \"selected Black\" area", () => hasArea(page, el("pie chart viewer"), "selected Black"));
      await session.step(62, "And pie chart viewer should not have a \"selected Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Caucasian"));
      await session.step(63, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(64, "Then no rows should be selected", () => noneSelected(page));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The selected-rows overlay needs a count aggregation and its own flag", async () => {
      await session.step(68, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(69, "Then pie chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(70, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"avg\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "avg"));
      await session.step(71, "Then pie chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(72, "And only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(73, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"count\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "count"));
      await session.step(74, "Then pie chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(75, "When user sets \"Show Selected Rows\" property of pie chart viewer to \"false\"", () => setProperty(page, "Show Selected Rows", el("pie chart viewer"), "false"));
      await session.step(76, "Then pie chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(77, "When user sets \"Show Selected Rows\" property of pie chart viewer to \"true\"", () => setProperty(page, "Show Selected Rows", el("pie chart viewer"), "true"));
      await session.step(78, "Then pie chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(79, "When user clears the row selection", () => clearSelection(page));
      await session.step(80, "Then no rows should be selected", () => noneSelected(page));
      await session.step(81, "And pie chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On Click Filter filters the table to the wedge's category, and a Category change releases that filter (core#982)", async () => {
      await session.step(86, "When user sets \"On Click\" property of pie chart viewer to \"Filter\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Filter"));
      await session.step(87, "And user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(88, "Then the filter should pass exactly the rows where \"RACE\" is \"Caucasian\"", () => filterIsExactlyCategory(page, "RACE", "Caucasian"));
      await session.step(89, "And 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(90, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(91, "And pie chart viewer should have a \"filtered Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "filtered Caucasian"));
      await session.step(92, "When user clicks on the \"empty space\" area of pie chart viewer", () => clickArea(page, "empty space", el("pie chart viewer")));
      await session.step(93, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(94, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(95, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(96, "When user sets \"Category\" property of pie chart viewer to \"SEX\"", () => setProperty(page, "Category", el("pie chart viewer"), "SEX"));
      await session.step(97, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(98, "And the \"slices\" reading of pie chart viewer should be 2", () => readingIs(page, "slices", el("pie chart viewer"), 2));
      await session.step(99, "When user sets \"Category\" property of pie chart viewer to \"RACE\"", () => setProperty(page, "Category", el("pie chart viewer"), "RACE"));
      await session.step(100, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control adds a second category to the filter", async () => {
      await session.step(104, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(105, "And user clicks on the \"slice Asian\" area of pie chart viewer holding Control", () => clickAreaHolding(page, "slice Asian", el("pie chart viewer"), "Control"));
      await session.step(106, "Then 911 rows should pass the filter", () => filterPasses(page, 911));
      await session.step(107, "And all rows where \"RACE\" is \"Asian\" should pass the filter", () => allOfFiltered(page, "RACE", "Asian"));
      await session.step(108, "And all rows where \"RACE\" is \"Caucasian\" should pass the filter", () => allOfFiltered(page, "RACE", "Caucasian"));
      await session.step(109, "And no rows where \"RACE\" is \"Black\" should pass the filter", () => noneOfFiltered(page, "RACE", "Black"));
      await session.step(110, "When user clicks on the \"empty space\" area of pie chart viewer", () => clickArea(page, "empty space", el("pie chart viewer")));
      await session.step(111, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the viewer releases the filter it held", async () => {
      await session.step(115, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(116, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(117, "When user clicks on close icon of pie chart viewer", () => clickOn(page, el("close icon of pie chart viewer")));
      await session.step(118, "Then pie chart viewer should be absent", () => shouldBe(page, el("pie chart viewer"), "absent"));
      await session.step(119, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(120, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Row Source","All"],["On Click","Filter"]]), [["Category","RACE"],["Row Source","All"],["On Click","Filter"]]);
      await session.step(124, "Then pie chart viewer should be visible", () => shouldBe(page, el("pie chart viewer"), "visible"));
      await session.step(125, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(126, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The wedge's filter intersects with a filter card and gives its rows back", async () => {
      await session.step(130, "When user adds a range filter on \"AGE\" from 30 to 50", () => addRangeFilter(page, "AGE", 30, 50));
      await session.step(131, "Then 494 rows should pass the filter", () => filterPasses(page, 494));
      await session.step(132, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(133, "Then 442 rows should pass the filter", () => filterPasses(page, 442));
      await session.step(134, "And no rows where \"RACE\" is \"Asian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Asian"));
      await session.step(135, "And no rows where \"RACE\" is \"Black\" should pass the filter", () => noneOfFiltered(page, "RACE", "Black"));
      await session.step(136, "When user clicks on the \"empty space\" area of pie chart viewer", () => clickArea(page, "empty space", el("pie chart viewer")));
      await session.step(137, "Then 494 rows should pass the filter", () => filterPasses(page, 494));
      await session.step(138, "And the filter should pass exactly the rows where \"AGE\" is between 30 and 50", () => filterIsExactly(page, "AGE", 30, 50));
      await session.step(139, "When user adds a range filter on \"AGE\" from 18 to 89", () => addRangeFilter(page, "AGE", 18, 89));
      await session.step(140, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(141, "When user sets \"On Click\" property of pie chart viewer to \"Select\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Select"));
      await session.step(142, "Then \"On Click\" property of pie chart viewer should be \"Select\"", () => propertyShouldBe(page, "On Click", el("pie chart viewer"), "Select"));
      await session.step(143, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("After On Click goes to Select and back the old wedge's filter does not come back on a filter change (core#982)", async () => {
      await session.step(147, "When user sets \"On Click\" property of pie chart viewer to \"Filter\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Filter"));
      await session.step(148, "And user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(149, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(150, "When user sets \"On Click\" property of pie chart viewer to \"Select\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Select"));
      await session.step(151, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(152, "When user sets \"On Click\" property of pie chart viewer to \"Filter\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Filter"));
      await session.step(153, "And user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
      await session.step(154, "Then the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
      await session.step(155, "When user adds a categorical filter on \"SEX\" keeping \"F, M\"", () => addCategoricalFilter(page, "SEX", "F, M"));
      await session.step(156, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(157, "When user sets \"On Click\" property of pie chart viewer to \"Select\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Select"));
      await session.step(158, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
