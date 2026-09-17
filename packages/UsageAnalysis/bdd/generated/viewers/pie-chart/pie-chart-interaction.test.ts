/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-interaction.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart]
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
import {clickOn, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addRangeFilter, allOfFiltered, clearSelection, filterIsExactly, filterIsExactlyCategory, filterPasses, filterPassesAll, noneOfFiltered, noneSelected, onlyOfAnySelected, onlyOfSelected, rowCount, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, clickArea, clickAreaHolding, hasArea, hasNoArea, moreHighlight, noErrors, propertyShouldBe, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart clicks — Select and Filter", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-interaction.feature", import.meta.url);
  test("Pie chart clicks — Select and Filter", {tag: ["@journey", "@viewers", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Row Source","All"]]));
    await session.step(23, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(24, "And the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(25, "And no rows should be selected", () => noneSelected(page));
    await run.scenario("On Click Select selects exactly the wedge's category", async () => {
      await session.step(28, "Then \"On Click\" property of pie chart viewer should be \"Select\"", () => propertyShouldBe(page, "On Click", el("pie chart viewer"), "Select"));
      await session.step(29, "And pie chart viewer should not have a \"selected Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Caucasian"));
      await session.step(30, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(31, "Then only rows where \"RACE\" is \"Caucasian\" should be selected", () => onlyOfSelected(page, "RACE", "Caucasian"));
      await session.step(32, "And 896 rows should be selected", () => selectedRowCount(page, 896));
      await session.step(33, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(34, "And pie chart viewer should have a \"selected Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "selected Caucasian"));
      await session.step(35, "And the \"pie\" area of pie chart viewer should contain the color \"#FF8C00\"", () => areaColor(page, "pie", el("pie chart viewer"), "#FF8C00"));
      await session.step(36, "And pie chart viewer should show more selection highlight than before", () => moreHighlight(page, el("pie chart viewer")));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Another wedge switches the selection and Escape clears it", async () => {
      await session.step(40, "When user clicks on the \"slice Asian\" area of pie chart viewer", () => clickArea(page, "slice Asian", el("pie chart viewer")));
      await session.step(41, "Then only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(42, "And 15 rows should be selected", () => selectedRowCount(page, 15));
      await session.step(43, "And pie chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(44, "And pie chart viewer should not have a \"selected Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Caucasian"));
      await session.step(45, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(46, "Then no rows should be selected", () => noneSelected(page));
      await session.step(47, "And pie chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control adds a category to the selection", async () => {
      await session.step(51, "When user clicks on the \"slice Asian\" area of pie chart viewer", () => clickArea(page, "slice Asian", el("pie chart viewer")));
      await session.step(52, "And user clicks on the \"slice Black\" area of pie chart viewer holding Control", () => clickAreaHolding(page, "slice Black", el("pie chart viewer"), "Control"));
      await session.step(53, "Then only rows where \"RACE\" is one of \"Asian, Black\" should be selected", () => onlyOfAnySelected(page, "RACE", "Asian, Black"));
      await session.step(54, "And 42 rows should be selected", () => selectedRowCount(page, 42));
      await session.step(55, "And pie chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(56, "And pie chart viewer should have a \"selected Black\" area", () => hasArea(page, el("pie chart viewer"), "selected Black"));
      await session.step(57, "And pie chart viewer should not have a \"selected Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Caucasian"));
      await session.step(58, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(59, "Then no rows should be selected", () => noneSelected(page));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The selected-rows overlay needs a count aggregation and its own flag", async () => {
      await session.step(63, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(64, "Then pie chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(65, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"avg\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "avg"));
      await session.step(66, "Then pie chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(67, "And only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(68, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"count\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "count"));
      await session.step(69, "Then pie chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(70, "When user sets \"Show Selected Rows\" property of pie chart viewer to \"false\"", () => setProperty(page, "Show Selected Rows", el("pie chart viewer"), "false"));
      await session.step(71, "Then pie chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(72, "When user sets \"Show Selected Rows\" property of pie chart viewer to \"true\"", () => setProperty(page, "Show Selected Rows", el("pie chart viewer"), "true"));
      await session.step(73, "Then pie chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(74, "When user clears the row selection", () => clearSelection(page));
      await session.step(75, "Then no rows should be selected", () => noneSelected(page));
      await session.step(76, "And pie chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("pie chart viewer"), "selected Asian"));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On Click Filter filters the table to the wedge's category", async () => {
      await session.step(80, "When user sets \"On Click\" property of pie chart viewer to \"Filter\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Filter"));
      await session.step(81, "And user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(82, "Then the filter should pass exactly the rows where \"RACE\" is \"Caucasian\"", () => filterIsExactlyCategory(page, "RACE", "Caucasian"));
      await session.step(83, "And 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(84, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(85, "And pie chart viewer should have a \"filtered Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "filtered Caucasian"));
      await session.step(86, "When user clicks on the \"empty space\" area of pie chart viewer", () => clickArea(page, "empty space", el("pie chart viewer")));
      await session.step(87, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control adds a second category to the filter", async () => {
      await session.step(91, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(92, "And user clicks on the \"slice Asian\" area of pie chart viewer holding Control", () => clickAreaHolding(page, "slice Asian", el("pie chart viewer"), "Control"));
      await session.step(93, "Then 911 rows should pass the filter", () => filterPasses(page, 911));
      await session.step(94, "And all rows where \"RACE\" is \"Asian\" should pass the filter", () => allOfFiltered(page, "RACE", "Asian"));
      await session.step(95, "And all rows where \"RACE\" is \"Caucasian\" should pass the filter", () => allOfFiltered(page, "RACE", "Caucasian"));
      await session.step(96, "And no rows where \"RACE\" is \"Black\" should pass the filter", () => noneOfFiltered(page, "RACE", "Black"));
      await session.step(97, "When user clicks on the \"empty space\" area of pie chart viewer", () => clickArea(page, "empty space", el("pie chart viewer")));
      await session.step(98, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(99, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the viewer releases the filter it held", async () => {
      await session.step(102, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(103, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(104, "When user clicks on close icon of pie chart viewer", () => clickOn(page, el("close icon of pie chart viewer")));
      await session.step(105, "Then pie chart viewer should be absent", () => shouldBe(page, el("pie chart viewer"), "absent"));
      await session.step(106, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(107, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Row Source","All"],["On Click","Filter"]]));
      await session.step(111, "Then pie chart viewer should be visible", () => shouldBe(page, el("pie chart viewer"), "visible"));
      await session.step(112, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(113, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The wedge's filter intersects with a filter card and gives its rows back", async () => {
      await session.step(117, "When user adds a range filter on \"AGE\" from 30 to 50", () => addRangeFilter(page, "AGE", 30, 50));
      await session.step(118, "Then 494 rows should pass the filter", () => filterPasses(page, 494));
      await session.step(119, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(120, "Then 442 rows should pass the filter", () => filterPasses(page, 442));
      await session.step(121, "And no rows where \"RACE\" is \"Asian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Asian"));
      await session.step(122, "And no rows where \"RACE\" is \"Black\" should pass the filter", () => noneOfFiltered(page, "RACE", "Black"));
      await session.step(123, "When user clicks on the \"empty space\" area of pie chart viewer", () => clickArea(page, "empty space", el("pie chart viewer")));
      await session.step(124, "Then 494 rows should pass the filter", () => filterPasses(page, 494));
      await session.step(125, "And the filter should pass exactly the rows where \"AGE\" is between 30 and 50", () => filterIsExactly(page, "AGE", 30, 50));
      await session.step(126, "When user adds a range filter on \"AGE\" from 18 to 89", () => addRangeFilter(page, "AGE", 18, 89));
      await session.step(127, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(128, "When user sets \"On Click\" property of pie chart viewer to \"Select\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Select"));
      await session.step(129, "Then \"On Click\" property of pie chart viewer should be \"Select\"", () => propertyShouldBe(page, "On Click", el("pie chart viewer"), "Select"));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
