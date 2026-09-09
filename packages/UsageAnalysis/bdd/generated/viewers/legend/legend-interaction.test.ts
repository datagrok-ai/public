/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/legend/legend-interaction.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.legend]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {categoricalColorIs, colorCategorical, colorOff, filterPassesAll, filterTo, filterToAnyOf, noColorCoding, noneSelected, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, clickLegendCross, clickLegendItem, clickLegendItemHolding, legendDocked, legendFewer, legendItemColor, legendItemsDiffer, legendLists, legendPlacedAsBefore, legendSameItems, legendSide, legendSlot, lessInk, moreInk, noErrors, repainted, setProperties, setProperty, showsFewerRows, showsMoreRows, showsRows, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Legend interaction", () => {
  const session = feature(test, "features/viewers/legend/legend-interaction.feature", import.meta.url);
  test("Legend interaction", {tag: ["@journey", "@viewers", "@realizes:viewers.legend"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]));
    await session.step(22, "Then legend of scatter plot viewer should be visible", () => shouldBe(page, el("legend of scatter plot viewer"), "visible"));
    await session.step(23, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(24, "And the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
    await run.scenario("The legend lists the categories of the color column", async () => {
      await session.step(27, "Then \"Caucasian\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Caucasian\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(28, "And \"Asian\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(29, "And \"Caucasian\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Caucasian\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a category filters the viewer, not the table", async () => {
      await session.step(33, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
      await session.step(34, "Then \"Black\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(35, "And \"Black\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(36, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(37, "And scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
      await session.step(38, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(39, "And no rows should be selected", () => noneSelected(page));
      await session.step(40, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
      await session.step(41, "Then \"Black\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(42, "And \"Black\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(43, "And scatter plot viewer should show more rows than before", () => showsMoreRows(page, el("scatter plot viewer")));
      await session.step(44, "And scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
    });
    await run.scenario("Control-click adds a category and the cross takes one back", async () => {
      await session.step(47, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
      await session.step(48, "And user clicks on \"Asian\" item in the legend of scatter plot viewer holding Control", () => clickLegendItemHolding(page, "Asian", el("scatter plot viewer"), "Control"));
      await session.step(49, "Then \"Black\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(50, "And \"Asian\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(51, "And scatter plot viewer should show more rows than before", () => showsMoreRows(page, el("scatter plot viewer")));
      await session.step(52, "When user clicks on the cross of \"Asian\" item in the legend of scatter plot viewer", () => clickLegendCross(page, "Asian", el("scatter plot viewer")));
      await session.step(53, "Then \"Asian\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(54, "And \"Asian\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(55, "And \"Black\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(56, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(57, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
      await session.step(58, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    });
    await run.scenario("A click never moves the legend, and a filter never moves it either", async () => {
      await session.step(61, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(62, "And user clicks on \"Caucasian\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Caucasian", el("scatter plot viewer")));
      await session.step(63, "Then the legend of scatter plot viewer should be placed as before", () => legendPlacedAsBefore(page, el("scatter plot viewer")));
      await session.step(64, "And the legend of scatter plot viewer should list the same items as before", () => legendSameItems(page, el("scatter plot viewer")));
      await session.step(65, "When user clicks on \"Caucasian\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Caucasian", el("scatter plot viewer")));
      await session.step(66, "And user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(67, "And user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(68, "Then the legend of scatter plot viewer should be placed as before", () => legendPlacedAsBefore(page, el("scatter plot viewer")));
      await session.step(69, "And the legend of scatter plot viewer should list the same items as before", () => legendSameItems(page, el("scatter plot viewer")));
      await session.step(70, "When user resets the filter", () => resetFilter(page));
    });
    await run.scenario("A filter that empties a category drops it from the legend", async () => {
      await session.step(73, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(74, "And user filters rows where \"RACE\" is one of \"Caucasian, Other\"", () => filterToAnyOf(page, "RACE", "Caucasian, Other"));
      await session.step(75, "Then the legend of scatter plot viewer should list fewer items than before", () => legendFewer(page, el("scatter plot viewer")));
      await session.step(76, "And \"Asian\" legend item in legend of scatter plot viewer should be absent", () => shouldBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "absent"));
      await session.step(77, "And \"Caucasian\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Caucasian\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(78, "When user resets the filter", () => resetFilter(page));
      await session.step(79, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    });
    await run.scenario("The category colors are the column's, and a recoloring reaches the items", async () => {
      await session.step(82, "Then the \"Caucasian\" and \"Asian\" items in the legend of scatter plot viewer should be colored differently", () => legendItemsDiffer(page, "Caucasian", "Asian", el("scatter plot viewer")));
      await session.step(83, "When user colors \"RACE\" column categorically:", () => colorCategorical(page, "RACE", [["Caucasian","#FF0000"],["Asian","#0000FF"]]));
      await session.step(86, "Then the \"Caucasian\" item in the legend of scatter plot viewer should be colored \"#FF0000\"", () => legendItemColor(page, "Caucasian", el("scatter plot viewer"), "#FF0000"));
      await session.step(87, "And the \"Asian\" item in the legend of scatter plot viewer should be colored \"#0000FF\"", () => legendItemColor(page, "Asian", el("scatter plot viewer"), "#0000FF"));
      await session.step(88, "And the categorical color of \"Caucasian\" in \"RACE\" column should be \"#FF0000\"", () => categoricalColorIs(page, "Caucasian", "RACE", "#FF0000"));
      await session.step(89, "And the \"view\" area of scatter plot viewer should contain the color \"#FF0000\"", () => areaColor(page, "view", el("scatter plot viewer"), "#FF0000"));
      await session.step(90, "When user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(91, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
      await session.step(92, "And the \"Caucasian\" item in the legend of scatter plot viewer should be colored \"#FF0000\"", () => legendItemColor(page, "Caucasian", el("scatter plot viewer"), "#FF0000"));
    });
    await run.scenario("Position and visibility move and hide the legend", async () => {
      await session.step(95, "When user sets \"Legend Position\" property of scatter plot viewer to \"Left\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "Left"));
      await session.step(96, "Then the legend of scatter plot viewer should be on the left", () => legendSide(page, el("scatter plot viewer"), "left"));
      await session.step(97, "And the legend of scatter plot viewer should be in the \"left\" slot", () => legendSlot(page, el("scatter plot viewer"), "left"));
      await session.step(98, "And the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
      await session.step(99, "When user sets \"Legend Position\" property of scatter plot viewer to \"Bottom\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "Bottom"));
      await session.step(100, "Then the legend of scatter plot viewer should be on the bottom", () => legendSide(page, el("scatter plot viewer"), "bottom"));
      await session.step(101, "When user sets \"Legend Visibility\" property of scatter plot viewer to \"Never\"", () => setProperty(page, "Legend Visibility", el("scatter plot viewer"), "Never"));
      await session.step(102, "Then legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
      await session.step(103, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Legend Visibility","Always"],["Legend Position","Right"]]));
      await session.step(106, "Then the legend of scatter plot viewer should be on the right", () => legendSide(page, el("scatter plot viewer"), "right"));
      await session.step(107, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A numerical color column shows a scale instead of items", async () => {
      await session.step(111, "When user sets \"Color\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "AGE"));
      await session.step(112, "Then the legend of scatter plot viewer should list 0 items", () => legendLists(page, el("scatter plot viewer"), 0));
      await session.step(113, "And scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(114, "When user sets \"Color\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE"));
      await session.step(115, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
