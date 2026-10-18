/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/legend/legend-interaction.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.legend]
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
import {shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, categoricalColorIs, colorCategorical, colorOff, filterPassesAll, filterTo, filterToAnyOf, noColorCoding, noneSelected, removeColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, clickLegendCross, clickLegendItem, clickLegendItemHolding, hasArea, hasNoArea, legendClearOfArea, legendDocked, legendFewer, legendInCorner, legendItemColor, legendItemsDiffer, legendLists, legendNarrower, legendPlacedAsBefore, legendSameItems, legendSide, legendSlot, legendWider, lessInk, moreInk, noErrors, repainted, setProperties, setProperty, showsFewerRows, showsMoreRows, showsRows, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Legend interaction", () => {
  const session = feature(test, "features/viewers/legend/legend-interaction.feature", import.meta.url);
  test("Legend interaction", {tag: ["@journey", "@viewers", "@realizes:viewers.legend", "@GROK-19041"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(28, "Then legend of scatter plot viewer should be visible", () => shouldBe(page, el("legend of scatter plot viewer"), "visible"));
    await session.step(29, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(30, "And the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
    await run.scenario("The legend lists the categories of the color column", async () => {
      await session.step(33, "Then \"Caucasian\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Caucasian\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(34, "And \"Asian\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(35, "And \"Caucasian\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Caucasian\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a category filters the viewer, not the table", async () => {
      await session.step(39, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
      await session.step(40, "Then \"Black\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(41, "And \"Black\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(42, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(43, "And scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
      await session.step(44, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(45, "And no rows should be selected", () => noneSelected(page));
      await session.step(46, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
      await session.step(47, "Then \"Black\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(48, "And \"Black\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(49, "And scatter plot viewer should show more rows than before", () => showsMoreRows(page, el("scatter plot viewer")));
      await session.step(50, "And scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
    });
    await run.scenario("Control-click adds a category and the cross takes one back", async () => {
      await session.step(53, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
      await session.step(54, "And user clicks on \"Asian\" item in the legend of scatter plot viewer holding Control", () => clickLegendItemHolding(page, "Asian", el("scatter plot viewer"), "Control"));
      await session.step(55, "Then \"Black\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(56, "And \"Asian\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(57, "And scatter plot viewer should show more rows than before", () => showsMoreRows(page, el("scatter plot viewer")));
      await session.step(58, "When user clicks on the cross of \"Asian\" item in the legend of scatter plot viewer", () => clickLegendCross(page, "Asian", el("scatter plot viewer")));
      await session.step(59, "Then \"Asian\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(60, "And \"Asian\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(61, "And \"Black\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(62, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(63, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
      await session.step(64, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    });
    await run.scenario("A click never moves the legend, and a filter never moves it either", async () => {
      await session.step(67, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(68, "And user clicks on \"Caucasian\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Caucasian", el("scatter plot viewer")));
      await session.step(69, "Then the legend of scatter plot viewer should be placed as before", () => legendPlacedAsBefore(page, el("scatter plot viewer")));
      await session.step(70, "And the legend of scatter plot viewer should list the same items as before", () => legendSameItems(page, el("scatter plot viewer")));
      await session.step(71, "When user clicks on \"Caucasian\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Caucasian", el("scatter plot viewer")));
      await session.step(72, "And user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(73, "And user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(74, "Then the legend of scatter plot viewer should be placed as before", () => legendPlacedAsBefore(page, el("scatter plot viewer")));
      await session.step(75, "And the legend of scatter plot viewer should list the same items as before", () => legendSameItems(page, el("scatter plot viewer")));
      await session.step(76, "When user resets the filter", () => resetFilter(page));
    });
    await run.scenario("A selection stays through a row filter and goes with the color column it was made on", async () => {
      await session.step(79, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
      await session.step(80, "Then \"Black\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(81, "And scatter plot viewer should show 25 rows", () => showsRows(page, el("scatter plot viewer"), 25));
      await session.step(82, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(83, "Then \"Black\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(84, "And scatter plot viewer should show 17 rows", () => showsRows(page, el("scatter plot viewer"), 17));
      await session.step(85, "When user resets the filter", () => resetFilter(page));
      await session.step(86, "Then scatter plot viewer should show 25 rows", () => showsRows(page, el("scatter plot viewer"), 25));
      await session.step(87, "When user sets \"Color\" property of scatter plot viewer to \"DIS_POP\"", () => setProperty(page, "Color", el("scatter plot viewer"), "DIS_POP"));
      await session.step(88, "Then the legend of scatter plot viewer should list 5 items", () => legendLists(page, el("scatter plot viewer"), 5));
      await session.step(89, "And \"Indigestion\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Indigestion\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(90, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(91, "When user sets \"Color\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE"));
      await session.step(92, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(93, "And \"Black\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(94, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter that empties a category drops it from the legend", async () => {
      await session.step(98, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(99, "And user filters rows where \"RACE\" is one of \"Caucasian, Other\"", () => filterToAnyOf(page, "RACE", "Caucasian, Other"));
      await session.step(100, "Then the legend of scatter plot viewer should list fewer items than before", () => legendFewer(page, el("scatter plot viewer")));
      await session.step(101, "And \"Asian\" legend item in legend of scatter plot viewer should be absent", () => shouldBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "absent"));
      await session.step(102, "And \"Caucasian\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Caucasian\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(103, "When user resets the filter", () => resetFilter(page));
      await session.step(104, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    });
    await run.scenario("The category colors are the column's, and a recoloring reaches the items", async () => {
      await session.step(107, "Then the \"Caucasian\" and \"Asian\" items in the legend of scatter plot viewer should be colored differently", () => legendItemsDiffer(page, "Caucasian", "Asian", el("scatter plot viewer")));
      await session.step(108, "When user colors \"RACE\" column categorically:", () => colorCategorical(page, "RACE", [["Caucasian","#FF0000"],["Asian","#0000FF"]]), [["Caucasian","#FF0000"],["Asian","#0000FF"]]);
      await session.step(111, "Then the \"Caucasian\" item in the legend of scatter plot viewer should be colored \"#FF0000\"", () => legendItemColor(page, "Caucasian", el("scatter plot viewer"), "#FF0000"));
      await session.step(112, "And the \"Asian\" item in the legend of scatter plot viewer should be colored \"#0000FF\"", () => legendItemColor(page, "Asian", el("scatter plot viewer"), "#0000FF"));
      await session.step(113, "And the categorical color of \"Caucasian\" in \"RACE\" column should be \"#FF0000\"", () => categoricalColorIs(page, "Caucasian", "RACE", "#FF0000"));
      await session.step(114, "And the \"view\" area of scatter plot viewer should contain the color \"#FF0000\"", () => areaColor(page, "view", el("scatter plot viewer"), "#FF0000"));
      await session.step(115, "When user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(116, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
      await session.step(117, "And the \"Caucasian\" item in the legend of scatter plot viewer should be colored \"#FF0000\"", () => legendItemColor(page, "Caucasian", el("scatter plot viewer"), "#FF0000"));
    });
    await run.scenario("Position and visibility move and hide the legend", async () => {
      await session.step(120, "When user sets \"Legend Position\" property of scatter plot viewer to \"Left\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "Left"));
      await session.step(121, "Then the legend of scatter plot viewer should be on the left", () => legendSide(page, el("scatter plot viewer"), "left"));
      await session.step(122, "And the legend of scatter plot viewer should be in the \"left\" slot", () => legendSlot(page, el("scatter plot viewer"), "left"));
      await session.step(123, "And the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
      await session.step(124, "When user sets \"Legend Position\" property of scatter plot viewer to \"Bottom\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "Bottom"));
      await session.step(125, "Then the legend of scatter plot viewer should be on the bottom", () => legendSide(page, el("scatter plot viewer"), "bottom"));
      await session.step(126, "When user sets \"Legend Visibility\" property of scatter plot viewer to \"Never\"", () => setProperty(page, "Legend Visibility", el("scatter plot viewer"), "Never"));
      await session.step(127, "Then legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
      await session.step(128, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Legend Visibility","Always"],["Legend Position","Right"]]), [["Legend Visibility","Always"],["Legend Position","Right"]]);
      await session.step(131, "Then the legend of scatter plot viewer should be on the right", () => legendSide(page, el("scatter plot viewer"), "right"));
      await session.step(132, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(133, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The docked legend is as wide as its labels need", async () => {
      await session.step(136, "When user adds a calculated column \"RACE described\" with formula \"\\\"race of the subject in the study: \\\" + ${RACE}\"", () => addCalculated(page, "RACE described", "\"race of the subject in the study: \" + ${RACE}"));
      await session.step(137, "And user sets \"Color\" property of scatter plot viewer to \"RACE described\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE described"));
      await session.step(138, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(139, "And the legend of scatter plot viewer should be wider than before", () => legendWider(page, el("scatter plot viewer")));
      await session.step(140, "And the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
      await session.step(141, "When user sets \"Color\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE"));
      await session.step(142, "Then the legend of scatter plot viewer should be narrower than before", () => legendNarrower(page, el("scatter plot viewer")));
      await session.step(143, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(144, "When user removes \"RACE described\" column", () => removeColumn(page, "RACE described"));
      await session.step(145, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A numerical color column shows a scale instead of items", async () => {
      await session.step(148, "When user sets \"Color\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "AGE"));
      await session.step(149, "Then the legend of scatter plot viewer should list 0 items", () => legendLists(page, el("scatter plot viewer"), 0));
      await session.step(150, "And scatter plot viewer should have a \"color scale\" area", () => hasArea(page, el("scatter plot viewer"), "color scale"));
      await session.step(151, "And scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(152, "When user sets \"Color\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE"));
      await session.step(153, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(154, "And scatter plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("scatter plot viewer"), "color scale"));
      await session.step(155, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An auto-placed legend leaves the corner the regression equation is drawn in", async () => {
      await session.step(159, "When user adds a calculated column \"CORNERS\" with formula \"if(${AGE} > 60, Mod(${AGE}, 2) * 100, 0)\"", () => addCalculated(page, "CORNERS", "if(${AGE} > 60, Mod(${AGE}, 2) * 100, 0)"));
      await session.step(160, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","AGE"],["Y","CORNERS"],["Legend Position","Auto"]]), [["X","AGE"],["Y","CORNERS"],["Legend Position","Auto"]]);
      await session.step(164, "Then the legend of scatter plot viewer should be in a corner", () => legendInCorner(page, el("scatter plot viewer")));
      await session.step(165, "And the legend of scatter plot viewer should be in the \"leftTop\" slot", () => legendSlot(page, el("scatter plot viewer"), "leftTop"));
      await session.step(166, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showRegressionLine","true"],["showRegressionLineEquation","true"]]), [["showRegressionLine","true"],["showRegressionLineEquation","true"]]);
      await session.step(169, "Then scatter plot viewer should have a \"regression stats\" area", () => hasArea(page, el("scatter plot viewer"), "regression stats"));
      await session.step(170, "And the legend of scatter plot viewer should not cover the \"regression stats\" area", () => legendClearOfArea(page, el("scatter plot viewer"), "regression stats"));
      await session.step(171, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(172, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["showRegressionLine","false"],["X","WEIGHT"],["Y","HEIGHT"],["Legend Position","Right"]]), [["showRegressionLine","false"],["X","WEIGHT"],["Y","HEIGHT"],["Legend Position","Right"]]);
      await session.step(177, "And user removes \"CORNERS\" column", () => removeColumn(page, "CORNERS"));
      await session.step(178, "Then the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
      await session.step(179, "And the legend of scatter plot viewer should be on the right", () => legendSide(page, el("scatter plot viewer"), "right"));
      await session.step(180, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
