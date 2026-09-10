/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-selector-and-legend.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, rightClickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {categoricalColorIs, colorOff, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaNotColor, hasArea, hasNoArea, legendItemColor, legendItemsDiffer, legendLists, legendSide, noErrors, pickColorSwatch, propertyShouldBe, readingIs, readingLower, repaintedBy, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart on-chart selector and legend", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-selector-and-legend.feature", import.meta.url);
  test("Pie chart on-chart selector and legend", {tag: ["@journey", "@viewers", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Legend Visibility","Always"]]));
    await session.step(21, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(22, "And legend of pie chart viewer should be visible", () => shouldBe(page, el("legend of pie chart viewer"), "visible"));
    await run.scenario("The legend lists exactly the category column's categories", async () => {
      await session.step(25, "Then the legend of pie chart viewer should list 4 items", () => legendLists(page, el("pie chart viewer"), 4));
      await session.step(26, "And \"Asian\" legend item in legend of pie chart viewer should be visible", () => shouldBe(page, el("\"Asian\" legend item in legend of pie chart viewer"), "visible"));
      await session.step(27, "And \"Black\" legend item in legend of pie chart viewer should be visible", () => shouldBe(page, el("\"Black\" legend item in legend of pie chart viewer"), "visible"));
      await session.step(28, "And \"Caucasian\" legend item in legend of pie chart viewer should be visible", () => shouldBe(page, el("\"Caucasian\" legend item in legend of pie chart viewer"), "visible"));
      await session.step(29, "And \"Other\" legend item in legend of pie chart viewer should be visible", () => shouldBe(page, el("\"Other\" legend item in legend of pie chart viewer"), "visible"));
      await session.step(30, "And the \"Asian\" and \"Caucasian\" items in the legend of pie chart viewer should be colored differently", () => legendItemsDiffer(page, "Asian", "Caucasian", el("pie chart viewer")));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Picking SEX in the on-chart selector re-splits the disc and the legend", async () => {
      await session.step(34, "When user hovers over pie chart viewer", () => hoverOver(page, el("pie chart viewer")));
      await session.step(35, "And user picks \"SEX\" in the \"category\" column selector of pie chart viewer", () => pickInColumnSelector(page, "SEX", "category", el("pie chart viewer")));
      await session.step(36, "Then \"Category\" property of pie chart viewer should be \"SEX\"", () => propertyShouldBe(page, "Category", el("pie chart viewer"), "SEX"));
      await session.step(37, "And the \"slices\" reading of pie chart viewer should be 2", () => readingIs(page, "slices", el("pie chart viewer"), 2));
      await session.step(38, "And pie chart viewer should have a \"slice F\" area", () => hasArea(page, el("pie chart viewer"), "slice F"));
      await session.step(39, "And pie chart viewer should have a \"slice M\" area", () => hasArea(page, el("pie chart viewer"), "slice M"));
      await session.step(40, "And the \"angle value of F\" reading of pie chart viewer should be 553", () => readingIs(page, "angle value of F", el("pie chart viewer"), 553));
      await session.step(41, "And the \"angle value of M\" reading of pie chart viewer should be 447", () => readingIs(page, "angle value of M", el("pie chart viewer"), 447));
      await session.step(42, "And the legend of pie chart viewer should list 2 items", () => legendLists(page, el("pie chart viewer"), 2));
      await session.step(43, "And pie chart viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pie chart viewer"), 1000));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Picking RACE back restores it", async () => {
      await session.step(47, "When user hovers over pie chart viewer", () => hoverOver(page, el("pie chart viewer")));
      await session.step(48, "And user picks \"RACE\" in the \"category\" column selector of pie chart viewer", () => pickInColumnSelector(page, "RACE", "category", el("pie chart viewer")));
      await session.step(49, "Then \"Category\" property of pie chart viewer should be \"RACE\"", () => propertyShouldBe(page, "Category", el("pie chart viewer"), "RACE"));
      await session.step(50, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(51, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 896", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 896));
      await session.step(52, "And pie chart viewer should not have a \"slice F\" area", () => hasNoArea(page, el("pie chart viewer"), "slice F"));
      await session.step(53, "And the legend of pie chart viewer should list 4 items", () => legendLists(page, el("pie chart viewer"), 4));
      await session.step(54, "And pie chart viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pie chart viewer"), 1000));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Position lays the legend out on each of the four sides", async () => {
      await session.step(58, "When user sets \"Legend Position\" property of pie chart viewer to \"Left\"", () => setProperty(page, "Legend Position", el("pie chart viewer"), "Left"));
      await session.step(59, "Then the legend of pie chart viewer should be on the left", () => legendSide(page, el("pie chart viewer"), "left"));
      await session.step(60, "And the \"pie radius\" reading of pie chart viewer should be lower than before", () => readingLower(page, "pie radius", el("pie chart viewer")));
      await session.step(61, "When user sets \"Legend Position\" property of pie chart viewer to \"Right\"", () => setProperty(page, "Legend Position", el("pie chart viewer"), "Right"));
      await session.step(62, "Then the legend of pie chart viewer should be on the right", () => legendSide(page, el("pie chart viewer"), "right"));
      await session.step(63, "When user sets \"Legend Position\" property of pie chart viewer to \"Top\"", () => setProperty(page, "Legend Position", el("pie chart viewer"), "Top"));
      await session.step(64, "Then the legend of pie chart viewer should be on the top", () => legendSide(page, el("pie chart viewer"), "top"));
      await session.step(65, "When user sets \"Legend Position\" property of pie chart viewer to \"Bottom\"", () => setProperty(page, "Legend Position", el("pie chart viewer"), "Bottom"));
      await session.step(66, "Then the legend of pie chart viewer should be on the bottom", () => legendSide(page, el("pie chart viewer"), "bottom"));
      await session.step(67, "And legend of pie chart viewer should be visible", () => shouldBe(page, el("legend of pie chart viewer"), "visible"));
      await session.step(68, "When user sets \"Legend Position\" property of pie chart viewer to \"Auto\"", () => setProperty(page, "Legend Position", el("pie chart viewer"), "Auto"));
      await session.step(69, "Then legend of pie chart viewer should be visible", () => shouldBe(page, el("legend of pie chart viewer"), "visible"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Visibility Never takes it away and Always brings it back", async () => {
      await session.step(73, "Then legend of pie chart viewer should be visible", () => shouldBe(page, el("legend of pie chart viewer"), "visible"));
      await session.step(74, "When user sets \"Legend Visibility\" property of pie chart viewer to \"Never\"", () => setProperty(page, "Legend Visibility", el("pie chart viewer"), "Never"));
      await session.step(75, "Then legend of pie chart viewer should be hidden", () => shouldBe(page, el("legend of pie chart viewer"), "hidden"));
      await session.step(76, "When user sets \"Legend Visibility\" property of pie chart viewer to \"Always\"", () => setProperty(page, "Legend Visibility", el("pie chart viewer"), "Always"));
      await session.step(77, "Then legend of pie chart viewer should be visible", () => shouldBe(page, el("legend of pie chart viewer"), "visible"));
      await session.step(78, "And the legend of pie chart viewer should list 4 items", () => legendLists(page, el("pie chart viewer"), 4));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A colour picked on a legend item reaches the wedge and the column", async () => {
      await session.step(82, "Given the \"pie\" area of pie chart viewer should not contain the color \"#9467BD\"", () => areaNotColor(page, "pie", el("pie chart viewer"), "#9467BD"));
      await session.step(83, "And the \"Asian\" item in the legend of pie chart viewer should be colored \"#1F77B4\"", () => legendItemColor(page, "Asian", el("pie chart viewer"), "#1F77B4"));
      await session.step(84, "When user right-clicks on \"Asian\" legend item in legend of pie chart viewer", () => rightClickOn(page, el("\"Asian\" legend item in legend of pie chart viewer")));
      await session.step(85, "Then \"Asian\" dialog should be visible", () => shouldBe(page, el("\"Asian\" dialog"), "visible"));
      await session.step(86, "When user picks the color \"#9467BD\" in the color picker dialog", () => pickColorSwatch(page, "#9467BD"));
      await session.step(87, "And user clicks on OK button in \"Asian\" dialog", () => clickOn(page, el("OK button in \"Asian\" dialog")));
      await session.step(88, "Then the categorical color of \"Asian\" in \"RACE\" column should be \"#9467BD\"", () => categoricalColorIs(page, "Asian", "RACE", "#9467BD"));
      await session.step(89, "And the \"Asian\" item in the legend of pie chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("pie chart viewer"), "#9467BD"));
      await session.step(90, "And the \"pie\" area of pie chart viewer should contain the color \"#9467BD\"", () => areaColor(page, "pie", el("pie chart viewer"), "#9467BD"));
      await session.step(91, "And the \"pie\" area of pie chart viewer should not contain the color \"#1F77B4\"", () => areaNotColor(page, "pie", el("pie chart viewer"), "#1F77B4"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Picking the default colour back takes the wedge back", async () => {
      await session.step(95, "When user right-clicks on \"Asian\" legend item in legend of pie chart viewer", () => rightClickOn(page, el("\"Asian\" legend item in legend of pie chart viewer")));
      await session.step(96, "And user picks the color \"#1F77B4\" in the color picker dialog", () => pickColorSwatch(page, "#1F77B4"));
      await session.step(97, "And user clicks on OK button in \"Asian\" dialog", () => clickOn(page, el("OK button in \"Asian\" dialog")));
      await session.step(98, "Then the categorical color of \"Asian\" in \"RACE\" column should be \"#1F77B4\"", () => categoricalColorIs(page, "Asian", "RACE", "#1F77B4"));
      await session.step(99, "And the \"Asian\" item in the legend of pie chart viewer should be colored \"#1F77B4\"", () => legendItemColor(page, "Asian", el("pie chart viewer"), "#1F77B4"));
      await session.step(100, "And the \"pie\" area of pie chart viewer should contain the color \"#1F77B4\"", () => areaColor(page, "pie", el("pie chart viewer"), "#1F77B4"));
      await session.step(101, "And the \"pie\" area of pie chart viewer should not contain the color \"#9467BD\"", () => areaNotColor(page, "pie", el("pie chart viewer"), "#9467BD"));
      await session.step(102, "When user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(103, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
      await session.step(104, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
