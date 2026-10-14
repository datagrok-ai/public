/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/color-coding/color-coding-pick-up-apply-and-viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
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
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, categoricalColorIs, colorCategorical, colorCodedAs, colorCodedCategorically, colorOff, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, closeContextMenu, legendItemColor, noErrors, pickFromAreaContextMenu, readingAsRemembered, readingsDiffer, readingsEqual, rememberReading, rightClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pick Up / Apply Coloring from the header menu, and a column's coloring kept by the viewers around it", () => {
  const session = feature(test, "features/viewers/color-coding/color-coding-pick-up-apply-and-viewers.feature", import.meta.url);
  test("Pick Up / Apply Coloring from the header menu, and a column's coloring kept by the viewers around it", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a calculated column \"Race_copy\" with formula \"${RACE}\"", () => addCalculated(page, "Race_copy", "${RACE}"));
    await session.step(19, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
    await session.step(20, "And \"Race_copy\" column should have no color coding", () => noColorCoding(page, "Race_copy"));
    await run.scenario("Apply Coloring is disabled until a coloring has been picked up", async () => {
      await session.step(23, "When user right-clicks on the \"header Race_copy\" area of grid", () => rightClickArea(page, "header Race_copy", el("grid")));
      await session.step(24, "And user hovers over \"Color Coding\" menu item", () => hoverOver(page, el("\"Color Coding\" menu item")));
      await session.step(25, "Then \"Color-Coding---Apply-Coloring\" menu item should be disabled", () => shouldBe(page, el("\"Color-Coding---Apply-Coloring\" menu item"), "disabled"));
      await session.step(26, "And \"Color-Coding---Pick-Up-Coloring\" menu item should be enabled", () => shouldBe(page, el("\"Color-Coding---Pick-Up-Coloring\" menu item"), "enabled"));
      await session.step(27, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(28, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pick Up Coloring on one column and Apply Coloring on another give it the same colours", async () => {
      await session.step(31, "When user colors \"RACE\" column categorically:", () => colorCategorical(page, "RACE", [["Caucasian","#3366CC"],["Other","#CC6699"]]), [["Caucasian","#3366CC"],["Other","#CC6699"]]);
      await session.step(34, "Then \"RACE\" column should be color-coded categorically", () => colorCodedCategorically(page, "RACE"));
      await session.step(35, "When user picks \"Color Coding > Pick Up Coloring\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Pick Up Coloring", "header RACE", el("grid")));
      await session.step(36, "And user right-clicks on the \"header Race_copy\" area of grid", () => rightClickArea(page, "header Race_copy", el("grid")));
      await session.step(37, "And user hovers over \"Color Coding\" menu item", () => hoverOver(page, el("\"Color Coding\" menu item")));
      await session.step(38, "Then \"Color-Coding---Apply-Coloring\" menu item should be enabled", () => shouldBe(page, el("\"Color-Coding---Apply-Coloring\" menu item"), "enabled"));
      await session.step(39, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(40, "And user picks \"Color Coding > Apply Coloring\" from the context menu of the \"header Race_copy\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Apply Coloring", "header Race_copy", el("grid")));
      await session.step(41, "Then \"Race_copy\" column should be color-coded categorically", () => colorCodedCategorically(page, "Race_copy"));
      await session.step(42, "And the categorical color of \"Caucasian\" in \"Race_copy\" column should be \"3366CC\"", () => categoricalColorIs(page, "Caucasian", "Race_copy", "3366CC"));
      await session.step(43, "And the categorical color of \"Other\" in \"Race_copy\" column should be \"CC6699\"", () => categoricalColorIs(page, "Other", "Race_copy", "CC6699"));
      await session.step(44, "And the \"color of cell 1 of Race_copy\" and \"color of cell 1 of RACE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of Race_copy", "color of cell 1 of RACE", el("grid")));
      await session.step(45, "And the \"color of cell 2 of Race_copy\" and \"color of cell 2 of RACE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 2 of Race_copy", "color of cell 2 of RACE", el("grid")));
      await session.step(46, "And the \"color of cell 1 of Race_copy\" and \"color of cell 2 of Race_copy\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of Race_copy", "color of cell 2 of Race_copy", el("grid")));
      await session.step(47, "When user picks \"Color Coding > Off\" from the context menu of the \"header Race_copy\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header Race_copy", el("grid")));
      await session.step(48, "And user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(49, "Then \"Race_copy\" column should have no color coding", () => noColorCoding(page, "Race_copy"));
      await session.step(50, "And \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A heat map added to the view leaves the grid's colour coding as it was", async () => {
      await session.step(54, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(55, "And user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(56, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(57, "And \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(58, "When user remembers the \"color of cell 1 of AGE\" reading of grid", () => rememberReading(page, "color of cell 1 of AGE", el("grid")));
      await session.step(59, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
      await session.step(60, "Then heat map viewer should be visible", () => shouldBe(page, el("heat map viewer"), "visible"));
      await session.step(61, "And \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(62, "And \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(63, "And the \"color of cell 1 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 1 of AGE", el("grid")));
      await session.step(64, "When user clicks on close icon of heat map viewer", () => clickOn(page, el("close icon of heat map viewer")));
      await session.step(65, "Then heat map viewer should be hidden", () => shouldBe(page, el("heat map viewer"), "hidden"));
      await session.step(66, "And \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(67, "And \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(68, "And the \"color of cell 1 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 1 of AGE", el("grid")));
      await session.step(69, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(70, "And user picks \"Color Coding > Off\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header SEX", el("grid")));
      await session.step(71, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The categorical colours set on the grid column are the scatter plot legend's", async () => {
      await session.step(75, "When user colors \"SEX\" column categorically:", () => colorCategorical(page, "SEX", [["M","#3366CC"],["F","#CC6699"]]), [["M","#3366CC"],["F","#CC6699"]]);
      await session.step(78, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","SEX"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Color","SEX"]]);
      await session.step(82, "Then the \"M\" item in the legend of scatter plot viewer should be colored \"#3366CC\"", () => legendItemColor(page, "M", el("scatter plot viewer"), "#3366CC"));
      await session.step(83, "And the \"F\" item in the legend of scatter plot viewer should be colored \"#CC6699\"", () => legendItemColor(page, "F", el("scatter plot viewer"), "#CC6699"));
      await session.step(84, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
      await session.step(85, "And user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(86, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
