/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-appearance.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnIncomplete, columnTag, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorCodedCategorically, colorConditional, colorOff, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaAtLeastTall, areaColor, areaNotColor, areasDiffer, areasSame, dragAreaBy, noErrors, notRepainted, pickFromAreaContextMenu, pointerAway, propertyShouldBe, readingAsRemembered, readingLower, readingReads, rememberReading, repainted, repaintedBy, setProperty, showsRows, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid cell appearance", () => {
  const session = feature(test, "features/viewers/grid/grid-appearance.feature", import.meta.url);
  test("Grid cell appearance", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(15, "And \"Color Coding\" property of grid should be \"Auto\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "Auto"));
    await session.step(16, "And \"HEIGHT\" column should have missing values", () => columnIncomplete(page, "HEIGHT"));
    await run.scenario("Linear coding from the header menu paints the column by value", async () => {
      await session.step(19, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(20, "Then \"AGE\" column should have tag \".color-coding-type\" equal to \"Linear\"", () => columnTag(page, "AGE", ".color-coding-type", "Linear"));
      await session.step(21, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(22, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(23, "Then the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(24, "And the \"cell 2 of AGE\" and \"cell 2 of USUBJID\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 2 of USUBJID", el("grid")));
      await session.step(25, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The grid-wide Color Coding overrides the column's", async () => {
      await session.step(28, "When user picks \"Grid Color Coding > None\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > None", "cell 8 of AGE", el("grid")));
      await session.step(29, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(30, "Then \"Color Coding\" property of grid should be \"None\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "None"));
      await session.step(31, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in the same colors", () => areasSame(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(32, "When user picks \"Grid Color Coding > All\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > All", "cell 8 of AGE", el("grid")));
      await session.step(33, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(34, "Then \"Color Coding\" property of grid should be \"All\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "All"));
      await session.step(35, "And the \"cell 2 of HEIGHT\" and \"cell 4 of HEIGHT\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of HEIGHT", "cell 4 of HEIGHT", el("grid")));
      await session.step(36, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(37, "When user picks \"Grid Color Coding > Auto\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > Auto", "cell 8 of AGE", el("grid")));
      await session.step(38, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(39, "Then \"Color Coding\" property of grid should be \"Auto\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "Auto"));
      await session.step(40, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(41, "And the \"cell 2 of HEIGHT\" and \"cell 4 of HEIGHT\" areas of grid should be painted in the same colors", () => areasSame(page, "cell 2 of HEIGHT", "cell 4 of HEIGHT", el("grid")));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Narrowing a coloured column keeps its colour and its text", async () => {
      await session.step(45, "When user picks \"Grid Color Coding > All\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > All", "cell 1 of AGE", el("grid")));
      await session.step(46, "And user remembers the \"color of cell 1 of AGE\" reading of grid", () => rememberReading(page, "color of cell 1 of AGE", el("grid")));
      await session.step(47, "And user drags the \"column resizer AGE\" area of grid by 30 pixels to the left", () => dragAreaBy(page, "column resizer AGE", el("grid"), 30, "left"));
      await session.step(48, "Then the \"column width of AGE\" reading of grid should be lower than before", () => readingLower(page, "column width of AGE", el("grid")));
      await session.step(49, "And the \"color of cell 1 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 1 of AGE", el("grid")));
      await session.step(50, "And the \"text of cell 1 of AGE\" reading of grid should be \"26\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26"));
      await session.step(51, "When user picks \"Column Sizing > Optimal\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Column Sizing > Optimal", "cell 1 of AGE", el("grid")));
      await session.step(52, "And user picks \"Grid Color Coding > Auto\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > Auto", "cell 1 of AGE", el("grid")));
      await session.step(53, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A custom format shows in the cell", async () => {
      await session.step(56, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(57, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(58, "When user enters \"0.00\" into Custom input in Format AGE dialog", () => enterInto(page, "0.00", el("Custom input in Format AGE dialog")));
      await session.step(59, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(60, "Then Format AGE dialog should be hidden", () => shouldBe(page, el("Format AGE dialog"), "hidden"));
      await session.step(61, "And \"AGE\" column should have tag \"format\" equal to \"0.00\"", () => columnTag(page, "AGE", "format", "0.00"));
      await session.step(62, "And the \"text of cell 1 of AGE\" reading of grid should be \"26.00\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26.00"));
      await session.step(63, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(64, "And user clears Custom input in Format AGE dialog", () => clearField(page, el("Custom input in Format AGE dialog")));
      await session.step(65, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(66, "Then the \"text of cell 1 of AGE\" reading of grid should be \"26\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26"));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Conditional coding paints the cells its ranges name", async () => {
      await session.step(70, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header HEIGHT", el("grid")));
      await session.step(71, "Then \"HEIGHT\" column should have tag \".color-coding-type\" equal to \"Conditional\"", () => columnTag(page, "HEIGHT", ".color-coding-type", "Conditional"));
      await session.step(72, "When user colors \"HEIGHT\" column conditionally:", () => colorConditional(page, "HEIGHT", [["<160","#0000FF"],[">180","#FF0000"]]));
      await session.step(75, "Then the \"cell 2 of HEIGHT\" area of grid should contain the color \"#0000FF\"", () => areaColor(page, "cell 2 of HEIGHT", el("grid"), "#0000FF"));
      await session.step(76, "And the \"cell 4 of HEIGHT\" area of grid should contain the color \"#FF0000\"", () => areaColor(page, "cell 4 of HEIGHT", el("grid"), "#FF0000"));
      await session.step(77, "And the \"cell 1 of HEIGHT\" area of grid should not contain the color \"#0000FF\"", () => areaNotColor(page, "cell 1 of HEIGHT", el("grid"), "#0000FF"));
      await session.step(78, "And the \"cell 1 of HEIGHT\" area of grid should not contain the color \"#FF0000\"", () => areaNotColor(page, "cell 1 of HEIGHT", el("grid"), "#FF0000"));
      await session.step(79, "When user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
      await session.step(80, "Then \"HEIGHT\" column should have no color coding", () => noColorCoding(page, "HEIGHT"));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Categorical coding tells the categories apart", async () => {
      await session.step(84, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(85, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(86, "And the \"cell 1 of SEX\" and \"cell 4 of SEX\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 1 of SEX", "cell 4 of SEX", el("grid")));
      await session.step(87, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(88, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row height and the missing-value colour paint through the grid's own properties", async () => {
      await session.step(92, "When user sets \"Row Height\" property of grid to \"48\"", () => setProperty(page, "Row Height", el("grid"), "48"));
      await session.step(93, "Then the \"cell 1 of AGE\" area of grid should be at least 40 pixels tall", () => areaAtLeastTall(page, "cell 1 of AGE", el("grid"), 40));
      await session.step(94, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(95, "When user sets \"Row Height\" property of grid to \"28\"", () => setProperty(page, "Row Height", el("grid"), "28"));
      await session.step(96, "And user makes row 298 current", () => makeRowCurrent(page, 298));
      await session.step(97, "And user makes row 300 current", () => makeRowCurrent(page, 300));
      await session.step(98, "And user sets \"Missing Value Color\" property of grid to \"#FFAAAA\"", () => setProperty(page, "Missing Value Color", el("grid"), "#FFAAAA"));
      await session.step(99, "Then the \"cell 298 of HEIGHT\" area of grid should contain the color \"#FFAAAA\"", () => areaColor(page, "cell 298 of HEIGHT", el("grid"), "#FFAAAA"));
      await session.step(100, "When user sets \"Missing Value Color\" property of grid to \"#00FF00\"", () => setProperty(page, "Missing Value Color", el("grid"), "#00FF00"));
      await session.step(101, "Then the \"cell 298 of HEIGHT\" area of grid should contain the color \"#00FF00\"", () => areaColor(page, "cell 298 of HEIGHT", el("grid"), "#00FF00"));
      await session.step(102, "And the \"cell 298 of HEIGHT\" area of grid should not contain the color \"#FFAAAA\"", () => areaNotColor(page, "cell 298 of HEIGHT", el("grid"), "#FFAAAA"));
      await session.step(103, "When user sets \"Missing Value Color\" property of grid to \"#FFFFFF\"", () => setProperty(page, "Missing Value Color", el("grid"), "#FFFFFF"));
      await session.step(104, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(105, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A bigger cell font repaints the grid and an idle grid does not", async () => {
      await session.step(108, "When user takes a snapshot of grid", () => takeSnapshot(page, el("grid")));
      await session.step(109, "Then grid should not have repainted", () => notRepainted(page, el("grid")));
      await session.step(110, "When user sets \"Default Cell Font\" property of grid to \"20px Roboto\"", () => setProperty(page, "Default Cell Font", el("grid"), "20px Roboto"));
      await session.step(111, "Then grid should have repainted by at least 3000 pixels", () => repaintedBy(page, el("grid"), 3000));
      await session.step(112, "When user sets \"Default Cell Font\" property of grid to \"12px Roboto\"", () => setProperty(page, "Default Cell Font", el("grid"), "12px Roboto"));
      await session.step(113, "Then grid should have repainted", () => repainted(page, el("grid")));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
