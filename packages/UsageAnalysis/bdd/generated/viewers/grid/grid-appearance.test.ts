/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-appearance.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, isExpanded, pressKey, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnIncomplete, columnTag, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorCodedCategorically, colorConditional, colorOff, noColorCoding, noneSelected, rowsRangeSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaAtLeastTall, areaColor, areaNotColor, areasDiffer, areasSame, clickArea, dragAreaBy, dragSelectionBetweenAreas, noErrors, pickColorSwatch, pickFromAreaContextMenu, pointerAway, propertyShouldBe, readingAsRemembered, readingDoesNotRead, readingLower, readingReads, readingsDiffer, rememberReading, repainted, repaintedBy, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid cell appearance", () => {
  const session = feature(test, "features/viewers/grid/grid-appearance.feature", import.meta.url);
  test("Grid cell appearance", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(25, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(26, "And \"Color Coding\" property of grid should be \"Auto\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "Auto"));
    await session.step(27, "And \"HEIGHT\" column should have missing values", () => columnIncomplete(page, "HEIGHT"));
    await run.scenario("Linear coding from the header menu paints the column by value", async () => {
      await session.step(30, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(31, "Then \"AGE\" column should have tag \".color-coding-type\" equal to \"Linear\"", () => columnTag(page, "AGE", ".color-coding-type", "Linear"));
      await session.step(32, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(33, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(34, "Then the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(35, "And the \"cell 2 of AGE\" and \"cell 2 of USUBJID\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 2 of USUBJID", el("grid")));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The grid-wide Color Coding overrides the column's", async () => {
      await session.step(39, "When user picks \"Grid Color Coding > None\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > None", "cell 8 of AGE", el("grid")));
      await session.step(40, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(41, "Then \"Color Coding\" property of grid should be \"None\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "None"));
      await session.step(42, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in the same colors", () => areasSame(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(43, "When user picks \"Grid Color Coding > All\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > All", "cell 8 of AGE", el("grid")));
      await session.step(44, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(45, "Then \"Color Coding\" property of grid should be \"All\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "All"));
      await session.step(46, "And the \"cell 2 of HEIGHT\" and \"cell 4 of HEIGHT\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of HEIGHT", "cell 4 of HEIGHT", el("grid")));
      await session.step(47, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(48, "When user picks \"Grid Color Coding > Auto\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > Auto", "cell 8 of AGE", el("grid")));
      await session.step(49, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(50, "Then \"Color Coding\" property of grid should be \"Auto\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "Auto"));
      await session.step(51, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(52, "And the \"cell 2 of HEIGHT\" and \"cell 4 of HEIGHT\" areas of grid should be painted in the same colors", () => areasSame(page, "cell 2 of HEIGHT", "cell 4 of HEIGHT", el("grid")));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Narrowing a coloured column keeps its colour", async () => {
      await session.step(56, "When user picks \"Grid Color Coding > All\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > All", "cell 1 of AGE", el("grid")));
      await session.step(57, "And user remembers the \"color of cell 1 of AGE\" reading of grid", () => rememberReading(page, "color of cell 1 of AGE", el("grid")));
      await session.step(58, "And user drags the \"column resizer AGE\" area of grid by 30 pixels to the left", () => dragAreaBy(page, "column resizer AGE", el("grid"), 30, "left"));
      await session.step(59, "Then the \"column width of AGE\" reading of grid should be lower than before", () => readingLower(page, "column width of AGE", el("grid")));
      await session.step(60, "And the \"color of cell 1 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 1 of AGE", el("grid")));
      await session.step(61, "When user picks \"Column Sizing > Optimal\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Column Sizing > Optimal", "cell 1 of AGE", el("grid")));
      await session.step(62, "And user picks \"Grid Color Coding > Auto\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > Auto", "cell 1 of AGE", el("grid")));
      await session.step(63, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A custom format shows in the cell", async () => {
      await session.step(66, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(67, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(68, "When user enters \"0.00\" into Custom input in Format AGE dialog", () => enterInto(page, "0.00", el("Custom input in Format AGE dialog")));
      await session.step(69, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(70, "Then Format AGE dialog should be hidden", () => shouldBe(page, el("Format AGE dialog"), "hidden"));
      await session.step(71, "And \"AGE\" column should have tag \"format\" equal to \"0.00\"", () => columnTag(page, "AGE", "format", "0.00"));
      await session.step(72, "And the \"text of cell 1 of AGE\" reading of grid should be \"26.00\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26.00"));
      await session.step(73, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
      await session.step(74, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(75, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
      await session.step(76, "Given Details accordion header in context panel is expanded", () => isExpanded(page, el("Details accordion header in context panel")));
      await session.step(77, "Then \"format\" table row in context panel should contain text \"0.00\"", () => shouldContainText(page, el("\"format\" table row in context panel"), "0.00"));
      await session.step(78, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(79, "And user clears Custom input in Format AGE dialog", () => clearField(page, el("Custom input in Format AGE dialog")));
      await session.step(80, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(81, "Then the \"text of cell 1 of AGE\" reading of grid should be \"26\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Conditional coding paints the cells its ranges name", async () => {
      await session.step(85, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header HEIGHT", el("grid")));
      await session.step(86, "Then \"HEIGHT\" column should have tag \".color-coding-type\" equal to \"Conditional\"", () => columnTag(page, "HEIGHT", ".color-coding-type", "Conditional"));
      await session.step(87, "When user colors \"HEIGHT\" column conditionally:", () => colorConditional(page, "HEIGHT", [["<160","#0000FF"],[">180","#FF0000"]]), [["<160","#0000FF"],[">180","#FF0000"]]);
      await session.step(90, "Then the \"cell 2 of HEIGHT\" area of grid should contain the color \"#0000FF\"", () => areaColor(page, "cell 2 of HEIGHT", el("grid"), "#0000FF"));
      await session.step(91, "And the \"cell 4 of HEIGHT\" area of grid should contain the color \"#FF0000\"", () => areaColor(page, "cell 4 of HEIGHT", el("grid"), "#FF0000"));
      await session.step(92, "And the \"cell 1 of HEIGHT\" area of grid should not contain the color \"#0000FF\"", () => areaNotColor(page, "cell 1 of HEIGHT", el("grid"), "#0000FF"));
      await session.step(93, "And the \"cell 1 of HEIGHT\" area of grid should not contain the color \"#FF0000\"", () => areaNotColor(page, "cell 1 of HEIGHT", el("grid"), "#FF0000"));
      await session.step(94, "When user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
      await session.step(95, "Then \"HEIGHT\" column should have no color coding", () => noColorCoding(page, "HEIGHT"));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Categorical coding tells the categories apart", async () => {
      await session.step(99, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(100, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(101, "And the \"cell 1 of SEX\" and \"cell 4 of SEX\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 1 of SEX", "cell 4 of SEX", el("grid")));
      await session.step(102, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(103, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(104, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row height and the missing-value colour paint through the grid's own properties", async () => {
      await session.step(107, "When user sets \"Row Height\" property of grid to \"48\"", () => setProperty(page, "Row Height", el("grid"), "48"));
      await session.step(108, "Then the \"cell 1 of AGE\" area of grid should be at least 40 pixels tall", () => areaAtLeastTall(page, "cell 1 of AGE", el("grid"), 40));
      await session.step(109, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(110, "When user sets \"Row Height\" property of grid to \"28\"", () => setProperty(page, "Row Height", el("grid"), "28"));
      await session.step(111, "And user makes row 298 current", () => makeRowCurrent(page, 298));
      await session.step(112, "And user makes row 300 current", () => makeRowCurrent(page, 300));
      await session.step(113, "And user sets \"Missing Value Color\" property of grid to \"#FFAAAA\"", () => setProperty(page, "Missing Value Color", el("grid"), "#FFAAAA"));
      await session.step(114, "Then the \"cell 298 of HEIGHT\" area of grid should contain the color \"#FFAAAA\"", () => areaColor(page, "cell 298 of HEIGHT", el("grid"), "#FFAAAA"));
      await session.step(115, "When user sets \"Missing Value Color\" property of grid to \"#00FF00\"", () => setProperty(page, "Missing Value Color", el("grid"), "#00FF00"));
      await session.step(116, "Then the \"cell 298 of HEIGHT\" area of grid should contain the color \"#00FF00\"", () => areaColor(page, "cell 298 of HEIGHT", el("grid"), "#00FF00"));
      await session.step(117, "And the \"cell 298 of HEIGHT\" area of grid should not contain the color \"#FFAAAA\"", () => areaNotColor(page, "cell 298 of HEIGHT", el("grid"), "#FFAAAA"));
      await session.step(118, "When user sets \"Missing Value Color\" property of grid to \"#FFFFFF\"", () => setProperty(page, "Missing Value Color", el("grid"), "#FFFFFF"));
      await session.step(119, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(120, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A bigger cell font repaints the grid", async () => {
      await session.step(123, "When user sets \"Default Cell Font\" property of grid to \"20px Roboto\"", () => setProperty(page, "Default Cell Font", el("grid"), "20px Roboto"));
      await session.step(124, "Then grid should have repainted by at least 3000 pixels", () => repaintedBy(page, el("grid"), 3000));
      await session.step(125, "When user sets \"Default Cell Font\" property of grid to \"12px Roboto\"", () => setProperty(page, "Default Cell Font", el("grid"), "12px Roboto"));
      await session.step(126, "Then grid should have repainted", () => repainted(page, el("grid")));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Selected Rows Color paints the selected rows and Escape takes it away", async () => {
      await session.step(130, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(131, "Then no rows should be selected", () => noneSelected(page));
      await session.step(132, "When user sets \"Selected Rows Color\" property of grid to \"#00FF00\"", () => setProperty(page, "Selected Rows Color", el("grid"), "#00FF00"));
      await session.step(133, "Then the \"cell 6 of USUBJID\" area of grid should not contain the color \"#00FF00\"", () => areaNotColor(page, "cell 6 of USUBJID", el("grid"), "#00FF00"));
      await session.step(134, "When user drags a selection box from the \"row header 5\" area to the \"row header 7\" area of grid", () => dragSelectionBetweenAreas(page, "row header 5", "row header 7", el("grid")));
      await session.step(135, "Then rows 5 to 7 should be selected", () => rowsRangeSelected(page, 5, 7));
      await session.step(136, "And the \"cell 6 of USUBJID\" area of grid should contain the color \"#00FF00\"", () => areaColor(page, "cell 6 of USUBJID", el("grid"), "#00FF00"));
      await session.step(137, "And the \"cell 9 of USUBJID\" area of grid should not contain the color \"#00FF00\"", () => areaNotColor(page, "cell 9 of USUBJID", el("grid"), "#00FF00"));
      await session.step(138, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(139, "Then no rows should be selected", () => noneSelected(page));
      await session.step(140, "And the \"cell 6 of USUBJID\" area of grid should not contain the color \"#00FF00\"", () => areaNotColor(page, "cell 6 of USUBJID", el("grid"), "#00FF00"));
      await session.step(141, "When user sets \"Selected Rows Color\" property of grid to \"819780688\"", () => setProperty(page, "Selected Rows Color", el("grid"), "819780688"));
      await session.step(142, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column's colour coding wins over the background its Style sets (GROK-18638)", async () => {
      await session.step(145, "When user clicks on the \"header SEVERITY\" area of grid", () => clickArea(page, "header SEVERITY", el("grid")));
      await session.step(146, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(147, "Then the context panel should show \"SEVERITY\"", () => contextPanelShows(page, "SEVERITY"));
      await session.step(148, "Given Style accordion header in context panel is expanded", () => isExpanded(page, el("Style accordion header in context panel")));
      await session.step(149, "And Content accordion header in context panel is expanded", () => isExpanded(page, el("Content accordion header in context panel")));
      await session.step(150, "When user clicks on editor of \"Back Color\" property in context panel", () => clickOn(page, el("editor of \"Back Color\" property in context panel")));
      await session.step(151, "And user picks the color \"#FFA500\" in the color picker", () => pickColorSwatch(page, "#FFA500"));
      await session.step(152, "Then the \"color of cell 1 of SEVERITY\" reading of grid should be \"#ffa500\"", () => readingReads(page, "color of cell 1 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(153, "And the \"color of cell 2 of SEVERITY\" reading of grid should be \"#ffa500\"", () => readingReads(page, "color of cell 2 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(154, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEVERITY\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEVERITY", el("grid")));
      await session.step(155, "Then the \"color of cell 1 of SEVERITY\" reading of grid should not be \"#ffa500\"", () => readingDoesNotRead(page, "color of cell 1 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(156, "And the \"color of cell 1 of SEVERITY\" and \"color of cell 3 of SEVERITY\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEVERITY", "color of cell 3 of SEVERITY", el("grid")));
      await session.step(157, "When user removes the coloring of \"SEVERITY\" column", () => colorOff(page, "SEVERITY"));
      await session.step(158, "Then the \"color of cell 1 of SEVERITY\" reading of grid should be \"#ffa500\"", () => readingReads(page, "color of cell 1 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(159, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
