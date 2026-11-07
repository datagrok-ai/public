/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-appearance.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#1200]
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
import {clearField, clickOn, enterInto, isExpanded, pressKey, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnIncomplete, columnTag, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorCodedCategorically, colorConditional, colorOff, noColorCoding, noneSelected, rowsRangeSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaAtLeastTall, areaColor, areaNotColor, areasDiffer, areasSame, clickArea, dragAreaBy, dragSelectionBetweenAreas, hoverArea, noErrors, pickColorSwatch, pickFromAreaContextMenu, pointerAway, propertyShouldBe, readingAsRemembered, readingDoesNotRead, readingLower, readingReads, readingsDiffer, rememberReading, repainted, repaintedBy, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid cell appearance", () => {
  const session = feature(test, "features/viewers/grid/grid-appearance.feature", import.meta.url);
  test("Grid cell appearance", {tag: ["@journey", "@viewers", "@realizes:viewers.grid", "@edge-cases", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 12, page);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(30, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(31, "And \"Color Coding\" property of grid should be \"Auto\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "Auto"));
    await session.step(32, "And \"HEIGHT\" column should have missing values", () => columnIncomplete(page, "HEIGHT"));
    await run.scenario("Linear coding from the header menu paints the column by value", async () => {
      await session.step(35, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(36, "Then \"AGE\" column should have tag \".color-coding-type\" equal to \"Linear\"", () => columnTag(page, "AGE", ".color-coding-type", "Linear"));
      await session.step(37, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(38, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(39, "Then the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(40, "And the \"cell 2 of AGE\" and \"cell 2 of USUBJID\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 2 of USUBJID", el("grid")));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The grid-wide Color Coding overrides the column's", async () => {
      await session.step(44, "When user picks \"Grid Color Coding > None\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > None", "cell 8 of AGE", el("grid")));
      await session.step(45, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(46, "Then \"Color Coding\" property of grid should be \"None\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "None"));
      await session.step(47, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in the same colors", () => areasSame(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(48, "When user picks \"Grid Color Coding > All\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > All", "cell 8 of AGE", el("grid")));
      await session.step(49, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(50, "Then \"Color Coding\" property of grid should be \"All\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "All"));
      await session.step(51, "And the \"cell 2 of HEIGHT\" and \"cell 4 of HEIGHT\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of HEIGHT", "cell 4 of HEIGHT", el("grid")));
      await session.step(52, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(53, "When user picks \"Grid Color Coding > Auto\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > Auto", "cell 8 of AGE", el("grid")));
      await session.step(54, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(55, "Then \"Color Coding\" property of grid should be \"Auto\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "Auto"));
      await session.step(56, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(57, "And the \"cell 2 of HEIGHT\" and \"cell 4 of HEIGHT\" areas of grid should be painted in the same colors", () => areasSame(page, "cell 2 of HEIGHT", "cell 4 of HEIGHT", el("grid")));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Narrowing a coloured column keeps its colour", async () => {
      await session.step(61, "When user picks \"Grid Color Coding > All\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > All", "cell 1 of AGE", el("grid")));
      await session.step(62, "And user remembers the \"color of cell 1 of AGE\" reading of grid", () => rememberReading(page, "color of cell 1 of AGE", el("grid")));
      await session.step(63, "And user drags the \"column resizer AGE\" area of grid by 30 pixels to the left", () => dragAreaBy(page, "column resizer AGE", el("grid"), 30, "left"));
      await session.step(64, "Then the \"column width of AGE\" reading of grid should be lower than before", () => readingLower(page, "column width of AGE", el("grid")));
      await session.step(65, "And the \"color of cell 1 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 1 of AGE", el("grid")));
      await session.step(66, "When user picks \"Column Sizing > Optimal\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Column Sizing > Optimal", "cell 1 of AGE", el("grid")));
      await session.step(67, "And user picks \"Grid Color Coding > Auto\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > Auto", "cell 1 of AGE", el("grid")));
      await session.step(68, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A custom format shows in the cell", async () => {
      await session.step(71, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(72, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(73, "When user enters \"0.00\" into Custom input in Format AGE dialog", () => enterInto(page, "0.00", el("Custom input in Format AGE dialog")));
      await session.step(74, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(75, "Then Format AGE dialog should be hidden", () => shouldBe(page, el("Format AGE dialog"), "hidden"));
      await session.step(76, "And \"AGE\" column should have tag \"format\" equal to \"0.00\"", () => columnTag(page, "AGE", "format", "0.00"));
      await session.step(77, "And the \"text of cell 1 of AGE\" reading of grid should be \"26.00\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26.00"));
      await session.step(78, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
      await session.step(79, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(80, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
      await session.step(81, "Given Details accordion header in context panel is expanded", () => isExpanded(page, el("Details accordion header in context panel")));
      await session.step(82, "Then \"format\" table row in context panel should contain text \"0.00\"", () => shouldContainText(page, el("\"format\" table row in context panel"), "0.00"));
      await session.step(83, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(84, "And user clears Custom input in Format AGE dialog", () => clearField(page, el("Custom input in Format AGE dialog")));
      await session.step(85, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(86, "Then the \"text of cell 1 of AGE\" reading of grid should be \"26\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26"));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Conditional coding paints the cells its ranges name", async () => {
      await session.step(90, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header HEIGHT", el("grid")));
      await session.step(91, "Then \"HEIGHT\" column should have tag \".color-coding-type\" equal to \"Conditional\"", () => columnTag(page, "HEIGHT", ".color-coding-type", "Conditional"));
      await session.step(92, "When user colors \"HEIGHT\" column conditionally:", () => colorConditional(page, "HEIGHT", [["<160","#0000FF"],[">180","#FF0000"]]), [["<160","#0000FF"],[">180","#FF0000"]]);
      await session.step(95, "Then the \"cell 2 of HEIGHT\" area of grid should contain the color \"#0000FF\"", () => areaColor(page, "cell 2 of HEIGHT", el("grid"), "#0000FF"));
      await session.step(96, "And the \"cell 4 of HEIGHT\" area of grid should contain the color \"#FF0000\"", () => areaColor(page, "cell 4 of HEIGHT", el("grid"), "#FF0000"));
      await session.step(97, "And the \"cell 1 of HEIGHT\" area of grid should not contain the color \"#0000FF\"", () => areaNotColor(page, "cell 1 of HEIGHT", el("grid"), "#0000FF"));
      await session.step(98, "And the \"cell 1 of HEIGHT\" area of grid should not contain the color \"#FF0000\"", () => areaNotColor(page, "cell 1 of HEIGHT", el("grid"), "#FF0000"));
      await session.step(99, "When user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
      await session.step(100, "Then \"HEIGHT\" column should have no color coding", () => noColorCoding(page, "HEIGHT"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Categorical coding tells the categories apart", async () => {
      await session.step(104, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(105, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(106, "And the \"cell 1 of SEX\" and \"cell 4 of SEX\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 1 of SEX", "cell 4 of SEX", el("grid")));
      await session.step(107, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(108, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(109, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row height and the missing-value colour paint through the grid's own properties", async () => {
      await session.step(112, "When user sets \"Row Height\" property of grid to \"48\"", () => setProperty(page, "Row Height", el("grid"), "48"));
      await session.step(113, "Then the \"cell 1 of AGE\" area of grid should be at least 40 pixels tall", () => areaAtLeastTall(page, "cell 1 of AGE", el("grid"), 40));
      await session.step(114, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(115, "When user sets \"Row Height\" property of grid to \"28\"", () => setProperty(page, "Row Height", el("grid"), "28"));
      await session.step(116, "And user makes row 298 current", () => makeRowCurrent(page, 298));
      await session.step(117, "And user makes row 300 current", () => makeRowCurrent(page, 300));
      await session.step(118, "And user sets \"Missing Value Color\" property of grid to \"#FFAAAA\"", () => setProperty(page, "Missing Value Color", el("grid"), "#FFAAAA"));
      await session.step(119, "Then the \"cell 298 of HEIGHT\" area of grid should contain the color \"#FFAAAA\"", () => areaColor(page, "cell 298 of HEIGHT", el("grid"), "#FFAAAA"));
      await session.step(120, "When user sets \"Missing Value Color\" property of grid to \"#00FF00\"", () => setProperty(page, "Missing Value Color", el("grid"), "#00FF00"));
      await session.step(121, "Then the \"cell 298 of HEIGHT\" area of grid should contain the color \"#00FF00\"", () => areaColor(page, "cell 298 of HEIGHT", el("grid"), "#00FF00"));
      await session.step(122, "And the \"cell 298 of HEIGHT\" area of grid should not contain the color \"#FFAAAA\"", () => areaNotColor(page, "cell 298 of HEIGHT", el("grid"), "#FFAAAA"));
      await session.step(123, "When user sets \"Missing Value Color\" property of grid to \"#FFFFFF\"", () => setProperty(page, "Missing Value Color", el("grid"), "#FFFFFF"));
      await session.step(124, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(125, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A bigger cell font repaints the grid", async () => {
      await session.step(128, "When user sets \"Default Cell Font\" property of grid to \"20px Roboto\"", () => setProperty(page, "Default Cell Font", el("grid"), "20px Roboto"));
      await session.step(129, "Then grid should have repainted by at least 3000 pixels", () => repaintedBy(page, el("grid"), 3000));
      await session.step(130, "When user sets \"Default Cell Font\" property of grid to \"12px Roboto\"", () => setProperty(page, "Default Cell Font", el("grid"), "12px Roboto"));
      await session.step(131, "Then grid should have repainted", () => repainted(page, el("grid")));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Selected Rows Color paints the selected rows and Escape takes it away", async () => {
      await session.step(135, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(136, "Then no rows should be selected", () => noneSelected(page));
      await session.step(137, "When user sets \"Selected Rows Color\" property of grid to \"#00FF00\"", () => setProperty(page, "Selected Rows Color", el("grid"), "#00FF00"));
      await session.step(138, "Then the \"cell 6 of USUBJID\" area of grid should not contain the color \"#00FF00\"", () => areaNotColor(page, "cell 6 of USUBJID", el("grid"), "#00FF00"));
      await session.step(139, "When user drags a selection box from the \"row header 5\" area to the \"row header 7\" area of grid", () => dragSelectionBetweenAreas(page, "row header 5", "row header 7", el("grid")));
      await session.step(140, "Then rows 5 to 7 should be selected", () => rowsRangeSelected(page, 5, 7));
      await session.step(141, "And the \"cell 6 of USUBJID\" area of grid should contain the color \"#00FF00\"", () => areaColor(page, "cell 6 of USUBJID", el("grid"), "#00FF00"));
      await session.step(142, "And the \"cell 9 of USUBJID\" area of grid should not contain the color \"#00FF00\"", () => areaNotColor(page, "cell 9 of USUBJID", el("grid"), "#00FF00"));
      await session.step(143, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(144, "Then no rows should be selected", () => noneSelected(page));
      await session.step(145, "And the \"cell 6 of USUBJID\" area of grid should not contain the color \"#00FF00\"", () => areaNotColor(page, "cell 6 of USUBJID", el("grid"), "#00FF00"));
      await session.step(146, "When user sets \"Selected Rows Color\" property of grid to \"819780688\"", () => setProperty(page, "Selected Rows Color", el("grid"), "819780688"));
      await session.step(147, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column's colour coding wins over the background its Style sets (GROK-18638)", async () => {
      await session.step(150, "When user clicks on the \"header SEVERITY\" area of grid", () => clickArea(page, "header SEVERITY", el("grid")));
      await session.step(151, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(152, "Then the context panel should show \"SEVERITY\"", () => contextPanelShows(page, "SEVERITY"));
      await session.step(153, "Given Style accordion header in context panel is expanded", () => isExpanded(page, el("Style accordion header in context panel")));
      await session.step(154, "And Content accordion header in context panel is expanded", () => isExpanded(page, el("Content accordion header in context panel")));
      await session.step(155, "When user clicks on editor of \"Back Color\" property in context panel", () => clickOn(page, el("editor of \"Back Color\" property in context panel")));
      await session.step(156, "And user picks the color \"#FFA500\" in the color picker", () => pickColorSwatch(page, "#FFA500"));
      await session.step(157, "Then the \"color of cell 1 of SEVERITY\" reading of grid should be \"#ffa500\"", () => readingReads(page, "color of cell 1 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(158, "And the \"color of cell 2 of SEVERITY\" reading of grid should be \"#ffa500\"", () => readingReads(page, "color of cell 2 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(159, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEVERITY\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEVERITY", el("grid")));
      await session.step(160, "Then the \"color of cell 1 of SEVERITY\" reading of grid should not be \"#ffa500\"", () => readingDoesNotRead(page, "color of cell 1 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(161, "And the \"color of cell 1 of SEVERITY\" and \"color of cell 3 of SEVERITY\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEVERITY", "color of cell 3 of SEVERITY", el("grid")));
      await session.step(162, "When user removes the coloring of \"SEVERITY\" column", () => colorOff(page, "SEVERITY"));
      await session.step(163, "Then the \"color of cell 1 of SEVERITY\" reading of grid should be \"#ffa500\"", () => readingReads(page, "color of cell 1 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(164, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Height typed in the context panel as 0 is clamped to 1 (core#1200) [typed=0, clamped=1]", async () => {
      await session.step(168, "Given user sets \"Row Height\" property of grid to \"28\"", () => setProperty(page, "Row Height", el("grid"), "28"));
      await session.step(169, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(170, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(171, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(172, "Given \"Rows\" category in context panel is expanded", () => isExpanded(page, el("\"Rows\" category in context panel")));
      await session.step(173, "When user enters \"0\" into \"Row Height\" property in context panel", () => enterInto(page, "0", el("\"Row Height\" property in context panel")));
      await session.step(174, "Then \"Row Height\" property of grid should be \"1\"", () => propertyShouldBe(page, "Row Height", el("grid"), "1"));
    });
    await run.scenario("Row Height typed in the context panel as 501 is clamped to 500 (core#1200) [typed=501, clamped=500]", async () => {
      await session.step(168, "Given user sets \"Row Height\" property of grid to \"28\"", () => setProperty(page, "Row Height", el("grid"), "28"));
      await session.step(169, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(170, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(171, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(172, "Given \"Rows\" category in context panel is expanded", () => isExpanded(page, el("\"Rows\" category in context panel")));
      await session.step(173, "When user enters \"501\" into \"Row Height\" property in context panel", () => enterInto(page, "501", el("\"Row Height\" property in context panel")));
      await session.step(174, "Then \"Row Height\" property of grid should be \"500\"", () => propertyShouldBe(page, "Row Height", el("grid"), "500"));
    });
    run.finish();
  });
});
