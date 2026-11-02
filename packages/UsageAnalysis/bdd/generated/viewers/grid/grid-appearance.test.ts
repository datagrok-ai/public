/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-appearance.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#1200, core#1198]
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
  test("Grid cell appearance", {tag: ["@journey", "@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:175", "@realizes:core#1200", "@known-failure:188", "@realizes:core#1198", "@known-failure:198"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 14, page);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(31, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(32, "And \"Color Coding\" property of grid should be \"Auto\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "Auto"));
    await session.step(33, "And \"HEIGHT\" column should have missing values", () => columnIncomplete(page, "HEIGHT"));
    await run.scenario("Linear coding from the header menu paints the column by value", async () => {
      await session.step(36, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(37, "Then \"AGE\" column should have tag \".color-coding-type\" equal to \"Linear\"", () => columnTag(page, "AGE", ".color-coding-type", "Linear"));
      await session.step(38, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(39, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(40, "Then the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(41, "And the \"cell 2 of AGE\" and \"cell 2 of USUBJID\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 2 of USUBJID", el("grid")));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The grid-wide Color Coding overrides the column's", async () => {
      await session.step(45, "When user picks \"Grid Color Coding > None\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > None", "cell 8 of AGE", el("grid")));
      await session.step(46, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(47, "Then \"Color Coding\" property of grid should be \"None\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "None"));
      await session.step(48, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in the same colors", () => areasSame(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(49, "When user picks \"Grid Color Coding > All\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > All", "cell 8 of AGE", el("grid")));
      await session.step(50, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(51, "Then \"Color Coding\" property of grid should be \"All\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "All"));
      await session.step(52, "And the \"cell 2 of HEIGHT\" and \"cell 4 of HEIGHT\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of HEIGHT", "cell 4 of HEIGHT", el("grid")));
      await session.step(53, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(54, "When user picks \"Grid Color Coding > Auto\" from the context menu of the \"cell 8 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > Auto", "cell 8 of AGE", el("grid")));
      await session.step(55, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(56, "Then \"Color Coding\" property of grid should be \"Auto\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "Auto"));
      await session.step(57, "And the \"cell 2 of AGE\" and \"cell 3 of AGE\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of AGE", "cell 3 of AGE", el("grid")));
      await session.step(58, "And the \"cell 2 of HEIGHT\" and \"cell 4 of HEIGHT\" areas of grid should be painted in the same colors", () => areasSame(page, "cell 2 of HEIGHT", "cell 4 of HEIGHT", el("grid")));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Narrowing a coloured column keeps its colour", async () => {
      await session.step(62, "When user picks \"Grid Color Coding > All\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > All", "cell 1 of AGE", el("grid")));
      await session.step(63, "And user remembers the \"color of cell 1 of AGE\" reading of grid", () => rememberReading(page, "color of cell 1 of AGE", el("grid")));
      await session.step(64, "And user drags the \"column resizer AGE\" area of grid by 30 pixels to the left", () => dragAreaBy(page, "column resizer AGE", el("grid"), 30, "left"));
      await session.step(65, "Then the \"column width of AGE\" reading of grid should be lower than before", () => readingLower(page, "column width of AGE", el("grid")));
      await session.step(66, "And the \"color of cell 1 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 1 of AGE", el("grid")));
      await session.step(67, "When user picks \"Column Sizing > Optimal\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Column Sizing > Optimal", "cell 1 of AGE", el("grid")));
      await session.step(68, "And user picks \"Grid Color Coding > Auto\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > Auto", "cell 1 of AGE", el("grid")));
      await session.step(69, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A custom format shows in the cell", async () => {
      await session.step(72, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(73, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(74, "When user enters \"0.00\" into Custom input in Format AGE dialog", () => enterInto(page, "0.00", el("Custom input in Format AGE dialog")));
      await session.step(75, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(76, "Then Format AGE dialog should be hidden", () => shouldBe(page, el("Format AGE dialog"), "hidden"));
      await session.step(77, "And \"AGE\" column should have tag \"format\" equal to \"0.00\"", () => columnTag(page, "AGE", "format", "0.00"));
      await session.step(78, "And the \"text of cell 1 of AGE\" reading of grid should be \"26.00\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26.00"));
      await session.step(79, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
      await session.step(80, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(81, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
      await session.step(82, "Given Details accordion header in context panel is expanded", () => isExpanded(page, el("Details accordion header in context panel")));
      await session.step(83, "Then \"format\" table row in context panel should contain text \"0.00\"", () => shouldContainText(page, el("\"format\" table row in context panel"), "0.00"));
      await session.step(84, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(85, "And user clears Custom input in Format AGE dialog", () => clearField(page, el("Custom input in Format AGE dialog")));
      await session.step(86, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(87, "Then the \"text of cell 1 of AGE\" reading of grid should be \"26\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26"));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Conditional coding paints the cells its ranges name", async () => {
      await session.step(91, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header HEIGHT", el("grid")));
      await session.step(92, "Then \"HEIGHT\" column should have tag \".color-coding-type\" equal to \"Conditional\"", () => columnTag(page, "HEIGHT", ".color-coding-type", "Conditional"));
      await session.step(93, "When user colors \"HEIGHT\" column conditionally:", () => colorConditional(page, "HEIGHT", [["<160","#0000FF"],[">180","#FF0000"]]), [["<160","#0000FF"],[">180","#FF0000"]]);
      await session.step(96, "Then the \"cell 2 of HEIGHT\" area of grid should contain the color \"#0000FF\"", () => areaColor(page, "cell 2 of HEIGHT", el("grid"), "#0000FF"));
      await session.step(97, "And the \"cell 4 of HEIGHT\" area of grid should contain the color \"#FF0000\"", () => areaColor(page, "cell 4 of HEIGHT", el("grid"), "#FF0000"));
      await session.step(98, "And the \"cell 1 of HEIGHT\" area of grid should not contain the color \"#0000FF\"", () => areaNotColor(page, "cell 1 of HEIGHT", el("grid"), "#0000FF"));
      await session.step(99, "And the \"cell 1 of HEIGHT\" area of grid should not contain the color \"#FF0000\"", () => areaNotColor(page, "cell 1 of HEIGHT", el("grid"), "#FF0000"));
      await session.step(100, "When user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
      await session.step(101, "Then \"HEIGHT\" column should have no color coding", () => noColorCoding(page, "HEIGHT"));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Categorical coding tells the categories apart", async () => {
      await session.step(105, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(106, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(107, "And the \"cell 1 of SEX\" and \"cell 4 of SEX\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 1 of SEX", "cell 4 of SEX", el("grid")));
      await session.step(108, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(109, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row height and the missing-value colour paint through the grid's own properties", async () => {
      await session.step(113, "When user sets \"Row Height\" property of grid to \"48\"", () => setProperty(page, "Row Height", el("grid"), "48"));
      await session.step(114, "Then the \"cell 1 of AGE\" area of grid should be at least 40 pixels tall", () => areaAtLeastTall(page, "cell 1 of AGE", el("grid"), 40));
      await session.step(115, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(116, "When user sets \"Row Height\" property of grid to \"28\"", () => setProperty(page, "Row Height", el("grid"), "28"));
      await session.step(117, "And user makes row 298 current", () => makeRowCurrent(page, 298));
      await session.step(118, "And user makes row 300 current", () => makeRowCurrent(page, 300));
      await session.step(119, "And user sets \"Missing Value Color\" property of grid to \"#FFAAAA\"", () => setProperty(page, "Missing Value Color", el("grid"), "#FFAAAA"));
      await session.step(120, "Then the \"cell 298 of HEIGHT\" area of grid should contain the color \"#FFAAAA\"", () => areaColor(page, "cell 298 of HEIGHT", el("grid"), "#FFAAAA"));
      await session.step(121, "When user sets \"Missing Value Color\" property of grid to \"#00FF00\"", () => setProperty(page, "Missing Value Color", el("grid"), "#00FF00"));
      await session.step(122, "Then the \"cell 298 of HEIGHT\" area of grid should contain the color \"#00FF00\"", () => areaColor(page, "cell 298 of HEIGHT", el("grid"), "#00FF00"));
      await session.step(123, "And the \"cell 298 of HEIGHT\" area of grid should not contain the color \"#FFAAAA\"", () => areaNotColor(page, "cell 298 of HEIGHT", el("grid"), "#FFAAAA"));
      await session.step(124, "When user sets \"Missing Value Color\" property of grid to \"#FFFFFF\"", () => setProperty(page, "Missing Value Color", el("grid"), "#FFFFFF"));
      await session.step(125, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(126, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A bigger cell font repaints the grid", async () => {
      await session.step(129, "When user sets \"Default Cell Font\" property of grid to \"20px Roboto\"", () => setProperty(page, "Default Cell Font", el("grid"), "20px Roboto"));
      await session.step(130, "Then grid should have repainted by at least 3000 pixels", () => repaintedBy(page, el("grid"), 3000));
      await session.step(131, "When user sets \"Default Cell Font\" property of grid to \"12px Roboto\"", () => setProperty(page, "Default Cell Font", el("grid"), "12px Roboto"));
      await session.step(132, "Then grid should have repainted", () => repainted(page, el("grid")));
      await session.step(133, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Selected Rows Color paints the selected rows and Escape takes it away", async () => {
      await session.step(136, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(137, "Then no rows should be selected", () => noneSelected(page));
      await session.step(138, "When user sets \"Selected Rows Color\" property of grid to \"#00FF00\"", () => setProperty(page, "Selected Rows Color", el("grid"), "#00FF00"));
      await session.step(139, "Then the \"cell 6 of USUBJID\" area of grid should not contain the color \"#00FF00\"", () => areaNotColor(page, "cell 6 of USUBJID", el("grid"), "#00FF00"));
      await session.step(140, "When user drags a selection box from the \"row header 5\" area to the \"row header 7\" area of grid", () => dragSelectionBetweenAreas(page, "row header 5", "row header 7", el("grid")));
      await session.step(141, "Then rows 5 to 7 should be selected", () => rowsRangeSelected(page, 5, 7));
      await session.step(142, "And the \"cell 6 of USUBJID\" area of grid should contain the color \"#00FF00\"", () => areaColor(page, "cell 6 of USUBJID", el("grid"), "#00FF00"));
      await session.step(143, "And the \"cell 9 of USUBJID\" area of grid should not contain the color \"#00FF00\"", () => areaNotColor(page, "cell 9 of USUBJID", el("grid"), "#00FF00"));
      await session.step(144, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(145, "Then no rows should be selected", () => noneSelected(page));
      await session.step(146, "And the \"cell 6 of USUBJID\" area of grid should not contain the color \"#00FF00\"", () => areaNotColor(page, "cell 6 of USUBJID", el("grid"), "#00FF00"));
      await session.step(147, "When user sets \"Selected Rows Color\" property of grid to \"819780688\"", () => setProperty(page, "Selected Rows Color", el("grid"), "819780688"));
      await session.step(148, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column's colour coding wins over the background its Style sets (GROK-18638)", async () => {
      await session.step(151, "When user clicks on the \"header SEVERITY\" area of grid", () => clickArea(page, "header SEVERITY", el("grid")));
      await session.step(152, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(153, "Then the context panel should show \"SEVERITY\"", () => contextPanelShows(page, "SEVERITY"));
      await session.step(154, "Given Style accordion header in context panel is expanded", () => isExpanded(page, el("Style accordion header in context panel")));
      await session.step(155, "And Content accordion header in context panel is expanded", () => isExpanded(page, el("Content accordion header in context panel")));
      await session.step(156, "When user clicks on editor of \"Back Color\" property in context panel", () => clickOn(page, el("editor of \"Back Color\" property in context panel")));
      await session.step(157, "And user picks the color \"#FFA500\" in the color picker", () => pickColorSwatch(page, "#FFA500"));
      await session.step(158, "Then the \"color of cell 1 of SEVERITY\" reading of grid should be \"#ffa500\"", () => readingReads(page, "color of cell 1 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(159, "And the \"color of cell 2 of SEVERITY\" reading of grid should be \"#ffa500\"", () => readingReads(page, "color of cell 2 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(160, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEVERITY\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEVERITY", el("grid")));
      await session.step(161, "Then the \"color of cell 1 of SEVERITY\" reading of grid should not be \"#ffa500\"", () => readingDoesNotRead(page, "color of cell 1 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(162, "And the \"color of cell 1 of SEVERITY\" and \"color of cell 3 of SEVERITY\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEVERITY", "color of cell 3 of SEVERITY", el("grid")));
      await session.step(163, "When user removes the coloring of \"SEVERITY\" column", () => colorOff(page, "SEVERITY"));
      await session.step(164, "Then the \"color of cell 1 of SEVERITY\" reading of grid should be \"#ffa500\"", () => readingReads(page, "color of cell 1 of SEVERITY", el("grid"), "#ffa500"));
      await session.step(165, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Height typed in the context panel as 0 is not taken [typed=0]", async () => {
      await session.step(169, "Given user sets \"Row Height\" property of grid to \"28\"", () => setProperty(page, "Row Height", el("grid"), "28"));
      await session.step(170, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(171, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(172, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(173, "Given \"Rows\" category in context panel is expanded", () => isExpanded(page, el("\"Rows\" category in context panel")));
      await session.step(174, "When user enters \"0\" into \"Row Height\" property in context panel", () => enterInto(page, "0", el("\"Row Height\" property in context panel")));
      await session.step(175, "Then \"Row Height\" property of grid should be \"28\"", () => propertyShouldBe(page, "Row Height", el("grid"), "28"));
    }, {knownFailure: true, line: 175});
    await run.scenario("Row Height typed in the context panel as 501 is not taken [typed=501]", async () => {
      await session.step(169, "Given user sets \"Row Height\" property of grid to \"28\"", () => setProperty(page, "Row Height", el("grid"), "28"));
      await session.step(170, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(171, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(172, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(173, "Given \"Rows\" category in context panel is expanded", () => isExpanded(page, el("\"Rows\" category in context panel")));
      await session.step(174, "When user enters \"501\" into \"Row Height\" property in context panel", () => enterInto(page, "501", el("\"Row Height\" property in context panel")));
      await session.step(175, "Then \"Row Height\" property of grid should be \"28\"", () => propertyShouldBe(page, "Row Height", el("grid"), "28"));
    }, {knownFailure: true, line: 175});
    await run.scenario("A malformed number format marks the format field invalid", async () => {
      await session.step(185, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(186, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(187, "When user enters \"0.0.0\" into Custom input in Format AGE dialog", () => enterInto(page, "0.0.0", el("Custom input in Format AGE dialog")));
      await session.step(188, "Then Custom input in Format AGE dialog should be invalid", () => shouldBe(page, el("Custom input in Format AGE dialog"), "invalid"));
    }, {knownFailure: true, line: 188});
    await run.scenario("A malformed number format is not taken and the old format stays", async () => {
      await session.step(192, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(193, "Then Format AGE dialog should be absent", () => shouldBe(page, el("Format AGE dialog"), "absent"));
      await session.step(194, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(195, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(196, "When user enters \"0.0.0\" into Custom input in Format AGE dialog", () => enterInto(page, "0.0.0", el("Custom input in Format AGE dialog")));
      await session.step(197, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(198, "Then \"AGE\" column should have tag \"format\" equal to \"\"", () => columnTag(page, "AGE", "format", ""));
    }, {knownFailure: true, line: 198});
    run.finish();
  });
});
