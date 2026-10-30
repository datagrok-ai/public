/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/sar-matrix-summary.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.sar-matrix]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {matricesBuilt, readingIsTextOf, textHoldsReading} from '../../bindings/sar-matrix.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKey, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldNotBe, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, currentRowValue, noneSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelOpen, openProject, switchView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clearCardSearch, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, clickAreaHolding, infoBalloonText, noErrors, pickFromAreaContextMenu, readingAsRemembered, readingIs, readingNotAsRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix Summary on series cut from the molecules", () => {
  const session = feature(test, "features/analyze/sar-matrix-summary.feature", import.meta.url);
  test("SAR Matrix Summary on series cut from the molecules", {tag: ["@journey", "@realizes:chem.cp.sar-matrix"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(22, "And user opens the \"SarMatrixDemo\" project", () => openProject(page, "SarMatrixDemo"));
    await session.step(23, "And SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
    await run.scenario("A key series opens its best measured compound on the SAR Matrix tab", async () => {
      await session.step(26, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(27, "Then the \"restored\" reading of SAR Matrix Viewer viewer should be \"true\"", () => readingReads(page, "restored", el("SAR Matrix Viewer viewer"), "true"));
      await session.step(28, "And the \"matrices\" reading of SAR Matrix Viewer viewer should be 134", () => readingIs(page, "matrices", el("SAR Matrix Viewer viewer"), 134));
      await session.step(29, "When user clicks on open icon of first key series row", () => clickOn(page, el("open icon of first key series row")));
      await session.step(30, "Then the \"tab\" reading of SAR Matrix Viewer viewer should be \"SAR Matrix\"", () => readingReads(page, "tab", el("SAR Matrix Viewer viewer"), "SAR Matrix"));
      await session.step(31, "And \"CMPD_CHEMBLID\" of the current row should be \"CHEMBL3422978\"", () => currentRowValue(page, "CMPD_CHEMBLID", "CHEMBL3422978"));
      await session.step(32, "And context panel should contain text \"SAR analysis\"", () => shouldContainText(page, el("context panel"), "SAR analysis"));
      await session.step(33, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(34, "Then the \"series\" reading of SAR Matrix Viewer viewer should be the text of name of first key series row", () => readingIsTextOf(page, "series", el("SAR Matrix Viewer viewer"), el("name of first key series row")));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An analog opens its prediction on the SAR Matrix tab", async () => {
      await session.step(38, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(39, "And user clicks on open icon of first analog worth making row", () => clickOn(page, el("open icon of first analog worth making row")));
      await session.step(40, "Then the \"tab\" reading of SAR Matrix Viewer viewer should be \"SAR Matrix\"", () => readingReads(page, "tab", el("SAR Matrix Viewer viewer"), "SAR Matrix"));
      await session.step(41, "And context panel should contain text \"Virtual analog\"", () => shouldContainText(page, el("context panel"), "Virtual analog"));
      await session.step(42, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(43, "Then name of first analog worth making row should contain the \"series\" reading of SAR Matrix Viewer viewer", () => textHoldsReading(page, el("name of first analog worth making row"), "series", el("SAR Matrix Viewer viewer")));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A level chip reads the Summary at one level", async () => {
      await session.step(47, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(48, "Then All level chip should be selected", () => shouldBe(page, el("All level chip"), "selected"));
      await session.step(49, "When user remembers the \"key series\" reading of SAR Matrix Viewer viewer", () => rememberReading(page, "key series", el("SAR Matrix Viewer viewer")));
      await session.step(50, "And user clicks on L2 level chip", () => clickOn(page, el("L2 level chip")));
      await session.step(51, "Then L2 level chip should be selected", () => shouldBe(page, el("L2 level chip"), "selected"));
      await session.step(52, "And All level chip should not be selected", () => shouldNotBe(page, el("All level chip"), "selected"));
      await session.step(53, "And the \"key series\" reading of SAR Matrix Viewer viewer should not be as remembered", () => readingNotAsRemembered(page, "key series", el("SAR Matrix Viewer viewer")));
      await session.step(54, "When user clicks on All level chip", () => clickOn(page, el("All level chip")));
      await session.step(55, "Then All level chip should be selected", () => shouldBe(page, el("All level chip"), "selected"));
      await session.step(56, "And the \"key series\" reading of SAR Matrix Viewer viewer should be as remembered", () => readingAsRemembered(page, "key series", el("SAR Matrix Viewer viewer")));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Make list collects without doubles, removes, hands a copy to the workspace and clears", async () => {
      await session.step(60, "When user clicks on Summary tab in SAR Matrix Viewer viewer", () => clickOn(page, el("Summary tab in SAR Matrix Viewer viewer")));
      await session.step(61, "And user clicks on cart of first analog worth making row", () => clickOn(page, el("cart of first analog worth making row")));
      await session.step(62, "Then an info balloon containing \"Added 1 to the Make list (1 total)\" should have been shown", () => infoBalloonText(page, "Added 1 to the Make list (1 total)"));
      await session.step(63, "When user clicks on cart of first analog worth making row", () => clickOn(page, el("cart of first analog worth making row")));
      await session.step(64, "Then an info balloon containing \"Added 0 to the Make list, 1 already there (1 total)\" should have been shown", () => infoBalloonText(page, "Added 0 to the Make list, 1 already there (1 total)"));
      await session.step(65, "And Add all to Make list button should have text \"Add all 14 to Make list\"", () => shouldHaveText(page, el("Add all to Make list button"), "Add all 14 to Make list"));
      await session.step(66, "When user clicks on Add all to Make list button", () => clickOn(page, el("Add all to Make list button")));
      await session.step(67, "Then an info balloon containing \"Added 13 to the Make list, 1 already there (14 total)\" should have been shown", () => infoBalloonText(page, "Added 13 to the Make list, 1 already there (14 total)"));
      await session.step(68, "When user clicks on Make list tab", () => clickOn(page, el("Make list tab")));
      await session.step(69, "Then Make list title should have text \"14 compounds\"", () => shouldHaveText(page, el("Make list title"), "14 compounds"));
      await session.step(70, "When user clicks on the \"cell 1 of Status\" area of Make list grid viewer", () => clickArea(page, "cell 1 of Status", el("Make list grid viewer")));
      await session.step(71, "And user clicks on Remove button", () => clickOn(page, el("Remove button")));
      await session.step(72, "Then Make list title should have text \"13 compounds\"", () => shouldHaveText(page, el("Make list title"), "13 compounds"));
      await session.step(73, "When user clicks on Add to workspace button", () => clickOn(page, el("Add to workspace button")));
      await session.step(74, "Then the \"SAR make list\" table view should open with 13 rows", () => tableViewOpened(page, "SAR make list", 13));
      await session.step(75, "When user switches to the \"SAR Matrix\" view", () => switchView(page, "SAR Matrix"));
      await session.step(76, "And user clicks on Make list tab", () => clickOn(page, el("Make list tab")));
      await session.step(77, "And user clicks on Clear button", () => clickOn(page, el("Clear button")));
      await session.step(78, "Then SAR Matrix note should contain text \"Nothing collected yet\"", () => shouldContainText(page, el("SAR Matrix note"), "Nothing collected yet"));
      await session.step(79, "And the \"make list\" reading of SAR Matrix Viewer viewer should be 0", () => readingIs(page, "make list", el("SAR Matrix Viewer viewer"), 0));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A measured cell opens its compound and a predicted cell its prediction", async () => {
      await session.step(83, "When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer", () => clickOn(page, el("SAR Matrix tab in SAR Matrix Viewer viewer")));
      await session.step(84, "And user clicks on first series card", () => clickOn(page, el("first series card")));
      await session.step(85, "Then the \"series\" reading of SAR Matrix Viewer viewer should be the text of name of first series card", () => readingIsTextOf(page, "series", el("SAR Matrix Viewer viewer"), el("name of first series card")));
      await session.step(86, "When user clicks on the \"cell 3 of c0\" area of grid viewer", () => clickArea(page, "cell 3 of c0", el("grid viewer")));
      await session.step(87, "Then \"CMPD_CHEMBLID\" of the current row should be \"CHEMBL3422763\"", () => currentRowValue(page, "CMPD_CHEMBLID", "CHEMBL3422763"));
      await session.step(88, "And context panel should contain text \"SAR analysis\"", () => shouldContainText(page, el("context panel"), "SAR analysis"));
      await session.step(89, "And context panel should not contain text \"Virtual analog\"", () => shouldNotContainText(page, el("context panel"), "Virtual analog"));
      await session.step(90, "When user clicks on the \"cell 1 of c0\" area of grid viewer", () => clickArea(page, "cell 1 of c0", el("grid viewer")));
      await session.step(91, "Then context panel should contain text \"Virtual analog\"", () => shouldContainText(page, el("context panel"), "Virtual analog"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl-click and Shift-click add measured compounds to the selection", async () => {
      await session.step(95, "When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer", () => clickOn(page, el("SAR Matrix tab in SAR Matrix Viewer viewer")));
      await session.step(96, "And user clicks on first series card", () => clickOn(page, el("first series card")));
      await session.step(97, "And user clears the row selection", () => clearSelection(page));
      await session.step(98, "And user clicks on the \"cell 3 of c0\" area of grid viewer holding Control", () => clickAreaHolding(page, "cell 3 of c0", el("grid viewer"), "Control"));
      await session.step(99, "Then 1 row should be selected", () => selectedRowCount(page, 1));
      await session.step(100, "When user clicks on the \"cell 4 of c0\" area of grid viewer holding Control", () => clickAreaHolding(page, "cell 4 of c0", el("grid viewer"), "Control"));
      await session.step(101, "Then 2 rows should be selected", () => selectedRowCount(page, 2));
      await session.step(102, "When user clicks on the \"cell 5 of c0\" area of grid viewer holding Shift", () => clickAreaHolding(page, "cell 5 of c0", el("grid viewer"), "Shift"));
      await session.step(103, "Then 3 rows should be selected", () => selectedRowCount(page, 3));
      await session.step(104, "When user clicks on the \"cell 4 of c0\" area of grid viewer holding Control", () => clickAreaHolding(page, "cell 4 of c0", el("grid viewer"), "Control"));
      await session.step(105, "Then 2 rows should be selected", () => selectedRowCount(page, 2));
      await session.step(106, "When user clears the row selection", () => clearSelection(page));
      await session.step(107, "Then no rows should be selected", () => noneSelected(page));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The matrix's cart, a cell's SAR analysis and the cell menu each add a predicted analog once", async () => {
      await session.step(111, "When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer", () => clickOn(page, el("SAR Matrix tab in SAR Matrix Viewer viewer")));
      await session.step(112, "And user clicks on first series card", () => clickOn(page, el("first series card")));
      await session.step(113, "And user clicks on the \"cell 1 of c0\" area of grid viewer", () => clickArea(page, "cell 1 of c0", el("grid viewer")));
      await session.step(114, "And user clicks on matrix cart", () => clickOn(page, el("matrix cart")));
      await session.step(115, "Then an info balloon containing \"Virtual analog added to the Make list (1 total)\" should have been shown", () => infoBalloonText(page, "Virtual analog added to the Make list (1 total)"));
      await session.step(116, "And the \"make list\" reading of SAR Matrix Viewer viewer should be 1", () => readingIs(page, "make list", el("SAR Matrix Viewer viewer"), 1));
      await session.step(117, "When user clicks on the \"cell 6 of c0\" area of grid viewer", () => clickArea(page, "cell 6 of c0", el("grid viewer")));
      await session.step(118, "And user clicks on Add to Make list button in context panel", () => clickOn(page, el("Add to Make list button in context panel")));
      await session.step(119, "Then an info balloon containing \"Virtual analog added to the Make list (2 total)\" should have been shown", () => infoBalloonText(page, "Virtual analog added to the Make list (2 total)"));
      await session.step(120, "And the \"make list\" reading of SAR Matrix Viewer viewer should be 2", () => readingIs(page, "make list", el("SAR Matrix Viewer viewer"), 2));
      await session.step(121, "When user picks \"SAR Matrix > Add this compound to make list\" from the context menu of the \"cell 9 of c0\" area of grid viewer", () => pickFromAreaContextMenu(page, "SAR Matrix > Add this compound to make list", "cell 9 of c0", el("grid viewer")));
      await session.step(122, "Then the \"make list\" reading of SAR Matrix Viewer viewer should be 3", () => readingIs(page, "make list", el("SAR Matrix Viewer viewer"), 3));
      await session.step(123, "When user picks \"SAR Matrix > Add this compound to make list\" from the context menu of the \"cell 1 of c0\" area of grid viewer", () => pickFromAreaContextMenu(page, "SAR Matrix > Add this compound to make list", "cell 1 of c0", el("grid viewer")));
      await session.step(124, "Then an info balloon containing \"Virtual analog is already on the Make list (3 total)\" should have been shown", () => infoBalloonText(page, "Virtual analog is already on the Make list (3 total)"));
      await session.step(125, "And the \"make list\" reading of SAR Matrix Viewer viewer should be 3", () => readingIs(page, "make list", el("SAR Matrix Viewer viewer"), 3));
      await session.step(126, "When user clicks on Make list tab", () => clickOn(page, el("Make list tab")));
      await session.step(127, "And user clicks on Clear button", () => clickOn(page, el("Clear button")));
      await session.step(128, "Then the \"make list\" reading of SAR Matrix Viewer viewer should be 0", () => readingIs(page, "make list", el("SAR Matrix Viewer viewer"), 0));
      await session.step(129, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Caption prints the substituents' molecular weight under the columns", async () => {
      await session.step(132, "When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer", () => clickOn(page, el("SAR Matrix tab in SAR Matrix Viewer viewer")));
      await session.step(133, "And user clicks on first series card", () => clickOn(page, el("first series card")));
      await session.step(134, "Then the \"column captions\" reading of SAR Matrix Viewer viewer should be \"μ 7.2, μ 5.6, μ 6.9, μ 7.5, μ 7.3, μ 7.9, μ 8.2, μ 5.0, μ 5.0, μ 5.2, μ 5.1, μ 7.7, μ 7.7, μ 6.1, μ 6.8, μ 5.7, μ 4.9, μ 5.2, μ 4.6, μ 4.6\"", () => readingReads(page, "column captions", el("SAR Matrix Viewer viewer"), "μ 7.2, μ 5.6, μ 6.9, μ 7.5, μ 7.3, μ 7.9, μ 8.2, μ 5.0, μ 5.0, μ 5.2, μ 5.1, μ 7.7, μ 7.7, μ 6.1, μ 6.8, μ 5.7, μ 4.9, μ 5.2, μ 4.6, μ 4.6"));
      await session.step(135, "When user selects \"Molecular weight\" in \"Caption\" input", () => selectIn(page, "Molecular weight", el("\"Caption\" input")));
      await session.step(136, "Then the \"column captions\" reading of SAR Matrix Viewer viewer should be \"MW 32, MW 17, MW 16, MW 20, MW 27, MW 36, MW 70, MW 31, MW 32, MW 18, MW 74, MW 46, MW 30, MW 60, MW 80, MW 52, MW 45, MW 46, MW 62, MW 46\"", () => readingReads(page, "column captions", el("SAR Matrix Viewer viewer"), "MW 32, MW 17, MW 16, MW 20, MW 27, MW 36, MW 70, MW 31, MW 32, MW 18, MW 74, MW 46, MW 30, MW 60, MW 80, MW 52, MW 45, MW 46, MW 62, MW 46"));
      await session.step(137, "When user selects \"Potency\" in \"Caption\" input", () => selectIn(page, "Potency", el("\"Caption\" input")));
      await session.step(138, "Then the \"column captions\" reading of SAR Matrix Viewer viewer should be \"μ 7.2, μ 5.6, μ 6.9, μ 7.5, μ 7.3, μ 7.9, μ 8.2, μ 5.0, μ 5.0, μ 5.2, μ 5.1, μ 7.7, μ 7.7, μ 6.1, μ 6.8, μ 5.7, μ 4.9, μ 5.2, μ 4.6, μ 4.6\"", () => readingReads(page, "column captions", el("SAR Matrix Viewer viewer"), "μ 7.2, μ 5.6, μ 6.9, μ 7.5, μ 7.3, μ 7.9, μ 8.2, μ 5.0, μ 5.0, μ 5.2, μ 5.1, μ 7.7, μ 7.7, μ 6.1, μ 6.8, μ 5.7, μ 4.9, μ 5.2, μ 4.6, μ 4.6"));
      await session.step(139, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Rank by reorders the series list without a rebuild, and the series filter counts what it keeps", async () => {
      await session.step(142, "When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer", () => clickOn(page, el("SAR Matrix tab in SAR Matrix Viewer viewer")));
      await session.step(143, "And user remembers the \"series list\" reading of SAR Matrix Viewer viewer", () => rememberReading(page, "series list", el("SAR Matrix Viewer viewer")));
      await session.step(144, "And user selects \"SAR discontinuity\" in \"Rank by\" input", () => selectIn(page, "SAR discontinuity", el("\"Rank by\" input")));
      await session.step(145, "Then the \"series list\" reading of SAR Matrix Viewer viewer should not be as remembered", () => readingNotAsRemembered(page, "series list", el("SAR Matrix Viewer viewer")));
      await session.step(146, "And the \"restored\" reading of SAR Matrix Viewer viewer should be \"true\"", () => readingReads(page, "restored", el("SAR Matrix Viewer viewer"), "true"));
      await session.step(147, "When user selects \"Potent compounds\" in \"Rank by\" input", () => selectIn(page, "Potent compounds", el("\"Rank by\" input")));
      await session.step(148, "Then the \"series list\" reading of SAR Matrix Viewer viewer should be as remembered", () => readingAsRemembered(page, "series list", el("SAR Matrix Viewer viewer")));
      await session.step(149, "When user clicks on series filter icon", () => clickOn(page, el("series filter icon")));
      await session.step(150, "And user hovers over \"Series\" filter card", () => hoverOver(page, el("\"Series\" filter card")));
      await session.step(151, "And user clicks on search icon of \"Series\" filter card", () => clickOn(page, el("search icon of \"Series\" filter card")));
      await session.step(152, "And user types \"Series 69\" into the search of the \"Series\" filter card", () => typeIntoCardSearch(page, "Series 69", "Series"));
      await session.step(153, "Then series match count should have text \"1 of 134 match\"", () => shouldHaveText(page, el("series match count"), "1 of 134 match"));
      await session.step(154, "And the \"series list\" reading of SAR Matrix Viewer viewer should be \"Series 69\"", () => readingReads(page, "series list", el("SAR Matrix Viewer viewer"), "Series 69"));
      await session.step(155, "When user clears the search of the \"Series\" filter card", () => clearCardSearch(page, "Series"));
      await session.step(156, "Then the \"series list\" reading of SAR Matrix Viewer viewer should be as remembered", () => readingAsRemembered(page, "series list", el("SAR Matrix Viewer viewer")));
      await session.step(157, "And series match count should have text \"\"", () => shouldHaveText(page, el("series match count"), ""));
      await session.step(158, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(159, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A transfer filter moves the pane to the transfer it keeps, and says when it keeps none", async () => {
      await session.step(162, "When user clicks on SAR Transfer tab in SAR Matrix Viewer viewer", () => clickOn(page, el("SAR Transfer tab in SAR Matrix Viewer viewer")));
      await session.step(163, "Then transfer title should have text \"Series 69 · Core 3 → Series 33 · Core 1\"", () => shouldHaveText(page, el("transfer title"), "Series 69 · Core 3 → Series 33 · Core 1"));
      await session.step(164, "When user clicks on transfer filter icon", () => clickOn(page, el("transfer filter icon")));
      await session.step(165, "And user hovers over \"Target\" filter card", () => hoverOver(page, el("\"Target\" filter card")));
      await session.step(166, "And user clicks on search icon of \"Target\" filter card", () => clickOn(page, el("search icon of \"Target\" filter card")));
      await session.step(167, "And user types \"Series 98\" into the search of the \"Target\" filter card", () => typeIntoCardSearch(page, "Series 98", "Target"));
      await session.step(168, "Then transfer title should have text \"Series 69 · Core 3 → Series 98 · Core 1\"", () => shouldHaveText(page, el("transfer title"), "Series 69 · Core 3 → Series 98 · Core 1"));
      await session.step(169, "When user types \"Series 5\" into the search of the \"Target\" filter card", () => typeIntoCardSearch(page, "Series 5", "Target"));
      await session.step(170, "Then SAR Matrix note should have text \"No transfer matches the filter — clear it to see the transfers again.\"", () => shouldHaveText(page, el("SAR Matrix note"), "No transfer matches the filter — clear it to see the transfers again."));
      await session.step(171, "When user clears the search of the \"Target\" filter card", () => clearCardSearch(page, "Target"));
      await session.step(172, "And user presses Escape", () => pressKey(page, "Escape"));
      await session.step(173, "And user selects \"Series 33 · Core 1 · R1\" in \"Transfer to\" input", () => selectIn(page, "Series 33 · Core 1 · R1", el("\"Transfer to\" input")));
      await session.step(174, "Then transfer title should have text \"Series 69 · Core 3 → Series 33 · Core 1\"", () => shouldHaveText(page, el("transfer title"), "Series 69 · Core 3 → Series 33 · Core 1"));
      await session.step(175, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
