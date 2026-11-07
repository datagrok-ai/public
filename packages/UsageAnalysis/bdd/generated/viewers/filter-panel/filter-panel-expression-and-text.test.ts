/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-expression-and-text.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
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
import {clearField, clickOn, dragSliderTo, hoverOver, pasteInto, pressKey, selectIn, shouldBe, shouldHaveText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactly, filterIsExactlyCategory, filterIsExactlyContains, filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickPanelMenu, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, loadLayout, noErrors, pickFromAreaContextMenu, readingAtLeast, readingHigherThanRemembered, readingReads, rememberReading, saveLayout} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Expression and text filter cards", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-expression-and-text.feature", import.meta.url);
  test("Expression and text filter cards", {tag: ["@journey", "@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(36, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(37, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
    await session.step(38, "Then \"Expression\" filter card should be visible", () => shouldBe(page, el("\"Expression\" filter card"), "visible"));
    await session.step(39, "And the \"type of Expression\" reading of filter panel should be \"expression\"", () => readingReads(page, "type of Expression", el("filter panel"), "expression"));
    await session.step(40, "And all rows should pass the filter", () => filterPassesAll(page));
    await run.scenario("Two rules from the form combine with OR, and the switch makes them AND", async () => {
      await session.step(43, "When user selects \"AGE\" in Column input in \"Expression\" filter card", () => selectIn(page, "AGE", el("Column input in \"Expression\" filter card")));
      await session.step(44, "And user selects \">\" in Operation input in \"Expression\" filter card", () => selectIn(page, ">", el("Operation input in \"Expression\" filter card")));
      await session.step(45, "And user types \"50\" into Value input in \"Expression\" filter card", () => typeInto(page, "50", el("Value input in \"Expression\" filter card")));
      await session.step(46, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(47, "Then 367 rows should pass the filter", () => filterPasses(page, 367));
      await session.step(48, "And the filter should pass exactly the rows where \"AGE\" is between 51 and 89", () => filterIsExactly(page, "AGE", 51, 89));
      await session.step(49, "And the \"categories of Expression\" reading of filter panel should be \"${AGE} > 50\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${AGE} > 50"));
      await session.step(50, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(51, "When user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
      await session.step(52, "And user selects \"<\" in Operation input in \"Expression\" filter card", () => selectIn(page, "<", el("Operation input in \"Expression\" filter card")));
      await session.step(53, "And user types \"160\" into Value input in \"Expression\" filter card", () => typeInto(page, "160", el("Value input in \"Expression\" filter card")));
      await session.step(54, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(55, "Then the \"categories of Expression\" reading of filter panel should be \"${AGE} > 50, ${HEIGHT} < 160\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${AGE} > 50, ${HEIGHT} < 160"));
      await session.step(56, "And mode of \"Expression\" filter card should have text \"OR\"", () => shouldHaveText(page, el("mode of \"Expression\" filter card"), "OR"));
      await session.step(57, "And 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(58, "When user clicks on mode of \"Expression\" filter card", () => clickOn(page, el("mode of \"Expression\" filter card")));
      await session.step(59, "Then mode of \"Expression\" filter card should have text \"AND\"", () => shouldHaveText(page, el("mode of \"Expression\" filter card"), "AND"));
      await session.step(60, "And 88 rows should pass the filter", () => filterPasses(page, 88));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A rule's own checkbox suspends it and keeps it in the list", async () => {
      await session.step(64, "When user clicks on the \"checkbox ${HEIGHT} < 160 of Expression\" area of filter panel", () => clickArea(page, "checkbox ${HEIGHT} < 160 of Expression", el("filter panel")));
      await session.step(65, "Then 367 rows should pass the filter", () => filterPasses(page, 367));
      await session.step(66, "And the \"categories of Expression\" reading of filter panel should be \"${AGE} > 50, ${HEIGHT} < 160\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${AGE} > 50, ${HEIGHT} < 160"));
      await session.step(67, "And the \"selected categories of Expression\" reading of filter panel should be \"${AGE} > 50\"", () => readingReads(page, "selected categories of Expression", el("filter panel"), "${AGE} > 50"));
      await session.step(68, "When user clicks on the \"checkbox ${HEIGHT} < 160 of Expression\" area of filter panel", () => clickArea(page, "checkbox ${HEIGHT} < 160 of Expression", el("filter panel")));
      await session.step(69, "Then 88 rows should pass the filter", () => filterPasses(page, 88));
      await session.step(70, "And the \"selected categories of Expression\" reading of filter panel should be \"${AGE} > 50, ${HEIGHT} < 160\"", () => readingReads(page, "selected categories of Expression", el("filter panel"), "${AGE} > 50, ${HEIGHT} < 160"));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Remove Query in a rule's menu takes the rule out", async () => {
      await session.step(74, "When user picks \"Remove Query\" from the context menu of the \"category ${AGE} > 50 of Expression\" area of filter panel", () => pickFromAreaContextMenu(page, "Remove Query", "category ${AGE} > 50 of Expression", el("filter panel")));
      await session.step(75, "Then 167 rows should pass the filter", () => filterPasses(page, 167));
      await session.step(76, "And the \"categories of Expression\" reading of filter panel should be \"${HEIGHT} < 160\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${HEIGHT} < 160"));
      await session.step(77, "And the filter should pass exactly the rows where \"HEIGHT\" is between 1 and 159.999", () => filterIsExactly(page, "HEIGHT", 1, 159.999));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Free-text mode takes a rule in the card's own syntax and the form comes back", async () => {
      await session.step(81, "When user hovers over \"Expression\" filter card", () => hoverOver(page, el("\"Expression\" filter card")));
      await session.step(82, "And user clicks on \"Switch to free-text mode\" icon in \"Expression\" filter card", () => clickOn(page, el("\"Switch to free-text mode\" icon in \"Expression\" filter card")));
      await session.step(83, "Then Column input in \"Expression\" filter card should be hidden", () => shouldBe(page, el("Column input in \"Expression\" filter card"), "hidden"));
      await session.step(84, "When user types \"${HEIGHT} < 150\" into the search of the \"Expression\" filter card", () => typeIntoCardSearch(page, "${HEIGHT} < 150", "Expression"));
      await session.step(85, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(86, "Then the \"categories of Expression\" reading of filter panel should be \"${HEIGHT} < 160, ${HEIGHT} < 150\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${HEIGHT} < 160, ${HEIGHT} < 150"));
      await session.step(87, "And 16 rows should pass the filter", () => filterPasses(page, 16));
      await session.step(88, "When user hovers over \"Expression\" filter card", () => hoverOver(page, el("\"Expression\" filter card")));
      await session.step(89, "And user clicks on \"Switch to free-text mode\" icon in \"Expression\" filter card", () => clickOn(page, el("\"Switch to free-text mode\" icon in \"Expression\" filter card")));
      await session.step(90, "Then Column input in \"Expression\" filter card should be visible", () => shouldBe(page, el("Column input in \"Expression\" filter card"), "visible"));
      await session.step(91, "And 16 rows should pass the filter", () => filterPasses(page, 16));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the rules back with the suspended one still unchecked", async () => {
      await session.step(95, "When user clicks on the \"checkbox ${HEIGHT} < 150 of Expression\" area of filter panel", () => clickArea(page, "checkbox ${HEIGHT} < 150 of Expression", el("filter panel")));
      await session.step(96, "Then 167 rows should pass the filter", () => filterPasses(page, 167));
      await session.step(97, "And the \"selected categories of Expression\" reading of filter panel should be \"${HEIGHT} < 160\"", () => readingReads(page, "selected categories of Expression", el("filter panel"), "${HEIGHT} < 160"));
      await session.step(98, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(99, "And user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(100, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(101, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(102, "When user loads the saved layout", () => loadLayout(page));
      await session.step(103, "Then \"Expression\" filter card should be visible", () => shouldBe(page, el("\"Expression\" filter card"), "visible"));
      await session.step(104, "And the \"categories of Expression\" reading of filter panel should be \"${HEIGHT} < 160, ${HEIGHT} < 150\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${HEIGHT} < 160, ${HEIGHT} < 150"));
      await session.step(105, "And the \"selected categories of Expression\" reading of filter panel should be \"${HEIGHT} < 160\"", () => readingReads(page, "selected categories of Expression", el("filter panel"), "${HEIGHT} < 160"));
      await session.step(106, "And 167 rows should pass the filter", () => filterPasses(page, 167));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The regex operation turns a pasted list into an alternation and drops a trailing comma", async () => {
      await session.step(110, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(111, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(112, "And user selects \"USUBJID\" in Column input in \"Expression\" filter card", () => selectIn(page, "USUBJID", el("Column input in \"Expression\" filter card")));
      await session.step(113, "And user selects \"regex\" in Operation input in \"Expression\" filter card", () => selectIn(page, "regex", el("Operation input in \"Expression\" filter card")));
      await session.step(114, "And user pastes \"5,15,25\" into Value input in \"Expression\" filter card", () => pasteInto(page, "5,15,25", el("Value input in \"Expression\" filter card")));
      await session.step(115, "Then Value input in \"Expression\" filter card should have the value \"5|15|25\"", () => shouldHaveValue(page, el("Value input in \"Expression\" filter card"), "5|15|25"));
      await session.step(116, "When user clears Value input in \"Expression\" filter card", () => clearField(page, el("Value input in \"Expression\" filter card")));
      await session.step(117, "Then Value input in \"Expression\" filter card should have the value \"\"", () => shouldHaveValue(page, el("Value input in \"Expression\" filter card"), ""));
      await session.step(118, "When user pastes \"5,15,25,\" into Value input in \"Expression\" filter card", () => pasteInto(page, "5,15,25,", el("Value input in \"Expression\" filter card")));
      await session.step(119, "Then Value input in \"Expression\" filter card should have the value \"5|15|25\"", () => shouldHaveValue(page, el("Value input in \"Expression\" filter card"), "5|15|25"));
      await session.step(120, "When user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(121, "Then the \"categories of Expression\" reading of filter panel should be \"${USUBJID} regex 5|15|25\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${USUBJID} regex 5|15|25"));
      await session.step(122, "And 357 rows should pass the filter", () => filterPasses(page, 357));
      await session.step(123, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The string and date operations keep the rows they name", async () => {
      await session.step(126, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(127, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(128, "And user selects \"SEX\" in Column input in \"Expression\" filter card", () => selectIn(page, "SEX", el("Column input in \"Expression\" filter card")));
      await session.step(129, "And user selects \"equals\" in Operation input in \"Expression\" filter card", () => selectIn(page, "equals", el("Operation input in \"Expression\" filter card")));
      await session.step(130, "And user types \"F\" into Value input in \"Expression\" filter card", () => typeInto(page, "F", el("Value input in \"Expression\" filter card")));
      await session.step(131, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(132, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(133, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
      await session.step(134, "And the \"categories of Expression\" reading of filter panel should be \"${SEX} equals F\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${SEX} equals F"));
      await session.step(135, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(136, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(137, "And user selects \"RACE\" in Column input in \"Expression\" filter card", () => selectIn(page, "RACE", el("Column input in \"Expression\" filter card")));
      await session.step(138, "And user selects \"contains\" in Operation input in \"Expression\" filter card", () => selectIn(page, "contains", el("Operation input in \"Expression\" filter card")));
      await session.step(139, "And user types \"an\" into Value input in \"Expression\" filter card", () => typeInto(page, "an", el("Value input in \"Expression\" filter card")));
      await session.step(140, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(141, "Then 911 rows should pass the filter", () => filterPasses(page, 911));
      await session.step(142, "And the filter should pass exactly the rows where \"RACE\" contains \"an\"", () => filterIsExactlyContains(page, "RACE", "an"));
      await session.step(143, "And the \"categories of Expression\" reading of filter panel should be \"${RACE} contains an\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${RACE} contains an"));
      await session.step(144, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(145, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(146, "And user selects \"STARTED\" in Column input in \"Expression\" filter card", () => selectIn(page, "STARTED", el("Column input in \"Expression\" filter card")));
      await session.step(147, "And user selects \"after\" in Operation input in \"Expression\" filter card", () => selectIn(page, "after", el("Operation input in \"Expression\" filter card")));
      await session.step(148, "And user types \"01/01/1991\" into Value input in \"Expression\" filter card", () => typeInto(page, "01/01/1991", el("Value input in \"Expression\" filter card")));
      await session.step(149, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(150, "Then 483 rows should pass the filter", () => filterPasses(page, 483));
      await session.step(151, "And the \"categories of Expression\" reading of filter panel should be \"${STARTED} after 01/01/1991\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${STARTED} after 01/01/1991"));
      await session.step(152, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(153, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(154, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The text card keeps the rows holding any of its terms, or all of them", async () => {
      await session.step(157, "When user opens beer dataset", () => openDataset(page, ds("beer")));
      await session.step(158, "And user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(159, "Then \"Aroma\" filter card should be visible", () => shouldBe(page, el("\"Aroma\" filter card"), "visible"));
      await session.step(160, "And the \"type of Aroma\" reading of filter panel should be \"text\"", () => readingReads(page, "type of Aroma", el("filter panel"), "text"));
      await session.step(161, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(162, "When user types \"malt\" into the search of the \"Aroma\" filter card", () => typeIntoCardSearch(page, "malt", "Aroma"));
      await session.step(163, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(164, "Then 92 rows should pass the filter", () => filterPasses(page, 92));
      await session.step(165, "And the \"categories of Aroma\" reading of filter panel should be \"malt\"", () => readingReads(page, "categories of Aroma", el("filter panel"), "malt"));
      await session.step(166, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(167, "When user types \"hop\" into the search of the \"Aroma\" filter card", () => typeIntoCardSearch(page, "hop", "Aroma"));
      await session.step(168, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(169, "Then the \"categories of Aroma\" reading of filter panel should be \"malt, hop\"", () => readingReads(page, "categories of Aroma", el("filter panel"), "malt, hop"));
      await session.step(170, "And mode of \"Aroma\" filter card should have text \"OR\"", () => shouldHaveText(page, el("mode of \"Aroma\" filter card"), "OR"));
      await session.step(171, "And 106 rows should pass the filter", () => filterPasses(page, 106));
      await session.step(172, "When user clicks on mode of \"Aroma\" filter card", () => clickOn(page, el("mode of \"Aroma\" filter card")));
      await session.step(173, "Then mode of \"Aroma\" filter card should have text \"AND\"", () => shouldHaveText(page, el("mode of \"Aroma\" filter card"), "AND"));
      await session.step(174, "And 88 rows should pass the filter", () => filterPasses(page, 88));
      await session.step(175, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("At fuzziness 0 a near miss keeps nothing, and a higher fuzziness brings rows in", async () => {
      await session.step(178, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(179, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(180, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(181, "When user types \"maltx\" into the search of the \"Aroma\" filter card", () => typeIntoCardSearch(page, "maltx", "Aroma"));
      await session.step(182, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(183, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(184, "And Fuzzyness input in \"Aroma\" filter card should have the value \"0\"", () => shouldHaveValue(page, el("Fuzzyness input in \"Aroma\" filter card"), "0"));
      await session.step(185, "When user drags the slider of Fuzzyness input in \"Aroma\" filter card to 0.5", () => dragSliderTo(page, el("Fuzzyness input in \"Aroma\" filter card"), 0.5));
      await session.step(186, "Then the \"rows shown\" reading of filter panel should be at least 1", () => readingAtLeast(page, "rows shown", el("filter panel"), 1));
      await session.step(187, "When user remembers the \"rows shown\" reading of filter panel", () => rememberReading(page, "rows shown", el("filter panel")));
      await session.step(188, "And user drags the slider of Fuzzyness input in \"Aroma\" filter card to 1", () => dragSliderTo(page, el("Fuzzyness input in \"Aroma\" filter card"), 1));
      await session.step(189, "Then the \"rows shown\" reading of filter panel should be higher than remembered", () => readingHigherThanRemembered(page, "rows shown", el("filter panel")));
      await session.step(190, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A free-text rule that does not parse is refused on Enter, and the rows stay", async () => {
      await session.step(194, "When user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(195, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(196, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(197, "And user hovers over \"Expression\" filter card", () => hoverOver(page, el("\"Expression\" filter card")));
      await session.step(198, "And user clicks on \"Switch to free-text mode\" icon in \"Expression\" filter card", () => clickOn(page, el("\"Switch to free-text mode\" icon in \"Expression\" filter card")));
      await session.step(199, "Then Column input in \"Expression\" filter card should be hidden", () => shouldBe(page, el("Column input in \"Expression\" filter card"), "hidden"));
      await session.step(200, "When user types \"${HEIGHT} <\" into the search of the \"Expression\" filter card", () => typeIntoCardSearch(page, "${HEIGHT} <", "Expression"));
      await session.step(201, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(202, "And user types \"${HEIGHT} < 150\" into the search of the \"Expression\" filter card", () => typeIntoCardSearch(page, "${HEIGHT} < 150", "Expression"));
      await session.step(203, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(204, "Then the \"categories of Expression\" reading of filter panel should be \"${HEIGHT} < 150\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${HEIGHT} < 150"));
      await session.step(205, "And 16 rows should pass the filter", () => filterPasses(page, 16));
      await session.step(206, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(207, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(208, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The is null and is not null operations of a numeric column keep the blank and the filled rows", async () => {
      await session.step(212, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(213, "And user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
      await session.step(214, "And user selects \"is not null\" in Operation input in \"Expression\" filter card", () => selectIn(page, "is not null", el("Operation input in \"Expression\" filter card")));
      await session.step(215, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(216, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
      await session.step(217, "When user picks \"Remove Query\" from the context menu of the \"category ${HEIGHT} is not null of Expression\" area of filter panel", () => pickFromAreaContextMenu(page, "Remove Query", "category ${HEIGHT} is not null of Expression", el("filter panel")));
      await session.step(218, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(219, "When user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
      await session.step(220, "And user selects \"is null\" in Operation input in \"Expression\" filter card", () => selectIn(page, "is null", el("Operation input in \"Expression\" filter card")));
      await session.step(221, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(222, "Then 128 rows should pass the filter", () => filterPasses(page, 128));
      await session.step(223, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(224, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(225, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
