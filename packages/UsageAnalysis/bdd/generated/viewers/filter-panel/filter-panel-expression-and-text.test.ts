/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-expression-and-text.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/nx.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {typeIntoCardSearch} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, dragSliderTo, hoverOver, pasteInto, pressKey, selectIn, shouldBe, shouldHaveText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactly, filterIsExactlyCategory, filterIsExactlyContains, filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickPanelMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, noErrors, pickFromAreaContextMenu, readingAtLeast, readingHigherThanRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Expression and text filter cards", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-expression-and-text.feature", import.meta.url);
  test("Expression and text filter cards", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(27, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(28, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
    await session.step(29, "Then \"Expression\" filter card should be visible", () => shouldBe(page, el("\"Expression\" filter card"), "visible"));
    await session.step(30, "And the \"type of Expression\" reading of filter panel should be \"expression\"", () => readingReads(page, "type of Expression", el("filter panel"), "expression"));
    await session.step(31, "And all rows should pass the filter", () => filterPassesAll(page));
    await run.scenario("Two rules from the form combine with OR, and the switch makes them AND", async () => {
      await session.step(34, "When user selects \"AGE\" in Column input in \"Expression\" filter card", () => selectIn(page, "AGE", el("Column input in \"Expression\" filter card")));
      await session.step(35, "And user selects \">\" in Operation input in \"Expression\" filter card", () => selectIn(page, ">", el("Operation input in \"Expression\" filter card")));
      await session.step(36, "And user types \"50\" into Value input in \"Expression\" filter card", () => typeInto(page, "50", el("Value input in \"Expression\" filter card")));
      await session.step(37, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(38, "Then 367 rows should pass the filter", () => filterPasses(page, 367));
      await session.step(39, "And the filter should pass exactly the rows where \"AGE\" is between 51 and 89", () => filterIsExactly(page, "AGE", 51, 89));
      await session.step(40, "And the \"categories of Expression\" reading of filter panel should be \"${AGE} > 50\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${AGE} > 50"));
      await session.step(41, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(42, "When user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
      await session.step(43, "And user selects \"<\" in Operation input in \"Expression\" filter card", () => selectIn(page, "<", el("Operation input in \"Expression\" filter card")));
      await session.step(44, "And user types \"160\" into Value input in \"Expression\" filter card", () => typeInto(page, "160", el("Value input in \"Expression\" filter card")));
      await session.step(45, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(46, "Then the \"categories of Expression\" reading of filter panel should be \"${AGE} > 50, ${HEIGHT} < 160\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${AGE} > 50, ${HEIGHT} < 160"));
      await session.step(47, "And mode of \"Expression\" filter card should have text \"OR\"", () => shouldHaveText(page, el("mode of \"Expression\" filter card"), "OR"));
      await session.step(48, "And 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(49, "When user clicks on mode of \"Expression\" filter card", () => clickOn(page, el("mode of \"Expression\" filter card")));
      await session.step(50, "Then mode of \"Expression\" filter card should have text \"AND\"", () => shouldHaveText(page, el("mode of \"Expression\" filter card"), "AND"));
      await session.step(51, "And 88 rows should pass the filter", () => filterPasses(page, 88));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A rule's own checkbox suspends it and keeps it in the list", async () => {
      await session.step(55, "When user clicks on the \"checkbox ${HEIGHT} < 160 of Expression\" area of filter panel", () => clickArea(page, "checkbox ${HEIGHT} < 160 of Expression", el("filter panel")));
      await session.step(56, "Then 367 rows should pass the filter", () => filterPasses(page, 367));
      await session.step(57, "And the \"categories of Expression\" reading of filter panel should be \"${AGE} > 50, ${HEIGHT} < 160\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${AGE} > 50, ${HEIGHT} < 160"));
      await session.step(58, "And the \"selected categories of Expression\" reading of filter panel should be \"${AGE} > 50\"", () => readingReads(page, "selected categories of Expression", el("filter panel"), "${AGE} > 50"));
      await session.step(59, "When user clicks on the \"checkbox ${HEIGHT} < 160 of Expression\" area of filter panel", () => clickArea(page, "checkbox ${HEIGHT} < 160 of Expression", el("filter panel")));
      await session.step(60, "Then 88 rows should pass the filter", () => filterPasses(page, 88));
      await session.step(61, "And the \"selected categories of Expression\" reading of filter panel should be \"${AGE} > 50, ${HEIGHT} < 160\"", () => readingReads(page, "selected categories of Expression", el("filter panel"), "${AGE} > 50, ${HEIGHT} < 160"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Remove Query in a rule's menu takes the rule out", async () => {
      await session.step(65, "When user picks \"Remove Query\" from the context menu of the \"category ${AGE} > 50 of Expression\" area of filter panel", () => pickFromAreaContextMenu(page, "Remove Query", "category ${AGE} > 50 of Expression", el("filter panel")));
      await session.step(66, "Then 167 rows should pass the filter", () => filterPasses(page, 167));
      await session.step(67, "And the \"categories of Expression\" reading of filter panel should be \"${HEIGHT} < 160\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${HEIGHT} < 160"));
      await session.step(68, "And the filter should pass exactly the rows where \"HEIGHT\" is between 1 and 159.999", () => filterIsExactly(page, "HEIGHT", 1, 159.999));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Free-text mode takes a rule in the card's own syntax and the form comes back", async () => {
      await session.step(72, "When user hovers over \"Expression\" filter card", () => hoverOver(page, el("\"Expression\" filter card")));
      await session.step(73, "And user clicks on \"Switch to free-text mode\" icon in \"Expression\" filter card", () => clickOn(page, el("\"Switch to free-text mode\" icon in \"Expression\" filter card")));
      await session.step(74, "Then Column input in \"Expression\" filter card should be hidden", () => shouldBe(page, el("Column input in \"Expression\" filter card"), "hidden"));
      await session.step(75, "When user types \"${HEIGHT} < 150\" into the search of the \"Expression\" filter card", () => typeIntoCardSearch(page, "${HEIGHT} < 150", "Expression"));
      await session.step(76, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(77, "Then the \"categories of Expression\" reading of filter panel should be \"${HEIGHT} < 160, ${HEIGHT} < 150\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${HEIGHT} < 160, ${HEIGHT} < 150"));
      await session.step(78, "And 16 rows should pass the filter", () => filterPasses(page, 16));
      await session.step(79, "When user hovers over \"Expression\" filter card", () => hoverOver(page, el("\"Expression\" filter card")));
      await session.step(80, "And user clicks on \"Switch to free-text mode\" icon in \"Expression\" filter card", () => clickOn(page, el("\"Switch to free-text mode\" icon in \"Expression\" filter card")));
      await session.step(81, "Then Column input in \"Expression\" filter card should be visible", () => shouldBe(page, el("Column input in \"Expression\" filter card"), "visible"));
      await session.step(82, "And 16 rows should pass the filter", () => filterPasses(page, 16));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The regex operation turns a pasted list into an alternation and drops a trailing comma", async () => {
      await session.step(86, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(87, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(88, "And user selects \"USUBJID\" in Column input in \"Expression\" filter card", () => selectIn(page, "USUBJID", el("Column input in \"Expression\" filter card")));
      await session.step(89, "And user selects \"regex\" in Operation input in \"Expression\" filter card", () => selectIn(page, "regex", el("Operation input in \"Expression\" filter card")));
      await session.step(90, "And user pastes \"5,15,25\" into Value input in \"Expression\" filter card", () => pasteInto(page, "5,15,25", el("Value input in \"Expression\" filter card")));
      await session.step(91, "Then Value input in \"Expression\" filter card should have the value \"5|15|25\"", () => shouldHaveValue(page, el("Value input in \"Expression\" filter card"), "5|15|25"));
      await session.step(92, "When user clears Value input in \"Expression\" filter card", () => clearField(page, el("Value input in \"Expression\" filter card")));
      await session.step(93, "Then Value input in \"Expression\" filter card should have the value \"\"", () => shouldHaveValue(page, el("Value input in \"Expression\" filter card"), ""));
      await session.step(94, "When user pastes \"5,15,25,\" into Value input in \"Expression\" filter card", () => pasteInto(page, "5,15,25,", el("Value input in \"Expression\" filter card")));
      await session.step(95, "Then Value input in \"Expression\" filter card should have the value \"5|15|25\"", () => shouldHaveValue(page, el("Value input in \"Expression\" filter card"), "5|15|25"));
      await session.step(96, "When user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(97, "Then the \"categories of Expression\" reading of filter panel should be \"${USUBJID} regex 5|15|25\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${USUBJID} regex 5|15|25"));
      await session.step(98, "And 357 rows should pass the filter", () => filterPasses(page, 357));
      await session.step(99, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The string and date operations keep the rows they name", async () => {
      await session.step(102, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(103, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(104, "And user selects \"SEX\" in Column input in \"Expression\" filter card", () => selectIn(page, "SEX", el("Column input in \"Expression\" filter card")));
      await session.step(105, "And user selects \"equals\" in Operation input in \"Expression\" filter card", () => selectIn(page, "equals", el("Operation input in \"Expression\" filter card")));
      await session.step(106, "And user types \"F\" into Value input in \"Expression\" filter card", () => typeInto(page, "F", el("Value input in \"Expression\" filter card")));
      await session.step(107, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(108, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(109, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
      await session.step(110, "And the \"categories of Expression\" reading of filter panel should be \"${SEX} equals F\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${SEX} equals F"));
      await session.step(111, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(112, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(113, "And user selects \"RACE\" in Column input in \"Expression\" filter card", () => selectIn(page, "RACE", el("Column input in \"Expression\" filter card")));
      await session.step(114, "And user selects \"contains\" in Operation input in \"Expression\" filter card", () => selectIn(page, "contains", el("Operation input in \"Expression\" filter card")));
      await session.step(115, "And user types \"an\" into Value input in \"Expression\" filter card", () => typeInto(page, "an", el("Value input in \"Expression\" filter card")));
      await session.step(116, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(117, "Then 911 rows should pass the filter", () => filterPasses(page, 911));
      await session.step(118, "And the filter should pass exactly the rows where \"RACE\" contains \"an\"", () => filterIsExactlyContains(page, "RACE", "an"));
      await session.step(119, "And the \"categories of Expression\" reading of filter panel should be \"${RACE} contains an\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${RACE} contains an"));
      await session.step(120, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(121, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(122, "And user selects \"STARTED\" in Column input in \"Expression\" filter card", () => selectIn(page, "STARTED", el("Column input in \"Expression\" filter card")));
      await session.step(123, "And user selects \"after\" in Operation input in \"Expression\" filter card", () => selectIn(page, "after", el("Operation input in \"Expression\" filter card")));
      await session.step(124, "And user types \"01/01/1991\" into Value input in \"Expression\" filter card", () => typeInto(page, "01/01/1991", el("Value input in \"Expression\" filter card")));
      await session.step(125, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(126, "Then 483 rows should pass the filter", () => filterPasses(page, 483));
      await session.step(127, "And the \"categories of Expression\" reading of filter panel should be \"${STARTED} after 01/01/1991\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${STARTED} after 01/01/1991"));
      await session.step(128, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(129, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The text card keeps the rows holding any of its terms, or all of them", async () => {
      await session.step(133, "When user opens beer dataset", () => openDataset(page, ds("beer")));
      await session.step(134, "And user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(135, "Then \"Aroma\" filter card should be visible", () => shouldBe(page, el("\"Aroma\" filter card"), "visible"));
      await session.step(136, "And the \"type of Aroma\" reading of filter panel should be \"text\"", () => readingReads(page, "type of Aroma", el("filter panel"), "text"));
      await session.step(137, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(138, "When user types \"malt\" into the search of the \"Aroma\" filter card", () => typeIntoCardSearch(page, "malt", "Aroma"));
      await session.step(139, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(140, "Then 92 rows should pass the filter", () => filterPasses(page, 92));
      await session.step(141, "And the \"categories of Aroma\" reading of filter panel should be \"malt\"", () => readingReads(page, "categories of Aroma", el("filter panel"), "malt"));
      await session.step(142, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(143, "When user types \"hop\" into the search of the \"Aroma\" filter card", () => typeIntoCardSearch(page, "hop", "Aroma"));
      await session.step(144, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(145, "Then the \"categories of Aroma\" reading of filter panel should be \"malt, hop\"", () => readingReads(page, "categories of Aroma", el("filter panel"), "malt, hop"));
      await session.step(146, "And mode of \"Aroma\" filter card should have text \"OR\"", () => shouldHaveText(page, el("mode of \"Aroma\" filter card"), "OR"));
      await session.step(147, "And 106 rows should pass the filter", () => filterPasses(page, 106));
      await session.step(148, "When user clicks on mode of \"Aroma\" filter card", () => clickOn(page, el("mode of \"Aroma\" filter card")));
      await session.step(149, "Then mode of \"Aroma\" filter card should have text \"AND\"", () => shouldHaveText(page, el("mode of \"Aroma\" filter card"), "AND"));
      await session.step(150, "And 88 rows should pass the filter", () => filterPasses(page, 88));
      await session.step(151, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("At fuzziness 0 a near miss keeps nothing, and a higher fuzziness brings rows in", async () => {
      await session.step(154, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(155, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(156, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(157, "When user types \"maltx\" into the search of the \"Aroma\" filter card", () => typeIntoCardSearch(page, "maltx", "Aroma"));
      await session.step(158, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(159, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(160, "And Fuzzyness input in \"Aroma\" filter card should have the value \"0\"", () => shouldHaveValue(page, el("Fuzzyness input in \"Aroma\" filter card"), "0"));
      await session.step(161, "When user drags the slider of Fuzzyness input in \"Aroma\" filter card to 0.5", () => dragSliderTo(page, el("Fuzzyness input in \"Aroma\" filter card"), 0.5));
      await session.step(162, "Then the \"rows shown\" reading of filter panel should be at least 1", () => readingAtLeast(page, "rows shown", el("filter panel"), 1));
      await session.step(163, "When user remembers the \"rows shown\" reading of filter panel", () => rememberReading(page, "rows shown", el("filter panel")));
      await session.step(164, "And user drags the slider of Fuzzyness input in \"Aroma\" filter card to 1", () => dragSliderTo(page, el("Fuzzyness input in \"Aroma\" filter card"), 1));
      await session.step(165, "Then the \"rows shown\" reading of filter panel should be higher than remembered", () => readingHigherThanRemembered(page, "rows shown", el("filter panel")));
      await session.step(166, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
