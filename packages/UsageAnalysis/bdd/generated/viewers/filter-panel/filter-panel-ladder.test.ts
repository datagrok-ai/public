/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-ladder.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {enterCardBound, pickCardIndicatorMenu} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, hoverOver, pressKey, shouldBe, shouldContainText, shouldHaveText, shouldNotContainText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addRangeFilter, filterIsExactlyCategory, filterPanelCount, filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, closeContextMenu, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickEmptySpace, dragAreaOntoWidget} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel core ladder", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-ladder.feature", import.meta.url);
  test("Filter panel core ladder", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 12, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(23, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(24, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(25, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
    await run.scenario("A card dragged from a grid header filters nothing", async () => {
      await session.step(28, "When user drags the \"header RACE\" area of grid onto the \"view\" area of filter panel", () => dragAreaOntoWidget(page, "header RACE", el("grid"), "view", el("filter panel")));
      await session.step(29, "Then \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(30, "And the \"cards\" reading of filter panel should be \"RACE\"", () => readingReads(page, "cards", el("filter panel"), "RACE"));
      await session.step(31, "And the \"type of RACE\" reading of filter panel should be \"categorical\"", () => readingReads(page, "type of RACE", el("filter panel"), "categorical"));
      await session.step(32, "And the \"categories of RACE\" reading of filter panel should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "categories of RACE", el("filter panel"), "Asian, Black, Caucasian, Other"));
      await session.step(33, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(34, "And the \"filters\" reading of filter panel should be 0", () => readingIs(page, "filters", el("filter panel"), 0));
      await session.step(35, "And the \"filtering of RACE\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of RACE", el("filter panel"), "false"));
      await session.step(36, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a category name keeps that category alone", async () => {
      await session.step(40, "When user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(41, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(42, "And the filter should pass exactly the rows where \"RACE\" is \"Black\"", () => filterIsExactlyCategory(page, "RACE", "Black"));
      await session.step(43, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(44, "And the \"rows shown\" reading of filter panel should be 27", () => readingIs(page, "rows shown", el("filter panel"), 27));
      await session.step(45, "And the \"filters\" reading of filter panel should be 1", () => readingIs(page, "filters", el("filter panel"), 1));
      await session.step(46, "And counter of filter panel should be visible", () => shouldBe(page, el("counter of filter panel"), "visible"));
      await session.step(47, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a checkbox adds a category to the ones kept", async () => {
      await session.step(51, "When user clicks on the \"checkbox Other of RACE\" area of filter panel", () => clickArea(page, "checkbox Other of RACE", el("filter panel")));
      await session.step(52, "Then 89 rows should pass the filter", () => filterPasses(page, 89));
      await session.step(53, "And the \"selected categories of RACE\" reading of filter panel should be \"Black, Other\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black, Other"));
      await session.step(54, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(55, "When user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(56, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(57, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second criterion composes with the first", async () => {
      await session.step(61, "When user adds a range filter on \"AGE\" from 30 to 60", () => addRangeFilter(page, "AGE", 30, 60));
      await session.step(62, "Then \"AGE\" filter card should be visible", () => shouldBe(page, el("\"AGE\" filter card"), "visible"));
      await session.step(63, "And 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(64, "And the \"min of AGE\" reading of filter panel should be 30", () => readingIs(page, "min of AGE", el("filter panel"), 30));
      await session.step(65, "And the \"max of AGE\" reading of filter panel should be 60", () => readingIs(page, "max of AGE", el("filter panel"), 60));
      await session.step(66, "And the \"filters\" reading of filter panel should be 2", () => readingIs(page, "filters", el("filter panel"), 2));
      await session.step(67, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The max field still answers after an end out of the column's range", async () => {
      await session.step(71, "When user picks \"Min / max\" from the indicator menu of the \"AGE\" filter card", () => pickCardIndicatorMenu(page, "Min / max", "AGE"));
      await session.step(72, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(73, "And user enters \"55\" into the max field of the \"AGE\" filter card", () => enterCardBound(page, "55", "max", "AGE"));
      await session.step(74, "Then 18 rows should pass the filter", () => filterPasses(page, 18));
      await session.step(75, "And the \"max of AGE\" reading of filter panel should be 55", () => readingIs(page, "max of AGE", el("filter panel"), 55));
      await session.step(76, "When user enters \"999\" into the max field of the \"AGE\" filter card", () => enterCardBound(page, "999", "max", "AGE"));
      await session.step(77, "Then 25 rows should pass the filter", () => filterPasses(page, 25));
      await session.step(78, "And the \"max of AGE\" reading of filter panel should be 89", () => readingIs(page, "max of AGE", el("filter panel"), 89));
      await session.step(79, "When user enters \"60\" into the max field of the \"AGE\" filter card", () => enterCardBound(page, "60", "max", "AGE"));
      await session.step(80, "Then 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(81, "And the \"max of AGE\" reading of filter panel should be 60", () => readingIs(page, "max of AGE", el("filter panel"), 60));
      await session.step(82, "And the \"min of AGE\" reading of filter panel should be 30", () => readingIs(page, "min of AGE", el("filter panel"), 30));
      await session.step(83, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The counter's tooltip names the cards that restrict rows", async () => {
      await session.step(87, "When user adds a card for \"SEX\" to the filter panel", () => addCardFor(page, "SEX"));
      await session.step(88, "Then the \"filters\" reading of filter panel should be 2", () => readingIs(page, "filters", el("filter panel"), 2));
      await session.step(89, "When user hovers over counter of filter panel", () => hoverOver(page, el("counter of filter panel")));
      await session.step(90, "Then tooltip should contain the text \"RACE\"", () => shouldContainText(page, el("tooltip"), "RACE"));
      await session.step(91, "And tooltip should contain the text \"Black\"", () => shouldContainText(page, el("tooltip"), "Black"));
      await session.step(92, "And tooltip should contain the text \"AGE\"", () => shouldContainText(page, el("tooltip"), "AGE"));
      await session.step(93, "And tooltip should contain the text \"[30,60]\"", () => shouldContainText(page, el("tooltip"), "[30,60]"));
      await session.step(94, "And tooltip should not contain the text \"Caucasian\"", () => shouldNotContainText(page, el("tooltip"), "Caucasian"));
      await session.step(95, "And tooltip should not contain the text \"Other\"", () => shouldNotContainText(page, el("tooltip"), "Other"));
      await session.step(96, "And tooltip should not contain the text \"SEX\"", () => shouldNotContainText(page, el("tooltip"), "SEX"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The master toggle stashes every card's state and gives it back", async () => {
      await session.step(100, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(101, "And user unchecks master of filter panel", () => uncheck(page, el("master of filter panel")));
      await session.step(102, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(103, "And the \"active\" reading of filter panel should be \"false\"", () => readingReads(page, "active", el("filter panel"), "false"));
      await session.step(104, "And the \"cards\" reading of filter panel should be \"SEX, AGE, RACE\"", () => readingReads(page, "cards", el("filter panel"), "SEX, AGE, RACE"));
      await session.step(105, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(106, "And \"AGE\" filter card should be disabled", () => shouldBe(page, el("\"AGE\" filter card"), "disabled"));
      await session.step(107, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(108, "When user checks master of filter panel", () => check(page, el("master of filter panel")));
      await session.step(109, "Then 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(110, "And counter of filter panel should be visible", () => shouldBe(page, el("counter of filter panel"), "visible"));
      await session.step(111, "And the \"active\" reading of filter panel should be \"true\"", () => readingReads(page, "active", el("filter panel"), "true"));
      await session.step(112, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The master toggle keeps a card switched off on its own switched off", async () => {
      await session.step(116, "When user hovers over \"RACE\" filter card", () => hoverOver(page, el("\"RACE\" filter card")));
      await session.step(117, "And user unchecks checkbox of \"RACE\" filter card", () => uncheck(page, el("checkbox of \"RACE\" filter card")));
      await session.step(118, "Then 708 rows should pass the filter", () => filterPasses(page, 708));
      await session.step(119, "And \"RACE\" filter card should be disabled", () => shouldBe(page, el("\"RACE\" filter card"), "disabled"));
      await session.step(120, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(121, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(122, "And user unchecks master of filter panel", () => uncheck(page, el("master of filter panel")));
      await session.step(123, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(124, "And the \"active\" reading of filter panel should be \"false\"", () => readingReads(page, "active", el("filter panel"), "false"));
      await session.step(125, "When user checks master of filter panel", () => check(page, el("master of filter panel")));
      await session.step(126, "Then 708 rows should pass the filter", () => filterPasses(page, 708));
      await session.step(127, "And the \"active\" reading of filter panel should be \"true\"", () => readingReads(page, "active", el("filter panel"), "true"));
      await session.step(128, "And \"RACE\" filter card should be disabled", () => shouldBe(page, el("\"RACE\" filter card"), "disabled"));
      await session.step(129, "And the \"enabled of RACE\" reading of filter panel should be \"false\"", () => readingReads(page, "enabled of RACE", el("filter panel"), "false"));
      await session.step(130, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(131, "When user hovers over \"RACE\" filter card", () => hoverOver(page, el("\"RACE\" filter card")));
      await session.step(132, "Then checkbox of \"RACE\" filter card should be unchecked", () => shouldBe(page, el("checkbox of \"RACE\" filter card"), "unchecked"));
      await session.step(133, "When user hovers over \"RACE\" filter card", () => hoverOver(page, el("\"RACE\" filter card")));
      await session.step(134, "And user checks checkbox of \"RACE\" filter card", () => check(page, el("checkbox of \"RACE\" filter card")));
      await session.step(135, "Then 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(136, "And \"RACE\" filter card should be enabled", () => shouldBe(page, el("\"RACE\" filter card"), "enabled"));
      await session.step(137, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(138, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(139, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Escape on the focused panel switches the cards off and gives them back unchanged", async () => {
      await session.step(142, "When user clicks on empty plot space of filter panel", () => clickEmptySpace(page, el("filter panel")));
      await session.step(143, "And user presses Escape", () => pressKey(page, "Escape"));
      await session.step(144, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(145, "And the \"active\" reading of filter panel should be \"false\"", () => readingReads(page, "active", el("filter panel"), "false"));
      await session.step(146, "And \"AGE\" filter card should be disabled", () => shouldBe(page, el("\"AGE\" filter card"), "disabled"));
      await session.step(147, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(148, "And the \"min of AGE\" reading of filter panel should be 30", () => readingIs(page, "min of AGE", el("filter panel"), 30));
      await session.step(149, "And the \"max of AGE\" reading of filter panel should be 60", () => readingIs(page, "max of AGE", el("filter panel"), 60));
      await session.step(150, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(151, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(152, "Then 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(153, "And the \"active\" reading of filter panel should be \"true\"", () => readingReads(page, "active", el("filter panel"), "true"));
      await session.step(154, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(155, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The header search hides cards and leaves the rows alone", async () => {
      await session.step(158, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(159, "And user clicks on search icon of filter panel", () => clickOn(page, el("search icon of filter panel")));
      await session.step(160, "And user types \"RACE\" into search of filter panel", () => typeInto(page, "RACE", el("search of filter panel")));
      await session.step(161, "Then \"AGE\" filter card should be hidden", () => shouldBe(page, el("\"AGE\" filter card"), "hidden"));
      await session.step(162, "And \"SEX\" filter card should be hidden", () => shouldBe(page, el("\"SEX\" filter card"), "hidden"));
      await session.step(163, "And \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(164, "And 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(165, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(166, "When user clears search of filter panel", () => clearField(page, el("search of filter panel")));
      await session.step(167, "Then \"AGE\" filter card should be visible", () => shouldBe(page, el("\"AGE\" filter card"), "visible"));
      await session.step(168, "And \"SEX\" filter card should be visible", () => shouldBe(page, el("\"SEX\" filter card"), "visible"));
      await session.step(169, "And 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(170, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The reset icon clears the criteria, keeps the cards and switches every one back on", async () => {
      await session.step(173, "When user hovers over \"AGE\" filter card", () => hoverOver(page, el("\"AGE\" filter card")));
      await session.step(174, "And user unchecks checkbox of \"AGE\" filter card", () => uncheck(page, el("checkbox of \"AGE\" filter card")));
      await session.step(175, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(176, "And \"AGE\" filter card should be disabled", () => shouldBe(page, el("\"AGE\" filter card"), "disabled"));
      await session.step(177, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(178, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(179, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(180, "And the \"filters\" reading of filter panel should be 0", () => readingIs(page, "filters", el("filter panel"), 0));
      await session.step(181, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(182, "And the \"cards\" reading of filter panel should be \"SEX, AGE, RACE\"", () => readingReads(page, "cards", el("filter panel"), "SEX, AGE, RACE"));
      await session.step(183, "And \"AGE\" filter card should be enabled", () => shouldBe(page, el("\"AGE\" filter card"), "enabled"));
      await session.step(184, "And \"RACE\" filter card should be enabled", () => shouldBe(page, el("\"RACE\" filter card"), "enabled"));
      await session.step(185, "And \"SEX\" filter card should be enabled", () => shouldBe(page, el("\"SEX\" filter card"), "enabled"));
      await session.step(186, "And the \"enabled of AGE\" reading of filter panel should be \"true\"", () => readingReads(page, "enabled of AGE", el("filter panel"), "true"));
      await session.step(187, "When user hovers over \"AGE\" filter card", () => hoverOver(page, el("\"AGE\" filter card")));
      await session.step(188, "Then checkbox of \"AGE\" filter card should be checked", () => shouldBe(page, el("checkbox of \"AGE\" filter card"), "checked"));
      await session.step(189, "When user hovers over \"RACE\" filter card", () => hoverOver(page, el("\"RACE\" filter card")));
      await session.step(190, "Then checkbox of \"RACE\" filter card should be checked", () => shouldBe(page, el("checkbox of \"RACE\" filter card"), "checked"));
      await session.step(191, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(192, "Then checkbox of \"SEX\" filter card should be checked", () => shouldBe(page, el("checkbox of \"SEX\" filter card"), "checked"));
      await session.step(193, "And the \"filtering of RACE\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of RACE", el("filter panel"), "false"));
      await session.step(194, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A card's checkbox suspends its criterion and keeps it", async () => {
      await session.step(197, "When user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(198, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(199, "When user hovers over \"RACE\" filter card", () => hoverOver(page, el("\"RACE\" filter card")));
      await session.step(200, "And user unchecks checkbox of \"RACE\" filter card", () => uncheck(page, el("checkbox of \"RACE\" filter card")));
      await session.step(201, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(202, "And \"RACE\" filter card should be disabled", () => shouldBe(page, el("\"RACE\" filter card"), "disabled"));
      await session.step(203, "And the \"enabled of RACE\" reading of filter panel should be \"false\"", () => readingReads(page, "enabled of RACE", el("filter panel"), "false"));
      await session.step(204, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(205, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(206, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(207, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(208, "Then \"RACE\" filter card should be enabled", () => shouldBe(page, el("\"RACE\" filter card"), "enabled"));
      await session.step(209, "And the \"enabled of RACE\" reading of filter panel should be \"true\"", () => readingReads(page, "enabled of RACE", el("filter panel"), "true"));
      await session.step(210, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(211, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
