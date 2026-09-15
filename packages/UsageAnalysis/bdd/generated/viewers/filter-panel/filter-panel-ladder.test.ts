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
import {addCardFor} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, hoverOver, shouldBe, shouldContainText, shouldHaveText, shouldNotContainText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addRangeFilter, filterIsExactlyCategory, filterPanelCount, filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel core ladder", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-ladder.feature", import.meta.url);
  test("Filter panel core ladder", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(16, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(17, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(18, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
    await run.scenario("A card added from the header picker filters nothing", async () => {
      await session.step(21, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(22, "Then \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(23, "And the \"cards\" reading of filter panel should be \"RACE\"", () => readingReads(page, "cards", el("filter panel"), "RACE"));
      await session.step(24, "And the \"type of RACE\" reading of filter panel should be \"categorical\"", () => readingReads(page, "type of RACE", el("filter panel"), "categorical"));
      await session.step(25, "And the \"categories of RACE\" reading of filter panel should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "categories of RACE", el("filter panel"), "Asian, Black, Caucasian, Other"));
      await session.step(26, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(27, "And the \"filters\" reading of filter panel should be 0", () => readingIs(page, "filters", el("filter panel"), 0));
      await session.step(28, "And the \"filtering of RACE\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of RACE", el("filter panel"), "false"));
      await session.step(29, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a category name keeps that category alone", async () => {
      await session.step(33, "When user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(34, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(35, "And the filter should pass exactly the rows where \"RACE\" is \"Black\"", () => filterIsExactlyCategory(page, "RACE", "Black"));
      await session.step(36, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(37, "And the \"rows shown\" reading of filter panel should be 27", () => readingIs(page, "rows shown", el("filter panel"), 27));
      await session.step(38, "And the \"filters\" reading of filter panel should be 1", () => readingIs(page, "filters", el("filter panel"), 1));
      await session.step(39, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a checkbox adds a category to the ones kept", async () => {
      await session.step(43, "When user clicks on the \"checkbox Other of RACE\" area of filter panel", () => clickArea(page, "checkbox Other of RACE", el("filter panel")));
      await session.step(44, "Then 89 rows should pass the filter", () => filterPasses(page, 89));
      await session.step(45, "And the \"selected categories of RACE\" reading of filter panel should be \"Black, Other\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black, Other"));
      await session.step(46, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(47, "When user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(48, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(49, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second criterion composes with the first", async () => {
      await session.step(53, "When user adds a range filter on \"AGE\" from 30 to 60", () => addRangeFilter(page, "AGE", 30, 60));
      await session.step(54, "Then \"AGE\" filter card should be visible", () => shouldBe(page, el("\"AGE\" filter card"), "visible"));
      await session.step(55, "And 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(56, "And the \"min of AGE\" reading of filter panel should be 30", () => readingIs(page, "min of AGE", el("filter panel"), 30));
      await session.step(57, "And the \"max of AGE\" reading of filter panel should be 60", () => readingIs(page, "max of AGE", el("filter panel"), 60));
      await session.step(58, "And the \"filters\" reading of filter panel should be 2", () => readingIs(page, "filters", el("filter panel"), 2));
      await session.step(59, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The counter's tooltip names the cards that restrict rows", async () => {
      await session.step(63, "When user adds a card for \"SEX\" to the filter panel", () => addCardFor(page, "SEX"));
      await session.step(64, "Then the \"filters\" reading of filter panel should be 2", () => readingIs(page, "filters", el("filter panel"), 2));
      await session.step(65, "When user hovers over counter of filter panel", () => hoverOver(page, el("counter of filter panel")));
      await session.step(66, "Then tooltip should contain the text \"RACE\"", () => shouldContainText(page, el("tooltip"), "RACE"));
      await session.step(67, "And tooltip should contain the text \"Black\"", () => shouldContainText(page, el("tooltip"), "Black"));
      await session.step(68, "And tooltip should contain the text \"AGE\"", () => shouldContainText(page, el("tooltip"), "AGE"));
      await session.step(69, "And tooltip should not contain the text \"SEX\"", () => shouldNotContainText(page, el("tooltip"), "SEX"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The master toggle stashes every card's state and gives it back", async () => {
      await session.step(73, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(74, "And user unchecks master of filter panel", () => uncheck(page, el("master of filter panel")));
      await session.step(75, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(76, "And the \"active\" reading of filter panel should be \"false\"", () => readingReads(page, "active", el("filter panel"), "false"));
      await session.step(77, "And the \"cards\" reading of filter panel should be \"SEX, AGE, RACE\"", () => readingReads(page, "cards", el("filter panel"), "SEX, AGE, RACE"));
      await session.step(78, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(79, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(80, "When user checks master of filter panel", () => check(page, el("master of filter panel")));
      await session.step(81, "Then 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(82, "And the \"active\" reading of filter panel should be \"true\"", () => readingReads(page, "active", el("filter panel"), "true"));
      await session.step(83, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The header search hides cards and leaves the rows alone", async () => {
      await session.step(87, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(88, "And user clicks on search icon of filter panel", () => clickOn(page, el("search icon of filter panel")));
      await session.step(89, "And user types \"RACE\" into search of filter panel", () => typeInto(page, "RACE", el("search of filter panel")));
      await session.step(90, "Then \"AGE\" filter card should be hidden", () => shouldBe(page, el("\"AGE\" filter card"), "hidden"));
      await session.step(91, "And \"SEX\" filter card should be hidden", () => shouldBe(page, el("\"SEX\" filter card"), "hidden"));
      await session.step(92, "And \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(93, "And 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(94, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(95, "When user clears search of filter panel", () => clearField(page, el("search of filter panel")));
      await session.step(96, "Then \"AGE\" filter card should be visible", () => shouldBe(page, el("\"AGE\" filter card"), "visible"));
      await session.step(97, "And \"SEX\" filter card should be visible", () => shouldBe(page, el("\"SEX\" filter card"), "visible"));
      await session.step(98, "And 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(99, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The reset icon clears the criteria and keeps the cards", async () => {
      await session.step(102, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(103, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(104, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(105, "And the \"filters\" reading of filter panel should be 0", () => readingIs(page, "filters", el("filter panel"), 0));
      await session.step(106, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(107, "And the \"cards\" reading of filter panel should be \"SEX, AGE, RACE\"", () => readingReads(page, "cards", el("filter panel"), "SEX, AGE, RACE"));
      await session.step(108, "And \"RACE\" filter card should be enabled", () => shouldBe(page, el("\"RACE\" filter card"), "enabled"));
      await session.step(109, "And the \"filtering of RACE\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of RACE", el("filter panel"), "false"));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A card's checkbox suspends its criterion and keeps it", async () => {
      await session.step(113, "When user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(114, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(115, "When user hovers over \"RACE\" filter card", () => hoverOver(page, el("\"RACE\" filter card")));
      await session.step(116, "And user unchecks checkbox of \"RACE\" filter card", () => uncheck(page, el("checkbox of \"RACE\" filter card")));
      await session.step(117, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(118, "And \"RACE\" filter card should be disabled", () => shouldBe(page, el("\"RACE\" filter card"), "disabled"));
      await session.step(119, "And the \"enabled of RACE\" reading of filter panel should be \"false\"", () => readingReads(page, "enabled of RACE", el("filter panel"), "false"));
      await session.step(120, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(121, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(122, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(123, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(124, "Then \"RACE\" filter card should be enabled", () => shouldBe(page, el("\"RACE\" filter card"), "enabled"));
      await session.step(125, "And the \"enabled of RACE\" reading of filter panel should be \"true\"", () => readingReads(page, "enabled of RACE", el("filter panel"), "true"));
      await session.step(126, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
