/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-categorical.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addCardFor, clearCardSearch, openCardIndicatorMenu, pickCardIndicatorMenu, typeIntoCardSearch} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldContainText, shouldHaveText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactlyCategory, filterIsExactlyContains, filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, closeContextMenu, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Categorical filter card", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-categorical.feature", import.meta.url);
  test("Categorical filter card", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(14, "And user adds a card for \"DIS_POP\" to the filter panel", () => addCardFor(page, "DIS_POP"));
    await session.step(15, "Then the \"type of DIS_POP\" reading of filter panel should be \"categorical\"", () => readingReads(page, "type of DIS_POP", el("filter panel"), "categorical"));
    await session.step(16, "And the \"categories of DIS_POP\" reading of filter panel should be \"AS, Indigestion, PsA, Psoriasis, RA, UC\"", () => readingReads(page, "categories of DIS_POP", el("filter panel"), "AS, Indigestion, PsA, Psoriasis, RA, UC"));
    await session.step(17, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(18, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
    await run.scenario("A name click keeps one category, a checkbox click adds another", async () => {
      await session.step(21, "When user clicks on the \"category RA of DIS_POP\" area of filter panel", () => clickArea(page, "category RA of DIS_POP", el("filter panel")));
      await session.step(22, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(23, "And the filter should pass exactly the rows where \"DIS_POP\" is \"RA\"", () => filterIsExactlyCategory(page, "DIS_POP", "RA"));
      await session.step(24, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"RA\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "RA"));
      await session.step(25, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(26, "When user clicks on the \"checkbox UC of DIS_POP\" area of filter panel", () => clickArea(page, "checkbox UC of DIS_POP", el("filter panel")));
      await session.step(27, "Then 557 rows should pass the filter", () => filterPasses(page, 557));
      await session.step(28, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"RA, UC\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "RA, UC"));
      await session.step(29, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The indicator menu selects, deselects and inverts every category", async () => {
      await session.step(33, "When user picks \"Deselect all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Deselect all", "DIS_POP"));
      await session.step(34, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(35, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), ""));
      await session.step(36, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(37, "When user picks \"Invert all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Invert all", "DIS_POP"));
      await session.step(38, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(39, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"AS, Indigestion, PsA, Psoriasis, RA, UC\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "AS, Indigestion, PsA, Psoriasis, RA, UC"));
      await session.step(40, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(41, "When user picks \"Deselect all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Deselect all", "DIS_POP"));
      await session.step(42, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(43, "When user picks \"Select all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Select all", "DIS_POP"));
      await session.step(44, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(45, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Radio mode keeps exactly one category and offers no batch operations", async () => {
      await session.step(49, "When user clicks on the \"category RA of DIS_POP\" area of filter panel", () => clickArea(page, "category RA of DIS_POP", el("filter panel")));
      await session.step(50, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(51, "When user picks \"Mode | Radio\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Mode | Radio", "DIS_POP"));
      await session.step(52, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(53, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(54, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"RA\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "RA"));
      await session.step(55, "When user opens the indicator menu of the \"DIS_POP\" filter card", () => openCardIndicatorMenu(page, "DIS_POP"));
      await session.step(56, "Then context menu should contain the text \"Mode\"", () => shouldContainText(page, el("context menu"), "Mode"));
      await session.step(57, "And context menu should not contain the text \"Select all\"", () => shouldNotContainText(page, el("context menu"), "Select all"));
      await session.step(58, "And context menu should not contain the text \"Invert all\"", () => shouldNotContainText(page, el("context menu"), "Invert all"));
      await session.step(59, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(60, "And user clicks on the \"category UC of DIS_POP\" area of filter panel", () => clickArea(page, "category UC of DIS_POP", el("filter panel")));
      await session.step(61, "Then 123 rows should pass the filter", () => filterPasses(page, 123));
      await session.step(62, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"UC\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "UC"));
      await session.step(63, "When user picks \"Mode | Multi-Select\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Mode | Multi-Select", "DIS_POP"));
      await session.step(64, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(65, "Then 123 rows should pass the filter", () => filterPasses(page, 123));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The in-card search narrows the table to the categories that match", async () => {
      await session.step(69, "When user picks \"Select all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Select all", "DIS_POP"));
      await session.step(70, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(71, "When user clicks on search icon of \"DIS_POP\" filter card", () => clickOn(page, el("search icon of \"DIS_POP\" filter card")));
      await session.step(72, "And user types \"Ps\" into the search of the \"DIS_POP\" filter card", () => typeIntoCardSearch(page, "Ps", "DIS_POP"));
      await session.step(73, "Then 242 rows should pass the filter", () => filterPasses(page, 242));
      await session.step(74, "And the filter should pass exactly the rows where \"DIS_POP\" contains \"Ps\"", () => filterIsExactlyContains(page, "DIS_POP", "Ps"));
      await session.step(75, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(76, "When user clears the search of the \"DIS_POP\" filter card", () => clearCardSearch(page, "DIS_POP"));
      await session.step(77, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(78, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A numeric card switches to categorical and back without filtering", async () => {
      await session.step(82, "When user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
      await session.step(83, "Then the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
      await session.step(84, "When user hovers over \"AGE\" filter card", () => hoverOver(page, el("\"AGE\" filter card")));
      await session.step(85, "And user clicks on \"Switch to categorical filter\" icon in \"AGE\" filter card", () => clickOn(page, el("\"Switch to categorical filter\" icon in \"AGE\" filter card")));
      await session.step(86, "Then the \"type of AGE\" reading of filter panel should be \"categorical\"", () => readingReads(page, "type of AGE", el("filter panel"), "categorical"));
      await session.step(87, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(88, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(89, "When user hovers over \"AGE\" filter card", () => hoverOver(page, el("\"AGE\" filter card")));
      await session.step(90, "And user clicks on \"Switch to histogram filter\" icon in \"AGE\" filter card", () => clickOn(page, el("\"Switch to histogram filter\" icon in \"AGE\" filter card")));
      await session.step(91, "Then the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
      await session.step(92, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(93, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
