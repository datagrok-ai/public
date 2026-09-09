/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-entry-points.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addCardFor, pickCounterMenu, pickPanelMenu} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPanelCount, filterPasses, filterPassesAll, openEmptyFilterPanel, openFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel entry points", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-entry-points.feature", import.meta.url);
  test("Filter panel entry points", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(13, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(14, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(15, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
    await run.scenario("Cards added from the header picker stack at the top", async () => {
      await session.step(18, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(19, "And user adds a card for \"SEX\" to the filter panel", () => addCardFor(page, "SEX"));
      await session.step(20, "Then the filter panel should have 2 filters", () => filterPanelCount(page, 2));
      await session.step(21, "And the \"cards\" reading of filter panel should be \"SEX, RACE\"", () => readingReads(page, "cards", el("filter panel"), "SEX, RACE"));
      await session.step(22, "And \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(23, "And \"SEX\" filter card should be visible", () => shouldBe(page, el("\"SEX\" filter card"), "visible"));
      await session.step(24, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(25, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A card added from the panel's own menu goes to the top too", async () => {
      await session.step(29, "When user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
      await session.step(30, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
      await session.step(31, "And the filter panel should have 3 filters", () => filterPanelCount(page, 3));
      await session.step(32, "And the \"cards\" reading of filter panel should be \"Flags, SEX, RACE\"", () => readingReads(page, "cards", el("filter panel"), "Flags, SEX, RACE"));
      await session.step(33, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(34, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A card's close icon removes it and releases its criterion", async () => {
      await session.step(38, "When user clicks on the \"category Caucasian of RACE\" area of filter panel", () => clickArea(page, "category Caucasian of RACE", el("filter panel")));
      await session.step(39, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(40, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(41, "When user hovers over \"RACE\" filter card", () => hoverOver(page, el("\"RACE\" filter card")));
      await session.step(42, "And user clicks on close of \"RACE\" filter card", () => clickOn(page, el("close of \"RACE\" filter card")));
      await session.step(43, "Then \"RACE\" filter card should be absent", () => shouldBe(page, el("\"RACE\" filter card"), "absent"));
      await session.step(44, "And the filter panel should have 2 filters", () => filterPanelCount(page, 2));
      await session.step(45, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(46, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Remove others keeps the cards that restrict rows, Remove All empties the panel", async () => {
      await session.step(50, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(51, "And user clicks on the \"category Caucasian of RACE\" area of filter panel", () => clickArea(page, "category Caucasian of RACE", el("filter panel")));
      await session.step(52, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(53, "And the filter panel should have 3 filters", () => filterPanelCount(page, 3));
      await session.step(54, "When user picks \"Remove others\" from the filter counter menu", () => pickCounterMenu(page, "Remove others"));
      await session.step(55, "Then \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(56, "And \"SEX\" filter card should be absent", () => shouldBe(page, el("\"SEX\" filter card"), "absent"));
      await session.step(57, "And \"Flags\" filter card should be absent", () => shouldBe(page, el("\"Flags\" filter card"), "absent"));
      await session.step(58, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(59, "And 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(60, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(61, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(62, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(63, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the panel releases its filtering and reopening restores the card", async () => {
      await session.step(67, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(68, "And user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(69, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(70, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(71, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(72, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(73, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(74, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(75, "And \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(76, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(77, "And 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(78, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(79, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
