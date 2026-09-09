/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addCardFor} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, hoverOver, shouldBe, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addRangeFilter, filterPanelCount, filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, loadLayout, noErrors, readingIs, readingReads, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel persistence", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-persistence.feature", import.meta.url);
  test("Filter panel persistence", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(13, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(14, "And all rows should pass the filter", () => filterPassesAll(page));
    await run.scenario("A saved layout brings the cards and their criteria back", async () => {
      await session.step(17, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(18, "And user clicks on the \"category Caucasian of RACE\" area of filter panel", () => clickArea(page, "category Caucasian of RACE", el("filter panel")));
      await session.step(19, "And user adds a range filter on \"AGE\" from 30 to 60", () => addRangeFilter(page, "AGE", 30, 60));
      await session.step(20, "Then 633 rows should pass the filter", () => filterPasses(page, 633));
      await session.step(21, "And the \"filters\" reading of filter panel should be 2", () => readingIs(page, "filters", el("filter panel"), 2));
      await session.step(22, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(23, "And user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(24, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(25, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(26, "And the \"filters\" reading of filter panel should be 0", () => readingIs(page, "filters", el("filter panel"), 0));
      await session.step(27, "When user loads the saved layout", () => loadLayout(page));
      await session.step(28, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(29, "And \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(30, "And \"AGE\" filter card should be visible", () => shouldBe(page, el("\"AGE\" filter card"), "visible"));
      await session.step(31, "And the \"selected categories of RACE\" reading of filter panel should be \"Caucasian\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Caucasian"));
      await session.step(32, "And the \"min of AGE\" reading of filter panel should be 30", () => readingIs(page, "min of AGE", el("filter panel"), 30));
      await session.step(33, "And the \"max of AGE\" reading of filter panel should be 60", () => readingIs(page, "max of AGE", el("filter panel"), 60));
      await session.step(34, "And 633 rows should pass the filter", () => filterPasses(page, 633));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout saved while the header search was open carries no residue", async () => {
      await session.step(38, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(39, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(40, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(41, "When user clicks on search icon of filter panel", () => clickOn(page, el("search icon of filter panel")));
      await session.step(42, "And user types \"RACE\" into search of filter panel", () => typeInto(page, "RACE", el("search of filter panel")));
      await session.step(43, "Then \"AGE\" filter card should be hidden", () => shouldBe(page, el("\"AGE\" filter card"), "hidden"));
      await session.step(44, "And \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(45, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(46, "And user clears search of filter panel", () => clearField(page, el("search of filter panel")));
      await session.step(47, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
      await session.step(48, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(49, "When user loads the saved layout", () => loadLayout(page));
      await session.step(50, "Then \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(51, "And \"AGE\" filter card should be visible", () => shouldBe(page, el("\"AGE\" filter card"), "visible"));
      await session.step(52, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(53, "And the \"filtering of RACE\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of RACE", el("filter panel"), "false"));
      await session.step(54, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project reopens with its filter panel and its criterion", async () => {
      await session.step(58, "When user clicks on the \"category Caucasian of RACE\" area of filter panel", () => clickArea(page, "category Caucasian of RACE", el("filter panel")));
      await session.step(59, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(60, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(61, "When user saves the current view as project \"bdd filter panel round trip\"", () => saveAsProject(page, "bdd filter panel round trip"));
      await session.step(62, "And user closes all views", () => closeAllViews(page));
      await session.step(63, "And user opens the \"bdd filter panel round trip\" project", () => openProject(page, "bdd filter panel round trip"));
      await session.step(64, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(65, "And \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(66, "And the \"selected categories of RACE\" reading of filter panel should be \"Caucasian\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Caucasian"));
      await session.step(67, "And 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
