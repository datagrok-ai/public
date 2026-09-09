/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/hierarchical-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clearCardSearch, clickHierarchicalRow, expandHierarchicalRow, hierarchicalHidesRow, hierarchicalListsRow, hierarchicalRowCounts, pickPanelMenu, typeIntoCardSearch} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {configureHierarchical, filterIsExactlyCategory, filterPanelCount, filterPasses, filterPassesAll, noneOfFiltered, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {loadLayout, noErrors, readingIs, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hierarchical filter card", () => {
  const session = feature(test, "features/viewers/filter-panel/hierarchical-filter.feature", import.meta.url);
  test("Hierarchical filter card", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(14, "And user picks \"Add Filter | Hierarchical\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Hierarchical"));
    await session.step(15, "And user configures the hierarchical filter with columns \"SEX, RACE\"", () => configureHierarchical(page, "SEX, RACE"));
    await session.step(16, "Then \"SEX / RACE\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE\" filter card"), "visible"));
    await session.step(17, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
    await session.step(18, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(19, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
    await run.scenario("A click on a branch keeps only its rows", async () => {
      await session.step(22, "Then the hierarchical filter card should list the \"F\" row", () => hierarchicalListsRow(page, "F"));
      await session.step(23, "And the \"F\" row of the hierarchical filter card should count 553 rows", () => hierarchicalRowCounts(page, "F", 553));
      await session.step(24, "When user clicks on the \"F\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "F"));
      await session.step(25, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(26, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
      await session.step(27, "And the \"filters\" reading of filter panel should be 1", () => readingIs(page, "filters", el("filter panel"), 1));
      await session.step(28, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(29, "And the \"M\" row of the hierarchical filter card should count 0 rows", () => hierarchicalRowCounts(page, "M", 0));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A branch opens to the next level and a click there narrows to the leaf", async () => {
      await session.step(33, "When user expands the \"F\" row of the hierarchical filter card", () => expandHierarchicalRow(page, "F"));
      await session.step(34, "Then the hierarchical filter card should list the \"F / Caucasian\" row", () => hierarchicalListsRow(page, "F / Caucasian"));
      await session.step(35, "And the hierarchical filter card should list the \"F / Other\" row", () => hierarchicalListsRow(page, "F / Other"));
      await session.step(36, "When user clicks on the \"F / Caucasian\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "F / Caucasian"));
      await session.step(37, "Then 480 rows should pass the filter", () => filterPasses(page, 480));
      await session.step(38, "And the \"F / Caucasian\" row of the hierarchical filter card should count 480 rows", () => hierarchicalRowCounts(page, "F / Caucasian", 480));
      await session.step(39, "And the \"F / Other\" row of the hierarchical filter card should count 0 rows", () => hierarchicalRowCounts(page, "F / Other", 0));
      await session.step(40, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A third level narrows inside the leaf of the second", async () => {
      await session.step(44, "When user configures the hierarchical filter with columns \"SEX, RACE, SEVERITY\"", () => configureHierarchical(page, "SEX, RACE, SEVERITY"));
      await session.step(45, "Then \"SEX / RACE / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE / SEVERITY\" filter card"), "visible"));
      await session.step(46, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(47, "When user expands the \"F\" row of the hierarchical filter card", () => expandHierarchicalRow(page, "F"));
      await session.step(48, "And user expands the \"F / Caucasian\" row of the hierarchical filter card", () => expandHierarchicalRow(page, "F / Caucasian"));
      await session.step(49, "Then the hierarchical filter card should list the \"F / Caucasian / None\" row", () => hierarchicalListsRow(page, "F / Caucasian / None"));
      await session.step(50, "When user clicks on the \"F / Caucasian / None\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "F / Caucasian / None"));
      await session.step(51, "Then 285 rows should pass the filter", () => filterPasses(page, 285));
      await session.step(52, "And the \"F / Caucasian / None\" row of the hierarchical filter card should count 285 rows", () => hierarchicalRowCounts(page, "F / Caucasian / None", 285));
      await session.step(53, "And the \"F / Caucasian / Low\" row of the hierarchical filter card should count 0 rows", () => hierarchicalRowCounts(page, "F / Caucasian / Low", 0));
      await session.step(54, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Reordering the levels rebuilds the roots and clears the criterion", async () => {
      await session.step(58, "When user configures the hierarchical filter with columns \"RACE, SEX\"", () => configureHierarchical(page, "RACE, SEX"));
      await session.step(59, "Then \"RACE / SEX\" filter card should be visible", () => shouldBe(page, el("\"RACE / SEX\" filter card"), "visible"));
      await session.step(60, "And the hierarchical filter card should list the \"Caucasian\" row", () => hierarchicalListsRow(page, "Caucasian"));
      await session.step(61, "And the hierarchical filter card should not list the \"F\" row", () => hierarchicalHidesRow(page, "F"));
      await session.step(62, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(63, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The card's search hides rows and leaves the table alone", async () => {
      await session.step(67, "When user clicks on the \"Asian\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "Asian"));
      await session.step(68, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(69, "When user hovers over \"RACE / SEX\" filter card", () => hoverOver(page, el("\"RACE / SEX\" filter card")));
      await session.step(70, "And user clicks on search icon of \"RACE / SEX\" filter card", () => clickOn(page, el("search icon of \"RACE / SEX\" filter card")));
      await session.step(71, "And user types \"Cau\" into the search of the \"RACE / SEX\" filter card", () => typeIntoCardSearch(page, "Cau", "RACE / SEX"));
      await session.step(72, "Then the hierarchical filter card should not list the \"Asian\" row", () => hierarchicalHidesRow(page, "Asian"));
      await session.step(73, "And the hierarchical filter card should list the \"Caucasian\" row", () => hierarchicalListsRow(page, "Caucasian"));
      await session.step(74, "And 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(75, "And no rows where \"RACE\" is \"Caucasian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Caucasian"));
      await session.step(76, "When user clears the search of the \"RACE / SEX\" filter card", () => clearCardSearch(page, "RACE / SEX"));
      await session.step(77, "Then the hierarchical filter card should list the \"Asian\" row", () => hierarchicalListsRow(page, "Asian"));
      await session.step(78, "And 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the levels and the criterion back", async () => {
      await session.step(82, "When user configures the hierarchical filter with columns \"SEX, RACE\"", () => configureHierarchical(page, "SEX, RACE"));
      await session.step(83, "And user clicks on the \"F\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "F"));
      await session.step(84, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(85, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(86, "And user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(87, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(88, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(89, "When user loads the saved layout", () => loadLayout(page));
      await session.step(90, "Then \"SEX / RACE\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE\" filter card"), "visible"));
      await session.step(91, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(92, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(93, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
