/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/hierarchical-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters, core#1223]
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
import {clickHierarchicalCheckbox, clickHierarchicalRow, expandHierarchicalRow, hierarchicalHidesRow, hierarchicalListsRow, hierarchicalRowCounts, hierarchicalRowState} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {configureHierarchical, filterIsExactlyCategory, filterPanelCount, filterPasses, filterPassesAll, noneOfFiltered, openEmptyFilterPanel, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clearCardSearch, pickPanelMenu, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {dragAreaToArea, loadLayout, noErrors, readingIs, readingReads, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hierarchical filter card", () => {
  const session = feature(test, "features/viewers/filter-panel/hierarchical-filter.feature", import.meta.url);
  test("Hierarchical filter card", {tag: ["@journey", "@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:165", "@realizes:core#1223"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(23, "And user picks \"Add Filter | Hierarchical\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Hierarchical"));
    await session.step(24, "And user configures the hierarchical filter with columns \"SEX, RACE\"", () => configureHierarchical(page, "SEX, RACE"));
    await session.step(25, "Then \"SEX / RACE\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE\" filter card"), "visible"));
    await session.step(26, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
    await session.step(27, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(28, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
    await run.scenario("A click on a branch keeps only its rows", async () => {
      await session.step(31, "Then the hierarchical filter card should list the \"F\" row", () => hierarchicalListsRow(page, "F"));
      await session.step(32, "And the \"F\" row of the hierarchical filter card should count 553 rows", () => hierarchicalRowCounts(page, "F", 553));
      await session.step(33, "When user clicks on the \"F\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "F"));
      await session.step(34, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(35, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
      await session.step(36, "And the \"filters\" reading of filter panel should be 1", () => readingIs(page, "filters", el("filter panel"), 1));
      await session.step(37, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(38, "And the \"M\" row of the hierarchical filter card should count 0 rows", () => hierarchicalRowCounts(page, "M", 0));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A branch opens to the next level and a click there narrows to the leaf", async () => {
      await session.step(42, "When user expands the \"F\" row of the hierarchical filter card", () => expandHierarchicalRow(page, "F"));
      await session.step(43, "Then the hierarchical filter card should list the \"F / Caucasian\" row", () => hierarchicalListsRow(page, "F / Caucasian"));
      await session.step(44, "And the hierarchical filter card should list the \"F / Other\" row", () => hierarchicalListsRow(page, "F / Other"));
      await session.step(45, "When user clicks on the \"F / Caucasian\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "F / Caucasian"));
      await session.step(46, "Then 480 rows should pass the filter", () => filterPasses(page, 480));
      await session.step(47, "And the \"F / Caucasian\" row of the hierarchical filter card should count 480 rows", () => hierarchicalRowCounts(page, "F / Caucasian", 480));
      await session.step(48, "And the \"F / Other\" row of the hierarchical filter card should count 0 rows", () => hierarchicalRowCounts(page, "F / Other", 0));
      await session.step(49, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A branch with one child switched off is partially checked", async () => {
      await session.step(53, "When user clicks on the \"F\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "F"));
      await session.step(54, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(55, "And the \"F\" row of the hierarchical filter card should be checked", () => hierarchicalRowState(page, "F", "checked"));
      await session.step(56, "And the \"F / Caucasian\" row of the hierarchical filter card should be checked", () => hierarchicalRowState(page, "F / Caucasian", "checked"));
      await session.step(57, "And the \"M\" row of the hierarchical filter card should be unchecked", () => hierarchicalRowState(page, "M", "unchecked"));
      await session.step(58, "When user clicks on the checkbox of the \"F / Caucasian\" row of the hierarchical filter card", () => clickHierarchicalCheckbox(page, "F / Caucasian"));
      await session.step(59, "Then 73 rows should pass the filter", () => filterPasses(page, 73));
      await session.step(60, "And the \"F / Caucasian\" row of the hierarchical filter card should be unchecked", () => hierarchicalRowState(page, "F / Caucasian", "unchecked"));
      await session.step(61, "And the \"F / Other\" row of the hierarchical filter card should be checked", () => hierarchicalRowState(page, "F / Other", "checked"));
      await session.step(62, "And the \"F\" row of the hierarchical filter card should be partially checked", () => hierarchicalRowState(page, "F", "partially"));
      await session.step(63, "And the \"M\" row of the hierarchical filter card should be unchecked", () => hierarchicalRowState(page, "M", "unchecked"));
      await session.step(64, "And no rows where \"RACE\" is \"Caucasian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Caucasian"));
      await session.step(65, "And no rows where \"SEX\" is \"M\" should pass the filter", () => noneOfFiltered(page, "SEX", "M"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A third level picked in the card's column list narrows inside the leaf of the second", async () => {
      await session.step(69, "When user hovers over \"SEX / RACE\" filter card", () => hoverOver(page, el("\"SEX / RACE\" filter card")));
      await session.step(70, "And user clicks on \"Select columns\" icon in \"SEX / RACE\" filter card", () => clickOn(page, el("\"Select columns\" icon in \"SEX / RACE\" filter card")));
      await session.step(71, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(72, "And the \"text of cell 1 of __name\" reading of Grid viewer in \"Select columns...\" dialog should be \"SEX\"", () => readingReads(page, "text of cell 1 of __name", el("Grid viewer in \"Select columns...\" dialog"), "SEX"));
      await session.step(73, "And the \"text of cell 2 of __name\" reading of Grid viewer in \"Select columns...\" dialog should be \"RACE\"", () => readingReads(page, "text of cell 2 of __name", el("Grid viewer in \"Select columns...\" dialog"), "RACE"));
      await session.step(74, "When user toggles the \"SEVERITY\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEVERITY", el("\"Select columns...\" dialog")));
      await session.step(75, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(76, "Then \"SEX / RACE / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE / SEVERITY\" filter card"), "visible"));
      await session.step(77, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(78, "When user expands the \"F\" row of the hierarchical filter card", () => expandHierarchicalRow(page, "F"));
      await session.step(79, "And user expands the \"F / Caucasian\" row of the hierarchical filter card", () => expandHierarchicalRow(page, "F / Caucasian"));
      await session.step(80, "Then the hierarchical filter card should list the \"F / Caucasian / None\" row", () => hierarchicalListsRow(page, "F / Caucasian / None"));
      await session.step(81, "When user clicks on the \"F / Caucasian / None\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "F / Caucasian / None"));
      await session.step(82, "Then 285 rows should pass the filter", () => filterPasses(page, 285));
      await session.step(83, "And the \"F / Caucasian / None\" row of the hierarchical filter card should count 285 rows", () => hierarchicalRowCounts(page, "F / Caucasian / None", 285));
      await session.step(84, "And the \"F / Caucasian / Low\" row of the hierarchical filter card should count 0 rows", () => hierarchicalRowCounts(page, "F / Caucasian / Low", 0));
      await session.step(85, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A grandchild switched off leaves its parent and grandparent partially checked", async () => {
      await session.step(89, "When user clicks on the \"F\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "F"));
      await session.step(90, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(91, "And the \"F / Caucasian\" row of the hierarchical filter card should be checked", () => hierarchicalRowState(page, "F / Caucasian", "checked"));
      await session.step(92, "And the \"F / Caucasian / None\" row of the hierarchical filter card should be checked", () => hierarchicalRowState(page, "F / Caucasian / None", "checked"));
      await session.step(93, "When user clicks on the checkbox of the \"F / Caucasian / None\" row of the hierarchical filter card", () => clickHierarchicalCheckbox(page, "F / Caucasian / None"));
      await session.step(94, "Then 268 rows should pass the filter", () => filterPasses(page, 268));
      await session.step(95, "And the \"F / Caucasian / None\" row of the hierarchical filter card should be unchecked", () => hierarchicalRowState(page, "F / Caucasian / None", "unchecked"));
      await session.step(96, "And the \"F / Caucasian / Low\" row of the hierarchical filter card should be checked", () => hierarchicalRowState(page, "F / Caucasian / Low", "checked"));
      await session.step(97, "And the \"F / Caucasian / High\" row of the hierarchical filter card should be checked", () => hierarchicalRowState(page, "F / Caucasian / High", "checked"));
      await session.step(98, "And the \"F / Caucasian / Medium\" row of the hierarchical filter card should be checked", () => hierarchicalRowState(page, "F / Caucasian / Medium", "checked"));
      await session.step(99, "And the hierarchical filter card should not list the \"F / Caucasian / Critical\" row", () => hierarchicalHidesRow(page, "F / Caucasian / Critical"));
      await session.step(100, "And the \"F / Caucasian\" row of the hierarchical filter card should be partially checked", () => hierarchicalRowState(page, "F / Caucasian", "partially"));
      await session.step(101, "And the \"F\" row of the hierarchical filter card should be partially checked", () => hierarchicalRowState(page, "F", "partially"));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Reordering the levels in the column list rebuilds the roots and clears the criterion", async () => {
      await session.step(105, "Then 268 rows should pass the filter", () => filterPasses(page, 268));
      await session.step(106, "When user hovers over \"SEX / RACE / SEVERITY\" filter card", () => hoverOver(page, el("\"SEX / RACE / SEVERITY\" filter card")));
      await session.step(107, "And user clicks on \"Select columns\" icon in \"SEX / RACE / SEVERITY\" filter card", () => clickOn(page, el("\"Select columns\" icon in \"SEX / RACE / SEVERITY\" filter card")));
      await session.step(108, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(109, "When user toggles the \"SEVERITY\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEVERITY", el("\"Select columns...\" dialog")));
      await session.step(110, "And user drags the \"cell 2 of __name\" area of Grid viewer in \"Select columns...\" dialog to the \"cell 1 of __name\" area", () => dragAreaToArea(page, "cell 2 of __name", el("Grid viewer in \"Select columns...\" dialog"), "cell 1 of __name"));
      await session.step(111, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(112, "Then \"RACE / SEX\" filter card should be visible", () => shouldBe(page, el("\"RACE / SEX\" filter card"), "visible"));
      await session.step(113, "And the hierarchical filter card should list the \"Caucasian\" row", () => hierarchicalListsRow(page, "Caucasian"));
      await session.step(114, "And the hierarchical filter card should not list the \"F\" row", () => hierarchicalHidesRow(page, "F"));
      await session.step(115, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(116, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(117, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The card's search hides rows and leaves the table alone", async () => {
      await session.step(120, "When user clicks on the \"Asian\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "Asian"));
      await session.step(121, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(122, "When user hovers over \"RACE / SEX\" filter card", () => hoverOver(page, el("\"RACE / SEX\" filter card")));
      await session.step(123, "And user clicks on search icon of \"RACE / SEX\" filter card", () => clickOn(page, el("search icon of \"RACE / SEX\" filter card")));
      await session.step(124, "And user types \"Cau\" into the search of the \"RACE / SEX\" filter card", () => typeIntoCardSearch(page, "Cau", "RACE / SEX"));
      await session.step(125, "Then the hierarchical filter card should not list the \"Asian\" row", () => hierarchicalHidesRow(page, "Asian"));
      await session.step(126, "And the hierarchical filter card should list the \"Caucasian\" row", () => hierarchicalListsRow(page, "Caucasian"));
      await session.step(127, "And 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(128, "And no rows where \"RACE\" is \"Caucasian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Caucasian"));
      await session.step(129, "When user clears the search of the \"RACE / SEX\" filter card", () => clearCardSearch(page, "RACE / SEX"));
      await session.step(130, "Then the hierarchical filter card should list the \"Asian\" row", () => hierarchicalListsRow(page, "Asian"));
      await session.step(131, "And 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the levels and the criterion back", async () => {
      await session.step(135, "When user configures the hierarchical filter with columns \"SEX, RACE\"", () => configureHierarchical(page, "SEX, RACE"));
      await session.step(136, "And user clicks on the \"F\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "F"));
      await session.step(137, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(138, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(139, "And user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(140, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(141, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(142, "When user loads the saved layout", () => loadLayout(page));
      await session.step(143, "Then \"SEX / RACE\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE\" filter card"), "visible"));
      await session.step(144, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(145, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(146, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(147, "When user configures the hierarchical filter with columns \"RACE, SEX\"", () => configureHierarchical(page, "RACE, SEX"));
      await session.step(148, "And user clicks on the \"Asian\" row of the hierarchical filter card", () => clickHierarchicalRow(page, "Asian"));
      await session.step(149, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(150, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(151, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(152, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(153, "When user loads the saved layout", () => loadLayout(page));
      await session.step(154, "Then \"SEX / RACE\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE\" filter card"), "visible"));
      await session.step(155, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(156, "And the \"F\" row of the hierarchical filter card should be checked", () => hierarchicalRowState(page, "F", "checked"));
      await session.step(157, "And the \"M\" row of the hierarchical filter card should be unchecked", () => hierarchicalRowState(page, "M", "unchecked"));
      await session.step(158, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the column of a hierarchical card's middle level drops that level", async () => {
      await session.step(162, "When user configures the hierarchical filter with columns \"SEX, RACE, SEVERITY\"", () => configureHierarchical(page, "SEX, RACE, SEVERITY"));
      await session.step(163, "Then \"SEX / RACE / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE / SEVERITY\" filter card"), "visible"));
      await session.step(164, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(165, "Then \"SEX / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / SEVERITY\" filter card"), "visible"));
    }, {knownFailure: true, line: 165});
    run.finish();
  });
});
