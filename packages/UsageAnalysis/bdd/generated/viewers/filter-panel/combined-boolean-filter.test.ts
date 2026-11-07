/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/combined-boolean-filter.feature
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
import {clickOn, enterInto, hoverOver, replaceCode, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, columnType, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, filterIsExactlyCategory, filterPanelCount, filterPasses, filterPassesAll, openEmptyFilterPanel, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, closeCurrentView, openDataset, openProject, openTableOf, saveAsProject, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickPanelMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, closeContextMenu, loadLayout, menuListsEnabled, noErrors, readingIs, readingReads, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {openViewerMenu, pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Combined boolean filter card", () => {
  const session = feature(test, "features/viewers/filter-panel/combined-boolean-filter.feature", import.meta.url);
  test("Combined boolean filter card", {tag: ["@journey", "@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(24, "And user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
    await session.step(25, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
    await session.step(26, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
    await session.step(27, "And the \"categories of Flags\" reading of filter panel should be \"CONTROL\"", () => readingReads(page, "categories of Flags", el("filter panel"), "CONTROL"));
    await session.step(28, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(29, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
    await run.scenario("A boolean column added later joins the card", async () => {
      await session.step(32, "When user adds a calculated column \"SEX_bool\" with formula \"${SEX} == \\\"F\\\"\"", () => addCalculated(page, "SEX_bool", "${SEX} == \"F\""));
      await session.step(33, "Then \"SEX_bool\" column should have type \"bool\"", () => columnType(page, "SEX_bool", "bool"));
      await session.step(34, "And the \"categories of Flags\" reading of filter panel should be \"CONTROL, SEX_bool\"", () => readingReads(page, "categories of Flags", el("filter panel"), "CONTROL, SEX_bool"));
      await session.step(35, "And the \"cards\" reading of filter panel should be \"Flags\"", () => readingReads(page, "cards", el("filter panel"), "Flags"));
      await session.step(36, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(37, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a flag keeps its true rows", async () => {
      await session.step(41, "When user clicks on the \"category CONTROL of Flags\" area of filter panel", () => clickArea(page, "category CONTROL of Flags", el("filter panel")));
      await session.step(42, "Then 6 rows should pass the filter", () => filterPasses(page, 6));
      await session.step(43, "And the filter should pass exactly the rows where \"CONTROL\" is \"true\"", () => filterIsExactlyCategory(page, "CONTROL", "true"));
      await session.step(44, "And the \"summary of Flags\" reading of filter panel should be \"CONTROL: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "CONTROL: True"));
      await session.step(45, "And the \"filters\" reading of filter panel should be 1", () => readingIs(page, "filters", el("filter panel"), 1));
      await session.step(46, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on another flag replaces the one the card kept", async () => {
      await session.step(50, "When user clicks on the \"category SEX_bool of Flags\" area of filter panel", () => clickArea(page, "category SEX_bool of Flags", el("filter panel")));
      await session.step(51, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(52, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
      await session.step(53, "And the \"summary of Flags\" reading of filter panel should be \"SEX_bool: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "SEX_bool: True"));
      await session.step(54, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the card and its flag back", async () => {
      await session.step(58, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(59, "And user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(60, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(61, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(62, "When user loads the saved layout", () => loadLayout(page));
      await session.step(63, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
      await session.step(64, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(65, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(66, "And the \"summary of Flags\" reading of filter panel should be \"SEX_bool: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "SEX_bool: True"));
      await session.step(67, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Two flags ticked at once combine with OR, and the card's switch makes them AND", async () => {
      await session.step(71, "Then the \"summary of Flags\" reading of filter panel should be \"SEX_bool: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "SEX_bool: True"));
      await session.step(72, "When user clicks on the \"true checkbox CONTROL of Flags\" area of filter panel", () => clickArea(page, "true checkbox CONTROL of Flags", el("filter panel")));
      await session.step(73, "Then the \"summary of Flags\" reading of filter panel should be \"CONTROL: True, SEX_bool: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "CONTROL: True, SEX_bool: True"));
      await session.step(74, "And mode of \"Flags\" filter card should have text \"OR\"", () => shouldHaveText(page, el("mode of \"Flags\" filter card"), "OR"));
      await session.step(75, "And 554 rows should pass the filter", () => filterPasses(page, 554));
      await session.step(76, "When user clicks on mode of \"Flags\" filter card", () => clickOn(page, el("mode of \"Flags\" filter card")));
      await session.step(77, "Then mode of \"Flags\" filter card should have text \"AND\"", () => shouldHaveText(page, el("mode of \"Flags\" filter card"), "AND"));
      await session.step(78, "And 5 rows should pass the filter", () => filterPasses(page, 5));
      await session.step(79, "When user clicks on mode of \"Flags\" filter card", () => clickOn(page, el("mode of \"Flags\" filter card")));
      await session.step(80, "Then mode of \"Flags\" filter card should have text \"OR\"", () => shouldHaveText(page, el("mode of \"Flags\" filter card"), "OR"));
      await session.step(81, "And 554 rows should pass the filter", () => filterPasses(page, 554));
      await session.step(82, "When user clicks on the \"true checkbox CONTROL of Flags\" area of filter panel", () => clickArea(page, "true checkbox CONTROL of Flags", el("filter panel")));
      await session.step(83, "Then the \"summary of Flags\" reading of filter panel should be \"SEX_bool: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "SEX_bool: True"));
      await session.step(84, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the card releases the rows and Add Filter brings it back once", async () => {
      await session.step(88, "When user hovers over \"Flags\" filter card", () => hoverOver(page, el("\"Flags\" filter card")));
      await session.step(89, "And user clicks on close of \"Flags\" filter card", () => clickOn(page, el("close of \"Flags\" filter card")));
      await session.step(90, "Then \"Flags\" filter card should be absent", () => shouldBe(page, el("\"Flags\" filter card"), "absent"));
      await session.step(91, "And the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(92, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(93, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(94, "When user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
      await session.step(95, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
      await session.step(96, "And the \"cards\" reading of filter panel should be \"Flags\"", () => readingReads(page, "cards", el("filter panel"), "Flags"));
      await session.step(97, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(98, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(99, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(100, "And user removes \"SEX_bool\" column", () => removeColumn(page, "SEX_bool"));
      await session.step(101, "Then the table should have 11 columns", () => columnCount(page, 11));
      await session.step(102, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project reopens with the card, and the card can be removed after the reopen", async () => {
      await session.step(106, "When user adds a calculated column \"SEX_bool\" with formula \"${SEX} == \\\"F\\\"\"", () => addCalculated(page, "SEX_bool", "${SEX} == \"F\""));
      await session.step(107, "And user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
      await session.step(108, "And user clicks on the \"category SEX_bool of Flags\" area of filter panel", () => clickArea(page, "category SEX_bool of Flags", el("filter panel")));
      await session.step(109, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(110, "When user saves the current view as project \"bdd combined boolean round trip\"", () => saveAsProject(page, "bdd combined boolean round trip"));
      await session.step(111, "And user closes all views", () => closeAllViews(page));
      await session.step(112, "And user opens the \"bdd combined boolean round trip\" project", () => openProject(page, "bdd combined boolean round trip"));
      await session.step(113, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(114, "And \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
      await session.step(115, "And the \"summary of Flags\" reading of filter panel should be \"SEX_bool: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "SEX_bool: True"));
      await session.step(116, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(117, "When user hovers over \"Flags\" filter card", () => hoverOver(page, el("\"Flags\" filter card")));
      await session.step(118, "And user clicks on close of \"Flags\" filter card", () => clickOn(page, el("close of \"Flags\" filter card")));
      await session.step(119, "Then \"Flags\" filter card should be absent", () => shouldBe(page, el("\"Flags\" filter card"), "absent"));
      await session.step(120, "And the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(121, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A panel that works its cards out puts one card on the two boolean columns, every time", async () => {
      await session.step(125, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(126, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(127, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(128, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(129, "And the \"cards\" reading of filter panel should be \"Flags, SEX, RACE, SEVERITY, DIS_POP, AGE, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "cards", el("filter panel"), "Flags, SEX, RACE, SEVERITY, DIS_POP, AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(130, "And the \"categories of Flags\" reading of filter panel should be \"CONTROL, SEX_bool\"", () => readingReads(page, "categories of Flags", el("filter panel"), "CONTROL, SEX_bool"));
      await session.step(131, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(132, "When user picks \"Remove All\" from the viewer menu of filter panel", () => pickFromViewerMenu(page, "Remove All", el("filter panel")));
      await session.step(133, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(134, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(135, "And user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(136, "Then the \"cards\" reading of filter panel should be \"Flags, SEX, RACE, SEVERITY, DIS_POP, AGE, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "cards", el("filter panel"), "Flags, SEX, RACE, SEVERITY, DIS_POP, AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(137, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(138, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Combined Boolean card can be added to a table with no boolean column, and a boolean column added later joins it", async () => {
      await session.step(142, "Given user opens a table \"no-flags\" with:", () => openTableOf(page, "no-flags", [["Name","Value"],["a","1"],["b","2"]]), [["Name","Value"],["a","1"],["b","2"]]);
      await session.step(146, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(147, "And user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
      await session.step(148, "Then the open menu should list \"Add Filter > Combined Boolean\" as enabled", () => menuListsEnabled(page, "Add Filter > Combined Boolean"));
      await session.step(149, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(150, "And user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
      await session.step(151, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
      await session.step(152, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(153, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(154, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(155, "When user enters \"BIG\" into Name input in \"Add New Column\" dialog", () => enterInto(page, "BIG", el("Name input in \"Add New Column\" dialog")));
      await session.step(156, "And user replaces the code of code editor in \"Add New Column\" dialog with \"${Value} > 1\"", () => replaceCode(page, el("code editor in \"Add New Column\" dialog"), "${Value} > 1"));
      await session.step(157, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(158, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(159, "And the table should have a column \"BIG\"", () => hasColumn(page, "BIG"));
      await session.step(160, "When user clicks on the \"category BIG of Flags\" area of filter panel", () => clickArea(page, "category BIG of Flags", el("filter panel")));
      await session.step(161, "Then 1 row should pass the filter", () => filterPasses(page, 1));
      await session.step(162, "When user closes the current view", () => closeCurrentView(page));
      await session.step(163, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(164, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(165, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
