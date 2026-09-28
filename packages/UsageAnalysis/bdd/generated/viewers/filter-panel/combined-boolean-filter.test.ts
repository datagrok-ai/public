/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/combined-boolean-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, columnType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, filterIsExactlyCategory, filterPanelCount, filterPasses, filterPassesAll, openEmptyFilterPanel, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickPanelMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, loadLayout, noErrors, readingIs, readingReads, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Combined boolean filter card", () => {
  const session = feature(test, "features/viewers/filter-panel/combined-boolean-filter.feature", import.meta.url);
  test("Combined boolean filter card", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(20, "And user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
    await session.step(21, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
    await session.step(22, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
    await session.step(23, "And the \"categories of Flags\" reading of filter panel should be \"CONTROL\"", () => readingReads(page, "categories of Flags", el("filter panel"), "CONTROL"));
    await session.step(24, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(25, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
    await run.scenario("A boolean column added later joins the card", async () => {
      await session.step(28, "When user adds a calculated column \"SEX_bool\" with formula \"${SEX} == \\\"F\\\"\"", () => addCalculated(page, "SEX_bool", "${SEX} == \"F\""));
      await session.step(29, "Then \"SEX_bool\" column should have type \"bool\"", () => columnType(page, "SEX_bool", "bool"));
      await session.step(30, "And the \"categories of Flags\" reading of filter panel should be \"CONTROL, SEX_bool\"", () => readingReads(page, "categories of Flags", el("filter panel"), "CONTROL, SEX_bool"));
      await session.step(31, "And the \"cards\" reading of filter panel should be \"Flags\"", () => readingReads(page, "cards", el("filter panel"), "Flags"));
      await session.step(32, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(33, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a flag keeps its true rows", async () => {
      await session.step(37, "When user clicks on the \"category CONTROL of Flags\" area of filter panel", () => clickArea(page, "category CONTROL of Flags", el("filter panel")));
      await session.step(38, "Then 6 rows should pass the filter", () => filterPasses(page, 6));
      await session.step(39, "And the filter should pass exactly the rows where \"CONTROL\" is \"true\"", () => filterIsExactlyCategory(page, "CONTROL", "true"));
      await session.step(40, "And the \"summary of Flags\" reading of filter panel should be \"CONTROL: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "CONTROL: True"));
      await session.step(41, "And the \"filters\" reading of filter panel should be 1", () => readingIs(page, "filters", el("filter panel"), 1));
      await session.step(42, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on another flag replaces the one the card kept", async () => {
      await session.step(46, "When user clicks on the \"category SEX_bool of Flags\" area of filter panel", () => clickArea(page, "category SEX_bool of Flags", el("filter panel")));
      await session.step(47, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(48, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
      await session.step(49, "And the \"summary of Flags\" reading of filter panel should be \"SEX_bool: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "SEX_bool: True"));
      await session.step(50, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the card and its flag back", async () => {
      await session.step(54, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(55, "And user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(56, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(57, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(58, "When user loads the saved layout", () => loadLayout(page));
      await session.step(59, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
      await session.step(60, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(61, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(62, "And the \"summary of Flags\" reading of filter panel should be \"SEX_bool: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "SEX_bool: True"));
      await session.step(63, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the card releases the rows and Add Filter brings it back once", async () => {
      await session.step(67, "When user hovers over \"Flags\" filter card", () => hoverOver(page, el("\"Flags\" filter card")));
      await session.step(68, "And user clicks on close of \"Flags\" filter card", () => clickOn(page, el("close of \"Flags\" filter card")));
      await session.step(69, "Then \"Flags\" filter card should be absent", () => shouldBe(page, el("\"Flags\" filter card"), "absent"));
      await session.step(70, "And the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(71, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(72, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(73, "When user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
      await session.step(74, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
      await session.step(75, "And the \"cards\" reading of filter panel should be \"Flags\"", () => readingReads(page, "cards", el("filter panel"), "Flags"));
      await session.step(76, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(77, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(78, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(79, "And user removes \"SEX_bool\" column", () => removeColumn(page, "SEX_bool"));
      await session.step(80, "Then the table should have 11 columns", () => columnCount(page, 11));
      await session.step(81, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project reopens with the card, and the card can be removed after the reopen", async () => {
      await session.step(85, "When user adds a calculated column \"SEX_bool\" with formula \"${SEX} == \\\"F\\\"\"", () => addCalculated(page, "SEX_bool", "${SEX} == \"F\""));
      await session.step(86, "And user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
      await session.step(87, "And user clicks on the \"category SEX_bool of Flags\" area of filter panel", () => clickArea(page, "category SEX_bool of Flags", el("filter panel")));
      await session.step(88, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(89, "When user saves the current view as project \"bdd combined boolean round trip\"", () => saveAsProject(page, "bdd combined boolean round trip"));
      await session.step(90, "And user closes all views", () => closeAllViews(page));
      await session.step(91, "And user opens the \"bdd combined boolean round trip\" project", () => openProject(page, "bdd combined boolean round trip"));
      await session.step(92, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(93, "And \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
      await session.step(94, "And the \"summary of Flags\" reading of filter panel should be \"SEX_bool: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "SEX_bool: True"));
      await session.step(95, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(96, "When user hovers over \"Flags\" filter card", () => hoverOver(page, el("\"Flags\" filter card")));
      await session.step(97, "And user clicks on close of \"Flags\" filter card", () => clickOn(page, el("close of \"Flags\" filter card")));
      await session.step(98, "Then \"Flags\" filter card should be absent", () => shouldBe(page, el("\"Flags\" filter card"), "absent"));
      await session.step(99, "And the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(100, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A panel that works its cards out puts one card on the two boolean columns, every time", async () => {
      await session.step(104, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(105, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(106, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(107, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(108, "And the \"cards\" reading of filter panel should be \"Flags, SEX, RACE, SEVERITY, DIS_POP, AGE, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "cards", el("filter panel"), "Flags, SEX, RACE, SEVERITY, DIS_POP, AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(109, "And the \"categories of Flags\" reading of filter panel should be \"CONTROL, SEX_bool\"", () => readingReads(page, "categories of Flags", el("filter panel"), "CONTROL, SEX_bool"));
      await session.step(110, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(111, "When user picks \"Remove All\" from the viewer menu of filter panel", () => pickFromViewerMenu(page, "Remove All", el("filter panel")));
      await session.step(112, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(113, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(114, "And user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(115, "Then the \"cards\" reading of filter panel should be \"Flags, SEX, RACE, SEVERITY, DIS_POP, AGE, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "cards", el("filter panel"), "Flags, SEX, RACE, SEVERITY, DIS_POP, AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(116, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(117, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
