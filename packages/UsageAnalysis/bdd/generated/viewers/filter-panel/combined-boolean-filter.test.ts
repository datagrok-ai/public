/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/combined-boolean-filter.feature
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
import {pickPanelMenu} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, columnType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, filterIsExactlyCategory, filterPanelCount, filterPasses, filterPassesAll, openEmptyFilterPanel, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, loadLayout, noErrors, readingIs, readingReads, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Combined boolean filter card", () => {
  const session = feature(test, "features/viewers/filter-panel/combined-boolean-filter.feature", import.meta.url);
  test("Combined boolean filter card", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(13, "And user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
    await session.step(14, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
    await session.step(15, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
    await session.step(16, "And the \"categories of Flags\" reading of filter panel should be \"CONTROL\"", () => readingReads(page, "categories of Flags", el("filter panel"), "CONTROL"));
    await session.step(17, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(18, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
    await run.scenario("A boolean column added later joins the card", async () => {
      await session.step(21, "When user adds a calculated column \"SEX_bool\" with formula \"${SEX} == \\\"F\\\"\"", () => addCalculated(page, "SEX_bool", "${SEX} == \"F\""));
      await session.step(22, "Then \"SEX_bool\" column should have type \"bool\"", () => columnType(page, "SEX_bool", "bool"));
      await session.step(23, "And the \"categories of Flags\" reading of filter panel should be \"CONTROL, SEX_bool\"", () => readingReads(page, "categories of Flags", el("filter panel"), "CONTROL, SEX_bool"));
      await session.step(24, "And the \"cards\" reading of filter panel should be \"Flags\"", () => readingReads(page, "cards", el("filter panel"), "Flags"));
      await session.step(25, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(26, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(27, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a flag keeps its true rows", async () => {
      await session.step(30, "When user clicks on the \"category CONTROL of Flags\" area of filter panel", () => clickArea(page, "category CONTROL of Flags", el("filter panel")));
      await session.step(31, "Then 6 rows should pass the filter", () => filterPasses(page, 6));
      await session.step(32, "And the filter should pass exactly the rows where \"CONTROL\" is \"true\"", () => filterIsExactlyCategory(page, "CONTROL", "true"));
      await session.step(33, "And the \"summary of Flags\" reading of filter panel should be \"CONTROL: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "CONTROL: True"));
      await session.step(34, "And the \"filters\" reading of filter panel should be 1", () => readingIs(page, "filters", el("filter panel"), 1));
      await session.step(35, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on another flag replaces the one the card kept", async () => {
      await session.step(39, "When user clicks on the \"category SEX_bool of Flags\" area of filter panel", () => clickArea(page, "category SEX_bool of Flags", el("filter panel")));
      await session.step(40, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(41, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
      await session.step(42, "And the \"summary of Flags\" reading of filter panel should be \"SEX_bool: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "SEX_bool: True"));
      await session.step(43, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the card and its flag back", async () => {
      await session.step(47, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(48, "And user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(49, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(50, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(51, "When user loads the saved layout", () => loadLayout(page));
      await session.step(52, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
      await session.step(53, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(54, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(55, "And the \"summary of Flags\" reading of filter panel should be \"SEX_bool: True\"", () => readingReads(page, "summary of Flags", el("filter panel"), "SEX_bool: True"));
      await session.step(56, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the card releases the rows and Add Filter brings it back once", async () => {
      await session.step(60, "When user hovers over \"Flags\" filter card", () => hoverOver(page, el("\"Flags\" filter card")));
      await session.step(61, "And user clicks on close of \"Flags\" filter card", () => clickOn(page, el("close of \"Flags\" filter card")));
      await session.step(62, "Then \"Flags\" filter card should be absent", () => shouldBe(page, el("\"Flags\" filter card"), "absent"));
      await session.step(63, "And the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(64, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(65, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(66, "When user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
      await session.step(67, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
      await session.step(68, "And the \"cards\" reading of filter panel should be \"Flags\"", () => readingReads(page, "cards", el("filter panel"), "Flags"));
      await session.step(69, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(70, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(71, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(72, "And user removes \"SEX_bool\" column", () => removeColumn(page, "SEX_bool"));
      await session.step(73, "Then the table should have 11 columns", () => columnCount(page, 11));
      await session.step(74, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
