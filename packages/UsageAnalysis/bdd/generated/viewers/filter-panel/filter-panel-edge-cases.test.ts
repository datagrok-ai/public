/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters, core#1134]
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
import {clickOn, hoverOver, pressKey, selectIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPanelCount, filterPasses, filterPassesAll, openEmptyFilterPanel, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, panelHasNoCardOfType, pickCardIndicatorMenu, pickPanelMenu, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, closeContextMenu, noErrors, pickFromAreaContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel edge cases", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-edge-cases.feature", import.meta.url);
  test("An opening bracket typed into a categorical card's search logs no error", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:29", "@realizes:core#1134"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await knownFailure(async () => {
      await session.step(23, "When user adds a card for \"DIS_POP\" to the filter panel", () => addCardFor(page, "DIS_POP"));
      await session.step(24, "And user hovers over \"DIS_POP\" filter card", () => hoverOver(page, el("\"DIS_POP\" filter card")));
      await session.step(25, "And user clicks on search icon of \"DIS_POP\" filter card", () => clickOn(page, el("search icon of \"DIS_POP\" filter card")));
      await session.step(26, "And user types \"Ps\" into the search of the \"DIS_POP\" filter card", () => typeIntoCardSearch(page, "Ps", "DIS_POP"));
      await session.step(27, "Then 242 rows should pass the filter", () => filterPasses(page, 242));
      await session.step(28, "When user types \"Ps(\" into the search of the \"DIS_POP\" filter card", () => typeIntoCardSearch(page, "Ps(", "DIS_POP"));
      await session.step(29, "Then no errors should have been logged", () => noErrors(page));
    }, 29);
  });
  test("Removing the column of a filtering card removes the card and releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(33, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
    await session.step(34, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
    await session.step(35, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
    await session.step(36, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(37, "Then the filter panel should have no \"RACE\" filter card", () => panelHasNoCardOfType(page, "RACE"));
    await session.step(38, "And the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(39, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
  });
  test("The missing-values menu of a numeric card keeps only the blank rows, then all of them again", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(44, "When user adds a card for \"HEIGHT\" to the filter panel", () => addCardFor(page, "HEIGHT"));
    await session.step(45, "And user picks \"Missing values | Show only missing value\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Show only missing value", "HEIGHT"));
    await session.step(46, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(47, "Then 128 rows should pass the filter", () => filterPasses(page, 128));
    await session.step(48, "When user picks \"Missing values | Filter out missing values\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Filter out missing values", "HEIGHT"));
    await session.step(49, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(50, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
    await session.step(51, "When user picks \"Missing values | Keep missing values\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Keep missing values", "HEIGHT"));
    await session.step(52, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(53, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
  });
  test("The expression card's \"is null\" keeps the blank rows and \"is not null\" the filled ones", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(58, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
    await session.step(59, "And user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
    await session.step(60, "And user selects \"is null\" in Operation input in \"Expression\" filter card", () => selectIn(page, "is null", el("Operation input in \"Expression\" filter card")));
    await session.step(61, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
    await session.step(62, "Then 128 rows should pass the filter", () => filterPasses(page, 128));
    await session.step(63, "When user picks \"Remove Query\" from the context menu of the \"category ${HEIGHT} is null of Expression\" area of filter panel", () => pickFromAreaContextMenu(page, "Remove Query", "category ${HEIGHT} is null of Expression", el("filter panel")));
    await session.step(64, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(65, "When user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
    await session.step(66, "And user selects \"is not null\" in Operation input in \"Expression\" filter card", () => selectIn(page, "is not null", el("Operation input in \"Expression\" filter card")));
    await session.step(67, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
    await session.step(68, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
    await session.step(69, "And no errors should have been logged", () => noErrors(page));
  });
  test("A free-text rule that does not parse is refused on Enter and the rows stay", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(73, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
    await session.step(74, "And user hovers over \"Expression\" filter card", () => hoverOver(page, el("\"Expression\" filter card")));
    await session.step(75, "And user clicks on \"Switch to free-text mode\" icon in \"Expression\" filter card", () => clickOn(page, el("\"Switch to free-text mode\" icon in \"Expression\" filter card")));
    await session.step(76, "And user types \"${HEIGHT} < 150\" into the search of the \"Expression\" filter card", () => typeIntoCardSearch(page, "${HEIGHT} < 150", "Expression"));
    await session.step(77, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(78, "Then 16 rows should pass the filter", () => filterPasses(page, 16));
    await session.step(79, "When user types \"${HEIGHT} <\" into the search of the \"Expression\" filter card", () => typeIntoCardSearch(page, "${HEIGHT} <", "Expression"));
    await session.step(80, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(81, "Then the \"categories of Expression\" reading of filter panel should be \"${HEIGHT} < 150\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${HEIGHT} < 150"));
    await session.step(82, "And 16 rows should pass the filter", () => filterPasses(page, 16));
    await session.step(83, "And no errors should have been logged", () => noErrors(page));
  });
});
