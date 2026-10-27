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
import {filterPanelCount, filterPasses, filterPassesAll, filterPassesFewer, openEmptyFilterPanel, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, enterCardBound, panelHasNoCardOfType, pickCardIndicatorMenu, pickPanelMenu, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, closeContextMenu, noErrors, pickFromAreaContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel edge cases", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-edge-cases.feature", import.meta.url);
  test("An opening bracket typed into a categorical card's search logs no error", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:31", "@realizes:core#1134"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await knownFailure(async () => {
      await session.step(25, "When user adds a card for \"DIS_POP\" to the filter panel", () => addCardFor(page, "DIS_POP"));
      await session.step(26, "And user hovers over \"DIS_POP\" filter card", () => hoverOver(page, el("\"DIS_POP\" filter card")));
      await session.step(27, "And user clicks on search icon of \"DIS_POP\" filter card", () => clickOn(page, el("search icon of \"DIS_POP\" filter card")));
      await session.step(28, "And user types \"Ps\" into the search of the \"DIS_POP\" filter card", () => typeIntoCardSearch(page, "Ps", "DIS_POP"));
      await session.step(29, "Then 242 rows should pass the filter", () => filterPasses(page, 242));
      await session.step(30, "When user types \"Ps(\" into the search of the \"DIS_POP\" filter card", () => typeIntoCardSearch(page, "Ps(", "DIS_POP"));
      await session.step(31, "Then no errors should have been logged", () => noErrors(page));
    }, 31);
  });
  test("Removing the column of a filtering card removes the card and releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(35, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
    await session.step(36, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
    await session.step(37, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
    await session.step(38, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(39, "Then the filter panel should have no \"RACE\" filter card", () => panelHasNoCardOfType(page, "RACE"));
    await session.step(40, "And the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(41, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
  test("The missing-values menu of a numeric card keeps only the blank rows, then all of them again", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(46, "When user adds a card for \"HEIGHT\" to the filter panel", () => addCardFor(page, "HEIGHT"));
    await session.step(47, "And user picks \"Missing values | Show only missing value\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Show only missing value", "HEIGHT"));
    await session.step(48, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(49, "Then 128 rows should pass the filter", () => filterPasses(page, 128));
    await session.step(50, "When user picks \"Missing values | Filter out missing values\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Filter out missing values", "HEIGHT"));
    await session.step(51, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(52, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
    await session.step(53, "When user picks \"Missing values | Keep missing values\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Keep missing values", "HEIGHT"));
    await session.step(54, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(55, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
  });
  test("The expression card's \"is null\" keeps the blank rows and \"is not null\" the filled ones", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(60, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
    await session.step(61, "And user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
    await session.step(62, "And user selects \"is null\" in Operation input in \"Expression\" filter card", () => selectIn(page, "is null", el("Operation input in \"Expression\" filter card")));
    await session.step(63, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
    await session.step(64, "Then 128 rows should pass the filter", () => filterPasses(page, 128));
    await session.step(65, "When user picks \"Remove Query\" from the context menu of the \"category ${HEIGHT} is null of Expression\" area of filter panel", () => pickFromAreaContextMenu(page, "Remove Query", "category ${HEIGHT} is null of Expression", el("filter panel")));
    await session.step(66, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(67, "When user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
    await session.step(68, "And user selects \"is not null\" in Operation input in \"Expression\" filter card", () => selectIn(page, "is not null", el("Operation input in \"Expression\" filter card")));
    await session.step(69, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
    await session.step(70, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
    await session.step(71, "And no errors should have been logged", () => noErrors(page));
  });
  test("A free-text rule that does not parse is refused on Enter and the rows stay", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(75, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
    await session.step(76, "And user hovers over \"Expression\" filter card", () => hoverOver(page, el("\"Expression\" filter card")));
    await session.step(77, "And user clicks on \"Switch to free-text mode\" icon in \"Expression\" filter card", () => clickOn(page, el("\"Switch to free-text mode\" icon in \"Expression\" filter card")));
    await session.step(78, "And user types \"${HEIGHT} < 150\" into the search of the \"Expression\" filter card", () => typeIntoCardSearch(page, "${HEIGHT} < 150", "Expression"));
    await session.step(79, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(80, "Then 16 rows should pass the filter", () => filterPasses(page, 16));
    await session.step(81, "When user types \"${HEIGHT} <\" into the search of the \"Expression\" filter card", () => typeIntoCardSearch(page, "${HEIGHT} <", "Expression"));
    await session.step(82, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(83, "Then the \"categories of Expression\" reading of filter panel should be \"${HEIGHT} < 150\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${HEIGHT} < 150"));
    await session.step(84, "And 16 rows should pass the filter", () => filterPasses(page, 16));
    await session.step(85, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the column of a filtering numeric card releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(89, "When user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
    await session.step(90, "And user picks \"Min / max\" from the indicator menu of the \"AGE\" filter card", () => pickCardIndicatorMenu(page, "Min / max", "AGE"));
    await session.step(91, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(92, "And user enters \"40\" into the max field of the \"AGE\" filter card", () => enterCardBound(page, "40", "max", "AGE"));
    await session.step(93, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
    await session.step(94, "When user removes \"AGE\" column", () => removeColumn(page, "AGE"));
    await session.step(95, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(96, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(97, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the column of a filtering flag in the Combined Boolean card releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:107", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await knownFailure(async () => {
      await session.step(103, "When user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
      await session.step(104, "And user clicks on the \"category CONTROL of Flags\" area of filter panel", () => clickArea(page, "category CONTROL of Flags", el("filter panel")));
      await session.step(105, "Then 6 rows should pass the filter", () => filterPasses(page, 6));
      await session.step(106, "When user removes \"CONTROL\" column", () => removeColumn(page, "CONTROL"));
      await session.step(107, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    }, 107);
  });
  test("Removing the column of a filtering Expression card rule releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:120", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await knownFailure(async () => {
      await session.step(114, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(115, "And user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
      await session.step(116, "And user selects \"is null\" in Operation input in \"Expression\" filter card", () => selectIn(page, "is null", el("Operation input in \"Expression\" filter card")));
      await session.step(117, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(118, "Then 128 rows should pass the filter", () => filterPasses(page, 128));
      await session.step(119, "When user removes \"HEIGHT\" column", () => removeColumn(page, "HEIGHT"));
      await session.step(120, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
    }, 120);
  });
});
