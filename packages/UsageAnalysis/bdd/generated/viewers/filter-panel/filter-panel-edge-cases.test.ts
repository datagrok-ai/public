/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters, core#1174, core#1134, core#1223, core#1224]
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
import {expandHierarchicalRow, hierarchicalListsRow, linkTables} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, pressKey, selectIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, configureHierarchical, filterIsExactlyCategory, filterPanelCount, filterPasses, filterPassesAll, filterPassesFewer, openEmptyFilterPanel, removeColumn, tableFilterCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, enterCardBound, forgetSavedState, panelHasNoCardOfType, pickCardIndicatorMenu, pickPanelMenu, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, closeContextMenu, loadLayout, menuListsEnabled, noErrors, pickFromAreaContextMenu, readingReads, saveLayout} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {openViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel edge cases", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-edge-cases.feature", import.meta.url);
  test("Removing the column of a filtering card removes the card and releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(23, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
    await session.step(24, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
    await session.step(25, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
    await session.step(26, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(27, "Then the filter panel should have no \"RACE\" filter card", () => panelHasNoCardOfType(page, "RACE"));
    await session.step(28, "And the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(29, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
  });
  test("The missing-values menu of a numeric card keeps only the blank rows, then all of them again", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(34, "When user adds a card for \"HEIGHT\" to the filter panel", () => addCardFor(page, "HEIGHT"));
    await session.step(35, "And user picks \"Missing values | Show only missing value\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Show only missing value", "HEIGHT"));
    await session.step(36, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(37, "Then 128 rows should pass the filter", () => filterPasses(page, 128));
    await session.step(38, "When user picks \"Missing values | Filter out missing values\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Filter out missing values", "HEIGHT"));
    await session.step(39, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(40, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
    await session.step(41, "When user picks \"Missing values | Keep missing values\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Keep missing values", "HEIGHT"));
    await session.step(42, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(43, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
  test("The expression card's \"is null\" keeps the blank rows and \"is not null\" the filled ones", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(48, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
    await session.step(49, "And user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
    await session.step(50, "And user selects \"is null\" in Operation input in \"Expression\" filter card", () => selectIn(page, "is null", el("Operation input in \"Expression\" filter card")));
    await session.step(51, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
    await session.step(52, "Then 128 rows should pass the filter", () => filterPasses(page, 128));
    await session.step(53, "When user picks \"Remove Query\" from the context menu of the \"category ${HEIGHT} is null of Expression\" area of filter panel", () => pickFromAreaContextMenu(page, "Remove Query", "category ${HEIGHT} is null of Expression", el("filter panel")));
    await session.step(54, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(55, "When user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
    await session.step(56, "And user selects \"is not null\" in Operation input in \"Expression\" filter card", () => selectIn(page, "is not null", el("Operation input in \"Expression\" filter card")));
    await session.step(57, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
    await session.step(58, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
  });
  test("A free-text rule that does not parse is refused on Enter and the rows stay", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(63, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
    await session.step(64, "And user hovers over \"Expression\" filter card", () => hoverOver(page, el("\"Expression\" filter card")));
    await session.step(65, "And user clicks on \"Switch to free-text mode\" icon in \"Expression\" filter card", () => clickOn(page, el("\"Switch to free-text mode\" icon in \"Expression\" filter card")));
    await session.step(66, "And user types \"${HEIGHT} < 150\" into the search of the \"Expression\" filter card", () => typeIntoCardSearch(page, "${HEIGHT} < 150", "Expression"));
    await session.step(67, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(68, "Then 16 rows should pass the filter", () => filterPasses(page, 16));
    await session.step(69, "When user types \"${HEIGHT} <\" into the search of the \"Expression\" filter card", () => typeIntoCardSearch(page, "${HEIGHT} <", "Expression"));
    await session.step(70, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(71, "Then the \"categories of Expression\" reading of filter panel should be \"${HEIGHT} < 150\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${HEIGHT} < 150"));
    await session.step(72, "And 16 rows should pass the filter", () => filterPasses(page, 16));
    await session.step(73, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the column of a filtering numeric card releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(77, "When user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
    await session.step(78, "And user picks \"Min / max\" from the indicator menu of the \"AGE\" filter card", () => pickCardIndicatorMenu(page, "Min / max", "AGE"));
    await session.step(79, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(80, "And user enters \"40\" into the max field of the \"AGE\" filter card", () => enterCardBound(page, "40", "max", "AGE"));
    await session.step(81, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
    await session.step(82, "When user removes \"AGE\" column", () => removeColumn(page, "AGE"));
    await session.step(83, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(84, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(85, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the column of a filtering flag in the Combined Boolean card releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:95", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await knownFailure(async () => {
      await session.step(91, "When user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
      await session.step(92, "And user clicks on the \"category CONTROL of Flags\" area of filter panel", () => clickArea(page, "category CONTROL of Flags", el("filter panel")));
      await session.step(93, "Then 6 rows should pass the filter", () => filterPasses(page, 6));
      await session.step(94, "When user picks \"Remove\" from the context menu of the \"header CONTROL\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header CONTROL", el("grid")));
      await session.step(95, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    }, 95);
  });
  test("Removing the column of a filtering Expression card rule releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:108", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await knownFailure(async () => {
      await session.step(102, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(103, "And user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
      await session.step(104, "And user selects \"is null\" in Operation input in \"Expression\" filter card", () => selectIn(page, "is null", el("Operation input in \"Expression\" filter card")));
      await session.step(105, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(106, "Then 128 rows should pass the filter", () => filterPasses(page, 128));
      await session.step(107, "When user picks \"Remove\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header HEIGHT", el("grid")));
      await session.step(108, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(109, "And no errors should have been logged", () => noErrors(page));
    }, 108);
  });
  test("Removing a column under a Numeric Columns rule recalculates the rule over the rest", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:122", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await knownFailure(async () => {
      await session.step(115, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(116, "And user selects \"Numeric Columns\" in Column input in \"Expression\" filter card", () => selectIn(page, "Numeric Columns", el("Column input in \"Expression\" filter card")));
      await session.step(117, "And user selects \">\" in Operation input in \"Expression\" filter card", () => selectIn(page, ">", el("Operation input in \"Expression\" filter card")));
      await session.step(118, "And user enters \"100\" into Value input in \"Expression\" filter card", () => enterInto(page, "100", el("Value input in \"Expression\" filter card")));
      await session.step(119, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(120, "Then 884 rows should pass the filter", () => filterPasses(page, 884));
      await session.step(121, "When user picks \"Remove\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header HEIGHT", el("grid")));
      await session.step(122, "Then 131 rows should pass the filter", () => filterPasses(page, 131));
      await session.step(123, "And no errors should have been logged", () => noErrors(page));
    }, 122);
  });
  test("The categorical search for \"a(b\" finds the category with that literal text [typed=a(b]", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@realizes:core#1134"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(127, "Given user opens a table \"special-categories\" with:", () => openTableOf(page, "special-categories", [["Name","Value"],["a(b","1"],["a.b","2"],["axb","3"],["a+b","4"],["aab","5"],["a[b","6"]]), [["Name","Value"],["a(b","1"],["a.b","2"],["axb","3"],["a+b","4"],["aab","5"],["a[b","6"]]);
    await session.step(135, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(136, "And user adds a card for \"Name\" to the filter panel", () => addCardFor(page, "Name"));
    await session.step(137, "And user hovers over \"Name\" filter card", () => hoverOver(page, el("\"Name\" filter card")));
    await session.step(138, "And user clicks on search icon of \"Name\" filter card", () => clickOn(page, el("search icon of \"Name\" filter card")));
    await session.step(139, "And user types \"a(b\" into the search of the \"Name\" filter card", () => typeIntoCardSearch(page, "a(b", "Name"));
    await session.step(140, "Then 1 row should pass the filter", () => filterPasses(page, 1));
    await session.step(141, "And the filter should pass exactly the rows where \"Name\" is \"a(b\"", () => filterIsExactlyCategory(page, "Name", "a(b"));
    await session.step(142, "And no errors should have been logged", () => noErrors(page));
  });
  test("The categorical search for \"a.b\" finds the category with that literal text [typed=a.b]", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@realizes:core#1134"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(127, "Given user opens a table \"special-categories\" with:", () => openTableOf(page, "special-categories", [["Name","Value"],["a(b","1"],["a.b","2"],["axb","3"],["a+b","4"],["aab","5"],["a[b","6"]]), [["Name","Value"],["a(b","1"],["a.b","2"],["axb","3"],["a+b","4"],["aab","5"],["a[b","6"]]);
    await session.step(135, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(136, "And user adds a card for \"Name\" to the filter panel", () => addCardFor(page, "Name"));
    await session.step(137, "And user hovers over \"Name\" filter card", () => hoverOver(page, el("\"Name\" filter card")));
    await session.step(138, "And user clicks on search icon of \"Name\" filter card", () => clickOn(page, el("search icon of \"Name\" filter card")));
    await session.step(139, "And user types \"a.b\" into the search of the \"Name\" filter card", () => typeIntoCardSearch(page, "a.b", "Name"));
    await session.step(140, "Then 1 row should pass the filter", () => filterPasses(page, 1));
    await session.step(141, "And the filter should pass exactly the rows where \"Name\" is \"a.b\"", () => filterIsExactlyCategory(page, "Name", "a.b"));
    await session.step(142, "And no errors should have been logged", () => noErrors(page));
  });
  test("The categorical search for \"a+b\" finds the category with that literal text [typed=a+b]", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@realizes:core#1134"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(127, "Given user opens a table \"special-categories\" with:", () => openTableOf(page, "special-categories", [["Name","Value"],["a(b","1"],["a.b","2"],["axb","3"],["a+b","4"],["aab","5"],["a[b","6"]]), [["Name","Value"],["a(b","1"],["a.b","2"],["axb","3"],["a+b","4"],["aab","5"],["a[b","6"]]);
    await session.step(135, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(136, "And user adds a card for \"Name\" to the filter panel", () => addCardFor(page, "Name"));
    await session.step(137, "And user hovers over \"Name\" filter card", () => hoverOver(page, el("\"Name\" filter card")));
    await session.step(138, "And user clicks on search icon of \"Name\" filter card", () => clickOn(page, el("search icon of \"Name\" filter card")));
    await session.step(139, "And user types \"a+b\" into the search of the \"Name\" filter card", () => typeIntoCardSearch(page, "a+b", "Name"));
    await session.step(140, "Then 1 row should pass the filter", () => filterPasses(page, 1));
    await session.step(141, "And the filter should pass exactly the rows where \"Name\" is \"a+b\"", () => filterIsExactlyCategory(page, "Name", "a+b"));
    await session.step(142, "And no errors should have been logged", () => noErrors(page));
  });
  test("The categorical search for \"a[b\" finds the category with that literal text [typed=a[b]", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@realizes:core#1134"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(127, "Given user opens a table \"special-categories\" with:", () => openTableOf(page, "special-categories", [["Name","Value"],["a(b","1"],["a.b","2"],["axb","3"],["a+b","4"],["aab","5"],["a[b","6"]]), [["Name","Value"],["a(b","1"],["a.b","2"],["axb","3"],["a+b","4"],["aab","5"],["a[b","6"]]);
    await session.step(135, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(136, "And user adds a card for \"Name\" to the filter panel", () => addCardFor(page, "Name"));
    await session.step(137, "And user hovers over \"Name\" filter card", () => hoverOver(page, el("\"Name\" filter card")));
    await session.step(138, "And user clicks on search icon of \"Name\" filter card", () => clickOn(page, el("search icon of \"Name\" filter card")));
    await session.step(139, "And user types \"a[b\" into the search of the \"Name\" filter card", () => typeIntoCardSearch(page, "a[b", "Name"));
    await session.step(140, "Then 1 row should pass the filter", () => filterPasses(page, 1));
    await session.step(141, "And the filter should pass exactly the rows where \"Name\" is \"a[b\"", () => filterIsExactlyCategory(page, "Name", "a[b"));
    await session.step(142, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the column of a hierarchical card's middle level drops that level", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:158", "@realizes:core#1223"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await knownFailure(async () => {
      await session.step(154, "When user picks \"Add Filter | Hierarchical\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Hierarchical"));
      await session.step(155, "And user configures the hierarchical filter with columns \"SEX, RACE, SEVERITY\"", () => configureHierarchical(page, "SEX, RACE, SEVERITY"));
      await session.step(156, "Then \"SEX / RACE / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE / SEVERITY\" filter card"), "visible"));
      await session.step(157, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(158, "Then \"SEX / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / SEVERITY\" filter card"), "visible"));
      await session.step(159, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(160, "When user expands the \"F\" row of the hierarchical filter card", () => expandHierarchicalRow(page, "F"));
      await session.step(161, "Then the hierarchical filter card should list the \"F / None\" row", () => hierarchicalListsRow(page, "F / None"));
      await session.step(162, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(163, "And no errors should have been logged", () => noErrors(page));
    }, 158);
  });
  test("A layout naming a removed level column brings the hierarchical card back without that level", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:174", "@realizes:core#1223"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await knownFailure(async () => {
      await session.step(168, "When user picks \"Add Filter | Hierarchical\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Hierarchical"));
      await session.step(169, "And user configures the hierarchical filter with columns \"SEX, RACE, SEVERITY\"", () => configureHierarchical(page, "SEX, RACE, SEVERITY"));
      await session.step(170, "Then \"SEX / RACE / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE / SEVERITY\" filter card"), "visible"));
      await session.step(171, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(172, "And user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(173, "And user loads the saved layout", () => loadLayout(page));
      await session.step(174, "Then \"SEX / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / SEVERITY\" filter card"), "visible"));
      await session.step(175, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(176, "When user expands the \"F\" row of the hierarchical filter card", () => expandHierarchicalRow(page, "F"));
      await session.step(177, "Then the hierarchical filter card should list the \"F / None\" row", () => hierarchicalListsRow(page, "F / None"));
      await session.step(178, "And no errors should have been logged", () => noErrors(page));
    }, 174);
  });
  test("Saving a filter state under a name already taken asks before overwriting it", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:196", "@realizes:core#1224"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await knownFailure(async () => {
      await session.step(183, "Given no saved filter state \"bdd edge filter state\" is kept, now or when the feature ends", () => forgetSavedState(page, "bdd edge filter state"));
      await session.step(184, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(185, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
      await session.step(186, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(187, "When user picks \"Save or Apply | Save...\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply | Save..."));
      await session.step(188, "And user types \"bdd edge filter state\" into Name input in \"Save filter preset\" dialog", () => typeInto(page, "bdd edge filter state", el("Name input in \"Save filter preset\" dialog")));
      await session.step(189, "And user clicks on OK button in \"Save filter preset\" dialog", () => clickOn(page, el("OK button in \"Save filter preset\" dialog")));
      await session.step(190, "Then \"Save filter preset\" dialog should be absent", () => shouldBe(page, el("\"Save filter preset\" dialog"), "absent"));
      await session.step(191, "When user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(192, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(193, "When user picks \"Save or Apply | Save...\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply | Save..."));
      await session.step(194, "And user types \"bdd edge filter state\" into Name input in \"Save filter preset\" dialog", () => typeInto(page, "bdd edge filter state", el("Name input in \"Save filter preset\" dialog")));
      await session.step(195, "And user clicks on OK button in \"Save filter preset\" dialog", () => clickOn(page, el("OK button in \"Save filter preset\" dialog")));
      await session.step(196, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(197, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
      await session.step(198, "And user picks \"Save or Apply | bdd edge filter state\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply | bdd edge filter state"));
      await session.step(199, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(200, "And the \"selected categories of RACE\" reading of filter panel should be \"Asian\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Asian"));
      await session.step(201, "And no errors should have been logged", () => noErrors(page));
    }, 196);
  });
  test("A Combined Boolean card added to a table with no boolean column picks up a boolean column added later", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(205, "Given user opens a table \"no-flags\" with:", () => openTableOf(page, "no-flags", [["Name","Value"],["a","1"],["b","2"]]), [["Name","Value"],["a","1"],["b","2"]]);
    await session.step(209, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(210, "And user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
    await session.step(211, "Then the open menu should list \"Add Filter > Combined Boolean\" as enabled", () => menuListsEnabled(page, "Add Filter > Combined Boolean"));
    await session.step(212, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(213, "And user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
    await session.step(214, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
    await session.step(215, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
    await session.step(216, "When user adds a calculated column \"BIG\" with formula \"${Value} > 1\"", () => addCalculated(page, "BIG", "${Value} > 1"));
    await session.step(217, "And user clicks on the \"category BIG of Flags\" area of filter panel", () => clickArea(page, "category BIG of Flags", el("filter panel")));
    await session.step(218, "Then 1 row should pass the filter", () => filterPasses(page, 1));
    await session.step(219, "And no errors should have been logged", () => noErrors(page));
  });
  test("Closing the view of a table whose filter narrows a linked table releases the linked table's filter", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(223, "Given user opens a table \"races\" with:", () => openTableOf(page, "races", [["RACE"],["Asian"],["Black"],["Caucasian"],["Other"]]), [["RACE"],["Asian"],["Black"],["Caucasian"],["Other"]]);
    await session.step(229, "And the \"races\" table is linked to the \"demog-1000\" table by \"RACE\" to \"RACE\" as \"filter to filter\"", () => linkTables(page, "races", "demog-1000", "RACE", "RACE", "filter to filter"));
    await session.step(230, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(231, "And user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
    await session.step(232, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
    await session.step(233, "Then 1 row should pass the filter", () => filterPasses(page, 1));
    await session.step(234, "And 15 rows of table \"demog-1000\" should pass the filter", () => tableFilterCount(page, 15, "demog-1000"));
    await session.step(235, "When user closes the current view", () => closeCurrentView(page));
    await session.step(236, "Then 1000 rows of table \"demog-1000\" should pass the filter", () => tableFilterCount(page, 1000, "demog-1000"));
    await session.step(237, "And no errors should have been logged", () => noErrors(page));
  });
});
