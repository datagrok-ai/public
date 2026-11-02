/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-range-handles.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters, public#1343, core#1174, core#1223, core#1224]
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
import {clickOn, enterInto, selectIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {configureHierarchical, filterPasses, filterPassesAll, openEmptyFilterPanel, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, forgetSavedState, pickPanelMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, dragAreaBy, dragAreaToArea, hoverArea, loadLayout, noErrors, pickFromAreaContextMenu, readingHigher, readingHigherThanRemembered, readingIs, readingLowerThanRemembered, readingReads, rememberReading, saveLayout} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("The range handles of a numeric filter card", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-range-handles.feature", import.meta.url);
  test("The min and max handles narrow the range, and dragged back past the ends they keep every row (public#1343)", {tag: ["@viewers", "@realizes:viewers.filters", "@realizes:public#1343"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(27, "And user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
    await session.step(28, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(29, "And the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
    await session.step(33, "When user remembers the \"min of AGE\" reading of filter panel", () => rememberReading(page, "min of AGE", el("filter panel")));
    await session.step(34, "And user hovers over the \"view\" area of histogram viewer in filter panel", () => hoverArea(page, "view", el("histogram viewer in filter panel")));
    await session.step(35, "And user drags the \"range min handle\" area of histogram viewer in filter panel to the \"bin 5\" area", () => dragAreaToArea(page, "range min handle", el("histogram viewer in filter panel"), "bin 5"));
    await session.step(36, "Then the \"min of AGE\" reading of filter panel should be higher than remembered", () => readingHigherThanRemembered(page, "min of AGE", el("filter panel")));
    await session.step(37, "And the \"range min\" reading of histogram viewer in filter panel should be higher than before", () => readingHigher(page, "range min", el("histogram viewer in filter panel")));
    await session.step(38, "And 776 rows should pass the filter", () => filterPasses(page, 776));
    await session.step(39, "When user remembers the \"max of AGE\" reading of filter panel", () => rememberReading(page, "max of AGE", el("filter panel")));
    await session.step(40, "And user hovers over the \"view\" area of histogram viewer in filter panel", () => hoverArea(page, "view", el("histogram viewer in filter panel")));
    await session.step(41, "And user drags the \"range max handle\" area of histogram viewer in filter panel to the \"bin 15\" area", () => dragAreaToArea(page, "range max handle", el("histogram viewer in filter panel"), "bin 15"));
    await session.step(42, "Then the \"max of AGE\" reading of filter panel should be lower than remembered", () => readingLowerThanRemembered(page, "max of AGE", el("filter panel")));
    await session.step(43, "And 730 rows should pass the filter", () => filterPasses(page, 730));
    await session.step(44, "And the \"filtering of AGE\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of AGE", el("filter panel"), "true"));
    await session.step(45, "When user hovers over the \"view\" area of histogram viewer in filter panel", () => hoverArea(page, "view", el("histogram viewer in filter panel")));
    await session.step(46, "And user drags the \"range min handle\" area of histogram viewer in filter panel by 400 pixels to the left", () => dragAreaBy(page, "range min handle", el("histogram viewer in filter panel"), 400, "left"));
    await session.step(47, "And user hovers over the \"view\" area of histogram viewer in filter panel", () => hoverArea(page, "view", el("histogram viewer in filter panel")));
    await session.step(48, "And user drags the \"range max handle\" area of histogram viewer in filter panel by 400 pixels to the right", () => dragAreaBy(page, "range max handle", el("histogram viewer in filter panel"), 400, "right"));
    await session.step(49, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(50, "And the \"min of AGE\" reading of filter panel should be 18", () => readingIs(page, "min of AGE", el("filter panel"), 18));
    await session.step(51, "And the \"max of AGE\" reading of filter panel should be 89", () => readingIs(page, "max of AGE", el("filter panel"), 89));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the column of a filtering flag in the Combined Boolean card releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@known-failure:60", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(27, "And user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
    await session.step(28, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(29, "And the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
    await knownFailure(async () => {
      await session.step(56, "When user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
      await session.step(57, "And user clicks on the \"category CONTROL of Flags\" area of filter panel", () => clickArea(page, "category CONTROL of Flags", el("filter panel")));
      await session.step(58, "Then 6 rows should pass the filter", () => filterPasses(page, 6));
      await session.step(59, "When user picks \"Remove\" from the context menu of the \"header CONTROL\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header CONTROL", el("grid")));
      await session.step(60, "Then all rows should pass the filter", () => filterPassesAll(page));
    }, 60);
  });
  test("Removing the column of a filtering Expression card rule releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@known-failure:70", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(27, "And user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
    await session.step(28, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(29, "And the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
    await knownFailure(async () => {
      await session.step(64, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(65, "And user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")));
      await session.step(66, "And user selects \"is null\" in Operation input in \"Expression\" filter card", () => selectIn(page, "is null", el("Operation input in \"Expression\" filter card")));
      await session.step(67, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(68, "Then 128 rows should pass the filter", () => filterPasses(page, 128));
      await session.step(69, "When user picks \"Remove\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header HEIGHT", el("grid")));
      await session.step(70, "Then all rows should pass the filter", () => filterPassesAll(page));
    }, 70);
  });
  test("Removing a column under a Numeric Columns rule works the rule out again over the rest", {tag: ["@viewers", "@realizes:viewers.filters", "@known-failure:81", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(27, "And user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
    await session.step(28, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(29, "And the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
    await knownFailure(async () => {
      await session.step(74, "When user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(75, "And user selects \"Numeric Columns\" in Column input in \"Expression\" filter card", () => selectIn(page, "Numeric Columns", el("Column input in \"Expression\" filter card")));
      await session.step(76, "And user selects \">\" in Operation input in \"Expression\" filter card", () => selectIn(page, ">", el("Operation input in \"Expression\" filter card")));
      await session.step(77, "And user enters \"100\" into Value input in \"Expression\" filter card", () => enterInto(page, "100", el("Value input in \"Expression\" filter card")));
      await session.step(78, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(79, "Then 884 rows should pass the filter", () => filterPasses(page, 884));
      await session.step(80, "When user picks \"Remove\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header HEIGHT", el("grid")));
      await session.step(81, "Then 131 rows should pass the filter", () => filterPasses(page, 131));
    }, 81);
  });
  test("Removing the column of a hierarchical card's middle level drops that level", {tag: ["@viewers", "@realizes:viewers.filters", "@known-failure:89", "@realizes:core#1223"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(27, "And user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
    await session.step(28, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(29, "And the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
    await knownFailure(async () => {
      await session.step(85, "When user picks \"Add Filter | Hierarchical\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Hierarchical"));
      await session.step(86, "And user configures the hierarchical filter with columns \"SEX, RACE, SEVERITY\"", () => configureHierarchical(page, "SEX, RACE, SEVERITY"));
      await session.step(87, "Then \"SEX / RACE / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE / SEVERITY\" filter card"), "visible"));
      await session.step(88, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(89, "Then \"SEX / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / SEVERITY\" filter card"), "visible"));
    }, 89);
  });
  test("A layout naming a removed level column brings the hierarchical card back without that level", {tag: ["@viewers", "@realizes:viewers.filters", "@known-failure:99", "@realizes:core#1223"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(27, "And user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
    await session.step(28, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(29, "And the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
    await knownFailure(async () => {
      await session.step(93, "When user picks \"Add Filter | Hierarchical\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Hierarchical"));
      await session.step(94, "And user configures the hierarchical filter with columns \"SEX, RACE, SEVERITY\"", () => configureHierarchical(page, "SEX, RACE, SEVERITY"));
      await session.step(95, "Then \"SEX / RACE / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE / SEVERITY\" filter card"), "visible"));
      await session.step(96, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(97, "And user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(98, "And user loads the saved layout", () => loadLayout(page));
      await session.step(99, "Then \"SEX / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / SEVERITY\" filter card"), "visible"));
    }, 99);
  });
  test("Saving a filter state under a name already taken asks before overwriting it", {tag: ["@viewers", "@realizes:viewers.filters", "@known-failure:116", "@realizes:core#1224"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(27, "And user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
    await session.step(28, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(29, "And the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
    await knownFailure(async () => {
      await session.step(103, "Given no saved filter state \"bdd edge filter state\" is kept, now or when the feature ends", () => forgetSavedState(page, "bdd edge filter state"));
      await session.step(104, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(105, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
      await session.step(106, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(107, "When user picks \"Save or Apply | Save...\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply | Save..."));
      await session.step(108, "And user types \"bdd edge filter state\" into Name input in \"Save filter preset\" dialog", () => typeInto(page, "bdd edge filter state", el("Name input in \"Save filter preset\" dialog")));
      await session.step(109, "And user clicks on OK button in \"Save filter preset\" dialog", () => clickOn(page, el("OK button in \"Save filter preset\" dialog")));
      await session.step(110, "Then \"Save filter preset\" dialog should be absent", () => shouldBe(page, el("\"Save filter preset\" dialog"), "absent"));
      await session.step(111, "When user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(112, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(113, "When user picks \"Save or Apply | Save...\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply | Save..."));
      await session.step(114, "And user types \"bdd edge filter state\" into Name input in \"Save filter preset\" dialog", () => typeInto(page, "bdd edge filter state", el("Name input in \"Save filter preset\" dialog")));
      await session.step(115, "And user clicks on OK button in \"Save filter preset\" dialog", () => clickOn(page, el("OK button in \"Save filter preset\" dialog")));
      await session.step(116, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    }, 116);
  });
});
