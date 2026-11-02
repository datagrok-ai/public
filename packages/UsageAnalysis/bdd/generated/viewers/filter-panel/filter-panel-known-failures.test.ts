/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters, core#1174, core#1223, core#1224]
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
import {clickOn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, addRangeFilter, configureHierarchical, filterPanelCount, filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {forgetSavedState, pickPanelMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, loadLayout, pickFromAreaContextMenu, readingReads, saveLayout} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel known failures", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-known-failures.feature", import.meta.url);
  test("Removing the column of a filtering flag in the Combined Boolean card releases its rows", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:29", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(18, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(19, "And all rows should pass the filter", () => filterPassesAll(page));
    await knownFailure(async () => {
      await session.step(23, "When user picks \"Add Filter | Combined Boolean\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Combined Boolean"));
      await session.step(24, "Then \"Flags\" filter card should be visible", () => shouldBe(page, el("\"Flags\" filter card"), "visible"));
      await session.step(25, "And the \"categories of Flags\" reading of filter panel should be \"CONTROL\"", () => readingReads(page, "categories of Flags", el("filter panel"), "CONTROL"));
      await session.step(26, "When user clicks on the \"category CONTROL of Flags\" area of filter panel", () => clickArea(page, "category CONTROL of Flags", el("filter panel")));
      await session.step(27, "Then 6 rows should pass the filter", () => filterPasses(page, 6));
      await session.step(28, "When user picks \"Remove\" from the context menu of the \"header CONTROL\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header CONTROL", el("grid")));
      await session.step(29, "Then all rows should pass the filter", () => filterPassesAll(page));
    }, 29);
  });
  test("A layout naming a removed level column brings the hierarchical card back without that level", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:39", "@realizes:core#1223"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(18, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(19, "And all rows should pass the filter", () => filterPassesAll(page));
    await knownFailure(async () => {
      await session.step(33, "When user picks \"Add Filter | Hierarchical\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Hierarchical"));
      await session.step(34, "And user configures the hierarchical filter with columns \"SEX, RACE, SEVERITY\"", () => configureHierarchical(page, "SEX, RACE, SEVERITY"));
      await session.step(35, "Then \"SEX / RACE / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE / SEVERITY\" filter card"), "visible"));
      await session.step(36, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(37, "And user picks \"Remove\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header RACE", el("grid")));
      await session.step(38, "And user loads the saved layout", () => loadLayout(page));
      await session.step(39, "Then \"SEX / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / SEVERITY\" filter card"), "visible"));
    }, 39);
  });
  test("Removing the column of a hierarchical card's middle level drops that level", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:47", "@realizes:core#1223"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(18, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(19, "And all rows should pass the filter", () => filterPassesAll(page));
    await knownFailure(async () => {
      await session.step(43, "When user picks \"Add Filter | Hierarchical\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Hierarchical"));
      await session.step(44, "And user configures the hierarchical filter with columns \"SEX, RACE, SEVERITY\"", () => configureHierarchical(page, "SEX, RACE, SEVERITY"));
      await session.step(45, "Then \"SEX / RACE / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / RACE / SEVERITY\" filter card"), "visible"));
      await session.step(46, "When user picks \"Remove\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header RACE", el("grid")));
      await session.step(47, "Then \"SEX / SEVERITY\" filter card should be visible", () => shouldBe(page, el("\"SEX / SEVERITY\" filter card"), "visible"));
    }, 47);
  });
  test("Saving a filter state under a name already taken asks before overwriting it", {tag: ["@viewers", "@realizes:viewers.filters", "@edge-cases", "@known-failure:64", "@realizes:core#1224"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(18, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(19, "And all rows should pass the filter", () => filterPassesAll(page));
    await knownFailure(async () => {
      await session.step(51, "Given no saved filter state \"bdd edge filter state\" is kept, now or when the feature ends", () => forgetSavedState(page, "bdd edge filter state"));
      await session.step(52, "When user adds a categorical filter on \"RACE\" keeping \"Asian\"", () => addCategoricalFilter(page, "RACE", "Asian"));
      await session.step(53, "And user adds a range filter on \"AGE\" from 30 to 60", () => addRangeFilter(page, "AGE", 30, 60));
      await session.step(54, "Then 8 rows should pass the filter", () => filterPasses(page, 8));
      await session.step(55, "When user picks \"Save or Apply | Save...\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply | Save..."));
      await session.step(56, "And user types \"bdd edge filter state\" into Name input in \"Save filter preset\" dialog", () => typeInto(page, "bdd edge filter state", el("Name input in \"Save filter preset\" dialog")));
      await session.step(57, "And user clicks on OK button in \"Save filter preset\" dialog", () => clickOn(page, el("OK button in \"Save filter preset\" dialog")));
      await session.step(58, "Then \"Save filter preset\" dialog should be absent", () => shouldBe(page, el("\"Save filter preset\" dialog"), "absent"));
      await session.step(59, "When user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(60, "Then the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(61, "When user picks \"Save or Apply | Save...\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply | Save..."));
      await session.step(62, "And user types \"bdd edge filter state\" into Name input in \"Save filter preset\" dialog", () => typeInto(page, "bdd edge filter state", el("Name input in \"Save filter preset\" dialog")));
      await session.step(63, "And user clicks on OK button in \"Save filter preset\" dialog", () => clickOn(page, el("OK button in \"Save filter preset\" dialog")));
      await session.step(64, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    }, 64);
  });
});
