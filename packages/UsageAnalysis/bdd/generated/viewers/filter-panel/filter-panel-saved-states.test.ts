/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-saved-states.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {forgetSavedState} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addRangeFilter, filterPanelCount, filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, pickPanelMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, closeContextMenu, menuDoesNotList, menuLists, noBalloons, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {openViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Saved filter panel states", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-saved-states.feature", import.meta.url);
  test("Saved filter panel states", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And no saved filter state \"bdd filter state\" is kept, now or when the feature ends", () => forgetSavedState(page, "bdd filter state"));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(19, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(20, "And all rows should pass the filter", () => filterPassesAll(page));
    await run.scenario("A state saved by name comes back from the Save or Apply menu", async () => {
      await session.step(23, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(24, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
      await session.step(25, "And user adds a range filter on \"AGE\" from 30 to 60", () => addRangeFilter(page, "AGE", 30, 60));
      await session.step(26, "Then 8 rows should pass the filter", () => filterPasses(page, 8));
      await session.step(27, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(28, "When user picks \"Save or Apply | Save...\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply | Save..."));
      await session.step(29, "And user types \"bdd filter state\" into Name input in \"Save filter preset\" dialog", () => typeInto(page, "bdd filter state", el("Name input in \"Save filter preset\" dialog")));
      await session.step(30, "And user clicks on OK button in \"Save filter preset\" dialog", () => clickOn(page, el("OK button in \"Save filter preset\" dialog")));
      await session.step(31, "Then \"Save filter preset\" dialog should be absent", () => shouldBe(page, el("\"Save filter preset\" dialog"), "absent"));
      await session.step(32, "When user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(33, "And user adds a range filter on \"AGE\" from 18 to 89", () => addRangeFilter(page, "AGE", 18, 89));
      await session.step(34, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(35, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(36, "And the \"min of AGE\" reading of filter panel should be 18", () => readingIs(page, "min of AGE", el("filter panel"), 18));
      await session.step(37, "When user picks \"Save or Apply | bdd filter state\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply | bdd filter state"));
      await session.step(38, "Then 8 rows should pass the filter", () => filterPasses(page, 8));
      await session.step(39, "And the \"selected categories of RACE\" reading of filter panel should be \"Asian\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Asian"));
      await session.step(40, "And the \"min of AGE\" reading of filter panel should be 30", () => readingIs(page, "min of AGE", el("filter panel"), 30));
      await session.step(41, "And the \"max of AGE\" reading of filter panel should be 60", () => readingIs(page, "max of AGE", el("filter panel"), 60));
      await session.step(42, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A table of another shape is not offered the state", async () => {
      await session.step(47, "When user opens beer dataset", () => openDataset(page, ds("beer")));
      await session.step(48, "And user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(49, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(50, "When user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
      await session.step(51, "Then the open menu should list \"Save or Apply > Save...\"", () => menuLists(page, "Save or Apply > Save..."));
      await session.step(52, "And the open menu should not list \"Save or Apply > bdd filter state\"", () => menuDoesNotList(page, "Save or Apply > bdd filter state"));
      await session.step(53, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(54, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
