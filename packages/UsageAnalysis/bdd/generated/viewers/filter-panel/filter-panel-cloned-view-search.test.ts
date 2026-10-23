/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-cloned-view-search.feature
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
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterIsExactlyContains, filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, clearCardSearch, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A search typed into a card of the original view is followed by the cloned view", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-cloned-view-search.feature", import.meta.url);
  test("The original's typed search is followed by the clone's card", {tag: ["@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(16, "And user adds a card for \"DIS_POP\" to the filter panel", () => addCardFor(page, "DIS_POP"));
    await session.step(17, "When user picks \"View > Layout > Clone View\" from the top menu", () => pickFromTopMenu(page, "View > Layout > Clone View"));
    await session.step(18, "Then the \"demog-1000 copy\" view should be current", () => viewIsCurrent(page, "demog-1000 copy"));
    await session.step(19, "And \"DIS_POP\" filter card should be visible", () => shouldBe(page, el("\"DIS_POP\" filter card"), "visible"));
    await session.step(20, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(23, "When user switches to the \"demog-1000\" view", () => switchView(page, "demog-1000"));
    await session.step(24, "And user hovers over \"DIS_POP\" filter card", () => hoverOver(page, el("\"DIS_POP\" filter card")));
    await session.step(25, "And user clicks on search icon of \"DIS_POP\" filter card", () => clickOn(page, el("search icon of \"DIS_POP\" filter card")));
    await session.step(26, "And user types \"Ps\" into the search of the \"DIS_POP\" filter card", () => typeIntoCardSearch(page, "Ps", "DIS_POP"));
    await session.step(27, "Then 242 rows should pass the filter", () => filterPasses(page, 242));
    await session.step(28, "And the filter should pass exactly the rows where \"DIS_POP\" contains \"Ps\"", () => filterIsExactlyContains(page, "DIS_POP", "Ps"));
    await session.step(29, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"PsA, Psoriasis\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "PsA, Psoriasis"));
    await session.step(30, "When user switches to the \"demog-1000 copy\" view", () => switchView(page, "demog-1000 copy"));
    await session.step(31, "Then the \"selected categories of DIS_POP\" reading of filter panel should be \"PsA, Psoriasis\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "PsA, Psoriasis"));
    await session.step(32, "And the \"filtering of DIS_POP\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of DIS_POP", el("filter panel"), "true"));
    await session.step(33, "And 242 rows should pass the filter", () => filterPasses(page, 242));
    await session.step(34, "When user switches to the \"demog-1000\" view", () => switchView(page, "demog-1000"));
    await session.step(35, "And user clears the search of the \"DIS_POP\" filter card", () => clearCardSearch(page, "DIS_POP"));
    await session.step(36, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(37, "When user switches to the \"demog-1000 copy\" view", () => switchView(page, "demog-1000 copy"));
    await session.step(38, "Then the \"filtering of DIS_POP\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of DIS_POP", el("filter panel"), "false"));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
});
