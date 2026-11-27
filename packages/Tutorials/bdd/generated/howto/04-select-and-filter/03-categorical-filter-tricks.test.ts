/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/03-categorical-filter-tricks.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKeyIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactlyContains, filterPasses, filterPassesAll, noneSelected, onlyOfSelected, openEmptyFilterPanel, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, pickCardIndicatorMenu, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/03-categorical-filter-tricks.feature", import.meta.url);
  test("Categorical filter tricks", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(21, "And user hovers over filter panel", () => hoverOver(page, el("filter panel")));
    await session.step(22, "And user adds a card for \"DIS_POP\" to the filter panel", () => addCardFor(page, "DIS_POP"));
    await session.step(23, "Then \"DIS_POP\" filter card should be visible", () => shouldBe(page, el("\"DIS_POP\" filter card"), "visible"));
    await session.step(24, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(25, "And no rows should be selected", () => noneSelected(page));
    await session.step(27, "When user clicks on the \"category RA of DIS_POP\" area of filter panel", () => clickArea(page, "category RA of DIS_POP", el("filter panel")), undefined, "Click a name to keep only that category");
    await session.step(28, "Then 2550 rows should pass the filter", () => filterPasses(page, 2550));
    await session.step(29, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"RA\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "RA"));
    await session.step(31, "When user presses ArrowDown in \"DIS_POP\" filter card", () => pressKeyIn(page, "ArrowDown", el("\"DIS_POP\" filter card")), undefined, "The arrow keys move to the next one");
    await session.step(32, "Then the \"selected categories of DIS_POP\" reading of filter panel should be \"UC\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "UC"));
    await session.step(33, "And 728 rows should pass the filter", () => filterPasses(page, 728));
    await session.step(35, "When user clicks on the \"checkbox RA of DIS_POP\" area of filter panel", () => clickArea(page, "checkbox RA of DIS_POP", el("filter panel")), undefined, "A checkbox adds a category");
    await session.step(36, "Then the \"selected categories of DIS_POP\" reading of filter panel should be \"RA, UC\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "RA, UC"));
    await session.step(37, "And 3278 rows should pass the filter", () => filterPasses(page, 3278));
    await session.step(39, "When user hovers over caption of \"DIS_POP\" filter card", () => hoverOver(page, el("caption of \"DIS_POP\" filter card")), undefined, "Invert all swaps what is kept and what is hidden");
    await session.step(40, "And user picks \"Invert all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Invert all", "DIS_POP"));
    await session.step(41, "Then the \"selected categories of DIS_POP\" reading of filter panel should be \"AS, Indigestion, PsA, Psoriasis\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "AS, Indigestion, PsA, Psoriasis"));
    await session.step(42, "And 2572 rows should pass the filter", () => filterPasses(page, 2572));
    await session.step(44, "When user hovers over \"DIS_POP\" filter card", () => hoverOver(page, el("\"DIS_POP\" filter card")), undefined, "Search the card to find categories");
    await session.step(45, "And user clicks on search icon of \"DIS_POP\" filter card", () => clickOn(page, el("search icon of \"DIS_POP\" filter card")));
    await session.step(46, "And user clicks on search of \"DIS_POP\" filter card", () => clickOn(page, el("search of \"DIS_POP\" filter card")));
    await session.step(47, "And user types \"Ps\" into the search of the \"DIS_POP\" filter card", () => typeIntoCardSearch(page, "Ps", "DIS_POP"));
    await session.step(48, "Then 1409 rows should pass the filter", () => filterPasses(page, 1409));
    await session.step(49, "And the filter should pass exactly the rows where \"DIS_POP\" contains \"Ps\"", () => filterIsExactlyContains(page, "DIS_POP", "Ps"));
    await session.step(51, "When user clicks on the \"count Psoriasis of DIS_POP\" area of filter panel", () => clickArea(page, "count Psoriasis of DIS_POP", el("filter panel")), undefined, "The numbers count rows; click one to select them");
    await session.step(52, "Then 1209 rows should be selected", () => selectedRowCount(page, 1209));
    await session.step(53, "And only rows where \"DIS_POP\" is \"Psoriasis\" should be selected", () => onlyOfSelected(page, "DIS_POP", "Psoriasis"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
  });
});
