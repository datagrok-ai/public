/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/06-text-search.feature
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
import {clickOn, hoverOver, pressKey, pressKeyIn, shouldBe, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactlyCategory, filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/06-text-search.feature", import.meta.url);
  test("Search text in a table", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens beer dataset", () => openDataset(page, ds("beer")));
    await session.step(18, "And the grid shows only the columns \"Styles, Origin, ABV, Aroma\"", () => gridShowsColumns(page, "Styles, Origin, ABV, Aroma"));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(23, "And user hovers over filter panel", () => hoverOver(page, el("filter panel")));
    await session.step(24, "And user adds a card for \"Aroma\" to the filter panel", () => addCardFor(page, "Aroma"));
    await session.step(25, "Then \"Aroma\" filter card should be visible", () => shouldBe(page, el("\"Aroma\" filter card"), "visible"));
    await session.step(26, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(28, "When user clicks on search of \"Aroma\" filter card", () => clickOn(page, el("search of \"Aroma\" filter card")), undefined, "A long text column gets a search box");
    await session.step(30, "And user types \"chocolate\" into the search of the \"Aroma\" filter card", () => typeIntoCardSearch(page, "chocolate", "Aroma"), undefined, "Type a word and press Enter");
    await session.step(31, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(32, "Then 25 rows should pass the filter", () => filterPasses(page, 25));
    await session.step(34, "When user types \"coffee\" into the search of the \"Aroma\" filter card", () => typeIntoCardSearch(page, "coffee", "Aroma"), undefined, "Add another word");
    await session.step(35, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(36, "Then the \"categories of Aroma\" reading of filter panel should be \"chocolate, coffee\"", () => readingReads(page, "categories of Aroma", el("filter panel"), "chocolate, coffee"));
    await session.step(37, "And 26 rows should pass the filter", () => filterPasses(page, 26));
    await session.step(39, "When user clicks on mode of \"Aroma\" filter card", () => clickOn(page, el("mode of \"Aroma\" filter card")), undefined, "OR keeps rows with either word; click for AND");
    await session.step(40, "Then mode of \"Aroma\" filter card should have text \"AND\"", () => shouldHaveText(page, el("mode of \"Aroma\" filter card"), "AND"));
    await session.step(41, "And 12 rows should pass the filter", () => filterPasses(page, 12));
    await session.step(43, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")), undefined, "Reset clears the filter");
    await session.step(44, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
    await session.step(45, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(47, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")), undefined, "Ctrl+F searches the whole table");
    await session.step(48, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await session.step(50, "When user types \"Origin = Belgium\" into table search", () => typeInto(page, "Origin = Belgium", el("table search")), undefined, "Type a condition and press Enter");
    await session.step(51, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(52, "Then 14 rows should pass the filter", () => filterPasses(page, 14));
    await session.step(53, "And the filter should pass exactly the rows where \"Origin\" is \"Belgium\"", () => filterIsExactlyCategory(page, "Origin", "Belgium"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
  });
});
