/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/11-filter-to-a-list-of-ids-and-copy-the-rows.feature
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
import {clickOn, hoverOver, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterIsExactlyAnyOf, filterPasses, filterPassesAll, openFilterPanel, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelClosed, clickViewTab, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, pasteIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {menuLists, noErrors, pickFromOpenMenu, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/11-filter-to-a-list-of-ids-and-copy-the-rows.feature", import.meta.url);
  test("Filter to a list of IDs and copy the rows into a new table", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(24, "And user opens the filter panel", () => openFilterPanel(page));
    await session.step(25, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(26, "When user hovers over caption of \"Structure\" filter card", () => hoverOver(page, el("caption of \"Structure\" filter card")));
    await session.step(28, "And user adds a card for \"Id\" to the filter panel", () => addCardFor(page, "Id"), undefined, "The + adds a filter for any column");
    await session.step(29, "Then \"Id\" filter card should be visible", () => shouldBe(page, el("\"Id\" filter card"), "visible"));
    await session.step(30, "When user hovers over \"Id\" filter card", () => hoverOver(page, el("\"Id\" filter card")));
    await session.step(32, "And user clicks on search icon of \"Id\" filter card", () => clickOn(page, el("search icon of \"Id\" filter card")), undefined, "The search icon opens the card's search box");
    await session.step(33, "And user clicks on search of \"Id\" filter card", () => clickOn(page, el("search of \"Id\" filter card")));
    await session.step(35, "And user pastes \"CAST-634783\\nCAST-634790\\nCAST-634812\\nCAST-634851\\nCAST-634880\" into the search of the \"Id\" filter card", () => pasteIntoCardSearch(page, "CAST-634783\\nCAST-634790\\nCAST-634812\\nCAST-634851\\nCAST-634880", "Id"), undefined, "Paste the list, one ID per line");
    await session.step(36, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
    await session.step(37, "And the filter should pass exactly the rows where \"Id\" is one of \"CAST-634783, CAST-634790, CAST-634812, CAST-634851, CAST-634880\"", () => filterIsExactlyAnyOf(page, "Id", "CAST-634783, CAST-634790, CAST-634812, CAST-634851, CAST-634880"));
    await session.step(38, "And grid should show 5 rows", () => showsRows(page, el("grid"), 5));
    await session.step(39, "When the video says \"Only the five listed compounds are left\"", () => videoSays(page, "Only the five listed compounds are left"));
    await session.step(41, "And user clicks on \"Filtered: 5\" text in status bar", () => clickOn(page, el("\"Filtered: 5\" text in status bar")), undefined, "The Filtered count's menu: keep, delete or copy out these rows");
    await session.step(42, "Then the open menu should list \"Keep Rows\"", () => menuLists(page, "Keep Rows"));
    await session.step(43, "And the open menu should list \"Delete Rows\"", () => menuLists(page, "Delete Rows"));
    await session.step(45, "When user picks \"Extract Rows\" from the open menu", () => pickFromOpenMenu(page, "Extract Rows"), undefined, "Extract Rows copies them into a new table");
    await session.step(46, "Then the table should have 5 rows", () => rowCount(page, 5));
    await session.step(47, "And every value of \"Id\" column should match \"^CAST-634(783|790|812|851|880)$\"", () => everyValueMatches(page, "Id", "^CAST-634(783|790|812|851|880)$"));
    await session.step(49, "When user hovers over \"Rows: 5\" text in status bar", () => hoverOver(page, el("\"Rows: 5\" text in status bar")), undefined, "The five compounds, now a table of their own");
    await session.step(50, "And user clicks on the tab of the \"spgi-100\" view", () => clickViewTab(page, "spgi-100"));
    await session.step(52, "And user hovers over \"Rows: 100\" text in status bar", () => hoverOver(page, el("\"Rows: 100\" text in status bar")), undefined, "The original table still has all 100 rows");
    await session.step(53, "Then the table should have 100 rows", () => rowCount(page, 100));
    await session.step(54, "And 5 rows should pass the filter", () => filterPasses(page, 5));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
  });
});
