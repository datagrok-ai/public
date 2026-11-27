/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/04-expression-filter.feature
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
import {clickOn, hoverOver, selectIn, shouldBe, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactly, filterPasses, filterPassesAll, openFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, pickFromOpenMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/04-expression-filter.feature", import.meta.url);
  test("Filter with an expression", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(16, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(19, "And user opens the filter panel", () => openFilterPanel(page));
    await session.step(20, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(21, "When user hovers over \"Filters\" dock panel", () => hoverOver(page, el("\"Filters\" dock panel")));
    await session.step(23, "And user clicks on \"font-icon-menu\" icon in \"Filters\" dock panel", () => clickOn(page, el("\"font-icon-menu\" icon in \"Filters\" dock panel")), undefined, "The panel's menu adds an expression filter");
    await session.step(24, "And user picks \"Add Filter | Expression\" from the open menu", () => pickFromOpenMenu(page, "Add Filter | Expression"));
    await session.step(25, "Then \"Expression\" filter card should be visible", () => shouldBe(page, el("\"Expression\" filter card"), "visible"));
    await session.step(26, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(28, "When user selects \"AGE\" in Column input in \"Expression\" filter card", () => selectIn(page, "AGE", el("Column input in \"Expression\" filter card")), undefined, "Pick a column, an operation and a value");
    await session.step(29, "And user selects \">\" in Operation input in \"Expression\" filter card", () => selectIn(page, ">", el("Operation input in \"Expression\" filter card")));
    await session.step(30, "And user types \"50\" into Value input in \"Expression\" filter card", () => typeInto(page, "50", el("Value input in \"Expression\" filter card")));
    await session.step(31, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
    await session.step(32, "Then 2176 rows should pass the filter", () => filterPasses(page, 2176));
    await session.step(33, "And the filter should pass exactly the rows where \"AGE\" is between 51 and 200", () => filterIsExactly(page, "AGE", 51, 200));
    await session.step(35, "When user selects \"HEIGHT\" in Column input in \"Expression\" filter card", () => selectIn(page, "HEIGHT", el("Column input in \"Expression\" filter card")), undefined, "Add another rule");
    await session.step(36, "And user selects \"<\" in Operation input in \"Expression\" filter card", () => selectIn(page, "<", el("Operation input in \"Expression\" filter card")));
    await session.step(37, "And user types \"160\" into Value input in \"Expression\" filter card", () => typeInto(page, "160", el("Value input in \"Expression\" filter card")));
    await session.step(38, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
    await session.step(39, "Then the \"categories of Expression\" reading of filter panel should be \"${AGE} > 50, ${HEIGHT} < 160\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${AGE} > 50, ${HEIGHT} < 160"));
    await session.step(40, "And 2628 rows should pass the filter", () => filterPasses(page, 2628));
    await session.step(42, "When user clicks on mode of \"Expression\" filter card", () => clickOn(page, el("mode of \"Expression\" filter card")), undefined, "OR keeps rows that match either rule; click it for AND");
    await session.step(43, "Then mode of \"Expression\" filter card should have text \"AND\"", () => shouldHaveText(page, el("mode of \"Expression\" filter card"), "AND"));
    await session.step(44, "And 499 rows should pass the filter", () => filterPasses(page, 499));
    await session.step(46, "When user clicks on the \"checkbox ${HEIGHT} < 160 of Expression\" area of filter panel", () => clickArea(page, "checkbox ${HEIGHT} < 160 of Expression", el("filter panel")), undefined, "A rule's checkbox turns it off");
    await session.step(47, "Then 2176 rows should pass the filter", () => filterPasses(page, 2176));
    await session.step(48, "And the \"selected categories of Expression\" reading of filter panel should be \"${AGE} > 50\"", () => readingReads(page, "selected categories of Expression", el("filter panel"), "${AGE} > 50"));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
