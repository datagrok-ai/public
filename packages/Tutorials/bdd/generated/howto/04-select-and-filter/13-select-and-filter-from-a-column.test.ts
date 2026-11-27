/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/13-select-and-filter-from-a-column.feature
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
import {clickOn, hoverOver, pressKey, pressKeyIn, shouldBe, shouldHaveText, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll, noneSelected, onlyBetweenSelected, onlyOfSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, dragAreaToArea, hoverArea, noErrors, readingHigher} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/13-select-and-filter-from-a-column.feature", import.meta.url);
  test("Select and filter from a column", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "Then no rows should be selected", () => noneSelected(page));
    await session.step(22, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(24, "When user clicks on the \"cell 4 of RACE\" area of grid", () => clickArea(page, "cell 4 of RACE", el("grid")), undefined, "Click a cell, then press Shift+Enter");
    await session.step(25, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
    await session.step(26, "Then 354 rows should be selected", () => selectedRowCount(page, 354));
    await session.step(27, "And only rows where \"RACE\" is \"Other\" should be selected", () => onlyOfSelected(page, "RACE", "Other"));
    await session.step(28, "When user hovers over \"Selected: 354\" text in status bar", () => hoverOver(page, el("\"Selected: 354\" text in status bar")));
    await session.step(29, "And the video says \"Every row with RACE = Other is selected\"", () => videoSays(page, "Every row with RACE = Other is selected"));
    await session.step(31, "When user presses Escape in grid", () => pressKeyIn(page, "Escape", el("grid")), undefined, "Esc clears the selection");
    await session.step(32, "Then no rows should be selected", () => noneSelected(page));
    await session.step(34, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")), undefined, "In the Search box (Ctrl+F), Shift+Enter selects the matches");
    await session.step(35, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await session.step(36, "When user types \"AGE > 70\" into table search", () => typeInto(page, "AGE > 70", el("table search")));
    await session.step(37, "And user presses Shift+Enter in table search", () => pressKeyIn(page, "Shift+Enter", el("table search")));
    await session.step(38, "Then 171 rows should be selected", () => selectedRowCount(page, 171));
    await session.step(39, "And only rows where \"AGE\" is between 71 and 200 should be selected", () => onlyBetweenSelected(page, "AGE", 71, 200));
    await session.step(40, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(41, "When the video says \"Enter would filter; Shift+Enter selects\"", () => videoSays(page, "Enter would filter; Shift+Enter selects"));
    await session.step(42, "And user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
    await session.step(44, "And user clicks on \"Column options\" icon in grid", () => clickOn(page, el("\"Column options\" icon in grid")), undefined, "The icon on a column header opens its options");
    await session.step(45, "Then column popup should be visible", () => shouldBe(page, el("column popup"), "visible"));
    await session.step(46, "And title of column popup should have text \"AGE\"", () => shouldHaveText(page, el("title of column popup"), "AGE"));
    await session.step(47, "When user hovers over the \"bin 1\" area of histogram viewer in column popup", () => hoverArea(page, "bin 1", el("histogram viewer in column popup")));
    await session.step(49, "And user drags the \"range min handle\" area of histogram viewer in column popup to the \"bin 10\" area", () => dragAreaToArea(page, "range min handle", el("histogram viewer in column popup"), "bin 10"), undefined, "Drag an end of the range to filter by this column");
    await session.step(50, "Then 1889 rows should pass the filter", () => filterPasses(page, 1889));
    await session.step(51, "And the \"range min\" reading of histogram viewer in column popup should be higher than before", () => readingHigher(page, "range min", el("histogram viewer in column popup")));
    await session.step(52, "And 171 rows should be selected", () => selectedRowCount(page, 171));
    await session.step(54, "When user clicks on \"Add filter\" action in column popup", () => clickOn(page, el("\"Add filter\" action in column popup")), undefined, "Add filter moves it to the filter panel");
    await session.step(55, "Then \"AGE\" filter card should be visible", () => shouldBe(page, el("\"AGE\" filter card"), "visible"));
    await session.step(56, "And 1889 rows should pass the filter", () => filterPasses(page, 1889));
    await session.step(58, "When user presses Escape", () => pressKey(page, "Escape"), undefined, "Esc closes the popup; the filter stays");
    await session.step(59, "Then column popup should be absent", () => shouldBe(page, el("column popup"), "absent"));
    await session.step(60, "And 1889 rows should pass the filter", () => filterPasses(page, 1889));
    await session.step(61, "And no errors should have been logged", () => noErrors(page));
  });
});
