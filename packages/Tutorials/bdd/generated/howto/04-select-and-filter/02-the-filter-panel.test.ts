/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/02-the-filter-panel.feature
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
import {clickOn, hoverOver, hoverText, shouldBe, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactlyCategory, filterPasses, filterPassesAll, openEmptyFilterPanel, passingFewer, passingMore, rememberPassing} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {addViewerWith, clickArea, dragAreaToArea, hoverArea, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/02-the-filter-panel.feature", import.meta.url);
  test("Filter rows with the filter panel", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(23, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(24, "And user hovers over filter panel", () => hoverOver(page, el("filter panel")));
    await session.step(25, "And user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
    await session.step(26, "And user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
    await session.step(27, "And user adds a card for \"SEX\" to the filter panel", () => addCardFor(page, "SEX"));
    await session.step(28, "And user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
    await session.step(29, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
    await session.step(31, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")), undefined, "The filter icon opens the filter panel");
    await session.step(32, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
    await session.step(33, "And the \"cards\" reading of filter panel should be \"SEX, AGE, RACE\"", () => readingReads(page, "cards", el("filter panel"), "SEX, AGE, RACE"));
    await session.step(34, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(36, "When user clicks on the \"category F of SEX\" area of filter panel", () => clickArea(page, "category F of SEX", el("filter panel")), undefined, "Click a category to keep only its rows");
    await session.step(37, "Then 3243 rows should pass the filter", () => filterPasses(page, 3243));
    await session.step(38, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
    await session.step(39, "When the video says \"The grid and the plot show only the rows that pass\"", () => videoSays(page, "The grid and the plot show only the rows that pass"));
    await session.step(40, "And user remembers how many rows pass the filter", () => rememberPassing(page));
    await session.step(42, "And user hovers over the \"bin 1\" area of histogram viewer in \"AGE\" filter card", () => hoverArea(page, "bin 1", el("histogram viewer in \"AGE\" filter card")), undefined, "Drag the ends of a range to narrow it");
    await session.step(43, "And user drags the \"range min handle\" area of histogram viewer in \"AGE\" filter card to the \"bin 7\" area", () => dragAreaToArea(page, "range min handle", el("histogram viewer in \"AGE\" filter card"), "bin 7"));
    await session.step(44, "Then fewer rows than remembered should pass the filter", () => passingFewer(page));
    await session.step(45, "When user hovers over the text \"Filtered\" in status bar", () => hoverText(page, "Filtered", el("status bar")));
    await session.step(46, "And the video says \"The status bar counts the rows that pass\"", () => videoSays(page, "The status bar counts the rows that pass"));
    await session.step(47, "And user remembers how many rows pass the filter", () => rememberPassing(page));
    await session.step(49, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")), undefined, "A card's checkbox turns that filter off");
    await session.step(50, "And user unchecks checkbox of \"SEX\" filter card", () => uncheck(page, el("checkbox of \"SEX\" filter card")));
    await session.step(51, "Then more rows than remembered should pass the filter", () => passingMore(page));
    await session.step(53, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")), undefined, "Reset clears every filter");
    await session.step(54, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
    await session.step(55, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
  });
});
