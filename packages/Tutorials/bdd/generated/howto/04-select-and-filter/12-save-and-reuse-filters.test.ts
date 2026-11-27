/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/12-save-and-reuse-filters.feature
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
import {clickOn, hoverOver, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, dialogCloses, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, forgetSavedState, pickPanelMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/12-save-and-reuse-filters.feature", import.meta.url);
  test("Save filters and apply them later", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And no saved filter state \"Asian women\" is kept, now or when the feature ends", () => forgetSavedState(page, "Asian women"));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(21, "And user hovers over filter panel", () => hoverOver(page, el("filter panel")));
    await session.step(22, "And user adds a card for \"SEX\" to the filter panel", () => addCardFor(page, "SEX"));
    await session.step(23, "And user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
    await session.step(24, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(26, "When user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")), undefined, "Set the filters you want to keep");
    await session.step(27, "Then 72 rows should pass the filter", () => filterPasses(page, 72));
    await session.step(28, "When user clicks on the \"category F of SEX\" area of filter panel", () => clickArea(page, "category F of SEX", el("filter panel")));
    await session.step(29, "Then 37 rows should pass the filter", () => filterPasses(page, 37));
    await session.step(31, "When user picks \"Save or Apply | Save...\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply | Save..."), undefined, "Save or Apply > Save... keeps these filters under a name");
    await session.step(32, "And user types \"Asian women\" into Name input in \"Save filter preset\" dialog", () => typeInto(page, "Asian women", el("Name input in \"Save filter preset\" dialog")));
    await session.step(33, "And user clicks on OK button in \"Save filter preset\" dialog", () => clickOn(page, el("OK button in \"Save filter preset\" dialog")));
    await session.step(34, "Then the \"Save filter preset\" dialog should close", () => dialogCloses(page, "Save filter preset"));
    await session.step(36, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")), undefined, "Reset clears the filters");
    await session.step(37, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
    await session.step(38, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(40, "When user picks \"Save or Apply > Asian women\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply > Asian women"), undefined, "Pick the name to bring them back");
    await session.step(41, "Then 37 rows should pass the filter", () => filterPasses(page, 37));
    await session.step(42, "And the \"selected categories of RACE\" reading of filter panel should be \"Asian\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Asian"));
    await session.step(43, "And the \"selected categories of SEX\" reading of filter panel should be \"F\"", () => readingReads(page, "selected categories of SEX", el("filter panel"), "F"));
    await session.step(44, "When the video says \"The preset is saved, so it is there next time too\"", () => videoSays(page, "The preset is saved, so it is there next time too"));
    await session.step(45, "Then no errors should have been logged", () => noErrors(page));
  });
});
