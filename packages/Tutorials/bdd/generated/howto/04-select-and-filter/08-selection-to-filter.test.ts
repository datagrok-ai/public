/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/08-selection-to-filter.feature
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
import {clickOn, enterInto, hoverOver, shouldBe, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterIsExactlyCategory, filterPasses, filterPassesAll, noneSelected, selectedPassFilter, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, dialogCloses, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridShowsColumns, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/08-selection-to-filter.feature", import.meta.url);
  test("Turn a selection into a filter", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(16, "And the grid shows only the columns \"USUBJID, AGE, SEX, RACE, HEIGHT, WEIGHT\"", () => gridShowsColumns(page, "USUBJID, AGE, SEX, RACE, HEIGHT, WEIGHT"));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "Then no rows should be selected", () => noneSelected(page));
    await session.step(22, "When user picks \"Select > Missing Values...\" from the top menu", () => pickFromTopMenu(page, "Select > Missing Values..."), undefined, "Select > Missing Values finds the empty cells");
    await session.step(23, "Then Select Missing Values dialog should be visible", () => shouldBe(page, el("Select Missing Values dialog"), "visible"));
    await session.step(24, "When user toggles the \"HEIGHT\" column in the column list of Select Missing Values dialog", () => toggleInColumnList(page, "HEIGHT", el("Select Missing Values dialog")));
    await session.step(25, "And user clicks on OK button in Select Missing Values dialog", () => clickOn(page, el("OK button in Select Missing Values dialog")));
    await session.step(26, "Then 751 rows should be selected", () => selectedRowCount(page, 751));
    await session.step(27, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(28, "When user hovers over \"Selected: 751\" text in status bar", () => hoverOver(page, el("\"Selected: 751\" text in status bar")));
    await session.step(29, "And the video says \"751 rows have no HEIGHT; they are selected\"", () => videoSays(page, "751 rows have no HEIGHT; they are selected"));
    await session.step(31, "When user picks \"Select > Selection to Filter\" from the top menu", () => pickFromTopMenu(page, "Select > Selection to Filter"), undefined, "Selection to Filter shows only the selected rows");
    await session.step(32, "Then 751 rows should pass the filter", () => filterPasses(page, 751));
    await session.step(33, "And every selected row should pass the filter", () => selectedPassFilter(page));
    await session.step(35, "When user picks \"Select > Filter to Column...\" from the top menu", () => pickFromTopMenu(page, "Select > Filter to Column..."), undefined, "Filter to Column keeps the filter as a column");
    await session.step(36, "Then Name input in \"Filter to Column\" dialog should have the value \"Filter\"", () => shouldHaveValue(page, el("Name input in \"Filter to Column\" dialog"), "Filter"));
    await session.step(37, "When user enters \"No HEIGHT\" into Name input in \"Filter to Column\" dialog", () => enterInto(page, "No HEIGHT", el("Name input in \"Filter to Column\" dialog")));
    await session.step(38, "And user clicks on OK button in \"Filter to Column\" dialog", () => clickOn(page, el("OK button in \"Filter to Column\" dialog")));
    await session.step(39, "Then the \"Filter to Column\" dialog should close", () => dialogCloses(page, "Filter to Column"));
    await session.step(40, "And \"No HEIGHT\" column should have type \"bool\"", () => columnType(page, "No HEIGHT", "bool"));
    await session.step(41, "And the filter should pass exactly the rows where \"No HEIGHT\" is \"true\"", () => filterIsExactlyCategory(page, "No HEIGHT", "true"));
    await session.step(42, "When the video says \"The new column is true in exactly these 751 rows\"", () => videoSays(page, "The new column is true in exactly these 751 rows"));
    await session.step(43, "Then no errors should have been logged", () => noErrors(page));
  });
});
