/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/10-remove-duplicate-rows.feature
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
import {clickOn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {noneSelected, rowCount, rowsRangeAllSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickAreaHolding, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/10-remove-duplicate-rows.feature", import.meta.url);
  test("Find and remove duplicate rows", {tag: ["@guide", "@tutorial", "@viewport:1280x720"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens iris dataset", () => openDataset(page, ds("iris")));
    await session.step(16, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(19, "Then no rows should be selected", () => noneSelected(page));
    await session.step(21, "When user picks \"Select > Duplicates...\" from the top menu", () => pickFromTopMenu(page, "Select > Duplicates..."), undefined, "Select > Duplicates finds rows that repeat");
    await session.step(22, "Then Select Duplicates dialog should be visible", () => shouldBe(page, el("Select Duplicates dialog"), "visible"));
    await session.step(24, "When user toggles the \"Sepal.Length\" column in the column list of Select Duplicates dialog", () => toggleInColumnList(page, "Sepal.Length", el("Select Duplicates dialog")), undefined, "Tick the columns that must match");
    await session.step(25, "And user toggles the \"Sepal.Width\" column in the column list of Select Duplicates dialog", () => toggleInColumnList(page, "Sepal.Width", el("Select Duplicates dialog")));
    await session.step(26, "And user toggles the \"Petal.Length\" column in the column list of Select Duplicates dialog", () => toggleInColumnList(page, "Petal.Length", el("Select Duplicates dialog")));
    await session.step(27, "And user toggles the \"Petal.Width\" column in the column list of Select Duplicates dialog", () => toggleInColumnList(page, "Petal.Width", el("Select Duplicates dialog")));
    await session.step(28, "And user toggles the \"Species\" column in the column list of Select Duplicates dialog", () => toggleInColumnList(page, "Species", el("Select Duplicates dialog")));
    await session.step(29, "And user clicks on OK button in Select Duplicates dialog", () => clickOn(page, el("OK button in Select Duplicates dialog")));
    await session.step(30, "Then Select Duplicates dialog should be hidden", () => shouldBe(page, el("Select Duplicates dialog"), "hidden"));
    await session.step(31, "And 2 rows should be selected", () => selectedRowCount(page, 2));
    await session.step(32, "And rows 102 to 102 should all be selected", () => rowsRangeAllSelected(page, 102, 102));
    await session.step(33, "And rows 143 to 143 should all be selected", () => rowsRangeAllSelected(page, 143, 143));
    await session.step(35, "When user picks \"Select > Selection to Top\" from the top menu", () => pickFromTopMenu(page, "Select > Selection to Top"), undefined, "Selection to Top brings them to the top of the grid");
    await session.step(36, "And the video says \"The same flower, entered twice\"", () => videoSays(page, "The same flower, entered twice"));
    await session.step(38, "And user clicks on the \"cell 102 of Species\" area of grid holding Control", () => clickAreaHolding(page, "cell 102 of Species", el("grid"), "Control"), undefined, "Ctrl+click one copy to keep it");
    await session.step(39, "Then 1 row should be selected", () => selectedRowCount(page, 1));
    await session.step(41, "When user picks \"Edit > Remove > Selected Rows\" from the top menu", () => pickFromTopMenu(page, "Edit > Remove > Selected Rows"), undefined, "Then remove the selected rows");
    await session.step(42, "Then the table should have 149 rows", () => rowCount(page, 149));
    await session.step(43, "And no rows should be selected", () => noneSelected(page));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
});
