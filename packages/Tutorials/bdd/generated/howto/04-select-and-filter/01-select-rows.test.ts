/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/01-select-rows.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {hoverOver, pressKeyIn, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {noneSelected, rowsRangeAllSelected, selectedPassFilter, selectedRowCount, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, clickAreaHolding, dragSelectionOverArea, noErrors, readingHigherThanRemembered, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/01-select-rows.feature", import.meta.url);
  test("Select rows", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "Then no rows should be selected", () => noneSelected(page));
    await session.step(24, "When user clicks on the \"cell 3 of AGE\" area of grid", () => clickArea(page, "cell 3 of AGE", el("grid")), undefined, "A click makes a row current, it does not select it");
    await session.step(25, "Then row 3 should be current", () => currentRowIs(page, 3));
    await session.step(26, "And no rows should be selected", () => noneSelected(page));
    await session.step(28, "When user clicks on the \"cell 8 of AGE\" area of grid holding Shift", () => clickAreaHolding(page, "cell 8 of AGE", el("grid"), "Shift"), undefined, "Shift+click selects every row in between");
    await session.step(29, "Then 6 rows should be selected", () => selectedRowCount(page, 6));
    await session.step(30, "And rows 3 to 8 should all be selected", () => rowsRangeAllSelected(page, 3, 8));
    await session.step(32, "When user clicks on the \"cell 11 of AGE\" area of grid holding Control", () => clickAreaHolding(page, "cell 11 of AGE", el("grid"), "Control"), undefined, "Ctrl+click adds a row, or takes it out");
    await session.step(33, "Then 7 rows should be selected", () => selectedRowCount(page, 7));
    await session.step(35, "When user hovers over \"Selected: 7\" text in status bar", () => hoverOver(page, el("\"Selected: 7\" text in status bar")), undefined, "The status bar counts the selected rows");
    await session.step(37, "When user picks \"Select > Invert\" from the top menu", () => pickFromTopMenu(page, "Select > Invert"), undefined, "The Select menu: all, none, invert and more");
    await session.step(38, "Then 5843 rows should be selected", () => selectedRowCount(page, 5843));
    await session.step(40, "When user presses Escape in grid", () => pressKeyIn(page, "Escape", el("grid")), undefined, "Esc clears the selection");
    await session.step(41, "Then no rows should be selected", () => noneSelected(page));
    await session.step(43, "When user remembers the \"rows selected\" reading of scatter plot viewer", () => rememberReading(page, "rows selected", el("scatter plot viewer")), undefined, "On a plot, Shift+drag a box around the points");
    await session.step(44, "And user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
    await session.step(45, "Then the \"rows selected\" reading of scatter plot viewer should be higher than remembered", () => readingHigherThanRemembered(page, "rows selected", el("scatter plot viewer")));
    await session.step(46, "And some rows should be selected", () => someSelected(page));
    await session.step(47, "When the video says \"The same rows are selected in the grid and in every viewer\"", () => videoSays(page, "The same rows are selected in the grid and in every viewer"));
    await session.step(48, "Then every selected row should pass the filter", () => selectedPassFilter(page));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
