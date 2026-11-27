/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/06-edit-a-table.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
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
import {closeVariablesPane} from '../../../bindings/howto-data.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, hoverOver, pressKey, pressKeyIn, selectIn, shouldBe, shouldHaveValue, typeInto, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {displayedInRow, hasColumn, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterPassesAll, rowCount, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, doubleClickArea, hoverArea, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/06-edit-a-table.feature", import.meta.url);
  test("Edit cells and rows, find and replace, edit many rows at once", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(25, "When user double-clicks on the \"cell 1 of price\" area of grid", () => doubleClickArea(page, "cell 1 of price", el("grid")), undefined, "Double-click a cell to edit it");
    await session.step(26, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
    await session.step(27, "And user types \"15900\" into cell editor", () => typeInto(page, "15900", el("cell editor")));
    await session.step(28, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(29, "Then the \"price\" cell of row 1 should be displayed as \"15900\"", () => displayedInRow(page, "price", 1, "15900"));
    await session.step(31, "When user picks \"Edit > Add Rows...\" from the top menu", () => pickFromTopMenu(page, "Edit > Add Rows..."), undefined, "Edit > Add Rows inserts empty rows at the current one");
    await session.step(32, "Then Add Rows dialog should be visible", () => shouldBe(page, el("Add Rows dialog"), "visible"));
    await session.step(33, "When user enters \"2\" into Count input in Add Rows dialog", () => enterInto(page, "2", el("Count input in Add Rows dialog")));
    await session.step(34, "And user clicks on OK button in Add Rows dialog", () => clickOn(page, el("OK button in Add Rows dialog")));
    await session.step(35, "Then Add Rows dialog should be hidden", () => shouldBe(page, el("Add Rows dialog"), "hidden"));
    await session.step(36, "And the table should have 32 rows", () => rowCount(page, 32));
    await session.step(37, "And the value of \"model\" column in row 1 should be \"\"", () => valueInRow(page, "model", 1, ""));
    await session.step(39, "When user clicks on the \"row header 1\" area of grid", () => clickArea(page, "row header 1", el("grid")), undefined, "Select rows by their headers; Shift+Delete removes them");
    await session.step(40, "And user clicks on the \"row header 2\" area of grid holding Shift", () => clickAreaHolding(page, "row header 2", el("grid"), "Shift"));
    await session.step(41, "Then 2 rows should be selected", () => selectedRowCount(page, 2));
    await session.step(42, "When user presses Shift+Delete", () => pressKey(page, "Shift+Delete"));
    await session.step(43, "Then the table should have 30 rows", () => rowCount(page, 30));
    await session.step(44, "And the value of \"model\" column in row 1 should be \"alfaromeo\"", () => valueInRow(page, "model", 1, "alfaromeo"));
    await session.step(46, "When user picks \"Remove\" from the context menu of the \"header symbol\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header symbol", el("grid")), undefined, "A column goes from its header's menu");
    await session.step(47, "Then the table should not have a column \"symbol\"", () => hasNoColumn(page, "symbol"));
    await session.step(49, "When user clicks on the \"cell 3 of model\" area of grid", () => clickArea(page, "cell 3 of model", el("grid")), undefined, "Ctrl+Z undoes the last change: symbol is back");
    await session.step(50, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
    await session.step(51, "Then the table should have a column \"symbol\"", () => hasColumn(page, "symbol"));
    await session.step(52, "When user hovers over the \"header symbol\" area of grid", () => hoverArea(page, "header symbol", el("grid")));
    await session.step(54, "When user presses Control+Y", () => pressKey(page, "Control+Y"), undefined, "Ctrl+Y does it again (both are in the Edit menu too)");
    await session.step(55, "Then the table should not have a column \"symbol\"", () => hasNoColumn(page, "symbol"));
    await session.step(57, "When user clicks on the \"cell 22 of model\" area of grid", () => clickArea(page, "cell 22 of model", el("grid")), undefined, "Ctrl+H opens Find and Replace");
    await session.step(58, "And user presses Control+H", () => pressKey(page, "Control+H"));
    await session.step(59, "Then Find and Replace dialog should be visible", () => shouldBe(page, el("Find and Replace dialog"), "visible"));
    await session.step(61, "When user enters \"toyota\" into \"Find what\" input in Find and Replace dialog", () => enterInto(page, "toyota", el("\"Find what\" input in Find and Replace dialog")), undefined, "What to find, and what to put instead");
    await session.step(62, "And user enters \"Toyota \" into \"Replace with\" input in Find and Replace dialog", () => enterInto(page, "Toyota ", el("\"Replace with\" input in Find and Replace dialog")));
    await session.step(64, "And user selects \"All Columns\" in \"Where columns\" input in Find and Replace dialog", () => selectIn(page, "All Columns", el("\"Where columns\" input in Find and Replace dialog")), undefined, "Search every column, or only some");
    await session.step(66, "And user checks \"Filter matching rows\" input in Find and Replace dialog", () => check(page, el("\"Filter matching rows\" input in Find and Replace dialog")), undefined, "Filter matching rows shows what will change");
    await session.step(67, "Then 4 rows should pass the filter", () => filterPasses(page, 4));
    await session.step(69, "When user clicks on OK button in Find and Replace dialog", () => clickOn(page, el("OK button in Find and Replace dialog")), undefined, "OK replaces every match");
    await session.step(70, "Then Find and Replace dialog should be hidden", () => shouldBe(page, el("Find and Replace dialog"), "hidden"));
    await session.step(71, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(72, "And the value of \"model\" column in row 22 should be \"Toyota 1\"", () => valueInRow(page, "model", 22, "Toyota 1"));
    await session.step(73, "And the value of \"model\" column in row 25 should be \"Toyota 4\"", () => valueInRow(page, "model", 25, "Toyota 4"));
    await session.step(75, "Given the variables pane is closed", () => closeVariablesPane(page));
    await session.step(77, "When user clicks on the \"row header 22\" area of grid", () => clickArea(page, "row header 22", el("grid")), undefined, "To set one value in many rows, select them");
    await session.step(78, "And user clicks on the \"row header 25\" area of grid holding Shift", () => clickAreaHolding(page, "row header 25", el("grid"), "Shift"));
    await session.step(79, "Then 4 rows should be selected", () => selectedRowCount(page, 4));
    await session.step(81, "When user picks \"Data > Batch Edit...\" from the top menu", () => pickFromTopMenu(page, "Data > Batch Edit..."), undefined, "Batch Edit: edit the selected rows together");
    await session.step(82, "Then \"Rows\" input should have value \"Selected\"", () => shouldHaveValue(page, el("\"Rows\" input"), "Selected"));
    await session.step(84, "When user enters \"1\" into \"turbo\" property", () => enterInto(page, "1", el("\"turbo\" property")), undefined, "A value typed once goes into every selected row");
    await session.step(85, "Then the \"turbo\" cell of row 22 should be displayed as \"1\"", () => displayedInRow(page, "turbo", 22, "1"));
    await session.step(86, "And the \"turbo\" cell of row 25 should be displayed as \"1\"", () => displayedInRow(page, "turbo", 25, "1"));
    await session.step(87, "And the \"turbo\" cell of row 26 should be displayed as \"0\"", () => displayedInRow(page, "turbo", 26, "0"));
    await session.step(88, "When the video says \"All four Toyotas now read 1 in turbo\" at the bottom", () => videoSaysAt(page, "All four Toyotas now read 1 in turbo", "bottom"));
    await session.step(89, "And user hovers over \"Batch Edit\" dock panel", () => hoverOver(page, el("\"Batch Edit\" dock panel")));
    await session.step(90, "And user clicks on the close icon of \"Batch Edit\" dock panel", () => clickOn(page, el("the close icon of \"Batch Edit\" dock panel")));
    await session.step(91, "Then \"Batch Edit\" dock panel should be absent", () => shouldBe(page, el("\"Batch Edit\" dock panel"), "absent"));
    await session.step(92, "And no errors should have been logged", () => noErrors(page));
  });
});
