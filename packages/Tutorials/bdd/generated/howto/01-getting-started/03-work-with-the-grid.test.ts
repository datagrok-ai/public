/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/01-getting-started/03-work-with-the-grid.feature
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
import {clickOn, hoverOver, pressKey, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {currentColumnIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {clickViewTab, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaTaller, clickArea, doubleClickArea, dragAreaBy, dragAreaToArea, hasArea, hasNoArea, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingHigher, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingIncludes, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Basics", () => {
  const session = feature(test, "features/howto/01-getting-started/03-work-with-the-grid.feature", import.meta.url);
  test("Work with the grid", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(20, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "When user clicks on the \"cell 1 of AGE\" area of grid", () => clickArea(page, "cell 1 of AGE", el("grid")), undefined, "Click a cell; the arrow keys move from there");
    await session.step(25, "And user presses ArrowDown", () => pressKey(page, "ArrowDown"));
    await session.step(26, "And user presses ArrowDown", () => pressKey(page, "ArrowDown"));
    await session.step(27, "And user presses ArrowRight", () => pressKey(page, "ArrowRight"));
    await session.step(28, "Then row 3 should be current", () => currentRowIs(page, 3));
    await session.step(30, "When user presses Control+End", () => pressKey(page, "Control+End"), undefined, "Ctrl+End jumps to the last row, Ctrl+Home back to the first");
    await session.step(31, "Then row 5850 should be current", () => currentRowIs(page, 5850));
    await session.step(32, "When user presses Control+Home", () => pressKey(page, "Control+Home"));
    await session.step(33, "Then row 1 should be current", () => currentRowIs(page, 1));
    await session.step(35, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")), undefined, "Double-click a header to sort by it");
    await session.step(36, "Then the \"sort column\" reading of grid should be \"AGE\"", () => readingReads(page, "sort column", el("grid"), "AGE"));
    await session.step(37, "And the \"sort direction\" reading of grid should be \"descending\"", () => readingReads(page, "sort direction", el("grid"), "descending"));
    await session.step(39, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")), undefined, "Again for the other direction, a third time to unsort");
    await session.step(40, "Then the \"sort direction\" reading of grid should be \"ascending\"", () => readingReads(page, "sort direction", el("grid"), "ascending"));
    await session.step(41, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
    await session.step(42, "Then the \"sort column\" reading of grid should be \"\"", () => readingReads(page, "sort column", el("grid"), ""));
    await session.step(44, "When user drags the \"column resizer RACE\" area of grid by 60 pixels to the right", () => dragAreaBy(page, "column resizer RACE", el("grid"), 60, "right"), undefined, "Drag a header's right edge to resize the column");
    await session.step(45, "Then the \"column width of RACE\" reading of grid should be higher than before", () => readingHigher(page, "column width of RACE", el("grid")));
    await session.step(47, "When user drags the \"header HEIGHT\" area of grid to the \"header AGE\" area", () => dragAreaToArea(page, "header HEIGHT", el("grid"), "header AGE"), undefined, "Drag a header to move the column");
    await session.step(48, "Then the \"column order\" reading of grid should include the text \"AGE, HEIGHT, SEX\"", () => readingIncludes(page, "column order", el("grid"), "AGE, HEIGHT, SEX"));
    await session.step(50, "When user picks \"Hide\" from the context menu of the \"header WEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Hide", "header WEIGHT", el("grid")), undefined, "Right-click a header to hide its column");
    await session.step(51, "Then grid should not have a \"header WEIGHT\" area", () => hasNoArea(page, el("grid"), "header WEIGHT"));
    await session.step(53, "When user picks \"Order or Hide Columns...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Order or Hide Columns...", "header AGE", el("grid")), undefined, "Order or Hide Columns lists them all: tick it to show it again");
    await session.step(54, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
    await session.step(55, "When user toggles the \"WEIGHT\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "WEIGHT", el("Order or Hide Columns dialog")));
    await session.step(56, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
    await session.step(57, "Then grid should have a \"header WEIGHT\" area", () => hasArea(page, el("grid"), "header WEIGHT"));
    await session.step(59, "When user drags the \"row resizer 3\" area of grid by 20 pixels to the down", () => dragAreaBy(page, "row resizer 3", el("grid"), 20, "down"), undefined, "Drag between row numbers to make rows taller");
    await session.step(60, "Then the \"row header 3\" area of grid should be taller than before", () => areaTaller(page, "row header 3", el("grid")));
    await session.step(62, "When user clicks on the tab of the \"spgi-100\" view", () => clickViewTab(page, "spgi-100"), undefined, "Now a wider table: click its tab");
    await session.step(63, "Then the \"spgi-100\" view should be current", () => viewIsCurrent(page, "spgi-100"));
    await session.step(65, "When user hovers over \"Columns: 88\" text in status bar", () => hoverOver(page, el("\"Columns: 88\" text in status bar")), undefined, "88 columns: more than fit on the screen");
    await session.step(67, "And user picks \"Pin > Pin Column\" from the context menu of the \"header Structure\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Column", "header Structure", el("grid")), undefined, "Pin a column to keep it in sight");
    await session.step(68, "Then \"Frozen Columns\" property of grid should be \"2\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "2"));
    await session.step(70, "When user clicks on \"Columns: 88\" text in status bar", () => clickOn(page, el("\"Columns: 88\" text in status bar")), undefined, "Columns in the status bar lists them all");
    await session.step(71, "Then column manager should be visible", () => shouldBe(page, el("column manager"), "visible"));
    await session.step(73, "When user types \"herg\" into \"Search\" input in column manager", () => typeInto(page, "herg", el("\"Search\" input in column manager")), undefined, "Type to find one by name");
    await session.step(74, "Then the \"rows shown\" reading of Grid viewer in column manager should be 1", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 1));
    await session.step(76, "When user clicks on the \"cell 36 of __name\" area of Grid viewer in column manager", () => clickArea(page, "cell 36 of __name", el("Grid viewer in column manager")), undefined, "Click it: the grid takes you there");
    await session.step(77, "Then the current column should be \"hERG Classification\"", () => currentColumnIs(page, "hERG Classification"));
    await session.step(78, "And no errors should have been logged", () => noErrors(page));
  });
});
