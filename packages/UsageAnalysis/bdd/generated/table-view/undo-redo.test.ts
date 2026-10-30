/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/undo-redo.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, hasColumn, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCodedAs, noColorCoding, tableColumns} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, openTableOf, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, noErrors, pickFromAreaContextMenu, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Ctrl+Z undoes and Ctrl+Shift+Z or Ctrl+Y redoes the table view's column and layout operations", () => {
  const session = feature(test, "features/table-view/undo-redo.feature", import.meta.url);
  test("Ctrl+Z undoes and Ctrl+Shift+Z or Ctrl+Y redoes the table view's column and layout operations", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens a table \"undo\" with:", () => openTableOf(page, "undo", [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cy","Oslo","3"]]), [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cy","Oslo","3"]]);
    await run.scenario("A column removed from its header menu comes back in its place, and redo removes it again", async () => {
      await session.step(28, "When user picks \"Remove\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header city", el("grid")));
      await session.step(29, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(30, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(31, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(32, "Then the table should have a column \"city\"", () => hasColumn(page, "city"));
      await session.step(33, "And the \"column order\" reading of grid should be \"name, city, weight\"", () => readingReads(page, "column order", el("grid"), "name, city, weight"));
      await session.step(34, "And the value of \"city\" column in row 2 should be \"Rome\"", () => valueInRow(page, "city", 2, "Rome"));
      await session.step(35, "When user presses Control+Shift+Z", () => pressKey(page, "Control+Shift+Z"));
      await session.step(36, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(37, "And the table should have 2 columns", () => columnCount(page, 2));
      await session.step(38, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(39, "Then the \"column order\" reading of grid should be \"name, city, weight\"", () => readingReads(page, "column order", el("grid"), "name, city, weight"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column renamed in Column Properties gets its name back, and Ctrl+Y renames it again", async () => {
      await session.step(43, "When user picks \"Column Properties...\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header city", el("grid")));
      await session.step(44, "And user enters \"town\" into \"New name:\" input in \"city\" dialog", () => enterInto(page, "town", el("\"New name:\" input in \"city\" dialog")));
      await session.step(45, "And user clicks on OK button in \"city\" dialog", () => clickOn(page, el("OK button in \"city\" dialog")));
      await session.step(46, "Then the table should have a column \"town\"", () => hasColumn(page, "town"));
      await session.step(47, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(48, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(49, "Then the table should have a column \"city\"", () => hasColumn(page, "city"));
      await session.step(50, "And the table should not have a column \"town\"", () => hasNoColumn(page, "town"));
      await session.step(51, "And the value of \"city\" column in row 2 should be \"Rome\"", () => valueInRow(page, "city", 2, "Rome"));
      await session.step(52, "When user presses Control+Y", () => pressKey(page, "Control+Y"));
      await session.step(53, "Then the table should have a column \"town\"", () => hasColumn(page, "town"));
      await session.step(54, "And the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(55, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(56, "Then the table should have a column \"city\"", () => hasColumn(page, "city"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color coding set from the header menu is taken off by Ctrl+Z and put back by redo", async () => {
      await session.step(60, "When user picks \"Color Coding > Linear\" from the context menu of the \"header weight\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header weight", el("grid")));
      await session.step(61, "Then \"weight\" column should be color-coded linearly", () => colorCodedAs(page, "weight", "linearly"));
      await session.step(62, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(63, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(64, "Then \"weight\" column should have no color coding", () => noColorCoding(page, "weight"));
      await session.step(65, "When user presses Control+Shift+Z", () => pressKey(page, "Control+Shift+Z"));
      await session.step(66, "Then \"weight\" column should be color-coded linearly", () => colorCodedAs(page, "weight", "linearly"));
      await session.step(67, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(68, "Then \"weight\" column should have no color coding", () => noColorCoding(page, "weight"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("View > Layout > Clear removes the viewers, Ctrl+Z brings them back and redo clears again", async () => {
      await session.step(72, "Given user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(73, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(74, "When user picks \"View > Layout > Clear\" from the top menu", () => pickFromTopMenu(page, "View > Layout > Clear"));
      await session.step(75, "Then the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
      await session.step(76, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(77, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(78, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(79, "Then the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
      await session.step(80, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(81, "When user presses Control+Shift+Z", () => pressKey(page, "Control+Shift+Z"));
      await session.step(82, "Then the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
      await session.step(83, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit > Undo and Edit > Redo do what the shortcuts do", async () => {
      await session.step(87, "When user picks \"Remove\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header city", el("grid")));
      await session.step(88, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(89, "When user picks \"Edit > Undo\" from the top menu", () => pickFromTopMenu(page, "Edit > Undo"));
      await session.step(90, "Then the table should have a column \"city\"", () => hasColumn(page, "city"));
      await session.step(91, "When user picks \"Edit > Redo\" from the top menu", () => pickFromTopMenu(page, "Edit > Redo"));
      await session.step(92, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(93, "When user picks \"Edit > Undo\" from the top menu", () => pickFromTopMenu(page, "Edit > Undo"));
      await session.step(94, "Then the \"column order\" reading of grid should be \"name, city, weight\"", () => readingReads(page, "column order", el("grid"), "name, city, weight"));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("While a dialog is open Ctrl+Z leaves the table alone, and after it closes Ctrl+Z works", async () => {
      await session.step(98, "When user picks \"Remove\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header city", el("grid")));
      await session.step(99, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(100, "When user clicks on the \"cell 1 of weight\" area of grid", () => clickArea(page, "cell 1 of weight", el("grid")));
      await session.step(101, "And user presses F2", () => pressKey(page, "F2"));
      await session.step(102, "Then \"weight\" dialog should be visible", () => shouldBe(page, el("\"weight\" dialog"), "visible"));
      await session.step(103, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(104, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(105, "And \"weight\" dialog should be visible", () => shouldBe(page, el("\"weight\" dialog"), "visible"));
      await session.step(106, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(107, "Then \"weight\" dialog should be hidden", () => shouldBe(page, el("\"weight\" dialog"), "hidden"));
      await session.step(108, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(109, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(110, "Then the table should have a column \"city\"", () => hasColumn(page, "city"));
      await session.step(111, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+Z in another table's view does not undo this table's change, and switching back brings the history", async () => {
      await session.step(114, "When user picks \"Remove\" from the context menu of the \"header city\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header city", el("grid")));
      await session.step(115, "Then the table should not have a column \"city\"", () => hasNoColumn(page, "city"));
      await session.step(116, "Given user opens a table \"other\" with:", () => openTableOf(page, "other", [["id","score"],["1","10"],["2","20"]]), [["id","score"],["1","10"],["2","20"]]);
      await session.step(120, "When user clicks on the \"cell 1 of id\" area of grid", () => clickArea(page, "cell 1 of id", el("grid")));
      await session.step(121, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(122, "Then table \"undo\" should have columns \"name, weight\"", () => tableColumns(page, "undo", "name, weight"));
      await session.step(123, "And table \"other\" should have columns \"id, score\"", () => tableColumns(page, "other", "id, score"));
      await session.step(124, "When user switches to the \"undo\" table view", () => switchTableView(page, "undo"));
      await session.step(125, "And user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(126, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(127, "Then the table should have a column \"city\"", () => hasColumn(page, "city"));
      await session.step(128, "And table \"other\" should have columns \"id, score\"", () => tableColumns(page, "other", "id, score"));
      await session.step(129, "When user switches to the \"other\" table view", () => switchTableView(page, "other"));
      await session.step(130, "And user closes the current view", () => closeCurrentView(page));
      await session.step(131, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
