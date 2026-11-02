/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/shortcuts.feature
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
import {clearField, enterInto, pressKey, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsSelected, currentColumnIs, noColumnsSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeConsole, closeCurrentView, currentViewType, openTableOf, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, clickAreaHolding, dragAreaToArea, noErrors, pickFromAreaContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The table view's keyboard shortcuts toggle the shell's panels, and Alt+C's Column Manager acts on the grid", () => {
  const session = feature(test, "features/table-view/shortcuts.feature", import.meta.url);
  test("The table view's keyboard shortcuts toggle the shell's panels, and Alt+C's Column Manager acts on the grid", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens a table \"bdd-colman\" with:", () => openTableOf(page, "bdd-colman", [["name","num","price","date"],["a","1","1.5","2020-01-01"],["b","2","2.5","2020-01-02"],["c","3","3.5","2020-01-03"]]), [["name","num","price","date"],["a","1","1.5","2020-01-01"],["b","2","2.5","2020-01-02"],["c","3","3.5","2020-01-03"]]);
    await run.scenario("Alt+X, Alt+V and the tilde toggle the Toolbox, the Variables panel and the Console", async () => {
      await session.step(32, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(33, "Then toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
      await session.step(34, "When user presses Alt+X", () => pressKey(page, "Alt+X"));
      await session.step(35, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
      await session.step(36, "And \"Toolbox\" status bar toggle should not be selected", () => shouldNotBe(page, el("\"Toolbox\" status bar toggle"), "selected"));
      await session.step(37, "When user presses Alt+X", () => pressKey(page, "Alt+X"));
      await session.step(38, "Then toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
      await session.step(39, "And \"Toolbox\" status bar toggle should be selected", () => shouldBe(page, el("\"Toolbox\" status bar toggle"), "selected"));
      await session.step(40, "And \"Variables\" dock panel should be absent", () => shouldBe(page, el("\"Variables\" dock panel"), "absent"));
      await session.step(41, "When user presses Alt+V", () => pressKey(page, "Alt+V"));
      await session.step(42, "Then \"Variables\" dock panel should be visible", () => shouldBe(page, el("\"Variables\" dock panel"), "visible"));
      await session.step(43, "And \"Variables\" status bar toggle should be selected", () => shouldBe(page, el("\"Variables\" status bar toggle"), "selected"));
      await session.step(44, "When user presses Alt+V", () => pressKey(page, "Alt+V"));
      await session.step(45, "Then \"Variables\" dock panel should be absent", () => shouldBe(page, el("\"Variables\" dock panel"), "absent"));
      await session.step(46, "And console should be absent", () => shouldBe(page, el("console"), "absent"));
      await session.step(47, "When user presses Backquote", () => pressKey(page, "Backquote"));
      await session.step(48, "Then console should be visible", () => shouldBe(page, el("console"), "visible"));
      await session.step(49, "And \"Console\" status bar toggle should be selected", () => shouldBe(page, el("\"Console\" status bar toggle"), "selected"));
      await session.step(50, "When user closes the console", () => closeConsole(page));
      await session.step(51, "Then console should be hidden", () => shouldBe(page, el("console"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Alt+S opens the Settings", async () => {
      await session.step(55, "When user presses Alt+S", () => pressKey(page, "Alt+S"));
      await session.step(56, "Then the current view should be a settings view", () => currentViewType(page, "settings"));
      await session.step(57, "And the \"Settings - General\" view should be current", () => viewIsCurrent(page, "Settings - General"));
      await session.step(58, "When user closes the current view", () => closeCurrentView(page));
      await session.step(59, "Then the \"bdd-colman\" view should be current", () => viewIsCurrent(page, "bdd-colman"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+Shift+V opens a second view on the same table", async () => {
      await session.step(63, "When user presses Control+Shift+V", () => pressKey(page, "Control+Shift+V"));
      await session.step(64, "Then the \"bdd-colman copy\" view should be current", () => viewIsCurrent(page, "bdd-colman copy"));
      await session.step(65, "And the open table views should be exactly \"bdd-colman, bdd-colman\"", () => openTableViewsExactly(page, "bdd-colman, bdd-colman"));
      await session.step(66, "When user closes the current view", () => closeCurrentView(page));
      await session.step(67, "Then the open table views should be exactly \"bdd-colman\"", () => openTableViewsExactly(page, "bdd-colman"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Alt+C docks the Columns pane with a row per column, a click on a row makes that column current, and Alt+C closes it", async () => {
      await session.step(71, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(72, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(73, "Then the current column should be \"name\"", () => currentColumnIs(page, "name"));
      await session.step(74, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(75, "Then column manager should be visible", () => shouldBe(page, el("column manager"), "visible"));
      await session.step(76, "And the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(77, "And the \"text of cell 1 of __name\" reading of Grid viewer in column manager should be \"name\"", () => readingReads(page, "text of cell 1 of __name", el("Grid viewer in column manager"), "name"));
      await session.step(78, "And the \"text of cell 4 of __name\" reading of Grid viewer in column manager should be \"date\"", () => readingReads(page, "text of cell 4 of __name", el("Grid viewer in column manager"), "date"));
      await session.step(79, "And the \"text of cell 1 of __t\" reading of Grid viewer in column manager should be \"s\"", () => readingReads(page, "text of cell 1 of __t", el("Grid viewer in column manager"), "s"));
      await session.step(80, "And the \"text of cell 2 of __t\" reading of Grid viewer in column manager should be \"i\"", () => readingReads(page, "text of cell 2 of __t", el("Grid viewer in column manager"), "i"));
      await session.step(81, "And the \"text of cell 3 of __t\" reading of Grid viewer in column manager should be \"d\"", () => readingReads(page, "text of cell 3 of __t", el("Grid viewer in column manager"), "d"));
      await session.step(82, "When user clicks on the \"cell 3 of __name\" area of Grid viewer in column manager", () => clickArea(page, "cell 3 of __name", el("Grid viewer in column manager")));
      await session.step(83, "Then the current column should be \"price\"", () => currentColumnIs(page, "price"));
      await session.step(84, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(85, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The search box narrows the list to the matching columns", async () => {
      await session.step(89, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(90, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(91, "Then the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(92, "When user enters \"pr\" into \"Search\" input in column manager", () => enterInto(page, "pr", el("\"Search\" input in column manager")));
      await session.step(93, "Then the \"rows shown\" reading of Grid viewer in column manager should be 1", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 1));
      await session.step(94, "And the \"text of cell 3 of __name\" reading of Grid viewer in column manager should be \"price\"", () => readingReads(page, "text of cell 3 of __name", el("Grid viewer in column manager"), "price"));
      await session.step(95, "When user enters \"a\" into \"Search\" input in column manager", () => enterInto(page, "a", el("\"Search\" input in column manager")));
      await session.step(96, "Then the \"rows shown\" reading of Grid viewer in column manager should be 2", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 2));
      await session.step(97, "When user clears \"Search\" input in column manager", () => clearField(page, el("\"Search\" input in column manager")));
      await session.step(98, "Then the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(99, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(100, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+click toggles a column's selection in the grid, and Esc clears it", async () => {
      await session.step(104, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(105, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(106, "Then no columns should be selected", () => noColumnsSelected(page));
      await session.step(107, "When user clicks on the \"cell 2 of __name\" area of Grid viewer in column manager holding Control", () => clickAreaHolding(page, "cell 2 of __name", el("Grid viewer in column manager"), "Control"));
      await session.step(108, "And user clicks on the \"cell 4 of __name\" area of Grid viewer in column manager holding Control", () => clickAreaHolding(page, "cell 4 of __name", el("Grid viewer in column manager"), "Control"));
      await session.step(109, "Then columns \"num, date\" should be selected", () => columnsSelected(page, "num, date"));
      await session.step(110, "When user clicks on the \"cell 4 of __name\" area of Grid viewer in column manager holding Control", () => clickAreaHolding(page, "cell 4 of __name", el("Grid viewer in column manager"), "Control"));
      await session.step(111, "Then columns \"num\" should be selected", () => columnsSelected(page, "num"));
      await session.step(112, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(113, "Then no columns should be selected", () => noColumnsSelected(page));
      await session.step(114, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(115, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Dragging a column's row reorders the columns of the grid", async () => {
      await session.step(119, "Then the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(120, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(121, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(122, "And user drags the \"cell 4 of __name\" area of Grid viewer in column manager to the \"cell 1 of __name\" area", () => dragAreaToArea(page, "cell 4 of __name", el("Grid viewer in column manager"), "cell 1 of __name"));
      await session.step(123, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(124, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(125, "And the \"column order\" reading of grid should be \"date, name, num, price\"", () => readingReads(page, "column order", el("grid"), "date, name, num, price"));
      await session.step(126, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(127, "Then the \"text of cell 1 of __name\" reading of Grid viewer in column manager should be \"date\"", () => readingReads(page, "text of cell 1 of __name", el("Grid viewer in column manager"), "date"));
      await session.step(128, "When user drags the \"cell 1 of __name\" area of Grid viewer in column manager to the \"cell 4 of __name\" area", () => dragAreaToArea(page, "cell 1 of __name", el("Grid viewer in column manager"), "cell 4 of __name"));
      await session.step(129, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(130, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(131, "And the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hide in the popup menu takes the column out of the grid, and Show brings it back", async () => {
      await session.step(135, "Then the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(136, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(137, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(138, "And user picks \"'num' column > Hide\" from the context menu of the \"cell 2 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "'num' column > Hide", "cell 2 of __name", el("Grid viewer in column manager")));
      await session.step(139, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(140, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(141, "And the \"column order\" reading of grid should be \"name, price, date\"", () => readingReads(page, "column order", el("grid"), "name, price, date"));
      await session.step(142, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(143, "Then the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(144, "When user picks \"'num' column > Show\" from the context menu of the \"cell 2 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "'num' column > Show", "cell 2 of __name", el("Grid viewer in column manager")));
      await session.step(145, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(146, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(147, "And the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(148, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Statistics in the popup menu adds a statistics column to the list", async () => {
      await session.step(151, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(152, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(153, "Then the \"column order\" reading of Grid viewer in column manager should be \"__t, __name\"", () => readingReads(page, "column order", el("Grid viewer in column manager"), "__t, __name"));
      await session.step(154, "When user picks \"Statistics > count\" from the context menu of the \"cell 1 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "Statistics > count", "cell 1 of __name", el("Grid viewer in column manager")));
      await session.step(155, "Then the \"column order\" reading of Grid viewer in column manager should be \"__t, __name, count\"", () => readingReads(page, "column order", el("Grid viewer in column manager"), "__t, __name, count"));
      await session.step(156, "And the \"text of cell 2 of count\" reading of Grid viewer in column manager should be \"3.00\"", () => readingReads(page, "text of cell 2 of count", el("Grid viewer in column manager"), "3.00"));
      await session.step(157, "When user picks \"Statistics > count\" from the context menu of the \"cell 1 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "Statistics > count", "cell 1 of __name", el("Grid viewer in column manager")));
      await session.step(158, "Then the \"column order\" reading of Grid viewer in column manager should be \"__t, __name\"", () => readingReads(page, "column order", el("Grid viewer in column manager"), "__t, __name"));
      await session.step(159, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(160, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(161, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
