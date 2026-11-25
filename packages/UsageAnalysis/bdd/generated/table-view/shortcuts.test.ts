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
import {clearField, clickOn, enterInto, pressKey, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsSelected, currentColumnIs, noColumnsSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeConsole, closeCurrentView, currentViewType, openTableOf, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, clickAreaHolding, dragAreaToArea, noErrors, pickFromAreaContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The table view's keyboard shortcuts toggle the shell's panels, and Alt+C's Column Manager acts on the grid", () => {
  const session = feature(test, "features/table-view/shortcuts.feature", import.meta.url);
  test("The table view's keyboard shortcuts toggle the shell's panels, and Alt+C's Column Manager acts on the grid", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens a table \"bdd-colman\" with:", () => openTableOf(page, "bdd-colman", [["name","num","price","date"],["a","1","1.5","2020-01-01"],["b","2","2.5","2020-01-02"],["c","3","3.5","2020-01-03"]]), [["name","num","price","date"],["a","1","1.5","2020-01-01"],["b","2","2.5","2020-01-02"],["c","3","3.5","2020-01-03"]]);
    await run.scenario("Alt+X, Alt+V and the tilde toggle the Toolbox, the Variables panel and the Console", async () => {
      await session.step(33, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(34, "Then toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
      await session.step(35, "When user presses Alt+X", () => pressKey(page, "Alt+X"));
      await session.step(36, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
      await session.step(37, "And \"Toolbox\" status bar toggle should not be selected", () => shouldNotBe(page, el("\"Toolbox\" status bar toggle"), "selected"));
      await session.step(38, "When user presses Alt+X", () => pressKey(page, "Alt+X"));
      await session.step(39, "Then toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
      await session.step(40, "And \"Toolbox\" status bar toggle should be selected", () => shouldBe(page, el("\"Toolbox\" status bar toggle"), "selected"));
      await session.step(41, "And \"Variables\" dock panel should be absent", () => shouldBe(page, el("\"Variables\" dock panel"), "absent"));
      await session.step(42, "When user presses Alt+V", () => pressKey(page, "Alt+V"));
      await session.step(43, "Then \"Variables\" dock panel should be visible", () => shouldBe(page, el("\"Variables\" dock panel"), "visible"));
      await session.step(44, "And \"Variables\" status bar toggle should be selected", () => shouldBe(page, el("\"Variables\" status bar toggle"), "selected"));
      await session.step(45, "When user presses Alt+V", () => pressKey(page, "Alt+V"));
      await session.step(46, "Then \"Variables\" dock panel should be absent", () => shouldBe(page, el("\"Variables\" dock panel"), "absent"));
      await session.step(47, "And console should be absent", () => shouldBe(page, el("console"), "absent"));
      await session.step(48, "When user presses Backquote", () => pressKey(page, "Backquote"));
      await session.step(49, "Then console should be visible", () => shouldBe(page, el("console"), "visible"));
      await session.step(50, "And \"Console\" status bar toggle should be selected", () => shouldBe(page, el("\"Console\" status bar toggle"), "selected"));
      await session.step(51, "When user closes the console", () => closeConsole(page));
      await session.step(52, "Then console should be hidden", () => shouldBe(page, el("console"), "hidden"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Alt+S opens the Settings", async () => {
      await session.step(56, "When user presses Alt+S", () => pressKey(page, "Alt+S"));
      await session.step(57, "Then the current view should be a settings view", () => currentViewType(page, "settings"));
      await session.step(58, "And the \"Settings - General\" view should be current", () => viewIsCurrent(page, "Settings - General"));
      await session.step(59, "When user closes the current view", () => closeCurrentView(page));
      await session.step(60, "Then the \"bdd-colman\" view should be current", () => viewIsCurrent(page, "bdd-colman"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+Shift+V opens a second view on the same table", async () => {
      await session.step(64, "When user presses Control+Shift+V", () => pressKey(page, "Control+Shift+V"));
      await session.step(65, "Then the \"bdd-colman copy\" view should be current", () => viewIsCurrent(page, "bdd-colman copy"));
      await session.step(66, "And the open table views should be exactly \"bdd-colman, bdd-colman\"", () => openTableViewsExactly(page, "bdd-colman, bdd-colman"));
      await session.step(67, "When user closes the current view", () => closeCurrentView(page));
      await session.step(68, "Then the open table views should be exactly \"bdd-colman\"", () => openTableViewsExactly(page, "bdd-colman"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+L docks the layout and viewer suggestions into the shell, and each panel's × closes it", async () => {
      await session.step(72, "When user presses Control+L", () => pressKey(page, "Control+L"));
      await session.step(73, "Then \"Layout suggestions\" dock panel should be visible", () => shouldBe(page, el("\"Layout suggestions\" dock panel"), "visible"));
      await session.step(74, "And \"Viewer suggestions\" dock panel should be visible", () => shouldBe(page, el("\"Viewer suggestions\" dock panel"), "visible"));
      await session.step(75, "When user clicks on the close icon of \"Layout suggestions\" dock panel", () => clickOn(page, el("the close icon of \"Layout suggestions\" dock panel")));
      await session.step(76, "And user clicks on the close icon of \"Viewer suggestions\" dock panel", () => clickOn(page, el("the close icon of \"Viewer suggestions\" dock panel")));
      await session.step(77, "Then \"Layout suggestions\" dock panel should be absent", () => shouldBe(page, el("\"Layout suggestions\" dock panel"), "absent"));
      await session.step(78, "And \"Viewer suggestions\" dock panel should be absent", () => shouldBe(page, el("\"Viewer suggestions\" dock panel"), "absent"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Alt+C docks the Columns pane with a row per column, a click on a row makes that column current, and Alt+C closes it", async () => {
      await session.step(82, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(83, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(84, "Then the current column should be \"name\"", () => currentColumnIs(page, "name"));
      await session.step(85, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(86, "Then column manager should be visible", () => shouldBe(page, el("column manager"), "visible"));
      await session.step(87, "And the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(88, "And the \"text of cell 1 of __name\" reading of Grid viewer in column manager should be \"name\"", () => readingReads(page, "text of cell 1 of __name", el("Grid viewer in column manager"), "name"));
      await session.step(89, "And the \"text of cell 4 of __name\" reading of Grid viewer in column manager should be \"date\"", () => readingReads(page, "text of cell 4 of __name", el("Grid viewer in column manager"), "date"));
      await session.step(90, "And the \"text of cell 1 of __t\" reading of Grid viewer in column manager should be \"s\"", () => readingReads(page, "text of cell 1 of __t", el("Grid viewer in column manager"), "s"));
      await session.step(91, "And the \"text of cell 2 of __t\" reading of Grid viewer in column manager should be \"i\"", () => readingReads(page, "text of cell 2 of __t", el("Grid viewer in column manager"), "i"));
      await session.step(92, "And the \"text of cell 3 of __t\" reading of Grid viewer in column manager should be \"d\"", () => readingReads(page, "text of cell 3 of __t", el("Grid viewer in column manager"), "d"));
      await session.step(93, "When user clicks on the \"cell 3 of __name\" area of Grid viewer in column manager", () => clickArea(page, "cell 3 of __name", el("Grid viewer in column manager")));
      await session.step(94, "Then the current column should be \"price\"", () => currentColumnIs(page, "price"));
      await session.step(95, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(96, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The search box narrows the list to the matching columns", async () => {
      await session.step(100, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(101, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(102, "Then the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(103, "When user enters \"pr\" into \"Search\" input in column manager", () => enterInto(page, "pr", el("\"Search\" input in column manager")));
      await session.step(104, "Then the \"rows shown\" reading of Grid viewer in column manager should be 1", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 1));
      await session.step(105, "And the \"text of cell 3 of __name\" reading of Grid viewer in column manager should be \"price\"", () => readingReads(page, "text of cell 3 of __name", el("Grid viewer in column manager"), "price"));
      await session.step(106, "When user enters \"a\" into \"Search\" input in column manager", () => enterInto(page, "a", el("\"Search\" input in column manager")));
      await session.step(107, "Then the \"rows shown\" reading of Grid viewer in column manager should be 2", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 2));
      await session.step(108, "When user clears \"Search\" input in column manager", () => clearField(page, el("\"Search\" input in column manager")));
      await session.step(109, "Then the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(110, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(111, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+click toggles a column's selection in the grid, and Esc clears it", async () => {
      await session.step(115, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(116, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(117, "Then no columns should be selected", () => noColumnsSelected(page));
      await session.step(118, "When user clicks on the \"cell 2 of __name\" area of Grid viewer in column manager holding Control", () => clickAreaHolding(page, "cell 2 of __name", el("Grid viewer in column manager"), "Control"));
      await session.step(119, "And user clicks on the \"cell 4 of __name\" area of Grid viewer in column manager holding Control", () => clickAreaHolding(page, "cell 4 of __name", el("Grid viewer in column manager"), "Control"));
      await session.step(120, "Then columns \"num, date\" should be selected", () => columnsSelected(page, "num, date"));
      await session.step(121, "When user clicks on the \"cell 4 of __name\" area of Grid viewer in column manager holding Control", () => clickAreaHolding(page, "cell 4 of __name", el("Grid viewer in column manager"), "Control"));
      await session.step(122, "Then columns \"num\" should be selected", () => columnsSelected(page, "num"));
      await session.step(123, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(124, "Then no columns should be selected", () => noColumnsSelected(page));
      await session.step(125, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(126, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Dragging a column's row reorders the columns of the grid", async () => {
      await session.step(130, "Then the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(131, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(132, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(133, "And user drags the \"cell 4 of __name\" area of Grid viewer in column manager to the \"cell 1 of __name\" area", () => dragAreaToArea(page, "cell 4 of __name", el("Grid viewer in column manager"), "cell 1 of __name"));
      await session.step(134, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(135, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(136, "And the \"column order\" reading of grid should be \"date, name, num, price\"", () => readingReads(page, "column order", el("grid"), "date, name, num, price"));
      await session.step(137, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(138, "Then the \"text of cell 1 of __name\" reading of Grid viewer in column manager should be \"date\"", () => readingReads(page, "text of cell 1 of __name", el("Grid viewer in column manager"), "date"));
      await session.step(139, "When user drags the \"cell 1 of __name\" area of Grid viewer in column manager to the \"cell 4 of __name\" area", () => dragAreaToArea(page, "cell 1 of __name", el("Grid viewer in column manager"), "cell 4 of __name"));
      await session.step(140, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(141, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(142, "And the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(143, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hide in the popup menu takes the column out of the grid, and Show brings it back", async () => {
      await session.step(146, "Then the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(147, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(148, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(149, "And user picks \"'num' column > Hide\" from the context menu of the \"cell 2 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "'num' column > Hide", "cell 2 of __name", el("Grid viewer in column manager")));
      await session.step(150, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(151, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(152, "And the \"column order\" reading of grid should be \"name, price, date\"", () => readingReads(page, "column order", el("grid"), "name, price, date"));
      await session.step(153, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(154, "Then the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(155, "When user picks \"'num' column > Show\" from the context menu of the \"cell 2 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "'num' column > Show", "cell 2 of __name", el("Grid viewer in column manager")));
      await session.step(156, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(157, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(158, "And the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(159, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Statistics in the popup menu adds a statistics column to the list", async () => {
      await session.step(162, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(163, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(164, "Then the \"column order\" reading of Grid viewer in column manager should be \"__t, __name\"", () => readingReads(page, "column order", el("Grid viewer in column manager"), "__t, __name"));
      await session.step(165, "When user picks \"Statistics > count\" from the context menu of the \"cell 1 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "Statistics > count", "cell 1 of __name", el("Grid viewer in column manager")));
      await session.step(166, "Then the \"column order\" reading of Grid viewer in column manager should be \"__t, __name, count\"", () => readingReads(page, "column order", el("Grid viewer in column manager"), "__t, __name, count"));
      await session.step(167, "And the \"text of cell 2 of count\" reading of Grid viewer in column manager should be \"3.00\"", () => readingReads(page, "text of cell 2 of count", el("Grid viewer in column manager"), "3.00"));
      await session.step(168, "When user picks \"Statistics > count\" from the context menu of the \"cell 1 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "Statistics > count", "cell 1 of __name", el("Grid viewer in column manager")));
      await session.step(169, "Then the \"column order\" reading of Grid viewer in column manager should be \"__t, __name\"", () => readingReads(page, "column order", el("Grid viewer in column manager"), "__t, __name"));
      await session.step(170, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(171, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(172, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
