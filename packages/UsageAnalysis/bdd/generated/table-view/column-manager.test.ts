/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/column-manager.feature
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
import {clearField, enterInto, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsSelected, currentColumnIs, noColumnsSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, dragAreaToArea, noErrors, pickFromAreaContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Column Manager lists the table's columns and acts on the table view's grid", () => {
  const session = feature(test, "features/table-view/column-manager.feature", import.meta.url);
  test("Column Manager lists the table's columns and acts on the table view's grid", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens a table \"bdd-colman\" with:", () => openTableOf(page, "bdd-colman", [["name","num","price","date"],["a","1","1.5","2020-01-01"],["b","2","2.5","2020-01-02"],["c","3","3.5","2020-01-03"]]), [["name","num","price","date"],["a","1","1.5","2020-01-01"],["b","2","2.5","2020-01-02"],["c","3","3.5","2020-01-03"]]);
    await run.scenario("Alt+C docks the Columns pane with a row per column, and Alt+C closes it", async () => {
      await session.step(28, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(29, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(30, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(31, "Then column manager should be visible", () => shouldBe(page, el("column manager"), "visible"));
      await session.step(32, "And the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(33, "And the \"text of cell 1 of __name\" reading of Grid viewer in column manager should be \"name\"", () => readingReads(page, "text of cell 1 of __name", el("Grid viewer in column manager"), "name"));
      await session.step(34, "And the \"text of cell 4 of __name\" reading of Grid viewer in column manager should be \"date\"", () => readingReads(page, "text of cell 4 of __name", el("Grid viewer in column manager"), "date"));
      await session.step(35, "And the \"text of cell 1 of __t\" reading of Grid viewer in column manager should be \"s\"", () => readingReads(page, "text of cell 1 of __t", el("Grid viewer in column manager"), "s"));
      await session.step(36, "And the \"text of cell 2 of __t\" reading of Grid viewer in column manager should be \"i\"", () => readingReads(page, "text of cell 2 of __t", el("Grid viewer in column manager"), "i"));
      await session.step(37, "And the \"text of cell 3 of __t\" reading of Grid viewer in column manager should be \"d\"", () => readingReads(page, "text of cell 3 of __t", el("Grid viewer in column manager"), "d"));
      await session.step(38, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(39, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The search box narrows the list to the matching columns", async () => {
      await session.step(43, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(44, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(45, "Then the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(46, "When user enters \"pr\" into \"Search\" input in column manager", () => enterInto(page, "pr", el("\"Search\" input in column manager")));
      await session.step(47, "Then the \"rows shown\" reading of Grid viewer in column manager should be 1", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 1));
      await session.step(48, "And the \"text of cell 3 of __name\" reading of Grid viewer in column manager should be \"price\"", () => readingReads(page, "text of cell 3 of __name", el("Grid viewer in column manager"), "price"));
      await session.step(49, "When user enters \"a\" into \"Search\" input in column manager", () => enterInto(page, "a", el("\"Search\" input in column manager")));
      await session.step(50, "Then the \"rows shown\" reading of Grid viewer in column manager should be 2", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 2));
      await session.step(51, "When user clears \"Search\" input in column manager", () => clearField(page, el("\"Search\" input in column manager")));
      await session.step(52, "Then the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(53, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(54, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a column's row makes it the current column", async () => {
      await session.step(58, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(59, "Then the current column should be \"name\"", () => currentColumnIs(page, "name"));
      await session.step(60, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(61, "And user clicks on the \"cell 3 of __name\" area of Grid viewer in column manager", () => clickArea(page, "cell 3 of __name", el("Grid viewer in column manager")));
      await session.step(62, "Then the current column should be \"price\"", () => currentColumnIs(page, "price"));
      await session.step(63, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(64, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+click toggles a column's selection in the grid, and Esc clears it", async () => {
      await session.step(68, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(69, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(70, "Then no columns should be selected", () => noColumnsSelected(page));
      await session.step(71, "When user clicks on the \"cell 2 of __name\" area of Grid viewer in column manager holding Control", () => clickAreaHolding(page, "cell 2 of __name", el("Grid viewer in column manager"), "Control"));
      await session.step(72, "And user clicks on the \"cell 4 of __name\" area of Grid viewer in column manager holding Control", () => clickAreaHolding(page, "cell 4 of __name", el("Grid viewer in column manager"), "Control"));
      await session.step(73, "Then columns \"num, date\" should be selected", () => columnsSelected(page, "num, date"));
      await session.step(74, "When user clicks on the \"cell 4 of __name\" area of Grid viewer in column manager holding Control", () => clickAreaHolding(page, "cell 4 of __name", el("Grid viewer in column manager"), "Control"));
      await session.step(75, "Then columns \"num\" should be selected", () => columnsSelected(page, "num"));
      await session.step(76, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(77, "Then no columns should be selected", () => noColumnsSelected(page));
      await session.step(78, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(79, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Dragging a column's row reorders the columns of the grid", async () => {
      await session.step(83, "Then the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(84, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(85, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(86, "And user drags the \"cell 4 of __name\" area of Grid viewer in column manager to the \"cell 1 of __name\" area", () => dragAreaToArea(page, "cell 4 of __name", el("Grid viewer in column manager"), "cell 1 of __name"));
      await session.step(87, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(88, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(89, "And the \"column order\" reading of grid should be \"date, name, num, price\"", () => readingReads(page, "column order", el("grid"), "date, name, num, price"));
      await session.step(90, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(91, "Then the \"text of cell 1 of __name\" reading of Grid viewer in column manager should be \"date\"", () => readingReads(page, "text of cell 1 of __name", el("Grid viewer in column manager"), "date"));
      await session.step(92, "When user drags the \"cell 1 of __name\" area of Grid viewer in column manager to the \"cell 4 of __name\" area", () => dragAreaToArea(page, "cell 1 of __name", el("Grid viewer in column manager"), "cell 4 of __name"));
      await session.step(93, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(94, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(95, "And the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hide in the popup menu takes the column out of the grid, and Show brings it back", async () => {
      await session.step(99, "Then the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(100, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(101, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(102, "And user picks \"'num' column > Hide\" from the context menu of the \"cell 2 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "'num' column > Hide", "cell 2 of __name", el("Grid viewer in column manager")));
      await session.step(103, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(104, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(105, "And the \"column order\" reading of grid should be \"name, price, date\"", () => readingReads(page, "column order", el("grid"), "name, price, date"));
      await session.step(106, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(107, "Then the \"rows shown\" reading of Grid viewer in column manager should be 4", () => readingIs(page, "rows shown", el("Grid viewer in column manager"), 4));
      await session.step(108, "When user picks \"'num' column > Show\" from the context menu of the \"cell 2 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "'num' column > Show", "cell 2 of __name", el("Grid viewer in column manager")));
      await session.step(109, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(110, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(111, "And the \"column order\" reading of grid should be \"name, num, price, date\"", () => readingReads(page, "column order", el("grid"), "name, num, price, date"));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Statistics in the popup menu adds a statistics column to the list", async () => {
      await session.step(115, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(116, "And user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(117, "Then the \"column order\" reading of Grid viewer in column manager should be \"__t, __name\"", () => readingReads(page, "column order", el("Grid viewer in column manager"), "__t, __name"));
      await session.step(118, "When user picks \"Statistics > count\" from the context menu of the \"cell 1 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "Statistics > count", "cell 1 of __name", el("Grid viewer in column manager")));
      await session.step(119, "Then the \"column order\" reading of Grid viewer in column manager should be \"__t, __name, count\"", () => readingReads(page, "column order", el("Grid viewer in column manager"), "__t, __name, count"));
      await session.step(120, "And the \"text of cell 2 of count\" reading of Grid viewer in column manager should be \"3.00\"", () => readingReads(page, "text of cell 2 of count", el("Grid viewer in column manager"), "3.00"));
      await session.step(121, "When user picks \"Statistics > count\" from the context menu of the \"cell 1 of __name\" area of Grid viewer in column manager", () => pickFromAreaContextMenu(page, "Statistics > count", "cell 1 of __name", el("Grid viewer in column manager")));
      await session.step(122, "Then the \"column order\" reading of Grid viewer in column manager should be \"__t, __name\"", () => readingReads(page, "column order", el("Grid viewer in column manager"), "__t, __name"));
      await session.step(123, "When user presses Alt+C", () => pressKey(page, "Alt+C"));
      await session.step(124, "Then column manager should be absent", () => shouldBe(page, el("column manager"), "absent"));
      await session.step(125, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
