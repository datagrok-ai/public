/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/edit-add-rows.feature
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
import {clickOn, enterInto, pressKey, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {makeRowCurrent, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Edit > Add Rows inserts empty rows at a position, and Ctrl+Z takes them out again", () => {
  const session = feature(test, "features/table-view/edit-add-rows.feature", import.meta.url);
  test("Count rows go in before the row at Position, the existing rows keep their values, and Ctrl+Z removes them", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens a table \"items\" with:", () => openTableOf(page, "items", [["name","qty"],["one","1"],["two","2"],["three","3"]]), [["name","qty"],["one","1"],["two","2"],["three","3"]]);
    await session.step(21, "When user picks \"Edit > Add Rows...\" from the top menu", () => pickFromTopMenu(page, "Edit > Add Rows..."));
    await session.step(22, "Then Add Rows dialog should be visible", () => shouldBe(page, el("Add Rows dialog"), "visible"));
    await session.step(23, "And Count input in Add Rows dialog should have value \"1\"", () => shouldHaveValue(page, el("Count input in Add Rows dialog"), "1"));
    await session.step(24, "When user enters \"1\" into Position input in Add Rows dialog", () => enterInto(page, "1", el("Position input in Add Rows dialog")));
    await session.step(25, "And user enters \"2\" into Count input in Add Rows dialog", () => enterInto(page, "2", el("Count input in Add Rows dialog")));
    await session.step(26, "And user clicks on OK button in Add Rows dialog", () => clickOn(page, el("OK button in Add Rows dialog")));
    await session.step(27, "Then Add Rows dialog should be hidden", () => shouldBe(page, el("Add Rows dialog"), "hidden"));
    await session.step(28, "And the table should have 5 rows", () => rowCount(page, 5));
    await session.step(29, "And the value of \"name\" column in row 1 should be \"one\"", () => valueInRow(page, "name", 1, "one"));
    await session.step(30, "And the value of \"name\" column in row 2 should be \"\"", () => valueInRow(page, "name", 2, ""));
    await session.step(31, "And the value of \"qty\" column in row 2 should be \"\"", () => valueInRow(page, "qty", 2, ""));
    await session.step(32, "And the value of \"name\" column in row 3 should be \"\"", () => valueInRow(page, "name", 3, ""));
    await session.step(33, "And the value of \"name\" column in row 4 should be \"two\"", () => valueInRow(page, "name", 4, "two"));
    await session.step(34, "And the value of \"qty\" column in row 4 should be \"2\"", () => valueInRow(page, "qty", 4, "2"));
    await session.step(35, "And the value of \"name\" column in row 5 should be \"three\"", () => valueInRow(page, "name", 5, "three"));
    await session.step(36, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
    await session.step(37, "Then the table should have 3 rows", () => rowCount(page, 3));
    await session.step(38, "And the value of \"name\" column in row 2 should be \"two\"", () => valueInRow(page, "name", 2, "two"));
    await session.step(39, "And the value of \"name\" column in row 3 should be \"three\"", () => valueInRow(page, "name", 3, "three"));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
  });
  test("Position starts at the current row", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens a table \"items\" with:", () => openTableOf(page, "items", [["name","qty"],["one","1"],["two","2"],["three","3"]]), [["name","qty"],["one","1"],["two","2"],["three","3"]]);
    await session.step(43, "Given user makes row 3 current", () => makeRowCurrent(page, 3));
    await session.step(44, "When user picks \"Edit > Add Rows...\" from the top menu", () => pickFromTopMenu(page, "Edit > Add Rows..."));
    await session.step(45, "Then Position input in Add Rows dialog should have value \"2\"", () => shouldHaveValue(page, el("Position input in Add Rows dialog"), "2"));
    await session.step(46, "When user clicks on OK button in Add Rows dialog", () => clickOn(page, el("OK button in Add Rows dialog")));
    await session.step(47, "Then the table should have 4 rows", () => rowCount(page, 4));
    await session.step(48, "And the value of \"name\" column in row 3 should be \"\"", () => valueInRow(page, "name", 3, ""));
    await session.step(49, "And the value of \"name\" column in row 4 should be \"three\"", () => valueInRow(page, "name", 4, "three"));
    await session.step(50, "And the value of \"name\" column in row 2 should be \"two\"", () => valueInRow(page, "name", 2, "two"));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("Position equal to the row count appends the rows at the end", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens a table \"items\" with:", () => openTableOf(page, "items", [["name","qty"],["one","1"],["two","2"],["three","3"]]), [["name","qty"],["one","1"],["two","2"],["three","3"]]);
    await session.step(54, "When user picks \"Edit > Add Rows...\" from the top menu", () => pickFromTopMenu(page, "Edit > Add Rows..."));
    await session.step(55, "And user enters \"3\" into Position input in Add Rows dialog", () => enterInto(page, "3", el("Position input in Add Rows dialog")));
    await session.step(56, "And user clicks on OK button in Add Rows dialog", () => clickOn(page, el("OK button in Add Rows dialog")));
    await session.step(57, "Then the table should have 4 rows", () => rowCount(page, 4));
    await session.step(58, "And the value of \"name\" column in row 3 should be \"three\"", () => valueInRow(page, "name", 3, "three"));
    await session.step(59, "And the value of \"name\" column in row 4 should be \"\"", () => valueInRow(page, "name", 4, ""));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
});
