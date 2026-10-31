/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/select-menu.feature
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
import {columnType, makeRowCurrent, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, onlyOfAnySelected, rowCount, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Select menu selects rows of the current table, Edit > Add Rows inserts empty rows, and Selection to Column keeps the selection", () => {
  const session = feature(test, "features/table-view/select-menu.feature", import.meta.url);
  test("The Select menu selects rows of the current table, Edit > Add Rows inserts empty rows, and Selection to Column keeps the selection", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]), [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]);
    await run.scenario("Duplicates selects every row whose key repeats another row's", async () => {
      await session.step(31, "When user picks \"Select > Duplicates...\" from the top menu", () => pickFromTopMenu(page, "Select > Duplicates..."));
      await session.step(32, "Then Select Duplicates dialog should be visible", () => shouldBe(page, el("Select Duplicates dialog"), "visible"));
      await session.step(33, "When user toggles the \"name\" column in the column list of Select Duplicates dialog", () => toggleInColumnList(page, "name", el("Select Duplicates dialog")));
      await session.step(34, "And user clicks on OK button in Select Duplicates dialog", () => clickOn(page, el("OK button in Select Duplicates dialog")));
      await session.step(35, "Then Select Duplicates dialog should be hidden", () => shouldBe(page, el("Select Duplicates dialog"), "hidden"));
      await session.step(36, "And only rows where \"name\" is one of \"a, b\" should be selected", () => onlyOfAnySelected(page, "name", "a, b"));
      await session.step(37, "And 4 rows should be selected", () => selectedRowCount(page, 4));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Missing Values selects the rows empty in the checked column", async () => {
      await session.step(41, "Given user clears the row selection", () => clearSelection(page));
      await session.step(42, "When user picks \"Select > Missing Values...\" from the top menu", () => pickFromTopMenu(page, "Select > Missing Values..."));
      await session.step(43, "Then Select Missing Values dialog should be visible", () => shouldBe(page, el("Select Missing Values dialog"), "visible"));
      await session.step(44, "When user toggles the \"v\" column in the column list of Select Missing Values dialog", () => toggleInColumnList(page, "v", el("Select Missing Values dialog")));
      await session.step(45, "And user clicks on OK button in Select Missing Values dialog", () => clickOn(page, el("OK button in Select Missing Values dialog")));
      await session.step(46, "Then Select Missing Values dialog should be hidden", () => shouldBe(page, el("Select Missing Values dialog"), "hidden"));
      await session.step(47, "And only rows where \"id\" is one of \"2, 4\" should be selected", () => onlyOfAnySelected(page, "id", "2, 4"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Random Rows selects as many rows as typed", async () => {
      await session.step(51, "When user picks \"Select > Random Rows...\" from the top menu", () => pickFromTopMenu(page, "Select > Random Rows..."));
      await session.step(52, "Then Select Random Rows dialog should be visible", () => shouldBe(page, el("Select Random Rows dialog"), "visible"));
      await session.step(53, "When user enters \"3\" into \"# rows\" input in Select Random Rows dialog", () => enterInto(page, "3", el("\"# rows\" input in Select Random Rows dialog")));
      await session.step(54, "And user clicks on OK button in Select Random Rows dialog", () => clickOn(page, el("OK button in Select Random Rows dialog")));
      await session.step(55, "Then Select Random Rows dialog should be hidden", () => shouldBe(page, el("Select Random Rows dialog"), "hidden"));
      await session.step(56, "And 3 rows should be selected", () => selectedRowCount(page, 3));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Add Rows inserts Count empty rows before the row at Position, which starts at the current row, and Ctrl+Z takes them out", async () => {
      await session.step(60, "Given user makes row 3 current", () => makeRowCurrent(page, 3));
      await session.step(61, "When user picks \"Edit > Add Rows...\" from the top menu", () => pickFromTopMenu(page, "Edit > Add Rows..."));
      await session.step(62, "Then Add Rows dialog should be visible", () => shouldBe(page, el("Add Rows dialog"), "visible"));
      await session.step(63, "And Count input in Add Rows dialog should have value \"1\"", () => shouldHaveValue(page, el("Count input in Add Rows dialog"), "1"));
      await session.step(64, "And Position input in Add Rows dialog should have value \"2\"", () => shouldHaveValue(page, el("Position input in Add Rows dialog"), "2"));
      await session.step(65, "When user enters \"1\" into Position input in Add Rows dialog", () => enterInto(page, "1", el("Position input in Add Rows dialog")));
      await session.step(66, "And user enters \"2\" into Count input in Add Rows dialog", () => enterInto(page, "2", el("Count input in Add Rows dialog")));
      await session.step(67, "And user clicks on OK button in Add Rows dialog", () => clickOn(page, el("OK button in Add Rows dialog")));
      await session.step(68, "Then Add Rows dialog should be hidden", () => shouldBe(page, el("Add Rows dialog"), "hidden"));
      await session.step(69, "And the table should have 7 rows", () => rowCount(page, 7));
      await session.step(70, "And the value of \"name\" column in row 1 should be \"a\"", () => valueInRow(page, "name", 1, "a"));
      await session.step(71, "And the value of \"name\" column in row 2 should be \"\"", () => valueInRow(page, "name", 2, ""));
      await session.step(72, "And the value of \"id\" column in row 3 should be \"\"", () => valueInRow(page, "id", 3, ""));
      await session.step(73, "And the value of \"id\" column in row 4 should be \"2\"", () => valueInRow(page, "id", 4, "2"));
      await session.step(74, "And the value of \"id\" column in row 7 should be \"5\"", () => valueInRow(page, "id", 7, "5"));
      await session.step(75, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(76, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(77, "Then the table should have 5 rows", () => rowCount(page, 5));
      await session.step(78, "And the value of \"id\" column in row 2 should be \"2\"", () => valueInRow(page, "id", 2, "2"));
      await session.step(79, "And the value of \"id\" column in row 3 should be \"3\"", () => valueInRow(page, "id", 3, "3"));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Selection to Column adds a boolean column true in the selected rows", async () => {
      await session.step(83, "Given user selects rows where \"name\" is \"b\"", () => selectWhereIs(page, "name", "b"));
      await session.step(84, "When user picks \"Select > Selection to Column...\" from the top menu", () => pickFromTopMenu(page, "Select > Selection to Column..."));
      await session.step(85, "Then Selection to Column dialog should be visible", () => shouldBe(page, el("Selection to Column dialog"), "visible"));
      await session.step(86, "And Name input in Selection to Column dialog should have value \"Selection\"", () => shouldHaveValue(page, el("Name input in Selection to Column dialog"), "Selection"));
      await session.step(87, "When user clicks on OK button in Selection to Column dialog", () => clickOn(page, el("OK button in Selection to Column dialog")));
      await session.step(88, "Then Selection to Column dialog should be hidden", () => shouldBe(page, el("Selection to Column dialog"), "hidden"));
      await session.step(89, "And a new column \"Selection\" should have been added", () => newColumnNamed(page, "Selection"));
      await session.step(90, "And \"Selection\" column should have type \"bool\"", () => columnType(page, "Selection", "bool"));
      await session.step(91, "And the value of \"Selection\" column in row 3 should be \"true\"", () => valueInRow(page, "Selection", 3, "true"));
      await session.step(92, "And the value of \"Selection\" column in row 4 should be \"true\"", () => valueInRow(page, "Selection", 4, "true"));
      await session.step(93, "And the value of \"Selection\" column in row 1 should be \"false\"", () => valueInRow(page, "Selection", 1, "false"));
      await session.step(94, "And the value of \"Selection\" column in row 5 should be \"false\"", () => valueInRow(page, "Selection", 5, "false"));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
