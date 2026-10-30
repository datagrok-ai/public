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
import {check, clickOn, enterInto, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {allRowsSelected, columnsSelected, noneSelected, onlyOfAnySelected, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Select menu selects rows and columns of the current table, and Selection to Column keeps the selection", () => {
  const session = feature(test, "features/table-view/select-menu.feature", import.meta.url);
  test("Duplicates selects every row whose key repeats another row's", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]), [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]);
    await session.step(29, "When user picks \"Select > Duplicates...\" from the top menu", () => pickFromTopMenu(page, "Select > Duplicates..."));
    await session.step(30, "Then Select Duplicates dialog should be visible", () => shouldBe(page, el("Select Duplicates dialog"), "visible"));
    await session.step(31, "When user toggles the \"name\" column in the column list of Select Duplicates dialog", () => toggleInColumnList(page, "name", el("Select Duplicates dialog")));
    await session.step(32, "And user clicks on OK button in Select Duplicates dialog", () => clickOn(page, el("OK button in Select Duplicates dialog")));
    await session.step(33, "Then Select Duplicates dialog should be hidden", () => shouldBe(page, el("Select Duplicates dialog"), "hidden"));
    await session.step(34, "And only rows where \"name\" is one of \"a, b\" should be selected", () => onlyOfAnySelected(page, "name", "a, b"));
    await session.step(35, "And 4 rows should be selected", () => selectedRowCount(page, 4));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
  });
  test("Superfluous rows only leaves one row of each duplicate group unselected", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]), [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]);
    await session.step(39, "When user picks \"Select > Duplicates...\" from the top menu", () => pickFromTopMenu(page, "Select > Duplicates..."));
    await session.step(40, "And user toggles the \"name\" column in the column list of Select Duplicates dialog", () => toggleInColumnList(page, "name", el("Select Duplicates dialog")));
    await session.step(41, "And user checks \"Superfluous rows only\" input in Select Duplicates dialog", () => check(page, el("\"Superfluous rows only\" input in Select Duplicates dialog")));
    await session.step(42, "And user clicks on OK button in Select Duplicates dialog", () => clickOn(page, el("OK button in Select Duplicates dialog")));
    await session.step(43, "Then only rows where \"id\" is one of \"2, 4\" should be selected", () => onlyOfAnySelected(page, "id", "2, 4"));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
  test("Missing Values selects the rows empty in the checked column", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]), [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]);
    await session.step(47, "When user picks \"Select > Missing Values...\" from the top menu", () => pickFromTopMenu(page, "Select > Missing Values..."));
    await session.step(48, "Then Select Missing Values dialog should be visible", () => shouldBe(page, el("Select Missing Values dialog"), "visible"));
    await session.step(49, "When user toggles the \"v\" column in the column list of Select Missing Values dialog", () => toggleInColumnList(page, "v", el("Select Missing Values dialog")));
    await session.step(50, "And user clicks on OK button in Select Missing Values dialog", () => clickOn(page, el("OK button in Select Missing Values dialog")));
    await session.step(51, "Then Select Missing Values dialog should be hidden", () => shouldBe(page, el("Select Missing Values dialog"), "hidden"));
    await session.step(52, "And only rows where \"id\" is one of \"2, 4\" should be selected", () => onlyOfAnySelected(page, "id", "2, 4"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
  test("Random Rows selects as many rows as typed, and a proportion link sets the count", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]), [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]);
    await session.step(56, "When user picks \"Select > Random Rows...\" from the top menu", () => pickFromTopMenu(page, "Select > Random Rows..."));
    await session.step(57, "Then Select Random Rows dialog should be visible", () => shouldBe(page, el("Select Random Rows dialog"), "visible"));
    await session.step(58, "When user enters \"3\" into \"# rows\" input in Select Random Rows dialog", () => enterInto(page, "3", el("\"# rows\" input in Select Random Rows dialog")));
    await session.step(59, "And user clicks on OK button in Select Random Rows dialog", () => clickOn(page, el("OK button in Select Random Rows dialog")));
    await session.step(60, "Then Select Random Rows dialog should be hidden", () => shouldBe(page, el("Select Random Rows dialog"), "hidden"));
    await session.step(61, "And 3 rows should be selected", () => selectedRowCount(page, 3));
    await session.step(62, "When user picks \"Select > Random Rows...\" from the top menu", () => pickFromTopMenu(page, "Select > Random Rows..."));
    await session.step(63, "And user clicks on \"100%\" label in Select Random Rows dialog", () => clickOn(page, el("\"100%\" label in Select Random Rows dialog")));
    await session.step(64, "Then \"# rows\" input in Select Random Rows dialog should have value \"5\"", () => shouldHaveValue(page, el("\"# rows\" input in Select Random Rows dialog"), "5"));
    await session.step(65, "When user clicks on OK button in Select Random Rows dialog", () => clickOn(page, el("OK button in Select Random Rows dialog")));
    await session.step(66, "Then all rows should be selected", () => allRowsSelected(page));
    await session.step(67, "And no errors should have been logged", () => noErrors(page));
  });
  test("All selects every row and every column, All Columns only the columns", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]), [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]);
    await session.step(70, "When user picks \"Select > All Columns\" from the top menu", () => pickFromTopMenu(page, "Select > All Columns"));
    await session.step(71, "Then columns \"id, name, v\" should be selected", () => columnsSelected(page, "id, name, v"));
    await session.step(72, "And no rows should be selected", () => noneSelected(page));
    await session.step(73, "When user picks \"Select > All\" from the top menu", () => pickFromTopMenu(page, "Select > All"));
    await session.step(74, "Then all rows should be selected", () => allRowsSelected(page));
    await session.step(75, "And columns \"id, name, v\" should be selected", () => columnsSelected(page, "id, name, v"));
    await session.step(76, "And no errors should have been logged", () => noErrors(page));
  });
  test("Selection to Column adds a boolean column true in the selected rows", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]), [["id","name","v"],["1","a","10"],["2","a",""],["3","b","30"],["4","b",""],["5","c","50"]]);
    await session.step(79, "Given user selects rows where \"name\" is \"b\"", () => selectWhereIs(page, "name", "b"));
    await session.step(80, "When user picks \"Select > Selection to Column...\" from the top menu", () => pickFromTopMenu(page, "Select > Selection to Column..."));
    await session.step(81, "Then Selection to Column dialog should be visible", () => shouldBe(page, el("Selection to Column dialog"), "visible"));
    await session.step(82, "And Name input in Selection to Column dialog should have value \"Selection\"", () => shouldHaveValue(page, el("Name input in Selection to Column dialog"), "Selection"));
    await session.step(83, "When user clicks on OK button in Selection to Column dialog", () => clickOn(page, el("OK button in Selection to Column dialog")));
    await session.step(84, "Then Selection to Column dialog should be hidden", () => shouldBe(page, el("Selection to Column dialog"), "hidden"));
    await session.step(85, "And a new column \"Selection\" should have been added", () => newColumnNamed(page, "Selection"));
    await session.step(86, "And \"Selection\" column should have type \"bool\"", () => columnType(page, "Selection", "bool"));
    await session.step(87, "And the value of \"Selection\" column in row 3 should be \"true\"", () => valueInRow(page, "Selection", 3, "true"));
    await session.step(88, "And the value of \"Selection\" column in row 4 should be \"true\"", () => valueInRow(page, "Selection", 4, "true"));
    await session.step(89, "And the value of \"Selection\" column in row 1 should be \"false\"", () => valueInRow(page, "Selection", 1, "false"));
    await session.step(90, "And the value of \"Selection\" column in row 5 should be \"false\"", () => valueInRow(page, "Selection", 5, "false"));
    await session.step(91, "And no errors should have been logged", () => noErrors(page));
  });
});
