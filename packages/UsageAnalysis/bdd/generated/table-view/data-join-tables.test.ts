/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/data-join-tables.feature
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
import {check, clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly, displayedInRowOfTable, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, selectedInTable, setCell, tableColumnIncomplete, tableColumns, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors, pickFromOpenMenu, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Data > Join Tables joins two open tables on their key columns into a new table view", () => {
  const session = feature(test, "features/table-view/data-join-tables.feature", import.meta.url);
  test("The dialog opens on both tables keyed by their shared column, and OK opens the inner join as a new table view", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens a table \"left\" with:", () => openTableOf(page, "left", [["id","name"],["1","a"],["2","b"],["3","c"]]), [["id","name"],["1","a"],["2","b"],["3","c"]]);
    await session.step(28, "And user opens a table \"right\" with:", () => openTableOf(page, "right", [["id","score"],["2","20"],["3","30"],["4","40"]]), [["id","score"],["2","20"],["3","30"],["4","40"]]);
    await session.step(35, "When user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
    await session.step(36, "Then Join Tables dialog should be visible", () => shouldBe(page, el("Join Tables dialog"), "visible"));
    await session.step(37, "And Key Columns input in Join Tables dialog should contain text \"id\"", () => shouldContainText(page, el("Key Columns input in Join Tables dialog"), "id"));
    await session.step(38, "And \"Join Type\" input in Join Tables dialog should have value \"inner\"", () => shouldHaveValue(page, el("\"Join Type\" input in Join Tables dialog"), "inner"));
    await session.step(39, "And grid viewer in \"Preview result columns\" dialog should show 2 rows", () => showsRows(page, el("grid viewer in \"Preview result columns\" dialog"), 2));
    await session.step(40, "And OK button in Join Tables dialog should be enabled", () => shouldBe(page, el("OK button in Join Tables dialog"), "enabled"));
    await session.step(41, "When user clicks on OK button in Join Tables dialog", () => clickOn(page, el("OK button in Join Tables dialog")));
    await session.step(42, "Then Join Tables dialog should be hidden", () => shouldBe(page, el("Join Tables dialog"), "hidden"));
    await session.step(43, "And the open table views should be exactly \"left, right, result\"", () => openTableViewsExactly(page, "left, right, result"));
    await session.step(44, "Given user switches to the \"result\" table view", () => switchTableView(page, "result"));
    await session.step(45, "Then the table should have 2 rows", () => rowCount(page, 2));
    await session.step(46, "And the table should have the columns \"id, name, right.id, score\"", () => columnsExactly(page, "id, name, right.id, score"));
    await session.step(47, "And the value of \"name\" column in row 1 should be \"b\"", () => valueInRow(page, "name", 1, "b"));
    await session.step(48, "And the value of \"score\" column in row 1 should be \"20\"", () => valueInRow(page, "score", 1, "20"));
    await session.step(49, "And the value of \"name\" column in row 2 should be \"c\"", () => valueInRow(page, "name", 2, "c"));
    await session.step(50, "And the value of \"score\" column in row 2 should be \"30\"", () => valueInRow(page, "score", 2, "30"));
    await session.step(51, "And table \"left\" should have 3 rows", () => tableRows(page, "left", 3));
    await session.step(52, "And table \"right\" should have 3 rows", () => tableRows(page, "right", 3));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
  test("The dialog counts the key rows without a match, and its menu selects them in that table", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens a table \"left\" with:", () => openTableOf(page, "left", [["id","name"],["1","a"],["2","b"],["3","c"]]), [["id","name"],["1","a"],["2","b"],["3","c"]]);
    await session.step(28, "And user opens a table \"right\" with:", () => openTableOf(page, "right", [["id","score"],["2","20"],["3","30"],["4","40"]]), [["id","score"],["2","20"],["3","30"],["4","40"]]);
    await session.step(56, "When user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
    await session.step(57, "Then Join Tables dialog should contain text \"1 / 3\"", () => shouldContainText(page, el("Join Tables dialog"), "1 / 3"));
    await session.step(58, "When user clicks on first \"context arrow down\" icon in Join Tables dialog", () => clickOn(page, el("first \"context arrow down\" icon in Join Tables dialog")));
    await session.step(59, "And user picks \"Select All Not Matching\" from the open menu", () => pickFromOpenMenu(page, "Select All Not Matching"));
    await session.step(60, "Then 1 row of table \"left\" should be selected", () => selectedInTable(page, 1, "left"));
    await session.step(61, "And 0 rows of table \"right\" should be selected", () => selectedInTable(page, 0, "right"));
    await session.step(62, "When user clicks on CANCEL button in Join Tables dialog", () => clickOn(page, el("CANCEL button in Join Tables dialog")));
    await session.step(63, "Then the open table views should be exactly \"left, right\"", () => openTableViewsExactly(page, "left, right"));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
  test("Join Type is passed on to the join: outer keeps every row of both tables", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens a table \"left\" with:", () => openTableOf(page, "left", [["id","name"],["1","a"],["2","b"],["3","c"]]), [["id","name"],["1","a"],["2","b"],["3","c"]]);
    await session.step(28, "And user opens a table \"right\" with:", () => openTableOf(page, "right", [["id","score"],["2","20"],["3","30"],["4","40"]]), [["id","score"],["2","20"],["3","30"],["4","40"]]);
    await session.step(67, "When user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
    await session.step(68, "And user selects \"outer\" in \"Join Type\" input in Join Tables dialog", () => selectIn(page, "outer", el("\"Join Type\" input in Join Tables dialog")));
    await session.step(69, "Then grid viewer in \"Preview result columns\" dialog should show 4 rows", () => showsRows(page, el("grid viewer in \"Preview result columns\" dialog"), 4));
    await session.step(70, "When user clicks on OK button in Join Tables dialog", () => clickOn(page, el("OK button in Join Tables dialog")));
    await session.step(71, "Then the open table views should be exactly \"left, right, result\"", () => openTableViewsExactly(page, "left, right, result"));
    await session.step(72, "Given user switches to the \"result\" table view", () => switchTableView(page, "result"));
    await session.step(73, "Then the table should have 4 rows", () => rowCount(page, 4));
    await session.step(74, "And the table should have the columns \"id, name, right.id, score\"", () => columnsExactly(page, "id, name, right.id, score"));
    await session.step(75, "And no errors should have been logged", () => noErrors(page));
  });
  test("Values decides which columns of the table go into the result", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens a table \"left\" with:", () => openTableOf(page, "left", [["id","name"],["1","a"],["2","b"],["3","c"]]), [["id","name"],["1","a"],["2","b"],["3","c"]]);
    await session.step(28, "And user opens a table \"right\" with:", () => openTableOf(page, "right", [["id","score"],["2","20"],["3","30"],["4","40"]]), [["id","score"],["2","20"],["3","30"],["4","40"]]);
    await session.step(78, "When user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
    await session.step(79, "And user clicks on Values input in Join Tables dialog", () => clickOn(page, el("Values input in Join Tables dialog")));
    await session.step(80, "And user toggles the \"name\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "name", el("\"Select columns...\" dialog")));
    await session.step(81, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(82, "And user clicks on OK button in Join Tables dialog", () => clickOn(page, el("OK button in Join Tables dialog")));
    await session.step(83, "Then the open table views should be exactly \"left, right, result\"", () => openTableViewsExactly(page, "left, right, result"));
    await session.step(84, "Given user switches to the \"result\" table view", () => switchTableView(page, "result"));
    await session.step(85, "Then the table should have 2 rows", () => rowCount(page, 2));
    await session.step(86, "And the table should have the columns \"id, right.id, score\"", () => columnsExactly(page, "id, right.id, score"));
    await session.step(87, "And no errors should have been logged", () => noErrors(page));
  });
  test("No key row matching says so and disables OK", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens a table \"left\" with:", () => openTableOf(page, "left", [["id","name"],["1","a"],["2","b"],["3","c"]]), [["id","name"],["1","a"],["2","b"],["3","c"]]);
    await session.step(28, "And user opens a table \"right\" with:", () => openTableOf(page, "right", [["id","score"],["2","20"],["3","30"],["4","40"]]), [["id","score"],["2","20"],["3","30"],["4","40"]]);
    await session.step(90, "Given user sets \"id\" column in row 1 to \"7\"", () => setCell(page, "id", 1, "7"));
    await session.step(91, "And user sets \"id\" column in row 2 to \"8\"", () => setCell(page, "id", 2, "8"));
    await session.step(92, "And user sets \"id\" column in row 3 to \"9\"", () => setCell(page, "id", 3, "9"));
    await session.step(93, "When user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
    await session.step(94, "Then Join Tables dialog should contain text \"None of key rows are matched\"", () => shouldContainText(page, el("Join Tables dialog"), "None of key rows are matched"));
    await session.step(95, "And OK button in Join Tables dialog should be disabled", () => shouldBe(page, el("OK button in Join Tables dialog"), "disabled"));
    await session.step(96, "When user clicks on CANCEL button in Join Tables dialog", () => clickOn(page, el("CANCEL button in Join Tables dialog")));
    await session.step(97, "Then the open table views should be exactly \"left, right\"", () => openTableViewsExactly(page, "left, right"));
    await session.step(98, "And no errors should have been logged", () => noErrors(page));
  });
  test("In-place is offered only for left and right, and a left join in place adds the columns to the left table", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens a table \"left\" with:", () => openTableOf(page, "left", [["id","name"],["1","a"],["2","b"],["3","c"]]), [["id","name"],["1","a"],["2","b"],["3","c"]]);
    await session.step(28, "And user opens a table \"right\" with:", () => openTableOf(page, "right", [["id","score"],["2","20"],["3","30"],["4","40"]]), [["id","score"],["2","20"],["3","30"],["4","40"]]);
    await session.step(101, "When user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
    await session.step(102, "Then \"In-place\" input in Join Tables dialog should be disabled", () => shouldBe(page, el("\"In-place\" input in Join Tables dialog"), "disabled"));
    await session.step(103, "When user selects \"left\" in \"Join Type\" input in Join Tables dialog", () => selectIn(page, "left", el("\"Join Type\" input in Join Tables dialog")));
    await session.step(104, "Then \"In-place\" input in Join Tables dialog should be enabled", () => shouldBe(page, el("\"In-place\" input in Join Tables dialog"), "enabled"));
    await session.step(105, "When user checks \"In-place\" input in Join Tables dialog", () => check(page, el("\"In-place\" input in Join Tables dialog")));
    await session.step(106, "And user clicks on OK button in Join Tables dialog", () => clickOn(page, el("OK button in Join Tables dialog")));
    await session.step(107, "Then Join Tables dialog should be hidden", () => shouldBe(page, el("Join Tables dialog"), "hidden"));
    await session.step(108, "And the open table views should be exactly \"left+right, right\"", () => openTableViewsExactly(page, "left+right, right"));
    await session.step(109, "And table \"left+right\" should have 3 rows", () => tableRows(page, "left+right", 3));
    await session.step(110, "And table \"left+right\" should have columns \"id, name, right.id, score\"", () => tableColumns(page, "left+right", "id, name, right.id, score"));
    await session.step(111, "And the \"score\" cell of row 2 of table \"left+right\" should be displayed as \"20\"", () => displayedInRowOfTable(page, "score", 2, "left+right", "20"));
    await session.step(112, "And table \"left+right\" should have missing values in \"score\" column", () => tableColumnIncomplete(page, "left+right", "score"));
    await session.step(113, "And no errors should have been logged", () => noErrors(page));
  });
});
