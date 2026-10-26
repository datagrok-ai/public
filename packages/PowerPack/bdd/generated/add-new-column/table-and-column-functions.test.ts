/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/table-and-column-functions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column, GROK-20931]
--- */
import {test} from '@playwright/test';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {acceptCompletion} from '../../bindings/add-new-column.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, holdsCode, pasteInto, pressKeyIn, shouldBe, shouldContainText, shouldNotBe, typeAtCaret, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, displayedInRow, everyValueMatches, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Add New Column: table and column arguments are offered as they are typed, and the table functions add their columns", () => {
  const session = feature(test, "features/add-new-column/table-and-column-functions.feature", import.meta.url);
  test("Add New Column: table and column arguments are offered as they are typed, and the table functions add their columns", {tag: ["@journey", "@realizes:powerpack.dialogs.add-new-column", "@realizes:GROK-20931"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(22, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(23, "And user switches to the \"demog\" table view", () => switchTableView(page, "demog"));
    await session.step(24, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
    await session.step(25, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
    await run.scenario("A column parameter offers the numeric columns of the table and takes ${col}", async () => {
      await session.step(28, "When user types \"CumSum(\" into formula editor", () => typeInto(page, "CumSum(", el("formula editor")));
      await session.step(29, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(30, "And \"AGE\" completion should be visible", () => shouldBe(page, el("\"AGE\" completion"), "visible"));
      await session.step(31, "And \"WEIGHT\" completion should be visible", () => shouldBe(page, el("\"WEIGHT\" completion"), "visible"));
      await session.step(32, "And \"SEX\" completion should not be visible", () => shouldNotBe(page, el("\"SEX\" completion"), "visible"));
      await session.step(33, "And formula hint should contain text \"CumSum(column:column\"", () => shouldContainText(page, el("formula hint"), "CumSum(column:column"));
      await session.step(34, "When user accepts the highlighted completion with Enter", () => acceptCompletion(page, "Enter"));
      await session.step(35, "Then formula editor should hold the code \"CumSum(${AGE}\"", () => holdsCode(page, el("formula editor"), "CumSum(${AGE}"));
      await session.step(36, "And completion list should be hidden", () => shouldBe(page, el("completion list"), "hidden"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A ColumnName parameter takes the name in quotes", async () => {
      await session.step(40, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(41, "And user types \"Avg(Column(\" into formula editor", () => typeInto(page, "Avg(Column(", el("formula editor")));
      await session.step(42, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(43, "And \"AGE\" completion should be visible", () => shouldBe(page, el("\"AGE\" completion"), "visible"));
      await session.step(44, "And \"SEX\" completion should be visible", () => shouldBe(page, el("\"SEX\" completion"), "visible"));
      await session.step(45, "And formula hint should contain text \"Column(columnName:ColumnName, tableName:TableName)\"", () => shouldContainText(page, el("formula hint"), "Column(columnName:ColumnName, tableName:TableName)"));
      await session.step(46, "When user accepts the highlighted completion with Enter", () => acceptCompletion(page, "Enter"));
      await session.step(47, "Then formula editor should hold the code \"Avg(Column(\\\"AGE\\\"\"", () => holdsCode(page, el("formula editor"), "Avg(Column(\"AGE\""));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A TableName parameter offers the open tables", async () => {
      await session.step(51, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(52, "And user types \"RowCount(\" into formula editor", () => typeInto(page, "RowCount(", el("formula editor")));
      await session.step(53, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(54, "And \"cars\" completion should be visible", () => shouldBe(page, el("\"cars\" completion"), "visible"));
      await session.step(55, "And \"demog\" completion should be visible", () => shouldBe(page, el("\"demog\" completion"), "visible"));
      await session.step(56, "When user accepts the highlighted completion with Enter", () => acceptCompletion(page, "Enter"));
      await session.step(57, "Then formula editor should hold the code \"RowCount(\\\"cars\\\"\"", () => holdsCode(page, el("formula editor"), "RowCount(\"cars\""));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The column list follows the table the call names", async () => {
      await session.step(61, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(62, "And user types \"Lookup(\" into formula editor", () => typeInto(page, "Lookup(", el("formula editor")));
      await session.step(63, "Then \"cars\" completion should be visible", () => shouldBe(page, el("\"cars\" completion"), "visible"));
      await session.step(64, "When user accepts the highlighted completion with Enter", () => acceptCompletion(page, "Enter"));
      await session.step(65, "And user types \", \" at the caret", () => typeAtCaret(page, ", "));
      await session.step(66, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(67, "And \"price\" completion should be visible", () => shouldBe(page, el("\"price\" completion"), "visible"));
      await session.step(68, "And \"model\" completion should be visible", () => shouldBe(page, el("\"model\" completion"), "visible"));
      await session.step(69, "And \"AGE\" completion should not be visible", () => shouldNotBe(page, el("\"AGE\" completion"), "visible"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The third argument of Value is a table again", async () => {
      await session.step(73, "When user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(74, "And user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(75, "And user types \"Value(\\\"AGE\\\", 0, \" into formula editor", () => typeInto(page, "Value(\"AGE\", 0, ", el("formula editor")));
      await session.step(76, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(77, "And \"cars\" completion should be visible", () => shouldBe(page, el("\"cars\" completion"), "visible"));
      await session.step(78, "And formula hint should contain text \"Value(columnName:ColumnName, row:int, tableName:TableName)\"", () => shouldContainText(page, el("formula hint"), "Value(columnName:ColumnName, row:int, tableName:TableName)"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("\"$\" brings the usual column list back, and a numeric parameter offers nothing", async () => {
      await session.step(82, "When user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(83, "And user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(84, "And user types \"Avg(Column(\" into formula editor", () => typeInto(page, "Avg(Column(", el("formula editor")));
      await session.step(85, "And user types \"$\" at the caret", () => typeAtCaret(page, "$"));
      await session.step(86, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(87, "And formula editor should hold the code \"Avg(Column($\"", () => holdsCode(page, el("formula editor"), "Avg(Column($"));
      await session.step(88, "When user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(89, "And user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(90, "And user types \"Sin(\" into formula editor", () => typeInto(page, "Sin(", el("formula editor")));
      await session.step(91, "Then completion list should be hidden", () => shouldBe(page, el("completion list"), "hidden"));
      await session.step(92, "And formula hint should contain text \"Sin(x:num)\"", () => shouldContainText(page, el("formula hint"), "Sin(x:num)"));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Column and RowCount of another table add their columns", async () => {
      await session.step(96, "When user presses Escape in formula editor", () => pressKeyIn(page, "Escape", el("formula editor")));
      await session.step(97, "And user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(98, "And user types \"mean_age\" into column name input", () => typeInto(page, "mean_age", el("column name input")));
      await session.step(99, "And user pastes \"Avg(Column(\\\"AGE\\\"))\" into formula editor", () => pasteInto(page, "Avg(Column(\"AGE\"))", el("formula editor")));
      await session.step(100, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(101, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(102, "And the table should have a column \"mean_age\"", () => hasColumn(page, "mean_age"));
      await session.step(103, "And \"mean_age\" column should have tag \"formula\" equal to \"Avg(Column(\\\"AGE\\\"))\"", () => columnTag(page, "mean_age", "formula", "Avg(Column(\"AGE\"))"));
      await session.step(104, "And the \"mean_age\" cell of row 1 should be displayed as \"45.69\"", () => displayedInRow(page, "mean_age", 1, "45.69"));
      await session.step(105, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(106, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(107, "When user types \"cars_rows\" into column name input", () => typeInto(page, "cars_rows", el("column name input")));
      await session.step(108, "And user pastes \"RowCount(\\\"cars\\\")\" into formula editor", () => pasteInto(page, "RowCount(\"cars\")", el("formula editor")));
      await session.step(109, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(110, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(111, "And the table should have a column \"cars_rows\"", () => hasColumn(page, "cars_rows"));
      await session.step(112, "And every value of \"cars_rows\" column should match \"^30$\"", () => everyValueMatches(page, "cars_rows", "^30$"));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
