/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/compare-columns/compare-columns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {topMenuListsEnabled} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, onlyOfAnySelected, removeColumn, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {call} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Data > Compare Columns counts the rows where two columns agree, selects the others and adds a match column", () => {
  const session = feature(test, "features/compare-columns/compare-columns.feature", import.meta.url);
  test("Data > Compare Columns counts the rows where two columns agree, selects the others and adds a match column", {tag: ["@journey", "@known-failure:44", "@known-failure:111"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens a table \"short\" with:", () => openTableOf(page, "short", [["x"],["1"],["2"]]), [["x"],["1"],["2"]]);
    await session.step(29, "And user opens a table \"twin\" with:", () => openTableOf(page, "twin", [["c"],["1"],["2"],["9"],["9"]]), [["c"],["1"],["2"],["9"],["9"]]);
    await session.step(35, "And user opens a table \"cmp\" with:", () => openTableOf(page, "cmp", [["a","b"],["1","1"],["2","3"],["3","3"],["4","5"]]), [["a","b"],["1","1"],["2","3"],["3","3"],["4","5"]]);
    await run.scenario("The Data menu offers Compare Columns", async () => {
      await session.step(44, "Then the top menu should list \"Data > Compare Columns...\" as enabled", () => topMenuListsEnabled(page, "Data > Compare Columns..."));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true, line: 44});
    await run.scenario("Two columns of one table: the counts, Select picks the mismatched rows, ADD COLUMN adds the match column", async () => {
      await session.step(48, "When user calls \"PowerPack:_compareColumns\" function", () => call(page, "PowerPack:_compareColumns"));
      await session.step(49, "Then \"Compare Columns\" dialog should be visible", () => shouldBe(page, el("\"Compare Columns\" dialog"), "visible"));
      await session.step(50, "When user selects \"cmp\" in Tables input in \"Compare Columns\" dialog", () => selectIn(page, "cmp", el("Tables input in \"Compare Columns\" dialog")));
      await session.step(51, "And user selects \"cmp\" in second choice input in \"Compare Columns\" dialog", () => selectIn(page, "cmp", el("second choice input in \"Compare Columns\" dialog")));
      await session.step(52, "And user selects \"a\" in Columns input in \"Compare Columns\" dialog", () => selectIn(page, "a", el("Columns input in \"Compare Columns\" dialog")));
      await session.step(53, "Then \"Compare Columns\" dialog should contain text \"Matched: 2\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "Matched: 2"));
      await session.step(54, "And \"Compare Columns\" dialog should contain text \"Mismatched: 2\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "Mismatched: 2"));
      await session.step(55, "When user clicks on \"Select\" link in \"Compare Columns\" dialog", () => clickOn(page, el("\"Select\" link in \"Compare Columns\" dialog")));
      await session.step(56, "Then only rows where \"a\" is one of \"2, 4\" should be selected", () => onlyOfAnySelected(page, "a", "2, 4"));
      await session.step(57, "When user clicks on \"Add Column\" button in \"Compare Columns\" dialog", () => clickOn(page, el("\"Add Column\" button in \"Compare Columns\" dialog")));
      await session.step(58, "Then the table should have the columns \"a, b, a == b\"", () => columnsExactly(page, "a, b, a == b"));
      await session.step(59, "And the value of \"a == b\" column in row 1 should be \"true\"", () => valueInRow(page, "a == b", 1, "true"));
      await session.step(60, "And the value of \"a == b\" column in row 2 should be \"false\"", () => valueInRow(page, "a == b", 2, "false"));
      await session.step(61, "And the value of \"a == b\" column in row 3 should be \"true\"", () => valueInRow(page, "a == b", 3, "true"));
      await session.step(62, "And the value of \"a == b\" column in row 4 should be \"false\"", () => valueInRow(page, "a == b", 4, "false"));
      await session.step(63, "And \"Compare Columns\" dialog should be visible", () => shouldBe(page, el("\"Compare Columns\" dialog"), "visible"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
      await session.step(65, "When user clears the row selection", () => clearSelection(page));
      await session.step(66, "Then 0 rows should be selected", () => selectedRowCount(page, 0));
      await session.step(67, "When user clicks on CANCEL button in \"Compare Columns\" dialog", () => clickOn(page, el("CANCEL button in \"Compare Columns\" dialog")));
      await session.step(68, "And user removes \"a == b\" column", () => removeColumn(page, "a == b"));
      await session.step(69, "Then the table should have the columns \"a, b\"", () => columnsExactly(page, "a, b"));
    });
    await run.scenario("A column compared with itself matches in every row", async () => {
      await session.step(72, "When user calls \"PowerPack:_compareColumns\" function", () => call(page, "PowerPack:_compareColumns"));
      await session.step(73, "And user selects \"cmp\" in Tables input in \"Compare Columns\" dialog", () => selectIn(page, "cmp", el("Tables input in \"Compare Columns\" dialog")));
      await session.step(74, "And user selects \"cmp\" in second choice input in \"Compare Columns\" dialog", () => selectIn(page, "cmp", el("second choice input in \"Compare Columns\" dialog")));
      await session.step(75, "And user selects \"b\" in Columns input in \"Compare Columns\" dialog", () => selectIn(page, "b", el("Columns input in \"Compare Columns\" dialog")));
      await session.step(76, "Then \"Compare Columns\" dialog should contain text \"Matched: 2\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "Matched: 2"));
      await session.step(77, "When user selects \"b\" in fourth choice input in \"Compare Columns\" dialog", () => selectIn(page, "b", el("fourth choice input in \"Compare Columns\" dialog")));
      await session.step(78, "Then \"Compare Columns\" dialog should contain text \"Matched: 4\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "Matched: 4"));
      await session.step(79, "And \"Compare Columns\" dialog should contain text \"Mismatched: 0\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "Mismatched: 0"));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
      await session.step(81, "When user clicks on CANCEL button in \"Compare Columns\" dialog", () => clickOn(page, el("CANCEL button in \"Compare Columns\" dialog")));
      await session.step(82, "Then \"Compare Columns\" dialog should be hidden", () => shouldBe(page, el("\"Compare Columns\" dialog"), "hidden"));
    });
    await run.scenario("Columns of tables with different row counts are not compared", async () => {
      await session.step(85, "When user calls \"PowerPack:_compareColumns\" function", () => call(page, "PowerPack:_compareColumns"));
      await session.step(86, "And user selects \"cmp\" in Tables input in \"Compare Columns\" dialog", () => selectIn(page, "cmp", el("Tables input in \"Compare Columns\" dialog")));
      await session.step(87, "And user selects \"short\" in second choice input in \"Compare Columns\" dialog", () => selectIn(page, "short", el("second choice input in \"Compare Columns\" dialog")));
      await session.step(88, "And user selects \"a\" in Columns input in \"Compare Columns\" dialog", () => selectIn(page, "a", el("Columns input in \"Compare Columns\" dialog")));
      await session.step(89, "And user selects \"x\" in fourth choice input in \"Compare Columns\" dialog", () => selectIn(page, "x", el("fourth choice input in \"Compare Columns\" dialog")));
      await session.step(90, "Then \"Compare Columns\" dialog should contain text \"Length mismatch\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "Length mismatch"));
      await session.step(91, "And \"Compare Columns\" dialog should contain text \"cmp: 4 rows\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "cmp: 4 rows"));
      await session.step(92, "And \"Compare Columns\" dialog should contain text \"short: 2 rows\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "short: 2 rows"));
      await session.step(93, "And \"Compare Columns\" dialog should not contain text \"Matched\"", () => shouldNotContainText(page, el("\"Compare Columns\" dialog"), "Matched"));
      await session.step(94, "And \"Add Column\" button in \"Compare Columns\" dialog should be absent", () => shouldBe(page, el("\"Add Column\" button in \"Compare Columns\" dialog"), "absent"));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
      await session.step(96, "When user clicks on CANCEL button in \"Compare Columns\" dialog", () => clickOn(page, el("CANCEL button in \"Compare Columns\" dialog")));
      await session.step(97, "Then \"Compare Columns\" dialog should be hidden", () => shouldBe(page, el("\"Compare Columns\" dialog"), "hidden"));
    });
    await run.scenario("A column of another table of the same length is compared, and ADD COLUMN adds the match column to the first table", async () => {
      await session.step(101, "When user calls \"PowerPack:_compareColumns\" function", () => call(page, "PowerPack:_compareColumns"));
      await session.step(102, "And user selects \"cmp\" in Tables input in \"Compare Columns\" dialog", () => selectIn(page, "cmp", el("Tables input in \"Compare Columns\" dialog")));
      await session.step(103, "And user selects \"twin\" in second choice input in \"Compare Columns\" dialog", () => selectIn(page, "twin", el("second choice input in \"Compare Columns\" dialog")));
      await session.step(104, "And user selects \"a\" in Columns input in \"Compare Columns\" dialog", () => selectIn(page, "a", el("Columns input in \"Compare Columns\" dialog")));
      await session.step(105, "And user selects \"c\" in fourth choice input in \"Compare Columns\" dialog", () => selectIn(page, "c", el("fourth choice input in \"Compare Columns\" dialog")));
      await session.step(106, "Then \"Compare Columns\" dialog should contain text \"Matched: 2\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "Matched: 2"));
      await session.step(107, "And \"Compare Columns\" dialog should contain text \"Mismatched: 2\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "Mismatched: 2"));
      await session.step(108, "When user clicks on \"Add Column\" button in \"Compare Columns\" dialog", () => clickOn(page, el("\"Add Column\" button in \"Compare Columns\" dialog")));
      await session.step(109, "Then an info balloon containing \"Column was added to cmp\" should have been shown", () => infoBalloonText(page, "Column was added to cmp"));
      await session.step(110, "And the table should have the columns \"a, b, a == c\"", () => columnsExactly(page, "a, b, a == c"));
      await session.step(111, "And the value of \"a == c\" column in row 1 should be \"true\"", () => valueInRow(page, "a == c", 1, "true"));
      await session.step(112, "And the value of \"a == c\" column in row 2 should be \"true\"", () => valueInRow(page, "a == c", 2, "true"));
      await session.step(113, "And the value of \"a == c\" column in row 3 should be \"false\"", () => valueInRow(page, "a == c", 3, "false"));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
      await session.step(115, "When user clicks on CANCEL button in \"Compare Columns\" dialog", () => clickOn(page, el("CANCEL button in \"Compare Columns\" dialog")));
      await session.step(116, "And user removes \"a == c\" column", () => removeColumn(page, "a == c"));
      await session.step(117, "Then the table should have the columns \"a, b\"", () => columnsExactly(page, "a, b"));
    }, {knownFailure: true, line: 111});
    run.finish();
  });
});
