/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/compare-columns/compare-columns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [core#1246, public#4475]
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
import {clickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {onlyOfAnySelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Data > Compare Columns counts the rows where two columns agree, selects the others and adds a match column", () => {
  const session = feature(test, "features/compare-columns/compare-columns.feature", import.meta.url);
  test("Data > Compare Columns counts the rows where two columns agree, selects the others and adds a match column", {tag: ["@journey", "@realizes:core#1246", "@realizes:public#4475"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens a table \"short\" with:", () => openTableOf(page, "short", [["x"],["1"],["2"]]), [["x"],["1"],["2"]]);
    await session.step(22, "And user opens a table \"twin\" with:", () => openTableOf(page, "twin", [["c"],["1"],["2"],["9"],["9"]]), [["c"],["1"],["2"],["9"],["9"]]);
    await session.step(28, "And user opens a table \"cmp\" with:", () => openTableOf(page, "cmp", [["a","b"],["1","1"],["2","3"],["3","3"],["4","5"]]), [["a","b"],["1","1"],["2","3"],["3","3"],["4","5"]]);
    await run.scenario("Two columns of one table: the counts, Select picks the mismatched rows, ADD COLUMN adds the match column", async () => {
      await session.step(37, "When user picks \"Data > Compare Columns...\" from the top menu", () => pickFromTopMenu(page, "Data > Compare Columns..."));
      await session.step(38, "Then \"Compare Columns\" dialog should be visible", () => shouldBe(page, el("\"Compare Columns\" dialog"), "visible"));
      await session.step(39, "When user selects \"cmp\" in Tables input in \"Compare Columns\" dialog", () => selectIn(page, "cmp", el("Tables input in \"Compare Columns\" dialog")));
      await session.step(40, "And user selects \"cmp\" in second choice input in \"Compare Columns\" dialog", () => selectIn(page, "cmp", el("second choice input in \"Compare Columns\" dialog")));
      await session.step(41, "And user selects \"a\" in Columns input in \"Compare Columns\" dialog", () => selectIn(page, "a", el("Columns input in \"Compare Columns\" dialog")));
      await session.step(42, "Then \"Compare Columns\" dialog should contain text \"Matched: 2\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "Matched: 2"));
      await session.step(43, "And \"Compare Columns\" dialog should contain text \"Mismatched: 2\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "Mismatched: 2"));
      await session.step(44, "When user clicks on \"Select\" link in \"Compare Columns\" dialog", () => clickOn(page, el("\"Select\" link in \"Compare Columns\" dialog")));
      await session.step(45, "Then only rows where \"a\" is one of \"2, 4\" should be selected", () => onlyOfAnySelected(page, "a", "2, 4"));
      await session.step(46, "When user clicks on \"Add Column\" button in \"Compare Columns\" dialog", () => clickOn(page, el("\"Add Column\" button in \"Compare Columns\" dialog")));
      await session.step(47, "Then a new column \"a == b\" should have been added", () => newColumnNamed(page, "a == b"));
      await session.step(48, "And the value of \"a == b\" column in row 1 should be \"true\"", () => valueInRow(page, "a == b", 1, "true"));
      await session.step(49, "And the value of \"a == b\" column in row 2 should be \"false\"", () => valueInRow(page, "a == b", 2, "false"));
      await session.step(50, "When user clicks on CANCEL button in \"Compare Columns\" dialog", () => clickOn(page, el("CANCEL button in \"Compare Columns\" dialog")));
    });
    await run.scenario("Columns of tables with different row counts are not compared", async () => {
      await session.step(53, "When user picks \"Data > Compare Columns...\" from the top menu", () => pickFromTopMenu(page, "Data > Compare Columns..."));
      await session.step(54, "And user selects \"cmp\" in Tables input in \"Compare Columns\" dialog", () => selectIn(page, "cmp", el("Tables input in \"Compare Columns\" dialog")));
      await session.step(55, "And user selects \"short\" in second choice input in \"Compare Columns\" dialog", () => selectIn(page, "short", el("second choice input in \"Compare Columns\" dialog")));
      await session.step(56, "And user selects \"a\" in Columns input in \"Compare Columns\" dialog", () => selectIn(page, "a", el("Columns input in \"Compare Columns\" dialog")));
      await session.step(57, "And user selects \"x\" in fourth choice input in \"Compare Columns\" dialog", () => selectIn(page, "x", el("fourth choice input in \"Compare Columns\" dialog")));
      await session.step(58, "Then \"Compare Columns\" dialog should contain text \"Length mismatch\"", () => shouldContainText(page, el("\"Compare Columns\" dialog"), "Length mismatch"));
      await session.step(59, "And \"Add Column\" button in \"Compare Columns\" dialog should be absent", () => shouldBe(page, el("\"Add Column\" button in \"Compare Columns\" dialog"), "absent"));
      await session.step(60, "When user clicks on CANCEL button in \"Compare Columns\" dialog", () => clickOn(page, el("CANCEL button in \"Compare Columns\" dialog")));
    });
    await run.scenario("A column of another table of the same length is compared, and ADD COLUMN adds the match column to the first table", async () => {
      await session.step(64, "When user picks \"Data > Compare Columns...\" from the top menu", () => pickFromTopMenu(page, "Data > Compare Columns..."));
      await session.step(65, "And user selects \"cmp\" in Tables input in \"Compare Columns\" dialog", () => selectIn(page, "cmp", el("Tables input in \"Compare Columns\" dialog")));
      await session.step(66, "And user selects \"twin\" in second choice input in \"Compare Columns\" dialog", () => selectIn(page, "twin", el("second choice input in \"Compare Columns\" dialog")));
      await session.step(67, "And user selects \"a\" in Columns input in \"Compare Columns\" dialog", () => selectIn(page, "a", el("Columns input in \"Compare Columns\" dialog")));
      await session.step(68, "And user selects \"c\" in fourth choice input in \"Compare Columns\" dialog", () => selectIn(page, "c", el("fourth choice input in \"Compare Columns\" dialog")));
      await session.step(69, "And user clicks on \"Add Column\" button in \"Compare Columns\" dialog", () => clickOn(page, el("\"Add Column\" button in \"Compare Columns\" dialog")));
      await session.step(70, "Then an info balloon containing \"Column was added to cmp\" should have been shown", () => infoBalloonText(page, "Column was added to cmp"));
      await session.step(71, "And the value of \"a == c\" column in row 1 should be \"true\"", () => valueInRow(page, "a == c", 1, "true"));
      await session.step(72, "And the value of \"a == c\" column in row 3 should be \"false\"", () => valueInRow(page, "a == c", 3, "false"));
    });
    run.finish();
  });
});
