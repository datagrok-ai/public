/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/append-tables.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column]
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
import {check, clickOn, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {addCalculated, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Append Tables keeps the formula of a calculated column and computes it for the appended rows", () => {
  const session = feature(test, "features/add-new-column/append-tables.feature", import.meta.url);
  test("Append Tables keeps the formula of a calculated column and computes it for the appended rows", {tag: ["@journey", "@realizes:powerpack.dialogs.add-new-column"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens a table \"first\" with:", () => openTableOf(page, "first", [["x"],["1"],["2"]]), [["x"],["1"],["2"]]);
    await session.step(20, "And user adds a calculated column \"dbl\" with formula \"${x} * 2\"", () => addCalculated(page, "dbl", "${x} * 2"));
    await session.step(21, "And user opens a table \"second\" with:", () => openTableOf(page, "second", [["x"],["10"]]), [["x"],["10"]]);
    await run.scenario("Both tables appended: the formula stays and covers the appended rows", async () => {
      await session.step(26, "When user picks \"Data > Append Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Append Tables..."));
      await session.step(27, "Then \"Append Tables\" dialog should be visible", () => shouldBe(page, el("\"Append Tables\" dialog"), "visible"));
      await session.step(28, "When user checks \"first\" checkbox in \"Append Tables\" dialog", () => check(page, el("\"first\" checkbox in \"Append Tables\" dialog")));
      await session.step(29, "And user checks \"second\" checkbox in \"Append Tables\" dialog", () => check(page, el("\"second\" checkbox in \"Append Tables\" dialog")));
      await session.step(30, "And user clicks on OK button in \"Append Tables\" dialog", () => clickOn(page, el("OK button in \"Append Tables\" dialog")));
      await session.step(31, "Then \"Append Tables\" dialog should be hidden", () => shouldBe(page, el("\"Append Tables\" dialog"), "hidden"));
      await session.step(32, "And the open table views should be exactly \"first, second, result\"", () => openTableViewsExactly(page, "first, second, result"));
      await session.step(33, "Given user switches to the \"result\" table view", () => switchTableView(page, "result"));
      await session.step(34, "Then the table should have 3 rows", () => rowCount(page, 3));
      await session.step(35, "And the table should have the columns \"x, dbl\"", () => columnsExactly(page, "x, dbl"));
      await session.step(36, "And \"dbl\" column should have tag \"formula\" equal to \"${x} * 2\"", () => columnTag(page, "dbl", "formula", "${x} * 2"));
      await session.step(37, "And the value of \"dbl\" column in row 1 should be \"2\"", () => valueInRow(page, "dbl", 1, "2"));
      await session.step(38, "And the value of \"dbl\" column in row 3 should be \"20\"", () => valueInRow(page, "dbl", 3, "20"));
      await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("OK stays disabled until a table is checked (core#784)", async () => {
      await session.step(43, "When user picks \"Data > Append Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Append Tables..."));
      await session.step(44, "Then \"Append Tables\" dialog should be visible", () => shouldBe(page, el("\"Append Tables\" dialog"), "visible"));
      await session.step(45, "And \"first\" checkbox in \"Append Tables\" dialog should be unchecked", () => shouldBe(page, el("\"first\" checkbox in \"Append Tables\" dialog"), "unchecked"));
      await session.step(46, "And \"second\" checkbox in \"Append Tables\" dialog should be unchecked", () => shouldBe(page, el("\"second\" checkbox in \"Append Tables\" dialog"), "unchecked"));
      await session.step(47, "And OK button in \"Append Tables\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Append Tables\" dialog"), "disabled"));
      await session.step(48, "When user checks \"second\" checkbox in \"Append Tables\" dialog", () => check(page, el("\"second\" checkbox in \"Append Tables\" dialog")));
      await session.step(49, "Then OK button in \"Append Tables\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Append Tables\" dialog"), "enabled"));
      await session.step(50, "When user unchecks \"second\" checkbox in \"Append Tables\" dialog", () => uncheck(page, el("\"second\" checkbox in \"Append Tables\" dialog")));
      await session.step(51, "Then OK button in \"Append Tables\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Append Tables\" dialog"), "disabled"));
      await session.step(52, "When user clicks on CANCEL button in \"Append Tables\" dialog", () => clickOn(page, el("CANCEL button in \"Append Tables\" dialog")));
      await session.step(53, "Then \"Append Tables\" dialog should be hidden", () => shouldBe(page, el("\"Append Tables\" dialog"), "hidden"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
