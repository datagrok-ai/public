/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/rows-that-fail.feature
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
import {check, clickOn, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {someValueContains, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Add New Column: the rows a formula fails on are left empty or filled, and their messages kept in a column", () => {
  const session = feature(test, "features/add-new-column/rows-that-fail.feature", import.meta.url);
  test("Add New Column: the rows a formula fails on are left empty or filled, and their messages kept in a column", {tag: ["@journey", "@realizes:powerpack.dialogs.add-new-column"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"dates\" with:", () => openTableOf(page, "dates", [["d","n"],["2024-01-05","12"],["soon","many"],["2024-03-07","30"]]), [["d","n"],["2024-01-05","12"],["soon","many"],["2024-03-07","30"]]);
    await run.scenario("Rows that fail are left empty, their messages go to an error column, and the dialog says so", async () => {
      await session.step(27, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(28, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(29, "When user types \"parsed\" into column name input", () => typeInto(page, "parsed", el("column name input")));
      await session.step(30, "And user types \"DateParse(${d})\" into formula editor", () => typeInto(page, "DateParse(${d})", el("formula editor")));
      await session.step(31, "And user clicks on \"Rows that fail\" icon in \"Add New Column\" dialog", () => clickOn(page, el("\"Rows that fail\" icon in \"Add New Column\" dialog")));
      await session.step(32, "Then \"If a row fails\" input should have value \"Leave empty\"", () => shouldHaveValue(page, el("\"If a row fails\" input"), "Leave empty"));
      await session.step(33, "When user checks \"Error column\" input", () => check(page, el("\"Error column\" input")));
      await session.step(34, "Then \"Column name\" input should have value \"parsed errors\"", () => shouldHaveValue(page, el("\"Column name\" input"), "parsed errors"));
      await session.step(35, "And \"Add New Column\" dialog should contain text \"1 of 3 preview rows failed and will be left empty.\"", () => shouldContainText(page, el("\"Add New Column\" dialog"), "1 of 3 preview rows failed and will be left empty."));
      await session.step(36, "And \"Add New Column\" dialog should contain text \"Messages are in \\\"parsed errors\\\".\"", () => shouldContainText(page, el("\"Add New Column\" dialog"), "Messages are in \"parsed errors\"."));
      await session.step(37, "When user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(38, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(39, "And a warning balloon containing \"1 of 3 rows failed and were left empty\" should have been shown", () => warningBalloonText(page, "1 of 3 rows failed and were left empty"));
      await session.step(40, "And the value of \"parsed\" column in row 2 should be \"\"", () => valueInRow(page, "parsed", 2, ""));
      await session.step(41, "And the value of \"parsed errors\" column in row 1 should be \"\"", () => valueInRow(page, "parsed errors", 1, ""));
      await session.step(42, "And some value of \"parsed errors\" column should contain \"soon\"", () => someValueContains(page, "parsed errors", "soon"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Use value fills the rows that fail, and Change... in the preview opens the same settings", async () => {
      await session.step(46, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(47, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(48, "When user types \"filled\" into column name input", () => typeInto(page, "filled", el("column name input")));
      await session.step(49, "And user types \"ParseInt(${n})\" into formula editor", () => typeInto(page, "ParseInt(${n})", el("formula editor")));
      await session.step(50, "And user clicks on \"Rows that fail\" icon in \"Add New Column\" dialog", () => clickOn(page, el("\"Rows that fail\" icon in \"Add New Column\" dialog")));
      await session.step(51, "And user selects \"Use value\" in \"If a row fails\" input", () => selectIn(page, "Use value", el("\"If a row fails\" input")));
      await session.step(52, "And user enters \"7\" into \"Value\" input", () => enterInto(page, "7", el("\"Value\" input")));
      await session.step(53, "And user clicks on formula editor", () => clickOn(page, el("formula editor")));
      await session.step(54, "Then \"Add New Column\" dialog should contain text \"1 of 3 preview rows failed and will be filled with 7.\"", () => shouldContainText(page, el("\"Add New Column\" dialog"), "1 of 3 preview rows failed and will be filled with 7."));
      await session.step(55, "When user clicks on \"Change...\" link in \"Add New Column\" dialog", () => clickOn(page, el("\"Change...\" link in \"Add New Column\" dialog")));
      await session.step(56, "Then \"If a row fails\" input should have value \"Use value\"", () => shouldHaveValue(page, el("\"If a row fails\" input"), "Use value"));
      await session.step(57, "When user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(58, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(59, "And the value of \"filled\" column in row 2 should be \"7\"", () => valueInRow(page, "filled", 2, "7"));
      await session.step(60, "And the value of \"filled\" column in row 3 should be \"30\"", () => valueInRow(page, "filled", 3, "30"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
