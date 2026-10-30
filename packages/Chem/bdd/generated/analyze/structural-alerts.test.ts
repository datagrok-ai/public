/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/structural-alerts.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, shouldBe, shouldContainText, shouldNotBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType, columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {removeColumn, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Chem > Analyze > Structural Alerts adds a boolean column per rule set ticked", () => {
  const session = feature(test, "features/analyze/structural-alerts.feature", import.meta.url);
  test("Chem > Analyze > Structural Alerts adds a boolean column per rule set ticked", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And user opens a table \"alerts\" with:", () => openTableOf(page, "alerts", [["name","smiles"],["ethanol","CCO"],["catechol","Oc1ccccc1O"],["quinone","O=C1C=CC(=O)C=C1"]]), [["name","smiles"],["ethanol","CCO"],["catechol","Oc1ccccc1O"],["quinone","O=C1C=CC(=O)C=C1"]]);
    await session.step(26, "And the semantic types of the current table are detected", () => detectTypes(page));
    await run.scenario("The dialog opens with PAINS ticked, and OK adds a PAINS column true for the molecules that raise one", async () => {
      await session.step(29, "When user picks \"Chem > Analyze > Structural Alerts...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Structural Alerts..."));
      await session.step(30, "Then \"Structural Alerts\" dialog should be visible", () => shouldBe(page, el("\"Structural Alerts\" dialog"), "visible"));
      await session.step(31, "And Molecules input in \"Structural Alerts\" dialog should contain text \"smiles\"", () => shouldContainText(page, el("Molecules input in \"Structural Alerts\" dialog"), "smiles"));
      await session.step(32, "And PAINS input in \"Structural Alerts\" dialog should be checked", () => shouldBe(page, el("PAINS input in \"Structural Alerts\" dialog"), "checked"));
      await session.step(33, "And BMS input in \"Structural Alerts\" dialog should not be checked", () => shouldNotBe(page, el("BMS input in \"Structural Alerts\" dialog"), "checked"));
      await session.step(34, "And Glaxo input in \"Structural Alerts\" dialog should not be checked", () => shouldNotBe(page, el("Glaxo input in \"Structural Alerts\" dialog"), "checked"));
      await session.step(35, "When user clicks on OK button in \"Structural Alerts\" dialog", () => clickOn(page, el("OK button in \"Structural Alerts\" dialog")));
      await session.step(36, "Then \"Structural Alerts\" dialog should be hidden", () => shouldBe(page, el("\"Structural Alerts\" dialog"), "hidden"));
      await session.step(37, "And a new column \"PAINS (smiles)\" should have been added", () => newColumnNamed(page, "PAINS (smiles)"));
      await session.step(38, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(39, "And \"PAINS (smiles)\" column should have type \"bool\"", () => columnType(page, "PAINS (smiles)", "bool"));
      await session.step(40, "And the value of \"PAINS (smiles)\" column in row 1 should be \"false\"", () => valueInRow(page, "PAINS (smiles)", 1, "false"));
      await session.step(41, "And the value of \"PAINS (smiles)\" column in row 2 should be \"true\"", () => valueInRow(page, "PAINS (smiles)", 2, "true"));
      await session.step(42, "And the value of \"PAINS (smiles)\" column in row 3 should be \"true\"", () => valueInRow(page, "PAINS (smiles)", 3, "true"));
      await session.step(43, "And the table should have 3 rows", () => rowCount(page, 3));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
      await session.step(45, "When user removes \"PAINS (smiles)\" column", () => removeColumn(page, "PAINS (smiles)"));
      await session.step(46, "Then the table should have the columns \"name, smiles\"", () => columnsExactly(page, "name, smiles"));
    });
    await run.scenario("Only the rule sets ticked get a column, each with its own alerts", async () => {
      await session.step(49, "When user picks \"Chem > Analyze > Structural Alerts...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Structural Alerts..."));
      await session.step(50, "And user unchecks PAINS input in \"Structural Alerts\" dialog", () => uncheck(page, el("PAINS input in \"Structural Alerts\" dialog")));
      await session.step(51, "And user checks BMS input in \"Structural Alerts\" dialog", () => check(page, el("BMS input in \"Structural Alerts\" dialog")));
      await session.step(52, "And user checks Glaxo input in \"Structural Alerts\" dialog", () => check(page, el("Glaxo input in \"Structural Alerts\" dialog")));
      await session.step(53, "And user clicks on OK button in \"Structural Alerts\" dialog", () => clickOn(page, el("OK button in \"Structural Alerts\" dialog")));
      await session.step(54, "Then a new column \"Glaxo (smiles)\" should have been added", () => newColumnNamed(page, "Glaxo (smiles)"));
      await session.step(55, "And a new column \"BMS (smiles)\" should have been added", () => newColumnNamed(page, "BMS (smiles)"));
      await session.step(56, "And 2 new columns should have been added", () => newColumnsCount(page, 2));
      await session.step(57, "And the value of \"BMS (smiles)\" column in row 1 should be \"false\"", () => valueInRow(page, "BMS (smiles)", 1, "false"));
      await session.step(58, "And the value of \"BMS (smiles)\" column in row 2 should be \"true\"", () => valueInRow(page, "BMS (smiles)", 2, "true"));
      await session.step(59, "And the value of \"BMS (smiles)\" column in row 3 should be \"true\"", () => valueInRow(page, "BMS (smiles)", 3, "true"));
      await session.step(60, "And the value of \"Glaxo (smiles)\" column in row 1 should be \"false\"", () => valueInRow(page, "Glaxo (smiles)", 1, "false"));
      await session.step(61, "And the value of \"Glaxo (smiles)\" column in row 2 should be \"false\"", () => valueInRow(page, "Glaxo (smiles)", 2, "false"));
      await session.step(62, "And the value of \"Glaxo (smiles)\" column in row 3 should be \"true\"", () => valueInRow(page, "Glaxo (smiles)", 3, "true"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
      await session.step(64, "When user removes \"BMS (smiles)\" column", () => removeColumn(page, "BMS (smiles)"));
      await session.step(65, "And user removes \"Glaxo (smiles)\" column", () => removeColumn(page, "Glaxo (smiles)"));
      await session.step(66, "Then the table should have the columns \"name, smiles\"", () => columnsExactly(page, "name, smiles"));
    });
    run.finish();
  });
});
