/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/pharmacophores.feature
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
import {columnSemType, columnType, columnsEqual, columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Chem > Analyze > Pharmacophores adds a highlighted molecule column and a boolean column per feature type", () => {
  const session = feature(test, "features/analyze/pharmacophores.feature", import.meta.url);
  test("Chem > Analyze > Pharmacophores adds a highlighted molecule column and a boolean column per feature type", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And user opens a table \"pharm\" with:", () => openTableOf(page, "pharm", [["name","smiles"],["ethanol","CCO"],["benzene","c1ccccc1"],["acid","CC(=O)O"],["ethylamine","CCN"]]), [["name","smiles"],["ethanol","CCO"],["benzene","c1ccccc1"],["acid","CC(=O)O"],["ethylamine","CCN"]]);
    await session.step(27, "And the semantic types of the current table are detected", () => detectTypes(page));
    await run.scenario("The dialog opens with four feature types ticked, and OK adds the highlight column and a column per type", async () => {
      await session.step(30, "When user picks \"Chem > Analyze > Pharmacophores...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Pharmacophores..."));
      await session.step(31, "Then \"Pharmacophores\" dialog should be visible", () => shouldBe(page, el("\"Pharmacophores\" dialog"), "visible"));
      await session.step(32, "And Molecules input in \"Pharmacophores\" dialog should contain text \"smiles\"", () => shouldContainText(page, el("Molecules input in \"Pharmacophores\" dialog"), "smiles"));
      await session.step(33, "And Donor input in \"Pharmacophores\" dialog should be checked", () => shouldBe(page, el("Donor input in \"Pharmacophores\" dialog"), "checked"));
      await session.step(34, "And Acceptor input in \"Pharmacophores\" dialog should be checked", () => shouldBe(page, el("Acceptor input in \"Pharmacophores\" dialog"), "checked"));
      await session.step(35, "And Hydrophobic input in \"Pharmacophores\" dialog should be checked", () => shouldBe(page, el("Hydrophobic input in \"Pharmacophores\" dialog"), "checked"));
      await session.step(36, "And Aromatic input in \"Pharmacophores\" dialog should be checked", () => shouldBe(page, el("Aromatic input in \"Pharmacophores\" dialog"), "checked"));
      await session.step(37, "And Positive input in \"Pharmacophores\" dialog should not be checked", () => shouldNotBe(page, el("Positive input in \"Pharmacophores\" dialog"), "checked"));
      await session.step(38, "And Negative input in \"Pharmacophores\" dialog should not be checked", () => shouldNotBe(page, el("Negative input in \"Pharmacophores\" dialog"), "checked"));
      await session.step(39, "And \"Halogen Bond\" input in \"Pharmacophores\" dialog should not be checked", () => shouldNotBe(page, el("\"Halogen Bond\" input in \"Pharmacophores\" dialog"), "checked"));
      await session.step(40, "When user clicks on OK button in \"Pharmacophores\" dialog", () => clickOn(page, el("OK button in \"Pharmacophores\" dialog")));
      await session.step(41, "Then \"Pharmacophores\" dialog should be hidden", () => shouldBe(page, el("\"Pharmacophores\" dialog"), "hidden"));
      await session.step(42, "And a new column \"Pharmacophore (smiles)\" should have been added", () => newColumnNamed(page, "Pharmacophore (smiles)"));
      await session.step(43, "And 5 new columns should have been added", () => newColumnsCount(page, 5));
      await session.step(44, "And \"Pharmacophore (smiles)\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Pharmacophore (smiles)", "Molecule"));
      await session.step(45, "And \"Pharmacophore (smiles)\" column should hold the same values as \"smiles\" column", () => columnsEqual(page, "Pharmacophore (smiles)", "smiles"));
      await session.step(46, "And \"Donor (smiles)\" column should have type \"bool\"", () => columnType(page, "Donor (smiles)", "bool"));
      await session.step(47, "And the value of \"Donor (smiles)\" column in row 1 should be \"true\"", () => valueInRow(page, "Donor (smiles)", 1, "true"));
      await session.step(48, "And the value of \"Donor (smiles)\" column in row 2 should be \"false\"", () => valueInRow(page, "Donor (smiles)", 2, "false"));
      await session.step(49, "And the value of \"Acceptor (smiles)\" column in row 1 should be \"true\"", () => valueInRow(page, "Acceptor (smiles)", 1, "true"));
      await session.step(50, "And the value of \"Acceptor (smiles)\" column in row 2 should be \"false\"", () => valueInRow(page, "Acceptor (smiles)", 2, "false"));
      await session.step(51, "And the value of \"Aromatic (smiles)\" column in row 1 should be \"false\"", () => valueInRow(page, "Aromatic (smiles)", 1, "false"));
      await session.step(52, "And the value of \"Aromatic (smiles)\" column in row 2 should be \"true\"", () => valueInRow(page, "Aromatic (smiles)", 2, "true"));
      await session.step(53, "And the value of \"Hydrophobic (smiles)\" column in row 2 should be \"true\"", () => valueInRow(page, "Hydrophobic (smiles)", 2, "true"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
      await session.step(55, "When user removes \"Pharmacophore (smiles)\" column", () => removeColumn(page, "Pharmacophore (smiles)"));
      await session.step(56, "And user removes \"Donor (smiles)\" column", () => removeColumn(page, "Donor (smiles)"));
      await session.step(57, "And user removes \"Acceptor (smiles)\" column", () => removeColumn(page, "Acceptor (smiles)"));
      await session.step(58, "And user removes \"Aromatic (smiles)\" column", () => removeColumn(page, "Aromatic (smiles)"));
      await session.step(59, "And user removes \"Hydrophobic (smiles)\" column", () => removeColumn(page, "Hydrophobic (smiles)"));
      await session.step(60, "Then the table should have the columns \"name, smiles\"", () => columnsExactly(page, "name, smiles"));
    });
    await run.scenario("Only the feature types ticked get a column", async () => {
      await session.step(63, "When user picks \"Chem > Analyze > Pharmacophores...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Pharmacophores..."));
      await session.step(64, "And user unchecks Donor input in \"Pharmacophores\" dialog", () => uncheck(page, el("Donor input in \"Pharmacophores\" dialog")));
      await session.step(65, "And user unchecks Acceptor input in \"Pharmacophores\" dialog", () => uncheck(page, el("Acceptor input in \"Pharmacophores\" dialog")));
      await session.step(66, "And user unchecks Hydrophobic input in \"Pharmacophores\" dialog", () => uncheck(page, el("Hydrophobic input in \"Pharmacophores\" dialog")));
      await session.step(67, "And user unchecks Aromatic input in \"Pharmacophores\" dialog", () => uncheck(page, el("Aromatic input in \"Pharmacophores\" dialog")));
      await session.step(68, "And user checks Positive input in \"Pharmacophores\" dialog", () => check(page, el("Positive input in \"Pharmacophores\" dialog")));
      await session.step(69, "And user checks Negative input in \"Pharmacophores\" dialog", () => check(page, el("Negative input in \"Pharmacophores\" dialog")));
      await session.step(70, "And user clicks on OK button in \"Pharmacophores\" dialog", () => clickOn(page, el("OK button in \"Pharmacophores\" dialog")));
      await session.step(71, "Then a new column \"Pharmacophore (smiles)\" should have been added", () => newColumnNamed(page, "Pharmacophore (smiles)"));
      await session.step(72, "And a new column \"Positive (smiles)\" should have been added", () => newColumnNamed(page, "Positive (smiles)"));
      await session.step(73, "And a new column \"Negative (smiles)\" should have been added", () => newColumnNamed(page, "Negative (smiles)"));
      await session.step(74, "And 3 new columns should have been added", () => newColumnsCount(page, 3));
      await session.step(75, "And the value of \"Positive (smiles)\" column in row 3 should be \"false\"", () => valueInRow(page, "Positive (smiles)", 3, "false"));
      await session.step(76, "And the value of \"Positive (smiles)\" column in row 4 should be \"true\"", () => valueInRow(page, "Positive (smiles)", 4, "true"));
      await session.step(77, "And the value of \"Negative (smiles)\" column in row 3 should be \"true\"", () => valueInRow(page, "Negative (smiles)", 3, "true"));
      await session.step(78, "And the value of \"Negative (smiles)\" column in row 4 should be \"false\"", () => valueInRow(page, "Negative (smiles)", 4, "false"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
      await session.step(80, "When user removes \"Pharmacophore (smiles)\" column", () => removeColumn(page, "Pharmacophore (smiles)"));
      await session.step(81, "And user removes \"Positive (smiles)\" column", () => removeColumn(page, "Positive (smiles)"));
      await session.step(82, "And user removes \"Negative (smiles)\" column", () => removeColumn(page, "Negative (smiles)"));
      await session.step(83, "Then the table should have the columns \"name, smiles\"", () => columnsExactly(page, "name, smiles"));
    });
    run.finish();
  });
});
