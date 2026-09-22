/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/calculate/identifiers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.calculate-identifiers]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {containerRunning} from '../../bindings/docker.js';
import {sameMolecules} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, distinctValues, everyValueBetween, everyValueContains, everyValueMatches, everyValueSameLength, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, tableColumns, tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Map Identifiers and Generate Conformers from the Calculate menu", () => {
  const session = feature(test, "features/calculate/identifiers.feature", import.meta.url);
  test("Map Identifiers and Generate Conformers from the Calculate menu", {tag: ["@journey", "@realizes:chem.cp.calculate-identifiers"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
    await run.scenario("Map Identifiers appends InChI keys for the chosen source", async () => {
      await session.step(20, "Given the \"chem-chem\" container is running", () => containerRunning(page, "chem-chem"));
      await session.step(21, "When user picks \"Chem > Calculate > Map Identifiers...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Map Identifiers..."));
      await session.step(22, "Then \"Map Identifiers\" dialog should be visible", () => shouldBe(page, el("\"Map Identifiers\" dialog"), "visible"));
      await session.step(23, "And Ids input in \"Map Identifiers\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Ids input in \"Map Identifiers\" dialog"), "canonical_smiles"));
      await session.step(24, "And \"From Source\" input in \"Map Identifiers\" dialog should have value \"smiles\"", () => shouldHaveValue(page, el("\"From Source\" input in \"Map Identifiers\" dialog"), "smiles"));
      await session.step(25, "When user selects \"inchi_key\" in \"To Source\" input in \"Map Identifiers\" dialog", () => selectIn(page, "inchi_key", el("\"To Source\" input in \"Map Identifiers\" dialog")));
      await session.step(26, "And user clicks on OK button in \"Map Identifiers\" dialog", () => clickOn(page, el("OK button in \"Map Identifiers\" dialog")));
      await session.step(27, "Then \"Map Identifiers\" dialog should be hidden", () => shouldBe(page, el("\"Map Identifiers\" dialog"), "hidden"));
      await session.step(28, "And the top menu command should have completed", () => commandCompleted(page));
      await session.step(29, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(30, "And a new column \"inchi_key\" should have been added", () => newColumnNamed(page, "inchi_key"));
      await session.step(31, "And \"inchi_key\" column should have no missing values", () => columnComplete(page, "inchi_key"));
      await session.step(32, "And every value of \"inchi_key\" column should match \"^[A-Z]{14}-[A-Z]{10}-[A-Z]$\"", () => everyValueMatches(page, "inchi_key", "^[A-Z]{14}-[A-Z]{10}-[A-Z]$"));
      await session.step(33, "And every value of \"inchi_key\" column should have the same length", () => everyValueSameLength(page, "inchi_key"));
      await session.step(34, "And the table should have 50 rows", () => rowCount(page, 50));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Map Identifiers to smiles returns the same molecules", async () => {
      await session.step(38, "Given the \"chem-chem\" container is running", () => containerRunning(page, "chem-chem"));
      await session.step(39, "When user picks \"Chem > Calculate > Map Identifiers...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Map Identifiers..."));
      await session.step(40, "And user selects \"smiles\" in \"To Source\" input in \"Map Identifiers\" dialog", () => selectIn(page, "smiles", el("\"To Source\" input in \"Map Identifiers\" dialog")));
      await session.step(41, "And user clicks on OK button in \"Map Identifiers\" dialog", () => clickOn(page, el("OK button in \"Map Identifiers\" dialog")));
      await session.step(42, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(43, "And a new column \"smiles\" should have been added", () => newColumnNamed(page, "smiles"));
      await session.step(44, "And \"smiles\" column should have no missing values", () => columnComplete(page, "smiles"));
      await session.step(45, "And every molecule of \"smiles\" column should be the same as in \"canonical_smiles\" column", () => sameMolecules(page, "smiles", "canonical_smiles"));
      await session.step(46, "And the table should have 50 rows", () => rowCount(page, 50));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Generate Conformers builds a conformer table for the dialog's own molecule", async () => {
      await session.step(50, "When user picks \"Chem > Calculate > Generate Conformers...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Generate Conformers..."));
      await session.step(51, "Then \"Generate Conformers\" dialog should be visible", () => shouldBe(page, el("\"Generate Conformers\" dialog"), "visible"));
      await session.step(52, "And \"Num conformers\" input in \"Generate Conformers\" dialog should have value \"50\"", () => shouldHaveValue(page, el("\"Num conformers\" input in \"Generate Conformers\" dialog"), "50"));
      await session.step(53, "And Optimize input in \"Generate Conformers\" dialog should be checked", () => shouldBe(page, el("Optimize input in \"Generate Conformers\" dialog"), "checked"));
      await session.step(54, "And \"RMS threshold\" input in \"Generate Conformers\" dialog should have value \"0.1\"", () => shouldHaveValue(page, el("\"RMS threshold\" input in \"Generate Conformers\" dialog"), "0.1"));
      await session.step(55, "And \"Max attempts\" input in \"Generate Conformers\" dialog should have value \"5000\"", () => shouldHaveValue(page, el("\"Max attempts\" input in \"Generate Conformers\" dialog"), "5000"));
      await session.step(56, "And \"Random seed\" input in \"Generate Conformers\" dialog should have value \"42\"", () => shouldHaveValue(page, el("\"Random seed\" input in \"Generate Conformers\" dialog"), "42"));
      await session.step(57, "When user clicks on OK button in \"Generate Conformers\" dialog", () => clickOn(page, el("OK button in \"Generate Conformers\" dialog")));
      await session.step(58, "Then \"Generate Conformers\" dialog should be hidden", () => shouldBe(page, el("\"Generate Conformers\" dialog"), "hidden"));
      await session.step(59, "And the top menu command should have completed", () => commandCompleted(page));
      await session.step(60, "And table \"conformers\" should be open", () => tableOpen(page, "conformers"));
      await session.step(61, "And table \"conformers\" should have columns \"smiles, molblock, conformer, energy, rmsd\"", () => tableColumns(page, "conformers", "smiles, molblock, conformer, energy, rmsd"));
      await session.step(62, "When user switches to the \"conformers\" table view", () => switchTableView(page, "conformers"));
      await session.step(63, "Then every value of \"smiles\" column should match \"^CCCC$\"", () => everyValueMatches(page, "smiles", "^CCCC$"));
      await session.step(64, "And the value of \"conformer\" column in row 1 should be \"1\"", () => valueInRow(page, "conformer", 1, "1"));
      await session.step(65, "And \"conformer\" column should have no missing values", () => columnComplete(page, "conformer"));
      await session.step(66, "And \"conformer\" column should have at least 2 distinct values", () => distinctValues(page, "conformer", 2));
      await session.step(67, "And every value of \"molblock\" column should contain \"M  END\"", () => everyValueContains(page, "molblock", "M  END"));
      await session.step(68, "And \"molblock\" column should have no missing values", () => columnComplete(page, "molblock"));
      await session.step(69, "And \"energy\" column should have no missing values", () => columnComplete(page, "energy"));
      await session.step(70, "And \"rmsd\" column should have no missing values", () => columnComplete(page, "rmsd"));
      await session.step(71, "And \"rmsd\" column should have at least 2 distinct values", () => distinctValues(page, "rmsd", 2));
      await session.step(72, "And every value of \"rmsd\" column should lie between 0 and 100", () => everyValueBetween(page, "rmsd", 0, 100));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
