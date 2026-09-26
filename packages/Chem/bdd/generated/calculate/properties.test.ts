/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/calculate/properties.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.calculate-properties]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, shouldBe, shouldContainText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnType, everyValueBetween, everyValueMatches, someValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, noNewColumn, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openDatasetRowsAs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Chemical properties, toxicity risks and InChI from the Calculate menu", () => {
  const session = feature(test, "features/calculate/properties.feature", import.meta.url);
  test("Chemical properties, toxicity risks and InChI from the Calculate menu", {tag: ["@journey", "@realizes:chem.cp.calculate-properties"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("Chemical Properties with no calculator ticked asks for one", async () => {
      await session.step(21, "When user picks \"Chem > Calculate > Chemical Properties...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Chemical Properties..."));
      await session.step(22, "Then \"Chemical Properties\" dialog should be visible", () => shouldBe(page, el("\"Chemical Properties\" dialog"), "visible"));
      await session.step(23, "And \"Chemical Properties (OCL)\" calculator in \"Chemical Properties\" dialog should not be checked", () => shouldNotBe(page, el("\"Chemical Properties (OCL)\" calculator in \"Chemical Properties\" dialog"), "checked"));
      await session.step(24, "When user clicks on OK button in \"Chemical Properties\" dialog", () => clickOn(page, el("OK button in \"Chemical Properties\" dialog")));
      await session.step(25, "Then a warning balloon containing \"Please select at least one calculation\" should have been shown", () => warningBalloonText(page, "Please select at least one calculation"));
      await session.step(26, "And no new column should have been added", () => noNewColumn(page));
      await session.step(27, "And the table should have 1000 rows", () => rowCount(page, 1000));
    });
    await run.scenario("Chemical Properties with MW alone adds an MW column", async () => {
      await session.step(30, "When user picks \"Chem > Calculate > Chemical Properties...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Chemical Properties..."));
      await session.step(31, "Then \"Chemical Properties\" dialog should be visible", () => shouldBe(page, el("\"Chemical Properties\" dialog"), "visible"));
      await session.step(32, "And Molecules input in \"Chemical Properties\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Chemical Properties\" dialog"), "canonical_smiles"));
      await session.step(33, "And \"Chemical Properties (OCL)\" calculator in \"Chemical Properties\" dialog should not be checked", () => shouldNotBe(page, el("\"Chemical Properties (OCL)\" calculator in \"Chemical Properties\" dialog"), "checked"));
      await session.step(34, "And MW input in \"Chemical Properties\" dialog should be checked", () => shouldBe(page, el("MW input in \"Chemical Properties\" dialog"), "checked"));
      await session.step(35, "And HBA input in \"Chemical Properties\" dialog should not be checked", () => shouldNotBe(page, el("HBA input in \"Chemical Properties\" dialog"), "checked"));
      await session.step(36, "And \"Molecule charge\" input in \"Chemical Properties\" dialog should not be checked", () => shouldNotBe(page, el("\"Molecule charge\" input in \"Chemical Properties\" dialog"), "checked"));
      await session.step(37, "When user checks \"Chemical Properties (OCL)\" calculator in \"Chemical Properties\" dialog", () => check(page, el("\"Chemical Properties (OCL)\" calculator in \"Chemical Properties\" dialog")));
      await session.step(38, "And user clicks on OK button in \"Chemical Properties\" dialog", () => clickOn(page, el("OK button in \"Chemical Properties\" dialog")));
      await session.step(39, "Then \"Chemical Properties\" dialog should be hidden", () => shouldBe(page, el("\"Chemical Properties\" dialog"), "hidden"));
      await session.step(40, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(41, "And a new column \"MW\" should have been added", () => newColumnNamed(page, "MW"));
      await session.step(42, "And every value of \"MW\" column should lie between 1 and 5000", () => everyValueBetween(page, "MW", 1, 5000));
      await session.step(43, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Chemical Properties with every property adds nine columns", async () => {
      await session.step(47, "When user picks \"Chem > Calculate > Chemical Properties...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Chemical Properties..."));
      await session.step(48, "And user checks \"Chemical Properties (OCL)\" calculator in \"Chemical Properties\" dialog", () => check(page, el("\"Chemical Properties (OCL)\" calculator in \"Chemical Properties\" dialog")));
      await session.step(49, "And user checks HBA input in \"Chemical Properties\" dialog", () => check(page, el("HBA input in \"Chemical Properties\" dialog")));
      await session.step(50, "And user checks HBD input in \"Chemical Properties\" dialog", () => check(page, el("HBD input in \"Chemical Properties\" dialog")));
      await session.step(51, "And user checks \"Log P\" input in \"Chemical Properties\" dialog", () => check(page, el("\"Log P\" input in \"Chemical Properties\" dialog")));
      await session.step(52, "And user checks \"Log S\" input in \"Chemical Properties\" dialog", () => check(page, el("\"Log S\" input in \"Chemical Properties\" dialog")));
      await session.step(53, "And user checks PSA input in \"Chemical Properties\" dialog", () => check(page, el("PSA input in \"Chemical Properties\" dialog")));
      await session.step(54, "And user checks \"Rotatable bonds\" input in \"Chemical Properties\" dialog", () => check(page, el("\"Rotatable bonds\" input in \"Chemical Properties\" dialog")));
      await session.step(55, "And user checks \"Stereo centers\" input in \"Chemical Properties\" dialog", () => check(page, el("\"Stereo centers\" input in \"Chemical Properties\" dialog")));
      await session.step(56, "And user checks \"Molecule charge\" input in \"Chemical Properties\" dialog", () => check(page, el("\"Molecule charge\" input in \"Chemical Properties\" dialog")));
      await session.step(57, "And user clicks on OK button in \"Chemical Properties\" dialog", () => clickOn(page, el("OK button in \"Chemical Properties\" dialog")));
      await session.step(58, "Then 9 new columns should have been added", () => newColumnsCount(page, 9));
      await session.step(59, "And a new column \"MW (2)\" should have been added", () => newColumnNamed(page, "MW (2)"));
      await session.step(60, "And every value of \"HBA\" column should lie between 0 and 1000", () => everyValueBetween(page, "HBA", 0, 1000));
      await session.step(61, "And every value of \"HBD\" column should lie between 0 and 1000", () => everyValueBetween(page, "HBD", 0, 1000));
      await session.step(62, "And every value of \"LogP\" column should lie between -50 and 50", () => everyValueBetween(page, "LogP", -50, 50));
      await session.step(63, "And every value of \"LogS\" column should lie between -50 and 50", () => everyValueBetween(page, "LogS", -50, 50));
      await session.step(64, "And every value of \"PSA\" column should lie between 0 and 1000", () => everyValueBetween(page, "PSA", 0, 1000));
      await session.step(65, "And every value of \"Rotatable bonds\" column should lie between 0 and 1000", () => everyValueBetween(page, "Rotatable bonds", 0, 1000));
      await session.step(66, "And every value of \"Stereo centers\" column should lie between 0 and 1000", () => everyValueBetween(page, "Stereo centers", 0, 1000));
      await session.step(67, "And every value of \"Molecule charge\" column should lie between -20 and 20", () => everyValueBetween(page, "Molecule charge", -20, 20));
      await session.step(68, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("To InchI adds a column of InChI strings", async () => {
      await session.step(72, "When user picks \"Chem > Calculate > To InchI...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > To InchI..."));
      await session.step(73, "Then Molecules input in \"To InchI\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"To InchI\" dialog"), "canonical_smiles"));
      await session.step(74, "When user clicks on OK button in \"To InchI\" dialog", () => clickOn(page, el("OK button in \"To InchI\" dialog")));
      await session.step(75, "Then a new column \"inchi\" should have been added", () => newColumnNamed(page, "inchi"));
      await session.step(76, "And \"inchi\" column should have no missing values", () => columnComplete(page, "inchi"));
      await session.step(77, "And every value of \"inchi\" column should match \"^InChI=1S/\"", () => everyValueMatches(page, "inchi", "^InChI=1S/"));
      await session.step(78, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("To InchI Keys adds a column of InChI keys", async () => {
      await session.step(82, "When user picks \"Chem > Calculate > To InchI Keys...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > To InchI Keys..."));
      await session.step(83, "And user clicks on OK button in \"To InchI Keys\" dialog", () => clickOn(page, el("OK button in \"To InchI Keys\" dialog")));
      await session.step(84, "Then a new column \"inchi_key\" should have been added", () => newColumnNamed(page, "inchi_key"));
      await session.step(85, "And \"inchi_key\" column should have no missing values", () => columnComplete(page, "inchi_key"));
      await session.step(86, "And every value of \"inchi_key\" column should match \"^[A-Z]{14}-[A-Z]{10}-[A-Z]$\"", () => everyValueMatches(page, "inchi_key", "^[A-Z]{14}-[A-Z]{10}-[A-Z]$"));
      await session.step(87, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Toxicity Risks opens with Mutagenicity alone and adds a column per ticked risk", async () => {
      await session.step(91, "Given user opens smiles dataset keeping the first 100 rows as \"smiles_100\"", () => openDatasetRowsAs(page, ds("smiles"), 100, "smiles_100"));
      await session.step(92, "When user picks \"Chem > Calculate > Toxicity Risks...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Toxicity Risks..."));
      await session.step(93, "Then Mutagenicity input in \"Toxicity Risks\" dialog should be checked", () => shouldBe(page, el("Mutagenicity input in \"Toxicity Risks\" dialog"), "checked"));
      await session.step(94, "And Tumorigenicity input in \"Toxicity Risks\" dialog should not be checked", () => shouldNotBe(page, el("Tumorigenicity input in \"Toxicity Risks\" dialog"), "checked"));
      await session.step(95, "And \"Irritating effects\" input in \"Toxicity Risks\" dialog should not be checked", () => shouldNotBe(page, el("\"Irritating effects\" input in \"Toxicity Risks\" dialog"), "checked"));
      await session.step(96, "And \"Reproductive effects\" input in \"Toxicity Risks\" dialog should not be checked", () => shouldNotBe(page, el("\"Reproductive effects\" input in \"Toxicity Risks\" dialog"), "checked"));
      await session.step(97, "When user checks Tumorigenicity input in \"Toxicity Risks\" dialog", () => check(page, el("Tumorigenicity input in \"Toxicity Risks\" dialog")));
      await session.step(98, "And user checks \"Irritating effects\" input in \"Toxicity Risks\" dialog", () => check(page, el("\"Irritating effects\" input in \"Toxicity Risks\" dialog")));
      await session.step(99, "And user checks \"Reproductive effects\" input in \"Toxicity Risks\" dialog", () => check(page, el("\"Reproductive effects\" input in \"Toxicity Risks\" dialog")));
      await session.step(100, "And user clicks on OK button in \"Toxicity Risks\" dialog", () => clickOn(page, el("OK button in \"Toxicity Risks\" dialog")));
      await session.step(101, "Then 4 new columns should have been added", () => newColumnsCount(page, 4));
      await session.step(102, "And a new column \"Mutagenicity\" should have been added", () => newColumnNamed(page, "Mutagenicity"));
      await session.step(103, "And \"Mutagenicity\" column should have type \"string\"", () => columnType(page, "Mutagenicity", "string"));
      await session.step(104, "And \"Mutagenicity\" column should have no missing values", () => columnComplete(page, "Mutagenicity"));
      await session.step(105, "And every value of \"Mutagenicity\" column should match \"^(Unknown|None|Low|High)$\"", () => everyValueMatches(page, "Mutagenicity", "^(Unknown|None|Low|High)$"));
      await session.step(106, "And every value of \"Tumorigenicity\" column should match \"^(Unknown|None|Low|High)$\"", () => everyValueMatches(page, "Tumorigenicity", "^(Unknown|None|Low|High)$"));
      await session.step(107, "And some value of \"Tumorigenicity\" column should match \"^(Low|High)$\"", () => someValueMatches(page, "Tumorigenicity", "^(Low|High)$"));
      await session.step(108, "And every value of \"Irritating effects\" column should match \"^(Unknown|None|Low|High)$\"", () => everyValueMatches(page, "Irritating effects", "^(Unknown|None|Low|High)$"));
      await session.step(109, "And some value of \"Irritating effects\" column should match \"^(Low|High)$\"", () => someValueMatches(page, "Irritating effects", "^(Low|High)$"));
      await session.step(110, "And every value of \"Reproductive effects\" column should match \"^(Unknown|None|Low|High)$\"", () => everyValueMatches(page, "Reproductive effects", "^(Unknown|None|Low|High)$"));
      await session.step(111, "And some value of \"Reproductive effects\" column should match \"^(Low|High)$\"", () => someValueMatches(page, "Reproductive effects", "^(Low|High)$"));
      await session.step(112, "And the table should have 100 rows", () => rowCount(page, 100));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
