/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/curate-mutate.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.transform-curate-mutate]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {changedMolecules, noMoreFragments, noSubstructure, sameMolecules, someMoleculesChanged} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, distinctValues} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Curate and Mutate from the Transform menu", () => {
  const session = feature(test, "features/transform/curate-mutate.feature", import.meta.url);
  test("Curate and Mutate from the Transform menu", {tag: ["@journey", "@realizes:chem.cp.transform-curate-mutate"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And user opens chem_standards dataset", () => openDataset(page, ds("chem_standards")));
    await run.scenario("Curate with its default options standardises the salts and leaves the parents alone", async () => {
      await session.step(23, "When user picks \"Chem > Transform > Curate...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Curate..."));
      await session.step(24, "Then \"Curate\" dialog should be visible", () => shouldBe(page, el("\"Curate\" dialog"), "visible"));
      await session.step(25, "And Molecules input in \"Curate\" dialog should contain text \"smiles\"", () => shouldContainText(page, el("Molecules input in \"Curate\" dialog"), "smiles"));
      await session.step(26, "And Normalization input in \"Curate\" dialog should be checked", () => shouldBe(page, el("Normalization input in \"Curate\" dialog"), "checked"));
      await session.step(27, "And Reionization input in \"Curate\" dialog should be checked", () => shouldBe(page, el("Reionization input in \"Curate\" dialog"), "checked"));
      await session.step(28, "And Neutralization input in \"Curate\" dialog should be checked", () => shouldBe(page, el("Neutralization input in \"Curate\" dialog"), "checked"));
      await session.step(29, "And \"Main fragment\" input in \"Curate\" dialog should be checked", () => shouldBe(page, el("\"Main fragment\" input in \"Curate\" dialog"), "checked"));
      await session.step(30, "And Kekulization input in \"Curate\" dialog should not be checked", () => shouldNotBe(page, el("Kekulization input in \"Curate\" dialog"), "checked"));
      await session.step(31, "And Tautomerization input in \"Curate\" dialog should not be checked", () => shouldNotBe(page, el("Tautomerization input in \"Curate\" dialog"), "checked"));
      await session.step(32, "When user clicks on OK button in \"Curate\" dialog", () => clickOn(page, el("OK button in \"Curate\" dialog")));
      await session.step(33, "Then \"Curate\" dialog should be hidden", () => shouldBe(page, el("\"Curate\" dialog"), "hidden"));
      await session.step(34, "And the top menu command should have completed", () => commandCompleted(page));
      await session.step(35, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(36, "And a new column \"curated_molecule\" should have been added", () => newColumnNamed(page, "curated_molecule"));
      await session.step(37, "And \"curated_molecule\" column should have no missing values", () => columnComplete(page, "curated_molecule"));
      await session.step(38, "And \"curated_molecule\" column should have semantic type \"Molecule\"", () => columnSemType(page, "curated_molecule", "Molecule"));
      await session.step(39, "And no molecule of \"curated_molecule\" column should have more fragments than in \"smiles\" column", () => noMoreFragments(page, "curated_molecule", "smiles"));
      await session.step(40, "And some but not all molecules of \"curated_molecule\" column should differ from \"smiles\" column", () => someMoleculesChanged(page, "curated_molecule", "smiles"));
      await session.step(41, "And 9 molecules of \"curated_molecule\" column should differ from \"smiles\" column", () => changedMolecules(page, 9, "curated_molecule", "smiles"));
      await session.step(42, "And no molecule of \"curated_molecule\" column should contain \"OC(=O)CCC(=O)O\"", () => noSubstructure(page, "curated_molecule", "OC(=O)CCC(=O)O"));
      await session.step(43, "And the table should have 14 rows", () => rowCount(page, 14));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second run with kekulization joins its own column", async () => {
      await session.step(47, "When user picks \"Chem > Transform > Curate...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Curate..."));
      await session.step(48, "And user checks Kekulization input in \"Curate\" dialog", () => check(page, el("Kekulization input in \"Curate\" dialog")));
      await session.step(49, "And user clicks on OK button in \"Curate\" dialog", () => clickOn(page, el("OK button in \"Curate\" dialog")));
      await session.step(50, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(51, "And a new column \"curated_molecule (2)\" should have been added", () => newColumnNamed(page, "curated_molecule (2)"));
      await session.step(52, "And \"curated_molecule (2)\" column should have no missing values", () => columnComplete(page, "curated_molecule (2)"));
      await session.step(53, "And every molecule of \"curated_molecule (2)\" column should be the same as in \"curated_molecule\" column", () => sameMolecules(page, "curated_molecule (2)", "curated_molecule"));
      await session.step(54, "And the table should have 14 rows", () => rowCount(page, 14));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Mutate returns a hundred molecules", async () => {
      await session.step(58, "When user picks \"Chem > Transform > Mutate...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Mutate..."));
      await session.step(59, "Then \"Mutate\" dialog should be visible", () => shouldBe(page, el("\"Mutate\" dialog"), "visible"));
      await session.step(60, "And Steps input in \"Mutate\" dialog should have value \"1\"", () => shouldHaveValue(page, el("Steps input in \"Mutate\" dialog"), "1"));
      await session.step(61, "And Randomize input in \"Mutate\" dialog should be checked", () => shouldBe(page, el("Randomize input in \"Mutate\" dialog"), "checked"));
      await session.step(62, "And \"Max random results\" input in \"Mutate\" dialog should have value \"100\"", () => shouldHaveValue(page, el("\"Max random results\" input in \"Mutate\" dialog"), "100"));
      await session.step(63, "When user clicks on OK button in \"Mutate\" dialog", () => clickOn(page, el("OK button in \"Mutate\" dialog")));
      await session.step(64, "Then \"Mutate\" dialog should be hidden", () => shouldBe(page, el("\"Mutate\" dialog"), "hidden"));
      await session.step(65, "And the top menu command should have completed", () => commandCompleted(page));
      await session.step(66, "And table \"mutations\" should be open", () => tableOpen(page, "mutations"));
      await session.step(67, "When user switches to the \"mutations\" table view", () => switchTableView(page, "mutations"));
      await session.step(68, "Then the table should have 100 rows", () => rowCount(page, 100));
      await session.step(69, "And \"mutations\" column should have no missing values", () => columnComplete(page, "mutations"));
      await session.step(70, "And \"mutations\" column should have semantic type \"Molecule\"", () => columnSemType(page, "mutations", "Molecule"));
      await session.step(71, "And \"mutations\" column should have at least 20 distinct values", () => distinctValues(page, "mutations", 20));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Mutate with two steps still returns a hundred molecules", async () => {
      await session.step(75, "When user picks \"Chem > Transform > Mutate...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Mutate..."));
      await session.step(76, "And user enters \"2\" into Steps input in \"Mutate\" dialog", () => enterInto(page, "2", el("Steps input in \"Mutate\" dialog")));
      await session.step(77, "And user clicks on OK button in \"Mutate\" dialog", () => clickOn(page, el("OK button in \"Mutate\" dialog")));
      await session.step(78, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(79, "And table \"mutations\" should be open", () => tableOpen(page, "mutations"));
      await session.step(80, "When user switches to the \"mutations\" table view", () => switchTableView(page, "mutations"));
      await session.step(81, "Then the table should have 100 rows", () => rowCount(page, 100));
      await session.step(82, "And \"mutations\" column should have no missing values", () => columnComplete(page, "mutations"));
      await session.step(83, "And \"mutations\" column should have at least 20 distinct values", () => distinctValues(page, "mutations", 20));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
