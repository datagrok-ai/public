/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/reactions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.transform-reactions]
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
import {moleculesMatching, someMoleculesChanged} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes, rowMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Reactions from the Transform menu", () => {
  const session = feature(test, "features/transform/reactions.feature", import.meta.url);
  test("Reactions from the Transform menu", {tag: ["@journey", "@realizes:chem.cp.transform-reactions"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("Remove Water and Salts adds the desalted molecules", async () => {
      await session.step(17, "When user picks \"Chem > Transform > Reactions > Remove Water and Salts...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Reactions > Remove Water and Salts..."));
      await session.step(18, "Then Molecules input in \"Remove Water and Salts\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Remove Water and Salts\" dialog"), "canonical_smiles"));
      await session.step(19, "When user clicks on OK button in \"Remove Water and Salts\" dialog", () => clickOn(page, el("OK button in \"Remove Water and Salts\" dialog")));
      await session.step(20, "Then a new column \"Desalted(canonical_smiles)\" should have been added", () => newColumnNamed(page, "Desalted(canonical_smiles)"));
      await session.step(21, "And \"Desalted(canonical_smiles)\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Desalted(canonical_smiles)", "Molecule"));
      await session.step(22, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Transformation with ester hydrolysis changes some molecules", async () => {
      await session.step(25, "When user picks \"Chem > Transform > Reactions > Transformation...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Reactions > Transformation..."));
      await session.step(26, "Then Molecules input in \"Run Reaction\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Run Reaction\" dialog"), "canonical_smiles"));
      await session.step(27, "And \"Remove salts and water\" input in \"Run Reaction\" dialog should be checked", () => shouldBe(page, el("\"Remove salts and water\" input in \"Run Reaction\" dialog"), "checked"));
      await session.step(28, "When user clicks on \"Ester Hydrolysis (Saponification)\" reaction in \"Run Reaction\" dialog", () => clickOn(page, el("\"Ester Hydrolysis (Saponification)\" reaction in \"Run Reaction\" dialog")));
      await session.step(29, "And user clicks on OK button in \"Run Reaction\" dialog", () => clickOn(page, el("OK button in \"Run Reaction\" dialog")));
      await session.step(30, "Then a new column \"Reacted(canonical_smiles)\" should have been added", () => newColumnNamed(page, "Reacted(canonical_smiles)"));
      await session.step(31, "And some but not all molecules of \"Reacted(canonical_smiles)\" column should differ from \"Desalted(canonical_smiles)\" column", () => someMoleculesChanged(page, "Reacted(canonical_smiles)", "Desalted(canonical_smiles)"));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Two-Component Reaction couples the acids with the amines", async () => {
      await session.step(35, "Given user opens a table \"amide_coupling\" with:", () => openTableOf(page, "amide_coupling", [["smiles1","smiles2"],["OC(=O)c1ccccc1","Nc1ccccc1"],["CC(=O)O","NCc1ccccc1"],["CCC(=O)O","NC1CCCCC1"],["CCCC(=O)O","CCCCN"],["OC(=O)C1CCCCC1","Cc1ccc(N)cc1"],["Cc1ccc(cc1)C(=O)O","Nc1ccc(Cl)cc1"],["OC(=O)c1ccc(Cl)cc1","NCCc1ccccc1"],["OC(=O)Cc1ccccc1","Nc1cccnc1"],["OC(=O)c1cccnc1","CCN"],["OC(=O)c1ccco1","COc1ccc(N)cc1"],["Cc1ccccc1","Nc1ccccc1"],["OC(=O)c1ccccc1","COc1ccccc1"]]), [["smiles1","smiles2"],["OC(=O)c1ccccc1","Nc1ccccc1"],["CC(=O)O","NCc1ccccc1"],["CCC(=O)O","NC1CCCCC1"],["CCCC(=O)O","CCCCN"],["OC(=O)C1CCCCC1","Cc1ccc(N)cc1"],["Cc1ccc(cc1)C(=O)O","Nc1ccc(Cl)cc1"],["OC(=O)c1ccc(Cl)cc1","NCCc1ccccc1"],["OC(=O)Cc1ccccc1","Nc1cccnc1"],["OC(=O)c1cccnc1","CCN"],["OC(=O)c1ccco1","COc1ccc(N)cc1"],["Cc1ccccc1","Nc1ccccc1"],["OC(=O)c1ccccc1","COc1ccccc1"]]);
      await session.step(49, "When user picks \"Chem > Transform > Reactions > Two-Component Reaction...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Reactions > Two-Component Reaction..."));
      await session.step(50, "Then \"Reactant 1\" input in \"Two-Component Reaction\" dialog should contain text \"smiles1\"", () => shouldContainText(page, el("\"Reactant 1\" input in \"Two-Component Reaction\" dialog"), "smiles1"));
      await session.step(51, "And \"Combination Mode\" input in \"Two-Component Reaction\" dialog should have value \"pairwise\"", () => shouldHaveValue(page, el("\"Combination Mode\" input in \"Two-Component Reaction\" dialog"), "pairwise"));
      await session.step(52, "When user selects \"smiles2\" in \"Reactant 2\" input in \"Two-Component Reaction\" dialog", () => selectIn(page, "smiles2", el("\"Reactant 2\" input in \"Two-Component Reaction\" dialog")));
      await session.step(53, "And user clicks on \"Amide Coupling\" reaction in \"Two-Component Reaction\" dialog", () => clickOn(page, el("\"Amide Coupling\" reaction in \"Two-Component Reaction\" dialog")));
      await session.step(54, "And user clicks on OK button in \"Two-Component Reaction\" dialog", () => clickOn(page, el("OK button in \"Two-Component Reaction\" dialog")));
      await session.step(55, "Then 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(56, "And a new column \"Product(smiles1+smiles2)\" should have been added", () => newColumnNamed(page, "Product(smiles1+smiles2)"));
      await session.step(57, "And 10 molecules of \"Product(smiles1+smiles2)\" column should contain \"[CX3](=O)[NX3;H1]\"", () => moleculesMatching(page, 10, "Product(smiles1+smiles2)", "[CX3](=O)[NX3;H1]"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Deprotect with its default fragment takes the Fmoc group off a protected amino acid", async () => {
      await session.step(61, "Given user opens a table \"fmoc_amino_acids\" with:", () => openTableOf(page, "fmoc_amino_acids", [["smiles"],["O=C(O)CNC(=O)OCC1c2ccccc2-c2ccccc21"]]), [["smiles"],["O=C(O)CNC(=O)OCC1c2ccccc2-c2ccccc21"]]);
      await session.step(64, "And the semantic types of the current table are detected", () => detectTypes(page));
      await session.step(65, "When user picks \"Chem > Transform > Reactions > Deprotect...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Reactions > Deprotect..."));
      await session.step(66, "Then \"Deprotect\" dialog should be visible", () => shouldBe(page, el("\"Deprotect\" dialog"), "visible"));
      await session.step(67, "And Molecules input in \"Deprotect\" dialog should contain text \"smiles\"", () => shouldContainText(page, el("Molecules input in \"Deprotect\" dialog"), "smiles"));
      await session.step(68, "When user clicks on OK button in \"Deprotect\" dialog", () => clickOn(page, el("OK button in \"Deprotect\" dialog")));
      await session.step(69, "Then a new column \"deprotected(smiles)\" should have been added", () => newColumnNamed(page, "deprotected(smiles)"));
      await session.step(70, "And the molecule in row 1 of \"deprotected(smiles)\" column should be \"NCC(=O)O\"", () => rowMolecule(page, 1, "deprotected(smiles)", "NCC(=O)O"));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
