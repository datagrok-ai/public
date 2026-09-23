/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/notation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.transform-notation-roundtrip]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {sameFlatMolecules, sameMolecules} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, everyValueContains, someValueDiffers} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Convert Notation and Recalculate Coordinates keep the molecules", () => {
  const session = feature(test, "features/transform/notation.feature", import.meta.url);
  test("Convert Notation and Recalculate Coordinates keep the molecules", {tag: ["@journey", "@realizes:chem.cp.transform-notation-roundtrip", "@known-failure", "@GROK-20956"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("Converting to molblock adds a column of V2000 molblocks", async () => {
      await session.step(21, "When user picks \"Chem > Transform > Convert Notation...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Convert Notation..."));
      await session.step(22, "Then Molecules input in \"Convert Notation\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Convert Notation\" dialog"), "canonical_smiles"));
      await session.step(23, "And \"Target Notation\" input in \"Convert Notation\" dialog should have value \"smiles\"", () => shouldHaveValue(page, el("\"Target Notation\" input in \"Convert Notation\" dialog"), "smiles"));
      await session.step(24, "And Overwrite input in \"Convert Notation\" dialog should not be checked", () => shouldNotBe(page, el("Overwrite input in \"Convert Notation\" dialog"), "checked"));
      await session.step(25, "And Join input in \"Convert Notation\" dialog should be checked", () => shouldBe(page, el("Join input in \"Convert Notation\" dialog"), "checked"));
      await session.step(26, "When user selects \"molblock\" in \"Target Notation\" input in \"Convert Notation\" dialog", () => selectIn(page, "molblock", el("\"Target Notation\" input in \"Convert Notation\" dialog")));
      await session.step(27, "And user clicks on OK button in \"Convert Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Notation\" dialog")));
      await session.step(28, "Then 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(29, "And a new column \"canonical_smiles_molblock\" should have been added", () => newColumnNamed(page, "canonical_smiles_molblock"));
      await session.step(30, "And \"canonical_smiles_molblock\" column should have no missing values", () => columnComplete(page, "canonical_smiles_molblock"));
      await session.step(31, "And every value of \"canonical_smiles_molblock\" column should contain \"V2000\"", () => everyValueContains(page, "canonical_smiles_molblock", "V2000"));
      await session.step(32, "And every value of \"canonical_smiles_molblock\" column should contain \"M  END\"", () => everyValueContains(page, "canonical_smiles_molblock", "M  END"));
      await session.step(33, "And \"canonical_smiles_molblock\" column should have semantic type \"Molecule\"", () => columnSemType(page, "canonical_smiles_molblock", "Molecule"));
      await session.step(34, "And every molecule of \"canonical_smiles_molblock\" column should be the same as in \"canonical_smiles\" column ignoring stereochemistry", () => sameFlatMolecules(page, "canonical_smiles_molblock", "canonical_smiles"));
      await session.step(35, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Converting the molblocks back to smiles gives the same molecules", async () => {
      await session.step(39, "When user picks \"Chem > Transform > Convert Notation...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Convert Notation..."));
      await session.step(40, "And user selects \"canonical_smiles_molblock\" in Molecules input in \"Convert Notation\" dialog", () => selectIn(page, "canonical_smiles_molblock", el("Molecules input in \"Convert Notation\" dialog")));
      await session.step(41, "And user selects \"smiles\" in \"Target Notation\" input in \"Convert Notation\" dialog", () => selectIn(page, "smiles", el("\"Target Notation\" input in \"Convert Notation\" dialog")));
      await session.step(42, "And user clicks on OK button in \"Convert Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Notation\" dialog")));
      await session.step(43, "Then a new column \"canonical_smiles_molblock_smiles\" should have been added", () => newColumnNamed(page, "canonical_smiles_molblock_smiles"));
      await session.step(44, "And every molecule of \"canonical_smiles_molblock_smiles\" column should be the same as in \"canonical_smiles\" column ignoring stereochemistry", () => sameFlatMolecules(page, "canonical_smiles_molblock_smiles", "canonical_smiles"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Recalculating coordinates with CoordGen moves the atoms and keeps the molecules", async () => {
      await session.step(48, "When user picks \"Chem > Transform > Recalculate Coordinates...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Recalculate Coordinates..."));
      await session.step(49, "Then Molecules input in \"Recalculate Coordinates\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Recalculate Coordinates\" dialog"), "canonical_smiles"));
      await session.step(50, "And Method input in \"Recalculate Coordinates\" dialog should have value \"OCL\"", () => shouldHaveValue(page, el("Method input in \"Recalculate Coordinates\" dialog"), "OCL"));
      await session.step(51, "And Join input in \"Recalculate Coordinates\" dialog should be checked", () => shouldBe(page, el("Join input in \"Recalculate Coordinates\" dialog"), "checked"));
      await session.step(52, "When user selects \"CoordGen\" in Method input in \"Recalculate Coordinates\" dialog", () => selectIn(page, "CoordGen", el("Method input in \"Recalculate Coordinates\" dialog")));
      await session.step(53, "And user clicks on OK button in \"Recalculate Coordinates\" dialog", () => clickOn(page, el("OK button in \"Recalculate Coordinates\" dialog")));
      await session.step(54, "Then a new column \"canonical_smiles_recalcCoords\" should have been added", () => newColumnNamed(page, "canonical_smiles_recalcCoords"));
      await session.step(55, "And every value of \"canonical_smiles_recalcCoords\" column should contain \"M  END\"", () => everyValueContains(page, "canonical_smiles_recalcCoords", "M  END"));
      await session.step(56, "And some value of \"canonical_smiles_recalcCoords\" column should differ from \"canonical_smiles_molblock\" column in the same row", () => someValueDiffers(page, "canonical_smiles_recalcCoords", "canonical_smiles_molblock"));
      await session.step(57, "And every molecule of \"canonical_smiles_recalcCoords\" column should be the same as in \"canonical_smiles\" column ignoring stereochemistry", () => sameFlatMolecules(page, "canonical_smiles_recalcCoords", "canonical_smiles"));
      await session.step(58, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Converted molecules keep their stereochemistry", async () => {
      await session.step(63, "Then every molecule of \"canonical_smiles_molblock\" column should be the same as in \"canonical_smiles\" column", () => sameMolecules(page, "canonical_smiles_molblock", "canonical_smiles"));
      await session.step(64, "And every molecule of \"canonical_smiles_recalcCoords\" column should be the same as in \"canonical_smiles\" column", () => sameMolecules(page, "canonical_smiles_recalcCoords", "canonical_smiles"));
    }, {knownFailure: true});
    run.finish();
  });
});
