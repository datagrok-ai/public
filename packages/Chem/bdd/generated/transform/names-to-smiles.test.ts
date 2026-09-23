/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/names-to-smiles.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.names-to-smiles]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {rowMolecule} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, setColumnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Names To Smiles over a column of compound names", () => {
  const session = feature(test, "features/transform/names-to-smiles.feature", import.meta.url);
  test("Names To Smiles over a column of compound names", {tag: ["@journey", "@realizes:chem.cp.names-to-smiles", "@known-failure", "@GROK-20955"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await run.scenario("The names become molecules in a canonical_smiles column", async () => {
      await session.step(15, "Given user opens a table \"compound_names\" with:", () => openTableOf(page, "compound_names", [["name","mol"],["aspirin","CCO"],["caffeine","CCC"]]), [["name","mol"],["aspirin","CCO"],["caffeine","CCC"]]);
      await session.step(19, "And user sets the semantic type of \"mol\" column to \"Molecule\"", () => setColumnSemType(page, "mol", "Molecule"));
      await session.step(20, "When user picks \"Chem > Transform > Names To Smiles...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Names To Smiles..."));
      await session.step(21, "Then \"Names To Smiles\" dialog should be visible", () => shouldBe(page, el("\"Names To Smiles\" dialog"), "visible"));
      await session.step(22, "When user clicks on OK button in \"Names To Smiles\" dialog", () => clickOn(page, el("OK button in \"Names To Smiles\" dialog")));
      await session.step(23, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(24, "And a new column \"canonical_smiles\" should have been added", () => newColumnNamed(page, "canonical_smiles"));
      await session.step(25, "And \"canonical_smiles\" column should have semantic type \"Molecule\"", () => columnSemType(page, "canonical_smiles", "Molecule"));
      await session.step(26, "And the molecule in row 1 of \"canonical_smiles\" column should be \"CC(=O)Oc1ccccc1C(=O)O\"", () => rowMolecule(page, 1, "canonical_smiles", "CC(=O)Oc1ccccc1C(=O)O"));
      await session.step(27, "And the molecule in row 2 of \"canonical_smiles\" column should be \"Cn1c(=O)c2c(ncn2C)n(C)c1=O\"", () => rowMolecule(page, 2, "canonical_smiles", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"));
      await session.step(28, "And the table should have 2 rows", () => rowCount(page, 2));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The names go into a column of their own beside an existing canonical_smiles", async () => {
      await session.step(33, "Given user opens a table \"named_molecules\" with:", () => openTableOf(page, "named_molecules", [["name","canonical_smiles"],["aspirin","CCO"],["caffeine","CCC"]]), [["name","canonical_smiles"],["aspirin","CCO"],["caffeine","CCC"]]);
      await session.step(37, "And user sets the semantic type of \"canonical_smiles\" column to \"Molecule\"", () => setColumnSemType(page, "canonical_smiles", "Molecule"));
      await session.step(38, "When user picks \"Chem > Transform > Names To Smiles...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Names To Smiles..."));
      await session.step(39, "And user clicks on OK button in \"Names To Smiles\" dialog", () => clickOn(page, el("OK button in \"Names To Smiles\" dialog")));
      await session.step(40, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(41, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
    }, {knownFailure: true});
    run.finish();
  });
});
