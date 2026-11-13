/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/elemental-analysis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.elemental-analysis]
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
import {check, clickOn, shouldBe, shouldContainText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, viewHoldsViewers} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Chem > Analyze commands that append columns to SMILES molecules", () => {
  const session = feature(test, "features/analyze/elemental-analysis.feature", import.meta.url);
  test("Chem > Analyze commands that append columns to SMILES molecules", {tag: ["@journey", "@realizes:chem.cp.elemental-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await run.scenario("The counts of every element are appended for SMILES molecules", async () => {
      await session.step(22, "Given user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
      await session.step(23, "When user picks \"Chem > Analyze > Elemental Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Elemental Analysis..."));
      await session.step(24, "Then \"Elemental Analysis\" dialog should be visible", () => shouldBe(page, el("\"Elemental Analysis\" dialog"), "visible"));
      await session.step(25, "And Molecules input in \"Elemental Analysis\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Elemental Analysis\" dialog"), "canonical_smiles"));
      await session.step(26, "And \"Radar Viewer\" input in \"Elemental Analysis\" dialog should not be checked", () => shouldNotBe(page, el("\"Radar Viewer\" input in \"Elemental Analysis\" dialog"), "checked"));
      await session.step(27, "And \"Radar Grid\" input in \"Elemental Analysis\" dialog should not be checked", () => shouldNotBe(page, el("\"Radar Grid\" input in \"Elemental Analysis\" dialog"), "checked"));
      await session.step(28, "When user clicks on OK button in \"Elemental Analysis\" dialog", () => clickOn(page, el("OK button in \"Elemental Analysis\" dialog")));
      await session.step(29, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(30, "And a new column \"C\" should have been added", () => newColumnNamed(page, "C"));
      await session.step(31, "And a new column \"N\" should have been added", () => newColumnNamed(page, "N"));
      await session.step(32, "And a new column \"O\" should have been added", () => newColumnNamed(page, "O"));
      await session.step(33, "And a new column \"Molecule Charge\" should have been added", () => newColumnNamed(page, "Molecule Charge"));
      await session.step(34, "And \"C\" column should have no missing values", () => columnComplete(page, "C"));
      await session.step(35, "And the table should have 50 rows", () => rowCount(page, 50));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Radar Viewer option adds a viewer", async () => {
      await session.step(39, "Given user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
      await session.step(40, "When user picks \"Chem > Analyze > Elemental Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Elemental Analysis..."));
      await session.step(41, "And user checks \"Radar Viewer\" input in \"Elemental Analysis\" dialog", () => check(page, el("\"Radar Viewer\" input in \"Elemental Analysis\" dialog")));
      await session.step(42, "And user clicks on OK button in \"Elemental Analysis\" dialog", () => clickOn(page, el("OK button in \"Elemental Analysis\" dialog")));
      await session.step(43, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(44, "And a new column \"C\" should have been added", () => newColumnNamed(page, "C"));
      await session.step(45, "And the current view should hold at least 2 viewers", () => viewHoldsViewers(page, 2));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Structural Alerts appends a boolean column for the rule set ticked", async () => {
      await session.step(49, "Given user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
      await session.step(50, "When user picks \"Chem > Analyze > Structural Alerts...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Structural Alerts..."));
      await session.step(51, "Then \"Structural Alerts\" dialog should be visible", () => shouldBe(page, el("\"Structural Alerts\" dialog"), "visible"));
      await session.step(52, "And Molecules input in \"Structural Alerts\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Structural Alerts\" dialog"), "canonical_smiles"));
      await session.step(53, "When user clicks on OK button in \"Structural Alerts\" dialog", () => clickOn(page, el("OK button in \"Structural Alerts\" dialog")));
      await session.step(54, "Then a new column \"PAINS (canonical_smiles)\" should have been added", () => newColumnNamed(page, "PAINS (canonical_smiles)"));
      await session.step(55, "And \"PAINS (canonical_smiles)\" column should have type \"bool\"", () => columnType(page, "PAINS (canonical_smiles)", "bool"));
    });
    await run.scenario("Pharmacophores appends the highlighted molecules and a boolean column per feature type ticked", async () => {
      await session.step(58, "Given user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
      await session.step(59, "When user picks \"Chem > Analyze > Pharmacophores...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Pharmacophores..."));
      await session.step(60, "Then \"Pharmacophores\" dialog should be visible", () => shouldBe(page, el("\"Pharmacophores\" dialog"), "visible"));
      await session.step(61, "And Molecules input in \"Pharmacophores\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Pharmacophores\" dialog"), "canonical_smiles"));
      await session.step(62, "When user clicks on OK button in \"Pharmacophores\" dialog", () => clickOn(page, el("OK button in \"Pharmacophores\" dialog")));
      await session.step(63, "Then 5 new columns should have been added", () => newColumnsCount(page, 5));
      await session.step(64, "And a new column \"Pharmacophore (canonical_smiles)\" should have been added", () => newColumnNamed(page, "Pharmacophore (canonical_smiles)"));
      await session.step(65, "And \"Pharmacophore (canonical_smiles)\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Pharmacophore (canonical_smiles)", "Molecule"));
    });
    run.finish();
  });
});
