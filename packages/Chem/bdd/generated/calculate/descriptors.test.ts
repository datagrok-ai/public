/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/calculate/descriptors.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.calculate-descriptors-docker]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {containerRunning} from '../../bindings/docker.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, expand, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, everyValueBetween} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Descriptors from the chem-chem service", () => {
  const session = feature(test, "features/calculate/descriptors.feature", import.meta.url);
  test("Descriptors from the chem-chem service", {tag: ["@journey", "@realizes:chem.cp.calculate-descriptors-docker"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(10, "And user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
    await run.scenario("The dialog offers the descriptor tree and appends the two chosen descriptors", async () => {
      await session.step(13, "Given the \"chem-chem\" container is running", () => containerRunning(page, "chem-chem"));
      await session.step(14, "When user picks \"Chem > Calculate > Descriptors (RDKit)...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Descriptors (RDKit)..."));
      await session.step(15, "Then \"Chemical Descriptors\" dialog should be visible", () => shouldBe(page, el("\"Chemical Descriptors\" dialog"), "visible"));
      await session.step(16, "And Molecules input in \"Chemical Descriptors\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Chemical Descriptors\" dialog"), "canonical_smiles"));
      await session.step(17, "When user clicks on None label in \"Chemical Descriptors\" dialog", () => clickOn(page, el("None label in \"Chemical Descriptors\" dialog")));
      await session.step(18, "Then \"Chemical Descriptors\" dialog should contain text \"0 checked\"", () => shouldContainText(page, el("\"Chemical Descriptors\" dialog"), "0 checked"));
      await session.step(19, "When user expands Descriptors tree node in \"Chemical Descriptors\" dialog", () => expand(page, el("Descriptors tree node in \"Chemical Descriptors\" dialog")));
      await session.step(20, "And user checks MolWt tree node in \"Chemical Descriptors\" dialog", () => check(page, el("MolWt tree node in \"Chemical Descriptors\" dialog")));
      await session.step(21, "And user expands Crippen tree node in \"Chemical Descriptors\" dialog", () => expand(page, el("Crippen tree node in \"Chemical Descriptors\" dialog")));
      await session.step(22, "And user checks MolLogP tree node in \"Chemical Descriptors\" dialog", () => check(page, el("MolLogP tree node in \"Chemical Descriptors\" dialog")));
      await session.step(23, "Then \"Chemical Descriptors\" dialog should contain text \"2 checked\"", () => shouldContainText(page, el("\"Chemical Descriptors\" dialog"), "2 checked"));
      await session.step(24, "When user clicks on OK button in \"Chemical Descriptors\" dialog", () => clickOn(page, el("OK button in \"Chemical Descriptors\" dialog")));
      await session.step(25, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(26, "And 2 new columns should have been added", () => newColumnsCount(page, 2));
      await session.step(27, "And a new column \"MolWt\" should have been added", () => newColumnNamed(page, "MolWt"));
      await session.step(28, "And a new column \"MolLogP\" should have been added", () => newColumnNamed(page, "MolLogP"));
      await session.step(29, "And \"MolWt\" column should have no missing values", () => columnComplete(page, "MolWt"));
      await session.step(30, "And every value of \"MolWt\" column should lie between 1 and 5000", () => everyValueBetween(page, "MolWt", 1, 5000));
      await session.step(31, "And \"MolLogP\" column should have no missing values", () => columnComplete(page, "MolLogP"));
      await session.step(32, "And the table should have 50 rows", () => rowCount(page, 50));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
