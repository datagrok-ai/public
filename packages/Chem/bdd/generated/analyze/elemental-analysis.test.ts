/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/elemental-analysis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.elemental-analysis]
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
import {columnComplete, columnType, everyValueBetween} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, viewHoldsViewers} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Elemental Analysis over SMILES, V2000 and V3000 molecules", () => {
  const session = feature(test, "features/analyze/elemental-analysis.feature", import.meta.url);
  test("Elemental Analysis over SMILES, V2000 and V3000 molecules", {tag: ["@journey", "@realizes:chem.cp.elemental-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And the package autostarts have completed", () => autostartsCompleted(page));
    await run.scenario("The counts of every element are appended for SMILES molecules", async () => {
      await session.step(12, "Given user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
      await session.step(13, "When user picks \"Chem > Analyze > Elemental Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Elemental Analysis..."));
      await session.step(14, "Then \"Elemental Analysis\" dialog should be visible", () => shouldBe(page, el("\"Elemental Analysis\" dialog"), "visible"));
      await session.step(15, "And Molecules input in \"Elemental Analysis\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Elemental Analysis\" dialog"), "canonical_smiles"));
      await session.step(16, "And \"Radar Viewer\" input in \"Elemental Analysis\" dialog should not be checked", () => shouldNotBe(page, el("\"Radar Viewer\" input in \"Elemental Analysis\" dialog"), "checked"));
      await session.step(17, "And \"Radar Grid\" input in \"Elemental Analysis\" dialog should not be checked", () => shouldNotBe(page, el("\"Radar Grid\" input in \"Elemental Analysis\" dialog"), "checked"));
      await session.step(18, "When user clicks on OK button in \"Elemental Analysis\" dialog", () => clickOn(page, el("OK button in \"Elemental Analysis\" dialog")));
      await session.step(19, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(20, "And a new column \"C\" should have been added", () => newColumnNamed(page, "C"));
      await session.step(21, "And a new column \"N\" should have been added", () => newColumnNamed(page, "N"));
      await session.step(22, "And a new column \"O\" should have been added", () => newColumnNamed(page, "O"));
      await session.step(23, "And a new column \"Molecule Charge\" should have been added", () => newColumnNamed(page, "Molecule Charge"));
      await session.step(24, "And \"C\" column should have no missing values", () => columnComplete(page, "C"));
      await session.step(25, "And every value of \"C\" column should lie between 1 and 200", () => everyValueBetween(page, "C", 1, 200));
      await session.step(26, "And \"C\" column should have type \"int\"", () => columnType(page, "C", "int"));
      await session.step(27, "And the table should have 50 rows", () => rowCount(page, 50));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The same counts come out of V2000 molecules", async () => {
      await session.step(31, "Given user opens mol1K.sdf dataset", () => openDataset(page, ds("mol1K.sdf")));
      await session.step(32, "When user picks \"Chem > Analyze > Elemental Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Elemental Analysis..."));
      await session.step(33, "And user clicks on OK button in \"Elemental Analysis\" dialog", () => clickOn(page, el("OK button in \"Elemental Analysis\" dialog")));
      await session.step(34, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(35, "And a new column \"C\" should have been added", () => newColumnNamed(page, "C"));
      await session.step(36, "And a new column \"Molecule Charge\" should have been added", () => newColumnNamed(page, "Molecule Charge"));
      await session.step(37, "And \"C\" column should have no missing values", () => columnComplete(page, "C"));
      await session.step(38, "And every value of \"C\" column should lie between 1 and 200", () => everyValueBetween(page, "C", 1, 200));
      await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The same counts come out of V3000 molecules", async () => {
      await session.step(43, "Given user opens ApprovedDrugs2015 dataset", () => openDataset(page, ds("ApprovedDrugs2015")));
      await session.step(44, "When user picks \"Chem > Analyze > Elemental Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Elemental Analysis..."));
      await session.step(45, "And user clicks on OK button in \"Elemental Analysis\" dialog", () => clickOn(page, el("OK button in \"Elemental Analysis\" dialog")));
      await session.step(46, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(47, "And a new column \"C\" should have been added", () => newColumnNamed(page, "C"));
      await session.step(48, "And a new column \"Molecule Charge\" should have been added", () => newColumnNamed(page, "Molecule Charge"));
      await session.step(49, "And \"C\" column should have no missing values", () => columnComplete(page, "C"));
      await session.step(50, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Radar Viewer option adds a viewer", async () => {
      await session.step(54, "Given user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
      await session.step(55, "When user picks \"Chem > Analyze > Elemental Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Elemental Analysis..."));
      await session.step(56, "And user checks \"Radar Viewer\" input in \"Elemental Analysis\" dialog", () => check(page, el("\"Radar Viewer\" input in \"Elemental Analysis\" dialog")));
      await session.step(57, "And user clicks on OK button in \"Elemental Analysis\" dialog", () => clickOn(page, el("OK button in \"Elemental Analysis\" dialog")));
      await session.step(58, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(59, "And a new column \"C\" should have been added", () => newColumnNamed(page, "C"));
      await session.step(60, "And the current view should hold at least 2 viewers", () => viewHoldsViewers(page, 2));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("SMARTS patterns are counted the same way", async () => {
      await session.step(64, "Given user opens ex-smarts dataset", () => openDataset(page, ds("ex-smarts")));
      await session.step(65, "When user picks \"Chem > Analyze > Elemental Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Elemental Analysis..."));
      await session.step(66, "Then \"Elemental Analysis\" dialog should be visible", () => shouldBe(page, el("\"Elemental Analysis\" dialog"), "visible"));
      await session.step(67, "When user clicks on OK button in \"Elemental Analysis\" dialog", () => clickOn(page, el("OK button in \"Elemental Analysis\" dialog")));
      await session.step(68, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(69, "And a new column \"Molecule Charge\" should have been added", () => newColumnNamed(page, "Molecule Charge"));
      await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
