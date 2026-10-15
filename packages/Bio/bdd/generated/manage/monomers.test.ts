/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/manage/monomers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.menu.manage.monomers, bio.menu.manage.match-with-library]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {monomerSketcherReady} from '../../bindings/library-files.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldHaveValue, shouldOffer} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, hasColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {closeCurrentView, openApp, openDataset, openTableOf, packageInstalled, switchTableView, viewHoldsViewers, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Monomers, and matching molecules against the library", () => {
  const session = feature(test, "features/manage/monomers.feature", import.meta.url);
  test("Monomers, and matching molecules against the library", {tag: ["@journey", "@realizes:bio.menu.manage.monomers", "@realizes:bio.menu.manage.match-with-library"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(20, "And the Bio package is initialized", () => bioInitialized(page));
    await run.scenario("Bio > Manage > Monomers lists every monomer of the libraries", async () => {
      await session.step(23, "When user picks \"Bio > Manage > Monomers\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomers"));
      await session.step(24, "Then the \"Manage Monomers\" view should be current", () => viewIsCurrent(page, "Manage Monomers"));
      await session.step(25, "And the current view should hold at least 1 viewer", () => viewHoldsViewers(page, 1));
      await session.step(26, "And the table should have a column \"Symbol\"", () => hasColumn(page, "Symbol"));
      await session.step(27, "And the table should have a column \"Polymer Type\"", () => hasColumn(page, "Polymer Type"));
      await session.step(28, "And \"Symbol\" column should have no missing values", () => columnComplete(page, "Symbol"));
      await session.step(29, "And \"Polymer Type\" column should have no missing values", () => columnComplete(page, "Polymer Type"));
      await session.step(30, "And the monomer sketcher of the Manage Monomers view should be ready", () => monomerSketcherReady(page));
      await session.step(31, "When user closes the current view", () => closeCurrentView(page));
    });
    await run.scenario("The Manage Monomer Libraries app opens the manager view", async () => {
      await session.step(34, "Given user opens the \"Manage Monomer Libraries\" app", () => openApp(page, "Manage Monomer Libraries"));
      await session.step(35, "Then the \"Manage Monomer Libraries\" view should be current", () => viewIsCurrent(page, "Manage Monomer Libraries"));
      await session.step(36, "And \"HELMCoreLibrary.json\" checkbox should be visible", () => shouldBe(page, el("\"HELMCoreLibrary.json\" checkbox"), "visible"));
      await session.step(37, "When user closes the current view", () => closeCurrentView(page));
    });
    await run.scenario("Match with Monomer Library offers the three polymer types", async () => {
      await session.step(40, "When user switches to the \"filter_HELM\" table view", () => switchTableView(page, "filter_HELM"));
      await session.step(41, "And user picks \"Bio > Manage > Match with Monomer Library...\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Match with Monomer Library..."));
      await session.step(42, "Then \"Match with Monomer Library\" dialog should be visible", () => shouldBe(page, el("\"Match with Monomer Library\" dialog"), "visible"));
      await session.step(43, "And \"Polymer Type\" input in \"Match with Monomer Library\" dialog should offer \"PEPTIDE, RNA, CHEM\"", () => shouldOffer(page, el("\"Polymer Type\" input in \"Match with Monomer Library\" dialog"), "PEPTIDE, RNA, CHEM"));
      await session.step(44, "And \"Polymer Type\" input in \"Match with Monomer Library\" dialog should have value \"PEPTIDE\"", () => shouldHaveValue(page, el("\"Polymer Type\" input in \"Match with Monomer Library\" dialog"), "PEPTIDE"));
      await session.step(45, "When user clicks on CANCEL button in \"Match with Monomer Library\" dialog", () => clickOn(page, el("CANCEL button in \"Match with Monomer Library\" dialog")));
      await session.step(46, "Then \"Match with Monomer Library\" dialog should be hidden", () => shouldBe(page, el("\"Match with Monomer Library\" dialog"), "hidden"));
    });
    await run.scenario("Match with Monomer Library names the monomer of each molecule", async () => {
      await session.step(50, "Given the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
      await session.step(51, "And user opens a table \"amino-acids\" with:", () => openTableOf(page, "amino-acids", [["molecule"],["C[C@H](N)C(=O)O"],["NCC(=O)O"]]), [["molecule"],["C[C@H](N)C(=O)O"],["NCC(=O)O"]]);
      await session.step(55, "And the semantic types of the current table are detected", () => detectTypes(page));
      await session.step(56, "Then \"molecule\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molecule", "Molecule"));
      await session.step(57, "When user picks \"Bio > Manage > Match with Monomer Library...\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Match with Monomer Library..."));
      await session.step(58, "Then \"Match with Monomer Library\" dialog should be visible", () => shouldBe(page, el("\"Match with Monomer Library\" dialog"), "visible"));
      await session.step(59, "When user clicks on OK button in \"Match with Monomer Library\" dialog", () => clickOn(page, el("OK button in \"Match with Monomer Library\" dialog")));
      await session.step(60, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(61, "And the table should have a column \"Matched monomer symbol\"", () => hasColumn(page, "Matched monomer symbol"));
      await session.step(62, "And the value of \"Matched monomer symbol\" column in row 1 should be \"A\"", () => valueInRow(page, "Matched monomer symbol", 1, "A"));
      await session.step(63, "And the value of \"Matched monomer symbol\" column in row 2 should be \"G\"", () => valueInRow(page, "Matched monomer symbol", 2, "G"));
      await session.step(64, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
