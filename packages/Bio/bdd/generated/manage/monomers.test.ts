/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/manage/monomers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.menu.manage.monomers, bio.menu.manage.match-with-library, bio.op.standardise_monomer_library]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {monomerSketcherReady, resultHoldsMonomer, standardiseLibrary} from '../../bindings/bio-b.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldHaveValue, shouldOffer} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {resultIsList} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {closeCurrentView, openApp, openDataset, switchTableView, viewHoldsViewers, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Monomers, and matching molecules against the library", () => {
  const session = feature(test, "features/manage/monomers.feature", import.meta.url);
  test("Monomers, and matching molecules against the library", {tag: ["@journey", "@realizes:bio.menu.manage.monomers", "@realizes:bio.menu.manage.match-with-library", "@realizes:bio.op.standardise_monomer_library"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(17, "And the Bio package is initialized", () => bioInitialized(page));
    await run.scenario("Bio > Manage > Monomers lists every monomer of the libraries", async () => {
      await session.step(20, "When user picks \"Bio > Manage > Monomers\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomers"));
      await session.step(21, "Then the \"Manage Monomers\" view should be current", () => viewIsCurrent(page, "Manage Monomers"));
      await session.step(22, "And the current view should hold at least 1 viewer", () => viewHoldsViewers(page, 1));
      await session.step(23, "And the table should have a column \"Symbol\"", () => hasColumn(page, "Symbol"));
      await session.step(24, "And the table should have a column \"Polymer Type\"", () => hasColumn(page, "Polymer Type"));
      await session.step(25, "And \"Symbol\" column should have no missing values", () => columnComplete(page, "Symbol"));
      await session.step(26, "And \"Polymer Type\" column should have no missing values", () => columnComplete(page, "Polymer Type"));
      await session.step(27, "And the monomer sketcher of the Manage Monomers view should be ready", () => monomerSketcherReady(page));
      await session.step(28, "When user closes the current view", () => closeCurrentView(page));
    });
    await run.scenario("The Manage Monomer Libraries app opens the manager view", async () => {
      await session.step(31, "Given user opens the \"Manage Monomer Libraries\" app", () => openApp(page, "Manage Monomer Libraries"));
      await session.step(32, "Then the \"Manage Monomer Libraries\" view should be current", () => viewIsCurrent(page, "Manage Monomer Libraries"));
      await session.step(33, "And \"HELMCoreLibrary.json\" checkbox should be visible", () => shouldBe(page, el("\"HELMCoreLibrary.json\" checkbox"), "visible"));
      await session.step(34, "When user closes the current view", () => closeCurrentView(page));
    });
    await run.scenario("Match with Monomer Library offers the three polymer types", async () => {
      await session.step(37, "When user switches to the \"filter_HELM\" table view", () => switchTableView(page, "filter_HELM"));
      await session.step(38, "And user picks \"Bio > Manage > Match with Monomer Library...\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Match with Monomer Library..."));
      await session.step(39, "Then \"Match with Monomer Library\" dialog should be visible", () => shouldBe(page, el("\"Match with Monomer Library\" dialog"), "visible"));
      await session.step(40, "And \"Polymer Type\" input in \"Match with Monomer Library\" dialog should offer \"PEPTIDE, RNA, CHEM\"", () => shouldOffer(page, el("\"Polymer Type\" input in \"Match with Monomer Library\" dialog"), "PEPTIDE, RNA, CHEM"));
      await session.step(41, "And \"Polymer Type\" input in \"Match with Monomer Library\" dialog should have value \"PEPTIDE\"", () => shouldHaveValue(page, el("\"Polymer Type\" input in \"Match with Monomer Library\" dialog"), "PEPTIDE"));
      await session.step(42, "When user clicks on CANCEL button in \"Match with Monomer Library\" dialog", () => clickOn(page, el("CANCEL button in \"Match with Monomer Library\" dialog")));
      await session.step(43, "Then \"Match with Monomer Library\" dialog should be hidden", () => shouldBe(page, el("\"Match with Monomer Library\" dialog"), "hidden"));
    });
    await run.scenario("Standardizing the core library keeps its monomers and their polymer types", async () => {
      await session.step(46, "When user standardises the \"HELMCoreLibrary.json\" monomer library", () => standardiseLibrary(page, "HELMCoreLibrary.json"));
      await session.step(47, "Then the result should be a list of 500 or more items", () => resultIsList(page, 500));
      await session.step(48, "And the result should hold \"PEPTIDE\" monomer \"A\"", () => resultHoldsMonomer(page, "PEPTIDE", "A"));
      await session.step(49, "And the result should hold \"PEPTIDE\" monomer \"meI\"", () => resultHoldsMonomer(page, "PEPTIDE", "meI"));
      await session.step(50, "And the result should hold \"RNA\" monomer \"A\"", () => resultHoldsMonomer(page, "RNA", "A"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
