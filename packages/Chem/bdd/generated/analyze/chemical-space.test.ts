/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/chemical-space.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.chemical-space]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe, shouldOffer} from '@datagrok-libraries/bdd/bindings/common/steps';
import {commandCompleted, newColumnMatching, newColumnsMatching, newestMatchingFilled, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openDatasetRows, sketcherIs, viewHoldsViewers} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Chemical Space over SMILES, V2000 and V3000 molecules", () => {
  const session = feature(test, "features/analyze/chemical-space.feature", import.meta.url);
  test("Chemical Space over SMILES, V2000 and V3000 molecules", {tag: ["@journey", "@realizes:chem.cp.chemical-space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await run.scenario("The dialog opens on the molecule column of smiles-50", async () => {
      await session.step(15, "Given user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
      await session.step(16, "When user picks \"Chem > Analyze > Chemical Space...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Chemical Space..."));
      await session.step(17, "Then \"Chem Space\" dialog should be visible", () => shouldBe(page, el("\"Chem Space\" dialog"), "visible"));
      await session.step(18, "And Column input in \"Chem Space\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Column input in \"Chem Space\" dialog"), "canonical_smiles"));
      await session.step(19, "And Method input in \"Chem Space\" dialog should have value \"UMAP\"", () => shouldHaveValue(page, el("Method input in \"Chem Space\" dialog"), "UMAP"));
      await session.step(20, "And Method input in \"Chem Space\" dialog should offer \"UMAP, t-SNE\"", () => shouldOffer(page, el("Method input in \"Chem Space\" dialog"), "UMAP, t-SNE"));
      await session.step(21, "And Similarity input in \"Chem Space\" dialog should have value \"Tanimoto\"", () => shouldHaveValue(page, el("Similarity input in \"Chem Space\" dialog"), "Tanimoto"));
      await session.step(22, "And Similarity input in \"Chem Space\" dialog should offer \"Tanimoto, Asymmetric, Cosine, Sokal\"", () => shouldOffer(page, el("Similarity input in \"Chem Space\" dialog"), "Tanimoto, Asymmetric, Cosine, Sokal"));
      await session.step(23, "And \"Plot embeddings\" input in \"Chem Space\" dialog should be checked", () => shouldBe(page, el("\"Plot embeddings\" input in \"Chem Space\" dialog"), "checked"));
      await session.step(24, "And \"Cluster MCS\" input in \"Chem Space\" dialog should not be checked", () => shouldNotBe(page, el("\"Cluster MCS\" input in \"Chem Space\" dialog"), "checked"));
      await session.step(25, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("UMAP on SMILES adds the embedding columns and plots them", async () => {
      await session.step(28, "When user clicks on OK button in \"Chem Space\" dialog", () => clickOn(page, el("OK button in \"Chem Space\" dialog")));
      await session.step(29, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(30, "And a new column matching \"^Embed_X_\" should have been added", () => newColumnMatching(page, "^Embed_X_"));
      await session.step(31, "And a new column matching \"^Embed_Y_\" should have been added", () => newColumnMatching(page, "^Embed_Y_"));
      await session.step(32, "And a new column matching \"^Cluster \" should have been added", () => newColumnMatching(page, "^Cluster "));
      await session.step(33, "And the newest column matching \"^Embed_X_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_X_"));
      await session.step(34, "And the newest column matching \"^Embed_Y_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_Y_"));
      await session.step(35, "And the current view should hold at least 2 viewers", () => viewHoldsViewers(page, 2));
      await session.step(36, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(37, "And the table should have 50 rows", () => rowCount(page, 50));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second run with t-SNE adds a second pair of embedding columns", async () => {
      await session.step(41, "When user picks \"Chem > Analyze > Chemical Space...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Chemical Space..."));
      await session.step(42, "And user selects \"t-SNE\" in Method input in \"Chem Space\" dialog", () => selectIn(page, "t-SNE", el("Method input in \"Chem Space\" dialog")));
      await session.step(43, "And user clicks on OK button in \"Chem Space\" dialog", () => clickOn(page, el("OK button in \"Chem Space\" dialog")));
      await session.step(44, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(45, "And 2 new columns matching \"^Embed_[XY]_\" should have been added", () => newColumnsMatching(page, 2, "^Embed_[XY]_"));
      await session.step(46, "And the newest column matching \"^Embed_X_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_X_"));
      await session.step(47, "And the table should have 50 rows", () => rowCount(page, 50));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("UMAP on V2000 molecules from an SDF", async () => {
      await session.step(51, "Given user opens mol1K.sdf dataset keeping the first 100 rows", () => openDatasetRows(page, ds("mol1K.sdf"), 100));
      await session.step(52, "When user picks \"Chem > Analyze > Chemical Space...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Chemical Space..."));
      await session.step(53, "Then \"Chem Space\" dialog should be visible", () => shouldBe(page, el("\"Chem Space\" dialog"), "visible"));
      await session.step(54, "When user clicks on OK button in \"Chem Space\" dialog", () => clickOn(page, el("OK button in \"Chem Space\" dialog")));
      await session.step(55, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(56, "And a new column matching \"^Embed_X_\" should have been added", () => newColumnMatching(page, "^Embed_X_"));
      await session.step(57, "And a new column matching \"^Embed_Y_\" should have been added", () => newColumnMatching(page, "^Embed_Y_"));
      await session.step(58, "And the newest column matching \"^Embed_X_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_X_"));
      await session.step(59, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("UMAP on V3000 molecules from an SDF", async () => {
      await session.step(63, "Given user opens ApprovedDrugs2015 dataset", () => openDataset(page, ds("ApprovedDrugs2015")));
      await session.step(64, "When user picks \"Chem > Analyze > Chemical Space...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Chemical Space..."));
      await session.step(65, "Then \"Chem Space\" dialog should be visible", () => shouldBe(page, el("\"Chem Space\" dialog"), "visible"));
      await session.step(66, "When user clicks on OK button in \"Chem Space\" dialog", () => clickOn(page, el("OK button in \"Chem Space\" dialog")));
      await session.step(67, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(68, "And a new column matching \"^Embed_X_\" should have been added", () => newColumnMatching(page, "^Embed_X_"));
      await session.step(69, "And a new column matching \"^Embed_Y_\" should have been added", () => newColumnMatching(page, "^Embed_Y_"));
      await session.step(70, "And the newest column matching \"^Embed_X_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_X_"));
      await session.step(71, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
