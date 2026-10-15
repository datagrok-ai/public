/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/calculate/clustering.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.calculate-clustering]
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
import {clickOn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openDatasetRowsAs, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("BitBIRCH clustering, Cluster MCS and the similarity matrix", () => {
  const session = feature(test, "features/calculate/clustering.feature", import.meta.url);
  test("BitBIRCH clustering, Cluster MCS and the similarity matrix", {tag: ["@journey", "@realizes:chem.cp.calculate-clustering"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("BitBIRCH clustering adds a cluster column", async () => {
      await session.step(19, "When user picks \"Chem > Calculate > BitBIRCH Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > BitBIRCH Clustering..."));
      await session.step(20, "Then \"BitBIRCH Clustering\" dialog should be visible", () => shouldBe(page, el("\"BitBIRCH Clustering\" dialog"), "visible"));
      await session.step(21, "And Molecules input in \"BitBIRCH Clustering\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"BitBIRCH Clustering\" dialog"), "canonical_smiles"));
      await session.step(22, "And Threshold input in \"BitBIRCH Clustering\" dialog should have value \"0.55\"", () => shouldHaveValue(page, el("Threshold input in \"BitBIRCH Clustering\" dialog"), "0.55"));
      await session.step(23, "And \"Fingerprint type\" input in \"BitBIRCH Clustering\" dialog should have value \"Morgan\"", () => shouldHaveValue(page, el("\"Fingerprint type\" input in \"BitBIRCH Clustering\" dialog"), "Morgan"));
      await session.step(24, "When user clicks on OK button in \"BitBIRCH Clustering\" dialog", () => clickOn(page, el("OK button in \"BitBIRCH Clustering\" dialog")));
      await session.step(25, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(26, "And a new column \"Cluster (BitBIRCH)\" should have been added", () => newColumnNamed(page, "Cluster (BitBIRCH)"));
      await session.step(27, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Similarity Matrix opens the matrix of the molecules", async () => {
      await session.step(30, "Given user opens smiles dataset keeping the first 50 rows as \"clustering_matrix_subset\"", () => openDatasetRowsAs(page, ds("smiles"), 50, "clustering_matrix_subset"));
      await session.step(31, "When user picks \"Chem > Calculate > Similarity Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Similarity Matrix..."));
      await session.step(32, "Then \"Similarity Matrix\" dialog should be visible", () => shouldBe(page, el("\"Similarity Matrix\" dialog"), "visible"));
      await session.step(33, "And Table input in \"Similarity Matrix\" dialog should have value \"clustering_matrix_subset\"", () => shouldHaveValue(page, el("Table input in \"Similarity Matrix\" dialog"), "clustering_matrix_subset"));
      await session.step(34, "And Molecules input in \"Similarity Matrix\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Similarity Matrix\" dialog"), "canonical_smiles"));
      await session.step(35, "And Symbols input in \"Similarity Matrix\" dialog should contain text \"molregno\"", () => shouldContainText(page, el("Symbols input in \"Similarity Matrix\" dialog"), "molregno"));
      await session.step(36, "When user clicks on OK button in \"Similarity Matrix\" dialog", () => clickOn(page, el("OK button in \"Similarity Matrix\" dialog")));
      await session.step(37, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(38, "And table \"canonical_smiles similarity matrix\" should be open", () => tableOpen(page, "canonical_smiles similarity matrix"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cluster MCS over the BitBIRCH clusters adds an MCS column", async () => {
      await session.step(42, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(43, "When user picks \"Chem > Calculate > BitBIRCH Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > BitBIRCH Clustering..."));
      await session.step(44, "And user clicks on OK button in \"BitBIRCH Clustering\" dialog", () => clickOn(page, el("OK button in \"BitBIRCH Clustering\" dialog")));
      await session.step(45, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(46, "And a new column \"Cluster (BitBIRCH)\" should have been added", () => newColumnNamed(page, "Cluster (BitBIRCH)"));
      await session.step(47, "When user picks \"Chem > Calculate > Cluster MCS...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Cluster MCS..."));
      await session.step(48, "And user clicks on OK button in \"Cluster MCS\" dialog", () => clickOn(page, el("OK button in \"Cluster MCS\" dialog")));
      await session.step(49, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(50, "And a new column \"Cluster MCS\" should have been added", () => newColumnNamed(page, "Cluster MCS"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
