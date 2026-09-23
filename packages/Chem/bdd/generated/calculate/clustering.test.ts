/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/calculate/clustering.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.calculate-clustering]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {matrixDiagonal} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, distinctValues, fewerDistinctThanRows, someValueDiffers} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, newColumnNamed, newestMatchingFilled, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openDatasetRowsAs, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("BitBIRCH clustering, Cluster MCS and the similarity matrix", () => {
  const session = feature(test, "features/calculate/clustering.feature", import.meta.url);
  test("BitBIRCH clustering, Cluster MCS and the similarity matrix", {tag: ["@journey", "@realizes:chem.cp.calculate-clustering"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("BitBIRCH clustering groups the molecules", async () => {
      await session.step(17, "When user picks \"Chem > Calculate > BitBIRCH Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > BitBIRCH Clustering..."));
      await session.step(18, "Then \"BitBIRCH Clustering\" dialog should be visible", () => shouldBe(page, el("\"BitBIRCH Clustering\" dialog"), "visible"));
      await session.step(19, "And Molecules input in \"BitBIRCH Clustering\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"BitBIRCH Clustering\" dialog"), "canonical_smiles"));
      await session.step(20, "And Threshold input in \"BitBIRCH Clustering\" dialog should have value \"0.55\"", () => shouldHaveValue(page, el("Threshold input in \"BitBIRCH Clustering\" dialog"), "0.55"));
      await session.step(21, "And \"Fingerprint type\" input in \"BitBIRCH Clustering\" dialog should have value \"Morgan\"", () => shouldHaveValue(page, el("\"Fingerprint type\" input in \"BitBIRCH Clustering\" dialog"), "Morgan"));
      await session.step(22, "When user clicks on OK button in \"BitBIRCH Clustering\" dialog", () => clickOn(page, el("OK button in \"BitBIRCH Clustering\" dialog")));
      await session.step(23, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(24, "And a new column \"Cluster (BitBIRCH)\" should have been added", () => newColumnNamed(page, "Cluster (BitBIRCH)"));
      await session.step(25, "And \"Cluster (BitBIRCH)\" column should have no missing values", () => columnComplete(page, "Cluster (BitBIRCH)"));
      await session.step(26, "And \"Cluster (BitBIRCH)\" column should have at least 2 distinct values", () => distinctValues(page, "Cluster (BitBIRCH)", 2));
      await session.step(27, "And \"Cluster (BitBIRCH)\" column should have fewer distinct values than the table has rows", () => fewerDistinctThanRows(page, "Cluster (BitBIRCH)"));
      await session.step(28, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The similarity matrix reads 1 down its diagonal", async () => {
      await session.step(32, "Given user opens smiles dataset keeping the first 50 rows as \"clustering_matrix_subset\"", () => openDatasetRowsAs(page, ds("smiles"), 50, "clustering_matrix_subset"));
      await session.step(33, "When user picks \"Chem > Calculate > Similarity Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Similarity Matrix..."));
      await session.step(34, "Then \"Similarity Matrix\" dialog should be visible", () => shouldBe(page, el("\"Similarity Matrix\" dialog"), "visible"));
      await session.step(35, "And Table input in \"Similarity Matrix\" dialog should have value \"clustering_matrix_subset\"", () => shouldHaveValue(page, el("Table input in \"Similarity Matrix\" dialog"), "clustering_matrix_subset"));
      await session.step(36, "And Molecules input in \"Similarity Matrix\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Similarity Matrix\" dialog"), "canonical_smiles"));
      await session.step(37, "And Symbols input in \"Similarity Matrix\" dialog should contain text \"molregno\"", () => shouldContainText(page, el("Symbols input in \"Similarity Matrix\" dialog"), "molregno"));
      await session.step(38, "When user clicks on OK button in \"Similarity Matrix\" dialog", () => clickOn(page, el("OK button in \"Similarity Matrix\" dialog")));
      await session.step(39, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(40, "And table \"canonical_smiles similarity matrix\" should be open", () => tableOpen(page, "canonical_smiles similarity matrix"));
      await session.step(41, "And the similarity columns of table \"canonical_smiles similarity matrix\" should be symmetric, read 1 on the diagonal and less somewhere off it", () => matrixDiagonal(page, "canonical_smiles similarity matrix", 1));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cluster MCS over the BitBIRCH clusters writes a structure on every row", async () => {
      await session.step(45, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(46, "When user picks \"Chem > Calculate > BitBIRCH Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > BitBIRCH Clustering..."));
      await session.step(47, "And user clicks on OK button in \"BitBIRCH Clustering\" dialog", () => clickOn(page, el("OK button in \"BitBIRCH Clustering\" dialog")));
      await session.step(48, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(49, "And a new column \"Cluster (BitBIRCH)\" should have been added", () => newColumnNamed(page, "Cluster (BitBIRCH)"));
      await session.step(50, "And \"Cluster (BitBIRCH)\" column should have at least 2 distinct values", () => distinctValues(page, "Cluster (BitBIRCH)", 2));
      await session.step(51, "And \"Cluster (BitBIRCH)\" column should have fewer distinct values than the table has rows", () => fewerDistinctThanRows(page, "Cluster (BitBIRCH)"));
      await session.step(52, "When user picks \"Chem > Calculate > Cluster MCS...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Cluster MCS..."));
      await session.step(53, "And user clicks on OK button in \"Cluster MCS\" dialog", () => clickOn(page, el("OK button in \"Cluster MCS\" dialog")));
      await session.step(54, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(55, "And a new column matching \"MCS|mcs|Scaffold\" should have been added", () => newColumnMatching(page, "MCS|mcs|Scaffold"));
      await session.step(56, "And the newest column matching \"MCS|mcs|Scaffold\" should have no missing values", () => newestMatchingFilled(page, "MCS|mcs|Scaffold"));
      await session.step(57, "And \"Cluster MCS\" column should have fewer distinct values than the table has rows", () => fewerDistinctThanRows(page, "Cluster MCS"));
      await session.step(58, "And the table should have 100 rows", () => rowCount(page, 100));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Butina clustering groups the same molecules its own way", async () => {
      await session.step(62, "Given user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
      await session.step(63, "When user picks \"Chem > Analyze > Butina Cluster...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Butina Cluster..."));
      await session.step(64, "Then \"Butina Molecules Clustering\" dialog should be visible", () => shouldBe(page, el("\"Butina Molecules Clustering\" dialog"), "visible"));
      await session.step(65, "And Molecules input in \"Butina Molecules Clustering\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Butina Molecules Clustering\" dialog"), "canonical_smiles"));
      await session.step(66, "And \"Distance Cutoff\" input in \"Butina Molecules Clustering\" dialog should have value \"0.4\"", () => shouldHaveValue(page, el("\"Distance Cutoff\" input in \"Butina Molecules Clustering\" dialog"), "0.4"));
      await session.step(67, "When user clicks on OK button in \"Butina Molecules Clustering\" dialog", () => clickOn(page, el("OK button in \"Butina Molecules Clustering\" dialog")));
      await session.step(68, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(69, "And a new column \"cluster (Butina)\" should have been added", () => newColumnNamed(page, "cluster (Butina)"));
      await session.step(70, "And \"cluster (Butina)\" column should have no missing values", () => columnComplete(page, "cluster (Butina)"));
      await session.step(71, "And \"cluster (Butina)\" column should have at least 2 distinct values", () => distinctValues(page, "cluster (Butina)", 2));
      await session.step(72, "And \"cluster (Butina)\" column should have fewer distinct values than the table has rows", () => fewerDistinctThanRows(page, "cluster (Butina)"));
      await session.step(73, "And the table should have 50 rows", () => rowCount(page, 50));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("BitBIRCH and Butina split the same table differently", async () => {
      await session.step(77, "When user picks \"Chem > Calculate > BitBIRCH Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > BitBIRCH Clustering..."));
      await session.step(78, "And user clicks on OK button in \"BitBIRCH Clustering\" dialog", () => clickOn(page, el("OK button in \"BitBIRCH Clustering\" dialog")));
      await session.step(79, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(80, "And a new column \"Cluster (BitBIRCH)\" should have been added", () => newColumnNamed(page, "Cluster (BitBIRCH)"));
      await session.step(81, "And \"Cluster (BitBIRCH)\" column should have at least 2 distinct values", () => distinctValues(page, "Cluster (BitBIRCH)", 2));
      await session.step(82, "And \"Cluster (BitBIRCH)\" column should have fewer distinct values than the table has rows", () => fewerDistinctThanRows(page, "Cluster (BitBIRCH)"));
      await session.step(83, "And some value of \"Cluster (BitBIRCH)\" column should differ from \"cluster (Butina)\" column in the same row", () => someValueDiffers(page, "Cluster (BitBIRCH)", "cluster (Butina)"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
