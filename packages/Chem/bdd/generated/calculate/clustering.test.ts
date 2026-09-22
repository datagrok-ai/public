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
import {columnComplete, distinctValues, fewerDistinctThanRows} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, newColumnNamed, newestMatchingFilled, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openDatasetRowsAs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("BitBIRCH clustering, Cluster MCS and the similarity matrix", () => {
  const session = feature(test, "features/calculate/clustering.feature", import.meta.url);
  test("BitBIRCH clustering, Cluster MCS and the similarity matrix", {tag: ["@journey", "@realizes:chem.cp.calculate-clustering"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(12, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("BitBIRCH clustering groups the molecules", async () => {
      await session.step(15, "When user picks \"Chem > Calculate > BitBIRCH Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > BitBIRCH Clustering..."));
      await session.step(16, "Then \"BitBIRCH Clustering\" dialog should be visible", () => shouldBe(page, el("\"BitBIRCH Clustering\" dialog"), "visible"));
      await session.step(17, "And Molecules input in \"BitBIRCH Clustering\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"BitBIRCH Clustering\" dialog"), "canonical_smiles"));
      await session.step(18, "And Threshold input in \"BitBIRCH Clustering\" dialog should have value \"0.55\"", () => shouldHaveValue(page, el("Threshold input in \"BitBIRCH Clustering\" dialog"), "0.55"));
      await session.step(19, "And \"Fingerprint type\" input in \"BitBIRCH Clustering\" dialog should have value \"Morgan\"", () => shouldHaveValue(page, el("\"Fingerprint type\" input in \"BitBIRCH Clustering\" dialog"), "Morgan"));
      await session.step(20, "When user clicks on OK button in \"BitBIRCH Clustering\" dialog", () => clickOn(page, el("OK button in \"BitBIRCH Clustering\" dialog")));
      await session.step(21, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(22, "And a new column \"Cluster (BitBIRCH)\" should have been added", () => newColumnNamed(page, "Cluster (BitBIRCH)"));
      await session.step(23, "And \"Cluster (BitBIRCH)\" column should have no missing values", () => columnComplete(page, "Cluster (BitBIRCH)"));
      await session.step(24, "And \"Cluster (BitBIRCH)\" column should have at least 2 distinct values", () => distinctValues(page, "Cluster (BitBIRCH)", 2));
      await session.step(25, "And \"Cluster (BitBIRCH)\" column should have fewer distinct values than the table has rows", () => fewerDistinctThanRows(page, "Cluster (BitBIRCH)"));
      await session.step(26, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(27, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cluster MCS writes a structure for every row", async () => {
      await session.step(30, "When user picks \"Chem > Calculate > Cluster MCS...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Cluster MCS..."));
      await session.step(31, "Then \"Cluster MCS\" dialog should be visible", () => shouldBe(page, el("\"Cluster MCS\" dialog"), "visible"));
      await session.step(32, "And Molecules input in \"Cluster MCS\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Cluster MCS\" dialog"), "canonical_smiles"));
      await session.step(33, "When user clicks on OK button in \"Cluster MCS\" dialog", () => clickOn(page, el("OK button in \"Cluster MCS\" dialog")));
      await session.step(34, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(35, "And a new column matching \"MCS|mcs|Scaffold\" should have been added", () => newColumnMatching(page, "MCS|mcs|Scaffold"));
      await session.step(36, "And the newest column matching \"MCS|mcs|Scaffold\" should have no missing values", () => newestMatchingFilled(page, "MCS|mcs|Scaffold"));
      await session.step(37, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The similarity matrix reads 1 down its diagonal", async () => {
      await session.step(41, "Given user opens smiles dataset keeping the first 50 rows as \"clustering_matrix_subset\"", () => openDatasetRowsAs(page, ds("smiles"), 50, "clustering_matrix_subset"));
      await session.step(42, "When user picks \"Chem > Calculate > Similarity Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Similarity Matrix..."));
      await session.step(43, "Then \"Similarity Matrix\" dialog should be visible", () => shouldBe(page, el("\"Similarity Matrix\" dialog"), "visible"));
      await session.step(44, "And Table input in \"Similarity Matrix\" dialog should have value \"clustering_matrix_subset\"", () => shouldHaveValue(page, el("Table input in \"Similarity Matrix\" dialog"), "clustering_matrix_subset"));
      await session.step(45, "And Molecules input in \"Similarity Matrix\" dialog should contain text \"canonical_smiles\"", () => shouldContainText(page, el("Molecules input in \"Similarity Matrix\" dialog"), "canonical_smiles"));
      await session.step(46, "And Symbols input in \"Similarity Matrix\" dialog should contain text \"molregno\"", () => shouldContainText(page, el("Symbols input in \"Similarity Matrix\" dialog"), "molregno"));
      await session.step(47, "When user clicks on OK button in \"Similarity Matrix\" dialog", () => clickOn(page, el("OK button in \"Similarity Matrix\" dialog")));
      await session.step(48, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(49, "And table \"canonical_smiles similarity matrix\" should be open", () => tableOpen(page, "canonical_smiles similarity matrix"));
      await session.step(50, "And the similarity columns of table \"canonical_smiles similarity matrix\" should be symmetric, read 1 on the diagonal and less somewhere off it", () => matrixDiagonal(page, "canonical_smiles similarity matrix", 1));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cluster MCS over real clusters writes a scaffold shared inside each of them", async () => {
      await session.step(54, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(55, "When user picks \"Chem > Calculate > BitBIRCH Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > BitBIRCH Clustering..."));
      await session.step(56, "And user clicks on OK button in \"BitBIRCH Clustering\" dialog", () => clickOn(page, el("OK button in \"BitBIRCH Clustering\" dialog")));
      await session.step(57, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(58, "And a new column \"Cluster (BitBIRCH)\" should have been added", () => newColumnNamed(page, "Cluster (BitBIRCH)"));
      await session.step(59, "And \"Cluster (BitBIRCH)\" column should have at least 2 distinct values", () => distinctValues(page, "Cluster (BitBIRCH)", 2));
      await session.step(60, "And \"Cluster (BitBIRCH)\" column should have fewer distinct values than the table has rows", () => fewerDistinctThanRows(page, "Cluster (BitBIRCH)"));
      await session.step(61, "When user picks \"Chem > Calculate > Cluster MCS...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Cluster MCS..."));
      await session.step(62, "And user clicks on OK button in \"Cluster MCS\" dialog", () => clickOn(page, el("OK button in \"Cluster MCS\" dialog")));
      await session.step(63, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(64, "And a new column matching \"MCS|mcs|Scaffold\" should have been added", () => newColumnMatching(page, "MCS|mcs|Scaffold"));
      await session.step(65, "And the newest column matching \"MCS|mcs|Scaffold\" should have no missing values", () => newestMatchingFilled(page, "MCS|mcs|Scaffold"));
      await session.step(66, "And \"Cluster MCS\" column should have fewer distinct values than the table has rows", () => fewerDistinctThanRows(page, "Cluster MCS"));
      await session.step(67, "And the table should have 100 rows", () => rowCount(page, 100));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
