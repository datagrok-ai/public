/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/clustering/bio-dialog.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [dendrogram.cp.hier-clustering-bio-sequence-path]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBe, shouldHaveText, shouldHaveValue, shouldHaveValueBetween, shouldOffer} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, currentRowIs, hasColumn, mouseOverRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newestMatchingDistinct, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {taskBarShown, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noBalloons, noErrors, readingIs, readingNotAsRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {noSuchReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hierarchical clustering from the Bio menu", () => {
  const session = feature(test, "features/clustering/bio-dialog.feature", import.meta.url);
  test("Hierarchical clustering from the Bio menu", {tag: ["@journey", "@realizes:dendrogram.cp.hier-clustering-bio-sequence-path", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await run.scenario("The dialog opens on the sequence column with every distance and linkage", async () => {
      await session.step(18, "Then \"sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "sequence", "Macromolecule"));
      await session.step(19, "And the \"cell type of sequence\" reading of grid should be \"sequence\"", () => readingReads(page, "cell type of sequence", el("grid"), "sequence"));
      await session.step(20, "When user picks \"Bio > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Hierarchical Clustering..."));
      await session.step(21, "Then \"Hierarchical Clustering\" dialog should be visible", () => shouldBe(page, el("\"Hierarchical Clustering\" dialog"), "visible"));
      await session.step(22, "And Table input in \"Hierarchical Clustering\" dialog should have value \"FASTA_PT_activity\"", () => shouldHaveValue(page, el("Table input in \"Hierarchical Clustering\" dialog"), "FASTA_PT_activity"));
      await session.step(23, "And editor of Features input in \"Hierarchical Clustering\" dialog should have text \"(1) sequence\"", () => shouldHaveText(page, el("editor of Features input in \"Hierarchical Clustering\" dialog"), "(1) sequence"));
      await session.step(24, "And Distance input in \"Hierarchical Clustering\" dialog should have value \"euclidean\"", () => shouldHaveValue(page, el("Distance input in \"Hierarchical Clustering\" dialog"), "euclidean"));
      await session.step(25, "And Linkage input in \"Hierarchical Clustering\" dialog should have value \"ward\"", () => shouldHaveValue(page, el("Linkage input in \"Hierarchical Clustering\" dialog"), "ward"));
      await session.step(26, "And Distance input in \"Hierarchical Clustering\" dialog should offer \"euclidean, manhattan\"", () => shouldOffer(page, el("Distance input in \"Hierarchical Clustering\" dialog"), "euclidean, manhattan"));
      await session.step(27, "And Linkage input in \"Hierarchical Clustering\" dialog should offer \"single, complete, average, weighted, centroid, median, ward\"", () => shouldOffer(page, el("Linkage input in \"Hierarchical Clustering\" dialog"), "single, complete, average, weighted, centroid, median, ward"));
      await session.step(28, "When user clicks on CANCEL button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("CANCEL button in \"Hierarchical Clustering\" dialog")));
      await session.step(29, "Then \"Hierarchical Clustering\" dialog should be hidden", () => shouldBe(page, el("\"Hierarchical Clustering\" dialog"), "hidden"));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Euclidean distance and ward linkage attach a tree with a leaf for every sequence", async () => {
      await session.step(33, "Given user watches the task bar", () => watchTaskBar(page));
      await session.step(34, "When user picks \"Bio > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Hierarchical Clustering..."));
      await session.step(35, "And user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
      await session.step(36, "Then \"Hierarchical Clustering\" dialog should be hidden", () => shouldBe(page, el("\"Hierarchical Clustering\" dialog"), "hidden"));
      await session.step(37, "And the \"tree leaves\" reading of grid should be 99", () => readingIs(page, "tree leaves", el("grid"), 99));
      await session.step(38, "And the task bar should have shown \"Creating dendrogram\"", () => taskBarShown(page, "Creating dendrogram"));
      await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The sequence tree follows the grid, and the grid follows the tree", async () => {
      await session.step(43, "When user clicks on the \"cell 3 of sequence_id\" area of grid", () => clickArea(page, "cell 3 of sequence_id", el("grid")));
      await session.step(44, "Then row 3 should be current", () => currentRowIs(page, 3));
      await session.step(45, "And the \"tree current node\" reading of grid should be \"2\"", () => readingReads(page, "tree current node", el("grid"), "2"));
      await session.step(46, "When user hovers over the \"leaf 10\" area of grid", () => hoverArea(page, "leaf 10", el("grid")));
      await session.step(47, "Then the \"tree mouse over node\" reading of grid should be \"10\"", () => readingReads(page, "tree mouse over node", el("grid"), "10"));
      await session.step(48, "And row 11 should be under the mouse", () => mouseOverRowIs(page, 11));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the tree, then manhattan distance and complete linkage attach a tree of another height", async () => {
      await session.step(52, "When user remembers the \"tree height\" reading of grid", () => rememberReading(page, "tree height", el("grid")));
      await session.step(53, "And user clicks on \"Remove Dendrogram\" icon", () => clickOn(page, el("\"Remove Dendrogram\" icon")));
      await session.step(54, "Then grid should not report a \"tree leaves\" reading", () => noSuchReading(page, el("grid"), "tree leaves"));
      await session.step(55, "When user picks \"Bio > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Hierarchical Clustering..."));
      await session.step(56, "And user selects \"manhattan\" in Distance input in \"Hierarchical Clustering\" dialog", () => selectIn(page, "manhattan", el("Distance input in \"Hierarchical Clustering\" dialog")));
      await session.step(57, "And user selects \"complete\" in Linkage input in \"Hierarchical Clustering\" dialog", () => selectIn(page, "complete", el("Linkage input in \"Hierarchical Clustering\" dialog")));
      await session.step(58, "Then Distance input in \"Hierarchical Clustering\" dialog should have value \"manhattan\"", () => shouldHaveValue(page, el("Distance input in \"Hierarchical Clustering\" dialog"), "manhattan"));
      await session.step(59, "And Linkage input in \"Hierarchical Clustering\" dialog should have value \"complete\"", () => shouldHaveValue(page, el("Linkage input in \"Hierarchical Clustering\" dialog"), "complete"));
      await session.step(60, "When user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
      await session.step(61, "Then the \"tree leaves\" reading of grid should be 99", () => readingIs(page, "tree leaves", el("grid"), 99));
      await session.step(62, "And the \"tree height\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "tree height", el("grid")));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Assign Clusters on the sequence tree adds a cluster number to every sequence", async () => {
      await session.step(67, "When user clicks on \"Remove Dendrogram\" icon", () => clickOn(page, el("\"Remove Dendrogram\" icon")));
      await session.step(68, "And user picks \"Bio > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Hierarchical Clustering..."));
      await session.step(69, "And user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
      await session.step(70, "Then the \"tree leaves\" reading of grid should be 99", () => readingIs(page, "tree leaves", el("grid"), 99));
      await session.step(71, "When user clicks on \"Assign Clusters\" icon", () => clickOn(page, el("\"Assign Clusters\" icon")));
      await session.step(72, "Then \"Assign Clusters\" dialog should be visible", () => shouldBe(page, el("\"Assign Clusters\" dialog"), "visible"));
      await session.step(73, "When user enters \"5\" into Clusters input in \"Assign Clusters\" dialog", () => enterInto(page, "5", el("Clusters input in \"Assign Clusters\" dialog")));
      await session.step(74, "Then Threshold input in \"Assign Clusters\" dialog should have a value between 11.05 and 11.07", () => shouldHaveValueBetween(page, el("Threshold input in \"Assign Clusters\" dialog"), 11.05, 11.07));
      await session.step(75, "When user clicks on Assign button in \"Assign Clusters\" dialog", () => clickOn(page, el("Assign button in \"Assign Clusters\" dialog")));
      await session.step(76, "Then the \"Assign Clusters\" dialog should close", () => dialogCloses(page, "Assign Clusters"));
      await session.step(77, "And the table should have a column \"Cluster (11.06)\"", () => hasColumn(page, "Cluster (11.06)"));
      await session.step(78, "And \"Cluster (11.06)\" column should have no missing values", () => columnComplete(page, "Cluster (11.06)"));
      await session.step(79, "And the newest column matching \"^Cluster \\(\" should have 5 distinct values", () => newestMatchingDistinct(page, "^Cluster \\(", 5));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A threshold above the tree's height keeps at least one cluster", async () => {
      await session.step(84, "When user clicks on \"Assign Clusters\" icon", () => clickOn(page, el("\"Assign Clusters\" icon")));
      await session.step(85, "And user enters \"100\" into Threshold input in \"Assign Clusters\" dialog", () => enterInto(page, "100", el("Threshold input in \"Assign Clusters\" dialog")));
      await session.step(86, "Then Clusters input in \"Assign Clusters\" dialog should have a value between 1 and 99", () => shouldHaveValueBetween(page, el("Clusters input in \"Assign Clusters\" dialog"), 1, 99));
    }, {knownFailure: true});
    run.finish();
  });
});
