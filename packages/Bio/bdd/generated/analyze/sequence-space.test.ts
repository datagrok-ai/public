/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/sequence-space.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.analyze.sequence-space]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, painted, propertyShouldBe, propertyShouldContain} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sequence Space", () => {
  const session = feature(test, "features/analyze/sequence-space.feature", import.meta.url);
  test("Sequence Space", {tag: ["@journey", "@realizes:bio.analyze.sequence-space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
    await session.step(10, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(11, "Then \"fasta\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "fasta", "Macromolecule"));
    await session.step(12, "And \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
    await run.scenario("The editor opens on the sequence column with the default engine", async () => {
      await session.step(15, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
      await session.step(16, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
      await session.step(17, "And editor of Column input in \"Sequence Space\" dialog should have text \"fasta\"", () => shouldHaveText(page, el("editor of Column input in \"Sequence Space\" dialog"), "fasta"));
      await session.step(18, "And Method input in \"Sequence Space\" dialog should have value \"UMAP\"", () => shouldHaveValue(page, el("Method input in \"Sequence Space\" dialog"), "UMAP"));
      await session.step(19, "And Similarity input in \"Sequence Space\" dialog should have value \"Hamming\"", () => shouldHaveValue(page, el("Similarity input in \"Sequence Space\" dialog"), "Hamming"));
      await session.step(20, "And \"Plot embeddings\" checkbox in \"Sequence Space\" dialog should be checked", () => shouldBe(page, el("\"Plot embeddings\" checkbox in \"Sequence Space\" dialog"), "checked"));
      await session.step(21, "And \"Cluster embeddings\" checkbox in \"Sequence Space\" dialog should be checked", () => shouldBe(page, el("\"Cluster embeddings\" checkbox in \"Sequence Space\" dialog"), "checked"));
    });
    await run.scenario("Running with the defaults appends the embeddings and docks the scatter plot", async () => {
      await session.step(24, "When user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
      await session.step(25, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(26, "And \"Sequence Space\" dialog should be hidden", () => shouldBe(page, el("\"Sequence Space\" dialog"), "hidden"));
      await session.step(27, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
      await session.step(28, "And a new column \"Embed_Y_1\" should have been added", () => newColumnNamed(page, "Embed_Y_1"));
      await session.step(29, "And a new column matching \"^Cluster \\(DBSCAN\\)\" should have been added", () => newColumnMatching(page, "^Cluster \\(DBSCAN\\)"));
      await session.step(30, "And \"Embed_X_1\" column should have no missing values", () => columnComplete(page, "Embed_X_1"));
      await session.step(31, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(32, "And title of scatter plot viewer should have text \"Sequence space\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Sequence space"));
      await session.step(33, "And \"X\" property of scatter plot viewer should be \"Embed_X_1\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "Embed_X_1"));
      await session.step(34, "And \"Y\" property of scatter plot viewer should be \"Embed_Y_1\"", () => propertyShouldBe(page, "Y", el("scatter plot viewer"), "Embed_Y_1"));
      await session.step(35, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The editor reopens over the first result and takes another method and metric", async () => {
      await session.step(40, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
      await session.step(41, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
      await session.step(42, "When user selects \"t-SNE\" in Method input in \"Sequence Space\" dialog", () => selectIn(page, "t-SNE", el("Method input in \"Sequence Space\" dialog")));
      await session.step(43, "And user selects \"Levenshtein\" in Similarity input in \"Sequence Space\" dialog", () => selectIn(page, "Levenshtein", el("Similarity input in \"Sequence Space\" dialog")));
      await session.step(44, "Then Method input in \"Sequence Space\" dialog should have value \"t-SNE\"", () => shouldHaveValue(page, el("Method input in \"Sequence Space\" dialog"), "t-SNE"));
      await session.step(45, "And Similarity input in \"Sequence Space\" dialog should have value \"Levenshtein\"", () => shouldHaveValue(page, el("Similarity input in \"Sequence Space\" dialog"), "Levenshtein"));
    });
    await run.scenario("The edited run is a second result computed with the edited settings", async () => {
      await session.step(48, "When user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
      await session.step(49, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(50, "And a new column \"Embed_X_2\" should have been added", () => newColumnNamed(page, "Embed_X_2"));
      await session.step(51, "And a new column \"Embed_Y_2\" should have been added", () => newColumnNamed(page, "Embed_Y_2"));
      await session.step(52, "And second scatter plot viewer should be visible", () => shouldBe(page, el("second scatter plot viewer"), "visible"));
      await session.step(53, "And \"X\" property of second scatter plot viewer should be \"Embed_X_2\"", () => propertyShouldBe(page, "X", el("second scatter plot viewer"), "Embed_X_2"));
      await session.step(54, "And \"Description\" property of second scatter plot viewer should contain \"method: t-SNE\"", () => propertyShouldContain(page, "Description", el("second scatter plot viewer"), "method: t-SNE"));
      await session.step(55, "And \"Description\" property of second scatter plot viewer should contain \"similarity: Levenshtein\"", () => propertyShouldContain(page, "Description", el("second scatter plot viewer"), "similarity: Levenshtein"));
      await session.step(56, "And second scatter plot viewer should be painted", () => painted(page, el("second scatter plot viewer")));
      await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
