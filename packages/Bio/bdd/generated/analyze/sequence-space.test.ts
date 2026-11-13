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
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, selectIn, shouldBe, shouldHaveText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnUnits, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {closeAllViews, dialogCloses, noProjectOnServer, openDataset, openProjectWithTable, openSaveDialog, projectsOnServer, saveFromDialog, savedAsSnapshot} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, painted, propertyShouldBe, propertyShouldContain} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sequence Space", () => {
  const session = feature(test, "features/analyze/sequence-space.feature", import.meta.url);
  test("Sequence Space", {tag: ["@journey", "@realizes:bio.analyze.sequence-space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
    await session.step(17, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(18, "Then \"fasta\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "fasta", "Macromolecule"));
    await session.step(19, "And \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
    await run.scenario("The editor opens on the sequence column with the default engine", async () => {
      await session.step(22, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
      await session.step(23, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
      await session.step(24, "And editor of Column input in \"Sequence Space\" dialog should have text \"fasta\"", () => shouldHaveText(page, el("editor of Column input in \"Sequence Space\" dialog"), "fasta"));
      await session.step(25, "And Method input in \"Sequence Space\" dialog should have value \"UMAP\"", () => shouldHaveValue(page, el("Method input in \"Sequence Space\" dialog"), "UMAP"));
      await session.step(26, "And Similarity input in \"Sequence Space\" dialog should have value \"Hamming\"", () => shouldHaveValue(page, el("Similarity input in \"Sequence Space\" dialog"), "Hamming"));
      await session.step(27, "And \"Plot embeddings\" checkbox in \"Sequence Space\" dialog should be checked", () => shouldBe(page, el("\"Plot embeddings\" checkbox in \"Sequence Space\" dialog"), "checked"));
      await session.step(28, "And \"Cluster embeddings\" checkbox in \"Sequence Space\" dialog should be checked", () => shouldBe(page, el("\"Cluster embeddings\" checkbox in \"Sequence Space\" dialog"), "checked"));
    });
    await run.scenario("Running with the defaults appends the embeddings and docks the scatter plot", async () => {
      await session.step(31, "When user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
      await session.step(32, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(33, "And \"Sequence Space\" dialog should be hidden", () => shouldBe(page, el("\"Sequence Space\" dialog"), "hidden"));
      await session.step(34, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
      await session.step(35, "And a new column \"Embed_Y_1\" should have been added", () => newColumnNamed(page, "Embed_Y_1"));
      await session.step(36, "And a new column matching \"^Cluster \\(DBSCAN\\)\" should have been added", () => newColumnMatching(page, "^Cluster \\(DBSCAN\\)"));
      await session.step(37, "And \"Embed_X_1\" column should have no missing values", () => columnComplete(page, "Embed_X_1"));
      await session.step(38, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(39, "And title of scatter plot viewer should have text \"Sequence space\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Sequence space"));
      await session.step(40, "And \"X\" property of scatter plot viewer should be \"Embed_X_1\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "Embed_X_1"));
      await session.step(41, "And \"Y\" property of scatter plot viewer should be \"Embed_Y_1\"", () => propertyShouldBe(page, "Y", el("scatter plot viewer"), "Embed_Y_1"));
      await session.step(42, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The editor reopens over the first result and takes another method and metric", async () => {
      await session.step(47, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
      await session.step(48, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
      await session.step(49, "When user selects \"t-SNE\" in Method input in \"Sequence Space\" dialog", () => selectIn(page, "t-SNE", el("Method input in \"Sequence Space\" dialog")));
      await session.step(50, "And user selects \"Levenshtein\" in Similarity input in \"Sequence Space\" dialog", () => selectIn(page, "Levenshtein", el("Similarity input in \"Sequence Space\" dialog")));
      await session.step(51, "Then Method input in \"Sequence Space\" dialog should have value \"t-SNE\"", () => shouldHaveValue(page, el("Method input in \"Sequence Space\" dialog"), "t-SNE"));
      await session.step(52, "And Similarity input in \"Sequence Space\" dialog should have value \"Levenshtein\"", () => shouldHaveValue(page, el("Similarity input in \"Sequence Space\" dialog"), "Levenshtein"));
    });
    await run.scenario("The edited run is a second result computed with the edited settings", async () => {
      await session.step(55, "When user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
      await session.step(56, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(57, "And a new column \"Embed_X_2\" should have been added", () => newColumnNamed(page, "Embed_X_2"));
      await session.step(58, "And a new column \"Embed_Y_2\" should have been added", () => newColumnNamed(page, "Embed_Y_2"));
      await session.step(59, "And a new column matching \"^Cluster \\(DBSCAN\\)\" should have been added", () => newColumnMatching(page, "^Cluster \\(DBSCAN\\)"));
      await session.step(60, "And second scatter plot viewer should be visible", () => shouldBe(page, el("second scatter plot viewer"), "visible"));
      await session.step(61, "And \"X\" property of second scatter plot viewer should be \"Embed_X_2\"", () => propertyShouldBe(page, "X", el("second scatter plot viewer"), "Embed_X_2"));
      await session.step(62, "And \"Description\" property of second scatter plot viewer should contain \"method: t-SNE\"", () => propertyShouldContain(page, "Description", el("second scatter plot viewer"), "method: t-SNE"));
      await session.step(63, "And \"Description\" property of second scatter plot viewer should contain \"similarity: Levenshtein\"", () => propertyShouldContain(page, "Description", el("second scatter plot viewer"), "similarity: Levenshtein"));
      await session.step(64, "And second scatter plot viewer should be painted", () => painted(page, el("second scatter plot viewer")));
      await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Both results come back from a saved project (GROK-19928)", async () => {
      await session.step(69, "Given no project named \"bdd-bio-seqspace-{run}\" is on the server", () => noProjectOnServer(page, session.text("bdd-bio-seqspace-{run}")));
      await session.step(70, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(71, "And user types \"bdd-bio-seqspace-{run}\" into text input in \"Save project\" dialog", () => typeInto(page, session.text("bdd-bio-seqspace-{run}"), el("text input in \"Save project\" dialog")));
      await session.step(72, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(73, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(74, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(75, "Then 1 project named \"bdd-bio-seqspace-{run}\" should be on the server", () => projectsOnServer(page, 1, session.text("bdd-bio-seqspace-{run}")));
      await session.step(76, "And the \"filter_FASTA\" table of the \"bdd-bio-seqspace-{run}\" project should be saved as a snapshot", () => savedAsSnapshot(page, "filter_FASTA", session.text("bdd-bio-seqspace-{run}")));
      await session.step(77, "When user closes all views", () => closeAllViews(page));
      await session.step(78, "And user opens the \"bdd-bio-seqspace-{run}\" project and waits for its table", () => openProjectWithTable(page, session.text("bdd-bio-seqspace-{run}")));
      await session.step(79, "Then \"fasta\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "fasta", "Macromolecule"));
      await session.step(80, "And \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
      await session.step(81, "And the table should have a column \"Embed_X_1\"", () => hasColumn(page, "Embed_X_1"));
      await session.step(82, "And the table should have a column \"Embed_Y_2\"", () => hasColumn(page, "Embed_Y_2"));
      await session.step(83, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(84, "And \"X\" property of scatter plot viewer should be \"Embed_X_1\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "Embed_X_1"));
      await session.step(85, "And \"Y\" property of scatter plot viewer should be \"Embed_Y_1\"", () => propertyShouldBe(page, "Y", el("scatter plot viewer"), "Embed_Y_1"));
      await session.step(86, "And \"X\" property of second scatter plot viewer should be \"Embed_X_2\"", () => propertyShouldBe(page, "X", el("second scatter plot viewer"), "Embed_X_2"));
      await session.step(87, "And \"Description\" property of second scatter plot viewer should contain \"method: t-SNE\"", () => propertyShouldContain(page, "Description", el("second scatter plot viewer"), "method: t-SNE"));
      await session.step(88, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(89, "And second scatter plot viewer should be painted", () => painted(page, el("second scatter plot viewer")));
      await session.step(90, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
