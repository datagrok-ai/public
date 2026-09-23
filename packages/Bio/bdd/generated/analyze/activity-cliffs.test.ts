/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/activity-cliffs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.analyze.activity-cliffs]
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
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, painted, propertyShouldBe, propertyShouldContain, readingAtLeast} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sequence Activity Cliffs", () => {
  const session = feature(test, "features/analyze/activity-cliffs.feature", import.meta.url);
  test("Sequence Activity Cliffs", {tag: ["@journey", "@realizes:bio.analyze.activity-cliffs"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens FASTA_sample dataset", () => openDataset(page, ds("FASTA_sample")));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(17, "Then \"Sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "Sequence", "Macromolecule"));
    await run.scenario("The editor opens with the sequence column, the first numeric column as the activity, and a cutoff", async () => {
      await session.step(20, "When user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
      await session.step(21, "Then \"Sequence Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "visible"));
      await session.step(22, "And editor of Column input in \"Sequence Activity Cliffs\" dialog should have text \"Sequence\"", () => shouldHaveText(page, el("editor of Column input in \"Sequence Activity Cliffs\" dialog"), "Sequence"));
      await session.step(23, "And Method input in \"Sequence Activity Cliffs\" dialog should have value \"UMAP\"", () => shouldHaveValue(page, el("Method input in \"Sequence Activity Cliffs\" dialog"), "UMAP"));
      await session.step(24, "And Similarity input in \"Sequence Activity Cliffs\" dialog should have value \"Hamming\"", () => shouldHaveValue(page, el("Similarity input in \"Sequence Activity Cliffs\" dialog"), "Hamming"));
      await session.step(25, "And \"Similarity cutoff\" input in \"Sequence Activity Cliffs\" dialog should have value \"80\"", () => shouldHaveValue(page, el("\"Similarity cutoff\" input in \"Sequence Activity Cliffs\" dialog"), "80"));
      await session.step(26, "And editor of Activities input in \"Sequence Activity Cliffs\" dialog should have text \"Length\"", () => shouldHaveText(page, el("editor of Activities input in \"Sequence Activity Cliffs\" dialog"), "Length"));
    });
    await run.scenario("Running on the Activity column docks a cliff scatter plot", async () => {
      await session.step(29, "When user selects \"Activity\" in Activities input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Activity", el("Activities input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(30, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
      await session.step(31, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(32, "And \"Sequence Activity Cliffs\" dialog should be hidden", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "hidden"));
      await session.step(33, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
      await session.step(34, "And a new column \"Embed_Y_1\" should have been added", () => newColumnNamed(page, "Embed_Y_1"));
      await session.step(35, "And a new column matching \"sali|SALI\" should have been added", () => newColumnMatching(page, "sali|SALI"));
      await session.step(36, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(37, "And title of scatter plot viewer should have text \"Activity cliffs\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Activity cliffs"));
      await session.step(38, "And \"X\" property of scatter plot viewer should be \"Embed_X_1\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "Embed_X_1"));
      await session.step(39, "And \"Y\" property of scatter plot viewer should be \"Embed_Y_1\"", () => propertyShouldBe(page, "Y", el("scatter plot viewer"), "Embed_Y_1"));
      await session.step(40, "And the \"cliffs\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "cliffs", el("scatter plot viewer"), 1));
      await session.step(41, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second analysis with other settings docks alongside the first", async () => {
      await session.step(46, "When user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
      await session.step(47, "Then \"Sequence Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "visible"));
      await session.step(48, "When user selects \"Activity\" in Activities input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Activity", el("Activities input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(49, "And user selects \"t-SNE\" in Method input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "t-SNE", el("Method input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(50, "And user selects \"Levenshtein\" in Similarity input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Levenshtein", el("Similarity input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(51, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
      await session.step(52, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(53, "And a new column \"Embed_X_2\" should have been added", () => newColumnNamed(page, "Embed_X_2"));
      await session.step(54, "And second scatter plot viewer should be visible", () => shouldBe(page, el("second scatter plot viewer"), "visible"));
      await session.step(55, "And \"X\" property of second scatter plot viewer should be \"Embed_X_2\"", () => propertyShouldBe(page, "X", el("second scatter plot viewer"), "Embed_X_2"));
      await session.step(56, "And \"Description\" property of second scatter plot viewer should contain \"method: t-SNE\"", () => propertyShouldContain(page, "Description", el("second scatter plot viewer"), "method: t-SNE"));
      await session.step(57, "And \"Description\" property of second scatter plot viewer should contain \"similarity: Levenshtein\"", () => propertyShouldContain(page, "Description", el("second scatter plot viewer"), "similarity: Levenshtein"));
      await session.step(58, "And the \"cliffs\" reading of second scatter plot viewer should be at least 1", () => readingAtLeast(page, "cliffs", el("second scatter plot viewer"), 1));
      await session.step(59, "And second scatter plot viewer should be painted", () => painted(page, el("second scatter plot viewer")));
      await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
