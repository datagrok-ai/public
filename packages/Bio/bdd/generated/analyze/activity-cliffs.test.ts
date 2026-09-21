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
import {noBalloons, noErrors, painted, propertyShouldBe, propertyShouldContain} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sequence Activity Cliffs", () => {
  const session = feature(test, "features/analyze/activity-cliffs.feature", import.meta.url);
  test("Sequence Activity Cliffs", {tag: ["@journey", "@realizes:bio.analyze.activity-cliffs"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And user opens FASTA_sample dataset", () => openDataset(page, ds("FASTA_sample")));
    await session.step(10, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(11, "Then \"Sequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "Sequence", "Macromolecule"));
    await run.scenario("The editor opens with the sequence column, an activity and a cutoff", async () => {
      await session.step(14, "When user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
      await session.step(15, "Then \"Sequence Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "visible"));
      await session.step(16, "And editor of Column input in \"Sequence Activity Cliffs\" dialog should have text \"Sequence\"", () => shouldHaveText(page, el("editor of Column input in \"Sequence Activity Cliffs\" dialog"), "Sequence"));
      await session.step(17, "And Method input in \"Sequence Activity Cliffs\" dialog should have value \"UMAP\"", () => shouldHaveValue(page, el("Method input in \"Sequence Activity Cliffs\" dialog"), "UMAP"));
      await session.step(18, "And Similarity input in \"Sequence Activity Cliffs\" dialog should have value \"Hamming\"", () => shouldHaveValue(page, el("Similarity input in \"Sequence Activity Cliffs\" dialog"), "Hamming"));
      await session.step(19, "And \"Similarity cutoff\" input in \"Sequence Activity Cliffs\" dialog should have value \"80\"", () => shouldHaveValue(page, el("\"Similarity cutoff\" input in \"Sequence Activity Cliffs\" dialog"), "80"));
    });
    await run.scenario("Running on the Activity column docks a cliff scatter plot", async () => {
      await session.step(22, "When user selects \"Activity\" in Activities input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Activity", el("Activities input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(23, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
      await session.step(24, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(25, "And \"Sequence Activity Cliffs\" dialog should be hidden", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "hidden"));
      await session.step(26, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
      await session.step(27, "And a new column \"Embed_Y_1\" should have been added", () => newColumnNamed(page, "Embed_Y_1"));
      await session.step(28, "And a new column matching \"sali|SALI\" should have been added", () => newColumnMatching(page, "sali|SALI"));
      await session.step(29, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(30, "And title of scatter plot viewer should have text \"Activity cliffs\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Activity cliffs"));
      await session.step(31, "And \"X\" property of scatter plot viewer should be \"Embed_X_1\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "Embed_X_1"));
      await session.step(32, "And \"Y\" property of scatter plot viewer should be \"Embed_Y_1\"", () => propertyShouldBe(page, "Y", el("scatter plot viewer"), "Embed_Y_1"));
      await session.step(33, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(34, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second analysis with other settings docks alongside the first", async () => {
      await session.step(38, "When user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
      await session.step(39, "Then \"Sequence Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "visible"));
      await session.step(40, "When user selects \"Activity\" in Activities input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Activity", el("Activities input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(41, "And user selects \"t-SNE\" in Method input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "t-SNE", el("Method input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(42, "And user selects \"Levenshtein\" in Similarity input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Levenshtein", el("Similarity input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(43, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
      await session.step(44, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(45, "And a new column \"Embed_X_2\" should have been added", () => newColumnNamed(page, "Embed_X_2"));
      await session.step(46, "And second scatter plot viewer should be visible", () => shouldBe(page, el("second scatter plot viewer"), "visible"));
      await session.step(47, "And \"X\" property of second scatter plot viewer should be \"Embed_X_2\"", () => propertyShouldBe(page, "X", el("second scatter plot viewer"), "Embed_X_2"));
      await session.step(48, "And \"Description\" property of second scatter plot viewer should contain \"method: t-SNE\"", () => propertyShouldContain(page, "Description", el("second scatter plot viewer"), "method: t-SNE"));
      await session.step(49, "And \"Description\" property of second scatter plot viewer should contain \"similarity: Levenshtein\"", () => propertyShouldContain(page, "Description", el("second scatter plot viewer"), "similarity: Levenshtein"));
      await session.step(50, "And second scatter plot viewer should be painted", () => painted(page, el("second scatter plot viewer")));
      await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
