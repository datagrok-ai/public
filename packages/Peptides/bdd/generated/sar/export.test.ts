/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sar/export.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {peptidesInitialized, sarReady} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {tableColumns} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {closeCurrentView, openDataset, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, pickFromContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Export peptide SAR results", () => {
  const session = feature(test, "features/sar/export.feature", import.meta.url);
  test("Export peptide SAR results", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(15, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(16, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(17, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(18, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(19, "And user enters \"93\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "93", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(20, "And Scaling input in \"Analyze Peptides\" dialog should have value \"none\"", () => shouldHaveValue(page, el("Scaling input in \"Analyze Peptides\" dialog"), "none"));
    await session.step(21, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(22, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(23, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(24, "Then no errors should have been logged", () => noErrors(page));
    await session.step(25, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("Both SAR viewers export the monomer-by-position count matrix [viewer=Sequence Variability Map]", async () => {
      await session.step(28, "When user picks \"Export > Export Invariant Map\" from the context menu of Sequence Variability Map viewer", () => pickFromContextMenu(page, "Export > Export Invariant Map", el("Sequence Variability Map viewer")));
      await session.step(29, "Then the \"Invariant Map\" view should be current", () => viewIsCurrent(page, "Invariant Map"));
      await session.step(30, "And table \"Invariant Map\" should have columns \"AAR, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17\"", () => tableColumns(page, "Invariant Map", "AAR, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17"));
      await session.step(31, "When user closes the current view", () => closeCurrentView(page));
      await session.step(32, "And user switches to the \"peptides\" table view", () => switchTableView(page, "peptides"));
      await session.step(33, "Then no errors should have been logged", () => noErrors(page));
      await session.step(34, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Both SAR viewers export the monomer-by-position count matrix [viewer=Most Potent Residues]", async () => {
      await session.step(28, "When user picks \"Export > Export Invariant Map\" from the context menu of Most Potent Residues viewer", () => pickFromContextMenu(page, "Export > Export Invariant Map", el("Most Potent Residues viewer")));
      await session.step(29, "Then the \"Invariant Map\" view should be current", () => viewIsCurrent(page, "Invariant Map"));
      await session.step(30, "And table \"Invariant Map\" should have columns \"AAR, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17\"", () => tableColumns(page, "Invariant Map", "AAR, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17"));
      await session.step(31, "When user closes the current view", () => closeCurrentView(page));
      await session.step(32, "And user switches to the \"peptides\" table view", () => switchTableView(page, "peptides"));
      await session.step(33, "Then no errors should have been logged", () => noErrors(page));
      await session.step(34, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The mutation-cliff export dialog makes a table of the pairs and their activities", async () => {
      await session.step(42, "When user picks \"Export > Export Mutation Cliffs...\" from the context menu of Sequence Variability Map viewer", () => pickFromContextMenu(page, "Export > Export Mutation Cliffs...", el("Sequence Variability Map viewer")));
      await session.step(43, "Then \"Export Mutation Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Export Mutation Cliffs\" dialog"), "visible"));
      await session.step(44, "And \"Extra columns\" input in \"Export Mutation Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Extra columns\" input in \"Export Mutation Cliffs\" dialog"), "visible"));
      await session.step(45, "When user clicks on OK button in \"Export Mutation Cliffs\" dialog", () => clickOn(page, el("OK button in \"Export Mutation Cliffs\" dialog")));
      await session.step(46, "Then \"Export Mutation Cliffs\" dialog should be hidden", () => shouldBe(page, el("\"Export Mutation Cliffs\" dialog"), "hidden"));
      await session.step(47, "And the \"Mutation Cliffs\" view should be current", () => viewIsCurrent(page, "Mutation Cliffs"));
      await session.step(48, "And table \"Mutation Cliffs\" should have columns \"Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta\"", () => tableColumns(page, "Mutation Cliffs", "Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta"));
      await session.step(49, "When user closes the current view", () => closeCurrentView(page));
      await session.step(50, "And user switches to the \"peptides\" table view", () => switchTableView(page, "peptides"));
      await session.step(51, "Then no errors should have been logged", () => noErrors(page));
      await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Selecting ID in the export dialog adds identifier columns for both members of a pair", async () => {
      await session.step(55, "When user picks \"Export > Export Mutation Cliffs...\" from the context menu of Sequence Variability Map viewer", () => pickFromContextMenu(page, "Export > Export Mutation Cliffs...", el("Sequence Variability Map viewer")));
      await session.step(56, "And user clicks on \"Extra columns\" input in \"Export Mutation Cliffs\" dialog", () => clickOn(page, el("\"Extra columns\" input in \"Export Mutation Cliffs\" dialog")));
      await session.step(57, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(58, "And the \"text of cell 2 of __name\" reading of grid in \"Select columns...\" dialog should be \"ID\"", () => readingReads(page, "text of cell 2 of __name", el("grid in \"Select columns...\" dialog"), "ID"));
      await session.step(59, "When user clicks on the \"cell 2 of x\" area of grid in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid in \"Select columns...\" dialog")));
      await session.step(60, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(61, "And user clicks on OK button in \"Export Mutation Cliffs\" dialog", () => clickOn(page, el("OK button in \"Export Mutation Cliffs\" dialog")));
      await session.step(62, "Then the \"Mutation Cliffs\" view should be current", () => viewIsCurrent(page, "Mutation Cliffs"));
      await session.step(63, "And table \"Mutation Cliffs\" should have columns \"Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta, Seq 1 ID, Seq 2 ID\"", () => tableColumns(page, "Mutation Cliffs", "Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta, Seq 1 ID, Seq 2 ID"));
      await session.step(64, "When user closes the current view", () => closeCurrentView(page));
      await session.step(65, "And user switches to the \"peptides\" table view", () => switchTableView(page, "peptides"));
      await session.step(66, "Then no errors should have been logged", () => noErrors(page));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
