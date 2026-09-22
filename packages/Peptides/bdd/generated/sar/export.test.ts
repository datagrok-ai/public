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
import {invariantCounts, mutationExtraValues, mutationPairs} from '../../bindings/exports.js';
import {peptidesInitialized, sarReady} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnType, columnUnits, joinedValues, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, tableColumns} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {closeCurrentView, openDataset, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, pickFromContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Export peptide SAR results", () => {
  const session = feature(test, "features/sar/export.feature", import.meta.url);
  test("Export peptide SAR results", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(13, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(14, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(15, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(16, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(17, "And user enters \"93\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "93", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(18, "And Scaling input in \"Analyze Peptides\" dialog should have value \"none\"", () => shouldHaveValue(page, el("Scaling input in \"Analyze Peptides\" dialog"), "none"));
    await session.step(19, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(20, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(21, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(22, "Then no errors should have been logged", () => noErrors(page));
    await session.step(23, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("Both SAR viewers export the full monomer-by-position count matrix [viewer=Sequence Variability Map]", async () => {
      await session.step(26, "When user picks \"Export > Export Invariant Map\" from the context menu of Sequence Variability Map viewer", () => pickFromContextMenu(page, "Export > Export Invariant Map", el("Sequence Variability Map viewer")));
      await session.step(27, "Then the \"Invariant Map\" view should be current", () => viewIsCurrent(page, "Invariant Map"));
      await session.step(28, "And table \"Invariant Map\" should have columns \"AAR, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17\"", () => tableColumns(page, "Invariant Map", "AAR, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17"));
      await session.step(29, "And the table should have 22 rows", () => rowCount(page, 22));
      await session.step(30, "And \"AAR\" column should have type \"string\"", () => columnType(page, "AAR", "string"));
      await session.step(31, "And the value of \"AAR\" column in row 14 should be \"NH2\"", () => valueInRow(page, "AAR", 14, "NH2"));
      await session.step(32, "And the value of \"1\" column in row 14 should be \"647\"", () => valueInRow(page, "1", 14, "647"));
      await session.step(33, "And the value of \"17\" column in row 3 should be \"647\"", () => valueInRow(page, "17", 3, "647"));
      await session.step(34, "And the value of \"2\" column in row 1 should be \"299\"", () => valueInRow(page, "2", 1, "299"));
      await session.step(35, "And the value of \"2\" column in row 12 should be \"9\"", () => valueInRow(page, "2", 12, "9"));
      await session.step(36, "And the value of \"10\" column in row 22 should be \"604\"", () => valueInRow(page, "10", 22, "604"));
      await session.step(37, "And the value of \"1\" column in row 1 should be \"0\"", () => valueInRow(page, "1", 1, "0"));
      await session.step(38, "And the invariant-map export should match the monomer counts of table \"peptides\"", () => invariantCounts(page, "peptides"));
      await session.step(39, "When user closes the current view", () => closeCurrentView(page));
      await session.step(40, "And user switches to the \"peptides\" table view", () => switchTableView(page, "peptides"));
      await session.step(41, "Then no errors should have been logged", () => noErrors(page));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Both SAR viewers export the full monomer-by-position count matrix [viewer=Most Potent Residues]", async () => {
      await session.step(26, "When user picks \"Export > Export Invariant Map\" from the context menu of Most Potent Residues viewer", () => pickFromContextMenu(page, "Export > Export Invariant Map", el("Most Potent Residues viewer")));
      await session.step(27, "Then the \"Invariant Map\" view should be current", () => viewIsCurrent(page, "Invariant Map"));
      await session.step(28, "And table \"Invariant Map\" should have columns \"AAR, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17\"", () => tableColumns(page, "Invariant Map", "AAR, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17"));
      await session.step(29, "And the table should have 22 rows", () => rowCount(page, 22));
      await session.step(30, "And \"AAR\" column should have type \"string\"", () => columnType(page, "AAR", "string"));
      await session.step(31, "And the value of \"AAR\" column in row 14 should be \"NH2\"", () => valueInRow(page, "AAR", 14, "NH2"));
      await session.step(32, "And the value of \"1\" column in row 14 should be \"647\"", () => valueInRow(page, "1", 14, "647"));
      await session.step(33, "And the value of \"17\" column in row 3 should be \"647\"", () => valueInRow(page, "17", 3, "647"));
      await session.step(34, "And the value of \"2\" column in row 1 should be \"299\"", () => valueInRow(page, "2", 1, "299"));
      await session.step(35, "And the value of \"2\" column in row 12 should be \"9\"", () => valueInRow(page, "2", 12, "9"));
      await session.step(36, "And the value of \"10\" column in row 22 should be \"604\"", () => valueInRow(page, "10", 22, "604"));
      await session.step(37, "And the value of \"1\" column in row 1 should be \"0\"", () => valueInRow(page, "1", 1, "0"));
      await session.step(38, "And the invariant-map export should match the monomer counts of table \"peptides\"", () => invariantCounts(page, "peptides"));
      await session.step(39, "When user closes the current view", () => closeCurrentView(page));
      await session.step(40, "And user switches to the \"peptides\" table view", () => switchTableView(page, "peptides"));
      await session.step(41, "Then no errors should have been logged", () => noErrors(page));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Mutation-cliff export preserves every pair and its activities", async () => {
      await session.step(50, "When user picks \"Export > Export Mutation Cliffs...\" from the context menu of Sequence Variability Map viewer", () => pickFromContextMenu(page, "Export > Export Mutation Cliffs...", el("Sequence Variability Map viewer")));
      await session.step(51, "Then \"Export Mutation Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Export Mutation Cliffs\" dialog"), "visible"));
      await session.step(52, "And \"Extra columns\" input in \"Export Mutation Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Extra columns\" input in \"Export Mutation Cliffs\" dialog"), "visible"));
      await session.step(53, "When user clicks on OK button in \"Export Mutation Cliffs\" dialog", () => clickOn(page, el("OK button in \"Export Mutation Cliffs\" dialog")));
      await session.step(54, "Then \"Export Mutation Cliffs\" dialog should be hidden", () => shouldBe(page, el("\"Export Mutation Cliffs\" dialog"), "hidden"));
      await session.step(55, "And the \"Mutation Cliffs\" view should be current", () => viewIsCurrent(page, "Mutation Cliffs"));
      await session.step(56, "And table \"Mutation Cliffs\" should have columns \"Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta\"", () => tableColumns(page, "Mutation Cliffs", "Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta"));
      await session.step(57, "And the table should have 6253 rows", () => rowCount(page, 6253));
      await session.step(58, "And \"Seq 1\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "Seq 1", "Macromolecule"));
      await session.step(59, "And \"Seq 2\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "Seq 2", "Macromolecule"));
      await session.step(60, "And \"Seq 1\" column should have units \"separator\"", () => columnUnits(page, "Seq 1", "separator"));
      await session.step(61, "And \"Mutation\" column should have semantic type \"MacromoleculeDifference\"", () => columnSemType(page, "Mutation", "MacromoleculeDifference"));
      await session.step(62, "And every value of \"Mutation\" column should be \"Seq 1\" and \"Seq 2\" of the same row joined by \"#\"", () => joinedValues(page, "Mutation", "Seq 1", "Seq 2", "#"));
      await session.step(63, "And the mutation-cliff export should contain every single-mutation pair from table \"peptides\"", () => mutationPairs(page, "peptides"));
      await session.step(64, "When user closes the current view", () => closeCurrentView(page));
      await session.step(65, "And user switches to the \"peptides\" table view", () => switchTableView(page, "peptides"));
      await session.step(66, "Then no errors should have been logged", () => noErrors(page));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Selecting ID exports identifiers for both members of every pair", async () => {
      await session.step(70, "When user picks \"Export > Export Mutation Cliffs...\" from the context menu of Sequence Variability Map viewer", () => pickFromContextMenu(page, "Export > Export Mutation Cliffs...", el("Sequence Variability Map viewer")));
      await session.step(71, "And user clicks on \"Extra columns\" input in \"Export Mutation Cliffs\" dialog", () => clickOn(page, el("\"Extra columns\" input in \"Export Mutation Cliffs\" dialog")));
      await session.step(72, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(73, "And the \"text of cell 2 of __name\" reading of grid in \"Select columns...\" dialog should be \"ID\"", () => readingReads(page, "text of cell 2 of __name", el("grid in \"Select columns...\" dialog"), "ID"));
      await session.step(74, "When user clicks on the \"cell 2 of x\" area of grid in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid in \"Select columns...\" dialog")));
      await session.step(75, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(76, "And user clicks on OK button in \"Export Mutation Cliffs\" dialog", () => clickOn(page, el("OK button in \"Export Mutation Cliffs\" dialog")));
      await session.step(77, "Then the \"Mutation Cliffs\" view should be current", () => viewIsCurrent(page, "Mutation Cliffs"));
      await session.step(78, "And table \"Mutation Cliffs\" should have columns \"Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta, Seq 1 ID, Seq 2 ID\"", () => tableColumns(page, "Mutation Cliffs", "Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta, Seq 1 ID, Seq 2 ID"));
      await session.step(79, "And the table should have 6253 rows", () => rowCount(page, 6253));
      await session.step(80, "And the mutation-cliff export should contain every single-mutation pair from table \"peptides\"", () => mutationPairs(page, "peptides"));
      await session.step(81, "And the mutation-cliff export should preserve \"ID\" values from table \"peptides\"", () => mutationExtraValues(page, "ID", "peptides"));
      await session.step(82, "When user closes the current view", () => closeCurrentView(page));
      await session.step(83, "And user switches to the \"peptides\" table view", () => switchTableView(page, "peptides"));
      await session.step(84, "Then no errors should have been logged", () => noErrors(page));
      await session.step(85, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
