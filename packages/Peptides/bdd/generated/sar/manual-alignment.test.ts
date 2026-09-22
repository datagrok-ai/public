/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sar/manual-alignment.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {alignedPositions} from '../../bindings/alignment.js';
import {peptidesInitialized, sarReady} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, hoverOver, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag, currentRowIs, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {currentColumnIs, onlyOfSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, readingIs, reportsNoError, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Manually align a peptide", () => {
  const session = feature(test, "features/sar/manual-alignment.feature", import.meta.url);
  test("Manually align a peptide", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(10, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(11, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(12, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(13, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(14, "And user enters \"93\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "93", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(15, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(16, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(17, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(18, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(19, "Then no errors should have been logged", () => noErrors(page));
    await session.step(20, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("A monomer cell opens the alignment editor for its row", async () => {
      await session.step(23, "Then \"2\" column should have semantic type \"Monomer\"", () => columnSemType(page, "2", "Monomer"));
      await session.step(24, "When user clicks on the \"cell 2 of 2\" area of grid", () => clickArea(page, "cell 2 of 2", el("grid")));
      await session.step(25, "Then row 2 should be current", () => currentRowIs(page, 2));
      await session.step(26, "And the current column should be \"2\"", () => currentColumnIs(page, "2"));
      await session.step(27, "And \"Manual Alignment\" pane in context panel should be visible", () => shouldBe(page, el("\"Manual Alignment\" pane in context panel"), "visible"));
      await session.step(28, "When user expands \"Manual Alignment\" pane in context panel", () => expand(page, el("\"Manual Alignment\" pane in context panel")));
      await session.step(29, "Then Sequence text area in \"Manual Alignment\" pane should have value \"NH2-M-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH\"", () => shouldHaveValue(page, el("Sequence text area in \"Manual Alignment\" pane"), "NH2-M-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"));
      await session.step(30, "And Apply button in \"Manual Alignment\" pane should be visible", () => shouldBe(page, el("Apply button in \"Manual Alignment\" pane"), "visible"));
      await session.step(31, "When user hovers over Sequence text area in \"Manual Alignment\" pane", () => hoverOver(page, el("Sequence text area in \"Manual Alignment\" pane")));
      await session.step(32, "Then Reset button in \"Manual Alignment\" pane should be visible", () => shouldBe(page, el("Reset button in \"Manual Alignment\" pane"), "visible"));
      await session.step(33, "Then no errors should have been logged", () => noErrors(page));
      await session.step(34, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Apply changes the intended position without shifting adjacent monomers", async () => {
      await session.step(37, "Then the \"count of cell M at 2\" reading of Sequence Variability Map viewer should be 9", () => readingIs(page, "count of cell M at 2", el("Sequence Variability Map viewer"), 9));
      await session.step(38, "And the \"count of cell V at 2\" reading of Sequence Variability Map viewer should be 3", () => readingIs(page, "count of cell V at 2", el("Sequence Variability Map viewer"), 3));
      await session.step(39, "When user enters \"NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH\" into Sequence text area in \"Manual Alignment\" pane", () => enterInto(page, "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH", el("Sequence text area in \"Manual Alignment\" pane")));
      await session.step(40, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
      await session.step(41, "When user clicks on Apply button in \"Manual Alignment\" pane", () => clickOn(page, el("Apply button in \"Manual Alignment\" pane")));
      await session.step(42, "Then the SAR analysis should be ready", () => sarReady(page));
      await session.step(43, "Then the value of \"AlignedSequence\" column in row 2 should be \"NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH\"", () => valueInRow(page, "AlignedSequence", 2, "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"));
      await session.step(44, "And the value of \"1\" column in row 2 should be \"NH2\"", () => valueInRow(page, "1", 2, "NH2"));
      await session.step(45, "And the value of \"2\" column in row 2 should be \"V\"", () => valueInRow(page, "2", 2, "V"));
      await session.step(46, "And the value of \"3\" column in row 2 should be \"A\"", () => valueInRow(page, "3", 2, "A"));
      await session.step(47, "And the value of \"17\" column in row 2 should be \"COOH\"", () => valueInRow(page, "17", 2, "COOH"));
      await session.step(48, "And the split monomer columns of row 2 should match its peptide sequence", () => alignedPositions(page, 2));
      await session.step(49, "And \"AlignedSequence\" column should have tag \"cell.renderer\" equal to \"sequence\"", () => columnTag(page, "AlignedSequence", "cell.renderer", "sequence"));
      await session.step(50, "And the open tableview should have 1 Sequence Variability Map viewer", () => viewerCount(page, 1, "Sequence Variability Map"));
      await session.step(51, "And the open tableview should have 1 Most Potent Residues viewer", () => viewerCount(page, 1, "Most Potent Residues"));
      await session.step(52, "And Sequence Variability Map viewer should report no error", () => reportsNoError(page, el("Sequence Variability Map viewer")));
      await session.step(53, "And the \"count of cell M at 2\" reading of Sequence Variability Map viewer should be 8", () => readingIs(page, "count of cell M at 2", el("Sequence Variability Map viewer"), 8));
      await session.step(54, "And the \"count of cell V at 2\" reading of Sequence Variability Map viewer should be 4", () => readingIs(page, "count of cell V at 2", el("Sequence Variability Map viewer"), 4));
      await session.step(55, "Then no errors should have been logged", () => noErrors(page));
      await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Reset keeps the applied sequence and its monomer columns", async () => {
      await session.step(59, "When user expands \"Manual Alignment\" pane in context panel", () => expand(page, el("\"Manual Alignment\" pane in context panel")));
      await session.step(60, "Then Sequence text area in \"Manual Alignment\" pane should have value \"NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH\"", () => shouldHaveValue(page, el("Sequence text area in \"Manual Alignment\" pane"), "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"));
      await session.step(61, "When user hovers over Sequence text area in \"Manual Alignment\" pane", () => hoverOver(page, el("Sequence text area in \"Manual Alignment\" pane")));
      await session.step(62, "When user clicks on Reset button in \"Manual Alignment\" pane", () => clickOn(page, el("Reset button in \"Manual Alignment\" pane")));
      await session.step(63, "Then Sequence text area in \"Manual Alignment\" pane should have value \"NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH\"", () => shouldHaveValue(page, el("Sequence text area in \"Manual Alignment\" pane"), "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"));
      await session.step(64, "And the value of \"AlignedSequence\" column in row 2 should be \"NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH\"", () => valueInRow(page, "AlignedSequence", 2, "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"));
      await session.step(65, "And the split monomer columns of row 2 should match its peptide sequence", () => alignedPositions(page, 2));
      await session.step(66, "Then no errors should have been logged", () => noErrors(page));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Reset discards an unsaved edit without changing the applied alignment", async () => {
      await session.step(70, "When user enters \"NH2-V-A-X-T-T-Y-K-N-Y-R-N-N-L-L--COOH\" into Sequence text area in \"Manual Alignment\" pane", () => enterInto(page, "NH2-V-A-X-T-T-Y-K-N-Y-R-N-N-L-L--COOH", el("Sequence text area in \"Manual Alignment\" pane")));
      await session.step(71, "Then Sequence text area in \"Manual Alignment\" pane should have value \"NH2-V-A-X-T-T-Y-K-N-Y-R-N-N-L-L--COOH\"", () => shouldHaveValue(page, el("Sequence text area in \"Manual Alignment\" pane"), "NH2-V-A-X-T-T-Y-K-N-Y-R-N-N-L-L--COOH"));
      await session.step(72, "When user clicks on Reset button in \"Manual Alignment\" pane", () => clickOn(page, el("Reset button in \"Manual Alignment\" pane")));
      await session.step(73, "Then Sequence text area in \"Manual Alignment\" pane should have value \"NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH\"", () => shouldHaveValue(page, el("Sequence text area in \"Manual Alignment\" pane"), "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"));
      await session.step(74, "And the value of \"AlignedSequence\" column in row 2 should be \"NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH\"", () => valueInRow(page, "AlignedSequence", 2, "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"));
      await session.step(75, "And the split monomer columns of row 2 should match its peptide sequence", () => alignedPositions(page, 2));
      await session.step(76, "Then no errors should have been logged", () => noErrors(page));
      await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The edited alignment remains selectable from its WebLogo header", async () => {
      await session.step(80, "When user clicks on the \"V at 2\" area of grid", () => clickArea(page, "V at 2", el("grid")));
      await session.step(81, "Then only rows where \"2\" is \"V\" should be selected", () => onlyOfSelected(page, "2", "V"));
      await session.step(82, "When user expands Distribution pane in context panel", () => expand(page, el("Distribution pane in context panel")));
      await session.step(83, "Then Distribution pane in context panel should contain text \"Mean difference\"", () => shouldContainText(page, el("Distribution pane in context panel"), "Mean difference"));
      await session.step(84, "When user expands Selection pane in context panel", () => expand(page, el("Selection pane in context panel")));
      await session.step(85, "Then Selection pane in context panel should not contain text \"No compounds selected\"", () => shouldNotContainText(page, el("Selection pane in context panel"), "No compounds selected"));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
      await session.step(87, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
