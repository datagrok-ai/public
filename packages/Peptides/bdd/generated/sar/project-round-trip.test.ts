/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sar/project-round-trip.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {analysisScaling} from '../../bindings/activity.js';
import {peptidesInitialized, sarReady, sarSetting} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, hasColumn, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, noneSelected, onlyOfSelected, rowCount, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaAtLeastTall, clickArea, noBalloons, noErrors, readingIs, readingReads, reportsNoError, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Save and reopen a peptide SAR analysis", () => {
  const session = feature(test, "features/sar/project-round-trip.feature", import.meta.url);
  test("Save and reopen a peptide SAR analysis", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(9, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(10, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(11, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(12, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(13, "And user enters \"93\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "93", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(14, "When user selects \"lg\" in Scaling input in \"Analyze Peptides\" dialog", () => selectIn(page, "lg", el("Scaling input in \"Analyze Peptides\" dialog")));
    await session.step(15, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(16, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(17, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(18, "Then no errors should have been logged", () => noErrors(page));
    await session.step(19, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("Reopening restores the data, settings, viewers and WebLogo headers", async () => {
      await session.step(22, "When user clicks on \"Invariant Map\" checkbox in Sequence Variability Map viewer", () => clickOn(page, el("\"Invariant Map\" checkbox in Sequence Variability Map viewer")));
      await session.step(23, "And user clicks on the \"cell M at 2\" area of Sequence Variability Map viewer", () => clickArea(page, "cell M at 2", el("Sequence Variability Map viewer")));
      await session.step(24, "Then 9 rows should be selected", () => selectedRowCount(page, 9));
      await session.step(25, "And only rows where \"2\" is \"M\" should be selected", () => onlyOfSelected(page, "2", "M"));
      await session.step(26, "When user saves the current view as project \"bdd-peptides-sar-roundtrip\"", () => saveAsProject(page, "bdd-peptides-sar-roundtrip"));
      await session.step(27, "And user closes all views", () => closeAllViews(page));
      await session.step(28, "And user opens the \"bdd-peptides-sar-roundtrip\" project", () => openProject(page, "bdd-peptides-sar-roundtrip"));
      await session.step(29, "Then the table should have 647 rows", () => rowCount(page, 647));
      await session.step(30, "And \"AlignedSequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AlignedSequence", "Macromolecule"));
      await session.step(31, "And \"2\" column should have semantic type \"Monomer\"", () => columnSemType(page, "2", "Monomer"));
      await session.step(32, "And the table should have a column \"17\"", () => hasColumn(page, "17"));
      await session.step(33, "And the table should not have a column \"18\"", () => hasNoColumn(page, "18"));
      await session.step(34, "And the SAR setting \"activityScaling\" should be \"lg\"", () => sarSetting(page, "activityScaling", "lg"));
      await session.step(35, "And the SAR activity column should use \"lg\" scaling", () => analysisScaling(page, "lg"));
      await session.step(36, "And the open tableview should have 1 Sequence Variability Map viewer", () => viewerCount(page, 1, "Sequence Variability Map"));
      await session.step(37, "And the open tableview should have 1 Most Potent Residues viewer", () => viewerCount(page, 1, "Most Potent Residues"));
      await session.step(38, "And the open tableview should have 1 MCL viewer", () => viewerCount(page, 1, "MCL"));
      await session.step(39, "And the open tableview should have 1 Logo Summary Table viewer", () => viewerCount(page, 1, "Logo Summary Table"));
      await session.step(40, "And the \"positions\" reading of Sequence Variability Map viewer should be 17", () => readingIs(page, "positions", el("Sequence Variability Map viewer"), 17));
      await session.step(41, "And the \"activity scaling\" reading of Sequence Variability Map viewer should be \"lg\"", () => readingReads(page, "activity scaling", el("Sequence Variability Map viewer"), "lg"));
      await session.step(42, "And the \"header 2\" area of grid should be at least 100 pixels tall", () => areaAtLeastTall(page, "header 2", el("grid"), 100));
      await session.step(43, "And Sequence Variability Map viewer should report no error", () => reportsNoError(page, el("Sequence Variability Map viewer")));
      await session.step(44, "Then no errors should have been logged", () => noErrors(page));
      await session.step(45, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The restored analysis responds to a new monomer selection", async () => {
      await session.step(48, "When user clears the row selection", () => clearSelection(page));
      await session.step(49, "Then no rows should be selected", () => noneSelected(page));
      await session.step(50, "When user clicks on \"Invariant Map\" checkbox in Sequence Variability Map viewer", () => clickOn(page, el("\"Invariant Map\" checkbox in Sequence Variability Map viewer")));
      await session.step(51, "And user clicks on the \"cell A at 2\" area of Sequence Variability Map viewer", () => clickArea(page, "cell A at 2", el("Sequence Variability Map viewer")));
      await session.step(52, "Then 299 rows should be selected", () => selectedRowCount(page, 299));
      await session.step(53, "And only rows where \"2\" is \"A\" should be selected", () => onlyOfSelected(page, "2", "A"));
      await session.step(54, "And Distribution pane in context panel should be present", () => shouldBe(page, el("Distribution pane in context panel"), "present"));
      await session.step(55, "And \"Mutation Cliffs pairs\" pane in context panel should be present", () => shouldBe(page, el("\"Mutation Cliffs pairs\" pane in context panel"), "present"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
      await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The restored WebLogo header still selects matching peptides", async () => {
      await session.step(60, "When user clears the row selection", () => clearSelection(page));
      await session.step(61, "And user clicks on the \"A at 2\" area of grid", () => clickArea(page, "A at 2", el("grid")));
      await session.step(62, "Then 299 rows should be selected", () => selectedRowCount(page, 299));
      await session.step(63, "And only rows where \"2\" is \"A\" should be selected", () => onlyOfSelected(page, "2", "A"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
      await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
