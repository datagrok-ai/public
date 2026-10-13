/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sar/weblogo-selection.feature
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
import {clickOn, collapse, enterInto, expand, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {allOfSelected, clearSelection, noneSelected, onlyOfSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaAtLeastTall, clickArea, clickAreaHolding, noBalloons, noErrors, painted, readingReads, reportsNoError, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Select peptides through WebLogo headers", () => {
  const session = feature(test, "features/sar/weblogo-selection.feature", import.meta.url);
  test("Select peptides through WebLogo headers", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(14, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(15, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(16, "When user picks \"Bio > Analyze > SAR...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > SAR..."));
    await session.step(17, "Then \"Analyze Peptides\" dialog should be visible", () => shouldBe(page, el("\"Analyze Peptides\" dialog"), "visible"));
    await session.step(18, "When user clicks on \"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog", () => clickOn(page, el("\"Adjust clustering parameters\" icon in \"Analyze Peptides\" dialog")));
    await session.step(19, "And user enters \"93\" into \"Similarity Threshold\" input in \"Analyze Peptides\" dialog", () => enterInto(page, "93", el("\"Similarity Threshold\" input in \"Analyze Peptides\" dialog")));
    await session.step(20, "Given user listens for \"peptides-sar-ready\" custom event", () => listenCustom(page, "peptides-sar-ready"));
    await session.step(21, "When user clicks on OK button in \"Analyze Peptides\" dialog", () => clickOn(page, el("OK button in \"Analyze Peptides\" dialog")));
    await session.step(22, "Then the SAR analysis should be ready", () => sarReady(page));
    await session.step(23, "And no rows should be selected", () => noneSelected(page));
    await session.step(24, "Then no errors should have been logged", () => noErrors(page));
    await session.step(25, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("A header glyph selects exactly its matching peptides", async () => {
      await session.step(28, "Then the \"header 2\" area of grid should be at least 100 pixels tall", () => areaAtLeastTall(page, "header 2", el("grid"), 100));
      await session.step(29, "And \"2\" column should have semantic type \"Monomer\"", () => columnSemType(page, "2", "Monomer"));
      await session.step(30, "When user clicks on the \"A at 2\" area of grid", () => clickArea(page, "A at 2", el("grid")));
      await session.step(31, "Then 299 rows should be selected", () => selectedRowCount(page, 299));
      await session.step(32, "And only rows where \"2\" is \"A\" should be selected", () => onlyOfSelected(page, "2", "A"));
      await session.step(33, "And context panel should contain text \"299 selected rows\"", () => shouldContainText(page, el("context panel"), "299 selected rows"));
      await session.step(34, "And context panel should contain text \"Selection Sources\"", () => shouldContainText(page, el("context panel"), "Selection Sources"));
      await session.step(35, "And context panel should contain text \"WebLogo\"", () => shouldContainText(page, el("context panel"), "WebLogo"));
      await session.step(36, "And context panel should contain text \"2:A\"", () => shouldContainText(page, el("context panel"), "2:A"));
      await session.step(37, "And the \"selected monomer-positions\" reading of Sequence Variability Map viewer should be \"\"", () => readingReads(page, "selected monomer-positions", el("Sequence Variability Map viewer"), ""));
      await session.step(38, "And Sequence Variability Map viewer should be painted", () => painted(page, el("Sequence Variability Map viewer")));
      await session.step(39, "And Most Potent Residues viewer should report no error", () => reportsNoError(page, el("Most Potent Residues viewer")));
      await session.step(40, "When user expands Distribution pane in context panel", () => expand(page, el("Distribution pane in context panel")));
      await session.step(41, "Then Distribution pane in context panel should contain text \"Mean difference\"", () => shouldContainText(page, el("Distribution pane in context panel"), "Mean difference"));
      await session.step(42, "And Distribution pane in context panel should contain text \"299 (\"", () => shouldContainText(page, el("Distribution pane in context panel"), "299 ("));
      await session.step(43, "And Distribution pane in context panel should not contain text \"No distribution\"", () => shouldNotContainText(page, el("Distribution pane in context panel"), "No distribution"));
      await session.step(44, "When user expands Selection pane in context panel", () => expand(page, el("Selection pane in context panel")));
      await session.step(45, "Then grid in Selection pane in context panel should show 299 rows", () => showsRows(page, el("grid in Selection pane in context panel"), 299));
      await session.step(46, "When user collapses Selection pane in context panel", () => collapse(page, el("Selection pane in context panel")));
      await session.step(47, "Then no errors should have been logged", () => noErrors(page));
      await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Shift adds the second glyph without losing the first selection", async () => {
      await session.step(51, "When user clicks on the \"Q at 4\" area of grid holding Shift", () => clickAreaHolding(page, "Q at 4", el("grid"), "Shift"));
      await session.step(52, "Then 314 rows should be selected", () => selectedRowCount(page, 314));
      await session.step(53, "And all rows where \"2\" is \"A\" should be selected", () => allOfSelected(page, "2", "A"));
      await session.step(54, "And all rows where \"4\" is \"Q\" should be selected", () => allOfSelected(page, "4", "Q"));
      await session.step(55, "And context panel should contain text \"314 selected rows\"", () => shouldContainText(page, el("context panel"), "314 selected rows"));
      await session.step(56, "And context panel should contain text \"2:A, 4:Q\"", () => shouldContainText(page, el("context panel"), "2:A, 4:Q"));
      await session.step(57, "And the \"selected monomer-positions\" reading of Sequence Variability Map viewer should be \"\"", () => readingReads(page, "selected monomer-positions", el("Sequence Variability Map viewer"), ""));
      await session.step(58, "When user expands Distribution pane in context panel", () => expand(page, el("Distribution pane in context panel")));
      await session.step(59, "Then Distribution pane in context panel should contain text \"Mean difference\"", () => shouldContainText(page, el("Distribution pane in context panel"), "Mean difference"));
      await session.step(60, "And Distribution pane in context panel should contain text \"314 (\"", () => shouldContainText(page, el("Distribution pane in context panel"), "314 ("));
      await session.step(61, "When user expands Selection pane in context panel", () => expand(page, el("Selection pane in context panel")));
      await session.step(62, "Then grid in Selection pane in context panel should show 314 rows", () => showsRows(page, el("grid in Selection pane in context panel"), 314));
      await session.step(63, "When user collapses Selection pane in context panel", () => collapse(page, el("Selection pane in context panel")));
      await session.step(64, "Then no errors should have been logged", () => noErrors(page));
      await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Control toggles off only the second glyph", async () => {
      await session.step(68, "When user clicks on the \"Q at 4\" area of grid holding Control", () => clickAreaHolding(page, "Q at 4", el("grid"), "Control"));
      await session.step(69, "Then 299 rows should be selected", () => selectedRowCount(page, 299));
      await session.step(70, "And only rows where \"2\" is \"A\" should be selected", () => onlyOfSelected(page, "2", "A"));
      await session.step(71, "And context panel should contain text \"2:A\"", () => shouldContainText(page, el("context panel"), "2:A"));
      await session.step(72, "And context panel should not contain text \"4:Q\"", () => shouldNotContainText(page, el("context panel"), "4:Q"));
      await session.step(73, "When user expands Distribution pane in context panel", () => expand(page, el("Distribution pane in context panel")));
      await session.step(74, "Then Distribution pane in context panel should contain text \"299 (\"", () => shouldContainText(page, el("Distribution pane in context panel"), "299 ("));
      await session.step(75, "And Distribution pane in context panel should not contain text \"314 (\"", () => shouldNotContainText(page, el("Distribution pane in context panel"), "314 ("));
      await session.step(76, "When user expands Selection pane in context panel", () => expand(page, el("Selection pane in context panel")));
      await session.step(77, "Then grid in Selection pane in context panel should show 299 rows", () => showsRows(page, el("grid in Selection pane in context panel"), 299));
      await session.step(78, "When user collapses Selection pane in context panel", () => collapse(page, el("Selection pane in context panel")));
      await session.step(79, "Then no errors should have been logged", () => noErrors(page));
      await session.step(80, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Clearing the selection clears both dependent panes", async () => {
      await session.step(83, "When user clears the row selection", () => clearSelection(page));
      await session.step(84, "Then no rows should be selected", () => noneSelected(page));
      await session.step(85, "And context panel should contain text \"0 selected rows\"", () => shouldContainText(page, el("context panel"), "0 selected rows"));
      await session.step(86, "And context panel should not contain text \"2:A\"", () => shouldNotContainText(page, el("context panel"), "2:A"));
      await session.step(87, "When user expands Distribution pane in context panel", () => expand(page, el("Distribution pane in context panel")));
      await session.step(88, "Then Distribution pane in context panel should contain text \"No distribution\"", () => shouldContainText(page, el("Distribution pane in context panel"), "No distribution"));
      await session.step(89, "When user expands Selection pane in context panel", () => expand(page, el("Selection pane in context panel")));
      await session.step(90, "Then Selection pane in context panel should contain text \"No compounds selected\"", () => shouldContainText(page, el("Selection pane in context panel"), "No compounds selected"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
      await session.step(92, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
