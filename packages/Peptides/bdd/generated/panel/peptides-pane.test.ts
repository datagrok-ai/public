/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/panel/peptides-pane.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {previewScaling} from '../../bindings/activity.js';
import {peptidesInitialized, sequenceSelection} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, collapse, expand, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, noneSelected, rowCount, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hasArea, noBalloons, noErrors, painted, readingBetween, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Peptide column information and SAR parameters", () => {
  const session = feature(test, "features/panel/peptides-pane.feature", import.meta.url);
  test("Peptide column information and SAR parameters", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(12, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(13, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(14, "Then the table should have 647 rows", () => rowCount(page, 647));
    await session.step(15, "And \"AlignedSequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AlignedSequence", "Macromolecule"));
    await session.step(16, "Then no errors should have been logged", () => noErrors(page));
    await session.step(17, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("The sequence renderer and column information are available from the header", async () => {
      await session.step(20, "Then \"AlignedSequence\" column should have tag \"cell.renderer\" equal to \"sequence\"", () => columnTag(page, "AlignedSequence", "cell.renderer", "sequence"));
      await session.step(21, "And the \"cell type of AlignedSequence\" reading of grid should be \"sequence\"", () => readingReads(page, "cell type of AlignedSequence", el("grid"), "sequence"));
      await session.step(22, "When user clicks on the \"header AlignedSequence\" area of grid", () => clickArea(page, "header AlignedSequence", el("grid")));
      await session.step(23, "Then the context panel should show \"AlignedSequence\"", () => contextPanelShows(page, "AlignedSequence"));
      await session.step(24, "And Details pane in context panel should be visible", () => shouldBe(page, el("Details pane in context panel"), "visible"));
      await session.step(25, "And Peptides pane in context panel should be visible", () => shouldBe(page, el("Peptides pane in context panel"), "visible"));
      await session.step(26, "When user expands Details pane in context panel", () => expand(page, el("Details pane in context panel")));
      await session.step(27, "Then Details pane in context panel should be expanded", () => shouldBe(page, el("Details pane in context panel"), "expanded"));
      await session.step(28, "And Details pane in context panel should contain text \"Data type\"", () => shouldContainText(page, el("Details pane in context panel"), "Data type"));
      await session.step(29, "And Details pane in context panel should contain text \"Semantic type\"", () => shouldContainText(page, el("Details pane in context panel"), "Semantic type"));
      await session.step(30, "Then no errors should have been logged", () => noErrors(page));
      await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Peptides pane previews the sequence and offers the SAR parameters", async () => {
      await session.step(34, "When user expands Peptides pane in context panel", () => expand(page, el("Peptides pane in context panel")));
      await session.step(35, "Then Peptides pane in context panel should be expanded", () => shouldBe(page, el("Peptides pane in context panel"), "expanded"));
      await session.step(36, "And \"Launch SAR\" button in Peptides pane should be visible", () => shouldBe(page, el("\"Launch SAR\" button in Peptides pane"), "visible"));
      await session.step(37, "And Activity input in Peptides pane should be visible", () => shouldBe(page, el("Activity input in Peptides pane"), "visible"));
      await session.step(38, "And editor of Activity input in Peptides pane should have text \"IC50\"", () => shouldHaveText(page, el("editor of Activity input in Peptides pane"), "IC50"));
      await session.step(39, "And Scaling input in Peptides pane should be enabled", () => shouldBe(page, el("Scaling input in Peptides pane"), "enabled"));
      await session.step(40, "And Scaling input in Peptides pane should have value \"none\"", () => shouldHaveValue(page, el("Scaling input in Peptides pane"), "none"));
      await session.step(41, "And Clusters input in Peptides pane should be visible", () => shouldBe(page, el("Clusters input in Peptides pane"), "visible"));
      await session.step(42, "And \"Generate clusters\" checkbox in Peptides pane should be checked", () => shouldBe(page, el("\"Generate clusters\" checkbox in Peptides pane"), "checked"));
      await session.step(43, "And WebLogo viewer in Peptides pane should be painted", () => painted(page, el("WebLogo viewer in Peptides pane")));
      await session.step(44, "And WebLogo viewer in Peptides pane should have a \"position 5\" area", () => hasArea(page, el("WebLogo viewer in Peptides pane"), "position 5"));
      await session.step(45, "And the \"rows shown\" reading of WebLogo viewer in Peptides pane should be 647", () => readingIs(page, "rows shown", el("WebLogo viewer in Peptides pane"), 647));
      await session.step(46, "Then no errors should have been logged", () => noErrors(page));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Changing the activity scale rebuilds the histogram with transformed values", async () => {
      await session.step(50, "Then the \"axis max\" reading of Histogram viewer in Peptides pane should be between 0.9 and 1.1", () => readingBetween(page, "axis max", el("Histogram viewer in Peptides pane"), 0.9, 1.1));
      await session.step(51, "And the peptide activity preview should use \"none\" scaling", () => previewScaling(page, "none"));
      await session.step(52, "When user selects \"lg\" in Scaling input in Peptides pane", () => selectIn(page, "lg", el("Scaling input in Peptides pane")));
      await session.step(53, "Then the \"axis max\" reading of Histogram viewer in Peptides pane should be between -7 and 0", () => readingBetween(page, "axis max", el("Histogram viewer in Peptides pane"), -7, 0));
      await session.step(54, "And Histogram viewer in Peptides pane should be painted", () => painted(page, el("Histogram viewer in Peptides pane")));
      await session.step(55, "And the peptide activity preview should use \"lg\" scaling", () => previewScaling(page, "lg"));
      await session.step(56, "When user selects \"none\" in Scaling input in Peptides pane", () => selectIn(page, "none", el("Scaling input in Peptides pane")));
      await session.step(57, "Then the \"axis max\" reading of Histogram viewer in Peptides pane should be between 0.9 and 1.1", () => readingBetween(page, "axis max", el("Histogram viewer in Peptides pane"), 0.9, 1.1));
      await session.step(58, "And the peptide activity preview should use \"none\" scaling", () => previewScaling(page, "none"));
      await session.step(59, "Then no errors should have been logged", () => noErrors(page));
      await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Clicking a WebLogo glyph selects exactly the peptides carrying that monomer", async () => {
      await session.step(63, "Given no rows should be selected", () => noneSelected(page));
      await session.step(64, "When user clicks on the \"monomer T at position 5\" area of WebLogo viewer in Peptides pane", () => clickArea(page, "monomer T at position 5", el("WebLogo viewer in Peptides pane")));
      await session.step(65, "Then 630 rows should be selected", () => selectedRowCount(page, 630));
      await session.step(66, "And only rows with \"T\" at position 5 of \"AlignedSequence\" column should be selected", () => sequenceSelection(page, "T", 5, "AlignedSequence"));
      await session.step(67, "When user clicks on the \"header AlignedSequence\" area of grid", () => clickArea(page, "header AlignedSequence", el("grid")));
      await session.step(68, "Then the context panel should show \"AlignedSequence\"", () => contextPanelShows(page, "AlignedSequence"));
      await session.step(69, "When user expands Peptides pane in context panel", () => expand(page, el("Peptides pane in context panel")));
      await session.step(70, "Then the \"rows selected\" reading of WebLogo viewer in Peptides pane should be 630", () => readingIs(page, "rows selected", el("WebLogo viewer in Peptides pane"), 630));
      await session.step(71, "When user clears the row selection", () => clearSelection(page));
      await session.step(72, "Then no rows should be selected", () => noneSelected(page));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
      await session.step(74, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Choosing a clusters column turns cluster generation off, and back", async () => {
      await session.step(77, "When user clicks on the \"header AlignedSequence\" area of grid", () => clickArea(page, "header AlignedSequence", el("grid")));
      await session.step(78, "Then the context panel should show \"AlignedSequence\"", () => contextPanelShows(page, "AlignedSequence"));
      await session.step(80, "When user collapses Details pane in context panel", () => collapse(page, el("Details pane in context panel")));
      await session.step(81, "And user expands Peptides pane in context panel", () => expand(page, el("Peptides pane in context panel")));
      await session.step(82, "Then \"Generate clusters\" checkbox in Peptides pane should be checked", () => shouldBe(page, el("\"Generate clusters\" checkbox in Peptides pane"), "checked"));
      await session.step(83, "When user selects \"ID\" in Clusters input in Peptides pane", () => selectIn(page, "ID", el("Clusters input in Peptides pane")));
      await session.step(84, "Then \"Generate clusters\" checkbox in Peptides pane should be unchecked", () => shouldBe(page, el("\"Generate clusters\" checkbox in Peptides pane"), "unchecked"));
      await session.step(85, "When user checks \"Generate clusters\" checkbox in Peptides pane", () => check(page, el("\"Generate clusters\" checkbox in Peptides pane")));
      await session.step(86, "Then \"Generate clusters\" checkbox in Peptides pane should be checked", () => shouldBe(page, el("\"Generate clusters\" checkbox in Peptides pane"), "checked"));
      await session.step(87, "And editor of Clusters input in Peptides pane should have text \"\"", () => shouldHaveText(page, el("editor of Clusters input in Peptides pane"), ""));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
      await session.step(89, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Without an analysis the Manual Alignment pane only says it needs one", async () => {
      await session.step(92, "When user picks \"Bio > Transform > Split to Monomers...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Split to Monomers..."));
      await session.step(93, "Then \"Split to Monomers\" dialog should be visible", () => shouldBe(page, el("\"Split to Monomers\" dialog"), "visible"));
      await session.step(94, "When user clicks on OK button in \"Split to Monomers\" dialog", () => clickOn(page, el("OK button in \"Split to Monomers\" dialog")));
      await session.step(95, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(96, "And \"2\" column should have semantic type \"Monomer\"", () => columnSemType(page, "2", "Monomer"));
      await session.step(97, "When user clicks on the \"cell 2 of 2\" area of grid", () => clickArea(page, "cell 2 of 2", el("grid")));
      await session.step(98, "Then \"Manual Alignment\" pane in context panel should be visible", () => shouldBe(page, el("\"Manual Alignment\" pane in context panel"), "visible"));
      await session.step(99, "When user expands \"Manual Alignment\" pane in context panel", () => expand(page, el("\"Manual Alignment\" pane in context panel")));
      await session.step(100, "Then \"Manual Alignment\" pane in context panel should contain text \"Manual alignment works with peptides analysis\"", () => shouldContainText(page, el("\"Manual Alignment\" pane in context panel"), "Manual alignment works with peptides analysis"));
      await session.step(101, "And Apply button in \"Manual Alignment\" pane should be absent", () => shouldBe(page, el("Apply button in \"Manual Alignment\" pane"), "absent"));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
      await session.step(103, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
