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
import {expand, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, noneSelected, rowCount, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hasArea, noBalloons, noErrors, painted, readingBetween, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Peptide column information and SAR parameters", () => {
  const session = feature(test, "features/panel/peptides-pane.feature", import.meta.url);
  test("Peptide column information and SAR parameters", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(9, "And user opens peptides dataset", () => openDataset(page, ds("peptides")));
    await session.step(10, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(11, "Then the table should have 647 rows", () => rowCount(page, 647));
    await session.step(12, "And \"AlignedSequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AlignedSequence", "Macromolecule"));
    await session.step(13, "Then no errors should have been logged", () => noErrors(page));
    await session.step(14, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("The sequence renderer and column information are available from the header", async () => {
      await session.step(17, "Then \"AlignedSequence\" column should have tag \"cell.renderer\" equal to \"sequence\"", () => columnTag(page, "AlignedSequence", "cell.renderer", "sequence"));
      await session.step(18, "And the \"cell type of AlignedSequence\" reading of grid should be \"sequence\"", () => readingReads(page, "cell type of AlignedSequence", el("grid"), "sequence"));
      await session.step(19, "When user clicks on the \"header AlignedSequence\" area of grid", () => clickArea(page, "header AlignedSequence", el("grid")));
      await session.step(20, "Then the context panel should show \"AlignedSequence\"", () => contextPanelShows(page, "AlignedSequence"));
      await session.step(21, "And Details pane in context panel should be visible", () => shouldBe(page, el("Details pane in context panel"), "visible"));
      await session.step(22, "And Peptides pane in context panel should be visible", () => shouldBe(page, el("Peptides pane in context panel"), "visible"));
      await session.step(23, "When user expands Details pane in context panel", () => expand(page, el("Details pane in context panel")));
      await session.step(24, "Then Details pane in context panel should be expanded", () => shouldBe(page, el("Details pane in context panel"), "expanded"));
      await session.step(25, "And Details pane in context panel should contain text \"Data type\"", () => shouldContainText(page, el("Details pane in context panel"), "Data type"));
      await session.step(26, "And Details pane in context panel should contain text \"Semantic type\"", () => shouldContainText(page, el("Details pane in context panel"), "Semantic type"));
      await session.step(27, "Then no errors should have been logged", () => noErrors(page));
      await session.step(28, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Peptides pane previews the sequence and offers the SAR parameters", async () => {
      await session.step(31, "When user expands Peptides pane in context panel", () => expand(page, el("Peptides pane in context panel")));
      await session.step(32, "Then Peptides pane in context panel should be expanded", () => shouldBe(page, el("Peptides pane in context panel"), "expanded"));
      await session.step(33, "And \"Launch SAR\" button in Peptides pane should be visible", () => shouldBe(page, el("\"Launch SAR\" button in Peptides pane"), "visible"));
      await session.step(34, "And Activity input in Peptides pane should be visible", () => shouldBe(page, el("Activity input in Peptides pane"), "visible"));
      await session.step(35, "And editor of Activity input in Peptides pane should have text \"IC50\"", () => shouldHaveText(page, el("editor of Activity input in Peptides pane"), "IC50"));
      await session.step(36, "And Scaling input in Peptides pane should be enabled", () => shouldBe(page, el("Scaling input in Peptides pane"), "enabled"));
      await session.step(37, "And Scaling input in Peptides pane should have value \"none\"", () => shouldHaveValue(page, el("Scaling input in Peptides pane"), "none"));
      await session.step(38, "And Clusters input in Peptides pane should be visible", () => shouldBe(page, el("Clusters input in Peptides pane"), "visible"));
      await session.step(39, "And \"Generate clusters\" checkbox in Peptides pane should be checked", () => shouldBe(page, el("\"Generate clusters\" checkbox in Peptides pane"), "checked"));
      await session.step(40, "And WebLogo viewer in Peptides pane should be painted", () => painted(page, el("WebLogo viewer in Peptides pane")));
      await session.step(41, "And WebLogo viewer in Peptides pane should have a \"position 5\" area", () => hasArea(page, el("WebLogo viewer in Peptides pane"), "position 5"));
      await session.step(42, "And the \"rows shown\" reading of WebLogo viewer in Peptides pane should be 647", () => readingIs(page, "rows shown", el("WebLogo viewer in Peptides pane"), 647));
      await session.step(43, "Then no errors should have been logged", () => noErrors(page));
      await session.step(44, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Changing the activity scale rebuilds the histogram with transformed values", async () => {
      await session.step(47, "Then the \"axis max\" reading of Histogram viewer in Peptides pane should be between 0.9 and 1.1", () => readingBetween(page, "axis max", el("Histogram viewer in Peptides pane"), 0.9, 1.1));
      await session.step(48, "And the peptide activity preview should use \"none\" scaling", () => previewScaling(page, "none"));
      await session.step(49, "When user selects \"lg\" in Scaling input in Peptides pane", () => selectIn(page, "lg", el("Scaling input in Peptides pane")));
      await session.step(50, "Then the \"axis max\" reading of Histogram viewer in Peptides pane should be between -7 and 0", () => readingBetween(page, "axis max", el("Histogram viewer in Peptides pane"), -7, 0));
      await session.step(51, "And Histogram viewer in Peptides pane should be painted", () => painted(page, el("Histogram viewer in Peptides pane")));
      await session.step(52, "And the peptide activity preview should use \"lg\" scaling", () => previewScaling(page, "lg"));
      await session.step(53, "When user selects \"none\" in Scaling input in Peptides pane", () => selectIn(page, "none", el("Scaling input in Peptides pane")));
      await session.step(54, "Then the \"axis max\" reading of Histogram viewer in Peptides pane should be between 0.9 and 1.1", () => readingBetween(page, "axis max", el("Histogram viewer in Peptides pane"), 0.9, 1.1));
      await session.step(55, "And the peptide activity preview should use \"none\" scaling", () => previewScaling(page, "none"));
      await session.step(56, "Then no errors should have been logged", () => noErrors(page));
      await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Clicking a WebLogo glyph selects exactly the peptides carrying that monomer", async () => {
      await session.step(60, "Given no rows should be selected", () => noneSelected(page));
      await session.step(61, "When user clicks on the \"monomer T at position 5\" area of WebLogo viewer in Peptides pane", () => clickArea(page, "monomer T at position 5", el("WebLogo viewer in Peptides pane")));
      await session.step(62, "Then 630 rows should be selected", () => selectedRowCount(page, 630));
      await session.step(63, "And only rows with \"T\" at position 5 of \"AlignedSequence\" column should be selected", () => sequenceSelection(page, "T", 5, "AlignedSequence"));
      await session.step(64, "When user clicks on the \"header AlignedSequence\" area of grid", () => clickArea(page, "header AlignedSequence", el("grid")));
      await session.step(65, "Then the context panel should show \"AlignedSequence\"", () => contextPanelShows(page, "AlignedSequence"));
      await session.step(66, "When user expands Peptides pane in context panel", () => expand(page, el("Peptides pane in context panel")));
      await session.step(67, "Then the \"rows selected\" reading of WebLogo viewer in Peptides pane should be 630", () => readingIs(page, "rows selected", el("WebLogo viewer in Peptides pane"), 630));
      await session.step(68, "When user clears the row selection", () => clearSelection(page));
      await session.step(69, "Then no rows should be selected", () => noneSelected(page));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
      await session.step(71, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
