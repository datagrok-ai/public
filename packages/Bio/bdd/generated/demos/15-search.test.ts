/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/15-search.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {dockedRightOfGrid, dockedUnderViewer} from '../../bindings/demo-layout.js';
import {subsetSpans, viewerGridColumnWidth} from '../../bindings/demo-search.js';
import {bioInitialized, sequenceFontSize} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, enterInto, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCategorical, filterIsExactlyContains, filterPanelHas, filterPasses, filterPassesAll, noneOfFiltered} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("4 · A peptide series · Search", () => {
  const session = feature(test, "features/demos/15-search.feature", import.meta.url);
  test("Find peptides by motif and by similarity", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And simple mode is off", () => simpleModeOff(page));
    await session.step(24, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(25, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(26, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(27, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(28, "And the grid shows only the columns \"cluster, sequence, activity\"", () => gridShowsColumns(page, "cluster, sequence, activity"));
    await session.step(29, "And user colors \"cluster\" column categorically:", () => colorCategorical(page, "cluster", [["0","#8e6cc0"],["1","#e9c46a"],["2","#2a9d8f"]]), [["0","#8e6cc0"],["1","#e9c46a"],["2","#2a9d8f"]]);
    await session.step(33, "And the \"cluster\" column of the grid is 60 px wide", () => gridColumnWidth(page, "cluster", 60));
    await session.step(34, "And the rows of the grid are 30 px tall", () => gridRowHeight(page, 30));
    await session.step(35, "And the monomers of the \"sequence\" column are drawn at 22 px", () => sequenceFontSize(page, "sequence", 22));
    await session.step(36, "And the \"sequence\" column of the grid is 500 px wide", () => gridColumnWidth(page, "sequence", 500));
    await session.step(37, "When user picks \"Bio > Search > Subsequence Search ...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Subsequence Search ..."));
    await session.step(38, "Then the filter panel should have a filter on \"sequence\" column", () => filterPanelHas(page, "sequence"));
    await session.step(39, "When user enters \"NEMTK\" into \"Substructure\" input in filters viewer", () => enterInto(page, "NEMTK", el("\"Substructure\" input in filters viewer")));
    await session.step(40, "Then 10 rows should pass the filter", () => filterPasses(page, 10));
    await session.step(41, "And the filter should pass exactly the rows where \"sequence\" contains \"NEMTK\"", () => filterIsExactlyContains(page, "sequence", "NEMTK"));
    await session.step(42, "And no rows where \"cluster\" is \"0\" should pass the filter", () => noneOfFiltered(page, "cluster", "0"));
    await session.step(43, "And no rows where \"cluster\" is \"2\" should pass the filter", () => noneOfFiltered(page, "cluster", "2"));
    await session.step(44, "When the video says \"Only the peptides with the motif remain: one family\"", () => videoSays(page, "Only the peptides with the motif remain: one family"));
    await session.step(45, "And user clicks on the \"cell 42 of sequence\" area of grid", () => clickArea(page, "cell 42 of sequence", el("grid")));
    await session.step(46, "Then row 42 should be current", () => currentRowIs(page, 42));
    await session.step(47, "When user picks \"Bio > Search > Similarity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Similarity Search"));
    await session.step(48, "Then \"Sequence Similarity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Similarity Search\" viewer"), "visible"));
    await session.step(49, "And the \"target row\" reading of \"Sequence Similarity Search\" viewer should be 41", () => readingIs(page, "target row", el("\"Sequence Similarity Search\" viewer"), 41));
    await session.step(50, "Given \"Sequence Similarity Search\" viewer is docked to the right of the grid, 50% of the view wide", () => dockedRightOfGrid(page, el("\"Sequence Similarity Search\" viewer"), 50));
    await session.step(51, "And the \"sequence\" column of the grid is 380 px wide", () => gridColumnWidth(page, "sequence", 380));
    await session.step(52, "And user sets \"limit\" property of \"Sequence Similarity Search\" viewer to \"6\"", () => setProperty(page, "limit", el("\"Sequence Similarity Search\" viewer"), "6"));
    await session.step(53, "Then the \"neighbours\" reading of \"Sequence Similarity Search\" viewer should be 7", () => readingIs(page, "neighbours", el("\"Sequence Similarity Search\" viewer"), 7));
    await session.step(54, "Given the \"similar (sequence)\" column of the \"Sequence Similarity Search\" viewer is 450 px wide", () => viewerGridColumnWidth(page, "similar (sequence)", "Sequence Similarity Search", 450));
    await session.step(55, "When the video says \"Its nearest neighbours: residues that differ stay in color\"", () => videoSays(page, "Its nearest neighbours: residues that differ stay in color"));
    await session.step(56, "And user clicks on the \"cell 55 of sequence\" area of grid", () => clickArea(page, "cell 55 of sequence", el("grid")));
    await session.step(57, "Then the \"target row\" reading of \"Sequence Similarity Search\" viewer should be 54", () => readingIs(page, "target row", el("\"Sequence Similarity Search\" viewer"), 54));
    await session.step(58, "And the \"neighbours\" reading of \"Sequence Similarity Search\" viewer should be 7", () => readingIs(page, "neighbours", el("\"Sequence Similarity Search\" viewer"), 7));
    await session.step(59, "When user clears \"Substructure\" input in filters viewer", () => clearField(page, el("\"Substructure\" input in filters viewer")));
    await session.step(60, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(61, "When user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"));
    await session.step(62, "Then \"Sequence Diversity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "visible"));
    await session.step(63, "Given \"Sequence Diversity Search\" viewer is docked under \"Sequence Similarity Search\" viewer", () => dockedUnderViewer(page, el("\"Sequence Diversity Search\" viewer"), el("\"Sequence Similarity Search\" viewer")));
    await session.step(64, "And user sets \"limit\" property of \"Sequence Diversity Search\" viewer to \"3\"", () => setProperty(page, "limit", el("\"Sequence Diversity Search\" viewer"), "3"));
    await session.step(65, "Then the \"subset size\" reading of \"Sequence Diversity Search\" viewer should be 3", () => readingIs(page, "subset size", el("\"Sequence Diversity Search\" viewer"), 3));
    await session.step(66, "And the subset of \"Sequence Diversity Search\" viewer should span 3 values of the \"cluster\" column", () => subsetSpans(page, el("\"Sequence Diversity Search\" viewer"), 3, "cluster"));
    await session.step(67, "When the video says \"The three most different peptides: one per family\"", () => videoSays(page, "The three most different peptides: one per family"));
    await session.step(68, "And user clicks on the \"row header 1\" area of third grid viewer", () => clickArea(page, "row header 1", el("third grid viewer")));
    await session.step(69, "And user clicks on the \"row header 2\" area of third grid viewer", () => clickArea(page, "row header 2", el("third grid viewer")));
    await session.step(70, "And user clicks on the \"row header 3\" area of third grid viewer", () => clickArea(page, "row header 3", el("third grid viewer")));
    await session.step(71, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(72, "And no errors should have been logged", () => noErrors(page));
  });
});
