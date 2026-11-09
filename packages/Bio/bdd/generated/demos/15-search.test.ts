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
import {subsetSpans, viewerGridColumnWidth} from '../../bindings/demo-search.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {enterInto, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterIsExactlyContains, filterPanelHas, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Search", () => {
  const session = feature(test, "features/demos/15-search.feature", import.meta.url);
  test("Find peptides by motif and by similarity", {tag: ["@guide", "@demo", "@viewport:1440x810"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(21, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the grid shows only the columns \"sequence, activity\"", () => gridShowsColumns(page, "sequence, activity"));
    await session.step(25, "And the \"sequence\" column of the grid is 300 px wide", () => gridColumnWidth(page, "sequence", 300));
    await session.step(26, "When user picks \"Bio > Search > Subsequence Search ...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Subsequence Search ..."));
    await session.step(27, "Then the filter panel should have a filter on \"sequence\" column", () => filterPanelHas(page, "sequence"));
    await session.step(28, "When user enters \"NEMTK\" into \"Substructure\" input in filters viewer", () => enterInto(page, "NEMTK", el("\"Substructure\" input in filters viewer")));
    await session.step(29, "Then 10 rows should pass the filter", () => filterPasses(page, 10));
    await session.step(30, "And the filter should pass exactly the rows where \"sequence\" contains \"NEMTK\"", () => filterIsExactlyContains(page, "sequence", "NEMTK"));
    await session.step(31, "When the video says \"Only the peptides that contain the motif remain\"", () => videoSays(page, "Only the peptides that contain the motif remain"));
    await session.step(32, "And user clicks on the \"cell 42 of sequence\" area of grid", () => clickArea(page, "cell 42 of sequence", el("grid")));
    await session.step(33, "Then row 42 should be current", () => currentRowIs(page, 42));
    await session.step(34, "When user picks \"Bio > Search > Similarity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Similarity Search"));
    await session.step(35, "Then \"Sequence Similarity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Similarity Search\" viewer"), "visible"));
    await session.step(36, "And the \"target row\" reading of \"Sequence Similarity Search\" viewer should be 41", () => readingIs(page, "target row", el("\"Sequence Similarity Search\" viewer"), 41));
    await session.step(37, "Given user sets \"limit\" property of \"Sequence Similarity Search\" viewer to \"6\"", () => setProperty(page, "limit", el("\"Sequence Similarity Search\" viewer"), "6"));
    await session.step(38, "Then the \"neighbours\" reading of \"Sequence Similarity Search\" viewer should be 7", () => readingIs(page, "neighbours", el("\"Sequence Similarity Search\" viewer"), 7));
    await session.step(39, "Given the \"similar (sequence)\" column of the \"Sequence Similarity Search\" viewer is 450 px wide", () => viewerGridColumnWidth(page, "similar (sequence)", "Sequence Similarity Search", 450));
    await session.step(40, "When the video says \"Its nearest neighbours: residues that differ stay in color\"", () => videoSays(page, "Its nearest neighbours: residues that differ stay in color"));
    await session.step(41, "And user clicks on the \"cell 55 of sequence\" area of grid", () => clickArea(page, "cell 55 of sequence", el("grid")));
    await session.step(42, "Then the \"target row\" reading of \"Sequence Similarity Search\" viewer should be 54", () => readingIs(page, "target row", el("\"Sequence Similarity Search\" viewer"), 54));
    await session.step(43, "And the \"neighbours\" reading of \"Sequence Similarity Search\" viewer should be 7", () => readingIs(page, "neighbours", el("\"Sequence Similarity Search\" viewer"), 7));
    await session.step(44, "When user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"));
    await session.step(45, "Then \"Sequence Diversity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "visible"));
    await session.step(46, "And the subset of \"Sequence Diversity Search\" viewer should span 3 values of the \"cluster\" column", () => subsetSpans(page, el("\"Sequence Diversity Search\" viewer"), 3, "cluster"));
    await session.step(47, "When the video says \"The most different peptides: all three families\"", () => videoSays(page, "The most different peptides: all three families"));
    await session.step(48, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
