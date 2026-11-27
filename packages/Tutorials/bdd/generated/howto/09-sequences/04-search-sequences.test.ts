/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/09-sequences/04-search-sequences.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bioReady, monomersDrawnAt, subsetSpans} from '../../../bindings/howto-sequences.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, enterInto, hoverOver, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCategorical, filterIsExactlyContains, filterPanelHas, filterPasses, filterPassesAll, noneOfFiltered} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Sequences", () => {
  const session = feature(test, "features/howto/09-sequences/04-search-sequences.feature", import.meta.url);
  test("Search sequences by motif, similarity and diversity", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(21, "And the Bio package is initialized", () => bioReady(page));
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the grid shows only the columns \"cluster, sequence, activity\"", () => gridShowsColumns(page, "cluster, sequence, activity"));
    await session.step(25, "And user colors \"cluster\" column categorically:", () => colorCategorical(page, "cluster", [["0","#8e6cc0"],["1","#e9c46a"],["2","#2a9d8f"]]), [["0","#8e6cc0"],["1","#e9c46a"],["2","#2a9d8f"]]);
    await session.step(29, "And the \"cluster\" column of the grid is 60 px wide", () => gridColumnWidth(page, "cluster", 60));
    await session.step(30, "And the rows of the grid are 30 px tall", () => gridRowHeight(page, 30));
    await session.step(31, "And the monomers of the \"sequence\" column are drawn at 22 px", () => monomersDrawnAt(page, "sequence", 22));
    await session.step(32, "And the \"sequence\" column of the grid is 420 px wide", () => gridColumnWidth(page, "sequence", 420));
    await session.step(33, "When the video says \"99 peptides; the first column, colored, is each one's family\"", () => videoSays(page, "99 peptides; the first column, colored, is each one's family"));
    await session.step(35, "And user picks \"Bio > Search > Subsequence Search ...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Subsequence Search ..."), undefined, "Bio > Search > Subsequence Search: find a motif");
    await session.step(36, "Then the filter panel should have a filter on \"sequence\" column", () => filterPanelHas(page, "sequence"));
    await session.step(38, "When user enters \"NEMTK\" into \"Substructure\" input in filters viewer", () => enterInto(page, "NEMTK", el("\"Substructure\" input in filters viewer")), undefined, "Type a motif into the sequence filter");
    await session.step(39, "Then 10 rows should pass the filter", () => filterPasses(page, 10));
    await session.step(40, "And the filter should pass exactly the rows where \"sequence\" contains \"NEMTK\"", () => filterIsExactlyContains(page, "sequence", "NEMTK"));
    await session.step(41, "And no rows where \"cluster\" is \"0\" should pass the filter", () => noneOfFiltered(page, "cluster", "0"));
    await session.step(42, "And no rows where \"cluster\" is \"2\" should pass the filter", () => noneOfFiltered(page, "cluster", "2"));
    await session.step(43, "When the video says \"Ten peptides hold the motif, all of one family\"", () => videoSays(page, "Ten peptides hold the motif, all of one family"));
    await session.step(44, "And user hovers over \"Filtered: 10\" text in status bar", () => hoverOver(page, el("\"Filtered: 10\" text in status bar")));
    await session.step(46, "And user clears \"Substructure\" input in filters viewer", () => clearField(page, el("\"Substructure\" input in filters viewer")), undefined, "Clear the motif: all 99 are back");
    await session.step(47, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(49, "When user clicks on the \"cell 5 of sequence\" area of grid", () => clickArea(page, "cell 5 of sequence", el("grid")), undefined, "Click a peptide");
    await session.step(50, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(52, "When user picks \"Bio > Search > Similarity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Similarity Search"), undefined, "Bio > Search > Similarity Search: its nearest neighbours");
    await session.step(53, "Then \"Sequence Similarity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Similarity Search\" viewer"), "visible"));
    await session.step(54, "And the \"target row\" reading of \"Sequence Similarity Search\" viewer should be 4", () => readingIs(page, "target row", el("\"Sequence Similarity Search\" viewer"), 4));
    await session.step(55, "And the \"neighbours\" reading of \"Sequence Similarity Search\" viewer should be 11", () => readingIs(page, "neighbours", el("\"Sequence Similarity Search\" viewer"), 11));
    await session.step(56, "When the video says \"Most similar first; what differs stays in color\"", () => videoSays(page, "Most similar first; what differs stays in color"));
    await session.step(58, "And user clicks on the \"cell 10 of sequence\" area of grid", () => clickArea(page, "cell 10 of sequence", el("grid")), undefined, "Click another peptide: the list follows");
    await session.step(59, "Then the \"target row\" reading of \"Sequence Similarity Search\" viewer should be 9", () => readingIs(page, "target row", el("\"Sequence Similarity Search\" viewer"), 9));
    await session.step(61, "When user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"), undefined, "Bio > Search > Diversity Search: the most different peptides");
    await session.step(62, "Then \"Sequence Diversity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "visible"));
    await session.step(63, "And the \"subset size\" reading of \"Sequence Diversity Search\" viewer should be 10", () => readingIs(page, "subset size", el("\"Sequence Diversity Search\" viewer"), 10));
    await session.step(64, "And the peptides \"Sequence Diversity Search\" viewer lists should come from 3 different \"cluster\" values", () => subsetSpans(page, el("\"Sequence Diversity Search\" viewer"), 3, "cluster"));
    await session.step(65, "When the video says \"Ten peptides unlike each other, from all three families\"", () => videoSays(page, "Ten peptides unlike each other, from all three families"));
    await session.step(67, "And user clicks on the \"row header 1\" area of third grid viewer", () => clickArea(page, "row header 1", el("third grid viewer")), undefined, "Click one to find it in the table");
    await session.step(68, "And user clicks on the \"row header 2\" area of third grid viewer", () => clickArea(page, "row header 2", el("third grid viewer")));
    await session.step(69, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(70, "And no errors should have been logged", () => noErrors(page));
  });
});
