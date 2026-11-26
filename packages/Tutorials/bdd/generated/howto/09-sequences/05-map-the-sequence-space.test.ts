/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/09-sequences/05-map-the-sequence-space.feature
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
import {bioReady} from '../../../bindings/howto-sequences.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldHaveText, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickFirstArea, hoverArea, noBalloons, noErrors, painted, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Sequences", () => {
  const session = feature(test, "features/howto/09-sequences/05-map-the-sequence-space.feature", import.meta.url);
  test("Map the sequence space and find the activity cliffs", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And simple mode is off", () => simpleModeOff(page));
    await session.step(22, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(23, "And the Bio package is initialized", () => bioReady(page));
    await session.step(24, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(25, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(26, "And the grid shows only the columns \"sequence, activity\"", () => gridShowsColumns(page, "sequence, activity"));
    await session.step(27, "And the \"sequence\" column of the grid is 260 px wide", () => gridColumnWidth(page, "sequence", 260));
    await session.step(29, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."), undefined, "Bio > Analyze > Sequence Space");
    await session.step(30, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
    await session.step(32, "When user selects \"Monomer chemical distance\" in Similarity input in \"Sequence Space\" dialog", () => selectIn(page, "Monomer chemical distance", el("Similarity input in \"Sequence Space\" dialog")), undefined, "Similarity: how two sequences are compared");
    await session.step(33, "Then Similarity input in \"Sequence Space\" dialog should have value \"Monomer chemical distance\"", () => shouldHaveValue(page, el("Similarity input in \"Sequence Space\" dialog"), "Monomer chemical distance"));
    await session.step(34, "When user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
    await session.step(35, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(36, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
    await session.step(37, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(38, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await session.step(39, "When the video says \"Each dot is a peptide; alike peptides land together, colored by cluster\"", () => videoSays(page, "Each dot is a peptide; alike peptides land together, colored by cluster"));
    await session.step(41, "And user hovers over the \"marker of row 5\" area of scatter plot viewer", () => hoverArea(page, "marker of row 5", el("scatter plot viewer")), undefined, "The pointer on a dot shows its peptide");
    await session.step(43, "And user clicks on the \"marker of row 15\" area of scatter plot viewer", () => clickArea(page, "marker of row 15", el("scatter plot viewer")), undefined, "Click a dot to find its peptide in the table");
    await session.step(44, "Then row 15 should be current", () => currentRowIs(page, 15));
    await session.step(45, "When the video says \"Activity cliffs: alike peptides whose activity differs a lot\"", () => videoSays(page, "Activity cliffs: alike peptides whose activity differs a lot"));
    await session.step(47, "And user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."), undefined, "Bio > Analyze > Activity Cliffs");
    await session.step(48, "Then \"Sequence Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "visible"));
    await session.step(49, "And editor of Activities input in \"Sequence Activity Cliffs\" dialog should have text \"activity\"", () => shouldHaveText(page, el("editor of Activities input in \"Sequence Activity Cliffs\" dialog"), "activity"));
    await session.step(51, "When user selects \"Monomer chemical distance\" in Similarity input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Monomer chemical distance", el("Similarity input in \"Sequence Activity Cliffs\" dialog")), undefined, "The sequences, their activity, and how to compare them");
    await session.step(52, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
    await session.step(53, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(54, "And second scatter plot viewer should be visible", () => shouldBe(page, el("second scatter plot viewer"), "visible"));
    await session.step(55, "And the \"cliffs\" reading of second scatter plot viewer should be 21", () => readingIs(page, "cliffs", el("second scatter plot viewer"), 21));
    await session.step(56, "And second scatter plot viewer should be painted", () => painted(page, el("second scatter plot viewer")));
    await session.step(57, "Given user sets \"zoomAndFilter\" property of second scatter plot viewer to \"no action\"", () => setProperty(page, "zoomAndFilter", el("second scatter plot viewer"), "no action"));
    await session.step(58, "When the video says \"Colored by activity; alike peptides linked, bigger dots on steeper cliffs\"", () => videoSays(page, "Colored by activity; alike peptides linked, bigger dots on steeper cliffs"));
    await session.step(60, "And user clicks on \"cliffs\" button", () => clickOn(page, el("\"cliffs\" button")), undefined, "This button lists the cliffs: 21 pairs");
    await session.step(61, "Then second grid viewer should be visible", () => shouldBe(page, el("second grid viewer"), "visible"));
    await session.step(63, "When user clicks on the first \"cell\" area of second grid viewer", () => clickFirstArea(page, "cell", el("second grid viewer")), undefined, "Click a pair: the plot zooms to it");
    await session.step(64, "Then row 70 should be current", () => currentRowIs(page, 70));
    await session.step(65, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(66, "When the video says \"A cliff: one residue apart, far apart in activity\"", () => videoSays(page, "A cliff: one residue apart, far apart in activity"));
    await session.step(67, "And user hovers over the \"cell 71 of activity\" area of grid", () => hoverArea(page, "cell 71 of activity", el("grid")));
    await session.step(68, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(69, "And no errors should have been logged", () => noErrors(page));
  });
});
