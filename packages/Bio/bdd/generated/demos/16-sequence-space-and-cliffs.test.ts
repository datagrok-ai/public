/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/16-sequence-space-and-cliffs.feature
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
import {dockedUnderGrid} from '../../bindings/demo-layout.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldHaveText, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickFirstArea, hoverArea, noBalloons, noErrors, painted, propertyShouldContain, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("4 · A peptide series · Sequence space and activity cliffs", () => {
  const session = feature(test, "features/demos/16-sequence-space-and-cliffs.feature", import.meta.url);
  test("Map the sequence space and find the activity cliffs", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And simple mode is off", () => simpleModeOff(page));
    await session.step(27, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(28, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(29, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(30, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(31, "And the grid shows only the columns \"sequence, activity\"", () => gridShowsColumns(page, "sequence, activity"));
    await session.step(32, "And the \"sequence\" column of the grid is 260 px wide", () => gridColumnWidth(page, "sequence", 260));
    await session.step(33, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
    await session.step(34, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
    await session.step(35, "When user selects \"Monomer chemical distance\" in Similarity input in \"Sequence Space\" dialog", () => selectIn(page, "Monomer chemical distance", el("Similarity input in \"Sequence Space\" dialog")));
    await session.step(36, "Then Similarity input in \"Sequence Space\" dialog should have value \"Monomer chemical distance\"", () => shouldHaveValue(page, el("Similarity input in \"Sequence Space\" dialog"), "Monomer chemical distance"));
    await session.step(37, "When user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
    await session.step(38, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(39, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
    await session.step(40, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(41, "And \"Description\" property of scatter plot viewer should contain \"similarity: Monomer chemical distance\"", () => propertyShouldContain(page, "Description", el("scatter plot viewer"), "similarity: Monomer chemical distance"));
    await session.step(42, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await session.step(43, "Given the grid shows only the columns \"sequence, activity\"", () => gridShowsColumns(page, "sequence, activity"));
    await session.step(44, "When the video says \"Placed by monomer chemistry: alike peptides land together\"", () => videoSays(page, "Placed by monomer chemistry: alike peptides land together"));
    await session.step(45, "And user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
    await session.step(46, "Then \"Sequence Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "visible"));
    await session.step(47, "And editor of Activities input in \"Sequence Activity Cliffs\" dialog should have text \"activity\"", () => shouldHaveText(page, el("editor of Activities input in \"Sequence Activity Cliffs\" dialog"), "activity"));
    await session.step(48, "When user selects \"Monomer chemical distance\" in Similarity input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Monomer chemical distance", el("Similarity input in \"Sequence Activity Cliffs\" dialog")));
    await session.step(49, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
    await session.step(50, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(51, "And second scatter plot viewer should be visible", () => shouldBe(page, el("second scatter plot viewer"), "visible"));
    await session.step(52, "And the \"cliffs\" reading of second scatter plot viewer should be 21", () => readingIs(page, "cliffs", el("second scatter plot viewer"), 21));
    await session.step(53, "And second scatter plot viewer should be painted", () => painted(page, el("second scatter plot viewer")));
    await session.step(54, "Given the grid shows only the columns \"sequence, activity\"", () => gridShowsColumns(page, "sequence, activity"));
    await session.step(55, "And user sets \"zoomAndFilter\" property of second scatter plot viewer to \"no action\"", () => setProperty(page, "zoomAndFilter", el("second scatter plot viewer"), "no action"));
    await session.step(56, "When user clicks on \"cliffs\" button", () => clickOn(page, el("\"cliffs\" button")));
    await session.step(57, "Then second grid viewer should be visible", () => shouldBe(page, el("second grid viewer"), "visible"));
    await session.step(58, "Given second grid viewer is docked under the grid, 35% of the view tall", () => dockedUnderGrid(page, el("second grid viewer"), 35));
    await session.step(59, "When user clicks on the first \"cell\" area of second grid viewer", () => clickFirstArea(page, "cell", el("second grid viewer")));
    await session.step(60, "Then row 70 should be current", () => currentRowIs(page, 70));
    await session.step(61, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(62, "And second scatter plot viewer should be painted", () => painted(page, el("second scatter plot viewer")));
    await session.step(63, "When the video says \"A cliff: one residue apart, far apart in activity\"", () => videoSays(page, "A cliff: one residue apart, far apart in activity"));
    await session.step(64, "And user hovers over the \"cell 71 of activity\" area of grid", () => hoverArea(page, "cell 71 of activity", el("grid")));
    await session.step(65, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
});
