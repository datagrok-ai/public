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
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldHaveText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickFirstArea, hoverArea, noBalloons, noErrors, oneTooltip, painted, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Sequence space and activity cliffs", () => {
  const session = feature(test, "features/demos/16-sequence-space-and-cliffs.feature", import.meta.url);
  test("Map the sequence space and find the activity cliffs", {tag: ["@guide", "@demo", "@viewport:1920x1080"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And simple mode is off", () => simpleModeOff(page));
    await session.step(23, "And user opens FASTA_PT_activity dataset", () => openDataset(page, ds("FASTA_PT_activity")));
    await session.step(24, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(25, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(26, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(27, "And the grid shows only the columns \"sequence, activity\"", () => gridShowsColumns(page, "sequence, activity"));
    await session.step(28, "And the \"sequence\" column of the grid is 260 px wide", () => gridColumnWidth(page, "sequence", 260));
    await session.step(29, "When the video says \"Sequence space: alike peptides land close together\"", () => videoSays(page, "Sequence space: alike peptides land close together"));
    await session.step(30, "And user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
    await session.step(31, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
    await session.step(32, "When user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
    await session.step(33, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(34, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
    await session.step(35, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(36, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await session.step(37, "Given the grid shows only the columns \"sequence, activity\"", () => gridShowsColumns(page, "sequence, activity"));
    await session.step(38, "When user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
    await session.step(39, "Then \"Sequence Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "visible"));
    await session.step(40, "And editor of Activities input in \"Sequence Activity Cliffs\" dialog should have text \"activity\"", () => shouldHaveText(page, el("editor of Activities input in \"Sequence Activity Cliffs\" dialog"), "activity"));
    await session.step(41, "When user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
    await session.step(42, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(43, "And second scatter plot viewer should be visible", () => shouldBe(page, el("second scatter plot viewer"), "visible"));
    await session.step(44, "And the \"cliffs\" reading of second scatter plot viewer should be 3", () => readingIs(page, "cliffs", el("second scatter plot viewer"), 3));
    await session.step(45, "And second scatter plot viewer should be painted", () => painted(page, el("second scatter plot viewer")));
    await session.step(46, "Given the grid shows only the columns \"sequence, activity\"", () => gridShowsColumns(page, "sequence, activity"));
    await session.step(47, "And user sets \"zoomAndFilter\" property of second scatter plot viewer to \"no action\"", () => setProperty(page, "zoomAndFilter", el("second scatter plot viewer"), "no action"));
    await session.step(48, "When user clicks on \"cliffs\" button", () => clickOn(page, el("\"cliffs\" button")));
    await session.step(49, "Then second grid viewer should be visible", () => shouldBe(page, el("second grid viewer"), "visible"));
    await session.step(50, "When user clicks on the first \"cell\" area of second grid viewer", () => clickFirstArea(page, "cell", el("second grid viewer")));
    await session.step(51, "Then row 70 should be current", () => currentRowIs(page, 70));
    await session.step(52, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(53, "And second scatter plot viewer should be painted", () => painted(page, el("second scatter plot viewer")));
    await session.step(54, "When the video says \"A cliff: one residue apart, far apart in activity\"", () => videoSays(page, "A cliff: one residue apart, far apart in activity"));
    await session.step(55, "And user hovers over the \"view\" area of second scatter plot viewer", () => hoverArea(page, "view", el("second scatter plot viewer")));
    await session.step(56, "Then exactly one tooltip should be shown", () => oneTooltip(page));
    await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(58, "And no errors should have been logged", () => noErrors(page));
  });
});
