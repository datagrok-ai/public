/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/07-molecules/05-activity-cliffs.feature
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
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, toggle, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff, toolboxPaneHidden, tooltipHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverFirstFreeArea, noErrors, readingIs, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragTopBorder, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Molecules", () => {
  const session = feature(test, "features/howto/07-molecules/05-activity-cliffs.feature", import.meta.url);
  test("Find activity cliffs", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And user opens activity_cliffs dataset", () => openDataset(page, ds("activity_cliffs")));
    await session.step(19, "And the grid shows only the columns \"chembl_tid, smiles, Activity\"", () => gridShowsColumns(page, "chembl_tid, smiles, Activity"));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(24, "When user picks \"Chem > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Activity Cliffs..."), undefined, "Chem > Analyze > Activity Cliffs");
    await session.step(25, "Then \"Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Activity Cliffs\" dialog"), "visible"));
    await session.step(27, "When user clicks on OK button in \"Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Activity Cliffs\" dialog")), undefined, "Check the molecule and activity columns, then OK");
    await session.step(28, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(29, "And the \"cliffs\" reading of scatter plot viewer should be 20", () => readingIs(page, "cliffs", el("scatter plot viewer"), 20));
    await session.step(30, "When the video says \"Lines join similar molecules whose activity differs a lot\"", () => videoSays(page, "Lines join similar molecules whose activity differs a lot"));
    await session.step(32, "And user toggles \"Show only cliffs\" input in scatter plot viewer", () => toggle(page, el("\"Show only cliffs\" input in scatter plot viewer")), undefined, "Show only cliffs hides the rest");
    await session.step(33, "Then the \"only cliffs\" reading of scatter plot viewer should be \"true\"", () => readingReads(page, "only cliffs", el("scatter plot viewer"), "true"));
    await session.step(35, "When user hovers over the first free \"line\" area of scatter plot viewer", () => hoverFirstFreeArea(page, "line", el("scatter plot viewer")), undefined, "Point at a line to see the pair");
    await session.step(36, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(37, "When the video says \"The tooltip shows both molecules and their activities\" at the top", () => videoSaysAt(page, "The tooltip shows both molecules and their activities", "top"));
    await session.step(38, "Given the tooltip is hidden", () => tooltipHidden(page));
    await session.step(40, "When user clicks on \"20 cliffs\" button in scatter plot viewer", () => clickOn(page, el("\"20 cliffs\" button in scatter plot viewer")), undefined, "The cliffs button lists every pair");
    await session.step(41, "Then the open tableview should have 2 grid viewers", () => viewerCount(page, 2, "grid"));
    await session.step(43, "When user drags the top border of second grid viewer by 250 pixels up", () => dragTopBorder(page, el("second grid viewer"), 250), undefined, "Drag its top border up to see more pairs");
    await session.step(44, "And the video says \"Each row is one pair: both molecules and how far their activities differ\"", () => videoSays(page, "Each row is one pair: both molecules and how far their activities differ"));
    await session.step(45, "Then the table should have 270 rows", () => rowCount(page, 270));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
});
