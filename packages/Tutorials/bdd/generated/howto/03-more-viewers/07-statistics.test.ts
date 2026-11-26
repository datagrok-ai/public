/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/07-statistics.feature
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
import {clickOn, hoverOver, pressKey, shouldBe, typeInto, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, hoverArea, noErrors, pickFromAreaContextMenu, readingIs, readingReads, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dockToViewEdge, dockedAtViewEdge, readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/07-statistics.feature", import.meta.url);
  test("Summarize every column with the Statistics viewer", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(19, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(20, "When the video says \"Statistics: a summary of every column at once\"", () => videoSays(page, "Statistics: a summary of every column at once"));
    await session.step(21, "And user clicks on statistics icon on toolbox", () => clickOn(page, el("statistics icon on toolbox")));
    await session.step(22, "Then statistics viewer should be added to the open tableview", () => viewerAdded(page, "statistics"));
    await session.step(24, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "Hide the toolbox for more room");
    await session.step(25, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(27, "When user docks statistics viewer to the bottom edge of the view", () => dockToViewEdge(page, el("statistics viewer"), "bottom"), undefined, "Drag its title bar to the bottom, so every statistic fits");
    await session.step(28, "Then statistics viewer should be docked along the bottom edge of the view", () => dockedAtViewEdge(page, el("statistics viewer"), "bottom"));
    await session.step(29, "When the video says \"HEIGHT is missing in 751 rows: see nulls\" at the top", () => videoSaysAt(page, "HEIGHT is missing in 751 rows: see nulls", "top"));
    await session.step(30, "And user hovers over the \"cell nulls of HEIGHT\" area of statistics viewer", () => hoverArea(page, "cell nulls of HEIGHT", el("statistics viewer")));
    await session.step(31, "Then the \"nulls of HEIGHT\" reading of statistics viewer should be \"751\"", () => readingReads(page, "nulls of HEIGHT", el("statistics viewer"), "751"));
    await session.step(33, "When user picks \"Statistics > sum\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Statistics > sum", "row AGE", el("statistics viewer")), undefined, "Right-click a row to add a statistic");
    await session.step(34, "Then the \"stats\" reading of statistics viewer should contain \"sum\"", () => readingContains(page, "stats", el("statistics viewer"), "sum"));
    await session.step(36, "When user picks \"Histograms > SEX\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Histograms > SEX", "row AGE", el("statistics viewer")), undefined, "or a histogram of each column, split by a category");
    await session.step(37, "Then the \"histogram columns\" reading of statistics viewer should be 1", () => readingIs(page, "histogram columns", el("statistics viewer"), 1));
    await session.step(38, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(40, "And user hovers over statistics viewer", () => hoverOver(page, el("statistics viewer")), undefined, "The settings hold the rest");
    await session.step(41, "And user clicks on settings icon of statistics viewer", () => clickOn(page, el("settings icon of statistics viewer")));
    await session.step(42, "Then the context panel should show \"Statistics\"", () => contextPanelShows(page, "Statistics"));
    await session.step(44, "When user types \"${AGE} > 60\" into \"Filter\" property in context panel", () => typeInto(page, "${AGE} > 60", el("\"Filter\" property in context panel")), undefined, "Data: a filter formula summarizes only some rows");
    await session.step(45, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(46, "Then the \"rows shown\" reading of statistics viewer should be 878", () => readingIs(page, "rows shown", el("statistics viewer"), 878));
    await session.step(47, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
});
