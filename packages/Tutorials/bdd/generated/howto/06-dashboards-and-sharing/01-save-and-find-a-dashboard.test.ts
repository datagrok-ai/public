/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/06-dashboards-and-sharing/01-save-and-find-a-dashboard.feature
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
import {clickOn, close, enterInto, expand, pressKey, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {favoriteOnServer, notFavorite} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, dialogCloses, doubleClickOpensProject, noProjectOnServer, openDataset, projectsOnServer, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableLeft} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {addViewerWith, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Dashboards", () => {
  const session = feature(test, "features/howto/06-dashboards-and-sharing/01-save-and-find-a-dashboard.feature", import.meta.url);
  test("Save a dashboard and find it again", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And no project named \"Demographics overview\" is on the server", () => noProjectOnServer(page, "Demographics overview"));
    await session.step(22, "And \"Demographics overview\" is not in favorites", () => notFavorite(page, "Demographics overview"));
    await session.step(23, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(24, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]);
    await session.step(28, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","HEIGHT"]]), [["valueColumnName","HEIGHT"]]);
    await session.step(30, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"]]), [["splitColumnName","RACE"]]);
    await session.step(32, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(33, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(34, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(35, "When the video says \"A dashboard is a table and its viewers, saved together\"", () => videoSays(page, "A dashboard is a table and its viewers, saved together"));
    await session.step(37, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")), undefined, "SAVE puts the table and its viewers on the server");
    await session.step(38, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(40, "When user enters \"Demographics overview\" into Name text input in \"Save project\" dialog", () => enterInto(page, "Demographics overview", el("Name text input in \"Save project\" dialog")), undefined, "Give it a name you will recognize");
    await session.step(41, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(42, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(43, "And 1 project named \"Demographics overview\" should be on the server", () => projectsOnServer(page, 1, "Demographics overview"));
    await session.step(44, "And \"Share Demographics overview\" dialog should be visible", () => shouldBe(page, el("\"Share Demographics overview\" dialog"), "visible"));
    await session.step(46, "When user clicks on CANCEL button in \"Share Demographics overview\" dialog", () => clickOn(page, el("CANCEL button in \"Share Demographics overview\" dialog")), undefined, "Next it offers to share. Until you do, only you can see it");
    await session.step(47, "Then the \"Share Demographics overview\" dialog should close", () => dialogCloses(page, "Share Demographics overview"));
    await session.step(49, "When user closes \"demog\" view tab", () => close(page, el("\"demog\" view tab")), undefined, "Close it for now. It is kept on the server");
    await session.step(50, "Then no table should be left in the workspace", () => noTableLeft(page));
    await session.step(51, "And browse panel should be visible", () => shouldBe(page, el("browse panel"), "visible"));
    await session.step(52, "When the video says \"Later, find it again in Browse\"", () => videoSays(page, "Later, find it again in Browse"));
    await session.step(54, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")), undefined, "Dashboards lists every dashboard you can open. Search by name");
    await session.step(55, "And user enters \"Demographics overview\" into gallery search", () => enterInto(page, "Demographics overview", el("gallery search")));
    await session.step(56, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
    await session.step(57, "And \"Demographics overview\" gallery card should be visible", () => shouldBe(page, el("\"Demographics overview\" gallery card"), "visible"));
    await session.step(59, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "The Context Panel shows details of what you click");
    await session.step(60, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(62, "When user clicks on \"Demographics overview\" gallery card", () => clickOn(page, el("\"Demographics overview\" gallery card")), undefined, "Click the card: who made it and when");
    await session.step(63, "Then the context panel should show \"Demographics overview\"", () => contextPanelShows(page, "Demographics overview"));
    await session.step(65, "When user expands \"My stuff\" tree node inside browse tree", () => expand(page, el("\"My stuff\" tree node inside browse tree")), undefined, "What you used lately is also under My stuff > Recent");
    await session.step(66, "And user expands \"My stuff > Recent\" tree node inside browse tree", () => expand(page, el("\"My stuff > Recent\" tree node inside browse tree")));
    await session.step(67, "Then \"My stuff > Recent > Demographics overview\" tree node inside browse tree should be visible", () => shouldBe(page, el("\"My stuff > Recent > Demographics overview\" tree node inside browse tree"), "visible"));
    await session.step(69, "When user picks \"Add To Favorites > Only for me\" from the context menu of \"Demographics overview\" gallery card", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el("\"Demographics overview\" gallery card")), undefined, "Keep the ones you open often in Favorites: only for you, or for a group");
    await session.step(70, "Then \"Demographics overview\" should be in favorites on the server", () => favoriteOnServer(page, "Demographics overview"));
    await session.step(72, "When user presses Escape", () => pressKey(page, "Escape"), undefined, "The menu stays open for more groups. Escape closes it");
    await session.step(74, "When user clicks on favorites sidebar tab", () => clickOn(page, el("favorites sidebar tab")), undefined, "The star on the left opens your favorites");
    await session.step(75, "Then favorites panel should be visible", () => shouldBe(page, el("favorites panel"), "visible"));
    await session.step(76, "And \"Demographics overview\" tree node in favorites panel should be visible", () => shouldBe(page, el("\"Demographics overview\" tree node in favorites panel"), "visible"));
    await session.step(78, "When user double-clicks on \"Demographics overview\" tree node in favorites panel and the project opens", () => doubleClickOpensProject(page, el("\"Demographics overview\" tree node in favorites panel")), undefined, "Double-click to open it, with all its viewers");
    await session.step(79, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(80, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(81, "And histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
    await session.step(82, "And bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
    await session.step(83, "And no errors should have been logged", () => noErrors(page));
  });
});
