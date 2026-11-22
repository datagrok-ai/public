/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/06-dashboards-and-sharing/01-save-and-find-a-dashboard.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, browsePanelOpen, closeAllViews, contextPanelClosed, dialogCloses, doubleClickOpensProject, noProjectOnServer, openDataset, projectsOnServer, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Dashboards", () => {
  const session = feature(test, "features/howto/06-dashboards-and-sharing/01-save-and-find-a-dashboard.feature", import.meta.url);
  test("Save a dashboard and find it again", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And no project named \"Demographics overview\" is on the server", () => noProjectOnServer(page, "Demographics overview"));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]);
    await session.step(21, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","HEIGHT"]]), [["valueColumnName","HEIGHT"]]);
    await session.step(23, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"]]), [["splitColumnName","RACE"]]);
    await session.step(25, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(26, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(27, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(28, "When the video says \"A dashboard is a table and its viewers, saved together\"", () => videoSays(page, "A dashboard is a table and its viewers, saved together"));
    await session.step(30, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")), undefined, "SAVE puts the table and its viewers on the server");
    await session.step(31, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(33, "When user enters \"Demographics overview\" into Name text input in \"Save project\" dialog", () => enterInto(page, "Demographics overview", el("Name text input in \"Save project\" dialog")), undefined, "Give it a name you will recognize");
    await session.step(34, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(35, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(36, "And 1 project named \"Demographics overview\" should be on the server", () => projectsOnServer(page, 1, "Demographics overview"));
    await session.step(37, "And \"Share Demographics overview\" dialog should be visible", () => shouldBe(page, el("\"Share Demographics overview\" dialog"), "visible"));
    await session.step(39, "When user clicks on CANCEL button in \"Share Demographics overview\" dialog", () => clickOn(page, el("CANCEL button in \"Share Demographics overview\" dialog")), undefined, "Next it offers to share. Until you do, only you can see it");
    await session.step(40, "Then the \"Share Demographics overview\" dialog should close", () => dialogCloses(page, "Share Demographics overview"));
    await session.step(41, "Given user closes all views", () => closeAllViews(page));
    await session.step(42, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "When the video says \"Later, find it again in Browse\"", () => videoSays(page, "Later, find it again in Browse"));
    await session.step(45, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")), undefined, "Dashboards lists every dashboard you can open");
    await session.step(46, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
    await session.step(48, "When user enters \"Demographics overview\" into gallery search", () => enterInto(page, "Demographics overview", el("gallery search")), undefined, "Search by name");
    await session.step(49, "Then \"Demographics overview\" gallery card should be visible", () => shouldBe(page, el("\"Demographics overview\" gallery card"), "visible"));
    await session.step(51, "When user double-clicks on \"Demographics overview\" gallery card and the project opens", () => doubleClickOpensProject(page, el("\"Demographics overview\" gallery card")), undefined, "Double-click to open it");
    await session.step(52, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(53, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(54, "And histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
    await session.step(55, "And bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
    await session.step(57, "When user expands \"My stuff\" tree node inside browse tree", () => expand(page, el("\"My stuff\" tree node inside browse tree")), undefined, "Recently opened dashboards are also under My stuff > Recent");
    await session.step(58, "And user expands \"My stuff > Recent\" tree node inside browse tree", () => expand(page, el("\"My stuff > Recent\" tree node inside browse tree")));
    await session.step(59, "Then \"My stuff > Recent > Demographics overview\" tree node inside browse tree should be visible", () => shouldBe(page, el("\"My stuff > Recent > Demographics overview\" tree node inside browse tree"), "visible"));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
});
