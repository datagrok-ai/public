/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-dock-layout-app.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, core#239]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openAppOnView} from '../../bindings/demo-oligo.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, pressKey, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, clickInMyStuffBucket, clickViewTab, closeAllViews, dialogCloses, dockViewBeside, doubleClickOpensProject, expandMyStuffBucket, includeSaveDialogView, openDataset, openSaveDialog, ownProjectGone, projectsOnServer, saveDialogViewExcluded, saveDialogViewIncluded, saveFromDialog, simpleModeOff, switchView, viewDockedBeside, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableLeft} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A dashboard saves an app view with its state and starts it only when opened for real", () => {
  const session = feature(test, "features/projects/projects-dock-layout-app.feature", import.meta.url);
  test("A dashboard saves an app view with its state and starts it only when opened for real", {tag: ["@journey", "@realizes:views.projects", "@realizes:core#239"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the user's own project \"BDDLayoutApp{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("BDDLayoutApp{time}")));
    await run.scenario("Usage Analysis is set to a sub-view and a filter, and docked beside a table", async () => {
      await session.step(20, "Given user opens cars dataset", () => openDataset(page, ds("cars")));
      await session.step(21, "And user opens the \"Usage Analysis\" app on its \"Overview\" view", () => openAppOnView(page, "Usage Analysis", "Overview"));
      await session.step(22, "When user clicks on Functions tab", () => clickOn(page, el("Functions tab")));
      await session.step(23, "Then the \"Functions\" view should be current", () => viewIsCurrent(page, "Functions"));
      await session.step(24, "When user clicks on \"Execution time\" button", () => clickOn(page, el("\"Execution time\" button")));
      await session.step(25, "Then \"Execution time\" button should be disabled", () => shouldBe(page, el("\"Execution time\" button"), "disabled"));
      await session.step(26, "When user types \"this month\" into Date text input", () => typeInto(page, "this month", el("Date text input")));
      await session.step(27, "And user clicks on Apply button", () => clickOn(page, el("Apply button")));
      await session.step(28, "Then Apply button should be disabled", () => shouldBe(page, el("Apply button"), "disabled"));
      await session.step(29, "When user drags the tab of the \"Functions\" view to the right side of the \"cars\" view", () => dockViewBeside(page, "Functions", "right", "cars"));
      await session.step(30, "Then the \"Functions\" view should be docked right-of the \"cars\" view", () => viewDockedBeside(page, "Functions", "right-of", "cars"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The dashboard reopens Usage Analysis where it was, on its sub-view and filter", async () => {
      await session.step(34, "When user clicks on the tab of the \"cars\" view", () => clickViewTab(page, "cars"));
      await session.step(35, "And user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(36, "Then the \"Functions\" view should be excluded in the Save project dialog", () => saveDialogViewExcluded(page, "Functions"));
      await session.step(37, "When user includes the \"Functions\" view in the Save project dialog", () => includeSaveDialogView(page, "Functions"));
      await session.step(38, "Then the \"Functions\" view should be included in the Save project dialog", () => saveDialogViewIncluded(page, "Functions"));
      await session.step(39, "When user types \"BDDLayoutApp{time}\" into text input in \"Save project\" dialog", () => typeInto(page, session.text("BDDLayoutApp{time}"), el("text input in \"Save project\" dialog")));
      await session.step(40, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(41, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(42, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(43, "Then 1 project named \"BDDLayoutApp{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLayoutApp{time}")));
      await session.step(44, "When user closes all views", () => closeAllViews(page));
      await session.step(45, "Then no table should be left in the workspace", () => noTableLeft(page));
      await session.step(46, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(47, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(48, "And user enters \"BDDLayoutApp{time}\" into gallery search", () => enterInto(page, session.text("BDDLayoutApp{time}"), el("gallery search")));
      await session.step(49, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(50, "And user double-clicks on BDDLayoutApp{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDLayoutApp{time} gallery card"))));
      await session.step(51, "Then the \"Functions\" view should be docked right-of the \"cars\" view", () => viewDockedBeside(page, "Functions", "right-of", "cars"));
      await session.step(52, "Given user switches to the \"Functions\" view", () => switchView(page, "Functions"));
      await session.step(53, "Then \"Execution time\" button should be disabled", () => shouldBe(page, el("\"Execution time\" button"), "disabled"));
      await session.step(54, "And Date text input should have value \"this month\"", () => shouldHaveValue(page, el("Date text input"), "this month"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Browse preview of the dashboard does not start the app, and pinning it does", async () => {
      await session.step(58, "When user closes all views", () => closeAllViews(page));
      await session.step(59, "Then no table should be left in the workspace", () => noTableLeft(page));
      await session.step(60, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(61, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
      await session.step(62, "When user expands the \"Dashboards\" bucket of My stuff", () => expandMyStuffBucket(page, "Dashboards"));
      await session.step(63, "And user clicks on \"BDDLayoutApp{time}\" in the \"Dashboards\" bucket of My stuff", () => clickInMyStuffBucket(page, session.text("BDDLayoutApp{time}"), "Dashboards"));
      await session.step(64, "Then \"cars\" view should be visible", () => shouldBe(page, el("\"cars\" view"), "visible"));
      await session.step(65, "And \"Functions\" view should be absent", () => shouldBe(page, el("\"Functions\" view"), "absent"));
      await session.step(66, "When user double-clicks on cars view", () => doubleClickOn(page, el("cars view")));
      await session.step(67, "Then \"Functions\" view should be visible", () => shouldBe(page, el("\"Functions\" view"), "visible"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
