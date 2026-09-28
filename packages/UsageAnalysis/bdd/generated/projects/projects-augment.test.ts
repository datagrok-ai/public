/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-augment.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/minimized-viewers.js';
import '../../bindings/projects-copies.js';
import '../../bindings/projects-derived.js';
import '../../bindings/projects-regressions.js';
import '../../bindings/projects-sources.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, collapse, doubleClickOn, dragTo, enterInto, expand, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, contextPanelOpen, contextPanelShows, currentViewType, dialogCloses, noProjectOnServer, projectsOnServer, savedWithDataSync, switchTableView, tableReloaded, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dashboardsPanelOpen, tableViewsOpen} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A table dragged onto an open project in the Dashboards panel joins it", () => {
  const session = feature(test, "features/projects/projects-augment.feature", import.meta.url);
  test("A table dragged onto an open project in the Dashboards panel joins it", {tag: ["@journey", "@serial", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And no project named \"BDDAugment{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDAugment{time}")));
    await run.scenario("demog is saved as a one-table project", async () => {
      await session.step(32, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(33, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(34, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(35, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(36, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(37, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(38, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(39, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(40, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(41, "When user enters \"BDDAugment{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDAugment{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(42, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(43, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(44, "And an info balloon containing \"BDDAugment{time}\\\" uploaded\" should have been shown", () => infoBalloonText(page, session.text("BDDAugment{time}\" uploaded")));
      await session.step(45, "And 1 project named \"BDDAugment{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDAugment{time}")));
      await session.step(46, "And \"Share BDDAugment{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDAugment{time}\" dialog")), "visible"));
      await session.step(47, "When user clicks on CANCEL button in \"Share BDDAugment{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDAugment{time}\" dialog"))));
      await session.step(48, "Then the \"Share BDDAugment{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDAugment{time}")));
    });
    await run.scenario("iris dragged onto the project's node moves into the project", async () => {
      await session.step(51, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(52, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(53, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(54, "When user double-clicks Files---Demo---iris.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---iris.csv tree node inside browse tree")));
      await session.step(55, "Then the \"iris\" view should be current", () => viewIsCurrent(page, "iris"));
      await session.step(56, "And the table should have 150 rows", () => rowCount(page, 150));
      await session.step(57, "Given the dashboards panel of the left sidebar is open", () => dashboardsPanelOpen(page));
      await session.step(58, "Then New-Dashboard---iris tree node inside browse tree should be visible", () => shouldBe(page, el("New-Dashboard---iris tree node inside browse tree"), "visible"));
      await session.step(59, "When user expands BDDAugment{time} tree node inside browse tree", () => expand(page, el(session.text("BDDAugment{time} tree node inside browse tree"))));
      await session.step(60, "Then BDDAugment{time}---demog tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDAugment{time}---demog tree node inside browse tree")), "visible"));
      await session.step(61, "When user collapses BDDAugment{time} tree node inside browse tree", () => collapse(page, el(session.text("BDDAugment{time} tree node inside browse tree"))));
      await session.step(62, "And user drags New-Dashboard---iris tree node inside browse tree to BDDAugment{time} tree node inside browse tree", () => dragTo(page, el("New-Dashboard---iris tree node inside browse tree"), el(session.text("BDDAugment{time} tree node inside browse tree"))));
      await session.step(63, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(64, "And Move entity dialog should contain text \"BDDAugment{time} project\"", () => shouldContainText(page, el("Move entity dialog"), session.text("BDDAugment{time} project")));
      await session.step(65, "And \"iris\" project table in Move entity dialog should be visible", () => shouldBe(page, el("\"iris\" project table in Move entity dialog"), "visible"));
      await session.step(66, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(67, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(68, "When user expands BDDAugment{time} tree node inside browse tree", () => expand(page, el(session.text("BDDAugment{time} tree node inside browse tree"))));
      await session.step(69, "Then BDDAugment{time}---iris tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDAugment{time}---iris tree node inside browse tree")), "visible"));
      await session.step(70, "And BDDAugment{time}---demog tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDAugment{time}---demog tree node inside browse tree")), "visible"));
      await session.step(71, "And New-Dashboard---iris tree node inside browse tree should be absent", () => shouldBe(page, el("New-Dashboard---iris tree node inside browse tree"), "absent"));
    });
    await run.scenario("The node's SAVE saves the project with both tables", async () => {
      await session.step(74, "When user clicks on Save button in BDDAugment{time} tree node inside browse tree", () => clickOn(page, el(session.text("Save button in BDDAugment{time} tree node inside browse tree"))));
      await session.step(75, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(76, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(78, "And Data sync switch in \"iris\" project table in \"Save project\" dialog should be unchecked", () => shouldBe(page, el("Data sync switch in \"iris\" project table in \"Save project\" dialog"), "unchecked"));
      await session.step(79, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(80, "When user checks Data sync switch in \"iris\" project table in \"Save project\" dialog", () => check(page, el("Data sync switch in \"iris\" project table in \"Save project\" dialog")));
      await session.step(81, "Then \"Creation script\" button in \"iris\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"iris\" project table in \"Save project\" dialog"), "visible"));
      await session.step(82, "And \"Save original project\" radio choice in \"Save project\" dialog should be checked", () => shouldBe(page, el("\"Save original project\" radio choice in \"Save project\" dialog"), "checked"));
      await session.step(83, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(84, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(85, "And the \"demog\" table of the \"BDDAugment{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDDAugment{time}")));
      await session.step(86, "And the \"iris\" table of the \"BDDAugment{time}\" project should be saved with data sync", () => savedWithDataSync(page, "iris", session.text("BDDAugment{time}")));
    });
    await run.scenario("Reopened from Dashboards, the project holds and opens both tables", async () => {
      await session.step(89, "When user closes all views", () => closeAllViews(page));
      await session.step(90, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(91, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(92, "And user enters \"BDDAugment{time}\" into gallery search", () => enterInto(page, session.text("BDDAugment{time}"), el("gallery search")));
      await session.step(93, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(94, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(95, "When user clicks on BDDAugment{time} gallery card", () => clickOn(page, el(session.text("BDDAugment{time} gallery card"))));
      await session.step(96, "Then the context panel should show \"BDDAugment{time}\"", () => contextPanelShows(page, session.text("BDDAugment{time}")));
      await session.step(97, "When user expands \"Content\" pane in context panel", () => expand(page, el("\"Content\" pane in context panel")));
      await session.step(98, "Then BDDAugment{time}---demog tree node in context panel should be visible", () => shouldBe(page, el(session.text("BDDAugment{time}---demog tree node in context panel")), "visible"));
      await session.step(99, "And BDDAugment{time}---iris tree node in context panel should be visible", () => shouldBe(page, el(session.text("BDDAugment{time}---iris tree node in context panel")), "visible"));
      await session.step(100, "When user double-clicks on BDDAugment{time} gallery card", () => doubleClickOn(page, el(session.text("BDDAugment{time} gallery card"))));
      await session.step(101, "Then the table views \"demog, iris\" should be open", () => tableViewsOpen(page, "demog, iris"));
      await session.step(102, "And table \"demog\" should have been reloaded by data sync with 5850 rows", () => tableReloaded(page, "demog", 5850));
      await session.step(103, "And table \"iris\" should have been reloaded by data sync with 150 rows", () => tableReloaded(page, "iris", 150));
      await session.step(104, "When user switches to the \"iris\" table view", () => switchTableView(page, "iris"));
      await session.step(105, "Then status bar should contain text \"Rows: 150\"", () => shouldContainText(page, el("status bar"), "Rows: 150"));
      await session.step(106, "And status bar should contain text \"Columns: 6\"", () => shouldContainText(page, el("status bar"), "Columns: 6"));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
      await session.step(108, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Delete Project removes the project", async () => {
      await session.step(111, "When user closes all views", () => closeAllViews(page));
      await session.step(112, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(113, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(114, "And user enters \"BDDAugment{time}\" into gallery search", () => enterInto(page, session.text("BDDAugment{time}"), el("gallery search")));
      await session.step(115, "And user picks \"Delete Project\" from the context menu of BDDAugment{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDAugment{time} gallery card"))));
      await session.step(116, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(117, "And \"Are you sure?\" dialog should contain text \"Delete project \\\"BDDAugment{time}\\\"?\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDAugment{time}\"?")));
      await session.step(118, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(120, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(121, "And 0 projects named \"BDDAugment{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDAugment{time}")));
      await session.step(122, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(123, "Then BDDAugment{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDAugment{time} gallery card")), "absent"));
      await session.step(124, "And no errors should have been logged", () => noErrors(page));
      await session.step(125, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
