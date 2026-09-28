/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-data-sync.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, enterInto, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, currentViewType, dialogCloses, loadedAsSnapshot, noProjectOnServer, openProjectWithTable, projectsOnServer, reloadedByDataSync, savedAsSnapshot, savedWithDataSync} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project saved with and without Data sync", () => {
  const session = feature(test, "features/projects/projects-data-sync.feature", import.meta.url);
  test("A project saved with and without Data sync", {tag: ["@journey", "@serial", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And no project named \"BDD-Sync-Orig-{time}\" is on the server", () => noProjectOnServer(page, session.text("BDD-Sync-Orig-{time}")));
    await session.step(20, "And no project named \"BDD-Sync-Copy-{time}\" is on the server", () => noProjectOnServer(page, session.text("BDD-Sync-Copy-{time}")));
    await run.scenario("A file opened from Browse is saved with Data sync on", async () => {
      await session.step(23, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(24, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(25, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(26, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(27, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(28, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(29, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(30, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(31, "When user enters \"BDD-Sync-Orig-{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDD-Sync-Orig-{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(32, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(33, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(34, "And an info balloon containing \"BDD-Sync-Orig-{time}\" should have been shown", () => infoBalloonText(page, session.text("BDD-Sync-Orig-{time}")));
      await session.step(35, "And 1 project named \"BDD-Sync-Orig-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-Sync-Orig-{time}")));
      await session.step(36, "And the \"demog\" table of the \"BDD-Sync-Orig-{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDD-Sync-Orig-{time}")));
      await session.step(37, "And \"Share BDD-Sync-Orig-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-Sync-Orig-{time}\" dialog")), "visible"));
      await session.step(38, "When user clicks on CANCEL button in \"Share BDD-Sync-Orig-{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDD-Sync-Orig-{time}\" dialog"))));
      await session.step(39, "Then the \"Share BDD-Sync-Orig-{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-Sync-Orig-{time}")));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A copy saved with Data sync off leaves the original synced", async () => {
      await session.step(43, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(44, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(45, "When user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(46, "Then Name text input in \"Save project\" dialog should have value \"Copy of BDD-Sync-Orig-{time}\"", () => shouldHaveValue(page, el("Name text input in \"Save project\" dialog"), session.text("Copy of BDD-Sync-Orig-{time}")));
      await session.step(47, "And choice input in \"demog\" project table in \"Save project\" dialog should have value \"Clone\"", () => shouldHaveValue(page, el("choice input in \"demog\" project table in \"Save project\" dialog"), "Clone"));
      await session.step(48, "When user enters \"BDD-Sync-Copy-{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDD-Sync-Copy-{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(49, "And user unchecks Data sync switch in \"demog\" project table in \"Save project\" dialog", () => uncheck(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog")));
      await session.step(50, "Then \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be hidden", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "hidden"));
      await session.step(51, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(52, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(53, "And 1 project named \"BDD-Sync-Copy-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-Sync-Copy-{time}")));
      await session.step(54, "And the \"demog\" table of the \"BDD-Sync-Copy-{time}\" project should be saved as a snapshot", () => savedAsSnapshot(page, "demog", session.text("BDD-Sync-Copy-{time}")));
      await session.step(55, "And the \"demog\" table of the \"BDD-Sync-Orig-{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDD-Sync-Orig-{time}")));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The copy opens from its snapshot", async () => {
      await session.step(59, "When user closes all views", () => closeAllViews(page));
      await session.step(60, "And user opens the \"BDD-Sync-Copy-{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDD-Sync-Copy-{time}")));
      await session.step(61, "Then the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(62, "And the table should have been loaded as a snapshot", () => loadedAsSnapshot(page));
      await session.step(63, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(64, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(65, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be unchecked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "unchecked"));
      await session.step(66, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(67, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Turning Data sync on for the copy makes it reread the file", async () => {
      await session.step(71, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(72, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(73, "When user checks Data sync switch in \"demog\" project table in \"Save project\" dialog", () => check(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog")));
      await session.step(74, "And user clicks on \"Creation script\" button in \"demog\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog")));
      await session.step(75, "Then \"demog\" project table in \"Save project\" dialog should contain text \"System:DemoFiles/demog.csv\"", () => shouldContainText(page, el("\"demog\" project table in \"Save project\" dialog"), "System:DemoFiles/demog.csv"));
      await session.step(76, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(77, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(78, "And the \"demog\" table of the \"BDD-Sync-Copy-{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDD-Sync-Copy-{time}")));
      await session.step(79, "When user closes all views", () => closeAllViews(page));
      await session.step(80, "And user opens the \"BDD-Sync-Copy-{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDD-Sync-Copy-{time}")));
      await session.step(81, "Then the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(82, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The original still opens by rereading the file", async () => {
      await session.step(86, "When user closes all views", () => closeAllViews(page));
      await session.step(87, "And user opens the \"BDD-Sync-Orig-{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDD-Sync-Orig-{time}")));
      await session.step(88, "Then the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(89, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
      await session.step(91, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
