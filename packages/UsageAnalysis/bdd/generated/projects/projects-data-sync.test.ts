/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-data-sync.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects]
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, enterInto, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, currentViewType, dialogCloses, noProjectOnServer, projectsOnServer, savedAsSnapshot, savedWithDataSync} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project saved with and without Data sync", () => {
  const session = feature(test, "features/projects/projects-data-sync.feature", import.meta.url);
  test("A project saved with and without Data sync", {tag: ["@journey", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And no project named \"BDD-Sync-Orig-{time}\" is on the server", () => noProjectOnServer(page, session.text("BDD-Sync-Orig-{time}")));
    await session.step(22, "And no project named \"BDD-Sync-Copy-{time}\" is on the server", () => noProjectOnServer(page, session.text("BDD-Sync-Copy-{time}")));
    await run.scenario("A file opened from Browse is saved with Data sync on", async () => {
      await session.step(25, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(26, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(27, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(28, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(29, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(30, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(31, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(32, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(33, "When user enters \"BDD-Sync-Orig-{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDD-Sync-Orig-{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(34, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(35, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(36, "And an info balloon containing \"BDD-Sync-Orig-{time}\" should have been shown", () => infoBalloonText(page, session.text("BDD-Sync-Orig-{time}")));
      await session.step(37, "And 1 project named \"BDD-Sync-Orig-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-Sync-Orig-{time}")));
      await session.step(38, "And the \"demog\" table of the \"BDD-Sync-Orig-{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDD-Sync-Orig-{time}")));
      await session.step(39, "And \"Share BDD-Sync-Orig-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-Sync-Orig-{time}\" dialog")), "visible"));
      await session.step(40, "When user clicks on CANCEL button in \"Share BDD-Sync-Orig-{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDD-Sync-Orig-{time}\" dialog"))));
      await session.step(41, "Then the \"Share BDD-Sync-Orig-{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-Sync-Orig-{time}")));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A copy saved with Data sync off leaves the original synced", async () => {
      await session.step(45, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(46, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(47, "When user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(48, "Then Name text input in \"Save project\" dialog should have value \"Copy of BDD-Sync-Orig-{time}\"", () => shouldHaveValue(page, el("Name text input in \"Save project\" dialog"), session.text("Copy of BDD-Sync-Orig-{time}")));
      await session.step(49, "And choice input in \"demog\" project table in \"Save project\" dialog should have value \"Clone\"", () => shouldHaveValue(page, el("choice input in \"demog\" project table in \"Save project\" dialog"), "Clone"));
      await session.step(50, "When user enters \"BDD-Sync-Copy-{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDD-Sync-Copy-{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(51, "And user unchecks Data sync switch in \"demog\" project table in \"Save project\" dialog", () => uncheck(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog")));
      await session.step(52, "Then Data sync switch in \"demog\" project table in \"Save project\" dialog should be unchecked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "unchecked"));
      await session.step(53, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be hidden", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "hidden"));
      await session.step(54, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(55, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(56, "And 1 project named \"BDD-Sync-Copy-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-Sync-Copy-{time}")));
      await session.step(57, "And the \"demog\" table of the \"BDD-Sync-Copy-{time}\" project should be saved as a snapshot", () => savedAsSnapshot(page, "demog", session.text("BDD-Sync-Copy-{time}")));
      await session.step(58, "And the \"demog\" table of the \"BDD-Sync-Orig-{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDD-Sync-Orig-{time}")));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The copy's Save dialog keeps Data sync off", async () => {
      await session.step(62, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(63, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(64, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be unchecked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "unchecked"));
      await session.step(65, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(66, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Turning Data sync on for the copy saves its creation script", async () => {
      await session.step(70, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(71, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(72, "When user checks Data sync switch in \"demog\" project table in \"Save project\" dialog", () => check(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog")));
      await session.step(73, "And user clicks on \"Creation script\" button in \"demog\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog")));
      await session.step(74, "Then \"demog\" project table in \"Save project\" dialog should contain text \"System:DemoFiles/demog.csv\"", () => shouldContainText(page, el("\"demog\" project table in \"Save project\" dialog"), "System:DemoFiles/demog.csv"));
      await session.step(75, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(76, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(77, "And the \"demog\" table of the \"BDD-Sync-Copy-{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDD-Sync-Copy-{time}")));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
      await session.step(79, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
