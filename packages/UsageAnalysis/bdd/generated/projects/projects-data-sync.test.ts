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
import '../../bindings/nx.js';
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
import {browsePanelOpen, currentViewType, dialogCloses, loadedAsSnapshot, noProjectOnServer, projectsOnServer, reloadedByDataSync, savedAsSnapshot, savedWithDataSync, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project saved with and without Data sync", () => {
  const session = feature(test, "features/projects/projects-data-sync.feature", import.meta.url);
  test("A project saved with and without Data sync", {tag: ["@journey", "@serial", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And no project named \"BDDSyncOrig{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSyncOrig{time}")));
    await session.step(19, "And no project named \"BDDSyncCopy{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSyncCopy{time}")));
    await run.scenario("A file opened from Browse is saved with Data sync on", async () => {
      await session.step(22, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(23, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(24, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(25, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(26, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(27, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(28, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(29, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(30, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(31, "When user enters \"BDDSyncOrig{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDSyncOrig{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(32, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(33, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(34, "And an info balloon containing 'Project \"BDDSyncOrig{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDSyncOrig{time}\" uploaded")));
      await session.step(35, "And 1 project named \"BDDSyncOrig{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDSyncOrig{time}")));
      await session.step(36, "And the \"demog\" table of the \"BDDSyncOrig{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDDSyncOrig{time}")));
      await session.step(37, "And \"Share BDDSyncOrig{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDSyncOrig{time}\" dialog")), "visible"));
      await session.step(38, "When user clicks on CANCEL button in \"Share BDDSyncOrig{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDSyncOrig{time}\" dialog"))));
      await session.step(39, "Then the \"Share BDDSyncOrig{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDSyncOrig{time}")));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A copy saved with Data sync off leaves the original synced, and is the open project", async () => {
      await session.step(43, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(44, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(45, "When user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(46, "Then Name text input in \"Save project\" dialog should have value \"Copy of BDDSyncOrig{time}\"", () => shouldHaveValue(page, el("Name text input in \"Save project\" dialog"), session.text("Copy of BDDSyncOrig{time}")));
      await session.step(47, "And choice input in \"demog\" project table in \"Save project\" dialog should have value \"Clone\"", () => shouldHaveValue(page, el("choice input in \"demog\" project table in \"Save project\" dialog"), "Clone"));
      await session.step(48, "When user enters \"BDDSyncCopy{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDSyncCopy{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(49, "And user unchecks Data sync switch in \"demog\" project table in \"Save project\" dialog", () => uncheck(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog")));
      await session.step(50, "Then \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be hidden", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "hidden"));
      await session.step(51, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(52, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(53, "And an info balloon containing 'Project \"BDDSyncCopy{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDSyncCopy{time}\" uploaded")));
      await session.step(54, "And 1 project named \"BDDSyncCopy{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDSyncCopy{time}")));
      await session.step(55, "And the \"demog\" table of the \"BDDSyncCopy{time}\" project should be saved as a snapshot", () => savedAsSnapshot(page, "demog", session.text("BDDSyncCopy{time}")));
      await session.step(56, "And the \"demog\" table of the \"BDDSyncOrig{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDDSyncOrig{time}")));
      await session.step(57, "When user clicks on Dashboards tab", () => clickOn(page, el("Dashboards tab")));
      await session.step(58, "Then BDDSyncCopy{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDSyncCopy{time} tree node inside browse tree")), "visible"));
      await session.step(59, "And BDDSyncOrig{time} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("BDDSyncOrig{time} tree node inside browse tree")), "absent"));
      await session.step(61, "When user clicks on Dashboards tab", () => clickOn(page, el("Dashboards tab")));
      await session.step(62, "Then \"New Dashboard\" tree node inside browse tree should be hidden", () => shouldBe(page, el("\"New Dashboard\" tree node inside browse tree"), "hidden"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The copy opens from its card as a snapshot, with no creation script", async () => {
      await session.step(66, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(67, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(68, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(69, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(70, "And user enters \"BDDSync\" into gallery search", () => enterInto(page, "BDDSync", el("gallery search")));
      await session.step(71, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(72, "Then BDDSyncOrig{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSyncOrig{time} gallery card")), "visible"));
      await session.step(73, "And BDDSyncCopy{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSyncCopy{time} gallery card")), "visible"));
      await session.step(74, "When user double-clicks on BDDSyncCopy{time} gallery card", () => doubleClickOn(page, el(session.text("BDDSyncCopy{time} gallery card"))));
      await session.step(75, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(76, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(77, "And the table should have been loaded as a snapshot", () => loadedAsSnapshot(page));
      await session.step(78, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(79, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(80, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be unchecked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "unchecked"));
      await session.step(81, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be hidden", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "hidden"));
      await session.step(82, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(83, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Data sync turned on for the copy, the copy rereads the file", async () => {
      await session.step(87, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(88, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(89, "When user checks Data sync switch in \"demog\" project table in \"Save project\" dialog", () => check(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog")));
      await session.step(90, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(91, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(92, "And the \"demog\" table of the \"BDDSyncCopy{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDDSyncCopy{time}")));
      await session.step(93, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(94, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(95, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(96, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(97, "And user enters \"BDDSyncCopy{time}\" into gallery search", () => enterInto(page, session.text("BDDSyncCopy{time}"), el("gallery search")));
      await session.step(98, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(99, "And user double-clicks on BDDSyncCopy{time} gallery card", () => doubleClickOn(page, el(session.text("BDDSyncCopy{time} gallery card"))));
      await session.step(100, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(101, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(102, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(103, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(104, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(105, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(106, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(107, "When user clicks on \"Creation script\" button in \"demog\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog")));
      await session.step(108, "Then \"demog\" project table in \"Save project\" dialog should contain text 'OpenFile(\"System:DemoFiles/demog.csv\")'", () => shouldContainText(page, el("\"demog\" project table in \"Save project\" dialog"), "OpenFile(\"System:DemoFiles/demog.csv\")"));
      await session.step(109, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(110, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(111, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The original still opens by rereading the file, with its creation script", async () => {
      await session.step(114, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(115, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(116, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(117, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(118, "And user enters \"BDDSyncOrig{time}\" into gallery search", () => enterInto(page, session.text("BDDSyncOrig{time}"), el("gallery search")));
      await session.step(119, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(120, "And user double-clicks on BDDSyncOrig{time} gallery card", () => doubleClickOn(page, el(session.text("BDDSyncOrig{time} gallery card"))));
      await session.step(121, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(122, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(123, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(124, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(125, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(126, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(127, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(128, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(129, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
      await session.step(131, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
