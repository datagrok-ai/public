/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-lifecycle-files.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, sharing.share-dialog]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe, shouldContainText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelShows, dialogCloses, noProjectOnServer, pickSharingUser, projectsOnServer, reloadedByDataSync, runningAccountSignedIn, savedWithDataSync, sharingPaneLists, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project built from a file, shared, renamed and deleted", () => {
  const session = feature(test, "features/projects/projects-lifecycle-files.feature", import.meta.url);
  test("A project built from a file, shared, renamed and deleted", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And no project named \"BDDLifeFiles{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeFiles{time}")));
    await session.step(27, "And no project named \"BDDLifeFilesRenamed{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeFilesRenamed{time}")));
    await run.scenario("The file is saved as a project with Data sync", async () => {
      await session.step(30, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(31, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(32, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(33, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(34, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(35, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(36, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(37, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(38, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(39, "When user enters \"BDDLifeFiles{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(40, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(41, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(42, "And an info balloon containing 'Project \"BDDLifeFiles{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDLifeFiles{time}\" uploaded")));
      await session.step(43, "And 1 project named \"BDDLifeFiles{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeFiles{time}")));
      await session.step(44, "And the \"demog\" table of the \"BDDLifeFiles{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDDLifeFiles{time}")));
      await session.step(45, "And \"Share BDDLifeFiles{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeFiles{time}\" dialog")), "visible"));
      await session.step(46, "When user clicks on CANCEL button in \"Share BDDLifeFiles{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(47, "Then the \"Share BDDLifeFiles{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeFiles{time}")));
      await session.step(48, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(49, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The project is shared to view and use, without notifications", async () => {
      await session.step(53, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(54, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(55, "And user enters \"BDDLifeFiles{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("gallery search")));
      await session.step(56, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(57, "And user picks \"Share...\" from the context menu of BDDLifeFiles{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(58, "Then \"Share BDDLifeFiles{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeFiles{time}\" dialog")), "visible"));
      await session.step(59, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(60, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDLifeFiles{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(61, "And user unchecks \"Send notifications\" input in \"Share BDDLifeFiles{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(62, "When user clicks on OK button in \"Share BDDLifeFiles{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(63, "Then the \"Share BDDLifeFiles{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeFiles{time}")));
      await session.step(64, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(65, "When user clicks on BDDLifeFiles{time} gallery card", () => clickOn(page, el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(66, "Then the context panel should show \"BDDLifeFiles{time}\"", () => contextPanelShows(page, session.text("BDDLifeFiles{time}")));
      await session.step(67, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(68, "And Sharing pane in context panel should contain text \"has special permissions\"", () => shouldContainText(page, el("Sharing pane in context panel"), "has special permissions"));
    });
    await run.scenario("With View and use the second account opens the project", async () => {
      await session.step(71, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(72, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(73, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(74, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(75, "And user enters \"BDDLifeFiles{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("gallery search")));
      await session.step(76, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(77, "And user double-clicks on BDDLifeFiles{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(78, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(79, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(80, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(81, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(82, "And user signs in as themselves again", () => signInAsSelf(page));
      await session.step(83, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(84, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(85, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(86, "And user enters \"BDDLifeFiles{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("gallery search")));
      await session.step(87, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(88, "Then BDDLifeFiles{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDLifeFiles{time} gallery card")), "visible"));
    });
    await run.scenario("The project is renamed from its card and opens under its new name", async () => {
      await session.step(91, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(92, "When user picks \"Rename...\" from the context menu of BDDLifeFiles{time} gallery card", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(93, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
      await session.step(94, "And Name input in Rename project dialog should have value \"BDDLifeFiles{time}\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), session.text("BDDLifeFiles{time}")));
      await session.step(95, "When user enters \"BDDLifeFilesRenamed{time}\" into Name input in Rename project dialog", () => enterInto(page, session.text("BDDLifeFilesRenamed{time}"), el("Name input in Rename project dialog")));
      await session.step(96, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(97, "Then the \"Rename project\" dialog should close", () => dialogCloses(page, "Rename project"));
      await session.step(98, "And 1 project named \"BDDLifeFilesRenamed{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeFilesRenamed{time}")));
      await session.step(99, "And 0 projects named \"BDDLifeFiles{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDLifeFiles{time}")));
      await session.step(100, "When user enters \"BDDLifeFilesRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFilesRenamed{time}"), el("gallery search")));
      await session.step(101, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(102, "And user double-clicks on BDDLifeFilesRenamed{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeFilesRenamed{time} gallery card"))));
      await session.step(103, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(104, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(105, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
      await session.step(107, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("After the rename, the second account finds the project under its new name and opens it", async () => {
      await session.step(110, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(111, "And user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(112, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(113, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(114, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(115, "And user enters \"BDDLifeFiles\" into gallery search", () => enterInto(page, "BDDLifeFiles", el("gallery search")));
      await session.step(116, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(117, "Then BDDLifeFilesRenamed{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDLifeFilesRenamed{time} gallery card")), "visible"));
      await session.step(118, "And BDDLifeFiles{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDLifeFiles{time} gallery card")), "absent"));
      await session.step(119, "When user double-clicks on BDDLifeFilesRenamed{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeFilesRenamed{time} gallery card"))));
      await session.step(120, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(121, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(122, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(123, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(124, "And user signs in as themselves again", () => signInAsSelf(page));
      await session.step(125, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    await run.scenario("The owner deletes the project from its card", async () => {
      await session.step(128, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(129, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(130, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(131, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(132, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(133, "And user enters \"BDDLifeFilesRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFilesRenamed{time}"), el("gallery search")));
      await session.step(134, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(135, "And user picks \"Delete Project\" from the context menu of BDDLifeFilesRenamed{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDLifeFilesRenamed{time} gallery card"))));
      await session.step(136, "Then \"Are you sure?\" dialog should contain text 'Delete project \"BDDLifeFilesRenamed{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDLifeFilesRenamed{time}\"?")));
      await session.step(137, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(138, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(139, "And 0 projects named \"BDDLifeFilesRenamed{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDLifeFilesRenamed{time}")));
      await session.step(140, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(141, "Then BDDLifeFilesRenamed{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDLifeFilesRenamed{time} gallery card")), "absent"));
    });
    run.finish();
  });
});
