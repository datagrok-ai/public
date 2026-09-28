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
import {projectSavedAgain, rememberProjectSaved} from '../../bindings/projects-sharing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, enterInto, followingShouldBe, isExpanded, shouldBe, shouldContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closePrivilegeTree, contextPanelShows, dialogCloses, noProjectOnServer, openSharingUserAccess, pickSharingUser, projectsOnServer, sharingPaneLists, sharingUserAccess, sharingUserShownAs, signInAsSecond, signInAsSelf, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A file-based project shared at two access levels and renamed", () => {
  const session = feature(test, "features/projects/projects-lifecycle-files.feature", import.meta.url);
  test("A file-based project shared at two access levels and renamed", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And no project named \"BDDLifeFiles{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeFiles{time}")));
    await session.step(24, "And no project named \"BDDLifeFilesRenamed{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeFilesRenamed{time}")));
    await run.scenario("The file is saved as a project", async () => {
      await session.step(27, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(28, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(29, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(30, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(31, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(32, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(33, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(34, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(35, "When user enters \"BDDLifeFiles{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(36, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(37, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(38, "And an info balloon containing 'Project \"BDDLifeFiles{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDLifeFiles{time}\" uploaded")));
      await session.step(39, "And 1 project named \"BDDLifeFiles{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeFiles{time}")));
      await session.step(40, "And \"Share BDDLifeFiles{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeFiles{time}\" dialog")), "visible"));
      await session.step(41, "When user clicks on CANCEL button in \"Share BDDLifeFiles{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(42, "Then the \"Share BDDLifeFiles{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeFiles{time}")));
      await session.step(43, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(44, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The owner shares it to view and use", async () => {
      await session.step(48, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(49, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(50, "And user enters \"BDDLifeFiles{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("gallery search")));
      await session.step(51, "And user picks \"Share...\" from the context menu of BDDLifeFiles{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(52, "Then \"Share BDDLifeFiles{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeFiles{time}\" dialog")), "visible"));
      await session.step(53, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(54, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDLifeFiles{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(55, "And user unchecks \"Send notifications\" input in \"Share BDDLifeFiles{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(56, "And user clicks on OK button in \"Share BDDLifeFiles{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(57, "Then the \"Share BDDLifeFiles{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeFiles{time}")));
      await session.step(58, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(59, "When user clicks on BDDLifeFiles{time} gallery card", () => clickOn(page, el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(60, "Then the context panel should show \"BDDLifeFiles{time}\"", () => contextPanelShows(page, session.text("BDDLifeFiles{time}")));
      await session.step(61, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(62, "And the sharing pane should show the sharing user as \"has special permissions\"", () => sharingUserShownAs(page, "has special permissions"));
    });
    await run.scenario("With View and use the second account opens it and can only save a copy", async () => {
      await session.step(65, "Given user signs in as the sharing user", () => signInAsSecond(page));
      await session.step(66, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(67, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(68, "And user enters \"BDDLifeFiles{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("gallery search")));
      await session.step(69, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(70, "And user double-clicks on BDDLifeFiles{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(71, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(72, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(74, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(75, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(76, "And \"Save original project\" radio choice in \"Save project\" dialog should be disabled", () => shouldBe(page, el("\"Save original project\" radio choice in \"Save project\" dialog"), "disabled"));
      await session.step(77, "And \"Save a copy\" radio choice in \"Save project\" dialog should be checked", () => shouldBe(page, el("\"Save a copy\" radio choice in \"Save project\" dialog"), "checked"));
      await session.step(78, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(79, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(80, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(81, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The owner raises the share to full access", async () => {
      await session.step(85, "Given user signs in as themselves again", () => signInAsSelf(page));
      await session.step(86, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(87, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(88, "And user enters \"BDDLifeFiles{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("gallery search")));
      await session.step(89, "And user picks \"Share...\" from the context menu of BDDLifeFiles{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(90, "Then \"Share BDDLifeFiles{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeFiles{time}\" dialog")), "visible"));
      await session.step(91, "And the access level of the sharing user in \"Share BDDLifeFiles{time}\" dialog should be \"View and use\"", () => sharingUserAccess(page, el(session.text("\"Share BDDLifeFiles{time}\" dialog")), "View and use"));
      await session.step(92, "When user opens the access level of the sharing user in \"Share BDDLifeFiles{time}\" dialog", () => openSharingUserAccess(page, el(session.text("\"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(93, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Full-access tree node inside privilege tree"],["View-and-use tree node inside privilege tree"],["Edit tree node inside privilege tree"],["Delete tree node inside privilege tree"],["Share tree node inside privilege tree"]]), [["Full-access tree node inside privilege tree"],["View-and-use tree node inside privilege tree"],["Edit tree node inside privilege tree"],["Delete tree node inside privilege tree"],["Share tree node inside privilege tree"]]);
      await session.step(99, "And Full-access tree node inside privilege tree should be unchecked", () => shouldBe(page, el("Full-access tree node inside privilege tree"), "unchecked"));
      await session.step(100, "When user checks Full-access tree node inside privilege tree", () => check(page, el("Full-access tree node inside privilege tree")));
      await session.step(101, "And user clicks outside the privilege tree", () => closePrivilegeTree(page));
      await session.step(102, "Then the access level of the sharing user in \"Share BDDLifeFiles{time}\" dialog should be \"Full access\"", () => sharingUserAccess(page, el(session.text("\"Share BDDLifeFiles{time}\" dialog")), "Full access"));
      await session.step(103, "When user clicks on OK button in \"Share BDDLifeFiles{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(104, "Then the \"Share BDDLifeFiles{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeFiles{time}")));
      await session.step(105, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
    });
    await run.scenario("The owner renames the project and it still opens", async () => {
      await session.step(108, "When user picks \"Rename...\" from the context menu of BDDLifeFiles{time} gallery card", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(109, "And user enters \"BDDLifeFilesRenamed{time}\" into Name input in Rename project dialog", () => enterInto(page, session.text("BDDLifeFilesRenamed{time}"), el("Name input in Rename project dialog")));
      await session.step(110, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(111, "Then the \"Rename project\" dialog should close", () => dialogCloses(page, "Rename project"));
      await session.step(112, "And 1 project named \"BDDLifeFilesRenamed{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeFilesRenamed{time}")));
      await session.step(113, "And 0 projects named \"BDDLifeFiles{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDLifeFiles{time}")));
      await session.step(114, "When user enters \"BDDLifeFilesRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFilesRenamed{time}"), el("gallery search")));
      await session.step(115, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(116, "And user double-clicks on BDDLifeFilesRenamed{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeFilesRenamed{time} gallery card"))));
      await session.step(117, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(118, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(119, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(120, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Full access the second account opens the renamed project and saves the original", async () => {
      await session.step(124, "Given user signs in as the sharing user", () => signInAsSecond(page));
      await session.step(125, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(126, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(127, "And user enters \"BDDLifeFiles\" into gallery search", () => enterInto(page, "BDDLifeFiles", el("gallery search")));
      await session.step(128, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(129, "Then BDDLifeFilesRenamed{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDLifeFilesRenamed{time} gallery card")), "visible"));
      await session.step(130, "When user double-clicks on BDDLifeFilesRenamed{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeFilesRenamed{time} gallery card"))));
      await session.step(131, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(132, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(133, "When user remembers when the \"BDDLifeFilesRenamed{time}\" project was saved", () => rememberProjectSaved(page, session.text("BDDLifeFilesRenamed{time}")));
      await session.step(134, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(135, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(136, "And \"Save original project\" radio choice in \"Save project\" dialog should be enabled", () => shouldBe(page, el("\"Save original project\" radio choice in \"Save project\" dialog"), "enabled"));
      await session.step(137, "And \"Save original project\" radio choice in \"Save project\" dialog should be checked", () => shouldBe(page, el("\"Save original project\" radio choice in \"Save project\" dialog"), "checked"));
      await session.step(138, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(139, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(140, "And an info balloon containing 'Project \"BDDLifeFilesRenamed{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDLifeFilesRenamed{time}\" uploaded")));
      await session.step(141, "And 1 project named \"BDDLifeFilesRenamed{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeFilesRenamed{time}")));
      await session.step(142, "And the \"BDDLifeFilesRenamed{time}\" project should have been saved again since remembered", () => projectSavedAgain(page, session.text("BDDLifeFilesRenamed{time}")));
      await session.step(143, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(144, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(145, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The owner deletes the project from its card", async () => {
      await session.step(148, "Given user signs in as themselves again", () => signInAsSelf(page));
      await session.step(149, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(150, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(151, "And user enters \"BDDLifeFilesRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFilesRenamed{time}"), el("gallery search")));
      await session.step(152, "And user picks \"Delete Project\" from the context menu of BDDLifeFilesRenamed{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDLifeFilesRenamed{time} gallery card"))));
      await session.step(153, "Then \"Are you sure?\" dialog should contain text 'Delete project \"BDDLifeFilesRenamed{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDLifeFilesRenamed{time}\"?")));
      await session.step(154, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(155, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(156, "And 0 projects named \"BDDLifeFilesRenamed{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDLifeFilesRenamed{time}")));
    });
    run.finish();
  });
});
