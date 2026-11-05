/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-lifecycle-files.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, sharing.share-dialog]
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
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe, shouldContainText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, doubleClickOpensProject, noProjectOnServer, pickSharingUser, projectsOnServer, reloadedByDataSync, runningAccountSignedIn, savedWithDataSync, sharingPaneLists, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project built from a file, shared and renamed", () => {
  const session = feature(test, "features/projects/projects-lifecycle-files.feature", import.meta.url);
  test("A project built from a file, shared and renamed", {tag: ["@journey", "@realizes:views.projects", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "And no project named \"BDDLifeFiles{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeFiles{time}")));
    await session.step(28, "And no project named \"BDDLifeFilesRenamed{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeFilesRenamed{time}")));
    await run.scenario("The file is saved as a project with Data sync", async () => {
      await session.step(31, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(32, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(33, "When user double-clicks Files---Demo---cars.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---cars.csv tree node inside browse tree")));
      await session.step(34, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
      await session.step(35, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(36, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(37, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(38, "And \"Creation script\" button in \"cars\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"cars\" project table in \"Save project\" dialog"), "visible"));
      await session.step(39, "And Data sync switch in \"cars\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"cars\" project table in \"Save project\" dialog"), "checked"));
      await session.step(40, "When user enters \"BDDLifeFiles{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(41, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(42, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(43, "And an info balloon containing 'Project \"BDDLifeFiles{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDLifeFiles{time}\" uploaded")));
      await session.step(44, "And 1 project named \"BDDLifeFiles{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeFiles{time}")));
      await session.step(45, "And the \"cars\" table of the \"BDDLifeFiles{time}\" project should be saved with data sync", () => savedWithDataSync(page, "cars", session.text("BDDLifeFiles{time}")));
      await session.step(46, "And \"Share BDDLifeFiles{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeFiles{time}\" dialog")), "visible"));
      await session.step(47, "When user clicks on CANCEL button in \"Share BDDLifeFiles{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(48, "Then the \"Share BDDLifeFiles{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeFiles{time}")));
      await session.step(49, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(50, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The project is shared to view and use, without notifications", async () => {
      await session.step(54, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(55, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(56, "And user enters \"BDDLifeFiles{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("gallery search")));
      await session.step(57, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(58, "And user picks \"Share...\" from the context menu of BDDLifeFiles{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(59, "Then \"Share BDDLifeFiles{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeFiles{time}\" dialog")), "visible"));
      await session.step(60, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(61, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDLifeFiles{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(62, "And user unchecks \"Send notifications\" input in \"Share BDDLifeFiles{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(63, "When user clicks on OK button in \"Share BDDLifeFiles{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDLifeFiles{time}\" dialog"))));
      await session.step(64, "Then the \"Share BDDLifeFiles{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeFiles{time}")));
      await session.step(65, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(66, "When user clicks on BDDLifeFiles{time} gallery card", () => clickOn(page, el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(67, "Then the context panel should show \"BDDLifeFiles{time}\"", () => contextPanelShows(page, session.text("BDDLifeFiles{time}")));
      await session.step(68, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(69, "And Sharing pane in context panel should contain text \"can view and use\"", () => shouldContainText(page, el("Sharing pane in context panel"), "can view and use"));
    });
    await run.scenario("With View and use the second account opens the project", async () => {
      await session.step(72, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(73, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(74, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(75, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(76, "And user enters \"BDDLifeFiles{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("gallery search")));
      await session.step(77, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(78, "And user double-clicks on BDDLifeFiles{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(79, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
      await session.step(80, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(81, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(82, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(83, "And user signs in as themselves again", () => signInAsSelf(page));
      await session.step(84, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(85, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(86, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(87, "And user enters \"BDDLifeFiles{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFiles{time}"), el("gallery search")));
      await session.step(88, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(89, "Then BDDLifeFiles{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDLifeFiles{time} gallery card")), "visible"));
    });
    await run.scenario("The project is renamed from its card and opens under its new name", async () => {
      await session.step(92, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(93, "When user picks \"Rename...\" from the context menu of BDDLifeFiles{time} gallery card", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDLifeFiles{time} gallery card"))));
      await session.step(94, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
      await session.step(95, "And Name input in Rename project dialog should have value \"BDDLifeFiles{time}\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), session.text("BDDLifeFiles{time}")));
      await session.step(96, "When user enters \"BDDLifeFilesRenamed{time}\" into Name input in Rename project dialog", () => enterInto(page, session.text("BDDLifeFilesRenamed{time}"), el("Name input in Rename project dialog")));
      await session.step(97, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(98, "Then the \"Rename project\" dialog should close", () => dialogCloses(page, "Rename project"));
      await session.step(99, "And 1 project named \"BDDLifeFilesRenamed{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeFilesRenamed{time}")));
      await session.step(100, "And 0 projects named \"BDDLifeFiles{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDLifeFiles{time}")));
      await session.step(101, "When user enters \"BDDLifeFilesRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeFilesRenamed{time}"), el("gallery search")));
      await session.step(102, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(103, "And user double-clicks on BDDLifeFilesRenamed{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDLifeFilesRenamed{time} gallery card"))));
      await session.step(104, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
      await session.step(105, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(106, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
      await session.step(108, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("After the rename, the second account finds the project under its new name and opens it", async () => {
      await session.step(111, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(112, "And user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(113, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(114, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(115, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(116, "And user enters \"BDDLifeFiles\" into gallery search", () => enterInto(page, "BDDLifeFiles", el("gallery search")));
      await session.step(117, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(118, "Then BDDLifeFilesRenamed{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDLifeFilesRenamed{time} gallery card")), "visible"));
      await session.step(119, "And BDDLifeFiles{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDLifeFiles{time} gallery card")), "absent"));
      await session.step(120, "When user double-clicks on BDDLifeFilesRenamed{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDLifeFilesRenamed{time} gallery card"))));
      await session.step(121, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
      await session.step(122, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(123, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(124, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(125, "And user signs in as themselves again", () => signInAsSelf(page));
      await session.step(126, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    run.finish();
  });
});
