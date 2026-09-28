/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-ui-smoke.feature
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
import {clipboardHoldsProjectForm} from '../../bindings/projects-sharing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, clipboardContains, doubleClickOn, downloadContains, enterInto, fileDownloaded, isExpanded, shouldBe, shouldContainText, shouldHaveValue, uncheck, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelShows, dialogCloses, galleryCountLower, noProjectOnServer, pickSharingUser, projectsOnServer, rememberGalleryCount, sharingPaneLists, sharingPaneListsNot, sharingUserShownAs, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project through its card in the Dashboards gallery", () => {
  const session = feature(test, "features/projects/projects-ui-smoke.feature", import.meta.url);
  test("A project through its card in the Dashboards gallery", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And no project named \"BDDSmoke{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSmoke{time}")));
    await session.step(25, "And no project named \"BDDSmokeRenamed{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSmokeRenamed{time}")));
    await run.scenario("A file opened from its context menu is saved as a project", async () => {
      await session.step(28, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(29, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(30, "When user picks \"Open\" from the context menu of Files---Demo---demog.csv tree node inside browse tree", () => pickFromContextMenu(page, "Open", el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(31, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(32, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(33, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(34, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(35, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(36, "When user enters \"BDDSmoke{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDSmoke{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(37, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(38, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(39, "And an info balloon containing 'Project \"BDDSmoke{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDSmoke{time}\" uploaded")));
      await session.step(40, "And 1 project named \"BDDSmoke{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDSmoke{time}")));
      await session.step(41, "And \"Share BDDSmoke{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDSmoke{time}\" dialog")), "visible"));
      await session.step(42, "When user clicks on CANCEL button in \"Share BDDSmoke{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(43, "Then the \"Share BDDSmoke{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDSmoke{time}")));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The card is found by the gallery search", async () => {
      await session.step(47, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(48, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(49, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(50, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(51, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
      await session.step(52, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(53, "And user enters \"BDDSmoke{time}\" into gallery search", () => enterInto(page, session.text("BDDSmoke{time}"), el("gallery search")));
      await session.step(54, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(55, "And BDDSmoke{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSmoke{time} gallery card")), "visible"));
    });
    await run.scenario("The project is shared from its card", async () => {
      await session.step(58, "When user clicks on BDDSmoke{time} gallery card", () => clickOn(page, el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(59, "Then the context panel should show \"BDDSmoke{time}\"", () => contextPanelShows(page, session.text("BDDSmoke{time}")));
      await session.step(60, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(61, "When user picks \"Share...\" from the context menu of BDDSmoke{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(62, "Then \"Share BDDSmoke{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDSmoke{time}\" dialog")), "visible"));
      await session.step(63, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(64, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDSmoke{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(65, "And user unchecks \"Send notifications\" input in \"Share BDDSmoke{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(66, "Then \"Send notifications\" input in \"Share BDDSmoke{time}\" dialog should be unchecked", () => shouldBe(page, el(session.text("\"Send notifications\" input in \"Share BDDSmoke{time}\" dialog")), "unchecked"));
      await session.step(67, "When user clicks on OK button in \"Share BDDSmoke{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(68, "Then the \"Share BDDSmoke{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDSmoke{time}")));
      await session.step(69, "When user clicks on BDDSmoke{time} gallery card", () => clickOn(page, el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(70, "Then the context panel should show \"BDDSmoke{time}\"", () => contextPanelShows(page, session.text("BDDSmoke{time}")));
      await session.step(71, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(72, "And the sharing pane should show the sharing user as \"has special permissions\"", () => sharingUserShownAs(page, "has special permissions"));
    });
    await run.scenario("The project is renamed from its card", async () => {
      await session.step(75, "When user picks \"Rename...\" from the context menu of BDDSmoke{time} gallery card", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(76, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
      await session.step(77, "And Name input in Rename project dialog should have value \"BDDSmoke{time}\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), session.text("BDDSmoke{time}")));
      await session.step(78, "When user enters \"BDDSmokeRenamed{time}\" into Name input in Rename project dialog", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("Name input in Rename project dialog")));
      await session.step(79, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(80, "Then the \"Rename project\" dialog should close", () => dialogCloses(page, "Rename project"));
      await session.step(81, "And 1 project named \"BDDSmokeRenamed{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDSmokeRenamed{time}")));
      await session.step(82, "And 0 projects named \"BDDSmoke{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDSmoke{time}")));
      await session.step(83, "When user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(84, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(85, "Then BDDSmokeRenamed{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSmokeRenamed{time} gallery card")), "visible"));
      await session.step(86, "And BDDSmoke{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDSmoke{time} gallery card")), "absent"));
    });
    await run.scenario("Copy puts the id, the grok name, the markup and the link on the clipboard", async () => {
      await session.step(89, "When user picks \"Copy > ID\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Copy > ID", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(90, "Then the clipboard should hold the id of the \"BDDSmokeRenamed{time}\" project", () => clipboardHoldsProjectForm(page, "id", session.text("BDDSmokeRenamed{time}")));
      await session.step(91, "When user picks \"Copy > Grok name\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Copy > Grok name", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(92, "Then the clipboard should hold the grok-name of the \"BDDSmokeRenamed{time}\" project", () => clipboardHoldsProjectForm(page, "grok-name", session.text("BDDSmokeRenamed{time}")));
      await session.step(93, "And the clipboard should contain text \":BDDSmokeRenamed{time}\"", () => clipboardContains(page, session.text(":BDDSmokeRenamed{time}")));
      await session.step(94, "When user picks \"Copy > Markup\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Copy > Markup", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(95, "Then the clipboard should hold the markup of the \"BDDSmokeRenamed{time}\" project", () => clipboardHoldsProjectForm(page, "markup", session.text("BDDSmokeRenamed{time}")));
      await session.step(96, "And the clipboard should contain text '.\"BDDSmokeRenamed{time}\"}'", () => clipboardContains(page, session.text(".\"BDDSmokeRenamed{time}\"}")));
      await session.step(97, "When user picks \"Copy > URL\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Copy > URL", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(98, "Then the clipboard should hold the url of the \"BDDSmokeRenamed{time}\" project", () => clipboardHoldsProjectForm(page, "url", session.text("BDDSmokeRenamed{time}")));
      await session.step(99, "And the clipboard should contain text \".BDDSmokeRenamed{time}\"", () => clipboardContains(page, session.text(".BDDSmokeRenamed{time}")));
    });
    await run.scenario("The project is added to favorites and taken out again", async () => {
      await session.step(102, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(103, "And \"My stuff > Favorites\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Favorites\" tree node inside browse tree")));
      await session.step(104, "Then \"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree")), "absent"));
      await session.step(105, "When user picks \"Add to favorites\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Add to favorites", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(106, "Then an info balloon containing \"to favorites\" should have been shown", () => infoBalloonText(page, "to favorites"));
      await session.step(107, "And \"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree")), "visible"));
      await session.step(108, "When user picks \"Remove from favorites\" from the context menu of \"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree", () => pickFromContextMenu(page, "Remove from favorites", el(session.text("\"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree"))));
      await session.step(109, "Then \"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree")), "absent"));
    });
    await run.scenario("The project is saved as a zip file", async () => {
      await session.step(112, "Given user watches downloads", () => watchDownloads(page));
      await session.step(113, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(114, "And user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(115, "And user picks \"Save as Zip\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Save as Zip", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(116, "Then a file \"BDDSmokeRenamed{time}.zip\" should have been downloaded", () => fileDownloaded(page, session.text("BDDSmokeRenamed{time}.zip")));
      await session.step(117, "And the downloaded file \"BDDSmokeRenamed{time}.zip\" should contain text \"demog\"", () => downloadContains(page, session.text("BDDSmokeRenamed{time}.zip"), "demog"));
    });
    await run.scenario("The card reopens the project with its data", async () => {
      await session.step(120, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(121, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(122, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(123, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(124, "And user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(125, "And user double-clicks on BDDSmokeRenamed{time} gallery card", () => doubleClickOn(page, el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(126, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(127, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(128, "And status bar should contain text \"5,850\"", () => shouldContainText(page, el("status bar"), "5,850"));
      await session.step(129, "And no errors should have been logged", () => noErrors(page));
      await session.step(130, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The project is deleted from its card", async () => {
      await session.step(133, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(134, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(135, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(136, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(137, "And user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(138, "Then BDDSmokeRenamed{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSmokeRenamed{time} gallery card")), "visible"));
      await session.step(139, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(140, "And user picks \"Delete Project\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(141, "Then \"Are you sure?\" dialog should contain text 'Delete project \"BDDSmokeRenamed{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDSmokeRenamed{time}\"?")));
      await session.step(142, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(143, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(144, "And 0 projects named \"BDDSmokeRenamed{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDSmokeRenamed{time}")));
      await session.step(145, "When user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(146, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(147, "And BDDSmokeRenamed{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDSmokeRenamed{time} gallery card")), "absent"));
    });
    run.finish();
  });
});
