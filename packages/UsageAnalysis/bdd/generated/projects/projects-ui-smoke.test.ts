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
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, clipboardContains, clipboardHas, doubleClickOn, downloadContains, enterInto, fileDownloaded, isExpanded, shouldBe, shouldContainText, shouldHaveValue, uncheck, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelShows, dialogCloses, galleryCountLower, noProjectOnServer, pickSharingUser, projectsOnServer, rememberGalleryCount, sharingPaneLists, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project through its card in the Dashboards gallery", () => {
  const session = feature(test, "features/projects/projects-ui-smoke.feature", import.meta.url);
  test("A project through its card in the Dashboards gallery", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "And no project named \"BDDSmoke{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSmoke{time}")));
    await session.step(28, "And no project named \"BDDSmokeRenamed{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSmokeRenamed{time}")));
    await run.scenario("A file opened from its context menu is saved as a project", async () => {
      await session.step(31, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(32, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(33, "When user picks \"Open\" from the context menu of Files---Demo---demog.csv tree node inside browse tree", () => pickFromContextMenu(page, "Open", el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(34, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(35, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(36, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(37, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(38, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(39, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(40, "When user enters \"BDDSmoke{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDSmoke{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(41, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(42, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(43, "And an info balloon containing 'Project \"BDDSmoke{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDSmoke{time}\" uploaded")));
      await session.step(44, "And 1 project named \"BDDSmoke{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDSmoke{time}")));
      await session.step(45, "And \"Share BDDSmoke{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDSmoke{time}\" dialog")), "visible"));
      await session.step(46, "When user clicks on CANCEL button in \"Share BDDSmoke{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(47, "Then the \"Share BDDSmoke{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDSmoke{time}")));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The card is found by the gallery search", async () => {
      await session.step(51, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(52, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(53, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(54, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(55, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
      await session.step(56, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(57, "And user enters \"BDDSmoke{time}\" into gallery search", () => enterInto(page, session.text("BDDSmoke{time}"), el("gallery search")));
      await session.step(58, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(59, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(60, "And BDDSmoke{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSmoke{time} gallery card")), "visible"));
    });
    await run.scenario("The project is shared from its card", async () => {
      await session.step(63, "When user clicks on BDDSmoke{time} gallery card", () => clickOn(page, el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(64, "Then the context panel should show \"BDDSmoke{time}\"", () => contextPanelShows(page, session.text("BDDSmoke{time}")));
      await session.step(65, "When user picks \"Share...\" from the context menu of BDDSmoke{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(66, "Then \"Share BDDSmoke{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDSmoke{time}\" dialog")), "visible"));
      await session.step(67, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(68, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDSmoke{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(69, "And user unchecks \"Send notifications\" input in \"Share BDDSmoke{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(70, "When user clicks on OK button in \"Share BDDSmoke{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(71, "Then the \"Share BDDSmoke{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDSmoke{time}")));
      await session.step(72, "When user clicks on BDDSmoke{time} gallery card", () => clickOn(page, el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(73, "Then the context panel should show \"BDDSmoke{time}\"", () => contextPanelShows(page, session.text("BDDSmoke{time}")));
      await session.step(74, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(75, "And Sharing pane in context panel should contain text \"has special permissions\"", () => shouldContainText(page, el("Sharing pane in context panel"), "has special permissions"));
    });
    await run.scenario("The project is renamed from its card", async () => {
      await session.step(78, "When user picks \"Rename...\" from the context menu of BDDSmoke{time} gallery card", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(79, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
      await session.step(80, "And Name input in Rename project dialog should have value \"BDDSmoke{time}\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), session.text("BDDSmoke{time}")));
      await session.step(81, "When user enters \"BDDSmokeRenamed{time}\" into Name input in Rename project dialog", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("Name input in Rename project dialog")));
      await session.step(82, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(83, "Then the \"Rename project\" dialog should close", () => dialogCloses(page, "Rename project"));
      await session.step(84, "And 1 project named \"BDDSmokeRenamed{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDSmokeRenamed{time}")));
      await session.step(85, "And 0 projects named \"BDDSmoke{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDSmoke{time}")));
      await session.step(86, "When user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(87, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(88, "Then BDDSmokeRenamed{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSmokeRenamed{time} gallery card")), "visible"));
    });
    await run.scenario("The Copy items put the Grok name, the markup and the URL on the clipboard", async () => {
      await session.step(91, "When user picks \"Copy > Grok name\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Copy > Grok name", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(92, "Then the clipboard should have text \"Admin:BDDSmokeRenamed{time}\"", () => clipboardHas(page, session.text("Admin:BDDSmokeRenamed{time}")));
      await session.step(93, "When user picks \"Copy > Markup\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Copy > Markup", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(94, "Then the clipboard should have text '#{x.Admin:BDDSmokeRenamed{time}.\"BDDSmokeRenamed{time}\"}'", () => clipboardHas(page, session.text("#{x.Admin:BDDSmokeRenamed{time}.\"BDDSmokeRenamed{time}\"}")));
      await session.step(95, "When user picks \"Copy > URL\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Copy > URL", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(96, "Then the clipboard should contain text \"/p/Admin.BDDSmokeRenamed{time}\"", () => clipboardContains(page, session.text("/p/Admin.BDDSmokeRenamed{time}")));
    });
    await run.scenario("The project is added to favorites and taken out again", async () => {
      await session.step(99, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(100, "And \"My stuff > Favorites\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Favorites\" tree node inside browse tree")));
      await session.step(101, "Then \"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree")), "absent"));
      await session.step(102, "When user picks \"Add to favorites\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Add to favorites", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(103, "Then an info balloon containing \"to favorites\" should have been shown", () => infoBalloonText(page, "to favorites"));
      await session.step(104, "And \"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree")), "visible"));
      await session.step(105, "When user picks \"Remove from favorites\" from the context menu of \"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree", () => pickFromContextMenu(page, "Remove from favorites", el(session.text("\"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree"))));
      await session.step(106, "Then \"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree")), "absent"));
    });
    await run.scenario("The project is saved as a zip file", async () => {
      await session.step(109, "Given user watches downloads", () => watchDownloads(page));
      await session.step(110, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(111, "And user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(112, "And user picks \"Save as Zip\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Save as Zip", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(113, "Then a file \"BDDSmokeRenamed{time}.zip\" should have been downloaded", () => fileDownloaded(page, session.text("BDDSmokeRenamed{time}.zip")));
      await session.step(114, "And the downloaded file \"BDDSmokeRenamed{time}.zip\" should contain text \"demog\"", () => downloadContains(page, session.text("BDDSmokeRenamed{time}.zip"), "demog"));
    });
    await run.scenario("The card reopens the project with its data", async () => {
      await session.step(117, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(118, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(119, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(120, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(121, "And user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(122, "And user double-clicks on BDDSmokeRenamed{time} gallery card", () => doubleClickOn(page, el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(123, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(124, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(125, "And status bar should contain text \"5,850\"", () => shouldContainText(page, el("status bar"), "5,850"));
      await session.step(126, "And no errors should have been logged", () => noErrors(page));
      await session.step(127, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The project is deleted from its card", async () => {
      await session.step(130, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(131, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(132, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(133, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(134, "And user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(135, "Then BDDSmokeRenamed{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSmokeRenamed{time} gallery card")), "visible"));
      await session.step(136, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(137, "And user picks \"Delete Project\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(138, "Then \"Are you sure?\" dialog should contain text 'Delete project \"BDDSmokeRenamed{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDSmokeRenamed{time}\"?")));
      await session.step(139, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(140, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(141, "And 0 projects named \"BDDSmokeRenamed{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDSmokeRenamed{time}")));
      await session.step(142, "When user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(143, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(144, "And BDDSmokeRenamed{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDSmokeRenamed{time} gallery card")), "absent"));
    });
    run.finish();
  });
});
