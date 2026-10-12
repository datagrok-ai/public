/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-ui-smoke.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, sharing.share-dialog]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, clipboardContains, clipboardHas, collapse, doubleClickOn, downloadContains, enterInto, expand, fileDownloaded, hoverOver, isExpanded, pressKey, rememberVisibleCount, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, typeAtCaret, uncheck, visibleCount, visibleMoreThanRemembered, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {favoriteOnServer, notFavorite, notFavoriteOnServer} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelShows, dialogCloses, galleryCountLower, noProjectOnServer, pickSharingUser, projectsOnServer, rememberGalleryCount, sharingPaneLists, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project through its card in the Dashboards gallery", () => {
  const session = feature(test, "features/projects/projects-ui-smoke.feature", import.meta.url);
  test("A project through its card in the Dashboards gallery", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(32, "And no project named \"BDDSmoke{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSmoke{time}")));
    await session.step(33, "And no project named \"BDDSmokeRenamed{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSmokeRenamed{time}")));
    await run.scenario("A file opened from its context menu is saved as a project", async () => {
      await session.step(36, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(37, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(38, "When user picks \"Open\" from the context menu of Files---Demo---demog.csv tree node inside browse tree", () => pickFromContextMenu(page, "Open", el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(39, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(40, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(41, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(42, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(43, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(44, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(45, "When user enters \"BDDSmoke{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDSmoke{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(46, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(47, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(48, "And an info balloon containing 'Project \"BDDSmoke{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDSmoke{time}\" uploaded")));
      await session.step(49, "And 1 project named \"BDDSmoke{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDSmoke{time}")));
      await session.step(50, "And \"Share BDDSmoke{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDSmoke{time}\" dialog")), "visible"));
      await session.step(51, "When user clicks on CANCEL button in \"Share BDDSmoke{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(52, "Then the \"Share BDDSmoke{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDSmoke{time}")));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The card is found by the gallery search", async () => {
      await session.step(56, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(57, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(58, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(59, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(60, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
      await session.step(61, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(62, "And user enters \"BDDSmoke{time}\" into gallery search", () => enterInto(page, session.text("BDDSmoke{time}"), el("gallery search")));
      await session.step(63, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(64, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(65, "And BDDSmoke{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSmoke{time} gallery card")), "visible"));
    });
    await run.scenario("A tag is added in the Details pane and the \"#tag\" search finds the project by it", async () => {
      await session.step(68, "When user clicks on BDDSmoke{time} gallery card", () => clickOn(page, el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(69, "Then the context panel should show \"BDDSmoke{time}\"", () => contextPanelShows(page, session.text("BDDSmoke{time}")));
      await session.step(70, "When user hovers over Tags details row in context panel", () => hoverOver(page, el("Tags details row in context panel")));
      await session.step(71, "And user hovers over edit icon of Tags details row in context panel", () => hoverOver(page, el("edit icon of Tags details row in context panel")));
      await session.step(72, "Then tooltip should contain text \"Edit tags and metadata\"", () => shouldContainText(page, el("tooltip"), "Edit tags and metadata"));
      await session.step(73, "When user clicks on edit icon of Tags details row in context panel", () => clickOn(page, el("edit icon of Tags details row in context panel")));
      await session.step(75, "And user clicks on value of Tags details row in context panel", () => clickOn(page, el("value of Tags details row in context panel")));
      await session.step(76, "And user types \"bddtag{time}\" at the caret", () => typeAtCaret(page, session.text("bddtag{time}")));
      await session.step(77, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(78, "And user presses Escape", () => pressKey(page, "Escape"));
      await session.step(79, "Then Details pane in context panel should contain text \"#bddtag{time}\"", () => shouldContainText(page, el("Details pane in context panel"), session.text("#bddtag{time}")));
      await session.step(80, "When user remembers the number of visible gallery card", () => rememberVisibleCount(page, el("gallery card")));
      await session.step(81, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(82, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(83, "Then there should be more visible gallery card than remembered", () => visibleMoreThanRemembered(page, el("gallery card")));
      await session.step(84, "When user enters \"#bddtag{time}\" into gallery search", () => enterInto(page, session.text("#bddtag{time}"), el("gallery search")));
      await session.step(85, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(86, "Then there should be 1 visible gallery card", () => visibleCount(page, 1, el("gallery card")));
      await session.step(87, "And BDDSmoke{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSmoke{time} gallery card")), "visible"));
      await session.step(88, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(89, "And user enters \"#bddnotag{time}\" into gallery search", () => enterInto(page, session.text("#bddnotag{time}"), el("gallery search")));
      await session.step(90, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(91, "Then gallery counter should have text \"0\"", () => shouldHaveText(page, el("gallery counter"), "0"));
      await session.step(92, "And there should be 0 visible gallery card", () => visibleCount(page, 0, el("gallery card")));
      await session.step(93, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(94, "And user enters \"BDDSmoke{time}\" into gallery search", () => enterInto(page, session.text("BDDSmoke{time}"), el("gallery search")));
      await session.step(95, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(96, "Then BDDSmoke{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSmoke{time} gallery card")), "visible"));
    });
    await run.scenario("The project is shared from its card", async () => {
      await session.step(99, "When user clicks on BDDSmoke{time} gallery card", () => clickOn(page, el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(100, "Then the context panel should show \"BDDSmoke{time}\"", () => contextPanelShows(page, session.text("BDDSmoke{time}")));
      await session.step(101, "When user picks \"Share...\" from the context menu of BDDSmoke{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(102, "Then \"Share BDDSmoke{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDSmoke{time}\" dialog")), "visible"));
      await session.step(103, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(104, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDSmoke{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(105, "And user unchecks \"Send notifications\" input in \"Share BDDSmoke{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(106, "When user clicks on OK button in \"Share BDDSmoke{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDSmoke{time}\" dialog"))));
      await session.step(107, "Then the \"Share BDDSmoke{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDSmoke{time}")));
      await session.step(108, "When user clicks on BDDSmoke{time} gallery card", () => clickOn(page, el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(109, "Then the context panel should show \"BDDSmoke{time}\"", () => contextPanelShows(page, session.text("BDDSmoke{time}")));
      await session.step(110, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(111, "And Sharing pane in context panel should contain text \"has special permissions\"", () => shouldContainText(page, el("Sharing pane in context panel"), "has special permissions"));
    });
    await run.scenario("The project is renamed from its card", async () => {
      await session.step(114, "When user picks \"Rename...\" from the context menu of BDDSmoke{time} gallery card", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDSmoke{time} gallery card"))));
      await session.step(115, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
      await session.step(116, "And Name input in Rename project dialog should have value \"BDDSmoke{time}\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), session.text("BDDSmoke{time}")));
      await session.step(117, "When user enters \"BDDSmokeRenamed{time}\" into Name input in Rename project dialog", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("Name input in Rename project dialog")));
      await session.step(118, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(119, "Then the \"Rename project\" dialog should close", () => dialogCloses(page, "Rename project"));
      await session.step(120, "And 1 project named \"BDDSmokeRenamed{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDSmokeRenamed{time}")));
      await session.step(121, "And 0 projects named \"BDDSmoke{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDSmoke{time}")));
      await session.step(122, "When user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(123, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(124, "Then BDDSmokeRenamed{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSmokeRenamed{time} gallery card")), "visible"));
    });
    await run.scenario("The Copy items put the Grok name, the markup and the URL on the clipboard", async () => {
      await session.step(127, "When user picks \"Copy > Grok name\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Copy > Grok name", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(128, "Then the clipboard should have text \"Admin:BDDSmokeRenamed{time}\"", () => clipboardHas(page, session.text("Admin:BDDSmokeRenamed{time}")));
      await session.step(129, "When user picks \"Copy > Markup\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Copy > Markup", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(130, "Then the clipboard should have text '#{x.Admin:BDDSmokeRenamed{time}.\"BDDSmokeRenamed{time}\"}'", () => clipboardHas(page, session.text("#{x.Admin:BDDSmokeRenamed{time}.\"BDDSmokeRenamed{time}\"}")));
      await session.step(131, "When user picks \"Copy > URL\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Copy > URL", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(132, "Then the clipboard should contain text \"/p/Admin.BDDSmokeRenamed{time}\"", () => clipboardContains(page, session.text("/p/Admin.BDDSmokeRenamed{time}")));
    });
    await run.scenario("The project is added to favorites and taken out again", async () => {
      await session.step(135, "Given \"BDDSmokeRenamed{time}\" is not in favorites", () => notFavorite(page, session.text("BDDSmokeRenamed{time}")));
      await session.step(136, "And \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(137, "And \"My stuff > Favorites\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Favorites\" tree node inside browse tree")));
      await session.step(138, "Then \"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree")), "absent"));
      await session.step(139, "When user picks \"Add To Favorites > Only for me\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(140, "Then \"BDDSmokeRenamed{time}\" should be in favorites on the server", () => favoriteOnServer(page, session.text("BDDSmokeRenamed{time}")));
      await session.step(141, "When user collapses \"My stuff > Favorites\" tree node inside browse tree", () => collapse(page, el("\"My stuff > Favorites\" tree node inside browse tree")));
      await session.step(142, "And user expands \"My stuff > Favorites\" tree node inside browse tree", () => expand(page, el("\"My stuff > Favorites\" tree node inside browse tree")));
      await session.step(143, "Then \"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree")), "visible"));
      await session.step(144, "When user picks \"Add To Favorites > Only for me\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(145, "Then \"BDDSmokeRenamed{time}\" should not be in favorites on the server", () => notFavoriteOnServer(page, session.text("BDDSmokeRenamed{time}")));
      await session.step(146, "When user collapses \"My stuff > Favorites\" tree node inside browse tree", () => collapse(page, el("\"My stuff > Favorites\" tree node inside browse tree")));
      await session.step(147, "And user expands \"My stuff > Favorites\" tree node inside browse tree", () => expand(page, el("\"My stuff > Favorites\" tree node inside browse tree")));
      await session.step(148, "Then \"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"My stuff > Favorites > BDDSmokeRenamed{time}\" tree node inside browse tree")), "absent"));
    });
    await run.scenario("The project is saved as a zip file", async () => {
      await session.step(151, "Given user watches downloads", () => watchDownloads(page));
      await session.step(152, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(153, "And user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(154, "And user picks \"Save as Zip\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Save as Zip", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(155, "Then a file \"BDDSmokeRenamed{time}.zip\" should have been downloaded", () => fileDownloaded(page, session.text("BDDSmokeRenamed{time}.zip")));
      await session.step(156, "And the downloaded file \"BDDSmokeRenamed{time}.zip\" should contain text \"demog\"", () => downloadContains(page, session.text("BDDSmokeRenamed{time}.zip"), "demog"));
    });
    await run.scenario("The card reopens the project with its data", async () => {
      await session.step(159, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(160, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(161, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(162, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(163, "And user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(164, "And user double-clicks on BDDSmokeRenamed{time} gallery card", () => doubleClickOn(page, el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(165, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(166, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(167, "And status bar should contain text \"5,850\"", () => shouldContainText(page, el("status bar"), "5,850"));
      await session.step(168, "And no errors should have been logged", () => noErrors(page));
      await session.step(169, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The project is deleted from its card", async () => {
      await session.step(172, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(173, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(174, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(175, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(176, "And user enters \"BDDSmokeRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDSmokeRenamed{time}"), el("gallery search")));
      await session.step(177, "Then BDDSmokeRenamed{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSmokeRenamed{time} gallery card")), "visible"));
      await session.step(178, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(179, "And user picks \"Delete Project\" from the context menu of BDDSmokeRenamed{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDSmokeRenamed{time} gallery card"))));
      await session.step(180, "Then \"Are you sure?\" dialog should contain text 'Delete project \"BDDSmokeRenamed{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDSmokeRenamed{time}\"?")));
      await session.step(181, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(182, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(183, "And 0 projects named \"BDDSmokeRenamed{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDSmokeRenamed{time}")));
      await session.step(184, "When user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(185, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(186, "And BDDSmokeRenamed{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDSmokeRenamed{time} gallery card")), "absent"));
    });
    run.finish();
  });
});
