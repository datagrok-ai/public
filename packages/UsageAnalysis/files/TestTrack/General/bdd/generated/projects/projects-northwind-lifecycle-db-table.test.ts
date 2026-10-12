/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-northwind-lifecycle-db-table.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, sharing.share-dialog]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, contextPanelShows, currentViewType, dialogCloses, galleryCountLower, noProjectOnServer, pickSharingUser, projectsOnServer, reloadedByDataSync, rememberGalleryCount, runningAccountSignedIn, savedWithDataSync, sharingPaneLists, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project of the NorthwindTest orders table: saved, reopened and shared", () => {
  const session = feature(test, "features/projects/projects-northwind-lifecycle-db-table.feature", import.meta.url);
  test("A project of the NorthwindTest orders table: saved, reopened and shared", {tag: ["@dev-only", "@journey", "@serial", "@realizes:views.projects", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(33, "Given user is logged in", () => loggedIn(page));
    await session.step(34, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(35, "And no project named \"BDDNwLifeDbTable{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDNwLifeDbTable{time}")));
    await run.scenario("Get All opens the database table into a table view", async () => {
      await session.step(38, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(39, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(40, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(41, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(42, "And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree")));
      await session.step(43, "And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
      await session.step(44, "When user picks \"Get All\" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => pickFromContextMenu(page, "Get All", el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
      await session.step(45, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(46, "And the \"orders\" view should be current", () => viewIsCurrent(page, "orders"));
      await session.step(47, "And the table should have 830 rows", () => rowCount(page, 830));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
      await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Saved with Data sync, the creation script reads the table from the database", async () => {
      await session.step(52, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(53, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(54, "And Data sync switch in \"orders\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"orders\" project table in \"Save project\" dialog"), "checked"));
      await session.step(55, "When user enters \"BDDNwLifeDbTable{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDNwLifeDbTable{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(56, "And user clicks on \"Creation script\" button in \"orders\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"orders\" project table in \"Save project\" dialog")));
      await session.step(57, "Then \"orders\" project table in \"Save project\" dialog should contain text \"DbQuery(Dbtests:PostgresTest, \\\"public.orders\\\"\"", () => shouldContainText(page, el("\"orders\" project table in \"Save project\" dialog"), "DbQuery(Dbtests:PostgresTest, \"public.orders\""));
      await session.step(58, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(59, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(60, "And an info balloon containing 'Project \"BDDNwLifeDbTable{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDNwLifeDbTable{time}\" uploaded")));
      await session.step(61, "And 1 project named \"BDDNwLifeDbTable{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDNwLifeDbTable{time}")));
      await session.step(62, "And the \"orders\" table of the \"BDDNwLifeDbTable{time}\" project should be saved with data sync", () => savedWithDataSync(page, "orders", session.text("BDDNwLifeDbTable{time}")));
      await session.step(63, "And \"Share BDDNwLifeDbTable{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDNwLifeDbTable{time}\" dialog")), "visible"));
      await session.step(64, "When user clicks on CANCEL button in \"Share BDDNwLifeDbTable{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDNwLifeDbTable{time}\" dialog"))));
      await session.step(65, "Then the \"Share BDDNwLifeDbTable{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDNwLifeDbTable{time}")));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The table project reopens from Dashboards by re-reading the table", async () => {
      await session.step(69, "When user closes all views", () => closeAllViews(page));
      await session.step(70, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(71, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(72, "And user enters \"BDDNwLifeDbTable{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbTable{time}"), el("gallery search")));
      await session.step(73, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(74, "And user double-clicks on BDDNwLifeDbTable{time} gallery card", () => doubleClickOn(page, el(session.text("BDDNwLifeDbTable{time} gallery card"))));
      await session.step(75, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(76, "And the \"orders\" view should be current", () => viewIsCurrent(page, "orders"));
      await session.step(77, "And the table should have 830 rows", () => rowCount(page, 830));
      await session.step(78, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
      await session.step(80, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(81, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(82, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(83, "And \"Creation script\" button in \"orders\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"orders\" project table in \"Save project\" dialog"), "visible"));
      await session.step(84, "And Data sync switch in \"orders\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"orders\" project table in \"Save project\" dialog"), "checked"));
      await session.step(85, "When user clicks on \"Creation script\" button in \"orders\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"orders\" project table in \"Save project\" dialog")));
      await session.step(86, "Then \"orders\" project table in \"Save project\" dialog should contain text \"DbQuery(Dbtests:PostgresTest, \\\"public.orders\\\"\"", () => shouldContainText(page, el("\"orders\" project table in \"Save project\" dialog"), "DbQuery(Dbtests:PostgresTest, \"public.orders\""));
      await session.step(87, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(88, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    });
    await run.scenario("Before the share, the second account does not see the NorthwindTest connection", async () => {
      await session.step(91, "When user closes all views", () => closeAllViews(page));
      await session.step(92, "And user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(93, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(94, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(95, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(96, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(97, "Then Databases---Postgres---Datagrok tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---Datagrok tree node inside browse tree"), "visible"));
      await session.step(98, "And Databases---Postgres---NorthwindTest tree node inside browse tree should be absent", () => shouldBe(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree"), "absent"));
      await session.step(99, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(100, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    await run.scenario("Only the project is shared with the second account", async () => {
      await session.step(103, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(104, "When user closes all views", () => closeAllViews(page));
      await session.step(105, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(106, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(107, "And user enters \"BDDNwLifeDbTable{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbTable{time}"), el("gallery search")));
      await session.step(108, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(109, "And user picks \"Share...\" from the context menu of BDDNwLifeDbTable{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDNwLifeDbTable{time} gallery card"))));
      await session.step(110, "Then \"Share BDDNwLifeDbTable{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDNwLifeDbTable{time}\" dialog")), "visible"));
      await session.step(111, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(112, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDNwLifeDbTable{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDNwLifeDbTable{time}\" dialog"))));
      await session.step(113, "And user clicks on OK button in \"Share BDDNwLifeDbTable{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDNwLifeDbTable{time}\" dialog"))));
      await session.step(114, "Then the \"Share BDDNwLifeDbTable{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDNwLifeDbTable{time}")));
      await session.step(115, "When user clicks on BDDNwLifeDbTable{time} gallery card", () => clickOn(page, el(session.text("BDDNwLifeDbTable{time} gallery card"))));
      await session.step(116, "Then the context panel should show \"BDDNwLifeDbTable{time}\"", () => contextPanelShows(page, session.text("BDDNwLifeDbTable{time}")));
      await session.step(117, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
    });
    await run.scenario("With the project shared, the second account sees the connection and opens the project with its rows", async () => {
      await session.step(120, "When user closes all views", () => closeAllViews(page));
      await session.step(121, "And user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(122, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(123, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(124, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(125, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(126, "Then Databases---Postgres---NorthwindTest tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree"), "visible"));
      await session.step(127, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(128, "And user enters \"BDDNwLifeDbTable{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbTable{time}"), el("gallery search")));
      await session.step(129, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(130, "And user double-clicks on BDDNwLifeDbTable{time} gallery card", () => doubleClickOn(page, el(session.text("BDDNwLifeDbTable{time} gallery card"))));
      await session.step(131, "Then the \"orders\" view should be current", () => viewIsCurrent(page, "orders"));
      await session.step(132, "And the table should have 830 rows", () => rowCount(page, 830));
      await session.step(133, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
      await session.step(134, "When user closes all views", () => closeAllViews(page));
      await session.step(135, "And user signs in as themselves again", () => signInAsSelf(page));
      await session.step(136, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    await run.scenario("Delete Project removes the table project", async () => {
      await session.step(139, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(140, "When user closes all views", () => closeAllViews(page));
      await session.step(141, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(142, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(143, "And user enters \"BDDNwLifeDbTable{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbTable{time}"), el("gallery search")));
      await session.step(144, "Then BDDNwLifeDbTable{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDNwLifeDbTable{time} gallery card")), "visible"));
      await session.step(145, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(146, "And user picks \"Delete Project\" from the context menu of BDDNwLifeDbTable{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDNwLifeDbTable{time} gallery card"))));
      await session.step(147, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(148, "And \"Are you sure?\" dialog should contain text \"Delete project \\\"BDDNwLifeDbTable{time}\\\"?\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDNwLifeDbTable{time}\"?")));
      await session.step(149, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(150, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(151, "And 0 projects named \"BDDNwLifeDbTable{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDNwLifeDbTable{time}")));
      await session.step(152, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(153, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(154, "And BDDNwLifeDbTable{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDNwLifeDbTable{time} gallery card")), "absent"));
      await session.step(155, "And no errors should have been logged", () => noErrors(page));
      await session.step(156, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
