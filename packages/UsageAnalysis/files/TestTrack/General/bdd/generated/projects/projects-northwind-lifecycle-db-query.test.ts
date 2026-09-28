/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-northwind-lifecycle-db-query.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects]
--- */
import {test} from '@playwright/test';
import '../../bindings/pkg-projects-copies.js';
import '../../bindings/pkg-projects-derived.js';
import '../../bindings/pkg-projects-regressions.js';
import '../../bindings/pkg-projects-sources.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {sharingUserUsesConnection} from '../../bindings/northwind.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, contextPanelOpen, contextPanelShows, currentViewType, dialogCloses, galleryCountLower, noProjectOnServer, pickSharingUser, projectsOnServer, reloadedByDataSync, rememberGalleryCount, savedWithDataSync, sharingPaneLists, sharingUserAccess, signInAsSecond, signInAsSelf, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project of the NorthwindTest query PostgresAll: saved, reopened, shared and renamed", () => {
  const session = feature(test, "features/projects/projects-northwind-lifecycle-db-query.feature", import.meta.url);
  test("A project of the NorthwindTest query PostgresAll: saved, reopened, shared and renamed", {tag: ["@dev-only", "@journey", "@serial", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(38, "Given user is logged in", () => loggedIn(page));
    await session.step(39, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(40, "And the sharing user may use the \"Dbtests:PostgresTest\" connection until the feature ends", () => sharingUserUsesConnection(page, "Dbtests:PostgresTest"));
    await session.step(41, "And no project named \"BDDNwLifeDbQuery{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDNwLifeDbQuery{time}")));
    await session.step(42, "And no project named \"BDDNwLifeDbQueryRenamed{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDNwLifeDbQueryRenamed{time}")));
    await run.scenario("The saved query runs from the tree into a table view", async () => {
      await session.step(45, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(46, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(47, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(48, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(49, "When user double-clicks on Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree", () => doubleClickOn(page, el("Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree")));
      await session.step(50, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(51, "And the \"PostgresAll\" view should be current", () => viewIsCurrent(page, "PostgresAll"));
      await session.step(52, "And the table should have 830 rows", () => rowCount(page, 830));
      await session.step(53, "Then no errors should have been logged", () => noErrors(page));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Saved with Data sync, the creation script calls the query", async () => {
      await session.step(57, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(58, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(59, "And Data sync switch in \"PostgresAll\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"PostgresAll\" project table in \"Save project\" dialog"), "checked"));
      await session.step(60, "When user enters \"BDDNwLifeDbQuery{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDNwLifeDbQuery{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(61, "And user clicks on \"Creation script\" button in \"PostgresAll\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"PostgresAll\" project table in \"Save project\" dialog")));
      await session.step(62, "Then creation script text in \"PostgresAll\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("creation script text in \"PostgresAll\" project table in \"Save project\" dialog"), "visible"));
      await session.step(63, "And creation script text in \"PostgresAll\" project table in \"Save project\" dialog should contain text \":PostgresAll()\"", () => shouldContainText(page, el("creation script text in \"PostgresAll\" project table in \"Save project\" dialog"), ":PostgresAll()"));
      await session.step(64, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(65, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(66, "And an info balloon containing 'Project \"BDDNwLifeDbQuery{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDNwLifeDbQuery{time}\" uploaded")));
      await session.step(67, "And 1 project named \"BDDNwLifeDbQuery{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDNwLifeDbQuery{time}")));
      await session.step(68, "And the \"PostgresAll\" table of the \"BDDNwLifeDbQuery{time}\" project should be saved with data sync", () => savedWithDataSync(page, "PostgresAll", session.text("BDDNwLifeDbQuery{time}")));
      await session.step(69, "And \"Share BDDNwLifeDbQuery{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDNwLifeDbQuery{time}\" dialog")), "visible"));
      await session.step(70, "When user clicks on CANCEL button in \"Share BDDNwLifeDbQuery{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDNwLifeDbQuery{time}\" dialog"))));
      await session.step(71, "Then the \"Share BDDNwLifeDbQuery{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDNwLifeDbQuery{time}")));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The query project reopens from Dashboards by re-running the query", async () => {
      await session.step(75, "When user closes all views", () => closeAllViews(page));
      await session.step(76, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(77, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(78, "And user enters \"BDDNwLifeDbQuery{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbQuery{time}"), el("gallery search")));
      await session.step(79, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(80, "And user double-clicks on BDDNwLifeDbQuery{time} gallery card", () => doubleClickOn(page, el(session.text("BDDNwLifeDbQuery{time} gallery card"))));
      await session.step(81, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(82, "And the \"PostgresAll\" view should be current", () => viewIsCurrent(page, "PostgresAll"));
      await session.step(83, "And the table should have 830 rows", () => rowCount(page, 830));
      await session.step(84, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
      await session.step(86, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The query project is shared with the second account", async () => {
      await session.step(89, "When user closes all views", () => closeAllViews(page));
      await session.step(90, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(91, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(92, "And user enters \"BDDNwLifeDbQuery{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbQuery{time}"), el("gallery search")));
      await session.step(93, "And user picks \"Share...\" from the context menu of BDDNwLifeDbQuery{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDNwLifeDbQuery{time} gallery card"))));
      await session.step(94, "Then \"Share BDDNwLifeDbQuery{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDNwLifeDbQuery{time}\" dialog")), "visible"));
      await session.step(96, "And \"Share BDDNwLifeDbQuery{time}\" dialog should contain text \"Full access\"", () => shouldContainText(page, el(session.text("\"Share BDDNwLifeDbQuery{time}\" dialog")), "Full access"));
      await session.step(97, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDNwLifeDbQuery{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDNwLifeDbQuery{time}\" dialog"))));
      await session.step(98, "Then share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(99, "When user clicks on OK button in \"Share BDDNwLifeDbQuery{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDNwLifeDbQuery{time}\" dialog"))));
      await session.step(100, "Then the \"Share BDDNwLifeDbQuery{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDNwLifeDbQuery{time}")));
      await session.step(101, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(102, "When user clicks on BDDNwLifeDbQuery{time} gallery card", () => clickOn(page, el(session.text("BDDNwLifeDbQuery{time} gallery card"))));
      await session.step(103, "Then the context panel should show \"BDDNwLifeDbQuery{time}\"", () => contextPanelShows(page, session.text("BDDNwLifeDbQuery{time}")));
      await session.step(104, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(105, "When user picks \"Share...\" from the context menu of BDDNwLifeDbQuery{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDNwLifeDbQuery{time} gallery card"))));
      await session.step(106, "Then the access level of the sharing user in \"Share BDDNwLifeDbQuery{time}\" dialog should be \"View and use\"", () => sharingUserAccess(page, el(session.text("\"Share BDDNwLifeDbQuery{time}\" dialog")), "View and use"));
      await session.step(107, "When user clicks on CANCEL button in \"Share BDDNwLifeDbQuery{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDNwLifeDbQuery{time}\" dialog"))));
      await session.step(108, "Then the \"Share BDDNwLifeDbQuery{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDNwLifeDbQuery{time}")));
      await session.step(109, "And no errors should have been logged", () => noErrors(page));
      await session.step(110, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The second account opens the shared query project", async () => {
      await session.step(113, "Given user signs in as the sharing user", () => signInAsSecond(page));
      await session.step(114, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(115, "And the sharing user may use the \"Dbtests:PostgresTest\" connection until the feature ends", () => sharingUserUsesConnection(page, "Dbtests:PostgresTest"));
      await session.step(116, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(117, "And user enters \"BDDNwLifeDbQuery{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbQuery{time}"), el("gallery search")));
      await session.step(118, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(119, "And user double-clicks on BDDNwLifeDbQuery{time} gallery card", () => doubleClickOn(page, el(session.text("BDDNwLifeDbQuery{time} gallery card"))));
      await session.step(120, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(121, "And the \"PostgresAll\" view should be current", () => viewIsCurrent(page, "PostgresAll"));
      await session.step(122, "And the table should have 830 rows", () => rowCount(page, 830));
      await session.step(123, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(124, "And no errors should have been logged", () => noErrors(page));
      await session.step(125, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Renamed, the query project opens under its new name", async () => {
      await session.step(128, "Given user signs in as themselves again", () => signInAsSelf(page));
      await session.step(129, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(130, "And the sharing user may use the \"Dbtests:PostgresTest\" connection until the feature ends", () => sharingUserUsesConnection(page, "Dbtests:PostgresTest"));
      await session.step(131, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(132, "And user enters \"BDDNwLifeDbQuery{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbQuery{time}"), el("gallery search")));
      await session.step(133, "And user picks \"Rename...\" from the context menu of BDDNwLifeDbQuery{time} gallery card", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDNwLifeDbQuery{time} gallery card"))));
      await session.step(134, "Then \"Rename project\" dialog should be visible", () => shouldBe(page, el("\"Rename project\" dialog"), "visible"));
      await session.step(135, "When user enters \"BDDNwLifeDbQueryRenamed{time}\" into Name input in \"Rename project\" dialog", () => enterInto(page, session.text("BDDNwLifeDbQueryRenamed{time}"), el("Name input in \"Rename project\" dialog")));
      await session.step(136, "And user clicks on OK button in \"Rename project\" dialog", () => clickOn(page, el("OK button in \"Rename project\" dialog")));
      await session.step(137, "Then the \"Rename project\" dialog should close", () => dialogCloses(page, "Rename project"));
      await session.step(138, "And 1 project named \"BDDNwLifeDbQueryRenamed{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDNwLifeDbQueryRenamed{time}")));
      await session.step(139, "And 0 projects named \"BDDNwLifeDbQuery{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDNwLifeDbQuery{time}")));
      await session.step(140, "When user enters \"BDDNwLifeDbQueryRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbQueryRenamed{time}"), el("gallery search")));
      await session.step(141, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(142, "And user double-clicks on BDDNwLifeDbQueryRenamed{time} gallery card", () => doubleClickOn(page, el(session.text("BDDNwLifeDbQueryRenamed{time} gallery card"))));
      await session.step(143, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(144, "And the table should have 830 rows", () => rowCount(page, 830));
      await session.step(145, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(146, "And no errors should have been logged", () => noErrors(page));
      await session.step(147, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Delete Project removes the query project", async () => {
      await session.step(150, "When user closes all views", () => closeAllViews(page));
      await session.step(151, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(152, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(153, "And user enters \"BDDNwLifeDbQueryRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbQueryRenamed{time}"), el("gallery search")));
      await session.step(154, "Then BDDNwLifeDbQueryRenamed{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDNwLifeDbQueryRenamed{time} gallery card")), "visible"));
      await session.step(155, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(156, "And user picks \"Delete Project\" from the context menu of BDDNwLifeDbQueryRenamed{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDNwLifeDbQueryRenamed{time} gallery card"))));
      await session.step(157, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(158, "And \"Are you sure?\" dialog should contain text \"Delete project \\\"BDDNwLifeDbQueryRenamed{time}\\\"?\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDNwLifeDbQueryRenamed{time}\"?")));
      await session.step(159, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(161, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(162, "And 0 projects named \"BDDNwLifeDbQueryRenamed{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDNwLifeDbQueryRenamed{time}")));
      await session.step(163, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(164, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(165, "And BDDNwLifeDbQueryRenamed{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDNwLifeDbQueryRenamed{time} gallery card")), "absent"));
      await session.step(166, "And no errors should have been logged", () => noErrors(page));
      await session.step(167, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
