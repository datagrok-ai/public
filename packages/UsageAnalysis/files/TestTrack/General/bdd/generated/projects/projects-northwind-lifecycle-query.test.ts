/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-northwind-lifecycle-query.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, sharing.share-dialog]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, doubleClickOn, enterInto, holdsCode, isExpanded, replaceCode, shouldBe, shouldContainText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, closeCurrentView, contextPanelShows, currentViewType, dialogCloses, noProjectOnServer, noQueryOnServer, pickSharingUser, projectsOnServer, queriesOnServer, reloadedByDataSync, savedWithDataSync, sharingPaneLists, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project of the user's own NorthwindTest query survives renaming the query", () => {
  const session = feature(test, "features/projects/projects-northwind-lifecycle-query.feature", import.meta.url);
  test("A project of the user's own NorthwindTest query survives renaming the query", {tag: ["@dev-only", "@journey", "@serial", "@realizes:views.projects", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(36, "Given user is logged in", () => loggedIn(page));
    await session.step(37, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(38, "And no project named \"BDDNwLifeQProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDNwLifeQProj{time}")));
    await session.step(39, "And no query named \"BDDNwLifeQ{time}, BDDNwLifeQRenamed{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDNwLifeQ{time}, BDDNwLifeQRenamed{time}")));
    await run.scenario("New SQL Query... on a table makes the query, saved under the feature's name", async () => {
      await session.step(42, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(43, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(44, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(45, "And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree")));
      await session.step(46, "And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
      await session.step(47, "When user picks \"New SQL Query...\" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => pickFromContextMenu(page, "New SQL Query...", el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
      await session.step(48, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(49, "And Name input should have value \"orders\"", () => shouldHaveValue(page, el("Name input"), "orders"));
      await session.step(50, "And code editor should hold the code \"select * from public.orders\"", () => holdsCode(page, el("code editor"), "select * from public.orders"));
      await session.step(51, "When user enters \"BDDNwLifeQ{time}\" into Name input", () => enterInto(page, session.text("BDDNwLifeQ{time}"), el("Name input")));
      await session.step(52, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(53, "Then 1 query named \"BDDNwLifeQ{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDNwLifeQ{time}")));
      await session.step(54, "When user closes the current view", () => closeCurrentView(page));
      await session.step(55, "Then no errors should have been logged", () => noErrors(page));
      await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The query runs from the tree into a table view", async () => {
      await session.step(59, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(60, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(61, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(62, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(63, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(64, "When user double-clicks on Databases---Postgres---NorthwindTest---BDDNwLifeQ{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("Databases---Postgres---NorthwindTest---BDDNwLifeQ{time} tree node inside browse tree"))));
      await session.step(65, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(66, "And the \"BDDNwLifeQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDNwLifeQ{time}")));
      await session.step(67, "And the table should have 830 rows", () => rowCount(page, 830));
      await session.step(68, "Then no errors should have been logged", () => noErrors(page));
      await session.step(69, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The query result is saved as a project with Data sync", async () => {
      await session.step(72, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(73, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(74, "And \"BDDNwLifeQ{time}\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDNwLifeQ{time}\" project table in \"Save project\" dialog")), "visible"));
      await session.step(75, "And Data sync switch in \"BDDNwLifeQ{time}\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el(session.text("Data sync switch in \"BDDNwLifeQ{time}\" project table in \"Save project\" dialog")), "checked"));
      await session.step(76, "When user enters \"BDDNwLifeQProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDNwLifeQProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(77, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(78, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(79, "And 1 project named \"BDDNwLifeQProj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDNwLifeQProj{time}")));
      await session.step(80, "And the \"BDDNwLifeQ{time}\" table of the \"BDDNwLifeQProj{time}\" project should be saved with data sync", () => savedWithDataSync(page, session.text("BDDNwLifeQ{time}"), session.text("BDDNwLifeQProj{time}")));
      await session.step(81, "And \"Share BDDNwLifeQProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDNwLifeQProj{time}\" dialog")), "visible"));
      await session.step(82, "When user clicks on CANCEL button in \"Share BDDNwLifeQProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDNwLifeQProj{time}\" dialog"))));
      await session.step(83, "Then the \"Share BDDNwLifeQProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDNwLifeQProj{time}")));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
      await session.step(85, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Only the project is shared with the second account, without notifications", async () => {
      await session.step(88, "When user closes all views", () => closeAllViews(page));
      await session.step(89, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(90, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(91, "And user enters \"BDDNwLifeQProj{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeQProj{time}"), el("gallery search")));
      await session.step(92, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(93, "And user picks \"Share...\" from the context menu of BDDNwLifeQProj{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDNwLifeQProj{time} gallery card"))));
      await session.step(94, "Then \"Share BDDNwLifeQProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDNwLifeQProj{time}\" dialog")), "visible"));
      await session.step(95, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(96, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDNwLifeQProj{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDNwLifeQProj{time}\" dialog"))));
      await session.step(97, "And user unchecks \"Send notifications\" input in \"Share BDDNwLifeQProj{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDNwLifeQProj{time}\" dialog"))));
      await session.step(98, "And user clicks on OK button in \"Share BDDNwLifeQProj{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDNwLifeQProj{time}\" dialog"))));
      await session.step(99, "Then the \"Share BDDNwLifeQProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDNwLifeQProj{time}")));
      await session.step(100, "When user clicks on BDDNwLifeQProj{time} gallery card", () => clickOn(page, el(session.text("BDDNwLifeQProj{time} gallery card"))));
      await session.step(101, "Then the context panel should show \"BDDNwLifeQProj{time}\"", () => contextPanelShows(page, session.text("BDDNwLifeQProj{time}")));
      await session.step(102, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
    });
    await run.scenario("The query is renamed in the Browse tree", async () => {
      await session.step(105, "When user closes all views", () => closeAllViews(page));
      await session.step(106, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(107, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(108, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(109, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(110, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(113, "When user collapses Databases---Postgres---NorthwindTest---BDDNwLifeQ{time} tree node inside browse tree", () => collapse(page, el(session.text("Databases---Postgres---NorthwindTest---BDDNwLifeQ{time} tree node inside browse tree"))));
      await session.step(114, "And user picks \"Rename...\" from the context menu of Databases---Postgres---NorthwindTest---BDDNwLifeQ{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Databases---Postgres---NorthwindTest---BDDNwLifeQ{time} tree node inside browse tree"))));
      await session.step(115, "Then \"Rename dataquery\" dialog should be visible", () => shouldBe(page, el("\"Rename dataquery\" dialog"), "visible"));
      await session.step(116, "When user enters \"BDDNwLifeQRenamed{time}\" into Name input in \"Rename dataquery\" dialog", () => enterInto(page, session.text("BDDNwLifeQRenamed{time}"), el("Name input in \"Rename dataquery\" dialog")));
      await session.step(117, "And user clicks on OK button in \"Rename dataquery\" dialog", () => clickOn(page, el("OK button in \"Rename dataquery\" dialog")));
      await session.step(118, "Then the \"Rename dataquery\" dialog should close", () => dialogCloses(page, "Rename dataquery"));
      await session.step(119, "And 1 query named \"BDDNwLifeQRenamed{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDNwLifeQRenamed{time}")));
      await session.step(120, "And 0 queries named \"BDDNwLifeQ{time}\" should be on the server", () => queriesOnServer(page, 0, session.text("BDDNwLifeQ{time}")));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
      await session.step(122, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The owner opens the project after the query was renamed (github-3550)", async () => {
      await session.step(125, "When user closes all views", () => closeAllViews(page));
      await session.step(126, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(127, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(128, "And user enters \"BDDNwLifeQProj{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeQProj{time}"), el("gallery search")));
      await session.step(129, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(130, "And user double-clicks on BDDNwLifeQProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDNwLifeQProj{time} gallery card"))));
      await session.step(131, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(132, "And the table should have 830 rows", () => rowCount(page, 830));
      await session.step(133, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(134, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
      await session.step(136, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(137, "When user closes all views", () => closeAllViews(page));
    });
    await run.scenario("With the query's SQL changed, the reopened project shows the new result", async () => {
      await session.step(140, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(141, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(142, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(143, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(144, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(145, "When user collapses Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree", () => collapse(page, el(session.text("Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree"))));
      await session.step(146, "And user picks \"Edit...\" from the context menu of Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree"))));
      await session.step(147, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(148, "When user replaces the code of code editor with \"select * from public.orders limit 100\"", () => replaceCode(page, el("code editor"), "select * from public.orders limit 100"));
      await session.step(149, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(150, "Then code editor should hold the code \"select * from public.orders limit 100\"", () => holdsCode(page, el("code editor"), "select * from public.orders limit 100"));
      await session.step(151, "When user closes the current view", () => closeCurrentView(page));
      await session.step(152, "And user reloads the page", () => reloadPage(page));
      await session.step(153, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(154, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(155, "And user enters \"BDDNwLifeQProj{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeQProj{time}"), el("gallery search")));
      await session.step(156, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(157, "And user double-clicks on BDDNwLifeQProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDNwLifeQProj{time} gallery card"))));
      await session.step(158, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(159, "And the table should have 100 rows", () => rowCount(page, 100));
      await session.step(160, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
      await session.step(161, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Delete Project removes the project, and Delete the renamed query", async () => {
      await session.step(164, "When user closes all views", () => closeAllViews(page));
      await session.step(165, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(166, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(167, "And user enters \"BDDNwLifeQProj{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeQProj{time}"), el("gallery search")));
      await session.step(168, "And user picks \"Delete Project\" from the context menu of BDDNwLifeQProj{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDNwLifeQProj{time} gallery card"))));
      await session.step(169, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(170, "And \"Are you sure?\" dialog should contain text \"Delete project \\\"BDDNwLifeQProj{time}\\\"?\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDNwLifeQProj{time}\"?")));
      await session.step(171, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(173, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(174, "And 0 projects named \"BDDNwLifeQProj{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDNwLifeQProj{time}")));
      await session.step(175, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(176, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(177, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(178, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(181, "When user collapses Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree", () => collapse(page, el(session.text("Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree"))));
      await session.step(182, "When user picks \"Delete\" from the context menu of Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete", el(session.text("Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree"))));
      await session.step(183, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(184, "And \"Are you sure?\" dialog should contain text \"Delete query \\\"BDDNwLifeQRenamed{time}\\\"?\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete query \"BDDNwLifeQRenamed{time}\"?")));
      await session.step(185, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(186, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(187, "And 0 queries named \"BDDNwLifeQRenamed{time}\" should be on the server", () => queriesOnServer(page, 0, session.text("BDDNwLifeQRenamed{time}")));
      await session.step(188, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(190, "Given Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(191, "Then Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree"), "visible"));
      await session.step(192, "And Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Databases---Postgres---NorthwindTest---BDDNwLifeQRenamed{time} tree node inside browse tree")), "absent"));
      await session.step(193, "And no errors should have been logged", () => noErrors(page));
      await session.step(194, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
