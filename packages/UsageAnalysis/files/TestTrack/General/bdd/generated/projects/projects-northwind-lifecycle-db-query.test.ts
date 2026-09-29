/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-northwind-lifecycle-db-query.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, currentViewType, dialogCloses, galleryCountLower, noProjectOnServer, projectsOnServer, reloadedByDataSync, rememberGalleryCount, savedWithDataSync, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project of the NorthwindTest query PostgresAll: saved and reopened with its creation script", () => {
  const session = feature(test, "features/projects/projects-northwind-lifecycle-db-query.feature", import.meta.url);
  test("A project of the NorthwindTest query PostgresAll: saved and reopened with its creation script", {tag: ["@dev-only", "@journey", "@serial", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "And no project named \"BDDNwLifeDbQuery{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDNwLifeDbQuery{time}")));
    await run.scenario("The saved query runs from the tree into a table view", async () => {
      await session.step(30, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(31, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(32, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(33, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(34, "When user double-clicks on Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree", () => doubleClickOn(page, el("Databases---Postgres---NorthwindTest---PostgresAll tree node inside browse tree")));
      await session.step(35, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(36, "And the \"PostgresAll\" view should be current", () => viewIsCurrent(page, "PostgresAll"));
      await session.step(37, "And the table should have 830 rows", () => rowCount(page, 830));
      await session.step(38, "Then no errors should have been logged", () => noErrors(page));
      await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Saved with Data sync, the creation script calls the query", async () => {
      await session.step(42, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(43, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(44, "And Data sync switch in \"PostgresAll\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"PostgresAll\" project table in \"Save project\" dialog"), "checked"));
      await session.step(45, "When user enters \"BDDNwLifeDbQuery{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDNwLifeDbQuery{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(46, "And user clicks on \"Creation script\" button in \"PostgresAll\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"PostgresAll\" project table in \"Save project\" dialog")));
      await session.step(47, "Then \"PostgresAll\" project table in \"Save project\" dialog should contain text \":PostgresAll()\"", () => shouldContainText(page, el("\"PostgresAll\" project table in \"Save project\" dialog"), ":PostgresAll()"));
      await session.step(48, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(49, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(50, "And an info balloon containing 'Project \"BDDNwLifeDbQuery{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDNwLifeDbQuery{time}\" uploaded")));
      await session.step(51, "And 1 project named \"BDDNwLifeDbQuery{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDNwLifeDbQuery{time}")));
      await session.step(52, "And the \"PostgresAll\" table of the \"BDDNwLifeDbQuery{time}\" project should be saved with data sync", () => savedWithDataSync(page, "PostgresAll", session.text("BDDNwLifeDbQuery{time}")));
      await session.step(53, "And \"Share BDDNwLifeDbQuery{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDNwLifeDbQuery{time}\" dialog")), "visible"));
      await session.step(54, "When user clicks on CANCEL button in \"Share BDDNwLifeDbQuery{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDNwLifeDbQuery{time}\" dialog"))));
      await session.step(55, "Then the \"Share BDDNwLifeDbQuery{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDNwLifeDbQuery{time}")));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The query project reopens from Dashboards by re-running the query", async () => {
      await session.step(59, "When user closes all views", () => closeAllViews(page));
      await session.step(60, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(61, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(62, "And user enters \"BDDNwLifeDbQuery{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbQuery{time}"), el("gallery search")));
      await session.step(63, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(64, "And user double-clicks on BDDNwLifeDbQuery{time} gallery card", () => doubleClickOn(page, el(session.text("BDDNwLifeDbQuery{time} gallery card"))));
      await session.step(65, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(66, "And the \"PostgresAll\" view should be current", () => viewIsCurrent(page, "PostgresAll"));
      await session.step(67, "And the table should have 830 rows", () => rowCount(page, 830));
      await session.step(68, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
      await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(71, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(72, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(73, "And \"Creation script\" button in \"PostgresAll\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"PostgresAll\" project table in \"Save project\" dialog"), "visible"));
      await session.step(74, "And Data sync switch in \"PostgresAll\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"PostgresAll\" project table in \"Save project\" dialog"), "checked"));
      await session.step(75, "When user clicks on \"Creation script\" button in \"PostgresAll\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"PostgresAll\" project table in \"Save project\" dialog")));
      await session.step(76, "Then \"PostgresAll\" project table in \"Save project\" dialog should contain text \":PostgresAll()\"", () => shouldContainText(page, el("\"PostgresAll\" project table in \"Save project\" dialog"), ":PostgresAll()"));
      await session.step(77, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(78, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    });
    await run.scenario("Delete Project removes the query project", async () => {
      await session.step(81, "When user closes all views", () => closeAllViews(page));
      await session.step(82, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(83, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(84, "And user enters \"BDDNwLifeDbQuery{time}\" into gallery search", () => enterInto(page, session.text("BDDNwLifeDbQuery{time}"), el("gallery search")));
      await session.step(85, "Then BDDNwLifeDbQuery{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDNwLifeDbQuery{time} gallery card")), "visible"));
      await session.step(86, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(87, "And user picks \"Delete Project\" from the context menu of BDDNwLifeDbQuery{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDNwLifeDbQuery{time} gallery card"))));
      await session.step(88, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(89, "And \"Are you sure?\" dialog should contain text \"Delete project \\\"BDDNwLifeDbQuery{time}\\\"?\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDNwLifeDbQuery{time}\"?")));
      await session.step(90, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(92, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(93, "And 0 projects named \"BDDNwLifeDbQuery{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDNwLifeDbQuery{time}")));
      await session.step(94, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(95, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(96, "And BDDNwLifeDbQuery{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDNwLifeDbQuery{time} gallery card")), "absent"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
      await session.step(98, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
