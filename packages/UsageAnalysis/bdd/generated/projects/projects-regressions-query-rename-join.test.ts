/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-query-rename-join.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, data.menu.join-tables, GROK-21026]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/minimized-viewers.js';
import '../../bindings/projects-copies.js';
import '../../bindings/projects-derived.js';
import '../../bindings/projects-sources.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noErrorsButPreviewNoise, tableViewTabs} from '../../bindings/projects-regressions.js';
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, doubleClickOn, enterInto, followingShouldBe, holdsCode, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickColumnIn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, closeAllViews, closeCurrentView, creationScriptHolds, currentViewType, dialogCloses, noProjectOnServer, noQueryOnServer, projectHoldsTables, queriesOnServer, rememberRows, savedWithDataSync, tableReloadedAsRemembered, tableRowsAsRemembered, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewsOpen} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a join over a query result after the query is renamed", () => {
  const session = feature(test, "features/projects/projects-regressions-query-rename-join.feature", import.meta.url);
  test("Projects regressions: a join over a query result after the query is renamed", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:data.menu.join-tables", "@known-failure", "@realizes:GROK-21026"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(39, "Given user is logged in", () => loggedIn(page));
    await session.step(40, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(41, "And no project named \"BDDRenJoinProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRenJoinProj{time}")));
    await session.step(42, "And no query named \"BDDRenJoinQ{time}, BDDRenJoinQR{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDRenJoinQ{time}, BDDRenJoinQR{time}")));
    await run.scenario("A join of entity_types with a query's result is saved with Data sync, and the query is renamed in the tree", async () => {
      await session.step(45, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(46, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(47, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(48, "And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas tree node inside browse tree")));
      await session.step(49, "And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree")));
      await session.step(50, "When user picks \"New SQL Query...\" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree", () => pickFromContextMenu(page, "New SQL Query...", el("Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree")));
      await session.step(51, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(52, "And code editor should hold the code \"select * from public.entity_types\"", () => holdsCode(page, el("code editor"), "select * from public.entity_types"));
      await session.step(53, "When user enters \"BDDRenJoinQ{time}\" into Name input", () => enterInto(page, session.text("BDDRenJoinQ{time}"), el("Name input")));
      await session.step(54, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(55, "Then 1 query named \"BDDRenJoinQ{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDRenJoinQ{time}")));
      await session.step(56, "When user closes the current view", () => closeCurrentView(page));
      await session.step(57, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(58, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(59, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(60, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(61, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(62, "When user double-clicks on Databases---Postgres---Datagrok---BDDRenJoinQ{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("Databases---Postgres---Datagrok---BDDRenJoinQ{time} tree node inside browse tree"))));
      await session.step(63, "Then the \"BDDRenJoinQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDRenJoinQ{time}")));
      await session.step(64, "When user remembers the row count of the table as \"rename join query\"", () => rememberRows(page, "rename join query"));
      await session.step(65, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(66, "And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas tree node inside browse tree")));
      await session.step(67, "And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree")));
      await session.step(68, "When user picks \"Get All\" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree", () => pickFromContextMenu(page, "Get All", el("Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree")));
      await session.step(69, "Then the \"entity_types\" view should be current", () => viewIsCurrent(page, "entity_types"));
      await session.step(70, "And table \"entity_types\" should have the \"rename join query\" row count", () => tableRowsAsRemembered(page, "entity_types", "rename join query"));
      await session.step(71, "When user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
      await session.step(72, "Then \"Join Tables\" dialog should be visible", () => shouldBe(page, el("\"Join Tables\" dialog"), "visible"));
      await session.step(73, "When user selects \"entity_types\" in join left table selector", () => selectIn(page, "entity_types", el("join left table selector")));
      await session.step(74, "And user selects \"BDDRenJoinQ{time}\" in join right table selector", () => selectIn(page, session.text("BDDRenJoinQ{time}"), el("join right table selector")));
      await session.step(75, "And user picks column \"id\" in join left key selector", () => pickColumnIn(page, "id", el("join left key selector")));
      await session.step(76, "And user picks column \"id\" in join right key selector", () => pickColumnIn(page, "id", el("join right key selector")));
      await session.step(77, "And user clicks on OK button in \"Join Tables\" dialog", () => clickOn(page, el("OK button in \"Join Tables\" dialog")));
      await session.step(78, "Then the \"Join Tables\" dialog should close", () => dialogCloses(page, "Join Tables"));
      await session.step(79, "And the \"result\" view should be current", () => viewIsCurrent(page, "result"));
      await session.step(80, "And table \"result\" should have the \"rename join query\" row count", () => tableRowsAsRemembered(page, "result", "rename join query"));
      await session.step(81, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(82, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(83, "And the following elements should be checked:", () => followingShouldBe(page, "checked", [[session.text("Data sync switch in \"BDDRenJoinQ{time}\" project table in \"Save project\" dialog")],["Data sync switch in \"entity_types\" project table in \"Save project\" dialog"],["Data sync switch in \"result\" project table in \"Save project\" dialog"]]), [[session.text("Data sync switch in \"BDDRenJoinQ{time}\" project table in \"Save project\" dialog")],["Data sync switch in \"entity_types\" project table in \"Save project\" dialog"],["Data sync switch in \"result\" project table in \"Save project\" dialog"]]);
      await session.step(87, "When user enters \"BDDRenJoinProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRenJoinProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(88, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(89, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(90, "And \"Share BDDRenJoinProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRenJoinProj{time}\" dialog")), "visible"));
      await session.step(91, "When user clicks on CANCEL button in \"Share BDDRenJoinProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRenJoinProj{time}\" dialog"))));
      await session.step(92, "Then the \"Share BDDRenJoinProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRenJoinProj{time}")));
      await session.step(93, "And the \"BDDRenJoinProj{time}\" project on the server should hold the tables \"BDDRenJoinQ{time}, entity_types, result\"", () => projectHoldsTables(page, session.text("BDDRenJoinProj{time}"), session.text("BDDRenJoinQ{time}, entity_types, result")));
      await session.step(94, "And the \"result\" table of the \"BDDRenJoinProj{time}\" project should be saved with data sync", () => savedWithDataSync(page, "result", session.text("BDDRenJoinProj{time}")));
      await session.step(95, "And the creation script of the \"result\" table of the \"BDDRenJoinProj{time}\" project on the server should contain 'JoinTables(\"entity_types\", \"BDDRenJoinQ{time}\"'", () => creationScriptHolds(page, "result", session.text("BDDRenJoinProj{time}"), session.text("JoinTables(\"entity_types\", \"BDDRenJoinQ{time}\"")));
      await session.step(96, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
      await session.step(97, "When user closes all views", () => closeAllViews(page));
      await session.step(98, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(99, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(100, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(101, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(102, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(105, "When user collapses Databases---Postgres---Datagrok---BDDRenJoinQ{time} tree node inside browse tree", () => collapse(page, el(session.text("Databases---Postgres---Datagrok---BDDRenJoinQ{time} tree node inside browse tree"))));
      await session.step(106, "And user picks \"Rename...\" from the context menu of Databases---Postgres---Datagrok---BDDRenJoinQ{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Databases---Postgres---Datagrok---BDDRenJoinQ{time} tree node inside browse tree"))));
      await session.step(107, "Then \"Rename dataquery\" dialog should be visible", () => shouldBe(page, el("\"Rename dataquery\" dialog"), "visible"));
      await session.step(108, "When user enters \"BDDRenJoinQR{time}\" into Name input in \"Rename dataquery\" dialog", () => enterInto(page, session.text("BDDRenJoinQR{time}"), el("Name input in \"Rename dataquery\" dialog")));
      await session.step(109, "And user clicks on OK button in \"Rename dataquery\" dialog", () => clickOn(page, el("OK button in \"Rename dataquery\" dialog")));
      await session.step(110, "Then the \"Rename dataquery\" dialog should close", () => dialogCloses(page, "Rename dataquery"));
      await session.step(111, "And 1 query named \"BDDRenJoinQR{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDRenJoinQR{time}")));
      await session.step(112, "And 0 queries named \"BDDRenJoinQ{time}\" should be on the server", () => queriesOnServer(page, 0, session.text("BDDRenJoinQ{time}")));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
      await session.step(114, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("After the rename, the join's creation script names the new query while the table keeps the old name", async () => {
      await session.step(117, "Then the creation script of the \"BDDRenJoinQ{time}\" table of the \"BDDRenJoinProj{time}\" project on the server should contain ':BDDRenJoinQR{time}()'", () => creationScriptHolds(page, session.text("BDDRenJoinQ{time}"), session.text("BDDRenJoinProj{time}"), session.text(":BDDRenJoinQR{time}()")));
      await session.step(118, "And the creation script of the \"result\" table of the \"BDDRenJoinProj{time}\" project on the server should contain 'JoinTables(\"entity_types\", \"BDDRenJoinQR{time}\"'", () => creationScriptHolds(page, "result", session.text("BDDRenJoinProj{time}"), session.text("JoinTables(\"entity_types\", \"BDDRenJoinQR{time}\"")));
      await session.step(119, "And the \"BDDRenJoinProj{time}\" project on the server should hold the tables \"BDDRenJoinQ{time}, entity_types, result\"", () => projectHoldsTables(page, session.text("BDDRenJoinProj{time}"), session.text("BDDRenJoinQ{time}, entity_types, result")));
    });
    await run.scenario("Reopened from the Dashboards gallery, the query's table and entity_types are re-read", async () => {
      await session.step(122, "When user closes all views", () => closeAllViews(page));
      await session.step(123, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(124, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(125, "And user enters \"BDDRenJoinProj{time}\" into gallery search", () => enterInto(page, session.text("BDDRenJoinProj{time}"), el("gallery search")));
      await session.step(126, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(127, "And user double-clicks on BDDRenJoinProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDRenJoinProj{time} gallery card"))));
      await session.step(129, "Then table \"BDDRenJoinQ{time}\" should have been reloaded by data sync with the \"rename join query\" row count", () => tableReloadedAsRemembered(page, session.text("BDDRenJoinQ{time}"), "rename join query"));
      await session.step(130, "And table \"entity_types\" should have been reloaded by data sync with the \"rename join query\" row count", () => tableReloadedAsRemembered(page, "entity_types", "rename join query"));
    });
    await run.scenario("The reopened project holds the join over the renamed query's result (GROK-21026)", async () => {
      await session.step(134, "Then table \"result\" should have been reloaded by data sync with the \"rename join query\" row count", () => tableReloadedAsRemembered(page, "result", "rename join query"));
      await session.step(135, "And the table views \"BDDRenJoinQ{time}, entity_types, result\" should be open", () => tableViewsOpen(page, session.text("BDDRenJoinQ{time}, entity_types, result")));
    }, {knownFailure: true});
    await run.scenario("The page is reloaded, dropping the open that never completed", async () => {
      await session.step(138, "When user reloads the page", () => reloadPage(page));
      await session.step(139, "Then the table view tabs should read \"\"", () => tableViewTabs(page, ""));
    });
    run.finish();
  });
});
