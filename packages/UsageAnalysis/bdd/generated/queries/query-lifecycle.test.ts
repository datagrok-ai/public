/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/query-lifecycle.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.queries]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, followingShouldBe, holdsCode, isExpanded, replaceCode, shouldBe, shouldContainText, shouldHaveValue, textAreaHolds} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeCurrentView, contextPanelOpen, contextPanelShows, currentViewType, dialogCloses, noQueryOnServer, queriesOnServer, refreshBrowse, standHasReachableConnection, toolboxPaneHidden, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A SQL query from creation to deletion", () => {
  const session = feature(test, "features/queries/query-lifecycle.feature", import.meta.url);
  test("A SQL query from creation to deletion", {tag: ["@journey", "@serial", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the stand has a reachable \"PostgresTest\" connection", () => standHasReachableConnection(page, "PostgresTest"));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And no query named \"BDD-Q-life-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Q-life-{time}")));
    await session.step(25, "And no query named \"BDD-Q-life-renamed-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Q-life-renamed-{time}")));
    await run.scenario("A new query is typed, run in its editor and on its own, and saved", async () => {
      await session.step(28, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(29, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(30, "When user picks \"New Query...\" from the context menu of Databases---Postgres---NorthwindTest tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(31, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(32, "When user enters \"BDD-Q-life-{time}\" into Name input", () => enterInto(page, session.text("BDD-Q-life-{time}"), el("Name input")));
      await session.step(33, "And user replaces the code of code editor with \"select * from products\"", () => replaceCode(page, el("code editor"), "select * from products"));
      await session.step(34, "And user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(35, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(36, "And the \"rows\" reading of grid should be 77", () => readingIs(page, "rows", el("grid"), 77));
      await session.step(37, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(38, "When user clicks on \"Run query...\" action in toolbox and the query it runs completes", () => clickRunsQuery(page, el("\"Run query...\" action in toolbox")));
      await session.step(39, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(40, "And the table should have 77 rows", () => rowCount(page, 77));
      await session.step(41, "When user closes the current view", () => closeCurrentView(page));
      await session.step(42, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(43, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(44, "Then 1 query named \"BDD-Q-life-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-Q-life-{time}")));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
      await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Edit renames the query and changes its SQL", async () => {
      await session.step(49, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
      await session.step(50, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(51, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(52, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(53, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(54, "When user picks \"Edit...\" from the context menu of Databases---Postgres---NorthwindTest---BDD-Q-life-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---NorthwindTest---BDD-Q-life-{time} tree node inside browse tree"))));
      await session.step(55, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(56, "And Name input should have value \"BDD-Q-life-{time}\"", () => shouldHaveValue(page, el("Name input"), session.text("BDD-Q-life-{time}")));
      await session.step(57, "And code editor should hold the code \"select * from products\"", () => holdsCode(page, el("code editor"), "select * from products"));
      await session.step(58, "When user enters \"BDD-Q-life-renamed-{time}\" into Name input", () => enterInto(page, session.text("BDD-Q-life-renamed-{time}"), el("Name input")));
      await session.step(59, "And user replaces the code of code editor with \"select * from orders\"", () => replaceCode(page, el("code editor"), "select * from orders"));
      await session.step(60, "And user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(61, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(62, "And the \"rows\" reading of grid should be 830", () => readingIs(page, "rows", el("grid"), 830));
      await session.step(63, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(64, "Then 1 query named \"BDD-Q-life-renamed-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-Q-life-renamed-{time}")));
      await session.step(65, "And 0 queries named \"BDD-Q-life-{time}\" should be on the server", () => queriesOnServer(page, 0, session.text("BDD-Q-life-{time}")));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The connection's queries view finds the query and its panes describe it", async () => {
      await session.step(70, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
      await session.step(71, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(72, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(73, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(74, "When user clicks on first Databases---Postgres---NorthwindTest tree node inside browse tree", () => clickOn(page, el("first Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(75, "Then the current view should be a queries view", () => currentViewType(page, "queries"));
      await session.step(76, "When user enters \"BDD-Q-life-renamed-{time}\" into gallery search", () => enterInto(page, session.text("BDD-Q-life-renamed-{time}"), el("gallery search")));
      await session.step(77, "Then \"BDD-Q-life-renamed-{time}\" gallery card should be visible", () => shouldBe(page, el(session.text("\"BDD-Q-life-renamed-{time}\" gallery card")), "visible"));
      await session.step(78, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(79, "When user clicks on \"BDD-Q-life-renamed-{time}\" gallery card", () => clickOn(page, el(session.text("\"BDD-Q-life-renamed-{time}\" gallery card"))));
      await session.step(80, "Then the context panel should show \"BDD-Q-life-renamed-{time}\"", () => contextPanelShows(page, session.text("BDD-Q-life-renamed-{time}")));
      await session.step(81, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Details section in context panel"],["Run section in context panel"],["Query section in context panel"],["Transformations section in context panel"],["Sharing section in context panel"],["Chats section in context panel"]]), [["Details section in context panel"],["Run section in context panel"],["Query section in context panel"],["Transformations section in context panel"],["Sharing section in context panel"],["Chats section in context panel"]]);
      await session.step(88, "Given Query section in context panel is expanded", () => isExpanded(page, el("Query section in context panel")));
      await session.step(89, "Then the text area of Query pane in context panel should hold \"select * from orders\"", () => textAreaHolds(page, el("Query pane in context panel"), "select * from orders"));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
      await session.step(91, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Delete asks first, then removes the query from the server and the tree", async () => {
      await session.step(94, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
      await session.step(95, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(96, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(97, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(98, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(99, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(100, "And user picks \"Delete\" from the context menu of Databases---Postgres---NorthwindTest---BDD-Q-life-renamed-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete", el(session.text("Databases---Postgres---NorthwindTest---BDD-Q-life-renamed-{time} tree node inside browse tree"))));
      await session.step(101, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(102, "And \"Are you sure?\" dialog should contain the text \"BDD-Q-life-renamed-{time}\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("BDD-Q-life-renamed-{time}")));
      await session.step(103, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(104, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(105, "And 0 queries named \"BDD-Q-life-renamed-{time}\" should be on the server", () => queriesOnServer(page, 0, session.text("BDD-Q-life-renamed-{time}")));
      await session.step(106, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(107, "Then Databases---Postgres---NorthwindTest---BDD-Q-life-renamed-{time} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Databases---Postgres---NorthwindTest---BDD-Q-life-renamed-{time} tree node inside browse tree")), "absent"));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
      await session.step(109, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
