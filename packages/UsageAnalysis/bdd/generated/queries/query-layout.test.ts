/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/query-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.queries]
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
import {clickOn, enterInto, isExpanded, replaceCode, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, currentViewType, dialogCloses, entityHasLayout, layoutsDeleted, noProjectOnServer, noQueryOnServer, openProjectWithTable, projectsOnServer, queriesOnServer, refreshBrowse, standHasReachableConnection, toolboxPaneHidden, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu, pointerAway, readingIs, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A query's layout", () => {
  const session = feature(test, "features/queries/query-layout.feature", import.meta.url);
  test("A query's layout", {tag: ["@journey", "@serial", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the stand has a reachable \"PostgresTest\" connection", () => standHasReachableConnection(page, "PostgresTest"));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And no query named \"BDD-Q-layout-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Q-layout-{time}")));
    await session.step(25, "And no project named \"BDD-Q-layout-proj-{time}\" is on the server", () => noProjectOnServer(page, session.text("BDD-Q-layout-proj-{time}")));
    await session.step(26, "And the layouts named \"BDD-Q-layout-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("BDD-Q-layout-{time}")));
    await run.scenario("The Layout tab waits for a run, then takes viewers from the toolbox", async () => {
      await session.step(29, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(30, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(31, "When user picks \"New Query...\" from the context menu of Databases---Postgres---NorthwindTest tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(32, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(33, "When user enters \"BDD-Q-layout-{time}\" into Name input", () => enterInto(page, session.text("BDD-Q-layout-{time}"), el("Name input")));
      await session.step(34, "And user replaces the code of code editor with \"select * from products\"", () => replaceCode(page, el("code editor"), "select * from products"));
      await session.step(35, "And user clicks on Layout tab", () => clickOn(page, el("Layout tab")));
      await session.step(36, "Then \"Run query to get data and edit layout\" text should be visible", () => shouldBe(page, el("\"Run query to get data and edit layout\" text"), "visible"));
      await session.step(37, "When user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(38, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(39, "And the \"rows\" reading of grid should be 77", () => readingIs(page, "rows", el("grid"), 77));
      await session.step(40, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(41, "When user moves the pointer away from play icon", () => pointerAway(page, el("play icon")));
      await session.step(42, "And user clicks on scatter plot icon in toolbox", () => clickOn(page, el("scatter plot icon in toolbox")));
      await session.step(43, "And user clicks on correlation plot icon in toolbox", () => clickOn(page, el("correlation plot icon in toolbox")));
      await session.step(44, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(45, "And correlation plot viewer should be visible", () => shouldBe(page, el("correlation plot viewer"), "visible"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Save keeps the query with its layout, and the saved query runs with its viewers", async () => {
      await session.step(50, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(51, "Then 1 query named \"BDD-Q-layout-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-Q-layout-{time}")));
      await session.step(52, "And the query \"BDD-Q-layout-{time}\" on the server should have a layout", () => entityHasLayout(page, "query", session.text("BDD-Q-layout-{time}")));
      await session.step(53, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
      await session.step(54, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(55, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(56, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(57, "And user picks \"Run\" from the context menu of Databases---Postgres---NorthwindTest---BDD-Q-layout-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Run", el(session.text("Databases---Postgres---NorthwindTest---BDD-Q-layout-{time} tree node inside browse tree"))));
      await session.step(58, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(59, "And the table should have 77 rows", () => rowCount(page, 77));
      await session.step(60, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(61, "And the open tableview should have 1 correlation plot viewer", () => viewerCount(page, 1, "correlation plot"));
    });
    await run.scenario("The result with one more viewer is saved as a project and opens with all three", async () => {
      await session.step(64, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(65, "When user clicks on histogram icon in toolbox", () => clickOn(page, el("histogram icon in toolbox")));
      await session.step(66, "Then the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
      await session.step(67, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(68, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(69, "When user enters \"BDD-Q-layout-proj-{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDD-Q-layout-proj-{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(70, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(71, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(72, "And 1 project named \"BDD-Q-layout-proj-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-Q-layout-proj-{time}")));
      await session.step(73, "And \"Share BDD-Q-layout-proj-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-Q-layout-proj-{time}\" dialog")), "visible"));
      await session.step(74, "When user clicks on CANCEL button in \"Share BDD-Q-layout-proj-{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDD-Q-layout-proj-{time}\" dialog"))));
      await session.step(75, "Then the \"Share BDD-Q-layout-proj-{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-Q-layout-proj-{time}")));
      await session.step(76, "When user closes all views", () => closeAllViews(page));
      await session.step(77, "And user opens the \"BDD-Q-layout-proj-{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDD-Q-layout-proj-{time}")));
      await session.step(78, "Then the table should have 77 rows", () => rowCount(page, 77));
      await session.step(79, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(80, "And the open tableview should have 1 correlation plot viewer", () => viewerCount(page, 1, "correlation plot"));
      await session.step(81, "And the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
      await session.step(83, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("REFRESH of the reopened result keeps its viewers", async () => {
      await session.step(86, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(87, "When user clicks on REFRESH button in toolbox and the query it runs completes", () => clickRunsQuery(page, el("REFRESH button in toolbox")));
      await session.step(88, "Then the table should have 77 rows", () => rowCount(page, 77));
      await session.step(89, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(90, "And the open tableview should have 1 correlation plot viewer", () => viewerCount(page, 1, "correlation plot"));
      await session.step(91, "And the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
      await session.step(93, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
