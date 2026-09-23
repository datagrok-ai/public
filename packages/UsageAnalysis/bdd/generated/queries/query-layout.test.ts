/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/query-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.queries]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, isExpanded, replaceCode, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, currentViewType, noProjectOnServer, noQueryOnServer, queriesOnServer, queryHasLayout, toolboxPaneHidden, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu, pointerAway, readingIs, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A query's layout", () => {
  const session = feature(test, "features/queries/query-layout.feature", import.meta.url);
  test("A query's layout", {tag: ["@journey", "@serial", "@realizes:views.queries", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And no query named \"BDD-Q-layout-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Q-layout-{time}")));
    await session.step(21, "And no project named \"BDD-Q-layout-{time}\" is on the server", () => noProjectOnServer(page, session.text("BDD-Q-layout-{time}")));
    await run.scenario("The Layout tab waits for a run, then takes viewers from the toolbox", async () => {
      await session.step(24, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(25, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(26, "When user picks \"New Query...\" from the context menu of Databases---Postgres---NorthwindTest tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(27, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(28, "When user enters \"BDD-Q-layout-{time}\" into Name input", () => enterInto(page, session.text("BDD-Q-layout-{time}"), el("Name input")));
      await session.step(29, "And user replaces the code of code editor with \"select * from products\"", () => replaceCode(page, el("code editor"), "select * from products"));
      await session.step(30, "And user clicks on Layout tab", () => clickOn(page, el("Layout tab")));
      await session.step(31, "Then \"Run query to get data and edit layout\" text should be visible", () => shouldBe(page, el("\"Run query to get data and edit layout\" text"), "visible"));
      await session.step(32, "When user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(33, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(34, "And the \"rows\" reading of grid should be 77", () => readingIs(page, "rows", el("grid"), 77));
      await session.step(35, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(36, "When user moves the pointer away from play icon", () => pointerAway(page, el("play icon")));
      await session.step(37, "And user clicks on scatter plot icon in toolbox", () => clickOn(page, el("scatter plot icon in toolbox")));
      await session.step(38, "And user clicks on correlation plot icon in toolbox", () => clickOn(page, el("correlation plot icon in toolbox")));
      await session.step(39, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(40, "And correlation plot viewer should be visible", () => shouldBe(page, el("correlation plot viewer"), "visible"));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Save keeps the query with its layout, and the saved query runs with its viewers", async () => {
      await session.step(50, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(51, "Then 1 query named \"BDD-Q-layout-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-Q-layout-{time}")));
      await session.step(52, "And the query \"BDD-Q-layout-{time}\" on the server should have a layout", () => queryHasLayout(page, session.text("BDD-Q-layout-{time}")));
      await session.step(53, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
      await session.step(54, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(55, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(56, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(58, "And user hovers over Databases---Postgres---NorthwindTest---BDD-Q-layout-{time} tree node inside browse tree", () => hoverOver(page, el(session.text("Databases---Postgres---NorthwindTest---BDD-Q-layout-{time} tree node inside browse tree"))));
      await session.step(59, "And user picks \"Run\" from the context menu of Databases---Postgres---NorthwindTest---BDD-Q-layout-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Run", el(session.text("Databases---Postgres---NorthwindTest---BDD-Q-layout-{time} tree node inside browse tree"))));
      await session.step(60, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(61, "And the table should have 77 rows", () => rowCount(page, 77));
      await session.step(62, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(63, "And the open tableview should have 1 correlation plot viewer", () => viewerCount(page, 1, "correlation plot"));
    }, {knownFailure: true});
    run.finish();
  });
});
