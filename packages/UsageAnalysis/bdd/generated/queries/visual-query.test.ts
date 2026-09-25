/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/visual-query.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.queries]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/nx.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addToBuilderRow, builderRowHolds} from '../../bindings/queries.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, currentViewType, noQueryOnServer, queriesOnServer, toolboxPaneHidden, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A visual query built on a table", () => {
  const session = feature(test, "features/queries/visual-query.feature", import.meta.url);
  test("A visual query built on a table", {tag: ["@journey", "@serial", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And no query named \"BDD-Q-vq-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Q-vq-{time}")));
    await run.scenario("The builder groups and aggregates a table", async () => {
      await session.step(23, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(24, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(25, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(26, "Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---NorthwindTest---Orders tree node inside browse tree"), "visible"));
      await session.step(27, "Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree")));
      await session.step(28, "And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
      await session.step(29, "When user picks \"New Visual Query...\" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---customers tree node inside browse tree", () => pickFromContextMenu(page, "New Visual Query...", el("Databases---Postgres---NorthwindTest---Schemas---public---customers tree node inside browse tree")));
      await session.step(30, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(31, "When user adds \"companyname\" to the \"Group-by\" row of the visual query", () => addToBuilderRow(page, "companyname", "Group-by"));
      await session.step(32, "Then the \"Group-by\" row of the visual query should hold \"companyname\"", () => builderRowHolds(page, "Group-by", "companyname"));
      await session.step(33, "When user adds \"customerid\" to the \"Aggregate\" row of the visual query", () => addToBuilderRow(page, "customerid", "Aggregate"));
      await session.step(34, "Then the \"Aggregate\" row of the visual query should hold \"values(customerid)\"", () => builderRowHolds(page, "Aggregate", "values(customerid)"));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
      await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The visual query is saved and runs as an ordinary query", async () => {
      await session.step(39, "When user enters \"BDD-Q-vq-{time}\" into Name input", () => enterInto(page, session.text("BDD-Q-vq-{time}"), el("Name input")));
      await session.step(40, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(41, "Then 1 query named \"BDD-Q-vq-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-Q-vq-{time}")));
      await session.step(42, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(43, "When user clicks on \"Run query...\" action in toolbox", () => clickOn(page, el("\"Run query...\" action in toolbox")));
      await session.step(44, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(45, "And the table should have 91 rows", () => rowCount(page, 91));
      await session.step(46, "And the table should have a column \"companyname\"", () => hasColumn(page, "companyname"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
      await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Reopening the saved query shows the builder with its rows", async () => {
      await session.step(53, "When user closes all views", () => closeAllViews(page));
      await session.step(54, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
      await session.step(55, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(56, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(57, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(58, "And user picks \"Edit...\" from the context menu of Databases---Postgres---NorthwindTest---BDD-Q-vq-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---NorthwindTest---BDD-Q-vq-{time} tree node inside browse tree"))));
      await session.step(59, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(60, "And the \"Group-by\" row of the visual query should hold \"companyname\"", () => builderRowHolds(page, "Group-by", "companyname"));
      await session.step(61, "And the \"Aggregate\" row of the visual query should hold \"values(customerid)\"", () => builderRowHolds(page, "Aggregate", "values(customerid)"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
