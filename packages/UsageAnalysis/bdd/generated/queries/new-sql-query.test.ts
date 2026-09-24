/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/new-sql-query.feature
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, holdsCode, isExpanded, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeCurrentView, currentViewType, noQueryOnServer, queriesOnServer, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A SQL query started from a table of the schema", () => {
  const session = feature(test, "features/queries/new-sql-query.feature", import.meta.url);
  test("The editor opens on the table's select, runs it and saves it", {tag: ["@serial", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And no query named \"BDD-Q-sql-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Q-sql-{time}")));
    await session.step(21, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(22, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(23, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(26, "Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---NorthwindTest---Orders tree node inside browse tree"), "visible"));
    await session.step(27, "Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree")));
    await session.step(28, "And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(29, "When user picks \"New SQL Query...\" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---products tree node inside browse tree", () => pickFromContextMenu(page, "New SQL Query...", el("Databases---Postgres---NorthwindTest---Schemas---public---products tree node inside browse tree")));
    await session.step(30, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
    await session.step(31, "And Name input should have value \"products\"", () => shouldHaveValue(page, el("Name input"), "products"));
    await session.step(32, "And code editor should hold the code \"select * from public.products\"", () => holdsCode(page, el("code editor"), "select * from public.products"));
    await session.step(33, "When user clicks on play icon", () => clickOn(page, el("play icon")));
    await session.step(34, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    await session.step(35, "And the \"rows\" reading of grid should be 77", () => readingIs(page, "rows", el("grid"), 77));
    await session.step(36, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(37, "When user clicks on \"Run query...\" action in toolbox", () => clickOn(page, el("\"Run query...\" action in toolbox")));
    await session.step(38, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(39, "And the table should have 77 rows", () => rowCount(page, 77));
    await session.step(40, "When user closes the current view", () => closeCurrentView(page));
    await session.step(41, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
    await session.step(42, "When user enters \"BDD-Q-sql-{time}\" into Name input", () => enterInto(page, session.text("BDD-Q-sql-{time}"), el("Name input")));
    await session.step(43, "And user clicks on Save button", () => clickOn(page, el("Save button")));
    await session.step(44, "Then 1 query named \"BDD-Q-sql-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-Q-sql-{time}")));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
