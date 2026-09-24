/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/query-transformations.feature
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
import {openAction} from '../../bindings/queries.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, isExpanded, replaceCode, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeCurrentView, currentViewType, dialogCloses, noQueryOnServer, queriesOnServer, queryNoTransformations, queryTransformations, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu, pointerAway, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Transformations saved with a query", () => {
  const session = feature(test, "features/queries/query-transformations.feature", import.meta.url);
  test("Transformations saved with a query", {tag: ["@journey", "@serial", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And no query named \"BDD-Q-tr-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Q-tr-{time}")));
    await run.scenario("A column added in the Transformations tab is in the result of the saved query", async () => {
      await session.step(23, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(24, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(25, "When user picks \"New Query...\" from the context menu of Databases---Postgres---NorthwindTest tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(26, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(27, "When user enters \"BDD-Q-tr-{time}\" into Name input", () => enterInto(page, session.text("BDD-Q-tr-{time}"), el("Name input")));
      await session.step(28, "And user replaces the code of code editor with \"select * from products\"", () => replaceCode(page, el("code editor"), "select * from products"));
      await session.step(29, "And user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(30, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(31, "And the \"rows\" reading of grid should be 77", () => readingIs(page, "rows", el("grid"), 77));
      await session.step(32, "When user moves the pointer away from play icon", () => pointerAway(page, el("play icon")));
      await session.step(33, "And user clicks on Transformations tab", () => clickOn(page, el("Transformations tab")));
      await session.step(34, "And user opens the \"Add New Column\" action of the transformations browser", () => openAction(page, "Add New Column"));
      await session.step(35, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(36, "When user enters \"doubled\" into Name input in \"Add New Column\" dialog", () => enterInto(page, "doubled", el("Name input in \"Add New Column\" dialog")));
      await session.step(37, "And user replaces the code of code editor in \"Add New Column\" dialog with \"${productid} * 2\"", () => replaceCode(page, el("code editor in \"Add New Column\" dialog"), "${productid} * 2"));
      await session.step(38, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(39, "Then the \"Add New Column\" dialog should close", () => dialogCloses(page, "Add New Column"));
      await session.step(40, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(41, "Then 1 query named \"BDD-Q-tr-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-Q-tr-{time}")));
      await session.step(44, "And the query \"BDD-Q-tr-{time}\" on the server should have transformations containing \"doubled\"", () => queryTransformations(page, session.text("BDD-Q-tr-{time}"), "doubled"));
      await session.step(45, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(46, "When user clicks on \"Run query...\" action in toolbox", () => clickOn(page, el("\"Run query...\" action in toolbox")));
      await session.step(47, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(48, "And the table should have 77 rows", () => rowCount(page, 77));
      await session.step(49, "And the table should have a column \"doubled\"", () => hasColumn(page, "doubled"));
      await session.step(50, "And the value of \"doubled\" column in row 1 should be \"2\"", () => valueInRow(page, "doubled", 1, "2"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
      await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Removing the step takes the column out of the saved query's result", async () => {
      await session.step(55, "When user closes the current view", () => closeCurrentView(page));
      await session.step(56, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(57, "When user clicks on Transformations tab", () => clickOn(page, el("Transformations tab")));
      await session.step(59, "And user hovers over last \"Remove step\" icon", () => hoverOver(page, el("last \"Remove step\" icon")));
      await session.step(60, "And user clicks on last \"Remove step\" icon", () => clickOn(page, el("last \"Remove step\" icon")));
      await session.step(61, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(62, "Then the query \"BDD-Q-tr-{time}\" on the server should not have transformations containing \"doubled\"", () => queryNoTransformations(page, session.text("BDD-Q-tr-{time}"), "doubled"));
      await session.step(63, "When user clicks on \"Run query...\" action in toolbox", () => clickOn(page, el("\"Run query...\" action in toolbox")));
      await session.step(64, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(65, "And the table should have 77 rows", () => rowCount(page, 77));
      await session.step(66, "And the table should not have a column \"doubled\"", () => hasNoColumn(page, "doubled"));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
      await session.step(68, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
