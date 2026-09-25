/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/query-post-process.feature
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
import {appendToEditor, clickOn, enterInto, isExpanded, replaceCode, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, currentViewType, noQueryOnServer, queriesOnServer, queryPostProcess, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu, pointerAway, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A query's post-process runs on its result", () => {
  const session = feature(test, "features/queries/query-post-process.feature", import.meta.url);
  test("A query's post-process runs on its result", {tag: ["@journey", "@serial", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And no query named \"BDD-Q-pp-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Q-pp-{time}")));
    await run.scenario("A line is typed into the Post-Process tab of a new query", async () => {
      await session.step(20, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(21, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(22, "When user picks \"New Query...\" from the context menu of Databases---Postgres---NorthwindTest tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(23, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(24, "When user enters \"BDD-Q-pp-{time}\" into Name input", () => enterInto(page, session.text("BDD-Q-pp-{time}"), el("Name input")));
      await session.step(25, "And user replaces the code of code editor with \"select * from products\"", () => replaceCode(page, el("code editor"), "select * from products"));
      await session.step(26, "And user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(27, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(28, "And the \"rows\" reading of grid should be 77", () => readingIs(page, "rows", el("grid"), 77));
      await session.step(29, "When user moves the pointer away from play icon", () => pointerAway(page, el("play icon")));
      await session.step(30, "And user clicks on Post-Process tab", () => clickOn(page, el("Post-Process tab")));
      await session.step(33, "And user appends \"grok.shell.info('PP' + result.rowCount);\" to code editor", () => appendToEditor(page, "grok.shell.info('PP' + result.rowCount);", el("code editor")));
      await session.step(34, "Then no errors should have been logged", () => noErrors(page));
      await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Save right after typing keeps the line", async () => {
      await session.step(40, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(41, "Then 1 query named \"BDD-Q-pp-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-Q-pp-{time}")));
      await session.step(42, "And the query \"BDD-Q-pp-{time}\" on the server should have a post-process containing \"grok.shell.info('PP' + result.rowCount);\"", () => queryPostProcess(page, session.text("BDD-Q-pp-{time}"), "grok.shell.info('PP' + result.rowCount);"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The saved query announces its row count when it is run from the tree", async () => {
      await session.step(46, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
      await session.step(47, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(48, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(49, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(50, "And user picks \"Run\" from the context menu of Databases---Postgres---NorthwindTest---BDD-Q-pp-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Run", el(session.text("Databases---Postgres---NorthwindTest---BDD-Q-pp-{time} tree node inside browse tree"))));
      await session.step(51, "Then the table should have 77 rows", () => rowCount(page, 77));
      await session.step(52, "And an info balloon containing \"PP77\" should have been shown", () => infoBalloonText(page, "PP77"));
    });
    run.finish();
  });
});
