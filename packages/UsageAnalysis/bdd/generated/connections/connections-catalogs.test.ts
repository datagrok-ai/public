/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-catalogs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {catalogCommentIs, catalogHasNoComment} from '../../bindings/connections.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, followingShouldBe, isExpanded, rightClickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The catalogs of an MS SQL connection", () => {
  const session = feature(test, "features/connections/connections-catalogs.feature", import.meta.url);
  test("The catalogs of an MS SQL connection", {tag: ["@connections", "@full-stand", "@serial", "@journey", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(29, "And Databases---MS-SQL tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---MS-SQL tree node inside browse tree")));
    await session.step(30, "And Databases---MS-SQL---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---MS-SQL---NorthwindTest tree node inside browse tree")));
    await run.scenario("The Catalogs group lists the connection's databases, and a catalog its tables", async () => {
      await session.step(33, "Given Databases---MS-SQL---NorthwindTest---Catalogs tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---MS-SQL---NorthwindTest---Catalogs tree node inside browse tree")));
      await session.step(34, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Databases---MS-SQL---NorthwindTest---Catalogs---northwind tree node inside browse tree"],["Databases---MS-SQL---NorthwindTest---Catalogs---tempdb tree node inside browse tree"],["Databases---MS-SQL---NorthwindTest---Catalogs---master tree node inside browse tree"]]), [["Databases---MS-SQL---NorthwindTest---Catalogs---northwind tree node inside browse tree"],["Databases---MS-SQL---NorthwindTest---Catalogs---tempdb tree node inside browse tree"],["Databases---MS-SQL---NorthwindTest---Catalogs---master tree node inside browse tree"]]);
      await session.step(38, "Given Databases---MS-SQL---NorthwindTest---Catalogs---northwind tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---MS-SQL---NorthwindTest---Catalogs---northwind tree node inside browse tree")));
      await session.step(39, "And Databases---MS-SQL---NorthwindTest---Catalogs---northwind---dbo tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---MS-SQL---NorthwindTest---Catalogs---northwind---dbo tree node inside browse tree")));
      await session.step(40, "Then Databases---MS-SQL---NorthwindTest---Catalogs---northwind---dbo---orders tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---MS-SQL---NorthwindTest---Catalogs---northwind---dbo---orders tree node inside browse tree"), "visible"));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A catalog's menu opens it as a table", async () => {
      await session.step(44, "When user right-clicks on Databases---MS-SQL---NorthwindTest---Catalogs---northwind tree node inside browse tree", () => rightClickOn(page, el("Databases---MS-SQL---NorthwindTest---Catalogs---northwind tree node inside browse tree")));
      await session.step(45, "Then the open menu should list \"Browse\"", () => menuLists(page, "Browse"));
      await session.step(46, "And the open menu should list \"Open as table\"", () => menuLists(page, "Open as table"));
      await session.step(47, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(48, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A comment on a catalog is saved on the server", async () => {
      await session.step(51, "Given the \"tempdb\" catalog of the \"MSSQLTest\" connection has no comment", () => catalogHasNoComment(page, "tempdb", "MSSQLTest"));
      await session.step(52, "When user clicks on Databases---MS-SQL---NorthwindTest---Catalogs---tempdb tree node inside browse tree", () => clickOn(page, el("Databases---MS-SQL---NorthwindTest---Catalogs---tempdb tree node inside browse tree")));
      await session.step(53, "Then the context panel should show \"tempdb\"", () => contextPanelShows(page, "tempdb"));
      await session.step(54, "Given \"Database meta\" section in context panel is expanded", () => isExpanded(page, el("\"Database meta\" section in context panel")));
      await session.step(55, "When user enters \"BDD comment {run}\" into Comment input in context panel", () => enterInto(page, session.text("BDD comment {run}"), el("Comment input in context panel")));
      await session.step(56, "And user clicks on SAVE button in context panel", () => clickOn(page, el("SAVE button in context panel")));
      await session.step(57, "Then the \"tempdb\" catalog of the \"MSSQLTest\" connection should have the comment \"BDD comment {run}\"", () => catalogCommentIs(page, "tempdb", "MSSQLTest", session.text("BDD comment {run}")));
      await session.step(58, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Clearing the comment in the pane clears it on the server", async () => {
      await session.step(65, "When user clears Comment input in context panel", () => clearField(page, el("Comment input in context panel")));
      await session.step(66, "And user clicks on SAVE button in context panel", () => clickOn(page, el("SAVE button in context panel")));
      await session.step(67, "Then the \"tempdb\" catalog of the \"MSSQLTest\" connection should have the comment \"\"", () => catalogCommentIs(page, "tempdb", "MSSQLTest", ""));
      await session.step(68, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true});
    run.finish();
  });
});
