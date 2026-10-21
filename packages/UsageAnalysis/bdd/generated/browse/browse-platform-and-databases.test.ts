/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-platform-and-databases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
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
import {clickOn, collapse, followingShouldBe, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, packageInstalled, standHasReachableConnection, standRunsService, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuDoesNotList, menuLists, noBalloons, noErrors, openContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Platform and Databases sections of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-platform-and-databases.feature", import.meta.url);
  test("The Platform section lists what an administrator manages", {tag: ["@browse", "@realizes:views.browse", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(45, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(46, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Platform---Plugins tree node inside browse tree"],["Platform---Credentials tree node inside browse tree"],["Platform---Functions tree node inside browse tree"],["Platform---Users tree node inside browse tree"],["Platform---Groups tree node inside browse tree"],["Platform---Roles tree node inside browse tree"],["Platform---Predictive-models tree node inside browse tree"],["Platform---Dockers tree node inside browse tree"]]), [["Platform---Plugins tree node inside browse tree"],["Platform---Credentials tree node inside browse tree"],["Platform---Functions tree node inside browse tree"],["Platform---Users tree node inside browse tree"],["Platform---Groups tree node inside browse tree"],["Platform---Roles tree node inside browse tree"],["Platform---Predictive-models tree node inside browse tree"],["Platform---Dockers tree node inside browse tree"]]);
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
    await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A Platform node opens the view it is named after", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(59, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(60, "When user clicks on Platform---Users tree node inside browse tree", () => clickOn(page, el("Platform---Users tree node inside browse tree")));
    await session.step(61, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
    await session.step(62, "When user clicks on Platform---Groups tree node inside browse tree", () => clickOn(page, el("Platform---Groups tree node inside browse tree")));
    await session.step(63, "Then the \"Groups\" view should be current", () => viewIsCurrent(page, "Groups"));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
    await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Databases section lists the connected providers", {tag: ["@browse", "@realizes:views.browse", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(69, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(70, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Databases---Postgres tree node inside browse tree"],["Databases---MySQL tree node inside browse tree"],["Databases---Oracle tree node inside browse tree"],["Databases---MariaDB tree node inside browse tree"],["Databases---MS-SQL tree node inside browse tree"]]), [["Databases---Postgres tree node inside browse tree"],["Databases---MySQL tree node inside browse tree"],["Databases---Oracle tree node inside browse tree"],["Databases---MariaDB tree node inside browse tree"],["Databases---MS-SQL tree node inside browse tree"]]);
    await session.step(76, "And no errors should have been logged", () => noErrors(page));
    await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A provider opens down to its connections", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(80, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(81, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(82, "Then Databases---Postgres---Datagrok tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---Datagrok tree node inside browse tree"), "visible"));
    await session.step(83, "And Databases---Postgres---CHEMBL tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---CHEMBL tree node inside browse tree"), "visible"));
    await session.step(84, "When user collapses Databases---Postgres tree node inside browse tree", () => collapse(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(85, "Then Databases---Postgres---Datagrok tree node inside browse tree should be hidden", () => shouldBe(page, el("Databases---Postgres---Datagrok tree node inside browse tree"), "hidden"));
    await session.step(86, "And Databases---Postgres tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres tree node inside browse tree"), "visible"));
    await session.step(87, "And no errors should have been logged", () => noErrors(page));
    await session.step(88, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A saved query opens from the tree with its details", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(93, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(94, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(95, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(96, "And the stand has a reachable \"PostgresTest\" connection", () => standHasReachableConnection(page, "PostgresTest"));
    await session.step(97, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(98, "When user clicks on Databases---Postgres---NorthwindTest---Orders tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Orders tree node inside browse tree")));
    await session.step(99, "Then the context panel should show \"Orders\"", () => contextPanelShows(page, "Orders"));
    await session.step(100, "And the page address should contain \"/func/Dbtests.PostgresOrders\"", () => urlShouldContain(page, "/func/Dbtests.PostgresOrders"));
    await session.step(101, "And \"Script\" accordion header in context panel should be visible", () => shouldBe(page, el("\"Script\" accordion header in context panel"), "visible"));
    await session.step(102, "And no errors should have been logged", () => noErrors(page));
    await session.step(103, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A table of a schema shows its details", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(106, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(107, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(108, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(109, "And the stand has a reachable \"PostgresTest\" connection", () => standHasReachableConnection(page, "PostgresTest"));
    await session.step(110, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(111, "And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree")));
    await session.step(112, "And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(113, "When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(114, "Then the context panel should show \"orders\"", () => contextPanelShows(page, "orders"));
    await session.step(115, "And no errors should have been logged", () => noErrors(page));
    await session.step(116, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Plugins opens the Plugins view [node=Plugins, view=Plugins]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(119, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(120, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(121, "When user clicks on Platform---Plugins tree node inside browse tree", () => clickOn(page, el("Platform---Plugins tree node inside browse tree")));
    await session.step(122, "Then the \"Plugins\" view should be current", () => viewIsCurrent(page, "Plugins"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Credentials opens the Credentials view [node=Credentials, view=Credentials]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(119, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(120, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(121, "When user clicks on Platform---Credentials tree node inside browse tree", () => clickOn(page, el("Platform---Credentials tree node inside browse tree")));
    await session.step(122, "Then the \"Credentials\" view should be current", () => viewIsCurrent(page, "Credentials"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Functions opens the Functions view [node=Functions, view=Functions]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(119, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(120, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(121, "When user clicks on Platform---Functions tree node inside browse tree", () => clickOn(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(122, "Then the \"Functions\" view should be current", () => viewIsCurrent(page, "Functions"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Roles opens the Roles view [node=Roles, view=Roles]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(119, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(120, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(121, "When user clicks on Platform---Roles tree node inside browse tree", () => clickOn(page, el("Platform---Roles tree node inside browse tree")));
    await session.step(122, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > MCP-Servers opens the MCP Servers view [node=MCP-Servers, view=MCP Servers]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(119, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(120, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(121, "When user clicks on Platform---MCP-Servers tree node inside browse tree", () => clickOn(page, el("Platform---MCP-Servers tree node inside browse tree")));
    await session.step(122, "Then the \"MCP Servers\" view should be current", () => viewIsCurrent(page, "MCP Servers"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Predictive-models opens the Models view [node=Predictive-models, view=Models]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(119, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(120, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(121, "When user clicks on Platform---Predictive-models tree node inside browse tree", () => clickOn(page, el("Platform---Predictive-models tree node inside browse tree")));
    await session.step(122, "Then the \"Models\" view should be current", () => viewIsCurrent(page, "Models"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Sync opens the Sync view [node=Sync, view=Sync]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(119, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(120, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(121, "When user clicks on Platform---Sync tree node inside browse tree", () => clickOn(page, el("Platform---Sync tree node inside browse tree")));
    await session.step(122, "Then the \"Sync\" view should be current", () => viewIsCurrent(page, "Sync"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Layouts opens the View layouts view [node=Layouts, view=View layouts]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(119, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(120, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(121, "When user clicks on Platform---Layouts tree node inside browse tree", () => clickOn(page, el("Platform---Layouts tree node inside browse tree")));
    await session.step(122, "Then the \"View layouts\" view should be current", () => viewIsCurrent(page, "View layouts"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > URL-Aliases opens the URL Aliases view [node=URL-Aliases, view=URL Aliases]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(119, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(120, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(121, "When user clicks on Platform---URL-Aliases tree node inside browse tree", () => clickOn(page, el("Platform---URL-Aliases tree node inside browse tree")));
    await session.step(122, "Then the \"URL Aliases\" view should be current", () => viewIsCurrent(page, "URL Aliases"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Settings opens the Settings view [node=Settings, view=Settings]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(119, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(120, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(121, "When user clicks on Platform---Settings tree node inside browse tree", () => clickOn(page, el("Platform---Settings tree node inside browse tree")));
    await session.step(122, "Then the \"Settings\" view should be current", () => viewIsCurrent(page, "Settings"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Sticky-Meta opens the Schemas view [node=Sticky-Meta, view=Schemas]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(119, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(120, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(121, "When user clicks on Platform---Sticky-Meta tree node inside browse tree", () => clickOn(page, el("Platform---Sticky-Meta tree node inside browse tree")));
    await session.step(122, "Then the \"Schemas\" view should be current", () => viewIsCurrent(page, "Schemas"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Notebooks opens the Notebooks view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(142, "Given the \"Notebooks\" package is installed", () => packageInstalled(page, "Notebooks"));
    await session.step(143, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(144, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(145, "When user clicks on Platform---Notebooks tree node inside browse tree", () => clickOn(page, el("Platform---Notebooks tree node inside browse tree")));
    await session.step(146, "Then the \"Notebooks\" view should be current", () => viewIsCurrent(page, "Notebooks"));
    await session.step(147, "And no errors should have been logged", () => noErrors(page));
    await session.step(148, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Plugins view lists the installed packages and offers no Delete for one", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(151, "Given the \"PowerPack\" package is installed", () => packageInstalled(page, "PowerPack"));
    await session.step(152, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(153, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(154, "When user clicks on Platform---Plugins tree node inside browse tree", () => clickOn(page, el("Platform---Plugins tree node inside browse tree")));
    await session.step(155, "Then the \"Plugins\" view should be current", () => viewIsCurrent(page, "Plugins"));
    await session.step(156, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"PowerPack\" gallery card"],["\"Tutorials\" gallery card"]]), [["\"PowerPack\" gallery card"],["\"Tutorials\" gallery card"]]);
    await session.step(159, "When user opens the context menu of \"PowerPack\" gallery card", () => openContextMenu(page, el("\"PowerPack\" gallery card")));
    await session.step(160, "Then the open menu should list \"Uninstall\"", () => menuLists(page, "Uninstall"));
    await session.step(161, "And the open menu should not list \"Delete\"", () => menuDoesNotList(page, "Delete"));
    await session.step(162, "And the open menu should not list \"Delete...\"", () => menuDoesNotList(page, "Delete..."));
    await session.step(163, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(164, "Then no errors should have been logged", () => noErrors(page));
    await session.step(165, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Dockers opens the Dockers view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(169, "Given the stand runs the \"Grok Spawner\" service", () => standRunsService(page, "Grok Spawner"));
    await session.step(170, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(171, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(172, "When user clicks on Platform---Dockers tree node inside browse tree", () => clickOn(page, el("Platform---Dockers tree node inside browse tree")));
    await session.step(173, "Then the \"Dockers\" view should be current", () => viewIsCurrent(page, "Dockers"));
    await session.step(174, "And no errors should have been logged", () => noErrors(page));
    await session.step(175, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
