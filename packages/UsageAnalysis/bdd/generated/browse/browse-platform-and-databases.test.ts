/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-platform-and-databases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
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
import {firstSignedIn, loggedIn, signBackInAsFirst, signInAsSecond, signedInAs} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, followingShouldBe, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, connectionOnServer, contextPanelOpen, contextPanelShows, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Platform and Databases sections of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-platform-and-databases.feature", import.meta.url);
  test("The Platform section lists what an administrator manages", {tag: ["@browse", "@realizes:views.browse", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(42, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(43, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Platform---Plugins tree node inside browse tree"],["Platform---Credentials tree node inside browse tree"],["Platform---Functions tree node inside browse tree"],["Platform---Users tree node inside browse tree"],["Platform---Groups tree node inside browse tree"],["Platform---Roles tree node inside browse tree"],["Platform---Predictive-models tree node inside browse tree"],["Platform---Dockers tree node inside browse tree"]]), [["Platform---Plugins tree node inside browse tree"],["Platform---Credentials tree node inside browse tree"],["Platform---Functions tree node inside browse tree"],["Platform---Users tree node inside browse tree"],["Platform---Groups tree node inside browse tree"],["Platform---Roles tree node inside browse tree"],["Platform---Predictive-models tree node inside browse tree"],["Platform---Dockers tree node inside browse tree"]]);
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
    await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A Platform node opens the view it is named after", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(56, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(57, "When user clicks on Platform---Users tree node inside browse tree", () => clickOn(page, el("Platform---Users tree node inside browse tree")));
    await session.step(58, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
    await session.step(59, "When user clicks on Platform---Groups tree node inside browse tree", () => clickOn(page, el("Platform---Groups tree node inside browse tree")));
    await session.step(60, "Then the \"Groups\" view should be current", () => viewIsCurrent(page, "Groups"));
    await session.step(61, "And no errors should have been logged", () => noErrors(page));
    await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Databases section lists the connected providers", {tag: ["@browse", "@realizes:views.browse", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(66, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(67, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Databases---Postgres tree node inside browse tree"],["Databases---MySQL tree node inside browse tree"],["Databases---Oracle tree node inside browse tree"],["Databases---MariaDB tree node inside browse tree"],["Databases---MS-SQL tree node inside browse tree"]]), [["Databases---Postgres tree node inside browse tree"],["Databases---MySQL tree node inside browse tree"],["Databases---Oracle tree node inside browse tree"],["Databases---MariaDB tree node inside browse tree"],["Databases---MS-SQL tree node inside browse tree"]]);
    await session.step(73, "And no errors should have been logged", () => noErrors(page));
    await session.step(74, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A provider opens down to its connections", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(77, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(78, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(79, "Then Databases---Postgres---Datagrok tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---Datagrok tree node inside browse tree"), "visible"));
    await session.step(80, "And Databases---Postgres---CHEMBL tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---CHEMBL tree node inside browse tree"), "visible"));
    await session.step(81, "When user collapses Databases---Postgres tree node inside browse tree", () => collapse(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(82, "Then Databases---Postgres---Datagrok tree node inside browse tree should be hidden", () => shouldBe(page, el("Databases---Postgres---Datagrok tree node inside browse tree"), "hidden"));
    await session.step(83, "And Databases---Postgres tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres tree node inside browse tree"), "visible"));
    await session.step(84, "And no errors should have been logged", () => noErrors(page));
    await session.step(85, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A saved query opens from the tree with its details", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(88, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(89, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(90, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(91, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(92, "When user clicks on Databases---Postgres---NorthwindTest---Orders tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Orders tree node inside browse tree")));
    await session.step(93, "Then the context panel should show \"Orders\"", () => contextPanelShows(page, "Orders"));
    await session.step(94, "And the page address should contain \"/func/Dbtests.PostgresOrders\"", () => urlShouldContain(page, "/func/Dbtests.PostgresOrders"));
    await session.step(95, "And \"Script\" accordion header in context panel should be visible", () => shouldBe(page, el("\"Script\" accordion header in context panel"), "visible"));
    await session.step(96, "And no errors should have been logged", () => noErrors(page));
    await session.step(97, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A table of a schema shows its details", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(100, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(101, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(102, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(103, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(104, "And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree")));
    await session.step(105, "And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(106, "When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(107, "Then the context panel should show \"orders\"", () => contextPanelShows(page, "orders"));
    await session.step(108, "And no errors should have been logged", () => noErrors(page));
    await session.step(109, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Plugins opens the Plugins view [node=Plugins, view=Plugins]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---Plugins tree node inside browse tree", () => clickOn(page, el("Platform---Plugins tree node inside browse tree")));
    await session.step(115, "Then the \"Plugins\" view should be current", () => viewIsCurrent(page, "Plugins"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Credentials opens the Credentials view [node=Credentials, view=Credentials]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---Credentials tree node inside browse tree", () => clickOn(page, el("Platform---Credentials tree node inside browse tree")));
    await session.step(115, "Then the \"Credentials\" view should be current", () => viewIsCurrent(page, "Credentials"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Functions opens the Functions view [node=Functions, view=Functions]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---Functions tree node inside browse tree", () => clickOn(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(115, "Then the \"Functions\" view should be current", () => viewIsCurrent(page, "Functions"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Roles opens the Roles view [node=Roles, view=Roles]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---Roles tree node inside browse tree", () => clickOn(page, el("Platform---Roles tree node inside browse tree")));
    await session.step(115, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Notebooks opens the Notebooks view [node=Notebooks, view=Notebooks]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---Notebooks tree node inside browse tree", () => clickOn(page, el("Platform---Notebooks tree node inside browse tree")));
    await session.step(115, "Then the \"Notebooks\" view should be current", () => viewIsCurrent(page, "Notebooks"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > MCP-Servers opens the MCP Servers view [node=MCP-Servers, view=MCP Servers]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---MCP-Servers tree node inside browse tree", () => clickOn(page, el("Platform---MCP-Servers tree node inside browse tree")));
    await session.step(115, "Then the \"MCP Servers\" view should be current", () => viewIsCurrent(page, "MCP Servers"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Predictive-models opens the Models view [node=Predictive-models, view=Models]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---Predictive-models tree node inside browse tree", () => clickOn(page, el("Platform---Predictive-models tree node inside browse tree")));
    await session.step(115, "Then the \"Models\" view should be current", () => viewIsCurrent(page, "Models"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Dockers opens the Dockers view [node=Dockers, view=Dockers]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---Dockers tree node inside browse tree", () => clickOn(page, el("Platform---Dockers tree node inside browse tree")));
    await session.step(115, "Then the \"Dockers\" view should be current", () => viewIsCurrent(page, "Dockers"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Sync opens the Sync view [node=Sync, view=Sync]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---Sync tree node inside browse tree", () => clickOn(page, el("Platform---Sync tree node inside browse tree")));
    await session.step(115, "Then the \"Sync\" view should be current", () => viewIsCurrent(page, "Sync"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Layouts opens the View layouts view [node=Layouts, view=View layouts]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---Layouts tree node inside browse tree", () => clickOn(page, el("Platform---Layouts tree node inside browse tree")));
    await session.step(115, "Then the \"View layouts\" view should be current", () => viewIsCurrent(page, "View layouts"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > URL-Aliases opens the URL Aliases view [node=URL-Aliases, view=URL Aliases]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---URL-Aliases tree node inside browse tree", () => clickOn(page, el("Platform---URL-Aliases tree node inside browse tree")));
    await session.step(115, "Then the \"URL Aliases\" view should be current", () => viewIsCurrent(page, "URL Aliases"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Settings opens the Settings view [node=Settings, view=Settings]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---Settings tree node inside browse tree", () => clickOn(page, el("Platform---Settings tree node inside browse tree")));
    await session.step(115, "Then the \"Settings\" view should be current", () => viewIsCurrent(page, "Settings"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Platform > Sticky-Meta opens the Schemas view [node=Sticky-Meta, view=Schemas]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(114, "When user clicks on Platform---Sticky-Meta tree node inside browse tree", () => clickOn(page, el("Platform---Sticky-Meta tree node inside browse tree")));
    await session.step(115, "Then the \"Schemas\" view should be current", () => viewIsCurrent(page, "Schemas"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Another user does not see a connection nobody shared with them", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(138, "Given a \"Postgres\" connection named \"BDD-Browse-Private-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Browse-Private-{run}")));
    await session.step(139, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(140, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(141, "Then Databases---Postgres---BDD-Browse-Private-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Browse-Private-{run} tree node inside browse tree")), "visible"));
    await session.step(142, "When user signs in as the second user", () => signInAsSecond(page));
    await session.step(143, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(144, "Then the second user should be signed in", () => signedInAs(page));
    await session.step(145, "When Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(146, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(147, "Then Databases---Postgres---Datagrok tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---Datagrok tree node inside browse tree"), "visible"));
    await session.step(148, "And Databases---Postgres---BDD-Browse-Private-{run} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Browse-Private-{run} tree node inside browse tree")), "absent"));
    await session.step(149, "And no errors should have been logged", () => noErrors(page));
    await session.step(150, "When user signs in again as the first user", () => signBackInAsFirst(page));
    await session.step(151, "Then the first user should be signed in", () => firstSignedIn(page));
  });
});
