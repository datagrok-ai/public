/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-sparql.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/grid.js';
import '../../bindings/nx.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {hiddenProvidersShown} from '../../bindings/connections.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, followingShouldBe, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, connectionsOnServer, dialogCloses, noConnectionOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noBalloons, noErrors, openContextMenu, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A SPARQL connection", () => {
  const session = feature(test, "features/connections/connections-sparql.feature", import.meta.url);
  test("Show more reveals Sparql among the hidden providers", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(20, "Then Databases---Sparql tree node inside browse tree should be hidden", () => shouldBe(page, el("Databases---Sparql tree node inside browse tree"), "hidden"));
    await session.step(21, "When user clicks on \"ellipsis-h\" icon inside Databases---Show-more tree node inside browse tree", () => clickOn(page, el("\"ellipsis-h\" icon inside Databases---Show-more tree node inside browse tree")));
    await session.step(22, "Then Databases---Sparql tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Sparql tree node inside browse tree"), "visible"));
    await session.step(23, "And Databases---Show-more tree node inside browse tree should be hidden", () => shouldBe(page, el("Databases---Show-more tree node inside browse tree"), "hidden"));
    await session.step(26, "When user clicks on Databases---Sparql tree node inside browse tree", () => clickOn(page, el("Databases---Sparql tree node inside browse tree")));
    await session.step(27, "And user opens the context menu of Databases---Sparql tree node inside browse tree", () => openContextMenu(page, el("Databases---Sparql tree node inside browse tree")));
    await session.step(28, "Then the open menu should list \"New connection...\"", () => menuLists(page, "New connection..."));
    await session.step(29, "And the open menu should list \"Browse connections\"", () => menuLists(page, "Browse connections"));
    await session.step(30, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(31, "Then no errors should have been logged", () => noErrors(page));
  });
  test("A Sparql connection is saved from its dialog and deleted", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(34, "Given no connection named \"BDD-Conn-Sparql-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-Sparql-{run}")));
    await session.step(35, "And the hidden providers of the Databases tree are shown", () => hiddenProvidersShown(page));
    await session.step(36, "When user clicks on Databases---Sparql tree node inside browse tree", () => clickOn(page, el("Databases---Sparql tree node inside browse tree")));
    await session.step(37, "And user picks \"New connection...\" from the context menu of Databases---Sparql tree node inside browse tree", () => pickFromContextMenu(page, "New connection...", el("Databases---Sparql tree node inside browse tree")));
    await session.step(38, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(39, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Endpoint input in \"Add new connection\" dialog"],["Prefixes input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]), [["Endpoint input in \"Add new connection\" dialog"],["Prefixes input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]);
    await session.step(43, "When user enters \"BDD-Conn-Sparql-{run}\" into Name input in \"Add new connection\" dialog", () => enterInto(page, session.text("BDD-Conn-Sparql-{run}"), el("Name input in \"Add new connection\" dialog")));
    await session.step(44, "And user enters \"http://data.ontotext.com/repositories/data-last\" into Endpoint input in \"Add new connection\" dialog", () => enterInto(page, "http://data.ontotext.com/repositories/data-last", el("Endpoint input in \"Add new connection\" dialog")));
    await session.step(45, "And user clicks on OK button in \"Add new connection\" dialog", () => clickOn(page, el("OK button in \"Add new connection\" dialog")));
    await session.step(46, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(47, "And 1 connection named \"BDD-Conn-Sparql-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Sparql-{run}")));
    await session.step(48, "Given Databases---Sparql tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Sparql tree node inside browse tree")));
    await session.step(49, "When user picks \"Delete...\" from the context menu of Databases---Sparql---BDD-Conn-Sparql-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el(session.text("Databases---Sparql---BDD-Conn-Sparql-{run} tree node inside browse tree"))));
    await session.step(50, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
    await session.step(51, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(52, "And 0 connections named \"BDD-Conn-Sparql-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Sparql-{run}")));
    await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
