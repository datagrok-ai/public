/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-delete.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
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
import {clickOn, isExpanded, rightClickOn, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, connectionOnServer, connectionsOnServer, dialogCloses, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Deleting a database connection", () => {
  const session = feature(test, "features/connections/connections-delete.feature", import.meta.url);
  test("DELETE removes the connection from the server and the tree", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(22, "Given a \"Postgres\" connection named \"BDD-Conn-Delete-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Delete-{run}")));
    await session.step(23, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(24, "When user right-clicks on Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree"))));
    await session.step(25, "And user picks \"Delete...\" from the open menu", () => pickFromOpenMenu(page, "Delete..."));
    await session.step(26, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(27, "And \"Are you sure?\" dialog should contain text \"Delete connection \\\"BDD-Conn-Delete-{run}\\\"?\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete connection \"BDD-Conn-Delete-{run}\"?")));
    await session.step(28, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
    await session.step(29, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(30, "And 0 connections named \"BDD-Conn-Delete-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Delete-{run}")));
    await session.step(31, "And Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree")), "absent"));
    await session.step(32, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
  });
  test("CANCEL keeps the connection", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(36, "Given a \"Postgres\" connection named \"BDD-Conn-Keep-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Keep-{run}")));
    await session.step(37, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(38, "When user right-clicks on Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree"))));
    await session.step(39, "And user picks \"Delete...\" from the open menu", () => pickFromOpenMenu(page, "Delete..."));
    await session.step(40, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(41, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
    await session.step(42, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(43, "And 1 connection named \"BDD-Conn-Keep-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Keep-{run}")));
    await session.step(44, "And Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree")), "visible"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
  test("A connection is deleted from the connections gallery", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(48, "Given a \"Postgres\" connection named \"BDD-Conn-Gallery-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Gallery-{run}")));
    await session.step(49, "When user right-clicks on Databases---Postgres tree node inside browse tree", () => rightClickOn(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(50, "And user picks \"Browse connections\" from the open menu", () => pickFromOpenMenu(page, "Browse connections"));
    await session.step(51, "Then the \"Postgres\" view should be current", () => viewIsCurrent(page, "Postgres"));
    await session.step(52, "When user types \"BDD-Conn-Gallery-{run}\" into gallery search", () => typeInto(page, session.text("BDD-Conn-Gallery-{run}"), el("gallery search")));
    await session.step(53, "Then \"BDD-Conn-Gallery-{run}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-Conn-Gallery-{run}\" link in gallery")), "visible"));
    await session.step(54, "When user picks \"Delete...\" from the context menu of \"BDD-Conn-Gallery-{run}\" link in gallery", () => pickFromContextMenu(page, "Delete...", el(session.text("\"BDD-Conn-Gallery-{run}\" link in gallery"))));
    await session.step(55, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(56, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
    await session.step(57, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(58, "And 0 connections named \"BDD-Conn-Gallery-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Gallery-{run}")));
    await session.step(59, "And \"BDD-Conn-Gallery-{run}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BDD-Conn-Gallery-{run}\" link in gallery")), "absent"));
    await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
