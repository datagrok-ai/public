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
import {clickOn, isExpanded, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, connectionOnServer, connectionsOnServer, dialogCloses, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Deleting a database connection", () => {
  const session = feature(test, "features/connections/connections-delete.feature", import.meta.url);
  test("DELETE removes the connection from the server and the tree", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(15, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(18, "Given a \"Postgres\" connection named \"BDD-Conn-Delete-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Delete-{run}")));
    await session.step(19, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(20, "When user picks \"Delete...\" from the context menu of Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el(session.text("Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree"))));
    await session.step(21, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(22, "And \"Are you sure?\" dialog should contain text \"Delete connection \\\"BDD-Conn-Delete-{run}\\\"?\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete connection \"BDD-Conn-Delete-{run}\"?")));
    await session.step(23, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
    await session.step(24, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(25, "And 0 connections named \"BDD-Conn-Delete-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Delete-{run}")));
    await session.step(26, "And Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree")), "absent"));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
  });
  test("CANCEL keeps the connection", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(15, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(31, "Given a \"Postgres\" connection named \"BDD-Conn-Keep-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Keep-{run}")));
    await session.step(32, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(33, "When user picks \"Delete...\" from the context menu of Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el(session.text("Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree"))));
    await session.step(34, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(35, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
    await session.step(36, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(37, "And 1 connection named \"BDD-Conn-Keep-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Keep-{run}")));
    await session.step(38, "And Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree")), "visible"));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("A connection is deleted from the connections gallery", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(15, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(42, "Given a \"Postgres\" connection named \"BDD-Conn-Gallery-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Gallery-{run}")));
    await session.step(43, "When user picks \"Browse connections\" from the context menu of Databases---Postgres tree node inside browse tree", () => pickFromContextMenu(page, "Browse connections", el("Databases---Postgres tree node inside browse tree")));
    await session.step(44, "Then the \"Postgres\" view should be current", () => viewIsCurrent(page, "Postgres"));
    await session.step(45, "When user types \"BDD-Conn-Gallery-{run}\" into gallery search", () => typeInto(page, session.text("BDD-Conn-Gallery-{run}"), el("gallery search")));
    await session.step(46, "Then \"BDD-Conn-Gallery-{run}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-Conn-Gallery-{run}\" link in gallery")), "visible"));
    await session.step(47, "When user picks \"Delete...\" from the context menu of \"BDD-Conn-Gallery-{run}\" link in gallery", () => pickFromContextMenu(page, "Delete...", el(session.text("\"BDD-Conn-Gallery-{run}\" link in gallery"))));
    await session.step(48, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(49, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
    await session.step(50, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(51, "And 0 connections named \"BDD-Conn-Gallery-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Gallery-{run}")));
    await session.step(52, "And \"BDD-Conn-Gallery-{run}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BDD-Conn-Gallery-{run}\" link in gallery")), "absent"));
    await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
