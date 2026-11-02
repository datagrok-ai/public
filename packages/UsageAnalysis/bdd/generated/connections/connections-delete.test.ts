/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-delete.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
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
import {clickOn, enterInto, isExpanded, replaceCode, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, connectionOnServer, connectionsOnServer, currentViewType, dialogCloses, noQueryOnServer, queriesOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Deleting a database connection", () => {
  const session = feature(test, "features/connections/connections-delete.feature", import.meta.url);
  test("DELETE removes the connection from the server and the tree", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(21, "Given a \"Postgres\" connection named \"BDD-Conn-Delete-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Delete-{run}")));
    await session.step(22, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(23, "When user picks \"Delete...\" from the context menu of Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el(session.text("Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree"))));
    await session.step(24, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(25, "And \"Are you sure?\" dialog should contain text \"Delete connection \\\"BDD-Conn-Delete-{run}\\\"?\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete connection \"BDD-Conn-Delete-{run}\"?")));
    await session.step(26, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
    await session.step(27, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(28, "And 0 connections named \"BDD-Conn-Delete-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Delete-{run}")));
    await session.step(29, "And Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Delete-{run} tree node inside browse tree")), "absent"));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
  test("CANCEL keeps the connection", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(34, "Given a \"Postgres\" connection named \"BDD-Conn-Keep-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Keep-{run}")));
    await session.step(35, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(36, "When user picks \"Delete...\" from the context menu of Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el(session.text("Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree"))));
    await session.step(37, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(38, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
    await session.step(39, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(40, "And 1 connection named \"BDD-Conn-Keep-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Keep-{run}")));
    await session.step(41, "And Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Keep-{run} tree node inside browse tree")), "visible"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
  test("A connection is deleted from the connections gallery", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(45, "Given a \"Postgres\" connection named \"BDD-Conn-Gallery-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Gallery-{run}")));
    await session.step(46, "When user picks \"Browse connections\" from the context menu of Databases---Postgres tree node inside browse tree", () => pickFromContextMenu(page, "Browse connections", el("Databases---Postgres tree node inside browse tree")));
    await session.step(47, "Then the \"Postgres\" view should be current", () => viewIsCurrent(page, "Postgres"));
    await session.step(48, "When user types \"BDD-Conn-Gallery-{run}\" into gallery search", () => typeInto(page, session.text("BDD-Conn-Gallery-{run}"), el("gallery search")));
    await session.step(49, "Then \"BDD-Conn-Gallery-{run}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-Conn-Gallery-{run}\" link in gallery")), "visible"));
    await session.step(50, "When user picks \"Delete...\" from the context menu of \"BDD-Conn-Gallery-{run}\" link in gallery", () => pickFromContextMenu(page, "Delete...", el(session.text("\"BDD-Conn-Gallery-{run}\" link in gallery"))));
    await session.step(51, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(52, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
    await session.step(53, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(54, "And 0 connections named \"BDD-Conn-Gallery-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Gallery-{run}")));
    await session.step(55, "And \"BDD-Conn-Gallery-{run}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BDD-Conn-Gallery-{run}\" link in gallery")), "absent"));
    await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A query saved under a connection is deleted with it", {tag: ["@connections", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(60, "Given no query named \"BDDConnCascadeQ{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDConnCascadeQ{time}")));
    await session.step(61, "And a \"Postgres\" connection named \"BDDConnCascade{time}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDDConnCascade{time}")));
    await session.step(62, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(63, "When user picks \"New Query...\" from the context menu of Databases---Postgres---BDDConnCascade{time} tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el(session.text("Databases---Postgres---BDDConnCascade{time} tree node inside browse tree"))));
    await session.step(64, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
    await session.step(65, "When user enters \"BDDConnCascadeQ{time}\" into Name input", () => enterInto(page, session.text("BDDConnCascadeQ{time}"), el("Name input")));
    await session.step(66, "And user replaces the code of code editor with \"select * from products\"", () => replaceCode(page, el("code editor"), "select * from products"));
    await session.step(67, "And user clicks on Save button", () => clickOn(page, el("Save button")));
    await session.step(68, "Then 1 query named \"BDDConnCascadeQ{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDConnCascadeQ{time}")));
    await session.step(69, "When user closes all views", () => closeAllViews(page));
    await session.step(70, "And user picks \"Delete...\" from the context menu of Databases---Postgres---BDDConnCascade{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el(session.text("Databases---Postgres---BDDConnCascade{time} tree node inside browse tree"))));
    await session.step(71, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(72, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
    await session.step(73, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(74, "And 0 connections named \"BDDConnCascade{time}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDDConnCascade{time}")));
    await session.step(75, "And 0 queries named \"BDDConnCascadeQ{time}\" should be on the server", () => queriesOnServer(page, 0, session.text("BDDConnCascadeQ{time}")));
    await session.step(76, "And Databases---Postgres---BDDConnCascade{time} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Databases---Postgres---BDDConnCascade{time} tree node inside browse tree")), "absent"));
    await session.step(77, "And no errors should have been logged", () => noErrors(page));
  });
});
