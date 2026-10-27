/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connection-edge-cases.feature
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
import {clickOn, enterInto, isExpanded, replaceCode, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, connectionOnServer, connectionsOnServer, currentViewType, dialogCloses, noConnectionOnServer, noQueryOnServer, queriesOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Connection edge cases: names the dialogs refuse, CANCEL, and a delete that takes the queries", () => {
  const session = feature(test, "features/connections/connection-edge-cases.feature", import.meta.url);
  test("A name of spaces keeps OK disabled, and CANCEL saves nothing", {tag: ["@connections", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(23, "Given no connection named \"BDDConnSpaces{time}\" is on the server", () => noConnectionOnServer(page, session.text("BDDConnSpaces{time}")));
    await session.step(24, "When user picks \"New connection...\" from the context menu of Databases---Postgres tree node inside browse tree", () => pickFromContextMenu(page, "New connection...", el("Databases---Postgres tree node inside browse tree")));
    await session.step(25, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(26, "When user enters \"   \" into Name input in \"Add new connection\" dialog", () => enterInto(page, "   ", el("Name input in \"Add new connection\" dialog")));
    await session.step(27, "And user enters \"db.datagrok.ai\" into Server input in \"Add new connection\" dialog", () => enterInto(page, "db.datagrok.ai", el("Server input in \"Add new connection\" dialog")));
    await session.step(28, "Then OK button in \"Add new connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "disabled"));
    await session.step(29, "When user enters \"BDDConnSpaces{time}\" into Name input in \"Add new connection\" dialog", () => enterInto(page, session.text("BDDConnSpaces{time}"), el("Name input in \"Add new connection\" dialog")));
    await session.step(30, "Then OK button in \"Add new connection\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "enabled"));
    await session.step(31, "When user clicks on CANCEL button in \"Add new connection\" dialog", () => clickOn(page, el("CANCEL button in \"Add new connection\" dialog")));
    await session.step(32, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(33, "And 0 connections named \"BDDConnSpaces{time}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDDConnSpaces{time}")));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
  });
  test("CANCEL in the Edit dialog keeps the name, and a name of spaces is refused there too", {tag: ["@connections", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(37, "Given a \"Postgres\" connection named \"BDDConnEditEdge{time}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDDConnEditEdge{time}")));
    await session.step(38, "And no connection named \"BDDConnEditedEdge{time}\" is on the server", () => noConnectionOnServer(page, session.text("BDDConnEditedEdge{time}")));
    await session.step(39, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(40, "When user picks \"Edit...\" from the context menu of Databases---Postgres---BDDConnEditEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---BDDConnEditEdge{time} tree node inside browse tree"))));
    await session.step(41, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
    await session.step(42, "When user enters \"   \" into Name input in \"Edit Connection\" dialog", () => enterInto(page, "   ", el("Name input in \"Edit Connection\" dialog")));
    await session.step(43, "Then OK button in \"Edit Connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Edit Connection\" dialog"), "disabled"));
    await session.step(44, "When user enters \"BDDConnEditedEdge{time}\" into Name input in \"Edit Connection\" dialog", () => enterInto(page, session.text("BDDConnEditedEdge{time}"), el("Name input in \"Edit Connection\" dialog")));
    await session.step(45, "Then OK button in \"Edit Connection\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Edit Connection\" dialog"), "enabled"));
    await session.step(46, "When user clicks on CANCEL button in \"Edit Connection\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Connection\" dialog")));
    await session.step(47, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
    await session.step(48, "And 1 connection named \"BDDConnEditEdge{time}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDDConnEditEdge{time}")));
    await session.step(49, "And 0 connections named \"BDDConnEditedEdge{time}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDDConnEditedEdge{time}")));
    await session.step(50, "And Databases---Postgres---BDDConnEditEdge{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDDConnEditEdge{time} tree node inside browse tree")), "visible"));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("Rename... refuses a name with no letter", {tag: ["@connections", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(54, "Given a \"Postgres\" connection named \"BDDConnRenameEdge{time}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDDConnRenameEdge{time}")));
    await session.step(55, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(56, "When user picks \"Rename...\" from the context menu of Databases---Postgres---BDDConnRenameEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Databases---Postgres---BDDConnRenameEdge{time} tree node inside browse tree"))));
    await session.step(57, "Then \"Rename dataconnection\" dialog should be visible", () => shouldBe(page, el("\"Rename dataconnection\" dialog"), "visible"));
    await session.step(58, "When user enters \"2026\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, "2026", el("Name input in \"Rename dataconnection\" dialog")));
    await session.step(59, "Then OK button in \"Rename dataconnection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Rename dataconnection\" dialog"), "disabled"));
    await session.step(60, "When user enters \"BDDConnRenamedEdge{time}\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, session.text("BDDConnRenamedEdge{time}"), el("Name input in \"Rename dataconnection\" dialog")));
    await session.step(61, "Then OK button in \"Rename dataconnection\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Rename dataconnection\" dialog"), "enabled"));
    await session.step(62, "When user clicks on CANCEL button in \"Rename dataconnection\" dialog", () => clickOn(page, el("CANCEL button in \"Rename dataconnection\" dialog")));
    await session.step(63, "Then the \"Rename dataconnection\" dialog should close", () => dialogCloses(page, "Rename dataconnection"));
    await session.step(64, "And 1 connection named \"BDDConnRenameEdge{time}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDDConnRenameEdge{time}")));
    await session.step(65, "And no errors should have been logged", () => noErrors(page));
  });
  test("Deleting a connection deletes the query saved under it", {tag: ["@connections", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(68, "Given no query named \"BDDConnCascadeQ{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDConnCascadeQ{time}")));
    await session.step(69, "And a \"Postgres\" connection named \"BDDConnCascade{time}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDDConnCascade{time}")));
    await session.step(70, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(71, "When user picks \"New Query...\" from the context menu of Databases---Postgres---BDDConnCascade{time} tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el(session.text("Databases---Postgres---BDDConnCascade{time} tree node inside browse tree"))));
    await session.step(72, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
    await session.step(73, "When user enters \"BDDConnCascadeQ{time}\" into Name input", () => enterInto(page, session.text("BDDConnCascadeQ{time}"), el("Name input")));
    await session.step(74, "And user replaces the code of code editor with \"select * from products\"", () => replaceCode(page, el("code editor"), "select * from products"));
    await session.step(75, "And user clicks on Save button", () => clickOn(page, el("Save button")));
    await session.step(76, "Then 1 query named \"BDDConnCascadeQ{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDConnCascadeQ{time}")));
    await session.step(77, "When user closes all views", () => closeAllViews(page));
    await session.step(78, "And user picks \"Delete...\" from the context menu of Databases---Postgres---BDDConnCascade{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el(session.text("Databases---Postgres---BDDConnCascade{time} tree node inside browse tree"))));
    await session.step(79, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(80, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
    await session.step(81, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(82, "And 0 connections named \"BDDConnCascade{time}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDDConnCascade{time}")));
    await session.step(83, "And 0 queries named \"BDDConnCascadeQ{time}\" should be on the server", () => queriesOnServer(page, 0, session.text("BDDConnCascadeQ{time}")));
    await session.step(84, "And Databases---Postgres---BDDConnCascade{time} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Databases---Postgres---BDDConnCascade{time} tree node inside browse tree")), "absent"));
    await session.step(85, "And no errors should have been logged", () => noErrors(page));
  });
});
