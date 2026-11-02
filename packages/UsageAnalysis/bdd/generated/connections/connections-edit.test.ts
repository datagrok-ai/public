/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-edit.feature
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
import {clickOn, enterInto, isExpanded, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, connectionOnServer, connectionsOnServer, dialogCloses, noConnectionOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Editing a database connection", () => {
  const session = feature(test, "features/connections/connections-edit.feature", import.meta.url);
  test("Editing a database connection", {tag: ["@connections", "@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(24, "And no connection named \"BDD-Conn-Edit-{run}, BDD-Conn-Edited-{run}, BDD-Conn-Renamed-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-Edit-{run}, BDD-Conn-Edited-{run}, BDD-Conn-Renamed-{run}")));
    await session.step(25, "And a \"Postgres\" connection named \"BDD-Conn-Edit-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Edit-{run}")));
    await session.step(26, "And a \"Postgres\" connection named \"BDDConnEditEdge{time}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDDConnEditEdge{time}")));
    await session.step(27, "And no connection named \"BDDConnEditedEdge{time}, BDDConnRenamedEdge{time}\" is on the server", () => noConnectionOnServer(page, session.text("BDDConnEditedEdge{time}, BDDConnRenamedEdge{time}")));
    await session.step(28, "And a \"Postgres\" connection named \"BDDConnRenameEdge{time}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDDConnRenameEdge{time}")));
    await session.step(29, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await run.scenario("The Edit dialog shows what the connection was saved with, refuses a name of spaces, and keeps the name on CANCEL", async () => {
      await session.step(32, "When user picks \"Edit...\" from the context menu of Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree"))));
      await session.step(33, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(34, "And Name input in \"Edit Connection\" dialog should have value \"BDD-Conn-Edit-{run}\"", () => shouldHaveValue(page, el("Name input in \"Edit Connection\" dialog"), session.text("BDD-Conn-Edit-{run}")));
      await session.step(35, "And Server input in \"Edit Connection\" dialog should have value \"db.datagrok.ai\"", () => shouldHaveValue(page, el("Server input in \"Edit Connection\" dialog"), "db.datagrok.ai"));
      await session.step(36, "And Db input in \"Edit Connection\" dialog should have value \"northwind\"", () => shouldHaveValue(page, el("Db input in \"Edit Connection\" dialog"), "northwind"));
      await session.step(37, "When user clicks on CANCEL button in \"Edit Connection\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Connection\" dialog")));
      await session.step(38, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
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
    await run.scenario("Renaming through the Edit dialog", async () => {
      await session.step(54, "When user picks \"Edit...\" from the context menu of Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree"))));
      await session.step(55, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(56, "When user enters \"BDD-Conn-Edited-{run}\" into Name input in \"Edit Connection\" dialog", () => enterInto(page, session.text("BDD-Conn-Edited-{run}"), el("Name input in \"Edit Connection\" dialog")));
      await session.step(57, "And user clicks on OK button in \"Edit Connection\" dialog", () => clickOn(page, el("OK button in \"Edit Connection\" dialog")));
      await session.step(58, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(59, "And 1 connection named \"BDD-Conn-Edited-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Edited-{run}")));
      await session.step(60, "And 0 connections named \"BDD-Conn-Edit-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Edit-{run}")));
      await session.step(61, "And Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree")), "visible"));
      await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Renaming through Rename..., which refuses a name with no letter and keeps the name on CANCEL", async () => {
      await session.step(65, "When user picks \"Rename...\" from the context menu of Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree"))));
      await session.step(66, "Then \"Rename dataconnection\" dialog should be visible", () => shouldBe(page, el("\"Rename dataconnection\" dialog"), "visible"));
      await session.step(67, "When user enters \"BDD-Conn-Renamed-{run}\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, session.text("BDD-Conn-Renamed-{run}"), el("Name input in \"Rename dataconnection\" dialog")));
      await session.step(68, "And user clicks on OK button in \"Rename dataconnection\" dialog", () => clickOn(page, el("OK button in \"Rename dataconnection\" dialog")));
      await session.step(69, "Then the \"Rename dataconnection\" dialog should close", () => dialogCloses(page, "Rename dataconnection"));
      await session.step(70, "And 1 connection named \"BDD-Conn-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Renamed-{run}")));
      await session.step(71, "And 0 connections named \"BDD-Conn-Edited-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Edited-{run}")));
      await session.step(72, "And Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree")), "visible"));
      await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(74, "When user picks \"Rename...\" from the context menu of Databases---Postgres---BDDConnRenameEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Databases---Postgres---BDDConnRenameEdge{time} tree node inside browse tree"))));
      await session.step(75, "Then \"Rename dataconnection\" dialog should be visible", () => shouldBe(page, el("\"Rename dataconnection\" dialog"), "visible"));
      await session.step(76, "When user enters \"2026\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, "2026", el("Name input in \"Rename dataconnection\" dialog")));
      await session.step(77, "Then OK button in \"Rename dataconnection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Rename dataconnection\" dialog"), "disabled"));
      await session.step(78, "When user enters \"BDDConnRenamedEdge{time}\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, session.text("BDDConnRenamedEdge{time}"), el("Name input in \"Rename dataconnection\" dialog")));
      await session.step(79, "Then OK button in \"Rename dataconnection\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Rename dataconnection\" dialog"), "enabled"));
      await session.step(80, "When user clicks on CANCEL button in \"Rename dataconnection\" dialog", () => clickOn(page, el("CANCEL button in \"Rename dataconnection\" dialog")));
      await session.step(81, "Then the \"Rename dataconnection\" dialog should close", () => dialogCloses(page, "Rename dataconnection"));
      await session.step(82, "And 1 connection named \"BDDConnRenameEdge{time}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDDConnRenameEdge{time}")));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clone... offers a copy under another name and saves nothing when cancelled", async () => {
      await session.step(86, "When user picks \"Clone...\" from the context menu of Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Clone...", el(session.text("Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree"))));
      await session.step(87, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(88, "And Name input in \"Edit Connection\" dialog should have value \"Copy of BDD-Conn-Renamed-{run}\"", () => shouldHaveValue(page, el("Name input in \"Edit Connection\" dialog"), session.text("Copy of BDD-Conn-Renamed-{run}")));
      await session.step(89, "And Server input in \"Edit Connection\" dialog should have value \"db.datagrok.ai\"", () => shouldHaveValue(page, el("Server input in \"Edit Connection\" dialog"), "db.datagrok.ai"));
      await session.step(90, "When user clicks on CANCEL button in \"Edit Connection\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Connection\" dialog")));
      await session.step(91, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(92, "And 0 connections named \"Copy of BDD-Conn-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("Copy of BDD-Conn-Renamed-{run}")));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
