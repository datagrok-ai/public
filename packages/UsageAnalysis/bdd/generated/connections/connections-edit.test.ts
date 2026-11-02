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
import {browsePanelOpen, connectionOnServer, connectionsOnServer, dialogCloses, noConnectionOnServer, refreshBrowse} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Editing a database connection", () => {
  const session = feature(test, "features/connections/connections-edit.feature", import.meta.url);
  test("Editing a database connection", {tag: ["@connections", "@journey", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(24, "And no connection named \"BDD-Conn-Edit-{run}, BDD-Conn-Edited-{run}, BDD-Conn-Renamed-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-Edit-{run}, BDD-Conn-Edited-{run}, BDD-Conn-Renamed-{run}")));
    await session.step(25, "And a \"Postgres\" connection named \"BDD-Conn-Edit-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Edit-{run}")));
    await session.step(26, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await run.scenario("The Edit dialog shows what the connection was saved with", async () => {
      await session.step(29, "When user picks \"Edit...\" from the context menu of Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree"))));
      await session.step(30, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(31, "And Name input in \"Edit Connection\" dialog should have value \"BDD-Conn-Edit-{run}\"", () => shouldHaveValue(page, el("Name input in \"Edit Connection\" dialog"), session.text("BDD-Conn-Edit-{run}")));
      await session.step(32, "And Server input in \"Edit Connection\" dialog should have value \"db.datagrok.ai\"", () => shouldHaveValue(page, el("Server input in \"Edit Connection\" dialog"), "db.datagrok.ai"));
      await session.step(33, "And Db input in \"Edit Connection\" dialog should have value \"northwind\"", () => shouldHaveValue(page, el("Db input in \"Edit Connection\" dialog"), "northwind"));
      await session.step(34, "When user clicks on CANCEL button in \"Edit Connection\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Connection\" dialog")));
      await session.step(35, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Renaming through the Edit dialog", async () => {
      await session.step(39, "When user picks \"Edit...\" from the context menu of Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree"))));
      await session.step(40, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(41, "When user enters \"BDD-Conn-Edited-{run}\" into Name input in \"Edit Connection\" dialog", () => enterInto(page, session.text("BDD-Conn-Edited-{run}"), el("Name input in \"Edit Connection\" dialog")));
      await session.step(42, "And user clicks on OK button in \"Edit Connection\" dialog", () => clickOn(page, el("OK button in \"Edit Connection\" dialog")));
      await session.step(43, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(44, "And 1 connection named \"BDD-Conn-Edited-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Edited-{run}")));
      await session.step(45, "And 0 connections named \"BDD-Conn-Edit-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Edit-{run}")));
      await session.step(46, "And Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree")), "visible"));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Renaming through Rename...", async () => {
      await session.step(50, "When user picks \"Rename...\" from the context menu of Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree"))));
      await session.step(51, "Then \"Rename dataconnection\" dialog should be visible", () => shouldBe(page, el("\"Rename dataconnection\" dialog"), "visible"));
      await session.step(52, "When user enters \"BDD-Conn-Renamed-{run}\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, session.text("BDD-Conn-Renamed-{run}"), el("Name input in \"Rename dataconnection\" dialog")));
      await session.step(53, "And user clicks on OK button in \"Rename dataconnection\" dialog", () => clickOn(page, el("OK button in \"Rename dataconnection\" dialog")));
      await session.step(54, "Then the \"Rename dataconnection\" dialog should close", () => dialogCloses(page, "Rename dataconnection"));
      await session.step(55, "And 1 connection named \"BDD-Conn-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Renamed-{run}")));
      await session.step(56, "And 0 connections named \"BDD-Conn-Edited-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Edited-{run}")));
      await session.step(57, "And Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree")), "visible"));
      await session.step(58, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Clone... offers a copy under another name and saves nothing when cancelled", async () => {
      await session.step(61, "When user picks \"Clone...\" from the context menu of Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Clone...", el(session.text("Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree"))));
      await session.step(62, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(63, "And Name input in \"Edit Connection\" dialog should have value \"Copy of BDD-Conn-Renamed-{run}\"", () => shouldHaveValue(page, el("Name input in \"Edit Connection\" dialog"), session.text("Copy of BDD-Conn-Renamed-{run}")));
      await session.step(64, "And Server input in \"Edit Connection\" dialog should have value \"db.datagrok.ai\"", () => shouldHaveValue(page, el("Server input in \"Edit Connection\" dialog"), "db.datagrok.ai"));
      await session.step(65, "When user clicks on CANCEL button in \"Edit Connection\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Connection\" dialog")));
      await session.step(66, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(67, "And 0 connections named \"Copy of BDD-Conn-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("Copy of BDD-Conn-Renamed-{run}")));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Edit dialog refuses a name of spaces and keeps the name on CANCEL", async () => {
      await session.step(72, "Given a \"Postgres\" connection named \"BDDConnEditEdge{time}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDDConnEditEdge{time}")));
      await session.step(73, "And no connection named \"BDDConnEditedEdge{time}\" is on the server", () => noConnectionOnServer(page, session.text("BDDConnEditedEdge{time}")));
      await session.step(74, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(75, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(76, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(77, "When user picks \"Edit...\" from the context menu of Databases---Postgres---BDDConnEditEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---BDDConnEditEdge{time} tree node inside browse tree"))));
      await session.step(78, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(79, "When user enters \"   \" into Name input in \"Edit Connection\" dialog", () => enterInto(page, "   ", el("Name input in \"Edit Connection\" dialog")));
      await session.step(80, "Then OK button in \"Edit Connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Edit Connection\" dialog"), "disabled"));
      await session.step(81, "When user enters \"BDDConnEditedEdge{time}\" into Name input in \"Edit Connection\" dialog", () => enterInto(page, session.text("BDDConnEditedEdge{time}"), el("Name input in \"Edit Connection\" dialog")));
      await session.step(82, "Then OK button in \"Edit Connection\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Edit Connection\" dialog"), "enabled"));
      await session.step(83, "When user clicks on CANCEL button in \"Edit Connection\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Connection\" dialog")));
      await session.step(84, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(85, "And 1 connection named \"BDDConnEditEdge{time}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDDConnEditEdge{time}")));
      await session.step(86, "And 0 connections named \"BDDConnEditedEdge{time}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDDConnEditedEdge{time}")));
      await session.step(87, "And Databases---Postgres---BDDConnEditEdge{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDDConnEditEdge{time} tree node inside browse tree")), "visible"));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Rename... refuses a name with no letter and keeps the name on CANCEL", async () => {
      await session.step(92, "Given a \"Postgres\" connection named \"BDDConnRenameEdge{time}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDDConnRenameEdge{time}")));
      await session.step(93, "And no connection named \"BDDConnRenamedEdge{time}\" is on the server", () => noConnectionOnServer(page, session.text("BDDConnRenamedEdge{time}")));
      await session.step(94, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(95, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(96, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(97, "When user picks \"Rename...\" from the context menu of Databases---Postgres---BDDConnRenameEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Databases---Postgres---BDDConnRenameEdge{time} tree node inside browse tree"))));
      await session.step(98, "Then \"Rename dataconnection\" dialog should be visible", () => shouldBe(page, el("\"Rename dataconnection\" dialog"), "visible"));
      await session.step(99, "When user enters \"2026\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, "2026", el("Name input in \"Rename dataconnection\" dialog")));
      await session.step(100, "Then OK button in \"Rename dataconnection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Rename dataconnection\" dialog"), "disabled"));
      await session.step(101, "When user enters \"BDDConnRenamedEdge{time}\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, session.text("BDDConnRenamedEdge{time}"), el("Name input in \"Rename dataconnection\" dialog")));
      await session.step(102, "Then OK button in \"Rename dataconnection\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Rename dataconnection\" dialog"), "enabled"));
      await session.step(103, "When user clicks on CANCEL button in \"Rename dataconnection\" dialog", () => clickOn(page, el("CANCEL button in \"Rename dataconnection\" dialog")));
      await session.step(104, "Then the \"Rename dataconnection\" dialog should close", () => dialogCloses(page, "Rename dataconnection"));
      await session.step(105, "And 1 connection named \"BDDConnRenameEdge{time}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDDConnRenameEdge{time}")));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
