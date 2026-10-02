/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-edit.feature
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
import {clickOn, enterInto, isExpanded, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, connectionOnServer, connectionsOnServer, dialogCloses, noConnectionOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Editing a database connection", () => {
  const session = feature(test, "features/connections/connections-edit.feature", import.meta.url);
  test("Editing a database connection", {tag: ["@connections", "@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(18, "And no connection named \"BDD-Conn-Edit-{run}, BDD-Conn-Edited-{run}, BDD-Conn-Renamed-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-Edit-{run}, BDD-Conn-Edited-{run}, BDD-Conn-Renamed-{run}")));
    await session.step(19, "And a \"Postgres\" connection named \"BDD-Conn-Edit-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Edit-{run}")));
    await session.step(20, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await run.scenario("The Edit dialog shows what the connection was saved with", async () => {
      await session.step(23, "When user picks \"Edit...\" from the context menu of Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree"))));
      await session.step(24, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(25, "And Name input in \"Edit Connection\" dialog should have value \"BDD-Conn-Edit-{run}\"", () => shouldHaveValue(page, el("Name input in \"Edit Connection\" dialog"), session.text("BDD-Conn-Edit-{run}")));
      await session.step(26, "And Server input in \"Edit Connection\" dialog should have value \"db.datagrok.ai\"", () => shouldHaveValue(page, el("Server input in \"Edit Connection\" dialog"), "db.datagrok.ai"));
      await session.step(27, "And Db input in \"Edit Connection\" dialog should have value \"northwind\"", () => shouldHaveValue(page, el("Db input in \"Edit Connection\" dialog"), "northwind"));
      await session.step(28, "When user clicks on CANCEL button in \"Edit Connection\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Connection\" dialog")));
      await session.step(29, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Renaming through the Edit dialog", async () => {
      await session.step(33, "When user picks \"Edit...\" from the context menu of Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree"))));
      await session.step(34, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(35, "When user enters \"BDD-Conn-Edited-{run}\" into Name input in \"Edit Connection\" dialog", () => enterInto(page, session.text("BDD-Conn-Edited-{run}"), el("Name input in \"Edit Connection\" dialog")));
      await session.step(36, "And user clicks on OK button in \"Edit Connection\" dialog", () => clickOn(page, el("OK button in \"Edit Connection\" dialog")));
      await session.step(37, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(38, "And 1 connection named \"BDD-Conn-Edited-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Edited-{run}")));
      await session.step(39, "And 0 connections named \"BDD-Conn-Edit-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Edit-{run}")));
      await session.step(40, "And Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree")), "visible"));
      await session.step(41, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Renaming through Rename...", async () => {
      await session.step(44, "When user picks \"Rename...\" from the context menu of Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Databases---Postgres---BDD-Conn-Edited-{run} tree node inside browse tree"))));
      await session.step(45, "Then \"Rename dataconnection\" dialog should be visible", () => shouldBe(page, el("\"Rename dataconnection\" dialog"), "visible"));
      await session.step(46, "When user enters \"BDD-Conn-Renamed-{run}\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, session.text("BDD-Conn-Renamed-{run}"), el("Name input in \"Rename dataconnection\" dialog")));
      await session.step(47, "And user clicks on OK button in \"Rename dataconnection\" dialog", () => clickOn(page, el("OK button in \"Rename dataconnection\" dialog")));
      await session.step(48, "Then the \"Rename dataconnection\" dialog should close", () => dialogCloses(page, "Rename dataconnection"));
      await session.step(49, "And 1 connection named \"BDD-Conn-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Renamed-{run}")));
      await session.step(50, "And 0 connections named \"BDD-Conn-Edited-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Edited-{run}")));
      await session.step(51, "And Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree")), "visible"));
      await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Clone... offers a copy under another name and saves nothing when cancelled", async () => {
      await session.step(55, "When user picks \"Clone...\" from the context menu of Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Clone...", el(session.text("Databases---Postgres---BDD-Conn-Renamed-{run} tree node inside browse tree"))));
      await session.step(56, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(57, "And Name input in \"Edit Connection\" dialog should have value \"Copy of BDD-Conn-Renamed-{run}\"", () => shouldHaveValue(page, el("Name input in \"Edit Connection\" dialog"), session.text("Copy of BDD-Conn-Renamed-{run}")));
      await session.step(58, "And Server input in \"Edit Connection\" dialog should have value \"db.datagrok.ai\"", () => shouldHaveValue(page, el("Server input in \"Edit Connection\" dialog"), "db.datagrok.ai"));
      await session.step(59, "When user clicks on CANCEL button in \"Edit Connection\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Connection\" dialog")));
      await session.step(60, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(61, "And 0 connections named \"Copy of BDD-Conn-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("Copy of BDD-Conn-Renamed-{run}")));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
