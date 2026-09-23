/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-edit.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {connectionTestEnded} from '../../bindings/connections.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, rightClickOn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {taskBarShown, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, connectionOnServer, connectionsOnServer, dialogCloses, noConnectionOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Editing a database connection", () => {
  const session = feature(test, "features/connections/connections-edit.feature", import.meta.url);
  test("Editing a database connection", {tag: ["@connections", "@journey", "@slow", "@full-stand"]}, async ({browser}) => {
    test.slow();
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(31, "And no connection named \"BDD-Conn-Edit-{run}, BDD-Conn-Edited-{run}, BDD-Conn-Renamed-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-Edit-{run}, BDD-Conn-Edited-{run}, BDD-Conn-Renamed-{run}")));
    await session.step(32, "And a \"Postgres\" connection named \"BDD-Conn-Edit-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Edit-{run}")));
    await session.step(33, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await run.scenario("The Edit dialog shows what the connection was saved with", async () => {
      await session.step(36, "When user right-clicks on Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree"))));
      await session.step(37, "And user picks \"Edit...\" from the open menu", () => pickFromOpenMenu(page, "Edit..."));
      await session.step(38, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(39, "And Name input in \"Edit Connection\" dialog should have value \"BDD-Conn-Edit-{run}\"", () => shouldHaveValue(page, el("Name input in \"Edit Connection\" dialog"), session.text("BDD-Conn-Edit-{run}")));
      await session.step(40, "And Server input in \"Edit Connection\" dialog should have value \"db.datagrok.ai\"", () => shouldHaveValue(page, el("Server input in \"Edit Connection\" dialog"), "db.datagrok.ai"));
      await session.step(41, "And Db input in \"Edit Connection\" dialog should have value \"northwind\"", () => shouldHaveValue(page, el("Db input in \"Edit Connection\" dialog"), "northwind"));
      await session.step(42, "When user clicks on CANCEL button in \"Edit Connection\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Connection\" dialog")));
      await session.step(43, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Renaming through the Edit dialog", async () => {
      await session.step(47, "When user right-clicks on Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---Postgres---BDD-Conn-Edit-{run} tree node inside browse tree"))));
      await session.step(48, "And user picks \"Edit...\" from the open menu", () => pickFromOpenMenu(page, "Edit..."));
      await session.step(49, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(50, "When user enters \"BDD-Conn-Edited-{run}\" into Name input in \"Edit Connection\" dialog", () => enterInto(page, session.text("BDD-Conn-Edited-{run}"), el("Name input in \"Edit Connection\" dialog")));
      await session.step(51, "And user clicks on OK button in \"Edit Connection\" dialog", () => clickOn(page, el("OK button in \"Edit Connection\" dialog")));
      await session.step(52, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(53, "And 1 connection named \"BDD-Conn-Edited-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Edited-{run}")));
      await session.step(54, "And 0 connections named \"BDD-Conn-Edit-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Edit-{run}")));
      await session.step(55, "And \"BDD-Conn-Edited-{run}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDD-Conn-Edited-{run}\" tree node inside browse tree")), "visible"));
      await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Renaming through Rename...", async () => {
      await session.step(59, "When user right-clicks on \"BDD-Conn-Edited-{run}\" tree node inside browse tree", () => rightClickOn(page, el(session.text("\"BDD-Conn-Edited-{run}\" tree node inside browse tree"))));
      await session.step(60, "And user picks \"Rename...\" from the open menu", () => pickFromOpenMenu(page, "Rename..."));
      await session.step(61, "Then \"Rename dataconnection\" dialog should be visible", () => shouldBe(page, el("\"Rename dataconnection\" dialog"), "visible"));
      await session.step(62, "When user enters \"BDD-Conn-Renamed-{run}\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, session.text("BDD-Conn-Renamed-{run}"), el("Name input in \"Rename dataconnection\" dialog")));
      await session.step(63, "And user clicks on OK button in \"Rename dataconnection\" dialog", () => clickOn(page, el("OK button in \"Rename dataconnection\" dialog")));
      await session.step(64, "Then the \"Rename dataconnection\" dialog should close", () => dialogCloses(page, "Rename dataconnection"));
      await session.step(65, "And 1 connection named \"BDD-Conn-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Renamed-{run}")));
      await session.step(66, "And 0 connections named \"BDD-Conn-Edited-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDD-Conn-Edited-{run}")));
      await session.step(67, "And \"BDD-Conn-Renamed-{run}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDD-Conn-Renamed-{run}\" tree node inside browse tree")), "visible"));
      await session.step(68, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Clone... offers a copy under another name and saves nothing when cancelled", async () => {
      await session.step(71, "When user right-clicks on \"BDD-Conn-Renamed-{run}\" tree node inside browse tree", () => rightClickOn(page, el(session.text("\"BDD-Conn-Renamed-{run}\" tree node inside browse tree"))));
      await session.step(72, "And user picks \"Clone...\" from the open menu", () => pickFromOpenMenu(page, "Clone..."));
      await session.step(73, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(74, "And Name input in \"Edit Connection\" dialog should have value \"Copy of BDD-Conn-Renamed-{run}\"", () => shouldHaveValue(page, el("Name input in \"Edit Connection\" dialog"), session.text("Copy of BDD-Conn-Renamed-{run}")));
      await session.step(75, "And Server input in \"Edit Connection\" dialog should have value \"db.datagrok.ai\"", () => shouldHaveValue(page, el("Server input in \"Edit Connection\" dialog"), "db.datagrok.ai"));
      await session.step(76, "When user clicks on CANCEL button in \"Edit Connection\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Connection\" dialog")));
      await session.step(77, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(78, "And 0 connections named \"Copy of BDD-Conn-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 0, session.text("Copy of BDD-Conn-Renamed-{run}")));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Wrong credentials make Test connection fail", async () => {
      await session.step(85, "When user right-clicks on \"BDD-Conn-Renamed-{run}\" tree node inside browse tree", () => rightClickOn(page, el(session.text("\"BDD-Conn-Renamed-{run}\" tree node inside browse tree"))));
      await session.step(86, "And user picks \"Edit...\" from the open menu", () => pickFromOpenMenu(page, "Edit..."));
      await session.step(87, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(88, "When user enters \"bdd_nobody\" into Login input in \"Edit Connection\" dialog", () => enterInto(page, "bdd_nobody", el("Login input in \"Edit Connection\" dialog")));
      await session.step(89, "And user enters \"wrong\" into Password input in \"Edit Connection\" dialog", () => enterInto(page, "wrong", el("Password input in \"Edit Connection\" dialog")));
      await session.step(90, "And user clicks on OK button in \"Edit Connection\" dialog", () => clickOn(page, el("OK button in \"Edit Connection\" dialog")));
      await session.step(91, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(92, "Given user watches the task bar", () => watchTaskBar(page));
      await session.step(93, "When user right-clicks on \"BDD-Conn-Renamed-{run}\" tree node inside browse tree", () => rightClickOn(page, el(session.text("\"BDD-Conn-Renamed-{run}\" tree node inside browse tree"))));
      await session.step(94, "And user picks \"Test connection\" from the open menu", () => pickFromOpenMenu(page, "Test connection"));
      await session.step(95, "Then the task bar should have shown \"Testing\"", () => taskBarShown(page, "Testing"));
      await session.step(96, "And the connection test should have ended on an error balloon containing \"password authentication failed\"", () => connectionTestEnded(page, "error", "password authentication failed"));
    });
    run.finish();
  });
});
