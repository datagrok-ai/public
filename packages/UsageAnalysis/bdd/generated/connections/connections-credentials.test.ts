/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-credentials.feature
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
import {clickOn, enterInto, enterSecret, isExpanded, rightClickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {taskBarShown, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, connectionOnServer, connectionsOnServer, dialogCloses, noConnectionOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Connections that log in with real credentials", () => {
  const session = feature(test, "features/connections/connections-credentials.feature", import.meta.url);
  test("The right credentials pass the Edit dialog's TEST", {tag: ["@connections", "@needs-credentials"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(24, "Given a \"Postgres\" connection named \"BDD-Conn-Creds-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Creds-{run}")));
    await session.step(25, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(26, "When user right-clicks on Databases---Postgres---BDD-Conn-Creds-{run} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---Postgres---BDD-Conn-Creds-{run} tree node inside browse tree"))));
    await session.step(27, "And user picks \"Edit...\" from the open menu", () => pickFromOpenMenu(page, "Edit..."));
    await session.step(28, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
    await session.step(29, "When user enters the DG_PG_LOGIN secret into Login input in \"Edit Connection\" dialog", () => enterSecret(page, "DG_PG_LOGIN", el("Login input in \"Edit Connection\" dialog")));
    await session.step(30, "And user enters the DG_PG_PASSWORD secret into Password input in \"Edit Connection\" dialog", () => enterSecret(page, "DG_PG_PASSWORD", el("Password input in \"Edit Connection\" dialog")));
    await session.step(31, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(32, "When user clicks on TEST button in \"Edit Connection\" dialog", () => clickOn(page, el("TEST button in \"Edit Connection\" dialog")));
    await session.step(33, "Then the task bar should have shown \"Testing\"", () => taskBarShown(page, "Testing"));
    await session.step(34, "And the connection test should have ended on an info balloon containing \"connected successfully\"", () => connectionTestEnded(page, "info", "connected successfully"));
    await session.step(35, "When user clicks on OK button in \"Edit Connection\" dialog", () => clickOn(page, el("OK button in \"Edit Connection\" dialog")));
    await session.step(36, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
    await session.step(37, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(38, "When user right-clicks on Databases---Postgres---BDD-Conn-Creds-{run} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---Postgres---BDD-Conn-Creds-{run} tree node inside browse tree"))));
    await session.step(39, "And user picks \"Test connection\" from the open menu", () => pickFromOpenMenu(page, "Test connection"));
    await session.step(40, "Then the connection test should have ended on an info balloon containing \"connected successfully\"", () => connectionTestEnded(page, "info", "connected successfully"));
  });
  test("The external provider's connection is created from the dialog and passes its test", {tag: ["@connections", "@needs-credentials"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(43, "Given no connection named \"BDD-Conn-Ext-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-Ext-{run}")));
    await session.step(44, "When user right-clicks on Databases---Postgres tree node inside browse tree", () => rightClickOn(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(45, "And user picks \"New connection...\" from the open menu", () => pickFromOpenMenu(page, "New connection..."));
    await session.step(46, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(47, "When user enters \"BDD-Conn-Ext-{run}\" into Name input in \"Add new connection\" dialog", () => enterInto(page, session.text("BDD-Conn-Ext-{run}"), el("Name input in \"Add new connection\" dialog")));
    await session.step(48, "And user enters \"db.datagrok.ai\" into Server input in \"Add new connection\" dialog", () => enterInto(page, "db.datagrok.ai", el("Server input in \"Add new connection\" dialog")));
    await session.step(49, "And user enters \"54327\" into Port input in \"Add new connection\" dialog", () => enterInto(page, "54327", el("Port input in \"Add new connection\" dialog")));
    await session.step(50, "And user enters \"test\" into Db input in \"Add new connection\" dialog", () => enterInto(page, "test", el("Db input in \"Add new connection\" dialog")));
    await session.step(51, "And user enters the DG_PG_EXT_LOGIN secret into Login input in \"Add new connection\" dialog", () => enterSecret(page, "DG_PG_EXT_LOGIN", el("Login input in \"Add new connection\" dialog")));
    await session.step(52, "And user enters the DG_PG_EXT_PASSWORD secret into Password input in \"Add new connection\" dialog", () => enterSecret(page, "DG_PG_EXT_PASSWORD", el("Password input in \"Add new connection\" dialog")));
    await session.step(53, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(54, "When user clicks on TEST button in \"Add new connection\" dialog", () => clickOn(page, el("TEST button in \"Add new connection\" dialog")));
    await session.step(55, "Then the task bar should have shown \"Testing\"", () => taskBarShown(page, "Testing"));
    await session.step(56, "And the connection test should have ended on an info balloon containing \"connected successfully\"", () => connectionTestEnded(page, "info", "connected successfully"));
    await session.step(57, "When user clicks on OK button in \"Add new connection\" dialog", () => clickOn(page, el("OK button in \"Add new connection\" dialog")));
    await session.step(58, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(59, "And 1 connection named \"BDD-Conn-Ext-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Ext-{run}")));
  });
});
