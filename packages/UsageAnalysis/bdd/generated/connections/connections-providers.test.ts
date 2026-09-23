/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-providers.feature
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
import {clickOn, enterInto, isExpanded, rightClickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {taskBarShown, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, connectionOnServer, connectionsOnServer, dialogCloses, noConnectionOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Editing connections of the other providers", () => {
  const session = feature(test, "features/connections/connections-providers.feature", import.meta.url);
  test("A Oracle connection is renamed, and its test without credentials fails [provider=Oracle, node=Oracle, short=Oracle]", {tag: ["@connections", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(26, "Given no connection named \"BDD-Conn-Oracle-{run}, BDD-Conn-Oracle-Renamed-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-Oracle-{run}, BDD-Conn-Oracle-Renamed-{run}")));
    await session.step(27, "And a \"Oracle\" connection named \"BDD-Conn-Oracle-{run}\" is on the server", () => connectionOnServer(page, "Oracle", session.text("BDD-Conn-Oracle-{run}")));
    await session.step(28, "And Databases---Oracle tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Oracle tree node inside browse tree")));
    await session.step(29, "When user right-clicks on Databases---Oracle---BDD-Conn-Oracle-{run} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---Oracle---BDD-Conn-Oracle-{run} tree node inside browse tree"))));
    await session.step(30, "And user picks \"Rename...\" from the open menu", () => pickFromOpenMenu(page, "Rename..."));
    await session.step(31, "Then \"Rename dataconnection\" dialog should be visible", () => shouldBe(page, el("\"Rename dataconnection\" dialog"), "visible"));
    await session.step(32, "When user enters \"BDD-Conn-Oracle-Renamed-{run}\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, session.text("BDD-Conn-Oracle-Renamed-{run}"), el("Name input in \"Rename dataconnection\" dialog")));
    await session.step(33, "And user clicks on OK button in \"Rename dataconnection\" dialog", () => clickOn(page, el("OK button in \"Rename dataconnection\" dialog")));
    await session.step(34, "Then the \"Rename dataconnection\" dialog should close", () => dialogCloses(page, "Rename dataconnection"));
    await session.step(35, "And 1 connection named \"BDD-Conn-Oracle-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Oracle-Renamed-{run}")));
    await session.step(36, "And \"BDD-Conn-Oracle-Renamed-{run}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDD-Conn-Oracle-Renamed-{run}\" tree node inside browse tree")), "visible"));
    await session.step(37, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(38, "When user right-clicks on \"BDD-Conn-Oracle-Renamed-{run}\" tree node inside browse tree", () => rightClickOn(page, el(session.text("\"BDD-Conn-Oracle-Renamed-{run}\" tree node inside browse tree"))));
    await session.step(39, "And user picks \"Test connection\" from the open menu", () => pickFromOpenMenu(page, "Test connection"));
    await session.step(40, "Then the task bar should have shown \"Testing\"", () => taskBarShown(page, "Testing"));
    await session.step(41, "And the connection test should have ended on an error balloon containing \"BDD-Conn-Oracle-Renamed-{run}\"", () => connectionTestEnded(page, "error", session.text("BDD-Conn-Oracle-Renamed-{run}")));
  });
  test("A MariaDB connection is renamed, and its test without credentials fails [provider=MariaDB, node=MariaDB, short=MariaDB]", {tag: ["@connections", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(26, "Given no connection named \"BDD-Conn-MariaDB-{run}, BDD-Conn-MariaDB-Renamed-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-MariaDB-{run}, BDD-Conn-MariaDB-Renamed-{run}")));
    await session.step(27, "And a \"MariaDB\" connection named \"BDD-Conn-MariaDB-{run}\" is on the server", () => connectionOnServer(page, "MariaDB", session.text("BDD-Conn-MariaDB-{run}")));
    await session.step(28, "And Databases---MariaDB tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---MariaDB tree node inside browse tree")));
    await session.step(29, "When user right-clicks on Databases---MariaDB---BDD-Conn-MariaDB-{run} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---MariaDB---BDD-Conn-MariaDB-{run} tree node inside browse tree"))));
    await session.step(30, "And user picks \"Rename...\" from the open menu", () => pickFromOpenMenu(page, "Rename..."));
    await session.step(31, "Then \"Rename dataconnection\" dialog should be visible", () => shouldBe(page, el("\"Rename dataconnection\" dialog"), "visible"));
    await session.step(32, "When user enters \"BDD-Conn-MariaDB-Renamed-{run}\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, session.text("BDD-Conn-MariaDB-Renamed-{run}"), el("Name input in \"Rename dataconnection\" dialog")));
    await session.step(33, "And user clicks on OK button in \"Rename dataconnection\" dialog", () => clickOn(page, el("OK button in \"Rename dataconnection\" dialog")));
    await session.step(34, "Then the \"Rename dataconnection\" dialog should close", () => dialogCloses(page, "Rename dataconnection"));
    await session.step(35, "And 1 connection named \"BDD-Conn-MariaDB-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-MariaDB-Renamed-{run}")));
    await session.step(36, "And \"BDD-Conn-MariaDB-Renamed-{run}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDD-Conn-MariaDB-Renamed-{run}\" tree node inside browse tree")), "visible"));
    await session.step(37, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(38, "When user right-clicks on \"BDD-Conn-MariaDB-Renamed-{run}\" tree node inside browse tree", () => rightClickOn(page, el(session.text("\"BDD-Conn-MariaDB-Renamed-{run}\" tree node inside browse tree"))));
    await session.step(39, "And user picks \"Test connection\" from the open menu", () => pickFromOpenMenu(page, "Test connection"));
    await session.step(40, "Then the task bar should have shown \"Testing\"", () => taskBarShown(page, "Testing"));
    await session.step(41, "And the connection test should have ended on an error balloon containing \"BDD-Conn-MariaDB-Renamed-{run}\"", () => connectionTestEnded(page, "error", session.text("BDD-Conn-MariaDB-Renamed-{run}")));
  });
  test("A MySQL connection is renamed, and its test without credentials fails [provider=MySQL, node=MySQL, short=MySQL]", {tag: ["@connections", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(26, "Given no connection named \"BDD-Conn-MySQL-{run}, BDD-Conn-MySQL-Renamed-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-MySQL-{run}, BDD-Conn-MySQL-Renamed-{run}")));
    await session.step(27, "And a \"MySQL\" connection named \"BDD-Conn-MySQL-{run}\" is on the server", () => connectionOnServer(page, "MySQL", session.text("BDD-Conn-MySQL-{run}")));
    await session.step(28, "And Databases---MySQL tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---MySQL tree node inside browse tree")));
    await session.step(29, "When user right-clicks on Databases---MySQL---BDD-Conn-MySQL-{run} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---MySQL---BDD-Conn-MySQL-{run} tree node inside browse tree"))));
    await session.step(30, "And user picks \"Rename...\" from the open menu", () => pickFromOpenMenu(page, "Rename..."));
    await session.step(31, "Then \"Rename dataconnection\" dialog should be visible", () => shouldBe(page, el("\"Rename dataconnection\" dialog"), "visible"));
    await session.step(32, "When user enters \"BDD-Conn-MySQL-Renamed-{run}\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, session.text("BDD-Conn-MySQL-Renamed-{run}"), el("Name input in \"Rename dataconnection\" dialog")));
    await session.step(33, "And user clicks on OK button in \"Rename dataconnection\" dialog", () => clickOn(page, el("OK button in \"Rename dataconnection\" dialog")));
    await session.step(34, "Then the \"Rename dataconnection\" dialog should close", () => dialogCloses(page, "Rename dataconnection"));
    await session.step(35, "And 1 connection named \"BDD-Conn-MySQL-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-MySQL-Renamed-{run}")));
    await session.step(36, "And \"BDD-Conn-MySQL-Renamed-{run}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDD-Conn-MySQL-Renamed-{run}\" tree node inside browse tree")), "visible"));
    await session.step(37, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(38, "When user right-clicks on \"BDD-Conn-MySQL-Renamed-{run}\" tree node inside browse tree", () => rightClickOn(page, el(session.text("\"BDD-Conn-MySQL-Renamed-{run}\" tree node inside browse tree"))));
    await session.step(39, "And user picks \"Test connection\" from the open menu", () => pickFromOpenMenu(page, "Test connection"));
    await session.step(40, "Then the task bar should have shown \"Testing\"", () => taskBarShown(page, "Testing"));
    await session.step(41, "And the connection test should have ended on an error balloon containing \"BDD-Conn-MySQL-Renamed-{run}\"", () => connectionTestEnded(page, "error", session.text("BDD-Conn-MySQL-Renamed-{run}")));
  });
  test("A MS SQL connection is renamed, and its test without credentials fails [provider=MS SQL, node=MS-SQL, short=MSSQL]", {tag: ["@connections", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(26, "Given no connection named \"BDD-Conn-MSSQL-{run}, BDD-Conn-MSSQL-Renamed-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-MSSQL-{run}, BDD-Conn-MSSQL-Renamed-{run}")));
    await session.step(27, "And a \"MS SQL\" connection named \"BDD-Conn-MSSQL-{run}\" is on the server", () => connectionOnServer(page, "MS SQL", session.text("BDD-Conn-MSSQL-{run}")));
    await session.step(28, "And Databases---MS-SQL tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---MS-SQL tree node inside browse tree")));
    await session.step(29, "When user right-clicks on Databases---MS-SQL---BDD-Conn-MSSQL-{run} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---MS-SQL---BDD-Conn-MSSQL-{run} tree node inside browse tree"))));
    await session.step(30, "And user picks \"Rename...\" from the open menu", () => pickFromOpenMenu(page, "Rename..."));
    await session.step(31, "Then \"Rename dataconnection\" dialog should be visible", () => shouldBe(page, el("\"Rename dataconnection\" dialog"), "visible"));
    await session.step(32, "When user enters \"BDD-Conn-MSSQL-Renamed-{run}\" into Name input in \"Rename dataconnection\" dialog", () => enterInto(page, session.text("BDD-Conn-MSSQL-Renamed-{run}"), el("Name input in \"Rename dataconnection\" dialog")));
    await session.step(33, "And user clicks on OK button in \"Rename dataconnection\" dialog", () => clickOn(page, el("OK button in \"Rename dataconnection\" dialog")));
    await session.step(34, "Then the \"Rename dataconnection\" dialog should close", () => dialogCloses(page, "Rename dataconnection"));
    await session.step(35, "And 1 connection named \"BDD-Conn-MSSQL-Renamed-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-MSSQL-Renamed-{run}")));
    await session.step(36, "And \"BDD-Conn-MSSQL-Renamed-{run}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDD-Conn-MSSQL-Renamed-{run}\" tree node inside browse tree")), "visible"));
    await session.step(37, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(38, "When user right-clicks on \"BDD-Conn-MSSQL-Renamed-{run}\" tree node inside browse tree", () => rightClickOn(page, el(session.text("\"BDD-Conn-MSSQL-Renamed-{run}\" tree node inside browse tree"))));
    await session.step(39, "And user picks \"Test connection\" from the open menu", () => pickFromOpenMenu(page, "Test connection"));
    await session.step(40, "Then the task bar should have shown \"Testing\"", () => taskBarShown(page, "Testing"));
    await session.step(41, "And the connection test should have ended on an error balloon containing \"BDD-Conn-MSSQL-Renamed-{run}\"", () => connectionTestEnded(page, "error", session.text("BDD-Conn-MSSQL-Renamed-{run}")));
  });
});
