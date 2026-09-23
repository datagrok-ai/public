/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-add.feature
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
import {clickOn, enterInto, enterSecret, followingShouldBe, isExpanded, rightClickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {taskBarShown, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, connectionDataSource, connectionsOnServer, dialogCloses, noConnectionOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Adding a database connection", () => {
  const session = feature(test, "features/connections/connections-add.feature", import.meta.url);
  test("The New connection dialog of Postgres asks for its fields [provider=Postgres, node=Postgres]", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(29, "When user right-clicks on Databases---Postgres tree node inside browse tree", () => rightClickOn(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(30, "And user picks \"New connection...\" from the open menu", () => pickFromOpenMenu(page, "New connection..."));
    await session.step(31, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(32, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]), [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]);
    await session.step(39, "And OK button in \"Add new connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "disabled"));
    await session.step(40, "When user clicks on CANCEL button in \"Add new connection\" dialog", () => clickOn(page, el("CANCEL button in \"Add new connection\" dialog")));
    await session.step(41, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
  test("The New connection dialog of MS SQL asks for its fields [provider=MS SQL, node=MS-SQL]", {tag: ["@connections", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(29, "When user right-clicks on Databases---MS-SQL tree node inside browse tree", () => rightClickOn(page, el("Databases---MS-SQL tree node inside browse tree")));
    await session.step(30, "And user picks \"New connection...\" from the open menu", () => pickFromOpenMenu(page, "New connection..."));
    await session.step(31, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(32, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]), [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]);
    await session.step(39, "And OK button in \"Add new connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "disabled"));
    await session.step(40, "When user clicks on CANCEL button in \"Add new connection\" dialog", () => clickOn(page, el("CANCEL button in \"Add new connection\" dialog")));
    await session.step(41, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
  test("The New connection dialog of Oracle asks for its fields [provider=Oracle, node=Oracle]", {tag: ["@connections", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(29, "When user right-clicks on Databases---Oracle tree node inside browse tree", () => rightClickOn(page, el("Databases---Oracle tree node inside browse tree")));
    await session.step(30, "And user picks \"New connection...\" from the open menu", () => pickFromOpenMenu(page, "New connection..."));
    await session.step(31, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(32, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]), [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]);
    await session.step(39, "And OK button in \"Add new connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "disabled"));
    await session.step(40, "When user clicks on CANCEL button in \"Add new connection\" dialog", () => clickOn(page, el("CANCEL button in \"Add new connection\" dialog")));
    await session.step(41, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
  test("The New connection dialog of MySQL asks for its fields [provider=MySQL, node=MySQL]", {tag: ["@connections", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(29, "When user right-clicks on Databases---MySQL tree node inside browse tree", () => rightClickOn(page, el("Databases---MySQL tree node inside browse tree")));
    await session.step(30, "And user picks \"New connection...\" from the open menu", () => pickFromOpenMenu(page, "New connection..."));
    await session.step(31, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(32, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]), [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]);
    await session.step(39, "And OK button in \"Add new connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "disabled"));
    await session.step(40, "When user clicks on CANCEL button in \"Add new connection\" dialog", () => clickOn(page, el("CANCEL button in \"Add new connection\" dialog")));
    await session.step(41, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
  test("The New connection dialog of MariaDB asks for its fields [provider=MariaDB, node=MariaDB]", {tag: ["@connections", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(29, "When user right-clicks on Databases---MariaDB tree node inside browse tree", () => rightClickOn(page, el("Databases---MariaDB tree node inside browse tree")));
    await session.step(30, "And user picks \"New connection...\" from the open menu", () => pickFromOpenMenu(page, "New connection..."));
    await session.step(31, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(32, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]), [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]);
    await session.step(39, "And OK button in \"Add new connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "disabled"));
    await session.step(40, "When user clicks on CANCEL button in \"Add new connection\" dialog", () => clickOn(page, el("CANCEL button in \"Add new connection\" dialog")));
    await session.step(41, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
  test("A connection string replaces the server fields", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(57, "When user right-clicks on Databases---Postgres tree node inside browse tree", () => rightClickOn(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(58, "And user picks \"New connection...\" from the open menu", () => pickFromOpenMenu(page, "New connection..."));
    await session.step(59, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(60, "When user selects \"Connection string\" in Configure input in \"Add new connection\" dialog", () => selectIn(page, "Connection string", el("Configure input in \"Add new connection\" dialog")));
    await session.step(61, "Then Conn-String input in \"Add new connection\" dialog should be visible", () => shouldBe(page, el("Conn-String input in \"Add new connection\" dialog"), "visible"));
    await session.step(62, "And Server input in \"Add new connection\" dialog should be hidden", () => shouldBe(page, el("Server input in \"Add new connection\" dialog"), "hidden"));
    await session.step(63, "When user clicks on CANCEL button in \"Add new connection\" dialog", () => clickOn(page, el("CANCEL button in \"Add new connection\" dialog")));
    await session.step(64, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(65, "And no errors should have been logged", () => noErrors(page));
  });
  test("A connection without a password is saved, and its TEST fails for the password", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(68, "Given no connection named \"BDD-Conn-Add-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-Add-{run}")));
    await session.step(69, "When user right-clicks on Databases---Postgres tree node inside browse tree", () => rightClickOn(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(70, "And user picks \"New connection...\" from the open menu", () => pickFromOpenMenu(page, "New connection..."));
    await session.step(71, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(72, "And OK button in \"Add new connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "disabled"));
    await session.step(73, "When user enters \"BDD-Conn-Add-{run}\" into Name input in \"Add new connection\" dialog", () => enterInto(page, session.text("BDD-Conn-Add-{run}"), el("Name input in \"Add new connection\" dialog")));
    await session.step(74, "Then OK button in \"Add new connection\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "enabled"));
    await session.step(75, "When user enters \"db.datagrok.ai\" into Server input in \"Add new connection\" dialog", () => enterInto(page, "db.datagrok.ai", el("Server input in \"Add new connection\" dialog")));
    await session.step(76, "And user enters \"54322\" into Port input in \"Add new connection\" dialog", () => enterInto(page, "54322", el("Port input in \"Add new connection\" dialog")));
    await session.step(77, "And user enters \"northwind\" into Db input in \"Add new connection\" dialog", () => enterInto(page, "northwind", el("Db input in \"Add new connection\" dialog")));
    await session.step(78, "And user enters \"datagrok\" into Login input in \"Add new connection\" dialog", () => enterInto(page, "datagrok", el("Login input in \"Add new connection\" dialog")));
    await session.step(79, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(80, "When user clicks on TEST button in \"Add new connection\" dialog", () => clickOn(page, el("TEST button in \"Add new connection\" dialog")));
    await session.step(81, "Then the task bar should have shown \"Testing\"", () => taskBarShown(page, "Testing"));
    await session.step(82, "And the connection test should have ended on an error balloon containing \"failed to connect\"", () => connectionTestEnded(page, "error", "failed to connect"));
    await session.step(83, "When user clicks on OK button in \"Add new connection\" dialog", () => clickOn(page, el("OK button in \"Add new connection\" dialog")));
    await session.step(84, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(85, "And 1 connection named \"BDD-Conn-Add-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Add-{run}")));
    await session.step(86, "And the \"BDD-Conn-Add-{run}\" connection on the server should have the data source \"Postgres\"", () => connectionDataSource(page, session.text("BDD-Conn-Add-{run}"), "Postgres"));
    await session.step(87, "Given Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(88, "Then Databases---Postgres---BDD-Conn-Add-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Add-{run} tree node inside browse tree")), "visible"));
  });
  test("A connection with the right credentials passes its TEST and is saved", {tag: ["@connections", "@needs-credentials"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(92, "Given no connection named \"BDD-Conn-Add-Ok-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-Add-Ok-{run}")));
    await session.step(93, "When user right-clicks on Databases---Postgres tree node inside browse tree", () => rightClickOn(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(94, "And user picks \"New connection...\" from the open menu", () => pickFromOpenMenu(page, "New connection..."));
    await session.step(95, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(96, "When user enters \"BDD-Conn-Add-Ok-{run}\" into Name input in \"Add new connection\" dialog", () => enterInto(page, session.text("BDD-Conn-Add-Ok-{run}"), el("Name input in \"Add new connection\" dialog")));
    await session.step(97, "And user enters \"db.datagrok.ai\" into Server input in \"Add new connection\" dialog", () => enterInto(page, "db.datagrok.ai", el("Server input in \"Add new connection\" dialog")));
    await session.step(98, "And user enters \"54322\" into Port input in \"Add new connection\" dialog", () => enterInto(page, "54322", el("Port input in \"Add new connection\" dialog")));
    await session.step(99, "And user enters \"northwind\" into Db input in \"Add new connection\" dialog", () => enterInto(page, "northwind", el("Db input in \"Add new connection\" dialog")));
    await session.step(100, "And user enters \"datagrok\" into Login input in \"Add new connection\" dialog", () => enterInto(page, "datagrok", el("Login input in \"Add new connection\" dialog")));
    await session.step(101, "And user enters the DG_PG_PASSWORD secret into Password input in \"Add new connection\" dialog", () => enterSecret(page, "DG_PG_PASSWORD", el("Password input in \"Add new connection\" dialog")));
    await session.step(102, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(103, "When user clicks on TEST button in \"Add new connection\" dialog", () => clickOn(page, el("TEST button in \"Add new connection\" dialog")));
    await session.step(104, "Then the task bar should have shown \"Testing\"", () => taskBarShown(page, "Testing"));
    await session.step(105, "And the connection test should have ended on an info balloon containing \"connected successfully\"", () => connectionTestEnded(page, "info", "connected successfully"));
    await session.step(106, "When user clicks on OK button in \"Add new connection\" dialog", () => clickOn(page, el("OK button in \"Add new connection\" dialog")));
    await session.step(107, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(108, "And 1 connection named \"BDD-Conn-Add-Ok-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Add-Ok-{run}")));
    await session.step(109, "And the \"BDD-Conn-Add-Ok-{run}\" connection on the server should have the data source \"Postgres\"", () => connectionDataSource(page, session.text("BDD-Conn-Add-Ok-{run}"), "Postgres"));
    await session.step(110, "Given Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(111, "Then Databases---Postgres---BDD-Conn-Add-Ok-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Add-Ok-{run} tree node inside browse tree")), "visible"));
    await session.step(112, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(113, "And no errors should have been logged", () => noErrors(page));
  });
});
