/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-add.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [core#1176]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
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
import {hiddenProvidersShown} from '../../bindings/connections.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, followingShouldBe, isExpanded, selectIn, shouldBe, shouldBecomeVisibleWithin} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, connectionDataSource, connectionsOnServer, dialogCloses, noConnectionOnServer, standHasReachableConnection} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Adding a database connection", () => {
  const session = feature(test, "features/connections/connections-add.feature", import.meta.url);
  test("The New connection dialog of Postgres asks for its fields [provider=Postgres, node=Postgres]", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(32, "And the hidden providers of the Databases tree are shown", () => hiddenProvidersShown(page));
    await session.step(35, "When user picks \"New connection...\" from the context menu of Databases---Postgres tree node inside browse tree", () => pickFromContextMenu(page, "New connection...", el("Databases---Postgres tree node inside browse tree")));
    await session.step(36, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(37, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]), [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]);
    await session.step(44, "And OK button in \"Add new connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "disabled"));
    await session.step(45, "When user clicks on CANCEL button in \"Add new connection\" dialog", () => clickOn(page, el("CANCEL button in \"Add new connection\" dialog")));
    await session.step(46, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
  test("The New connection dialog of MS SQL asks for its fields [provider=MS SQL, node=MS-SQL]", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(32, "And the hidden providers of the Databases tree are shown", () => hiddenProvidersShown(page));
    await session.step(35, "When user picks \"New connection...\" from the context menu of Databases---MS-SQL tree node inside browse tree", () => pickFromContextMenu(page, "New connection...", el("Databases---MS-SQL tree node inside browse tree")));
    await session.step(36, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(37, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]), [["Name input in \"Add new connection\" dialog"],["Server input in \"Add new connection\" dialog"],["Db input in \"Add new connection\" dialog"],["Login input in \"Add new connection\" dialog"],["Password input in \"Add new connection\" dialog"],["TEST button in \"Add new connection\" dialog"]]);
    await session.step(44, "And OK button in \"Add new connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "disabled"));
    await session.step(45, "When user clicks on CANCEL button in \"Add new connection\" dialog", () => clickOn(page, el("CANCEL button in \"Add new connection\" dialog")));
    await session.step(46, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
  test("A connection string replaces the server fields", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(32, "And the hidden providers of the Databases tree are shown", () => hiddenProvidersShown(page));
    await session.step(55, "When user picks \"New connection...\" from the context menu of Databases---Postgres tree node inside browse tree", () => pickFromContextMenu(page, "New connection...", el("Databases---Postgres tree node inside browse tree")));
    await session.step(56, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(57, "When user selects \"Connection string\" in Configure input in \"Add new connection\" dialog", () => selectIn(page, "Connection string", el("Configure input in \"Add new connection\" dialog")));
    await session.step(58, "Then Conn-String input in \"Add new connection\" dialog should be visible", () => shouldBe(page, el("Conn-String input in \"Add new connection\" dialog"), "visible"));
    await session.step(59, "And Server input in \"Add new connection\" dialog should be hidden", () => shouldBe(page, el("Server input in \"Add new connection\" dialog"), "hidden"));
    await session.step(60, "When user clicks on CANCEL button in \"Add new connection\" dialog", () => clickOn(page, el("CANCEL button in \"Add new connection\" dialog")));
    await session.step(61, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(62, "And no errors should have been logged", () => noErrors(page));
  });
  test("A named connection is saved under its provider", {tag: ["@connections"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(32, "And the hidden providers of the Databases tree are shown", () => hiddenProvidersShown(page));
    await session.step(65, "Given no connection named \"BDD-Conn-Add-{run}\" is on the server", () => noConnectionOnServer(page, session.text("BDD-Conn-Add-{run}")));
    await session.step(66, "When user picks \"New connection...\" from the context menu of Databases---Postgres tree node inside browse tree", () => pickFromContextMenu(page, "New connection...", el("Databases---Postgres tree node inside browse tree")));
    await session.step(67, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(68, "And OK button in \"Add new connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "disabled"));
    await session.step(69, "When user enters \"BDD-Conn-Add-{run}\" into Name input in \"Add new connection\" dialog", () => enterInto(page, session.text("BDD-Conn-Add-{run}"), el("Name input in \"Add new connection\" dialog")));
    await session.step(70, "Then OK button in \"Add new connection\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "enabled"));
    await session.step(71, "When user enters \"db.datagrok.ai\" into Server input in \"Add new connection\" dialog", () => enterInto(page, "db.datagrok.ai", el("Server input in \"Add new connection\" dialog")));
    await session.step(72, "And user enters \"54322\" into Port input in \"Add new connection\" dialog", () => enterInto(page, "54322", el("Port input in \"Add new connection\" dialog")));
    await session.step(73, "And user enters \"northwind\" into Db input in \"Add new connection\" dialog", () => enterInto(page, "northwind", el("Db input in \"Add new connection\" dialog")));
    await session.step(74, "And user enters \"datagrok\" into Login input in \"Add new connection\" dialog", () => enterInto(page, "datagrok", el("Login input in \"Add new connection\" dialog")));
    await session.step(75, "And user clicks on OK button in \"Add new connection\" dialog", () => clickOn(page, el("OK button in \"Add new connection\" dialog")));
    await session.step(76, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(77, "And 1 connection named \"BDD-Conn-Add-{run}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDD-Conn-Add-{run}")));
    await session.step(78, "And the \"BDD-Conn-Add-{run}\" connection on the server should have the data source \"Postgres\"", () => connectionDataSource(page, session.text("BDD-Conn-Add-{run}"), "Postgres"));
    await session.step(79, "Given Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(80, "Then Databases---Postgres---BDD-Conn-Add-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Databases---Postgres---BDD-Conn-Add-{run} tree node inside browse tree")), "visible"));
    await session.step(81, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(82, "And no errors should have been logged", () => noErrors(page));
  });
  test("A name of spaces keeps OK disabled, and CANCEL with a name typed saves nothing", {tag: ["@connections", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(32, "And the hidden providers of the Databases tree are shown", () => hiddenProvidersShown(page));
    await session.step(86, "Given no connection named \"BDDConnSpaces{time}\" is on the server", () => noConnectionOnServer(page, session.text("BDDConnSpaces{time}")));
    await session.step(87, "When user picks \"New connection...\" from the context menu of Databases---Postgres tree node inside browse tree", () => pickFromContextMenu(page, "New connection...", el("Databases---Postgres tree node inside browse tree")));
    await session.step(88, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(89, "When user enters \"BDDConnSpaces{time}\" into Name input in \"Add new connection\" dialog", () => enterInto(page, session.text("BDDConnSpaces{time}"), el("Name input in \"Add new connection\" dialog")));
    await session.step(90, "And user enters \"db.datagrok.ai\" into Server input in \"Add new connection\" dialog", () => enterInto(page, "db.datagrok.ai", el("Server input in \"Add new connection\" dialog")));
    await session.step(91, "Then OK button in \"Add new connection\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "enabled"));
    await session.step(92, "When user enters \"   \" into Name input in \"Add new connection\" dialog", () => enterInto(page, "   ", el("Name input in \"Add new connection\" dialog")));
    await session.step(93, "Then OK button in \"Add new connection\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "disabled"));
    await session.step(94, "When user enters \"BDDConnSpaces{time}\" into Name input in \"Add new connection\" dialog", () => enterInto(page, session.text("BDDConnSpaces{time}"), el("Name input in \"Add new connection\" dialog")));
    await session.step(95, "Then OK button in \"Add new connection\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Add new connection\" dialog"), "enabled"));
    await session.step(96, "When user clicks on CANCEL button in \"Add new connection\" dialog", () => clickOn(page, el("CANCEL button in \"Add new connection\" dialog")));
    await session.step(97, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(98, "And 0 connections named \"BDDConnSpaces{time}\" should be on the server", () => connectionsOnServer(page, 0, session.text("BDDConnSpaces{time}")));
    await session.step(99, "And no errors should have been logged", () => noErrors(page));
  });
  test("A clone keeps the credentials and lists its schemas", {tag: ["@connections", "@edge-cases", "@known-failure:115", "@realizes:core#1176"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(32, "And the hidden providers of the Databases tree are shown", () => hiddenProvidersShown(page));
    await knownFailure(async () => {
      await session.step(105, "Given no connection named \"BDDConnClone{time}\" is on the server", () => noConnectionOnServer(page, session.text("BDDConnClone{time}")));
      await session.step(106, "And the stand has a reachable \"PostgresTest\" connection", () => standHasReachableConnection(page, "PostgresTest"));
      await session.step(107, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(108, "When user picks \"Clone...\" from the context menu of Databases---Postgres---NorthwindTest tree node inside browse tree", () => pickFromContextMenu(page, "Clone...", el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(109, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
      await session.step(110, "When user enters \"BDDConnClone{time}\" into Name input in \"Edit Connection\" dialog", () => enterInto(page, session.text("BDDConnClone{time}"), el("Name input in \"Edit Connection\" dialog")));
      await session.step(111, "And user clicks on OK button in \"Edit Connection\" dialog", () => clickOn(page, el("OK button in \"Edit Connection\" dialog")));
      await session.step(112, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
      await session.step(113, "And 1 connection named \"BDDConnClone{time}\" should be on the server", () => connectionsOnServer(page, 1, session.text("BDDConnClone{time}")));
      await session.step(114, "Given Databases---Postgres---BDDConnClone{time} tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("Databases---Postgres---BDDConnClone{time} tree node inside browse tree"))));
      await session.step(115, "Then Databases---Postgres---BDDConnClone{time}---Schemas tree node inside browse tree should become visible within 60 seconds", () => shouldBecomeVisibleWithin(page, el(session.text("Databases---Postgres---BDDConnClone{time}---Schemas tree node inside browse tree")), 60));
    }, 115);
  });
});
