/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-external-provider.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noTableOnConnection} from '../../bindings/connections.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, replaceCode, rightClickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeCurrentView, currentViewType, noQueryOnServer, queriesOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromOpenMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Writing to an external Postgres through queries", () => {
  const session = feature(test, "features/connections/connections-external-provider.feature", import.meta.url);
  test("Writing to an external Postgres through queries", {tag: ["@connections", "@full-stand", "@serial", "@needs-credentials", "@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "And no table \"bdd_tmp_{time}\" is in the database of the \"PostgreSQLDBTests\" connection", () => noTableOnConnection(page, session.text("bdd_tmp_{time}"), "PostgreSQLDBTests"));
    await session.step(34, "And no query named \"BDD-Conn-Ext-Create-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Conn-Ext-Create-{time}")));
    await session.step(35, "And no query named \"BDD-Conn-Ext-Insert-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Conn-Ext-Insert-{time}")));
    await session.step(36, "And no query named \"BDD-Conn-Ext-Update-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Conn-Ext-Update-{time}")));
    await session.step(37, "And no query named \"BDD-Conn-Ext-Select-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-Conn-Ext-Select-{time}")));
    await session.step(38, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(39, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await run.scenario("A CREATE TABLE query is saved and run", async () => {
      await session.step(42, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(43, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(44, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(45, "When user right-clicks on first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree", () => rightClickOn(page, el("first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree")));
      await session.step(46, "And user picks \"New Query...\" from the open menu", () => pickFromOpenMenu(page, "New Query..."));
      await session.step(47, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(48, "When user enters \"BDD-Conn-Ext-Create-{time}\" into Name input", () => enterInto(page, session.text("BDD-Conn-Ext-Create-{time}"), el("Name input")));
      await session.step(49, "And user replaces the code of code editor with \"create table bdd_tmp_{time} (id int, name varchar(50))\"", () => replaceCode(page, el("code editor"), session.text("create table bdd_tmp_{time} (id int, name varchar(50))")));
      await session.step(50, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(51, "Then 1 query named \"BDD-Conn-Ext-Create-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-Conn-Ext-Create-{time}")));
      await session.step(52, "When user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(53, "Then no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An INSERT query writes a row the SELECT query reads back", async () => {
      await session.step(57, "When user closes the current view", () => closeCurrentView(page));
      await session.step(58, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(59, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(60, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(61, "When user right-clicks on first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree", () => rightClickOn(page, el("first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree")));
      await session.step(62, "And user picks \"New Query...\" from the open menu", () => pickFromOpenMenu(page, "New Query..."));
      await session.step(63, "And user enters \"BDD-Conn-Ext-Insert-{time}\" into Name input", () => enterInto(page, session.text("BDD-Conn-Ext-Insert-{time}"), el("Name input")));
      await session.step(64, "And user replaces the code of code editor with \"insert into bdd_tmp_{time} (id, name) values (1, 'test')\"", () => replaceCode(page, el("code editor"), session.text("insert into bdd_tmp_{time} (id, name) values (1, 'test')")));
      await session.step(65, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(66, "Then 1 query named \"BDD-Conn-Ext-Insert-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-Conn-Ext-Insert-{time}")));
      await session.step(67, "When user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(68, "Then no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(69, "When user closes the current view", () => closeCurrentView(page));
      await session.step(70, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(71, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(72, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(73, "When user right-clicks on first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree", () => rightClickOn(page, el("first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree")));
      await session.step(74, "And user picks \"New Query...\" from the open menu", () => pickFromOpenMenu(page, "New Query..."));
      await session.step(75, "And user enters \"BDD-Conn-Ext-Select-{time}\" into Name input", () => enterInto(page, session.text("BDD-Conn-Ext-Select-{time}"), el("Name input")));
      await session.step(76, "And user replaces the code of code editor with \"select id, name from bdd_tmp_{time}\"", () => replaceCode(page, el("code editor"), session.text("select id, name from bdd_tmp_{time}")));
      await session.step(77, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(78, "And user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(79, "Then the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
      await session.step(80, "And the \"text of cell 1 of name\" reading of grid should be \"test\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "test"));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An UPDATE query changes the row", async () => {
      await session.step(84, "When user closes the current view", () => closeCurrentView(page));
      await session.step(85, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(86, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(87, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(88, "When user right-clicks on first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree", () => rightClickOn(page, el("first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree")));
      await session.step(89, "And user picks \"New Query...\" from the open menu", () => pickFromOpenMenu(page, "New Query..."));
      await session.step(90, "And user enters \"BDD-Conn-Ext-Update-{time}\" into Name input", () => enterInto(page, session.text("BDD-Conn-Ext-Update-{time}"), el("Name input")));
      await session.step(91, "And user replaces the code of code editor with \"update bdd_tmp_{time} set name = 'bdd' where id = 1\"", () => replaceCode(page, el("code editor"), session.text("update bdd_tmp_{time} set name = 'bdd' where id = 1")));
      await session.step(92, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(93, "Then 1 query named \"BDD-Conn-Ext-Update-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-Conn-Ext-Update-{time}")));
      await session.step(94, "When user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(95, "Then no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(96, "When user closes the current view", () => closeCurrentView(page));
      await session.step(97, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(98, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(99, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(100, "When user right-clicks on first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree", () => rightClickOn(page, el("first Databases---Postgres---PostgreSQLDBTests tree node inside browse tree")));
      await session.step(101, "And user picks \"New Query...\" from the open menu", () => pickFromOpenMenu(page, "New Query..."));
      await session.step(102, "And user replaces the code of code editor with \"select id, name from bdd_tmp_{time}\"", () => replaceCode(page, el("code editor"), session.text("select id, name from bdd_tmp_{time}")));
      await session.step(103, "And user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(104, "Then the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
      await session.step(105, "And the \"text of cell 1 of name\" reading of grid should be \"bdd\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "bdd"));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
