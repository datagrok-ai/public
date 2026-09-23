/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-identifiers.feature
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
import {clickOn, enterInto, enterSecret, isExpanded, rightClickOn, selectIn, shouldBe, shouldOffer} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, connectionOnServer, contextPanelOpen, dialogCloses, reloadPage, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloon, noBalloons, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Identifiers configured on a connection", () => {
  const session = feature(test, "features/connections/connections-identifiers.feature", import.meta.url);
  test("Identifiers configured on a connection", {tag: ["@connections", "@needs-credentials", "@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(29, "And a \"Postgres\" connection named \"BDD-Conn-Ident-{time}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Ident-{time}")));
    await session.step(30, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(31, "When user right-clicks on Databases---Postgres---BDD-Conn-Ident-{time} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---Postgres---BDD-Conn-Ident-{time} tree node inside browse tree"))));
    await session.step(32, "And user picks \"Edit...\" from the open menu", () => pickFromOpenMenu(page, "Edit..."));
    await session.step(33, "And user enters the DG_PG_LOGIN secret into Login input in \"Edit Connection\" dialog", () => enterSecret(page, "DG_PG_LOGIN", el("Login input in \"Edit Connection\" dialog")));
    await session.step(34, "And user enters the DG_PG_PASSWORD secret into Password input in \"Edit Connection\" dialog", () => enterSecret(page, "DG_PG_PASSWORD", el("Password input in \"Edit Connection\" dialog")));
    await session.step(35, "And user clicks on OK button in \"Edit Connection\" dialog", () => clickOn(page, el("OK button in \"Edit Connection\" dialog")));
    await session.step(36, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
    await run.scenario("The primary schema opens the identifiers view", async () => {
      await session.step(39, "When user right-clicks on Databases---Postgres---BDD-Conn-Ident-{time} tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---Postgres---BDD-Conn-Ident-{time} tree node inside browse tree"))));
      await session.step(40, "And user picks \"Configure Identifiers...\" from the open menu", () => pickFromOpenMenu(page, "Configure Identifiers..."));
      await session.step(41, "Then \"Select primary schema for Identifiers Configuration\" dialog should be visible", () => shouldBe(page, el("\"Select primary schema for Identifiers Configuration\" dialog"), "visible"));
      await session.step(42, "And Schema input in \"Select primary schema for Identifiers Configuration\" dialog should offer \"public, pg_catalog, information_schema\"", () => shouldOffer(page, el("Schema input in \"Select primary schema for Identifiers Configuration\" dialog"), "public, pg_catalog, information_schema"));
      await session.step(43, "When user selects \"public\" in Schema input in \"Select primary schema for Identifiers Configuration\" dialog", () => selectIn(page, "public", el("Schema input in \"Select primary schema for Identifiers Configuration\" dialog")));
      await session.step(44, "And user clicks on OK button in \"Select primary schema for Identifiers Configuration\" dialog", () => clickOn(page, el("OK button in \"Select primary schema for Identifiers Configuration\" dialog")));
      await session.step(45, "Then the \"Select primary schema for Identifiers Configuration\" dialog should close", () => dialogCloses(page, "Select primary schema for Identifiers Configuration"));
      await session.step(46, "And \"Add a new identifier\" icon should be visible", () => shouldBe(page, el("\"Add a new identifier\" icon"), "visible"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An identifier is added and saved", async () => {
      await session.step(50, "When user clicks on \"Add a new identifier\" icon", () => clickOn(page, el("\"Add a new identifier\" icon")));
      await session.step(51, "Then \"Add Identifier\" dialog should be visible", () => shouldBe(page, el("\"Add Identifier\" dialog"), "visible"));
      await session.step(52, "When user enters \"CUSTOMER_ID\" into \"Semantic Type\" input in \"Add Identifier\" dialog", () => enterInto(page, "CUSTOMER_ID", el("\"Semantic Type\" input in \"Add Identifier\" dialog")));
      await session.step(53, "And user selects \"customers\" in Table input in \"Add Identifier\" dialog", () => selectIn(page, "customers", el("Table input in \"Add Identifier\" dialog")));
      await session.step(54, "And user selects \"customerid\" in Column input in \"Add Identifier\" dialog", () => selectIn(page, "customerid", el("Column input in \"Add Identifier\" dialog")));
      await session.step(55, "And user enters \"[A-Z]{5}\" into \"Match Regexp\" input in \"Add Identifier\" dialog", () => enterInto(page, "[A-Z]{5}", el("\"Match Regexp\" input in \"Add Identifier\" dialog")));
      await session.step(56, "And user clicks on Add button in \"Add Identifier\" dialog", () => clickOn(page, el("Add button in \"Add Identifier\" dialog")));
      await session.step(57, "Then the \"Add Identifier\" dialog should close", () => dialogCloses(page, "Add Identifier"));
      await session.step(58, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(59, "Then an info balloon should have been shown", () => infoBalloon(page));
      await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The table's column comes with the identifier's semantic type", async () => {
      await session.step(63, "When user reloads the page", () => reloadPage(page));
      await session.step(64, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(65, "And Databases---Postgres---BDD-Conn-Ident-{time} tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("Databases---Postgres---BDD-Conn-Ident-{time} tree node inside browse tree"))));
      await session.step(66, "And Databases---Postgres---BDD-Conn-Ident-{time}---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("Databases---Postgres---BDD-Conn-Ident-{time}---Schemas tree node inside browse tree"))));
      await session.step(67, "And Databases---Postgres---BDD-Conn-Ident-{time}---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("Databases---Postgres---BDD-Conn-Ident-{time}---Schemas---public tree node inside browse tree"))));
      await session.step(68, "When user right-clicks on Databases---Postgres---BDD-Conn-Ident-{time}---Schemas---public---customers tree node inside browse tree", () => rightClickOn(page, el(session.text("Databases---Postgres---BDD-Conn-Ident-{time}---Schemas---public---customers tree node inside browse tree"))));
      await session.step(69, "And user picks \"Get All\" from the open menu", () => pickFromOpenMenu(page, "Get All"));
      await session.step(70, "Then the \"customers\" view should be current", () => viewIsCurrent(page, "customers"));
      await session.step(71, "And the table should have 91 rows", () => rowCount(page, 91));
      await session.step(72, "And \"customerid\" column should have semantic type \"CUSTOMER_ID\"", () => columnSemType(page, "customerid", "CUSTOMER_ID"));
      await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
