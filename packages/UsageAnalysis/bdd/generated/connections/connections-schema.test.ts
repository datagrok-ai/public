/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-schema.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {followingShouldBe, isExpanded, rightClickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noBalloons, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The schemas of a connection and the schema view", () => {
  const session = feature(test, "features/connections/connections-schema.feature", import.meta.url);
  test("The schemas of a connection and the schema view", {tag: ["@connections", "@journey", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(28, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(29, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(30, "And Databases---Postgres---CHEMBL tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---CHEMBL tree node inside browse tree")));
    await run.scenario("The Schemas group lists the connection's schemas", async () => {
      await session.step(33, "Given \"Postgres-Chembl-Schemas\" tree group inside browse tree is expanded", () => isExpanded(page, el("\"Postgres-Chembl-Schemas\" tree group inside browse tree")));
      await session.step(34, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Databases---Postgres---CHEMBL---Schemas---public tree node inside browse tree"],["Databases---Postgres---CHEMBL---Schemas---information-schema tree node inside browse tree"]]), [["Databases---Postgres---CHEMBL---Schemas---public tree node inside browse tree"],["Databases---Postgres---CHEMBL---Schemas---information-schema tree node inside browse tree"]]);
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A schema's menu offers the schema view and the table actions", async () => {
      await session.step(40, "When user right-clicks on Databases---Postgres---CHEMBL---Schemas---public tree node inside browse tree", () => rightClickOn(page, el("Databases---Postgres---CHEMBL---Schemas---public tree node inside browse tree")));
      await session.step(41, "Then the open menu should list \"Browse\"", () => menuLists(page, "Browse"));
      await session.step(42, "And the open menu should list \"Open as table\"", () => menuLists(page, "Open as table"));
      await session.step(43, "And the open menu should list \"New Table...\"", () => menuLists(page, "New Table..."));
      await session.step(44, "And the open menu should list \"Import Table...\"", () => menuLists(page, "Import Table..."));
      await session.step(45, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(46, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Browse opens the schema view with a box per table", async () => {
      await session.step(49, "When user right-clicks on Databases---Postgres---CHEMBL---Schemas---public tree node inside browse tree", () => rightClickOn(page, el("Databases---Postgres---CHEMBL---Schemas---public tree node inside browse tree")));
      await session.step(50, "And user picks \"Browse\" from the open menu", () => pickFromOpenMenu(page, "Browse"));
      await session.step(51, "Then the \"Schema: public\" view should be current", () => viewIsCurrent(page, "Schema: public"));
      await session.step(52, "And \"activities\" schema table should be visible", () => shouldBe(page, el("\"activities\" schema table"), "visible"));
      await session.step(53, "And \"molecule_dictionary\" schema table should be visible", () => shouldBe(page, el("\"molecule_dictionary\" schema table"), "visible"));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A table box's menu is the table's menu, and Get Top 100 reads the table", async () => {
      await session.step(58, "When user right-clicks on \"activities\" schema table", () => rightClickOn(page, el("\"activities\" schema table")));
      await session.step(59, "Then the open menu should list \"Get All\"", () => menuLists(page, "Get All"));
      await session.step(60, "And the open menu should list \"Get Top 100\"", () => menuLists(page, "Get Top 100"));
      await session.step(61, "And the open menu should list \"New SQL Query...\"", () => menuLists(page, "New SQL Query..."));
      await session.step(62, "And the open menu should list \"New Visual Query...\"", () => menuLists(page, "New Visual Query..."));
      await session.step(63, "When user picks \"Get Top 100\" from the open menu", () => pickFromOpenMenu(page, "Get Top 100"));
      await session.step(64, "Then the \"activities\" view should be current", () => viewIsCurrent(page, "activities"));
      await session.step(65, "And the table should have 100 rows", () => rowCount(page, 100));
      await session.step(66, "And the table should have a column \"activity_id\"", () => hasColumn(page, "activity_id"));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
