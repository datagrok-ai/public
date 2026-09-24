/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/get-all-top-100.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.queries]
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
import {isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Get All and Get Top 100 on a table of the schema", () => {
  const session = feature(test, "features/queries/get-all-top-100.feature", import.meta.url);
  test("Get All opens the whole table, Get Top 100 its first hundred rows", {tag: ["@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(19, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(20, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(23, "Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---NorthwindTest---Orders tree node inside browse tree"), "visible"));
    await session.step(24, "Given Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree")));
    await session.step(25, "And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(26, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(27, "When user picks \"Get All\" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => pickFromContextMenu(page, "Get All", el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(28, "Then the task bar should have finished \"orders\"", () => taskBarFinished(page, "orders"));
    await session.step(29, "And table \"orders\" should have 830 rows", () => tableRows(page, "orders", 830));
    await session.step(30, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(32, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(33, "When user picks \"Get Top 100\" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => pickFromContextMenu(page, "Get Top 100", el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(34, "Then the task bar should have finished \"orders\"", () => taskBarFinished(page, "orders"));
    await session.step(35, "And table \"orders (2)\" should have 100 rows", () => tableRows(page, "orders (2)", 100));
    await session.step(36, "And table \"orders\" should have 830 rows", () => tableRows(page, "orders", 830));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
