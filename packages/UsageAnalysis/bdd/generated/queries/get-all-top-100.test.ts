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
import '../../bindings/queries.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {collapse, hoverOver, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Get All and Get Top 100 on a table of the schema", () => {
  const session = feature(test, "features/queries/get-all-top-100.feature", import.meta.url);
  test("Postgres opens the whole table and its first hundred rows [provider=Postgres, other=PostgresDart]", {tag: ["@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(21, "When user collapses Databases---PostgresDart tree node inside browse tree", () => collapse(page, el("Databases---PostgresDart tree node inside browse tree")));
    await session.step(22, "Given Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(23, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(26, "Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---NorthwindTest---Orders tree node inside browse tree"), "visible"));
    await session.step(28, "Given Databases---Postgres---NorthwindTest schemas node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest schemas node inside browse tree")));
    await session.step(29, "And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(30, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(32, "When user hovers over Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => hoverOver(page, el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(33, "And user picks \"Get All\" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => pickFromContextMenu(page, "Get All", el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(34, "Then the task bar should have finished \"orders\"", () => taskBarFinished(page, "orders"));
    await session.step(35, "And table \"orders\" should have 830 rows", () => tableRows(page, "orders", 830));
    await session.step(36, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(37, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(38, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(40, "When user hovers over Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => hoverOver(page, el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(41, "And user picks \"Get Top 100\" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => pickFromContextMenu(page, "Get Top 100", el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(42, "Then the task bar should have finished \"orders\"", () => taskBarFinished(page, "orders"));
    await session.step(43, "And table \"orders (2)\" should have 100 rows", () => tableRows(page, "orders (2)", 100));
    await session.step(44, "And table \"orders\" should have 830 rows", () => tableRows(page, "orders", 830));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("PostgresDart opens the whole table and its first hundred rows [provider=PostgresDart, other=Postgres]", {tag: ["@realizes:views.queries", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(21, "When user collapses Databases---Postgres tree node inside browse tree", () => collapse(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(22, "Given Databases---PostgresDart tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---PostgresDart tree node inside browse tree")));
    await session.step(23, "And Databases---PostgresDart---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---PostgresDart---NorthwindTest tree node inside browse tree")));
    await session.step(26, "Then Databases---PostgresDart---NorthwindTest---Orders tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---PostgresDart---NorthwindTest---Orders tree node inside browse tree"), "visible"));
    await session.step(28, "Given Databases---PostgresDart---NorthwindTest schemas node inside browse tree is expanded", () => isExpanded(page, el("Databases---PostgresDart---NorthwindTest schemas node inside browse tree")));
    await session.step(29, "And Databases---PostgresDart---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---PostgresDart---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(30, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(32, "When user hovers over Databases---PostgresDart---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => hoverOver(page, el("Databases---PostgresDart---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(33, "And user picks \"Get All\" from the context menu of Databases---PostgresDart---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => pickFromContextMenu(page, "Get All", el("Databases---PostgresDart---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(34, "Then the task bar should have finished \"orders\"", () => taskBarFinished(page, "orders"));
    await session.step(35, "And table \"orders\" should have 830 rows", () => tableRows(page, "orders", 830));
    await session.step(36, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(37, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(38, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(40, "When user hovers over Databases---PostgresDart---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => hoverOver(page, el("Databases---PostgresDart---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(41, "And user picks \"Get Top 100\" from the context menu of Databases---PostgresDart---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => pickFromContextMenu(page, "Get Top 100", el("Databases---PostgresDart---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(42, "Then the task bar should have finished \"orders\"", () => taskBarFinished(page, "orders"));
    await session.step(43, "And table \"orders (2)\" should have 100 rows", () => tableRows(page, "orders (2)", 100));
    await session.step(44, "And table \"orders\" should have 830 rows", () => tableRows(page, "orders", 830));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
