/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/columns-inspect.feature
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
import {clickEveryColumn, everyColumnShown} from '../../bindings/queries.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {collapse, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Every column of a schema in the context panel", () => {
  const session = feature(test, "features/queries/columns-inspect.feature", import.meta.url);
  test("Every column of every table of PostgresDart's public schema is shown on click [provider=PostgresDart, other=Postgres]", {tag: ["@realizes:views.queries", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(19, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(20, "When user collapses Databases---Postgres tree node inside browse tree", () => collapse(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(21, "Given Databases---PostgresDart tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---PostgresDart tree node inside browse tree")));
    await session.step(22, "And Databases---PostgresDart---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---PostgresDart---NorthwindTest tree node inside browse tree")));
    await session.step(25, "Then Databases---PostgresDart---NorthwindTest---Orders tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---PostgresDart---NorthwindTest---Orders tree node inside browse tree"), "visible"));
    await session.step(27, "Given Databases---PostgresDart---NorthwindTest schemas node inside browse tree is expanded", () => isExpanded(page, el("Databases---PostgresDart---NorthwindTest schemas node inside browse tree")));
    await session.step(28, "And Databases---PostgresDart---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---PostgresDart---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(29, "When user clicks every column of every table under Databases---PostgresDart---NorthwindTest---Schemas---public tree node inside browse tree", () => clickEveryColumn(page, el("Databases---PostgresDart---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(30, "Then every clicked column should have been shown in the context panel with \"General, Actions, Inspect, Database meta\"", () => everyColumnShown(page, "General, Actions, Inspect, Database meta"));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
    await session.step(32, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Every column of every table of Postgres's public schema is shown on click [provider=Postgres, other=PostgresDart]", {tag: ["@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(19, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(20, "When user collapses Databases---PostgresDart tree node inside browse tree", () => collapse(page, el("Databases---PostgresDart tree node inside browse tree")));
    await session.step(21, "Given Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(22, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(25, "Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---Postgres---NorthwindTest---Orders tree node inside browse tree"), "visible"));
    await session.step(27, "Given Databases---Postgres---NorthwindTest schemas node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest schemas node inside browse tree")));
    await session.step(28, "And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(29, "When user clicks every column of every table under Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree", () => clickEveryColumn(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(30, "Then every clicked column should have been shown in the context panel with \"General, Actions, Inspect, Database meta\"", () => everyColumnShown(page, "General, Actions, Inspect, Database meta"));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
    await session.step(32, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
