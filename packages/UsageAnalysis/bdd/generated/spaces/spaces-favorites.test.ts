/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-favorites.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, noSpaceOnServer, spaceOnServer, spacesOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A space in favorites", () => {
  const session = feature(test, "features/spaces/spaces-favorites.feature", import.meta.url);
  test("A space in favorites", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And no space named \"BDD-Fav\" is on the server", () => noSpaceOnServer(page, "BDD-Fav"));
    await run.scenario("A space is added to favorites", async () => {
      await session.step(22, "Given a space named \"BDD-Fav\" is on the server", () => spaceOnServer(page, "BDD-Fav"));
      await session.step(23, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(24, "And \"My stuff > Favorites > BDD-Fav\" tree node inside browse tree should be absent", () => shouldBe(page, el("\"My stuff > Favorites > BDD-Fav\" tree node inside browse tree"), "absent"));
      await session.step(25, "When user picks \"Add To Favorites > Only for me\" from the context menu of Spaces---BDD-Fav tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el("Spaces---BDD-Fav tree node inside browse tree")));
      await session.step(26, "Then \"My stuff > Favorites > BDD-Fav\" tree node inside browse tree should be present", () => shouldBe(page, el("\"My stuff > Favorites > BDD-Fav\" tree node inside browse tree"), "present"));
    });
    await run.scenario("A space is removed from favorites", async () => {
      await session.step(29, "When user picks \"Add To Favorites > Only for me\" from the context menu of Spaces---BDD-Fav tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el("Spaces---BDD-Fav tree node inside browse tree")));
      await session.step(30, "Then \"My stuff > Favorites > BDD-Fav\" tree node inside browse tree should be absent", () => shouldBe(page, el("\"My stuff > Favorites > BDD-Fav\" tree node inside browse tree"), "absent"));
      await session.step(31, "And \"My stuff > Favorites\" tree node inside browse tree should be present", () => shouldBe(page, el("\"My stuff > Favorites\" tree node inside browse tree"), "present"));
      await session.step(32, "And 1 space named \"BDD-Fav\" should be on the server", () => spacesOnServer(page, 1, "BDD-Fav"));
      await session.step(33, "And Spaces---BDD-Fav tree node inside browse tree should be visible", () => shouldBe(page, el("Spaces---BDD-Fav tree node inside browse tree"), "visible"));
    });
    run.finish();
  });
});
