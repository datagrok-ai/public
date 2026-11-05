/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-favorites.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space, core#840]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {collapse, expand, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {favoriteOnServer, notFavorite, notFavoriteOnServer} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, childSpaceOnServer, noSpaceOnServer, spaceOnServer, spacesOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A space in favorites", () => {
  const session = feature(test, "features/spaces/spaces-favorites.feature", import.meta.url);
  test("A space in favorites", {tag: ["@journey", "@spaces", "@realizes:views.space", "@realizes:core#840"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And no space named \"BDD-Fav\" is on the server", () => noSpaceOnServer(page, "BDD-Fav"));
    await session.step(20, "And no space named \"BDD-FavChild\" is on the server", () => noSpaceOnServer(page, "BDD-FavChild"));
    await run.scenario("A space is added to favorites", async () => {
      await session.step(24, "Given a space named \"BDD-Fav\" is on the server", () => spaceOnServer(page, "BDD-Fav"));
      await session.step(25, "And a space named \"BDD-FavChild\" under \"BDD-Fav\" is on the server", () => childSpaceOnServer(page, "BDD-FavChild", "BDD-Fav"));
      await session.step(26, "And \"BDD-FavChild\" is not in favorites", () => notFavorite(page, "BDD-FavChild"));
      await session.step(27, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(28, "And \"My stuff > Favorites > BDD-Fav\" tree node inside browse tree should be absent", () => shouldBe(page, el("\"My stuff > Favorites > BDD-Fav\" tree node inside browse tree"), "absent"));
      await session.step(29, "When user picks \"Add To Favorites > Only for me\" from the context menu of Spaces---BDD-Fav tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el("Spaces---BDD-Fav tree node inside browse tree")));
      await session.step(30, "Then \"My stuff > Favorites > BDD-Fav\" tree node inside browse tree should be present", () => shouldBe(page, el("\"My stuff > Favorites > BDD-Fav\" tree node inside browse tree"), "present"));
      await session.step(31, "Given My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
      await session.step(32, "And My-stuff---Favorites tree node inside browse tree is expanded", () => isExpanded(page, el("My-stuff---Favorites tree node inside browse tree")));
      await session.step(33, "When user expands \"My stuff > Favorites > BDD-Fav\" tree node inside browse tree", () => expand(page, el("\"My stuff > Favorites > BDD-Fav\" tree node inside browse tree")));
      await session.step(34, "Then My-stuff---Favorites---BDD-Fav---BDD-FavChild tree node inside browse tree should be visible", () => shouldBe(page, el("My-stuff---Favorites---BDD-Fav---BDD-FavChild tree node inside browse tree"), "visible"));
      await session.step(35, "Given Spaces---BDD-Fav tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces---BDD-Fav tree node inside browse tree")));
      await session.step(36, "When user picks \"Add To Favorites > Only for me\" from the context menu of Spaces---BDD-Fav---BDD-FavChild tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el("Spaces---BDD-Fav---BDD-FavChild tree node inside browse tree")));
      await session.step(37, "Then \"BDD-FavChild\" should be in favorites on the server", () => favoriteOnServer(page, "BDD-FavChild"));
      await session.step(38, "When user collapses \"My stuff > Favorites\" tree node inside browse tree", () => collapse(page, el("\"My stuff > Favorites\" tree node inside browse tree")));
      await session.step(39, "And user expands \"My stuff > Favorites\" tree node inside browse tree", () => expand(page, el("\"My stuff > Favorites\" tree node inside browse tree")));
      await session.step(40, "Then \"My stuff > Favorites > BDD-Fav\" tree node inside browse tree should be present", () => shouldBe(page, el("\"My stuff > Favorites > BDD-Fav\" tree node inside browse tree"), "present"));
      await session.step(41, "And My-stuff---Favorites---BDD-Fav---BDD-FavChild tree node inside browse tree should be visible", () => shouldBe(page, el("My-stuff---Favorites---BDD-Fav---BDD-FavChild tree node inside browse tree"), "visible"));
      await session.step(42, "And My-stuff---Favorites---BDD-FavChild tree node inside browse tree should be visible", () => shouldBe(page, el("My-stuff---Favorites---BDD-FavChild tree node inside browse tree"), "visible"));
      await session.step(43, "When user picks \"Add To Favorites > Only for me\" from the context menu of Spaces---BDD-Fav---BDD-FavChild tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el("Spaces---BDD-Fav---BDD-FavChild tree node inside browse tree")));
      await session.step(44, "Then \"BDD-FavChild\" should not be in favorites on the server", () => notFavoriteOnServer(page, "BDD-FavChild"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A space is removed from favorites", async () => {
      await session.step(48, "When user picks \"Add To Favorites > Only for me\" from the context menu of Spaces---BDD-Fav tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el("Spaces---BDD-Fav tree node inside browse tree")));
      await session.step(49, "Then \"My stuff > Favorites > BDD-Fav\" tree node inside browse tree should be absent", () => shouldBe(page, el("\"My stuff > Favorites > BDD-Fav\" tree node inside browse tree"), "absent"));
      await session.step(50, "And \"My stuff > Favorites\" tree node inside browse tree should be present", () => shouldBe(page, el("\"My stuff > Favorites\" tree node inside browse tree"), "present"));
      await session.step(51, "And 1 space named \"BDD-Fav\" should be on the server", () => spacesOnServer(page, 1, "BDD-Fav"));
      await session.step(52, "And Spaces---BDD-Fav tree node inside browse tree should be visible", () => shouldBe(page, el("Spaces---BDD-Fav tree node inside browse tree"), "visible"));
    });
    run.finish();
  });
});
