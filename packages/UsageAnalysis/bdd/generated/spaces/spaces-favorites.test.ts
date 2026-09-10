/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-favorites.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space]
--- */
import {test} from '@playwright/test';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noSpaceOnServer, spacesOnServer} from '../../bindings/spaces.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A space in favorites", () => {
  const session = feature(test, "features/spaces/spaces-favorites.feature", import.meta.url);
  test("A space in favorites", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(14, "And no space named \"BDD-Fav\" is on the server", () => noSpaceOnServer(page, "BDD-Fav"));
    await run.scenario("A space is added to favorites", async () => {
      await session.step(17, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(18, "And user enters \"BDD-Fav\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Fav", el("Name input in Create Space dialog")));
      await session.step(19, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(20, "Then 1 space named \"BDD-Fav\" should be on the server", () => spacesOnServer(page, 1, "BDD-Fav"));
      await session.step(21, "And \"My stuff > Favorites > BDD-Fav\" tree node inside browse tree should be absent", () => shouldBe(page, el("\"My stuff > Favorites > BDD-Fav\" tree node inside browse tree"), "absent"));
      await session.step(22, "When user picks \"Add to favorites\" from the context menu of BDD-Fav tree node inside browse tree", () => pickFromContextMenu(page, "Add to favorites", el("BDD-Fav tree node inside browse tree")));
      await session.step(23, "Then \"My stuff > Favorites > BDD-Fav\" tree node inside browse tree should be present", () => shouldBe(page, el("\"My stuff > Favorites > BDD-Fav\" tree node inside browse tree"), "present"));
    });
    await run.scenario("A space is removed from favorites", async () => {
      await session.step(26, "When user picks \"Remove from favorites\" from the context menu of BDD-Fav tree node inside browse tree", () => pickFromContextMenu(page, "Remove from favorites", el("BDD-Fav tree node inside browse tree")));
      await session.step(27, "Then \"My stuff > Favorites > BDD-Fav\" tree node inside browse tree should be absent", () => shouldBe(page, el("\"My stuff > Favorites > BDD-Fav\" tree node inside browse tree"), "absent"));
      await session.step(28, "And \"My stuff > Favorites\" tree node inside browse tree should be present", () => shouldBe(page, el("\"My stuff > Favorites\" tree node inside browse tree"), "present"));
      await session.step(29, "And 1 space named \"BDD-Fav\" should be on the server", () => spacesOnServer(page, 1, "BDD-Fav"));
      await session.step(30, "And BDD-Fav tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Fav tree node inside browse tree"), "visible"));
    });
    run.finish();
  });
});
