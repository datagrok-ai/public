/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-my-stuff.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {collapse, followingShouldBe, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The My stuff section of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-my-stuff.feature", import.meta.url);
  test("My stuff gathers the user's own things", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(28, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(31, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["My-stuff---Recent tree node inside browse tree"],["My-stuff---Favorites tree node inside browse tree"],["My-stuff---Shared-with-me tree node inside browse tree"]]));
    await session.step(35, "And no errors should have been logged", () => noErrors(page));
    await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The section closes on its own twistie", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(28, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(39, "Given My-stuff---Recent tree node inside browse tree should be visible", () => shouldBe(page, el("My-stuff---Recent tree node inside browse tree"), "visible"));
    await session.step(40, "When user collapses My stuff tree node inside browse tree", () => collapse(page, el("My stuff tree node inside browse tree")));
    await session.step(41, "Then My-stuff---Recent tree node inside browse tree should be hidden", () => shouldBe(page, el("My-stuff---Recent tree node inside browse tree"), "hidden"));
    await session.step(42, "And My-stuff---Favorites tree node inside browse tree should be hidden", () => shouldBe(page, el("My-stuff---Favorites tree node inside browse tree"), "hidden"));
    await session.step(43, "And My stuff tree node inside browse tree should be visible", () => shouldBe(page, el("My stuff tree node inside browse tree"), "visible"));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
    await session.step(45, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
