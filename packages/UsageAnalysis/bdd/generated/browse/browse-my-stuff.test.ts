/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-my-stuff.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {collapse, expand, followingShouldBe, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {favoriteOnServer, notFavorite, notFavoriteOnServer, recentOnServer} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, closeAllViews, connectionOnServer, noProjectOnServer, openDataset, openProject, saveAsProject, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The My stuff section of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-my-stuff.feature", import.meta.url);
  test("My stuff gathers the user's own things", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(33, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["My-stuff---Recent tree node inside browse tree"],["My-stuff---Favorites tree node inside browse tree"],["My-stuff---Shared-with-me tree node inside browse tree"]]), [["My-stuff---Recent tree node inside browse tree"],["My-stuff---Favorites tree node inside browse tree"],["My-stuff---Shared-with-me tree node inside browse tree"]]);
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The section closes on its own twistie", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(41, "Given My-stuff---Recent tree node inside browse tree should be visible", () => shouldBe(page, el("My-stuff---Recent tree node inside browse tree"), "visible"));
    await session.step(42, "When user collapses My stuff tree node inside browse tree", () => collapse(page, el("My stuff tree node inside browse tree")));
    await session.step(43, "Then My-stuff---Recent tree node inside browse tree should be hidden", () => shouldBe(page, el("My-stuff---Recent tree node inside browse tree"), "hidden"));
    await session.step(44, "And My-stuff---Favorites tree node inside browse tree should be hidden", () => shouldBe(page, el("My-stuff---Favorites tree node inside browse tree"), "hidden"));
    await session.step(45, "And My stuff tree node inside browse tree should be visible", () => shouldBe(page, el("My stuff tree node inside browse tree"), "visible"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A project just opened is in Recent", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(52, "Given no project named \"bdd-browse-recent-{run}\" is on the server", () => noProjectOnServer(page, session.text("bdd-browse-recent-{run}")));
    await session.step(53, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(54, "When user saves the current view as project \"bdd-browse-recent-{run}\"", () => saveAsProject(page, session.text("bdd-browse-recent-{run}")));
    await session.step(55, "And user closes all views", () => closeAllViews(page));
    await session.step(56, "And user opens the \"bdd-browse-recent-{run}\" project", () => openProject(page, session.text("bdd-browse-recent-{run}")));
    await session.step(57, "Then \"bdd-browse-recent-{run}\" should be among the recently used entities on the server", () => recentOnServer(page, session.text("bdd-browse-recent-{run}")));
    await session.step(59, "And user closes all views", () => closeAllViews(page));
    await session.step(60, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(61, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(62, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(63, "When user collapses My-stuff---Recent tree node inside browse tree", () => collapse(page, el("My-stuff---Recent tree node inside browse tree")));
    await session.step(64, "And user expands My-stuff---Recent tree node inside browse tree", () => expand(page, el("My-stuff---Recent tree node inside browse tree")));
    await session.step(65, "Then My-stuff---Recent---bdd-browse-recent-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("My-stuff---Recent---bdd-browse-recent-{run} tree node inside browse tree")), "visible"));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
    await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An entity added to favorites from its menu is in Favorites, and out after Remove", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(70, "Given a \"Postgres\" connection named \"BDD-Browse-Fav-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Browse-Fav-{run}")));
    await session.step(71, "And \"BDD-Browse-Fav-{run}\" is not in favorites", () => notFavorite(page, session.text("BDD-Browse-Fav-{run}")));
    await session.step(72, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(73, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(74, "When user picks \"Add To Favorites\" from the context menu of Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites", el(session.text("Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree"))));
    await session.step(75, "Then \"BDD-Browse-Fav-{run}\" should be in favorites on the server", () => favoriteOnServer(page, session.text("BDD-Browse-Fav-{run}")));
    await session.step(76, "When user collapses My-stuff---Favorites tree node inside browse tree", () => collapse(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(77, "And user expands My-stuff---Favorites tree node inside browse tree", () => expand(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(78, "Then My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree")), "visible"));
    await session.step(79, "When user picks \"Remove from favorites\" from the context menu of Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Remove from favorites", el(session.text("Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree"))));
    await session.step(80, "Then \"BDD-Browse-Fav-{run}\" should not be in favorites on the server", () => notFavoriteOnServer(page, session.text("BDD-Browse-Fav-{run}")));
    await session.step(81, "When user collapses My-stuff---Favorites tree node inside browse tree", () => collapse(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(82, "And user expands My-stuff---Favorites tree node inside browse tree", () => expand(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(83, "Then My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree")), "absent"));
    await session.step(84, "And no errors should have been logged", () => noErrors(page));
    await session.step(85, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The home share is added to favorites from its menu", {tag: ["@browse", "@realizes:views.browse", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(89, "Given \"My files\" is not in favorites", () => notFavorite(page, "My files"));
    await session.step(90, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(91, "When user picks \"Add To Favorites\" from the context menu of Files---My-files tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites", el("Files---My-files tree node inside browse tree")));
    await session.step(92, "Then \"My files\" should be in favorites on the server", () => favoriteOnServer(page, "My files"));
    await session.step(93, "When user collapses My-stuff---Favorites tree node inside browse tree", () => collapse(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(94, "And user expands My-stuff---Favorites tree node inside browse tree", () => expand(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(95, "Then My-stuff---Favorites---My-files tree node inside browse tree should be visible", () => shouldBe(page, el("My-stuff---Favorites---My-files tree node inside browse tree"), "visible"));
    await session.step(96, "And no errors should have been logged", () => noErrors(page));
    await session.step(97, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
