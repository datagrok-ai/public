/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-my-stuff.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
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
import {clickOn, collapse, dragTo, expand, followingShouldBe, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {favoriteOnServer, inFavorites, notFavorite, notFavoriteOnServer, recentOnServer} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, closeAllViews, connectionOnServer, noProjectOnServer, openDataset, openProject, saveAsProject, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The My stuff section of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-my-stuff.feature", import.meta.url);
  test("My stuff gathers the user's own things", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(39, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["My-stuff---Recent tree node inside browse tree"],["My-stuff---Favorites tree node inside browse tree"],["My-stuff---Shared-with-me tree node inside browse tree"]]), [["My-stuff---Recent tree node inside browse tree"],["My-stuff---Favorites tree node inside browse tree"],["My-stuff---Shared-with-me tree node inside browse tree"]]);
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
    await session.step(44, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The section closes on its own twistie", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(47, "Given My-stuff---Recent tree node inside browse tree should be visible", () => shouldBe(page, el("My-stuff---Recent tree node inside browse tree"), "visible"));
    await session.step(48, "When user collapses My stuff tree node inside browse tree", () => collapse(page, el("My stuff tree node inside browse tree")));
    await session.step(49, "Then My-stuff---Recent tree node inside browse tree should be hidden", () => shouldBe(page, el("My-stuff---Recent tree node inside browse tree"), "hidden"));
    await session.step(50, "And My-stuff---Favorites tree node inside browse tree should be hidden", () => shouldBe(page, el("My-stuff---Favorites tree node inside browse tree"), "hidden"));
    await session.step(51, "And My stuff tree node inside browse tree should be visible", () => shouldBe(page, el("My stuff tree node inside browse tree"), "visible"));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
    await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A project just opened is in Recent", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(58, "Given no project named \"bdd-browse-recent-{run}\" is on the server", () => noProjectOnServer(page, session.text("bdd-browse-recent-{run}")));
    await session.step(59, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(60, "When user saves the current view as project \"bdd-browse-recent-{run}\"", () => saveAsProject(page, session.text("bdd-browse-recent-{run}")));
    await session.step(61, "And user closes all views", () => closeAllViews(page));
    await session.step(62, "And user opens the \"bdd-browse-recent-{run}\" project", () => openProject(page, session.text("bdd-browse-recent-{run}")));
    await session.step(63, "Then \"bdd-browse-recent-{run}\" should be among the recently used entities on the server", () => recentOnServer(page, session.text("bdd-browse-recent-{run}")));
    await session.step(65, "And user closes all views", () => closeAllViews(page));
    await session.step(66, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(67, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(68, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(69, "When user collapses My-stuff---Recent tree node inside browse tree", () => collapse(page, el("My-stuff---Recent tree node inside browse tree")));
    await session.step(70, "And user expands My-stuff---Recent tree node inside browse tree", () => expand(page, el("My-stuff---Recent tree node inside browse tree")));
    await session.step(71, "Then My-stuff---Recent---bdd-browse-recent-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("My-stuff---Recent---bdd-browse-recent-{run} tree node inside browse tree")), "visible"));
    await session.step(72, "And no errors should have been logged", () => noErrors(page));
    await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An entity added to favorites from its menu is in Favorites, and out when picked again", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(76, "Given a \"Postgres\" connection named \"BDD-Browse-Fav-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Browse-Fav-{run}")));
    await session.step(77, "And \"BDD-Browse-Fav-{run}\" is not in favorites", () => notFavorite(page, session.text("BDD-Browse-Fav-{run}")));
    await session.step(78, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(79, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(80, "When user picks \"Add To Favorites > Only for me\" from the context menu of Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el(session.text("Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree"))));
    await session.step(81, "Then \"BDD-Browse-Fav-{run}\" should be in favorites on the server", () => favoriteOnServer(page, session.text("BDD-Browse-Fav-{run}")));
    await session.step(82, "When user collapses My-stuff---Favorites tree node inside browse tree", () => collapse(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(83, "And user expands My-stuff---Favorites tree node inside browse tree", () => expand(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(84, "Then My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree")), "visible"));
    await session.step(85, "When user picks \"Add To Favorites > Only for me\" from the context menu of Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el(session.text("Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree"))));
    await session.step(86, "Then \"BDD-Browse-Fav-{run}\" should not be in favorites on the server", () => notFavoriteOnServer(page, session.text("BDD-Browse-Fav-{run}")));
    await session.step(87, "When user collapses My-stuff---Favorites tree node inside browse tree", () => collapse(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(88, "And user expands My-stuff---Favorites tree node inside browse tree", () => expand(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(89, "Then My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree")), "absent"));
    await session.step(90, "And no errors should have been logged", () => noErrors(page));
    await session.step(91, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The home share is added to favorites from its menu", {tag: ["@browse", "@realizes:views.browse", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(95, "Given \"My files\" is not in favorites", () => notFavorite(page, "My files"));
    await session.step(96, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(97, "When user picks \"Add To Favorites > Only for me\" from the context menu of Files---My-files tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el("Files---My-files tree node inside browse tree")));
    await session.step(98, "Then \"My files\" should be in favorites on the server", () => favoriteOnServer(page, "My files"));
    await session.step(99, "When user collapses My-stuff---Favorites tree node inside browse tree", () => collapse(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(100, "And user expands My-stuff---Favorites tree node inside browse tree", () => expand(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(101, "Then My-stuff---Favorites---My-files tree node inside browse tree should be visible", () => shouldBe(page, el("My-stuff---Favorites---My-files tree node inside browse tree"), "visible"));
    await session.step(102, "And no errors should have been logged", () => noErrors(page));
    await session.step(103, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Favorites star of the sidebar opens the panel with the account's favorites", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(106, "Given a \"Postgres\" connection named \"BDD-Browse-FavPanel-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Browse-FavPanel-{run}")));
    await session.step(107, "And the connection \"BDD-Browse-FavPanel-{run}\" is in favorites", () => inFavorites(page, "connection", session.text("BDD-Browse-FavPanel-{run}")));
    await session.step(108, "And favorites panel should be absent", () => shouldBe(page, el("favorites panel"), "absent"));
    await session.step(109, "When user clicks on favorites sidebar tab", () => clickOn(page, el("favorites sidebar tab")));
    await session.step(110, "Then favorites panel should be visible", () => shouldBe(page, el("favorites panel"), "visible"));
    await session.step(111, "And \"BDD-Browse-FavPanel-{run}\" tree node in favorites panel should be visible", () => shouldBe(page, el(session.text("\"BDD-Browse-FavPanel-{run}\" tree node in favorites panel")), "visible"));
    await session.step(112, "When user clicks on favorites sidebar tab", () => clickOn(page, el("favorites sidebar tab")));
    await session.step(113, "Then favorites panel should be absent", () => shouldBe(page, el("favorites panel"), "absent"));
    await session.step(114, "And no errors should have been logged", () => noErrors(page));
    await session.step(115, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A node dragged from the tree onto the Favorites panel becomes a favorite", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(118, "Given a \"Postgres\" connection named \"BDD-Browse-FavDrag-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Browse-FavDrag-{run}")));
    await session.step(119, "And \"BDD-Browse-FavDrag-{run}\" is not in favorites", () => notFavorite(page, session.text("BDD-Browse-FavDrag-{run}")));
    await session.step(120, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(121, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(122, "When user clicks on favorites sidebar tab", () => clickOn(page, el("favorites sidebar tab")));
    await session.step(123, "Then favorites panel should be visible", () => shouldBe(page, el("favorites panel"), "visible"));
    await session.step(124, "And \"BDD-Browse-FavDrag-{run}\" tree node in favorites panel should be absent", () => shouldBe(page, el(session.text("\"BDD-Browse-FavDrag-{run}\" tree node in favorites panel")), "absent"));
    await session.step(125, "When user drags Databases---Postgres---BDD-Browse-FavDrag-{run} tree node inside browse tree to favorites panel", () => dragTo(page, el(session.text("Databases---Postgres---BDD-Browse-FavDrag-{run} tree node inside browse tree")), el("favorites panel")));
    await session.step(126, "Then \"BDD-Browse-FavDrag-{run}\" should be in favorites on the server", () => favoriteOnServer(page, session.text("BDD-Browse-FavDrag-{run}")));
    await session.step(127, "And \"BDD-Browse-FavDrag-{run}\" tree node in favorites panel should be visible", () => shouldBe(page, el(session.text("\"BDD-Browse-FavDrag-{run}\" tree node in favorites panel")), "visible"));
    await session.step(128, "When user clicks on favorites sidebar tab", () => clickOn(page, el("favorites sidebar tab")));
    await session.step(129, "Then favorites panel should be absent", () => shouldBe(page, el("favorites panel"), "absent"));
    await session.step(130, "And no errors should have been logged", () => noErrors(page));
    await session.step(131, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
