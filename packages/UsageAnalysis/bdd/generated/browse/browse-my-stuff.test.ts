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
import {favoriteOnServer, inFavorites, notFavorite, notFavoriteOnServer, openAddress, recentOnServer} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, closeAllViews, connectionOnServer, contextPanelShows, currentViewShows, currentViewType, noProjectOnServer, openDataset, openProject, saveAsProject, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hasArea, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The My stuff section of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-my-stuff.feature", import.meta.url);
  test("My stuff gathers the user's own things", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(39, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(42, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["My-stuff---Recent tree node inside browse tree"],["My-stuff---Favorites tree node inside browse tree"],["My-stuff---Shared-with-me tree node inside browse tree"]]), [["My-stuff---Recent tree node inside browse tree"],["My-stuff---Favorites tree node inside browse tree"],["My-stuff---Shared-with-me tree node inside browse tree"]]);
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The section closes on its own twistie", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(39, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(50, "Given My-stuff---Recent tree node inside browse tree should be visible", () => shouldBe(page, el("My-stuff---Recent tree node inside browse tree"), "visible"));
    await session.step(51, "When user collapses My stuff tree node inside browse tree", () => collapse(page, el("My stuff tree node inside browse tree")));
    await session.step(52, "Then My-stuff---Recent tree node inside browse tree should be hidden", () => shouldBe(page, el("My-stuff---Recent tree node inside browse tree"), "hidden"));
    await session.step(53, "And My-stuff---Favorites tree node inside browse tree should be hidden", () => shouldBe(page, el("My-stuff---Favorites tree node inside browse tree"), "hidden"));
    await session.step(54, "And My stuff tree node inside browse tree should be visible", () => shouldBe(page, el("My stuff tree node inside browse tree"), "visible"));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
    await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A project just opened is in Recent", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(39, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(61, "Given no project named \"bdd-browse-recent-{run}\" is on the server", () => noProjectOnServer(page, session.text("bdd-browse-recent-{run}")));
    await session.step(62, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(63, "When user saves the current view as project \"bdd-browse-recent-{run}\"", () => saveAsProject(page, session.text("bdd-browse-recent-{run}")));
    await session.step(64, "And user closes all views", () => closeAllViews(page));
    await session.step(65, "And user opens the \"bdd-browse-recent-{run}\" project", () => openProject(page, session.text("bdd-browse-recent-{run}")));
    await session.step(66, "Then \"bdd-browse-recent-{run}\" should be among the recently used entities on the server", () => recentOnServer(page, session.text("bdd-browse-recent-{run}")));
    await session.step(68, "And user closes all views", () => closeAllViews(page));
    await session.step(69, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(70, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(71, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(72, "When user collapses My-stuff---Recent tree node inside browse tree", () => collapse(page, el("My-stuff---Recent tree node inside browse tree")));
    await session.step(73, "And user expands My-stuff---Recent tree node inside browse tree", () => expand(page, el("My-stuff---Recent tree node inside browse tree")));
    await session.step(74, "Then My-stuff---Recent---bdd-browse-recent-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("My-stuff---Recent---bdd-browse-recent-{run} tree node inside browse tree")), "visible"));
    await session.step(75, "When user clicks on My-stuff---Recent tree node inside browse tree", () => clickOn(page, el("My-stuff---Recent tree node inside browse tree")));
    await session.step(76, "Then the current view should show the text \"bdd-browse-recent-{run}\"", () => currentViewShows(page, session.text("bdd-browse-recent-{run}")));
    await session.step(77, "And calendar viewer should have a \"busiest day\" area", () => hasArea(page, el("calendar viewer"), "busiest day"));
    await session.step(78, "And \"Show all days\" text should be hidden", () => shouldBe(page, el("\"Show all days\" text"), "hidden"));
    await session.step(79, "When user clicks on the \"busiest day\" area of calendar viewer", () => clickArea(page, "busiest day", el("calendar viewer")));
    await session.step(80, "Then \"Show all days\" text should be visible", () => shouldBe(page, el("\"Show all days\" text"), "visible"));
    await session.step(81, "When user clicks on \"Show all days\" text", () => clickOn(page, el("\"Show all days\" text")));
    await session.step(82, "Then \"Show all days\" text should be hidden", () => shouldBe(page, el("\"Show all days\" text"), "hidden"));
    await session.step(83, "When user closes all views", () => closeAllViews(page));
    await session.step(84, "And user opens the address \"/recent\"", () => openAddress(page, "/recent"));
    await session.step(85, "Then the current view should be a RecentActivity view", () => currentViewType(page, "RecentActivity"));
    await session.step(86, "And the current view should show the text \"bdd-browse-recent-{run}\"", () => currentViewShows(page, session.text("bdd-browse-recent-{run}")));
    await session.step(87, "When user clicks on \"bdd-browse-recent-{run}\" link", () => clickOn(page, el(session.text("\"bdd-browse-recent-{run}\" link"))));
    await session.step(88, "Then the context panel should show \"bdd-browse-recent-{run}\"", () => contextPanelShows(page, session.text("bdd-browse-recent-{run}")));
    await session.step(89, "And no errors should have been logged", () => noErrors(page));
    await session.step(90, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An entity added to favorites from its menu is in Favorites, and out when picked again", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(39, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(93, "Given a \"Postgres\" connection named \"BDD-Browse-Fav-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Browse-Fav-{run}")));
    await session.step(94, "And \"BDD-Browse-Fav-{run}\" is not in favorites", () => notFavorite(page, session.text("BDD-Browse-Fav-{run}")));
    await session.step(95, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(96, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(97, "When user picks \"Add To Favorites > Only for me\" from the context menu of Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el(session.text("Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree"))));
    await session.step(98, "Then \"BDD-Browse-Fav-{run}\" should be in favorites on the server", () => favoriteOnServer(page, session.text("BDD-Browse-Fav-{run}")));
    await session.step(99, "When user collapses My-stuff---Favorites tree node inside browse tree", () => collapse(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(100, "And user expands My-stuff---Favorites tree node inside browse tree", () => expand(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(101, "Then My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree")), "visible"));
    await session.step(102, "When user picks \"Add To Favorites > Only for me\" from the context menu of Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el(session.text("Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree"))));
    await session.step(103, "Then \"BDD-Browse-Fav-{run}\" should not be in favorites on the server", () => notFavoriteOnServer(page, session.text("BDD-Browse-Fav-{run}")));
    await session.step(104, "When user collapses My-stuff---Favorites tree node inside browse tree", () => collapse(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(105, "And user expands My-stuff---Favorites tree node inside browse tree", () => expand(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(106, "Then My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree")), "absent"));
    await session.step(107, "And no errors should have been logged", () => noErrors(page));
    await session.step(108, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The home share is added to favorites from its menu", {tag: ["@browse", "@realizes:views.browse", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(39, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(112, "Given \"My files\" is not in favorites", () => notFavorite(page, "My files"));
    await session.step(113, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(114, "When user picks \"Add To Favorites > Only for me\" from the context menu of Files---My-files tree node inside browse tree", () => pickFromContextMenu(page, "Add To Favorites > Only for me", el("Files---My-files tree node inside browse tree")));
    await session.step(115, "Then \"My files\" should be in favorites on the server", () => favoriteOnServer(page, "My files"));
    await session.step(116, "When user collapses My-stuff---Favorites tree node inside browse tree", () => collapse(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(117, "And user expands My-stuff---Favorites tree node inside browse tree", () => expand(page, el("My-stuff---Favorites tree node inside browse tree")));
    await session.step(118, "Then My-stuff---Favorites---My-files tree node inside browse tree should be visible", () => shouldBe(page, el("My-stuff---Favorites---My-files tree node inside browse tree"), "visible"));
    await session.step(119, "And no errors should have been logged", () => noErrors(page));
    await session.step(120, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Favorites star of the sidebar opens the panel with the account's favorites", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(39, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(123, "Given a \"Postgres\" connection named \"BDD-Browse-FavPanel-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Browse-FavPanel-{run}")));
    await session.step(124, "And the connection \"BDD-Browse-FavPanel-{run}\" is in favorites", () => inFavorites(page, "connection", session.text("BDD-Browse-FavPanel-{run}")));
    await session.step(125, "And favorites panel should be absent", () => shouldBe(page, el("favorites panel"), "absent"));
    await session.step(126, "When user clicks on favorites sidebar tab", () => clickOn(page, el("favorites sidebar tab")));
    await session.step(127, "Then favorites panel should be visible", () => shouldBe(page, el("favorites panel"), "visible"));
    await session.step(128, "And \"BDD-Browse-FavPanel-{run}\" tree node in favorites panel should be visible", () => shouldBe(page, el(session.text("\"BDD-Browse-FavPanel-{run}\" tree node in favorites panel")), "visible"));
    await session.step(129, "When user clicks on favorites sidebar tab", () => clickOn(page, el("favorites sidebar tab")));
    await session.step(130, "Then favorites panel should be absent", () => shouldBe(page, el("favorites panel"), "absent"));
    await session.step(131, "And no errors should have been logged", () => noErrors(page));
    await session.step(132, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A node dragged from the tree onto the Favorites panel becomes a favorite", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(39, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(135, "Given a \"Postgres\" connection named \"BDD-Browse-FavDrag-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Browse-FavDrag-{run}")));
    await session.step(136, "And \"BDD-Browse-FavDrag-{run}\" is not in favorites", () => notFavorite(page, session.text("BDD-Browse-FavDrag-{run}")));
    await session.step(137, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(138, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(139, "When user clicks on favorites sidebar tab", () => clickOn(page, el("favorites sidebar tab")));
    await session.step(140, "Then favorites panel should be visible", () => shouldBe(page, el("favorites panel"), "visible"));
    await session.step(141, "And \"BDD-Browse-FavDrag-{run}\" tree node in favorites panel should be absent", () => shouldBe(page, el(session.text("\"BDD-Browse-FavDrag-{run}\" tree node in favorites panel")), "absent"));
    await session.step(142, "When user drags Databases---Postgres---BDD-Browse-FavDrag-{run} tree node inside browse tree to favorites panel", () => dragTo(page, el(session.text("Databases---Postgres---BDD-Browse-FavDrag-{run} tree node inside browse tree")), el("favorites panel")));
    await session.step(143, "Then \"BDD-Browse-FavDrag-{run}\" should be in favorites on the server", () => favoriteOnServer(page, session.text("BDD-Browse-FavDrag-{run}")));
    await session.step(144, "And \"BDD-Browse-FavDrag-{run}\" tree node in favorites panel should be visible", () => shouldBe(page, el(session.text("\"BDD-Browse-FavDrag-{run}\" tree node in favorites panel")), "visible"));
    await session.step(145, "When user clicks on favorites sidebar tab", () => clickOn(page, el("favorites sidebar tab")));
    await session.step(146, "Then favorites panel should be absent", () => shouldBe(page, el("favorites panel"), "absent"));
    await session.step(147, "And no errors should have been logged", () => noErrors(page));
    await session.step(148, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
