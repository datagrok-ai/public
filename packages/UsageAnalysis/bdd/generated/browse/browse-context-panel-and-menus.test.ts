/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-context-panel-and-menus.feature
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
import {clickOn, hoverOver, isExpanded, pressKey, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {favoriteOnServer, notFavorite, notFavoriteOnServer} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, connectionOnServer, contextPanelOpen, contextPanelShows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuDoesNotList, menuLists, noBalloons, noErrors, openContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The context panel and the context menus of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-context-panel-and-menus.feature", import.meta.url);
  test("F4 hides the context panel and shows it again", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(35, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(36, "When user presses F4", () => pressKey(page, "F4"));
    await session.step(37, "Then context panel should be hidden", () => shouldBe(page, el("context panel"), "hidden"));
    await session.step(38, "When user presses F4", () => pressKey(page, "F4"));
    await session.step(39, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
    await session.step(41, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The context panel follows the node that was clicked", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(44, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(45, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(46, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(47, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(48, "Then the context panel should show \"demog.csv\"", () => contextPanelShows(page, "demog.csv"));
    await session.step(51, "When Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(52, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(53, "And user clicks on Databases---Postgres---Datagrok tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(54, "Then the context panel should show \"Datagrok\"", () => contextPanelShows(page, "Datagrok"));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
    await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A connection offers Browse, the query commands, Edit, Rename, Clone, Delete and Clear cache", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(59, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(60, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(61, "When user opens the context menu of Databases---Postgres---Datagrok tree node inside browse tree", () => openContextMenu(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(62, "Then the open menu should list \"Browse\"", () => menuLists(page, "Browse"));
    await session.step(63, "And the open menu should list \"New Query...\"", () => menuLists(page, "New Query..."));
    await session.step(64, "And the open menu should list \"Edit...\"", () => menuLists(page, "Edit..."));
    await session.step(65, "And the open menu should list \"Rename...\"", () => menuLists(page, "Rename..."));
    await session.step(66, "And the open menu should list \"Clone...\"", () => menuLists(page, "Clone..."));
    await session.step(67, "And the open menu should list \"Delete...\"", () => menuLists(page, "Delete..."));
    await session.step(68, "And the open menu should list \"Clear cache\"", () => menuLists(page, "Clear cache"));
    await session.step(70, "And the open menu should list \"Add to favorites\"", () => menuLists(page, "Add to favorites"));
    await session.step(71, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(72, "Then no errors should have been logged", () => noErrors(page));
    await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A file is not offered the commands of a connection", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(76, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(77, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(78, "When user opens the context menu of Files---Demo---demog.csv tree node inside browse tree", () => openContextMenu(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(79, "Then the open menu should list \"Open\"", () => menuLists(page, "Open"));
    await session.step(80, "And the open menu should list \"Download\"", () => menuLists(page, "Download"));
    await session.step(81, "And the open menu should not list \"New Query...\"", () => menuDoesNotList(page, "New Query..."));
    await session.step(82, "And the open menu should not list \"Clear cache\"", () => menuDoesNotList(page, "Clear cache"));
    await session.step(84, "And the open menu should not list \"Add to favorites\"", () => menuDoesNotList(page, "Add to favorites"));
    await session.step(85, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(86, "Then no errors should have been logged", () => noErrors(page));
    await session.step(87, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Back and Forward walk the panel's history", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(90, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(91, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(92, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(93, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(94, "Then the context panel should show \"demog.csv\"", () => contextPanelShows(page, "demog.csv"));
    await session.step(95, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(96, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(97, "When user clicks on Databases---Postgres---Datagrok tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(98, "Then the context panel should show \"Datagrok\"", () => contextPanelShows(page, "Datagrok"));
    await session.step(100, "When user hovers over \"Expand all\" icon", () => hoverOver(page, el("\"Expand all\" icon")));
    await session.step(101, "And user clicks on \"Back\" icon", () => clickOn(page, el("\"Back\" icon")));
    await session.step(102, "Then context panel should contain text \"demog.csv\"", () => shouldContainText(page, el("context panel"), "demog.csv"));
    await session.step(103, "And context panel should not contain text \"Datagrok\"", () => shouldNotContainText(page, el("context panel"), "Datagrok"));
    await session.step(104, "When user hovers over \"Expand all\" icon", () => hoverOver(page, el("\"Expand all\" icon")));
    await session.step(105, "And user clicks on \"Forward\" icon", () => clickOn(page, el("\"Forward\" icon")));
    await session.step(106, "Then context panel should contain text \"Datagrok\"", () => shouldContainText(page, el("context panel"), "Datagrok"));
    await session.step(107, "And context panel should not contain text \"demog.csv\"", () => shouldNotContainText(page, el("context panel"), "demog.csv"));
    await session.step(108, "And no errors should have been logged", () => noErrors(page));
    await session.step(109, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Collapse all and Expand all fold every pane of the panel", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(113, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(114, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(115, "When user clicks on Databases---Postgres---Datagrok tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(116, "Then the context panel should show \"Datagrok\"", () => contextPanelShows(page, "Datagrok"));
    await session.step(117, "When user clicks on \"Expand all\" icon", () => clickOn(page, el("\"Expand all\" icon")));
    await session.step(118, "Then \"Details\" accordion header in context panel should be expanded", () => shouldBe(page, el("\"Details\" accordion header in context panel"), "expanded"));
    await session.step(119, "When user clicks on \"Collapse all\" icon", () => clickOn(page, el("\"Collapse all\" icon")));
    await session.step(120, "Then \"Details\" accordion header in context panel should be collapsed", () => shouldBe(page, el("\"Details\" accordion header in context panel"), "collapsed"));
    await session.step(121, "When user clicks on \"Expand all\" icon", () => clickOn(page, el("\"Expand all\" icon")));
    await session.step(122, "Then \"Details\" accordion header in context panel should be expanded", () => shouldBe(page, el("\"Details\" accordion header in context panel"), "expanded"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The star beside a connection's name adds it to favorites and takes it out", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(127, "Given a \"Postgres\" connection named \"BDD-Browse-Star-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Browse-Star-{run}")));
    await session.step(128, "And \"BDD-Browse-Star-{run}\" is not in favorites", () => notFavorite(page, session.text("BDD-Browse-Star-{run}")));
    await session.step(129, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(130, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(131, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(132, "When user clicks on Databases---Postgres---BDD-Browse-Star-{run} tree node inside browse tree", () => clickOn(page, el(session.text("Databases---Postgres---BDD-Browse-Star-{run} tree node inside browse tree"))));
    await session.step(133, "Then the context panel should show \"BDD-Browse-Star-{run}\"", () => contextPanelShows(page, session.text("BDD-Browse-Star-{run}")));
    await session.step(134, "When user clicks on favorite star in context panel", () => clickOn(page, el("favorite star in context panel")));
    await session.step(135, "Then \"BDD-Browse-Star-{run}\" should be in favorites on the server", () => favoriteOnServer(page, session.text("BDD-Browse-Star-{run}")));
    await session.step(136, "When user clicks on favorite star in context panel", () => clickOn(page, el("favorite star in context panel")));
    await session.step(137, "Then \"BDD-Browse-Star-{run}\" should not be in favorites on the server", () => notFavoriteOnServer(page, session.text("BDD-Browse-Star-{run}")));
    await session.step(138, "And no errors should have been logged", () => noErrors(page));
    await session.step(139, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A file has no favorite star, a connection has one", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(142, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(143, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(144, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(145, "When user clicks on Databases---Postgres---Datagrok tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(146, "Then the context panel should show \"Datagrok\"", () => contextPanelShows(page, "Datagrok"));
    await session.step(147, "And favorite star in context panel should be visible", () => shouldBe(page, el("favorite star in context panel"), "visible"));
    await session.step(148, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(149, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(150, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(151, "Then the context panel should show \"demog.csv\"", () => contextPanelShows(page, "demog.csv"));
    await session.step(152, "And favorite star in context panel should be absent", () => shouldBe(page, el("favorite star in context panel"), "absent"));
    await session.step(153, "And no errors should have been logged", () => noErrors(page));
    await session.step(154, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
