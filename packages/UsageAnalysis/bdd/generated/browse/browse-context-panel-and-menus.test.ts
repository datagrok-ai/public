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
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, isExpanded, pressKey, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuDoesNotList, menuLists, noBalloons, noErrors, openContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The context panel and the context menus of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-context-panel-and-menus.feature", import.meta.url);
  test("F4 hides the context panel and shows it again", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(34, "When user presses F4", () => pressKey(page, "F4"));
    await session.step(35, "Then context panel should be hidden", () => shouldBe(page, el("context panel"), "hidden"));
    await session.step(36, "When user presses F4", () => pressKey(page, "F4"));
    await session.step(37, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
    await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The context panel follows the node that was clicked", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(42, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(43, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(44, "When user clicks on Tutorials tree node inside browse tree", () => clickOn(page, el("Tutorials tree node inside browse tree")));
    await session.step(45, "Then the context panel should show \"Tutorials\"", () => contextPanelShows(page, "Tutorials"));
    await session.step(48, "When Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(49, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(50, "And user clicks on Databases---Postgres---Datagrok tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(51, "Then the context panel should show \"Datagrok\"", () => contextPanelShows(page, "Datagrok"));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
    await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A connection offers Browse, the query commands, Edit, Rename, Clone, Delete and Clear cache", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(56, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(57, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(58, "When user opens the context menu of Databases---Postgres---Datagrok tree node inside browse tree", () => openContextMenu(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(59, "Then the open menu should list \"Browse\"", () => menuLists(page, "Browse"));
    await session.step(60, "And the open menu should list \"New Query...\"", () => menuLists(page, "New Query..."));
    await session.step(61, "And the open menu should list \"Edit...\"", () => menuLists(page, "Edit..."));
    await session.step(62, "And the open menu should list \"Rename...\"", () => menuLists(page, "Rename..."));
    await session.step(63, "And the open menu should list \"Clone...\"", () => menuLists(page, "Clone..."));
    await session.step(64, "And the open menu should list \"Delete...\"", () => menuLists(page, "Delete..."));
    await session.step(65, "And the open menu should list \"Clear cache\"", () => menuLists(page, "Clear cache"));
    await session.step(67, "And the open menu should list \"Add to favorites\"", () => menuLists(page, "Add to favorites"));
    await session.step(68, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(69, "Then no errors should have been logged", () => noErrors(page));
    await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A file is not offered the commands of a connection", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(73, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(74, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(75, "When user opens the context menu of Files---Demo---demog.csv tree node inside browse tree", () => openContextMenu(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(76, "Then the open menu should list \"Open\"", () => menuLists(page, "Open"));
    await session.step(77, "And the open menu should list \"Download\"", () => menuLists(page, "Download"));
    await session.step(78, "And the open menu should not list \"New Query...\"", () => menuDoesNotList(page, "New Query..."));
    await session.step(79, "And the open menu should not list \"Clear cache\"", () => menuDoesNotList(page, "Clear cache"));
    await session.step(81, "And the open menu should not list \"Add to favorites\"", () => menuDoesNotList(page, "Add to favorites"));
    await session.step(82, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(83, "Then no errors should have been logged", () => noErrors(page));
    await session.step(84, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Back and Forward walk the panel's history", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(87, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(88, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(89, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(90, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(91, "Then the context panel should show \"demog.csv\"", () => contextPanelShows(page, "demog.csv"));
    await session.step(92, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(93, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(94, "When user clicks on Databases---Postgres---Datagrok tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(95, "Then the context panel should show \"Datagrok\"", () => contextPanelShows(page, "Datagrok"));
    await session.step(97, "When user hovers over \"Expand all\" icon", () => hoverOver(page, el("\"Expand all\" icon")));
    await session.step(98, "And user clicks on \"Back\" icon", () => clickOn(page, el("\"Back\" icon")));
    await session.step(99, "Then context panel should contain text \"demog.csv\"", () => shouldContainText(page, el("context panel"), "demog.csv"));
    await session.step(100, "And context panel should not contain text \"Datagrok\"", () => shouldNotContainText(page, el("context panel"), "Datagrok"));
    await session.step(101, "When user hovers over \"Expand all\" icon", () => hoverOver(page, el("\"Expand all\" icon")));
    await session.step(102, "And user clicks on \"Forward\" icon", () => clickOn(page, el("\"Forward\" icon")));
    await session.step(103, "Then context panel should contain text \"Datagrok\"", () => shouldContainText(page, el("context panel"), "Datagrok"));
    await session.step(104, "And context panel should not contain text \"demog.csv\"", () => shouldNotContainText(page, el("context panel"), "demog.csv"));
    await session.step(105, "And no errors should have been logged", () => noErrors(page));
    await session.step(106, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Collapse all and Expand all fold every pane of the panel", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(109, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(110, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(111, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(112, "When user clicks on Databases---Postgres---Datagrok tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(113, "Then the context panel should show \"Datagrok\"", () => contextPanelShows(page, "Datagrok"));
    await session.step(114, "When user clicks on \"Expand all\" icon", () => clickOn(page, el("\"Expand all\" icon")));
    await session.step(115, "Then \"Details\" accordion header in context panel should be expanded", () => shouldBe(page, el("\"Details\" accordion header in context panel"), "expanded"));
    await session.step(116, "When user clicks on \"Collapse all\" icon", () => clickOn(page, el("\"Collapse all\" icon")));
    await session.step(117, "Then \"Details\" accordion header in context panel should be collapsed", () => shouldBe(page, el("\"Details\" accordion header in context panel"), "collapsed"));
    await session.step(118, "When user clicks on \"Expand all\" icon", () => clickOn(page, el("\"Expand all\" icon")));
    await session.step(119, "Then \"Details\" accordion header in context panel should be expanded", () => shouldBe(page, el("\"Details\" accordion header in context panel"), "expanded"));
    await session.step(120, "And no errors should have been logged", () => noErrors(page));
    await session.step(121, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
