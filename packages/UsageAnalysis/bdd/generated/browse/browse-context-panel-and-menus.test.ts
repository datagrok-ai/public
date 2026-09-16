/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-context-panel-and-menus.feature
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
import {clickOn, isExpanded, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuDoesNotList, menuLists, noBalloons, noErrors, openContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The context panel and the context menus of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-context-panel-and-menus.feature", import.meta.url);
  test("F4 hides the context panel and shows it again", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(32, "When user presses F4", () => pressKey(page, "F4"));
    await session.step(33, "Then context panel should be hidden", () => shouldBe(page, el("context panel"), "hidden"));
    await session.step(34, "When user presses F4", () => pressKey(page, "F4"));
    await session.step(35, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
    await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The context panel follows the node that was clicked", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(40, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(41, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(42, "When user clicks on Tutorials tree node inside browse tree", () => clickOn(page, el("Tutorials tree node inside browse tree")));
    await session.step(43, "Then the context panel should show \"Tutorials\"", () => contextPanelShows(page, "Tutorials"));
    await session.step(46, "When Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(47, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(48, "And user clicks on Databases---Postgres---Datagrok tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(49, "Then the context panel should show \"Datagrok\"", () => contextPanelShows(page, "Datagrok"));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
    await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A connection offers the commands of a connection", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(54, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(55, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(56, "When user opens the context menu of Databases---Postgres---Datagrok tree node inside browse tree", () => openContextMenu(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(57, "Then the open menu should list \"Browse\"", () => menuLists(page, "Browse"));
    await session.step(58, "And the open menu should list \"New Query...\"", () => menuLists(page, "New Query..."));
    await session.step(59, "And the open menu should list \"Edit...\"", () => menuLists(page, "Edit..."));
    await session.step(60, "And the open menu should list \"Rename...\"", () => menuLists(page, "Rename..."));
    await session.step(61, "And the open menu should list \"Clone...\"", () => menuLists(page, "Clone..."));
    await session.step(62, "And the open menu should list \"Delete...\"", () => menuLists(page, "Delete..."));
    await session.step(63, "And the open menu should list \"Clear cache\"", () => menuLists(page, "Clear cache"));
    await session.step(65, "And the open menu should list \"Add to favorites\"", () => menuLists(page, "Add to favorites"));
    await session.step(66, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(67, "Then no errors should have been logged", () => noErrors(page));
    await session.step(68, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A file is not offered the commands of a connection", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(71, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(72, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(73, "When user opens the context menu of Files---Demo---demog.csv tree node inside browse tree", () => openContextMenu(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(74, "Then the open menu should list \"Open\"", () => menuLists(page, "Open"));
    await session.step(75, "And the open menu should list \"Download\"", () => menuLists(page, "Download"));
    await session.step(76, "And the open menu should not list \"New Query...\"", () => menuDoesNotList(page, "New Query..."));
    await session.step(77, "And the open menu should not list \"Clear cache\"", () => menuDoesNotList(page, "Clear cache"));
    await session.step(79, "And the open menu should not list \"Add to favorites\"", () => menuDoesNotList(page, "Add to favorites"));
    await session.step(80, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(81, "Then no errors should have been logged", () => noErrors(page));
    await session.step(82, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
