/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-apps-and-dashboards.feature
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
import {clearField, clickOn, collapse, doubleClickOn, expand, hoverOver, isExpanded, rightClickOn, shouldBe, shouldBecomeVisibleWithin, shouldContainText, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, currentViewShows, packageInstalled, scriptOnServer, urlShouldContain, viewHoldsViewers, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Apps and Dashboards sections of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-apps-and-dashboards.feature", import.meta.url);
  test("The Apps section lists the installed applications", {tag: ["@browse", "@realizes:views.browse", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(35, "Given Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(36, "Then Apps---Tutorials tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Tutorials tree node inside browse tree"), "visible"));
    await session.step(37, "And Misc tree node inside browse tree should be visible", () => shouldBe(page, el("Misc tree node inside browse tree"), "visible"));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
    await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An application opens from the tree and becomes the current object", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(42, "Given Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(43, "When user clicks on Apps---Tutorials tree node inside browse tree", () => clickOn(page, el("Apps---Tutorials tree node inside browse tree")));
    await session.step(44, "Then the \"Tutorials\" view should be current", () => viewIsCurrent(page, "Tutorials"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Model Hub opens from the Compute group, and its Tools menu opens Compare Runs", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(50, "Given the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(51, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(52, "And Apps---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute tree node inside browse tree")));
    await session.step(53, "When user clicks on Apps---Compute---Model-Hub tree node inside browse tree", () => clickOn(page, el("Apps---Compute---Model-Hub tree node inside browse tree")));
    await session.step(54, "Then the \"Model Hub\" view should be current", () => viewIsCurrent(page, "Model Hub"));
    await session.step(55, "When user picks \"Tools > Compare Runs...\" from the top menu", () => pickFromTopMenu(page, "Tools > Compare Runs..."));
    await session.step(56, "Then the \"Run Comparison\" view should be current", () => viewIsCurrent(page, "Run Comparison"));
    await session.step(57, "And the current view should show the text \"Add at least two runs (or tables) to compare\"", () => currentViewShows(page, "Add at least two runs (or tables) to compare"));
    await session.step(58, "And no errors should have been logged", () => noErrors(page));
    await session.step(59, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Dashboards node opens the list of projects", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(62, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(63, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
    await session.step(64, "And the page address should contain \"/projects\"", () => urlShouldContain(page, "/projects"));
    await session.step(67, "And card of gallery should be visible", () => shouldBe(page, el("card of gallery"), "visible"));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
    await session.step(69, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The application list comes back within five seconds of opening Apps", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(72, "Given the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(73, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(74, "And user collapses Apps tree node inside browse tree", () => collapse(page, el("Apps tree node inside browse tree")));
    await session.step(75, "And Apps---Compute tree node inside browse tree should be hidden", () => shouldBe(page, el("Apps---Compute tree node inside browse tree"), "hidden"));
    await session.step(76, "When user expands Apps tree node inside browse tree", () => expand(page, el("Apps tree node inside browse tree")));
    await session.step(77, "Then Apps---Compute tree node inside browse tree should become visible within 5 seconds", () => shouldBecomeVisibleWithin(page, el("Apps---Compute tree node inside browse tree"), 5));
    await session.step(78, "And Apps---Chem tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Chem tree node inside browse tree"), "visible"));
    await session.step(79, "And no errors should have been logged", () => noErrors(page));
  });
  test("An application's tooltip says what it does and which package it comes from", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(82, "Given the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(83, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(84, "And Apps---Chem tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem tree node inside browse tree")));
    await session.step(85, "And Apps---Chem---Reactions tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem---Reactions tree node inside browse tree")));
    await session.step(86, "When user hovers over Apps---Chem---Reactions---Reaction-Enumerator tree node inside browse tree", () => hoverOver(page, el("Apps---Chem---Reactions---Reaction-Enumerator tree node inside browse tree")));
    await session.step(87, "Then tooltip should contain text \"Reaction Enumerator\"", () => shouldContainText(page, el("tooltip"), "Reaction Enumerator"));
    await session.step(88, "And tooltip should contain text \"Forward-reaction library enumeration\"", () => shouldContainText(page, el("tooltip"), "Forward-reaction library enumeration"));
    await session.step(89, "And tooltip should contain text \"Package\"", () => shouldContainText(page, el("tooltip"), "Package"));
    await session.step(90, "And tooltip should contain text \"Chem\"", () => shouldContainText(page, el("tooltip"), "Chem"));
    await session.step(91, "And no errors should have been logged", () => noErrors(page));
  });
  test("A dashboard the Chem package ships opens from the Dashboards gallery with its viewers", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(94, "Given the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(95, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(96, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
    await session.step(97, "When user types \"chemical_space_demo\" into gallery search", () => typeInto(page, "chemical_space_demo", el("gallery search")));
    await session.step(98, "And user double-clicks on \"ChemicalSpaceDemo\" gallery card", () => doubleClickOn(page, el("\"ChemicalSpaceDemo\" gallery card")));
    await session.step(99, "Then the current view should hold at least 2 viewers", () => viewHoldsViewers(page, 2));
    await session.step(100, "And no errors should have been logged", () => noErrors(page));
    await session.step(101, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The context panel follows from one dashboard to the next", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(105, "Given the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(106, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(107, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(108, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
    await session.step(109, "When user types \"chemical_space_demo\" into gallery search", () => typeInto(page, "chemical_space_demo", el("gallery search")));
    await session.step(110, "And user clicks on \"ChemicalSpaceDemo\" gallery card", () => clickOn(page, el("\"ChemicalSpaceDemo\" gallery card")));
    await session.step(111, "Then the context panel should show \"chemical_space_demo\"", () => contextPanelShows(page, "chemical_space_demo"));
    await session.step(112, "When user types \"demo_activity_cliffs\" into gallery search", () => typeInto(page, "demo_activity_cliffs", el("gallery search")));
    await session.step(113, "And user clicks on \"DemoActivityCliffs\" gallery card", () => clickOn(page, el("\"DemoActivityCliffs\" gallery card")));
    await session.step(114, "Then the context panel should show \"demo_activity_cliffs\"", () => contextPanelShows(page, "demo_activity_cliffs"));
    await session.step(115, "And context panel should not contain text \"chemical_space_demo\"", () => shouldNotContainText(page, el("context panel"), "chemical_space_demo"));
    await session.step(116, "When user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(117, "Then no errors should have been logged", () => noErrors(page));
    await session.step(118, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Model Hub catalog lists the model", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(122, "Given the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(123, "And a script \"BddBrowseModel\" is on the server:", () => scriptOnServer(page, "BddBrowseModel", "//language: javascript\n//meta.role: model\n//description: A model a BDD feature saved\n//input: int x = 1\n//output: int result\nresult = x + 1;"));
    await session.step(132, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(133, "And Apps---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute tree node inside browse tree")));
    await session.step(136, "When user clicks on Apps---Compute---Model-Hub tree node inside browse tree", () => clickOn(page, el("Apps---Compute---Model-Hub tree node inside browse tree")));
    await session.step(137, "Then the \"Model Hub\" view should be current", () => viewIsCurrent(page, "Model Hub"));
    await session.step(138, "And \"BddBrowseModel\" link should be visible", () => shouldBe(page, el("\"BddBrowseModel\" link"), "visible"));
    await session.step(139, "And no errors should have been logged", () => noErrors(page));
  });
  test("Uncategorized opens to the model, a hover explains it and a click previews it", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(122, "Given the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(123, "And a script \"BddBrowseModel\" is on the server:", () => scriptOnServer(page, "BddBrowseModel", "//language: javascript\n//meta.role: model\n//description: A model a BDD feature saved\n//input: int x = 1\n//output: int result\nresult = x + 1;"));
    await session.step(132, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(133, "And Apps---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute tree node inside browse tree")));
    await session.step(142, "Given Apps---Compute---Model-Hub tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute---Model-Hub tree node inside browse tree")));
    await session.step(143, "When user expands Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree", () => expand(page, el("Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree")));
    await session.step(144, "Then Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree"), "visible"));
    await session.step(145, "When user hovers over Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree", () => hoverOver(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree")));
    await session.step(146, "Then tooltip should contain text \"A model a BDD feature saved\"", () => shouldContainText(page, el("tooltip"), "A model a BDD feature saved"));
    await session.step(147, "When user clicks on Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree", () => clickOn(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree")));
    await session.step(148, "Then the \"BddBrowseModel preview\" view should be current", () => viewIsCurrent(page, "BddBrowseModel preview"));
    await session.step(149, "And no errors should have been logged", () => noErrors(page));
    await session.step(150, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A double click keeps the model's view and its menu offers Run", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(122, "Given the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(123, "And a script \"BddBrowseModel\" is on the server:", () => scriptOnServer(page, "BddBrowseModel", "//language: javascript\n//meta.role: model\n//description: A model a BDD feature saved\n//input: int x = 1\n//output: int result\nresult = x + 1;"));
    await session.step(132, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(133, "And Apps---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute tree node inside browse tree")));
    await session.step(153, "Given Apps---Compute---Model-Hub tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute---Model-Hub tree node inside browse tree")));
    await session.step(154, "And Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree")));
    await session.step(155, "When user double-clicks on Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree", () => doubleClickOn(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree")));
    await session.step(156, "Then the \"BddBrowseModel preview\" view should be current", () => viewIsCurrent(page, "BddBrowseModel preview"));
    await session.step(157, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(158, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
    await session.step(159, "And \"BddBrowseModel preview\" view should be present", () => shouldBe(page, el("\"BddBrowseModel preview\" view"), "present"));
    await session.step(161, "Given Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree")));
    await session.step(162, "And Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree"), "visible"));
    await session.step(163, "When user right-clicks on Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree", () => rightClickOn(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree")));
    await session.step(164, "Then the open menu should list \"Run...\"", () => menuLists(page, "Run..."));
    await session.step(165, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(166, "Then no errors should have been logged", () => noErrors(page));
    await session.step(167, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
