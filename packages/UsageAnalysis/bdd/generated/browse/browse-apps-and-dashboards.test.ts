/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-apps-and-dashboards.feature
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
import {clearField, clickOn, collapse, doubleClickOn, expand, hoverOver, isExpanded, rightClickOn, shouldBe, shouldBecomeVisibleWithin, shouldContainText, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, contextPanelOpen, contextPanelShows, openDataset, openProject, packageInstalled, saveAsProject, scriptOnServer, urlShouldContain, viewHoldsViewers, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, closeContextMenu, menuLists, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Apps and Dashboards sections of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-apps-and-dashboards.feature", import.meta.url);
  test("The Apps section lists the installed applications", {tag: ["@browse", "@realizes:views.browse", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "Given Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(34, "Then Tutorials tree node inside browse tree should be visible", () => shouldBe(page, el("Tutorials tree node inside browse tree"), "visible"));
    await session.step(35, "And Misc tree node inside browse tree should be visible", () => shouldBe(page, el("Misc tree node inside browse tree"), "visible"));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
    await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An application opens from the tree and becomes the current object", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(40, "Given Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(41, "When user clicks on Tutorials tree node inside browse tree", () => clickOn(page, el("Tutorials tree node inside browse tree")));
    await session.step(42, "Then the \"Tutorials\" view should be current", () => viewIsCurrent(page, "Tutorials"));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
    await session.step(44, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Model Hub opens from the Compute group", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(48, "Given the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(49, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(50, "And Apps---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute tree node inside browse tree")));
    await session.step(51, "When user clicks on Apps---Compute---Model-Hub tree node inside browse tree", () => clickOn(page, el("Apps---Compute---Model-Hub tree node inside browse tree")));
    await session.step(52, "Then the \"Model Hub\" view should be current", () => viewIsCurrent(page, "Model Hub"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
    await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Dashboards node opens the list of projects", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(57, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(58, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
    await session.step(59, "And the page address should contain \"/projects\"", () => urlShouldContain(page, "/projects"));
    await session.step(62, "And card of gallery should be visible", () => shouldBe(page, el("card of gallery"), "visible"));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
    await session.step(64, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A dashboard saved from a table view opens again with its viewers", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(67, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(68, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(69, "When user saves the current view as project \"BDD-Browse-Dash-{run}\"", () => saveAsProject(page, session.text("BDD-Browse-Dash-{run}")));
    await session.step(70, "And user closes all views", () => closeAllViews(page));
    await session.step(71, "And user opens the \"BDD-Browse-Dash-{run}\" project", () => openProject(page, session.text("BDD-Browse-Dash-{run}")));
    await session.step(72, "Then the current view should hold at least 2 viewers", () => viewHoldsViewers(page, 2));
    await session.step(73, "And no errors should have been logged", () => noErrors(page));
    await session.step(74, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The application list comes back within five seconds of opening Apps", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(77, "Given the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(78, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(79, "And user collapses Apps tree node inside browse tree", () => collapse(page, el("Apps tree node inside browse tree")));
    await session.step(80, "And Apps---Compute tree node inside browse tree should be hidden", () => shouldBe(page, el("Apps---Compute tree node inside browse tree"), "hidden"));
    await session.step(81, "When user expands Apps tree node inside browse tree", () => expand(page, el("Apps tree node inside browse tree")));
    await session.step(82, "Then Apps---Compute tree node inside browse tree should become visible within 5 seconds", () => shouldBecomeVisibleWithin(page, el("Apps---Compute tree node inside browse tree"), 5));
    await session.step(83, "And Apps---Chem tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Chem tree node inside browse tree"), "visible"));
    await session.step(84, "And no errors should have been logged", () => noErrors(page));
  });
  test("An application's tooltip says what it does and which package it comes from", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(87, "Given the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(88, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(89, "And Apps---Chem tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem tree node inside browse tree")));
    await session.step(90, "And Apps---Chem---Reactions tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem---Reactions tree node inside browse tree")));
    await session.step(91, "When user hovers over Apps---Chem---Reactions---Reaction-Enumerator tree node inside browse tree", () => hoverOver(page, el("Apps---Chem---Reactions---Reaction-Enumerator tree node inside browse tree")));
    await session.step(92, "Then tooltip should contain text \"Reaction Enumerator\"", () => shouldContainText(page, el("tooltip"), "Reaction Enumerator"));
    await session.step(93, "And tooltip should contain text \"Forward-reaction library enumeration\"", () => shouldContainText(page, el("tooltip"), "Forward-reaction library enumeration"));
    await session.step(94, "And tooltip should contain text \"Package\"", () => shouldContainText(page, el("tooltip"), "Package"));
    await session.step(95, "And tooltip should contain text \"Chem\"", () => shouldContainText(page, el("tooltip"), "Chem"));
    await session.step(96, "And no errors should have been logged", () => noErrors(page));
  });
  test("A dashboard the Chem package ships opens from the Dashboards gallery with its viewers", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(99, "Given the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(100, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(101, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
    await session.step(102, "When user types \"chemical_space_demo\" into gallery search", () => typeInto(page, "chemical_space_demo", el("gallery search")));
    await session.step(103, "And user double-clicks on \"ChemicalSpaceDemo\" gallery card", () => doubleClickOn(page, el("\"ChemicalSpaceDemo\" gallery card")));
    await session.step(104, "Then the current view should hold at least 2 viewers", () => viewHoldsViewers(page, 2));
    await session.step(105, "And no errors should have been logged", () => noErrors(page));
    await session.step(106, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The context panel follows from one dashboard to the next", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(110, "Given the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(111, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(112, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(113, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
    await session.step(114, "When user types \"chemical_space_demo\" into gallery search", () => typeInto(page, "chemical_space_demo", el("gallery search")));
    await session.step(115, "And user clicks on \"ChemicalSpaceDemo\" gallery card", () => clickOn(page, el("\"ChemicalSpaceDemo\" gallery card")));
    await session.step(116, "Then the context panel should show \"chemical_space_demo\"", () => contextPanelShows(page, "chemical_space_demo"));
    await session.step(117, "When user types \"demo_activity_cliffs\" into gallery search", () => typeInto(page, "demo_activity_cliffs", el("gallery search")));
    await session.step(118, "And user clicks on \"DemoActivityCliffs\" gallery card", () => clickOn(page, el("\"DemoActivityCliffs\" gallery card")));
    await session.step(119, "Then the context panel should show \"demo_activity_cliffs\"", () => contextPanelShows(page, "demo_activity_cliffs"));
    await session.step(120, "And context panel should not contain text \"chemical_space_demo\"", () => shouldNotContainText(page, el("context panel"), "chemical_space_demo"));
    await session.step(121, "When user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(122, "Then no errors should have been logged", () => noErrors(page));
    await session.step(123, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Model Hub catalog lists the model", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(127, "Given the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(128, "And a script \"BddBrowseModel\" is on the server:", () => scriptOnServer(page, "BddBrowseModel", "//language: javascript\n//meta.role: model\n//description: A model a BDD feature saved\n//input: int x = 1\n//output: int result\nresult = x + 1;"));
    await session.step(137, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(138, "And Apps---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute tree node inside browse tree")));
    await session.step(141, "When user clicks on Apps---Compute---Model-Hub tree node inside browse tree", () => clickOn(page, el("Apps---Compute---Model-Hub tree node inside browse tree")));
    await session.step(142, "Then the \"Model Hub\" view should be current", () => viewIsCurrent(page, "Model Hub"));
    await session.step(143, "And \"BddBrowseModel\" link should be visible", () => shouldBe(page, el("\"BddBrowseModel\" link"), "visible"));
    await session.step(144, "And no errors should have been logged", () => noErrors(page));
  });
  test("Uncategorized opens to the model, a hover explains it and a click previews it", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(127, "Given the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(128, "And a script \"BddBrowseModel\" is on the server:", () => scriptOnServer(page, "BddBrowseModel", "//language: javascript\n//meta.role: model\n//description: A model a BDD feature saved\n//input: int x = 1\n//output: int result\nresult = x + 1;"));
    await session.step(137, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(138, "And Apps---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute tree node inside browse tree")));
    await session.step(147, "Given Apps---Compute---Model-Hub tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute---Model-Hub tree node inside browse tree")));
    await session.step(148, "When user expands Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree", () => expand(page, el("Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree")));
    await session.step(149, "Then Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree"), "visible"));
    await session.step(150, "When user hovers over Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree", () => hoverOver(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree")));
    await session.step(151, "Then tooltip should contain text \"A model a BDD feature saved\"", () => shouldContainText(page, el("tooltip"), "A model a BDD feature saved"));
    await session.step(152, "When user clicks on Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree", () => clickOn(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree")));
    await session.step(153, "Then the \"BddBrowseModel preview\" view should be current", () => viewIsCurrent(page, "BddBrowseModel preview"));
    await session.step(154, "And no errors should have been logged", () => noErrors(page));
    await session.step(155, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A double click keeps the model's view and its menu offers Run", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(127, "Given the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
    await session.step(128, "And a script \"BddBrowseModel\" is on the server:", () => scriptOnServer(page, "BddBrowseModel", "//language: javascript\n//meta.role: model\n//description: A model a BDD feature saved\n//input: int x = 1\n//output: int result\nresult = x + 1;"));
    await session.step(137, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(138, "And Apps---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute tree node inside browse tree")));
    await session.step(158, "Given Apps---Compute---Model-Hub tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute---Model-Hub tree node inside browse tree")));
    await session.step(159, "And Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree")));
    await session.step(160, "When user double-clicks on Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree", () => doubleClickOn(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree")));
    await session.step(161, "Then the \"BddBrowseModel preview\" view should be current", () => viewIsCurrent(page, "BddBrowseModel preview"));
    await session.step(162, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(163, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
    await session.step(164, "And \"BddBrowseModel preview\" view should be present", () => shouldBe(page, el("\"BddBrowseModel preview\" view"), "present"));
    await session.step(166, "Given Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree")));
    await session.step(167, "And Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree"), "visible"));
    await session.step(168, "When user right-clicks on Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree", () => rightClickOn(page, el("Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree")));
    await session.step(169, "Then the open menu should list \"Run...\"", () => menuLists(page, "Run..."));
    await session.step(170, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(171, "Then no errors should have been logged", () => noErrors(page));
    await session.step(172, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
