/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-url-parameters.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, views.functions, GROK-20930, GROK-20929]
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
import {appendToEditor, clickOn, clipboardContains, collapse, doubleClickOn, enterInto, hoverOver, isExpanded, replaceCode, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeCurrentView, currentViewType, dialogCloses, noProjectOnServer, noQueryOnServer, projectsOnServer, queriesOnServer, simpleModeOff, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, menuLists, noErrors, pickFromContextMenu, pickFromOpenMenu, pointerAway, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A dashboard on a query with a parameter: Toolbox > Source, URL parameters and saved values", () => {
  const session = feature(test, "features/projects/projects-url-parameters.feature", import.meta.url);
  test("A dashboard on a query with a parameter: Toolbox > Source, URL parameters and saved values", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:views.functions", "@known-failure", "@realizes:GROK-20930", "@realizes:GROK-20929"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 16, page);
    await session.step(35, "Given user is logged in", () => loggedIn(page));
    await session.step(36, "And simple mode is off", () => simpleModeOff(page));
    await session.step(37, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(38, "And no project named \"BDDUrlParamProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlParamProj{time}")));
    await session.step(39, "And no project named \"BDDUrlParamCopy{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlParamCopy{time}")));
    await session.step(40, "And no query named \"BDDUrlParamQuery{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDUrlParamQuery{time}")));
    await run.scenario("A query with a typeName parameter is made on entity_types", async () => {
      await session.step(43, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(44, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(45, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(46, "And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas tree node inside browse tree")));
      await session.step(47, "And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree")));
      await session.step(48, "When user picks \"New SQL Query...\" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree", () => pickFromContextMenu(page, "New SQL Query...", el("Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree")));
      await session.step(49, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(50, "When user enters \"BDDUrlParamQuery{time}\" into Name input", () => enterInto(page, session.text("BDDUrlParamQuery{time}"), el("Name input")));
      await session.step(51, "And user replaces the code of code editor with '--input: string typeName = \"Project\"'", () => replaceCode(page, el("code editor"), "--input: string typeName = \"Project\""));
      await session.step(52, "And user appends \"select name from entity_types where name = @typeName\" to code editor", () => appendToEditor(page, "select name from entity_types where name = @typeName", el("code editor")));
      await session.step(53, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(54, "Then 1 query named \"BDDUrlParamQuery{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDUrlParamQuery{time}")));
      await session.step(55, "When user closes the current view", () => closeCurrentView(page));
      await session.step(56, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The query runs from the tree with its default value", async () => {
      await session.step(59, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(60, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(61, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(62, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(63, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(64, "When user double-clicks on Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree"))));
      await session.step(65, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(66, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
      await session.step(67, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
    });
    await run.scenario("Before the save, the Source link runs the query", async () => {
      await session.step(70, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(71, "Then Source pane in toolbox should be visible", () => shouldBe(page, el("Source pane in toolbox"), "visible"));
      await session.step(72, "When user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(73, "Then tooltip should contain text \"/func/\"", () => shouldContainText(page, el("tooltip"), "/func/"));
      await session.step(74, "And tooltip should contain text \"BDDUrlParamQuery{time}\"", () => shouldContainText(page, el("tooltip"), session.text("BDDUrlParamQuery{time}")));
      await session.step(75, "And tooltip should contain text \"typeName=%22Project%22\"", () => shouldContainText(page, el("tooltip"), "typeName=%22Project%22"));
      await session.step(76, "And tooltip should contain text \"run=true\"", () => shouldContainText(page, el("tooltip"), "run=true"));
    });
    await run.scenario("The Save dialog exposes typeName under the table's URL Parameters", async () => {
      await session.step(79, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(80, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(81, "And \"Creation script\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el(session.text("\"Creation script\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog")), "visible"));
      await session.step(82, "And Data sync switch in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el(session.text("Data sync switch in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog")), "checked"));
      await session.step(83, "When user enters \"BDDUrlParamProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(84, "And user clicks on \"URL Parameters\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog", () => clickOn(page, el(session.text("\"URL Parameters\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog"))));
      await session.step(85, "Then \"URL alias\" text input in \"Save project\" dialog should have value \"typeName\"", () => shouldHaveValue(page, el("\"URL alias\" text input in \"Save project\" dialog"), "typeName"));
    });
    await run.scenario("The Share link line follows the name typed into the dialog", async () => {
      await session.step(89, "Then \"Save project\" dialog should contain text \".BDDUrlParamProj{time}?typeName=\"", () => shouldContainText(page, el("\"Save project\" dialog"), session.text(".BDDUrlParamProj{time}?typeName=")));
    }, {knownFailure: true});
    await run.scenario("The alias is changed to \"type\" and the dashboard is saved", async () => {
      await session.step(92, "When user enters \"type\" into \"URL alias\" text input in \"Save project\" dialog", () => enterInto(page, "type", el("\"URL alias\" text input in \"Save project\" dialog")));
      await session.step(93, "Then \"Save project\" dialog should contain text \"?type=Project\"", () => shouldContainText(page, el("\"Save project\" dialog"), "?type=Project"));
      await session.step(94, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(95, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(96, "And an info balloon containing 'Project \"BDDUrlParamProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDUrlParamProj{time}\" uploaded")));
      await session.step(97, "And \"Share BDDUrlParamProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDUrlParamProj{time}\" dialog")), "visible"));
      await session.step(98, "When user clicks on CANCEL button in \"Share BDDUrlParamProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDUrlParamProj{time}\" dialog"))));
      await session.step(99, "Then the \"Share BDDUrlParamProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDUrlParamProj{time}")));
    });
    await run.scenario("After the save, the Source link opens the dashboard without run=true", async () => {
      await session.step(102, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(103, "And user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(104, "Then tooltip should contain text \"/p/\"", () => shouldContainText(page, el("tooltip"), "/p/"));
      await session.step(105, "And tooltip should contain text \".BDDUrlParamProj{time}?\"", () => shouldContainText(page, el("tooltip"), session.text(".BDDUrlParamProj{time}?")));
      await session.step(106, "And tooltip should not contain text \"run=true\"", () => shouldNotContainText(page, el("tooltip"), "run=true"));
    });
    await run.scenario("Right after the save, the Source link already carries the alias", async () => {
      await session.step(110, "Then tooltip should contain text \".BDDUrlParamProj{time}?type=Project\"", () => shouldContainText(page, el("tooltip"), session.text(".BDDUrlParamProj{time}?type=Project")));
    }, {knownFailure: true});
    await run.scenario("Right after the save, Source offers the sliders icon", async () => {
      await session.step(114, "Then sliders-h icon in Source pane in toolbox should be visible", () => shouldBe(page, el("sliders-h icon in Source pane in toolbox"), "visible"));
      await session.step(115, "When user hovers over sliders-h icon in Source pane in toolbox", () => hoverOver(page, el("sliders-h icon in Source pane in toolbox")));
      await session.step(116, "Then tooltip should contain text \"Choose which parameters the dashboard link carries\"", () => shouldContainText(page, el("tooltip"), "Choose which parameters the dashboard link carries"));
    }, {knownFailure: true});
    await run.scenario("Reopened, Source offers the sliders icon, and the link follows a new value before REFRESH", async () => {
      await session.step(119, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(120, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(121, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(122, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(123, "And user enters \"BDDUrlParamProj{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("gallery search")));
      await session.step(124, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(125, "And user double-clicks on BDDUrlParamProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(126, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(127, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
      await session.step(128, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(129, "Then sliders-h icon in Source pane in toolbox should be visible", () => shouldBe(page, el("sliders-h icon in Source pane in toolbox"), "visible"));
      await session.step(130, "When user enters \"Script\" into \"Type Name\" input in Source pane in toolbox", () => enterInto(page, "Script", el("\"Type Name\" input in Source pane in toolbox")));
      await session.step(131, "And user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(132, "Then tooltip should contain text \"Copy link:\"", () => shouldContainText(page, el("tooltip"), "Copy link:"));
      await session.step(133, "And tooltip should contain text \".BDDUrlParamProj{time}?type=Script\"", () => shouldContainText(page, el("tooltip"), session.text(".BDDUrlParamProj{time}?type=Script")));
      await session.step(134, "When user clicks on REFRESH button in Source pane in toolbox", () => clickOn(page, el("REFRESH button in Source pane in toolbox")));
      await session.step(135, "Then the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
      await session.step(136, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
      await session.step(137, "When user clicks on copy icon in Source pane in toolbox", () => clickOn(page, el("copy icon in Source pane in toolbox")));
      await session.step(138, "Then check icon in Source pane in toolbox should be visible", () => shouldBe(page, el("check icon in Source pane in toolbox"), "visible"));
      await session.step(139, "And the clipboard should contain text \".BDDUrlParamProj{time}?type=Script\"", () => clipboardContains(page, session.text(".BDDUrlParamProj{time}?type=Script")));
    });
    await run.scenario("The new value is saved into the project (GROK-20006)", async () => {
      await session.step(142, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(143, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(144, "When user clicks on \"Creation script\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog", () => clickOn(page, el(session.text("\"Creation script\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog"))));
      await session.step(145, "Then \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog should contain text '\"Script\"'", () => shouldContainText(page, el(session.text("\"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog")), "\"Script\""));
      await session.step(146, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(147, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(148, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(149, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(150, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(151, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(152, "And user enters \"BDDUrlParamProj{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("gallery search")));
      await session.step(153, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(154, "And user double-clicks on BDDUrlParamProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(155, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(156, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
      await session.step(157, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
    });
    await run.scenario("A copy saved with a third value leaves the original's value alone", async () => {
      await session.step(160, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(161, "When user enters \"Project\" into \"Type Name\" input in Source pane in toolbox", () => enterInto(page, "Project", el("\"Type Name\" input in Source pane in toolbox")));
      await session.step(162, "And user clicks on REFRESH button in Source pane in toolbox", () => clickOn(page, el("REFRESH button in Source pane in toolbox")));
      await session.step(163, "Then the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
      await session.step(164, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(165, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(166, "When user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(167, "And user enters \"BDDUrlParamCopy{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDUrlParamCopy{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(168, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(169, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(170, "And an info balloon containing 'Project \"BDDUrlParamCopy{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDUrlParamCopy{time}\" uploaded")));
      await session.step(171, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(172, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(173, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(174, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(175, "And user enters \"BDDUrlParamCopy{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamCopy{time}"), el("gallery search")));
      await session.step(176, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(177, "And user double-clicks on BDDUrlParamCopy{time} gallery card", () => doubleClickOn(page, el(session.text("BDDUrlParamCopy{time} gallery card"))));
      await session.step(178, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(179, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
      await session.step(180, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(181, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(182, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(183, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(184, "And user enters \"BDDUrlParamProj{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("gallery search")));
      await session.step(185, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(186, "And user double-clicks on BDDUrlParamProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(187, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(188, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
    });
    await run.scenario("The sliders icon takes the parameter out of the link, and the save applies it", async () => {
      await session.step(191, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(192, "When user clicks on sliders-h icon in Source pane in toolbox", () => clickOn(page, el("sliders-h icon in Source pane in toolbox")));
      await session.step(193, "Then the open menu should list \"typeName\"", () => menuLists(page, "typeName"));
      await session.step(194, "When user picks \"typeName\" from the open menu", () => pickFromOpenMenu(page, "typeName"));
      await session.step(195, "Then \"Save dashboard to apply changes\" text in toolbox should be visible", () => shouldBe(page, el("\"Save dashboard to apply changes\" text in toolbox"), "visible"));
    });
    await run.scenario("With the parameter taken out, the Source link no longer carries it", async () => {
      await session.step(199, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(200, "And user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(201, "Then tooltip should contain text \"/p/\"", () => shouldContainText(page, el("tooltip"), "/p/"));
      await session.step(202, "And tooltip should not contain text \"?type=\"", () => shouldNotContainText(page, el("tooltip"), "?type="));
    }, {knownFailure: true});
    await run.scenario("The save applies the change and the hint goes", async () => {
      await session.step(205, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(206, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(207, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(208, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(209, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(210, "And sliders-h icon in Source pane in toolbox should be visible", () => shouldBe(page, el("sliders-h icon in Source pane in toolbox"), "visible"));
      await session.step(211, "And \"Save dashboard to apply changes\" text in toolbox should be hidden", () => shouldBe(page, el("\"Save dashboard to apply changes\" text in toolbox"), "hidden"));
    });
    await run.scenario("The query and the projects are deleted", async () => {
      await session.step(214, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(215, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(216, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(217, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(218, "And user enters \"BDDUrlParam\" into gallery search", () => enterInto(page, "BDDUrlParam", el("gallery search")));
      await session.step(219, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(220, "And user picks \"Delete Project\" from the context menu of BDDUrlParamProj{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(221, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(222, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(223, "When user picks \"Delete Project\" from the context menu of BDDUrlParamCopy{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDUrlParamCopy{time} gallery card"))));
      await session.step(224, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(225, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(226, "And 0 projects named \"BDDUrlParamProj{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDUrlParamProj{time}")));
      await session.step(227, "And 0 projects named \"BDDUrlParamCopy{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDUrlParamCopy{time}")));
      await session.step(228, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(229, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(230, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(231, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(232, "When user collapses Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree", () => collapse(page, el(session.text("Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree"))));
      await session.step(233, "And user picks \"Delete\" from the context menu of Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete", el(session.text("Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree"))));
      await session.step(234, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(235, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(236, "And 0 queries named \"BDDUrlParamQuery{time}\" should be on the server", () => queriesOnServer(page, 0, session.text("BDDUrlParamQuery{time}")));
    });
    run.finish();
  });
});
