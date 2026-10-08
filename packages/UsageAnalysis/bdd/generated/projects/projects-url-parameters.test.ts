/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-url-parameters.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, views.functions, GROK-20930, GROK-21030, GROK-20929]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {appendToEditor, clickOn, clipboardContains, collapse, doubleClickOn, enterInto, hoverOver, isExpanded, pressKey, replaceCode, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, closeCurrentView, contextPanelShows, currentViewType, dialogCloses, noProjectOnServer, noQueryOnServer, projectsOnServer, queriesOnServer, simpleModeOff, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, menuLists, noBalloons, noErrors, pickFromContextMenu, pickFromOpenMenu, pointerAway, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A dashboard on a query with a parameter: Toolbox > Source, URL parameters and saved values", () => {
  const session = feature(test, "features/projects/projects-url-parameters.feature", import.meta.url);
  test("A dashboard on a query with a parameter: Toolbox > Source, URL parameters and saved values", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:views.functions", "@known-failure", "@realizes:GROK-20930", "@realizes:GROK-21030", "@realizes:GROK-20929"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 18, page);
    await session.step(42, "Given user is logged in", () => loggedIn(page));
    await session.step(43, "And simple mode is off", () => simpleModeOff(page));
    await session.step(44, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(45, "And no project named \"BDDUrlParamProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlParamProj{time}")));
    await session.step(46, "And no project named \"BDDUrlParamCopy{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlParamCopy{time}")));
    await session.step(47, "And no query named \"BDDUrlParamQuery{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDUrlParamQuery{time}")));
    await run.scenario("A query with a typeName parameter is made on entity_types", async () => {
      await session.step(50, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(51, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(52, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(53, "And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas tree node inside browse tree")));
      await session.step(54, "And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree")));
      await session.step(55, "When user picks \"New SQL Query...\" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree", () => pickFromContextMenu(page, "New SQL Query...", el("Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree")));
      await session.step(56, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(57, "When user enters \"BDDUrlParamQuery{time}\" into Name input", () => enterInto(page, session.text("BDDUrlParamQuery{time}"), el("Name input")));
      await session.step(58, "And user replaces the code of code editor with '--input: string typeName = \"Project\"'", () => replaceCode(page, el("code editor"), "--input: string typeName = \"Project\""));
      await session.step(59, "And user appends \"select name from entity_types where name = @typeName\" to code editor", () => appendToEditor(page, "select name from entity_types where name = @typeName", el("code editor")));
      await session.step(60, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(61, "Then 1 query named \"BDDUrlParamQuery{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDUrlParamQuery{time}")));
      await session.step(62, "When user closes the current view", () => closeCurrentView(page));
      await session.step(63, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The query runs from the tree with its default value", async () => {
      await session.step(66, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(67, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(68, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(69, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(70, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(71, "When user double-clicks on Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree"))));
      await session.step(72, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(73, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
      await session.step(74, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
    });
    await run.scenario("Before the save, the Source link runs the query", async () => {
      await session.step(77, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(78, "Then Source pane in toolbox should be visible", () => shouldBe(page, el("Source pane in toolbox"), "visible"));
      await session.step(79, "When user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(80, "Then tooltip should contain text \"/func/\"", () => shouldContainText(page, el("tooltip"), "/func/"));
      await session.step(81, "And tooltip should contain text \"BDDUrlParamQuery{time}\"", () => shouldContainText(page, el("tooltip"), session.text("BDDUrlParamQuery{time}")));
      await session.step(82, "And tooltip should contain text \"typeName=%22Project%22\"", () => shouldContainText(page, el("tooltip"), "typeName=%22Project%22"));
      await session.step(83, "And tooltip should contain text \"run=true\"", () => shouldContainText(page, el("tooltip"), "run=true"));
    });
    await run.scenario("The Save dialog exposes typeName under the table's URL Parameters", async () => {
      await session.step(86, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(87, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(88, "And \"Creation script\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el(session.text("\"Creation script\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog")), "visible"));
      await session.step(89, "And Data sync switch in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el(session.text("Data sync switch in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog")), "checked"));
      await session.step(90, "When user enters \"BDDUrlParamProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(91, "And user clicks on \"URL Parameters\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog", () => clickOn(page, el(session.text("\"URL Parameters\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog"))));
      await session.step(92, "Then \"URL alias\" text input in \"Save project\" dialog should have value \"typeName\"", () => shouldHaveValue(page, el("\"URL alias\" text input in \"Save project\" dialog"), "typeName"));
    });
    await run.scenario("The Share link line follows the name typed into the dialog", async () => {
      await session.step(96, "Then \"Save project\" dialog should contain text \".BDDUrlParamProj{time}?typeName=\"", () => shouldContainText(page, el("\"Save project\" dialog"), session.text(".BDDUrlParamProj{time}?typeName=")));
    }, {knownFailure: true});
    await run.scenario("The alias is changed to \"type\" and the dashboard is saved", async () => {
      await session.step(99, "When user enters \"type\" into \"URL alias\" text input in \"Save project\" dialog", () => enterInto(page, "type", el("\"URL alias\" text input in \"Save project\" dialog")));
      await session.step(100, "Then \"Save project\" dialog should contain text \"?type=Project\"", () => shouldContainText(page, el("\"Save project\" dialog"), "?type=Project"));
      await session.step(101, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(102, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(103, "And an info balloon containing 'Project \"BDDUrlParamProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDUrlParamProj{time}\" uploaded")));
      await session.step(104, "And \"Share BDDUrlParamProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDUrlParamProj{time}\" dialog")), "visible"));
      await session.step(105, "When user clicks on CANCEL button in \"Share BDDUrlParamProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDUrlParamProj{time}\" dialog"))));
      await session.step(106, "Then the \"Share BDDUrlParamProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDUrlParamProj{time}")));
    });
    await run.scenario("After the save, the Source link opens the dashboard without run=true", async () => {
      await session.step(109, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(110, "And user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(111, "Then tooltip should contain text \"/p/\"", () => shouldContainText(page, el("tooltip"), "/p/"));
      await session.step(112, "And tooltip should contain text \".BDDUrlParamProj{time}?\"", () => shouldContainText(page, el("tooltip"), session.text(".BDDUrlParamProj{time}?")));
      await session.step(113, "And tooltip should not contain text \"run=true\"", () => shouldNotContainText(page, el("tooltip"), "run=true"));
    });
    await run.scenario("Right after the save, the Source link already carries the alias", async () => {
      await session.step(117, "Then tooltip should contain text \".BDDUrlParamProj{time}?type=Project\"", () => shouldContainText(page, el("tooltip"), session.text(".BDDUrlParamProj{time}?type=Project")));
    }, {knownFailure: true});
    await run.scenario("Right after the save, Source offers the sliders icon", async () => {
      await session.step(121, "Then sliders-h icon in Source pane in toolbox should be visible", () => shouldBe(page, el("sliders-h icon in Source pane in toolbox"), "visible"));
      await session.step(122, "When user hovers over sliders-h icon in Source pane in toolbox", () => hoverOver(page, el("sliders-h icon in Source pane in toolbox")));
      await session.step(123, "Then tooltip should contain text \"Choose which parameters the dashboard link carries\"", () => shouldContainText(page, el("tooltip"), "Choose which parameters the dashboard link carries"));
    }, {knownFailure: true});
    await run.scenario("The project's address takes the alias, ignores the parameter's own name, and falls back to the saved value", async () => {
      await session.step(126, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(127, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(128, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(129, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(130, "And user enters \"BDDUrlParamProj{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("gallery search")));
      await session.step(131, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(132, "And user clicks on BDDUrlParamProj{time} gallery card", () => clickOn(page, el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(133, "Then the context panel should show \"BDDUrlParamProj{time}\"", () => contextPanelShows(page, session.text("BDDUrlParamProj{time}")));
      await session.step(134, "When user clicks on \"Links...\" link in Details pane in context panel", () => clickOn(page, el("\"Links...\" link in Details pane in context panel")));
      await session.step(135, "Then \"Links to BDDUrlParamProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Links to BDDUrlParamProj{time}\" dialog")), "visible"));
      await session.step(136, "And \"Grok name\" input in \"Links to BDDUrlParamProj{time}\" dialog should have value \"Admin:BDDUrlParamProj{time}\"", () => shouldHaveValue(page, el(session.text("\"Grok name\" input in \"Links to BDDUrlParamProj{time}\" dialog")), session.text("Admin:BDDUrlParamProj{time}")));
      await session.step(137, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(138, "Then the \"Links to BDDUrlParamProj{time}\" dialog should close", () => dialogCloses(page, session.text("Links to BDDUrlParamProj{time}")));
      await session.step(139, "When user opens the address \"/p/Admin.BDDUrlParamProj{time}?type=Script\"", () => openAddress(page, session.text("/p/Admin.BDDUrlParamProj{time}?type=Script")));
      await session.step(140, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(141, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
      await session.step(142, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(143, "Then \"Type Name\" input in Source pane in toolbox should have value \"Script\"", () => shouldHaveValue(page, el("\"Type Name\" input in Source pane in toolbox"), "Script"));
      await session.step(144, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(145, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(146, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(147, "When user opens the address \"/p/Admin.BDDUrlParamProj{time}\"", () => openAddress(page, session.text("/p/Admin.BDDUrlParamProj{time}")));
      await session.step(148, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(149, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
      await session.step(150, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(151, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(152, "When user opens the address \"/p/Admin.BDDUrlParamProj{time}?typeName=Script\"", () => openAddress(page, session.text("/p/Admin.BDDUrlParamProj{time}?typeName=Script")));
      await session.step(153, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(154, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
    });
    await run.scenario("Reopened, Source offers the sliders icon, and the link follows a new value before REFRESH", async () => {
      await session.step(157, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(158, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(159, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(160, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(161, "And user enters \"BDDUrlParamProj{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("gallery search")));
      await session.step(162, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(163, "And user double-clicks on BDDUrlParamProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(164, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(165, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
      await session.step(166, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(167, "Then sliders-h icon in Source pane in toolbox should be visible", () => shouldBe(page, el("sliders-h icon in Source pane in toolbox"), "visible"));
      await session.step(168, "When user enters \"Script\" into \"Type Name\" input in Source pane in toolbox", () => enterInto(page, "Script", el("\"Type Name\" input in Source pane in toolbox")));
      await session.step(169, "And user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(170, "Then tooltip should contain text \"Copy link:\"", () => shouldContainText(page, el("tooltip"), "Copy link:"));
      await session.step(171, "And tooltip should contain text \".BDDUrlParamProj{time}?type=Script\"", () => shouldContainText(page, el("tooltip"), session.text(".BDDUrlParamProj{time}?type=Script")));
      await session.step(172, "When user clicks on REFRESH button in Source pane in toolbox", () => clickOn(page, el("REFRESH button in Source pane in toolbox")));
      await session.step(173, "Then the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
      await session.step(174, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
      await session.step(175, "When user clicks on copy icon in Source pane in toolbox", () => clickOn(page, el("copy icon in Source pane in toolbox")));
      await session.step(176, "Then check icon in Source pane in toolbox should be visible", () => shouldBe(page, el("check icon in Source pane in toolbox"), "visible"));
      await session.step(177, "And the clipboard should contain text \".BDDUrlParamProj{time}?type=Script\"", () => clipboardContains(page, session.text(".BDDUrlParamProj{time}?type=Script")));
    });
    await run.scenario("The new value is saved into the project (GROK-20006)", async () => {
      await session.step(180, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(181, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(182, "When user clicks on \"Creation script\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog", () => clickOn(page, el(session.text("\"Creation script\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog"))));
      await session.step(183, "Then \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog should contain text '\"Script\"'", () => shouldContainText(page, el(session.text("\"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog")), "\"Script\""));
      await session.step(184, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(185, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(186, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(187, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(188, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(189, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(190, "And user enters \"BDDUrlParamProj{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("gallery search")));
      await session.step(191, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(192, "And user double-clicks on BDDUrlParamProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(193, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(194, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
      await session.step(195, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
    });
    await run.scenario("A copy saved with a third value leaves the original's value alone", async () => {
      await session.step(198, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(199, "When user enters \"Project\" into \"Type Name\" input in Source pane in toolbox", () => enterInto(page, "Project", el("\"Type Name\" input in Source pane in toolbox")));
      await session.step(200, "And user clicks on REFRESH button in Source pane in toolbox", () => clickOn(page, el("REFRESH button in Source pane in toolbox")));
      await session.step(201, "Then the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
      await session.step(202, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(203, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(204, "When user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(205, "And user enters \"BDDUrlParamCopy{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDUrlParamCopy{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(206, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(207, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(208, "And an info balloon containing 'Project \"BDDUrlParamCopy{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDUrlParamCopy{time}\" uploaded")));
      await session.step(209, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(210, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(211, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(212, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(213, "And user enters \"BDDUrlParamCopy{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamCopy{time}"), el("gallery search")));
      await session.step(214, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(215, "And user double-clicks on BDDUrlParamCopy{time} gallery card", () => doubleClickOn(page, el(session.text("BDDUrlParamCopy{time} gallery card"))));
      await session.step(216, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(217, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
      await session.step(218, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(219, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(220, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(221, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(222, "And user enters \"BDDUrlParamProj{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("gallery search")));
      await session.step(223, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(224, "And user double-clicks on BDDUrlParamProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(225, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(226, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
    });
    await run.scenario("The sliders icon takes the parameter out of the link, and the save applies it", async () => {
      await session.step(229, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(230, "When user clicks on sliders-h icon in Source pane in toolbox", () => clickOn(page, el("sliders-h icon in Source pane in toolbox")));
      await session.step(231, "Then the open menu should list \"typeName\"", () => menuLists(page, "typeName"));
      await session.step(232, "When user picks \"typeName\" from the open menu", () => pickFromOpenMenu(page, "typeName"));
      await session.step(233, "Then \"Save dashboard to apply changes\" text in toolbox should be visible", () => shouldBe(page, el("\"Save dashboard to apply changes\" text in toolbox"), "visible"));
    });
    await run.scenario("The save applies the change and the hint goes", async () => {
      await session.step(236, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(237, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(238, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(239, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(240, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(241, "And sliders-h icon in Source pane in toolbox should be visible", () => shouldBe(page, el("sliders-h icon in Source pane in toolbox"), "visible"));
      await session.step(242, "And \"Save dashboard to apply changes\" text in toolbox should be hidden", () => shouldBe(page, el("\"Save dashboard to apply changes\" text in toolbox"), "hidden"));
    });
    await run.scenario("Saved with the parameter taken out, the Source link no longer carries it", async () => {
      await session.step(246, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(247, "And user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(248, "Then tooltip should contain text \"/p/\"", () => shouldContainText(page, el("tooltip"), "/p/"));
      await session.step(249, "And tooltip should not contain text \"?type=\"", () => shouldNotContainText(page, el("tooltip"), "?type="));
    }, {knownFailure: true});
    await run.scenario("With the parameter taken out, the address no longer takes it", async () => {
      await session.step(252, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(253, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(254, "When user opens the address \"/p/Admin.BDDUrlParamProj{time}?type=Project\"", () => openAddress(page, session.text("/p/Admin.BDDUrlParamProj{time}?type=Project")));
      await session.step(255, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(256, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
    });
    await run.scenario("The query and the projects are deleted", async () => {
      await session.step(259, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(260, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(261, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(262, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(263, "And user enters \"BDDUrlParam\" into gallery search", () => enterInto(page, "BDDUrlParam", el("gallery search")));
      await session.step(264, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(265, "And user picks \"Delete Project\" from the context menu of BDDUrlParamProj{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(266, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(267, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(268, "When user picks \"Delete Project\" from the context menu of BDDUrlParamCopy{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDUrlParamCopy{time} gallery card"))));
      await session.step(269, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(270, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(271, "And 0 projects named \"BDDUrlParamProj{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDUrlParamProj{time}")));
      await session.step(272, "And 0 projects named \"BDDUrlParamCopy{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDUrlParamCopy{time}")));
      await session.step(273, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(274, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(275, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(276, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(277, "When user collapses Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree", () => collapse(page, el(session.text("Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree"))));
      await session.step(278, "And user picks \"Delete\" from the context menu of Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete", el(session.text("Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree"))));
      await session.step(279, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(280, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(281, "And 0 queries named \"BDDUrlParamQuery{time}\" should be on the server", () => queriesOnServer(page, 0, session.text("BDDUrlParamQuery{time}")));
    });
    run.finish();
  });
});
