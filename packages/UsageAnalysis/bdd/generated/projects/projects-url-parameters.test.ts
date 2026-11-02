/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-url-parameters.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, views.functions, core#569, core#566, core#630, core#568, GROK-20006, core#1127]
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
import {clickOn, clipboardContains, enterInto, hoverOver, isExpanded, pressKey, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {dblclickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, doubleClickOpensProject, noProjectOnServer, queryTextOnServer, simpleModeOff, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorOrWarningBalloonMatching, infoBalloonText, menuLists, noBalloons, pickFromContextMenu, pickFromOpenMenu, pointerAway, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A dashboard on a query with a parameter: Toolbox > Source, URL parameters and saved values", () => {
  const session = feature(test, "features/projects/projects-url-parameters.feature", import.meta.url);
  test("A dashboard on a query with a parameter: Toolbox > Source, URL parameters and saved values", {tag: ["@journey", "@realizes:views.projects", "@realizes:views.functions", "@realizes:core#569", "@realizes:core#566", "@realizes:core#630", "@realizes:core#568", "@realizes:GROK-20006", "@realizes:core#1127", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 17, page);
    await session.step(46, "Given user is logged in", () => loggedIn(page));
    await session.step(47, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(48, "And simple mode is off", () => simpleModeOff(page));
    await session.step(49, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(50, "And no project named \"BDDUrlParamProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlParamProj{time}")));
    await session.step(51, "And no project named \"BDDUrlParamCopy{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlParamCopy{time}")));
    await session.step(52, "And a query \"BDDUrlParamQuery{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDUrlParamQuery{time}"), "System:Datagrok", "--input: string typeName = \"Project\"\nselect name from entity_types where name = @typeName"));
    await run.scenario("The query runs from the tree with its default value", async () => {
      await session.step(59, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(60, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(61, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(62, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(63, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(64, "When user double-clicks on Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree and the query it runs completes", () => dblclickRunsQuery(page, el(session.text("Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree"))));
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
    });
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
      await session.step(103, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(104, "And user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(105, "Then tooltip should contain text \"/p/\"", () => shouldContainText(page, el("tooltip"), "/p/"));
      await session.step(106, "And tooltip should contain text \".BDDUrlParamProj{time}?\"", () => shouldContainText(page, el("tooltip"), session.text(".BDDUrlParamProj{time}?")));
      await session.step(107, "And tooltip should not contain text \"run=true\"", () => shouldNotContainText(page, el("tooltip"), "run=true"));
    });
    await run.scenario("Right after the save, the Source link already carries the alias", async () => {
      await session.step(111, "Then tooltip should contain text \".BDDUrlParamProj{time}?type=Project\"", () => shouldContainText(page, el("tooltip"), session.text(".BDDUrlParamProj{time}?type=Project")));
    });
    await run.scenario("Right after the save, Source offers the sliders icon", async () => {
      await session.step(115, "Then sliders-h icon in Source pane in toolbox should be visible", () => shouldBe(page, el("sliders-h icon in Source pane in toolbox"), "visible"));
      await session.step(116, "When user hovers over sliders-h icon in Source pane in toolbox", () => hoverOver(page, el("sliders-h icon in Source pane in toolbox")));
      await session.step(117, "Then tooltip should contain text \"Choose which parameters the dashboard link carries\"", () => shouldContainText(page, el("tooltip"), "Choose which parameters the dashboard link carries"));
    });
    await run.scenario("The project's address takes the alias, ignores the parameter's own name, and falls back to the saved value", async () => {
      await session.step(120, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(121, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(122, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(123, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(124, "And user enters \"BDDUrlParamProj{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("gallery search")));
      await session.step(125, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(126, "And user clicks on BDDUrlParamProj{time} gallery card", () => clickOn(page, el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(127, "Then the context panel should show \"BDDUrlParamProj{time}\"", () => contextPanelShows(page, session.text("BDDUrlParamProj{time}")));
      await session.step(128, "When user clicks on \"Links...\" link in Details pane in context panel", () => clickOn(page, el("\"Links...\" link in Details pane in context panel")));
      await session.step(129, "Then \"Links to BDDUrlParamProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Links to BDDUrlParamProj{time}\" dialog")), "visible"));
      await session.step(130, "And \"Grok name\" input in \"Links to BDDUrlParamProj{time}\" dialog should have value \"Admin:BDDUrlParamProj{time}\"", () => shouldHaveValue(page, el(session.text("\"Grok name\" input in \"Links to BDDUrlParamProj{time}\" dialog")), session.text("Admin:BDDUrlParamProj{time}")));
      await session.step(131, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(132, "Then the \"Links to BDDUrlParamProj{time}\" dialog should close", () => dialogCloses(page, session.text("Links to BDDUrlParamProj{time}")));
      await session.step(133, "When user opens the address \"/p/Admin.BDDUrlParamProj{time}?type=Script\"", () => openAddress(page, session.text("/p/Admin.BDDUrlParamProj{time}?type=Script")));
      await session.step(134, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(135, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
      await session.step(136, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(137, "Then \"Type Name\" input in Source pane in toolbox should have value \"Script\"", () => shouldHaveValue(page, el("\"Type Name\" input in Source pane in toolbox"), "Script"));
      await session.step(138, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(139, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(140, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(141, "When user opens the address \"/p/Admin.BDDUrlParamProj{time}\"", () => openAddress(page, session.text("/p/Admin.BDDUrlParamProj{time}")));
      await session.step(142, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(143, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
      await session.step(144, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(145, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(146, "When user opens the address \"/p/Admin.BDDUrlParamProj{time}?typeName=Script\"", () => openAddress(page, session.text("/p/Admin.BDDUrlParamProj{time}?typeName=Script")));
      await session.step(147, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(148, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
    });
    await run.scenario("Reopened, Source offers the sliders icon, and the link follows a new value before REFRESH", async () => {
      await session.step(151, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(152, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(153, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(154, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(155, "And user enters \"BDDUrlParamProj{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("gallery search")));
      await session.step(156, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(157, "And user double-clicks on BDDUrlParamProj{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(158, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(159, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
      await session.step(160, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(161, "Then sliders-h icon in Source pane in toolbox should be visible", () => shouldBe(page, el("sliders-h icon in Source pane in toolbox"), "visible"));
      await session.step(162, "When user enters \"Script\" into \"Type Name\" input in Source pane in toolbox", () => enterInto(page, "Script", el("\"Type Name\" input in Source pane in toolbox")));
      await session.step(163, "And user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(164, "Then tooltip should contain text \"Copy link:\"", () => shouldContainText(page, el("tooltip"), "Copy link:"));
      await session.step(165, "And tooltip should contain text \".BDDUrlParamProj{time}?type=Script\"", () => shouldContainText(page, el("tooltip"), session.text(".BDDUrlParamProj{time}?type=Script")));
      await session.step(166, "When user clicks on REFRESH button in Source pane in toolbox", () => clickOn(page, el("REFRESH button in Source pane in toolbox")));
      await session.step(167, "Then the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
      await session.step(168, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
      await session.step(169, "When user clicks on copy icon in Source pane in toolbox", () => clickOn(page, el("copy icon in Source pane in toolbox")));
      await session.step(170, "Then check icon in Source pane in toolbox should be visible", () => shouldBe(page, el("check icon in Source pane in toolbox"), "visible"));
      await session.step(171, "And the clipboard should contain text \".BDDUrlParamProj{time}?type=Script\"", () => clipboardContains(page, session.text(".BDDUrlParamProj{time}?type=Script")));
    });
    await run.scenario("The new value is saved into the project (GROK-20006)", async () => {
      await session.step(175, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(176, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(177, "When user clicks on \"Creation script\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog", () => clickOn(page, el(session.text("\"Creation script\" button in \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog"))));
      await session.step(178, "Then \"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog should contain text '\"Script\"'", () => shouldContainText(page, el(session.text("\"BDDUrlParamQuery{time}\" project table in \"Save project\" dialog")), "\"Script\""));
      await session.step(179, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(180, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(181, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(182, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(183, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(184, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(185, "And user enters \"BDDUrlParamProj{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("gallery search")));
      await session.step(186, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(187, "And user double-clicks on BDDUrlParamProj{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(188, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(189, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
      await session.step(190, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
    });
    await run.scenario("A copy saved with a third value leaves the original's value alone", async () => {
      await session.step(193, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(194, "When user enters \"Project\" into \"Type Name\" input in Source pane in toolbox", () => enterInto(page, "Project", el("\"Type Name\" input in Source pane in toolbox")));
      await session.step(195, "And user clicks on REFRESH button in Source pane in toolbox", () => clickOn(page, el("REFRESH button in Source pane in toolbox")));
      await session.step(196, "Then the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
      await session.step(197, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(198, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(199, "When user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(200, "And user enters \"BDDUrlParamCopy{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDUrlParamCopy{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(201, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(202, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(203, "And an info balloon containing 'Project \"BDDUrlParamCopy{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDUrlParamCopy{time}\" uploaded")));
      await session.step(204, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(205, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(206, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(207, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(208, "And user enters \"BDDUrlParamCopy{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamCopy{time}"), el("gallery search")));
      await session.step(209, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(210, "And user double-clicks on BDDUrlParamCopy{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDUrlParamCopy{time} gallery card"))));
      await session.step(211, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(212, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
      await session.step(213, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(214, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(215, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(216, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(217, "And user enters \"BDDUrlParamProj{time}\" into gallery search", () => enterInto(page, session.text("BDDUrlParamProj{time}"), el("gallery search")));
      await session.step(218, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(219, "And user double-clicks on BDDUrlParamProj{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDUrlParamProj{time} gallery card"))));
      await session.step(220, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(221, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
    });
    await run.scenario("The sliders icon takes the parameter out of the link, and the save applies it", async () => {
      await session.step(224, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(225, "When user clicks on sliders-h icon in Source pane in toolbox", () => clickOn(page, el("sliders-h icon in Source pane in toolbox")));
      await session.step(226, "Then the open menu should list \"typeName\"", () => menuLists(page, "typeName"));
      await session.step(227, "When user picks \"typeName\" from the open menu", () => pickFromOpenMenu(page, "typeName"));
      await session.step(228, "Then \"Save dashboard to apply changes\" text in toolbox should be visible", () => shouldBe(page, el("\"Save dashboard to apply changes\" text in toolbox"), "visible"));
    });
    await run.scenario("The save applies the change and the hint goes", async () => {
      await session.step(231, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(232, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(233, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(234, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(235, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(236, "And sliders-h icon in Source pane in toolbox should be visible", () => shouldBe(page, el("sliders-h icon in Source pane in toolbox"), "visible"));
      await session.step(237, "And \"Save dashboard to apply changes\" text in toolbox should be hidden", () => shouldBe(page, el("\"Save dashboard to apply changes\" text in toolbox"), "hidden"));
    });
    await run.scenario("Saved with the parameter taken out, the Source link no longer carries it", async () => {
      await session.step(241, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(242, "And user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(243, "Then tooltip should contain text \"/p/\"", () => shouldContainText(page, el("tooltip"), "/p/"));
      await session.step(244, "And tooltip should not contain text \"?type=\"", () => shouldNotContainText(page, el("tooltip"), "?type="));
    });
    await run.scenario("With the parameter taken out, the address no longer takes it", async () => {
      await session.step(247, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(248, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(249, "When user opens the address \"/p/Admin.BDDUrlParamProj{time}?type=Project\"", () => openAddress(page, session.text("/p/Admin.BDDUrlParamProj{time}?type=Project")));
      await session.step(250, "Then the \"BDDUrlParamQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlParamQuery{time}")));
      await session.step(251, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
    });
    await run.scenario("A value of the wrong type in the address is refused with a balloon naming the parameter and the value", async () => {
      await session.step(255, "Given no project named \"BDDUrlEdgeP{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlEdgeP{time}")));
      await session.step(256, "And a query \"BDDUrlEdgeQ{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDUrlEdgeQ{time}"), "System:Datagrok", "--input: int lim = 2\n--input: string typeName = \"Project\"\nselect name from entity_types where name <> @typeName order by name limit @lim"));
      await session.step(262, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(263, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(264, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(265, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(266, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(267, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(268, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(269, "When user double-clicks on Databases---Postgres---Datagrok---BDDUrlEdgeQ{time} tree node inside browse tree and the query it runs completes", () => dblclickRunsQuery(page, el(session.text("Databases---Postgres---Datagrok---BDDUrlEdgeQ{time} tree node inside browse tree"))));
      await session.step(270, "Then the \"BDDUrlEdgeQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlEdgeQ{time}")));
      await session.step(271, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(272, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(273, "When user enters \"BDDUrlEdgeP{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDUrlEdgeP{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(274, "And user clicks on \"URL Parameters\" button in \"BDDUrlEdgeQ{time}\" project table in \"Save project\" dialog", () => clickOn(page, el(session.text("\"URL Parameters\" button in \"BDDUrlEdgeQ{time}\" project table in \"Save project\" dialog"))));
      await session.step(275, "Then first \"URL alias\" text input in \"Save project\" dialog should have value \"lim\"", () => shouldHaveValue(page, el("first \"URL alias\" text input in \"Save project\" dialog"), "lim"));
      await session.step(276, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(277, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(278, "And \"Share BDDUrlEdgeP{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDUrlEdgeP{time}\" dialog")), "visible"));
      await session.step(279, "When user clicks on CANCEL button in \"Share BDDUrlEdgeP{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDUrlEdgeP{time}\" dialog"))));
      await session.step(280, "Then the \"Share BDDUrlEdgeP{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDUrlEdgeP{time}")));
      await session.step(281, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(282, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(283, "When user opens the address \"/p/{login}.BDDUrlEdgeP{time}?lim=abc\"", () => openAddress(page, session.text("/p/{login}.BDDUrlEdgeP{time}?lim=abc")));
      await session.step(284, "Then an error or warning balloon matching \"lim.*abc|abc.*lim\" should have been shown", () => errorOrWarningBalloonMatching(page, "lim.*abc|abc.*lim"));
    });
    run.finish();
  });
});
