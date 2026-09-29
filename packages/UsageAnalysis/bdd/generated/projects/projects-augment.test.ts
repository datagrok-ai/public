/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-augment.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, file.menu.save.tables-as-project]
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
import {check, clickOn, collapse, doubleClickOn, dragTo, enterInto, hoverOver, isExpanded, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, contextPanelShows, dialogCloses, noProjectOnServer, projectsOnServer, saveDialogDataSync, savedWithDataSync, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu, pointerAway} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project saved with a table left out, then augmented by dropping the table onto it", () => {
  const session = feature(test, "features/projects/projects-augment.feature", import.meta.url);
  test("A project saved with a table left out, then augmented by dropping the table onto it", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:file.menu.save.tables-as-project", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And simple mode is off", () => simpleModeOff(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(28, "And no project named \"BDDAugment{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDAugment{time}")));
    await run.scenario("The cross icon of the Save dialog leaves iris out", async () => {
      await session.step(31, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(32, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(33, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(34, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(35, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
      await session.step(36, "And user double-clicks Files---Demo---iris.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---iris.csv tree node inside browse tree")));
      await session.step(37, "Then the \"iris\" view should be current", () => viewIsCurrent(page, "iris"));
      await session.step(38, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(39, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(40, "And \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(41, "And \"iris\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"iris\" project table in \"Save project\" dialog"), "visible"));
      await session.step(42, "And times icon in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("times icon in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(43, "And times icon in \"iris\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("times icon in \"iris\" project table in \"Save project\" dialog"), "visible"));
      await session.step(44, "When user hovers over times icon in \"iris\" project table in \"Save project\" dialog", () => hoverOver(page, el("times icon in \"iris\" project table in \"Save project\" dialog")));
      await session.step(45, "Then tooltip should contain text \"Exclude table from the project.\"", () => shouldContainText(page, el("tooltip"), "Exclude table from the project."));
      await session.step(46, "When user clicks on times icon in \"iris\" project table in \"Save project\" dialog", () => clickOn(page, el("times icon in \"iris\" project table in \"Save project\" dialog")));
      await session.step(47, "Then plus icon in \"iris\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("plus icon in \"iris\" project table in \"Save project\" dialog"), "visible"));
      await session.step(48, "And times icon in \"iris\" project table in \"Save project\" dialog should be hidden", () => shouldBe(page, el("times icon in \"iris\" project table in \"Save project\" dialog"), "hidden"));
    });
    await run.scenario("The plus icon that replaces the cross has its tooltip", async () => {
      await session.step(52, "When user moves the pointer away from \"iris\" project table in \"Save project\" dialog", () => pointerAway(page, el("\"iris\" project table in \"Save project\" dialog")));
      await session.step(53, "And user hovers over plus icon in \"iris\" project table in \"Save project\" dialog", () => hoverOver(page, el("plus icon in \"iris\" project table in \"Save project\" dialog")));
      await session.step(54, "Then tooltip should contain text \"Include table to the project.\"", () => shouldContainText(page, el("tooltip"), "Include table to the project."));
    }, {knownFailure: true});
    await run.scenario("The plus icon brings iris back, the cross leaves it out again, and the project is saved", async () => {
      await session.step(57, "When user clicks on plus icon in \"iris\" project table in \"Save project\" dialog", () => clickOn(page, el("plus icon in \"iris\" project table in \"Save project\" dialog")));
      await session.step(58, "Then times icon in \"iris\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("times icon in \"iris\" project table in \"Save project\" dialog"), "visible"));
      await session.step(59, "And plus icon in \"iris\" project table in \"Save project\" dialog should be hidden", () => shouldBe(page, el("plus icon in \"iris\" project table in \"Save project\" dialog"), "hidden"));
      await session.step(60, "When user clicks on times icon in \"iris\" project table in \"Save project\" dialog", () => clickOn(page, el("times icon in \"iris\" project table in \"Save project\" dialog")));
      await session.step(61, "Then plus icon in \"iris\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("plus icon in \"iris\" project table in \"Save project\" dialog"), "visible"));
      await session.step(62, "When user enters \"BDDAugment{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDAugment{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(63, "Then \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(64, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(65, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(66, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(67, "And an info balloon containing 'Project \"BDDAugment{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDAugment{time}\" uploaded")));
      await session.step(68, "And \"Share BDDAugment{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDAugment{time}\" dialog")), "visible"));
      await session.step(69, "When user clicks on CANCEL button in \"Share BDDAugment{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDAugment{time}\" dialog"))));
      await session.step(70, "Then the \"Share BDDAugment{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDAugment{time}")));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The project holds demog only", async () => {
      await session.step(74, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(75, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(76, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(77, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(78, "And user enters \"BDDAugment{time}\" into gallery search", () => enterInto(page, session.text("BDDAugment{time}"), el("gallery search")));
      await session.step(79, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(80, "And user clicks on BDDAugment{time} gallery card", () => clickOn(page, el(session.text("BDDAugment{time} gallery card"))));
      await session.step(81, "Then the context panel should show \"BDDAugment{time}\"", () => contextPanelShows(page, session.text("BDDAugment{time}")));
      await session.step(82, "Given Content pane in context panel is expanded", () => isExpanded(page, el("Content pane in context panel")));
      await session.step(83, "Then Content pane in context panel should contain text \"demog\"", () => shouldContainText(page, el("Content pane in context panel"), "demog"));
      await session.step(84, "And Content pane in context panel should not contain text \"iris\"", () => shouldNotContainText(page, el("Content pane in context panel"), "iris"));
      await session.step(85, "Given user watches the task bar", () => watchTaskBar(page));
      await session.step(86, "When user double-clicks on BDDAugment{time} gallery card", () => doubleClickOn(page, el(session.text("BDDAugment{time} gallery card"))));
      await session.step(87, "Then the task bar should have finished \"Opening project\"", () => taskBarFinished(page, "Opening project"));
      await session.step(88, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(89, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(90, "And iris view should be absent", () => shouldBe(page, el("iris view"), "absent"));
    });
    await run.scenario("iris dragged onto the project's node in the Dashboards panel moves into the project", async () => {
      await session.step(93, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(94, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(95, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(96, "When user double-clicks Files---Demo---iris.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---iris.csv tree node inside browse tree")));
      await session.step(97, "Then the \"iris\" view should be current", () => viewIsCurrent(page, "iris"));
      await session.step(98, "When user clicks on Dashboards tab", () => clickOn(page, el("Dashboards tab")));
      await session.step(99, "Then \"New Dashboard > iris\" tree node inside browse tree should be visible", () => shouldBe(page, el("\"New Dashboard > iris\" tree node inside browse tree"), "visible"));
      await session.step(100, "Given BDDAugment{time} tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("BDDAugment{time} tree node inside browse tree"))));
      await session.step(101, "Then \"BDDAugment{time} > demog\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDDAugment{time} > demog\" tree node inside browse tree")), "visible"));
      await session.step(102, "When user collapses BDDAugment{time} tree node inside browse tree", () => collapse(page, el(session.text("BDDAugment{time} tree node inside browse tree"))));
      await session.step(103, "And user drags \"New Dashboard > iris\" tree node inside browse tree to BDDAugment{time} tree node inside browse tree", () => dragTo(page, el("\"New Dashboard > iris\" tree node inside browse tree"), el(session.text("BDDAugment{time} tree node inside browse tree"))));
      await session.step(104, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(105, "And Move entity dialog should contain text \"BDDAugment{time} project\"", () => shouldContainText(page, el("Move entity dialog"), session.text("BDDAugment{time} project")));
      await session.step(106, "And Move entity dialog should contain text \"iris\"", () => shouldContainText(page, el("Move entity dialog"), "iris"));
      await session.step(107, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(108, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(109, "Given BDDAugment{time} tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("BDDAugment{time} tree node inside browse tree"))));
      await session.step(110, "Then \"BDDAugment{time} > iris\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDDAugment{time} > iris\" tree node inside browse tree")), "visible"));
    });
    await run.scenario("The node's SAVE saves the project with iris, its Data sync switched on", async () => {
      await session.step(113, "When user clicks on Save button in BDDAugment{time} tree node inside browse tree", () => clickOn(page, el(session.text("Save button in BDDAugment{time} tree node inside browse tree"))));
      await session.step(114, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(115, "And \"iris\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"iris\" project table in \"Save project\" dialog"), "visible"));
      await session.step(116, "And \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(117, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(118, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(119, "And Data sync switch in \"iris\" project table in \"Save project\" dialog should be unchecked", () => shouldBe(page, el("Data sync switch in \"iris\" project table in \"Save project\" dialog"), "unchecked"));
      await session.step(120, "When user checks Data sync switch in \"iris\" project table in \"Save project\" dialog", () => check(page, el("Data sync switch in \"iris\" project table in \"Save project\" dialog")));
      await session.step(121, "And user clicks on \"Creation script\" button in \"iris\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"iris\" project table in \"Save project\" dialog")));
      await session.step(122, "Then \"iris\" project table in \"Save project\" dialog should contain text 'OpenFile(\"System:DemoFiles/iris.csv\")'", () => shouldContainText(page, el("\"iris\" project table in \"Save project\" dialog"), "OpenFile(\"System:DemoFiles/iris.csv\")"));
      await session.step(123, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(124, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(125, "And the \"iris\" table of the \"BDDAugment{time}\" project should be saved with data sync", () => savedWithDataSync(page, "iris", session.text("BDDAugment{time}")));
      await session.step(126, "And the \"demog\" table of the \"BDDAugment{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDDAugment{time}")));
      await session.step(128, "When user clicks on Dashboards tab", () => clickOn(page, el("Dashboards tab")));
      await session.step(129, "Then \"New Dashboard\" tree node inside browse tree should be hidden", () => shouldBe(page, el("\"New Dashboard\" tree node inside browse tree"), "hidden"));
    });
    await run.scenario("Reopened from Dashboards, the project holds and opens both tables", async () => {
      await session.step(132, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(133, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(134, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(135, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(136, "And user enters \"BDDAugment{time}\" into gallery search", () => enterInto(page, session.text("BDDAugment{time}"), el("gallery search")));
      await session.step(137, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(138, "And user clicks on BDDAugment{time} gallery card", () => clickOn(page, el(session.text("BDDAugment{time} gallery card"))));
      await session.step(139, "Then the context panel should show \"BDDAugment{time}\"", () => contextPanelShows(page, session.text("BDDAugment{time}")));
      await session.step(140, "Given Content pane in context panel is expanded", () => isExpanded(page, el("Content pane in context panel")));
      await session.step(141, "Then Content pane in context panel should contain text \"demog\"", () => shouldContainText(page, el("Content pane in context panel"), "demog"));
      await session.step(142, "And Content pane in context panel should contain text \"iris\"", () => shouldContainText(page, el("Content pane in context panel"), "iris"));
      await session.step(143, "When user double-clicks on BDDAugment{time} gallery card", () => doubleClickOn(page, el(session.text("BDDAugment{time} gallery card"))));
      await session.step(144, "Then the \"demog\" table view should open with 5850 rows", () => tableViewOpened(page, "demog", 5850));
      await session.step(145, "And the \"iris\" table view should open with 150 rows", () => tableViewOpened(page, "iris", 150));
      await session.step(146, "When user clicks on iris view", () => clickOn(page, el("iris view")));
      await session.step(147, "Then the \"iris\" view should be current", () => viewIsCurrent(page, "iris"));
      await session.step(148, "And status bar should contain text \"Rows: 150\"", () => shouldContainText(page, el("status bar"), "Rows: 150"));
      await session.step(149, "And the table should have 6 columns", () => columnCount(page, 6));
      await session.step(150, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(151, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(152, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(153, "And the Save project dialog should save the tables \"demog, iris\" with data sync", () => saveDialogDataSync(page, "demog, iris"));
      await session.step(154, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(155, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    });
    await run.scenario("The project is deleted from its card", async () => {
      await session.step(158, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(159, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(160, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(161, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(162, "And user enters \"BDDAugment{time}\" into gallery search", () => enterInto(page, session.text("BDDAugment{time}"), el("gallery search")));
      await session.step(163, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(164, "And user picks \"Delete Project\" from the context menu of BDDAugment{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDAugment{time} gallery card"))));
      await session.step(165, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(166, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(167, "And 0 projects named \"BDDAugment{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDAugment{time}")));
      await session.step(168, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(169, "Then BDDAugment{time} gallery card should be absent", () => shouldBe(page, el(session.text("BDDAugment{time} gallery card")), "absent"));
    });
    run.finish();
  });
});
