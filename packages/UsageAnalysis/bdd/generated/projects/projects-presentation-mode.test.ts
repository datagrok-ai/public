/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-presentation-mode.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects]
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
import {clickOn, doubleClickOn, enterInto, hoverOver, isExpanded, pressKey, shouldBe, shouldBeSwitchedOn, shouldContainText, switchOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, dialogCloses, noProjectOnServer, projectsOnServer, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Presentation mode switch of a saved project's Save dialog", () => {
  const session = feature(test, "features/projects/projects-presentation-mode.feature", import.meta.url);
  test("The Presentation mode switch of a saved project's Save dialog", {tag: ["@journey", "@serial", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "And no project named \"BDDPresentProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDPresentProj{time}")));
    await session.step(34, "And no project named \"BDDPresentPlain{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDPresentPlain{time}")));
    await session.step(35, "And no project named \"BDDPresentNew{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDPresentNew{time}")));
    await run.scenario("demog with a scatter plot is saved as a project", async () => {
      await session.step(38, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(39, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(40, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(41, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(42, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(43, "When user clicks on \"scatter plot\" icon in toolbox", () => clickOn(page, el("\"scatter plot\" icon in toolbox")));
      await session.step(44, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(45, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(46, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(47, "When user enters \"BDDPresentProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDPresentProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(48, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(49, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(50, "And an info balloon containing 'Project \"BDDPresentProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDPresentProj{time}\" uploaded")));
      await session.step(51, "And \"Share BDDPresentProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDPresentProj{time}\" dialog")), "visible"));
      await session.step(52, "When user clicks on CANCEL button in \"Share BDDPresentProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDPresentProj{time}\" dialog"))));
      await session.step(53, "Then the \"Share BDDPresentProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDPresentProj{time}")));
    });
    await run.scenario("A second project is saved in design mode", async () => {
      await session.step(56, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(57, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(58, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(59, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(60, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(61, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(62, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(63, "When user enters \"BDDPresentPlain{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDPresentPlain{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(64, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(65, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(66, "And an info balloon containing 'Project \"BDDPresentPlain{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDPresentPlain{time}\" uploaded")));
      await session.step(67, "When user clicks on CANCEL button in \"Share BDDPresentPlain{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDPresentPlain{time}\" dialog"))));
      await session.step(68, "Then the \"Share BDDPresentPlain{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDPresentPlain{time}")));
    });
    await run.scenario("The project opens in design mode", async () => {
      await session.step(71, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(72, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(73, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(74, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(75, "And user enters \"BDDPresent\" into gallery search", () => enterInto(page, "BDDPresent", el("gallery search")));
      await session.step(76, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(77, "And user double-clicks on BDDPresentProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDPresentProj{time} gallery card"))));
      await session.step(78, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(79, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(80, "And \"Data\" menu item should be visible", () => shouldBe(page, el("\"Data\" menu item"), "visible"));
      await session.step(81, "And browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
      await session.step(82, "And status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(83, "And toolbox tab should be visible", () => shouldBe(page, el("toolbox tab"), "visible"));
    });
    await run.scenario("The Save dialog offers Presentation mode, and its tooltip tells how to get back", async () => {
      await session.step(86, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(87, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(88, "And \"Presentation mode\" input in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Presentation mode\" input in \"Save project\" dialog"), "visible"));
      await session.step(89, "When user hovers over \"Presentation mode\" input in \"Save project\" dialog", () => hoverOver(page, el("\"Presentation mode\" input in \"Save project\" dialog")));
      await session.step(90, "Then tooltip should contain text \"visualization\"", () => shouldContainText(page, el("tooltip"), "visualization"));
      await session.step(91, "And tooltip should contain text \"Design mode\"", () => shouldContainText(page, el("tooltip"), "Design mode"));
      await session.step(92, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(93, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    });
    await run.scenario("Presentation mode is turned on and saved into the project", async () => {
      await session.step(96, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(97, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(98, "When user switches on \"Presentation mode\" input in \"Save project\" dialog", () => switchOn(page, el("\"Presentation mode\" input in \"Save project\" dialog")));
      await session.step(99, "Then \"Presentation mode\" input in \"Save project\" dialog should be switched on", () => shouldBeSwitchedOn(page, el("\"Presentation mode\" input in \"Save project\" dialog")));
      await session.step(100, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(101, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(102, "And an info balloon containing 'Project \"BDDPresentProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDPresentProj{time}\" uploaded")));
      await session.step(103, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(104, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    });
    await run.scenario("The project reopens in presentation mode", async () => {
      await session.step(107, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(108, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(109, "And user enters \"BDDPresent\" into gallery search", () => enterInto(page, "BDDPresent", el("gallery search")));
      await session.step(110, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(111, "And user double-clicks on BDDPresentProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDPresentProj{time} gallery card"))));
      await session.step(112, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(113, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(114, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(115, "And browse tab should be hidden", () => shouldBe(page, el("browse tab"), "hidden"));
      await session.step(116, "And toolbox tab should be hidden", () => shouldBe(page, el("toolbox tab"), "hidden"));
      await session.step(117, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(118, "And an info balloon containing \"Press F7 to go back to the design mode\" should have been shown", () => infoBalloonText(page, "Press F7 to go back to the design mode"));
    });
    await run.scenario("back to design mode and F7 switch between the modes", async () => {
      await session.step(121, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")));
      await session.step(122, "Then browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
      await session.step(123, "And status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(124, "And \"back to design mode\" link should be absent", () => shouldBe(page, el("\"back to design mode\" link"), "absent"));
      await session.step(125, "And an info balloon containing \"Press F7 to go back to the presentation mode\" should have been shown", () => infoBalloonText(page, "Press F7 to go back to the presentation mode"));
      await session.step(126, "When user presses F7", () => pressKey(page, "F7"));
      await session.step(127, "Then status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(128, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(129, "When user presses F7", () => pressKey(page, "F7"));
      await session.step(130, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(131, "And \"back to design mode\" link should be absent", () => shouldBe(page, el("\"back to design mode\" link"), "absent"));
    });
    await run.scenario("A new project saved from the Dashboards panel with the switch on opens in presentation mode", async () => {
      await session.step(134, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(135, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(136, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(137, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(138, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(139, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(140, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(141, "When user clicks on Dashboards tab", () => clickOn(page, el("Dashboards tab")));
      await session.step(142, "Then \"New Dashboard > demog\" tree node inside browse tree should be visible", () => shouldBe(page, el("\"New Dashboard > demog\" tree node inside browse tree"), "visible"));
      await session.step(143, "When user clicks on Save button in \"New Dashboard\" tree node inside browse tree", () => clickOn(page, el("Save button in \"New Dashboard\" tree node inside browse tree")));
      await session.step(144, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(145, "When user enters \"BDDPresentNew{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDPresentNew{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(146, "And user switches on \"Presentation mode\" input in \"Save project\" dialog", () => switchOn(page, el("\"Presentation mode\" input in \"Save project\" dialog")));
      await session.step(147, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(148, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(149, "And an info balloon containing 'Project \"BDDPresentNew{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDPresentNew{time}\" uploaded")));
      await session.step(150, "When user clicks on OK button in \"Share BDDPresentNew{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDPresentNew{time}\" dialog"))));
      await session.step(151, "Then the \"Share BDDPresentNew{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDPresentNew{time}")));
      await session.step(153, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(154, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(155, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")));
      await session.step(156, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(158, "When user clicks on Dashboards tab", () => clickOn(page, el("Dashboards tab")));
      await session.step(159, "Then \"New Dashboard\" tree node inside browse tree should be hidden", () => shouldBe(page, el("\"New Dashboard\" tree node inside browse tree"), "hidden"));
      await session.step(160, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(161, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(162, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(163, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(164, "And user enters \"BDDPresent\" into gallery search", () => enterInto(page, "BDDPresent", el("gallery search")));
      await session.step(165, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(166, "Then BDDPresentNew{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDPresentNew{time} gallery card")), "visible"));
      await session.step(167, "When user double-clicks on BDDPresentNew{time} gallery card", () => doubleClickOn(page, el(session.text("BDDPresentNew{time} gallery card"))));
      await session.step(168, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(169, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(170, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(171, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")));
      await session.step(172, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
    });
    await run.scenario("The project's address opens it in presentation mode, and ?mode=presentation any project", async () => {
      await session.step(175, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(176, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(177, "When user opens the address \"/p/Admin.BDDPresentProj{time}\"", () => openAddress(page, session.text("/p/Admin.BDDPresentProj{time}")));
      await session.step(178, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(179, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(180, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(181, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")));
      await session.step(182, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(183, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(184, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(185, "When user opens the address \"/p/Admin.BDDPresentPlain{time}\"", () => openAddress(page, session.text("/p/Admin.BDDPresentPlain{time}")));
      await session.step(186, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(187, "And grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(188, "And scatter plot viewer should be absent", () => shouldBe(page, el("scatter plot viewer"), "absent"));
      await session.step(189, "And status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(190, "And \"back to design mode\" link should be absent", () => shouldBe(page, el("\"back to design mode\" link"), "absent"));
      await session.step(191, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(192, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(193, "When user opens the address \"/p/Admin.BDDPresentPlain{time}?mode=presentation\"", () => openAddress(page, session.text("/p/Admin.BDDPresentPlain{time}?mode=presentation")));
      await session.step(194, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(195, "And scatter plot viewer should be absent", () => shouldBe(page, el("scatter plot viewer"), "absent"));
      await session.step(196, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(197, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(198, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")));
      await session.step(199, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
    });
    await run.scenario("The projects are deleted from their cards", async () => {
      await session.step(202, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(203, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(204, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(205, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(206, "And user enters \"BDDPresent\" into gallery search", () => enterInto(page, "BDDPresent", el("gallery search")));
      await session.step(207, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(208, "And user picks \"Delete Project\" from the context menu of BDDPresentProj{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDPresentProj{time} gallery card"))));
      await session.step(209, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(210, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(211, "When user picks \"Delete Project\" from the context menu of BDDPresentPlain{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDPresentPlain{time} gallery card"))));
      await session.step(212, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(213, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(214, "When user picks \"Delete Project\" from the context menu of BDDPresentNew{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDPresentNew{time} gallery card"))));
      await session.step(215, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(216, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(217, "And 0 projects named \"BDDPresentProj{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDPresentProj{time}")));
      await session.step(218, "And 0 projects named \"BDDPresentPlain{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDPresentPlain{time}")));
      await session.step(219, "And 0 projects named \"BDDPresentNew{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDPresentNew{time}")));
    });
    run.finish();
  });
});
