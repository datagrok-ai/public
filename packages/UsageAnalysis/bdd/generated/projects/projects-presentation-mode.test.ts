/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-presentation-mode.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects]
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
import {clickOn, doubleClickOn, enterInto, hoverOver, isExpanded, pressKey, shouldBe, shouldBeSwitchedOn, shouldContainText, switchOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, dialogCloses, doubleClickOpensProject, noProjectOnServer, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Presentation mode switch of a saved project's Save dialog", () => {
  const session = feature(test, "features/projects/projects-presentation-mode.feature", import.meta.url);
  test("The Presentation mode switch of a saved project's Save dialog", {tag: ["@journey", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(34, "And no project named \"BDDPresentProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDPresentProj{time}")));
    await session.step(35, "And no project named \"BDDPresentPlain{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDPresentPlain{time}")));
    await session.step(36, "And no project named \"BDDPresentNew{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDPresentNew{time}")));
    await run.scenario("demog with a scatter plot is saved as a project", async () => {
      await session.step(39, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(40, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(41, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(42, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(43, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(44, "When user clicks on \"scatter plot\" icon in toolbox", () => clickOn(page, el("\"scatter plot\" icon in toolbox")));
      await session.step(45, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(46, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(47, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(48, "When user enters \"BDDPresentProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDPresentProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(49, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(50, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(51, "And an info balloon containing 'Project \"BDDPresentProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDPresentProj{time}\" uploaded")));
      await session.step(52, "And \"Share BDDPresentProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDPresentProj{time}\" dialog")), "visible"));
      await session.step(53, "When user clicks on CANCEL button in \"Share BDDPresentProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDPresentProj{time}\" dialog"))));
      await session.step(54, "Then the \"Share BDDPresentProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDPresentProj{time}")));
    });
    await run.scenario("A second project is saved in design mode", async () => {
      await session.step(57, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(58, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(59, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(60, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(61, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(62, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(63, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(64, "When user enters \"BDDPresentPlain{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDPresentPlain{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(65, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(66, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(67, "And an info balloon containing 'Project \"BDDPresentPlain{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDPresentPlain{time}\" uploaded")));
      await session.step(68, "When user clicks on CANCEL button in \"Share BDDPresentPlain{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDPresentPlain{time}\" dialog"))));
      await session.step(69, "Then the \"Share BDDPresentPlain{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDPresentPlain{time}")));
    });
    await run.scenario("The project opens in design mode", async () => {
      await session.step(72, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(73, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(74, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(75, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(76, "And user enters \"BDDPresent\" into gallery search", () => enterInto(page, "BDDPresent", el("gallery search")));
      await session.step(77, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(78, "And user double-clicks on BDDPresentProj{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDPresentProj{time} gallery card"))));
      await session.step(79, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(80, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(81, "And \"Data\" menu item should be visible", () => shouldBe(page, el("\"Data\" menu item"), "visible"));
      await session.step(82, "And browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
      await session.step(83, "And status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(84, "And toolbox tab should be visible", () => shouldBe(page, el("toolbox tab"), "visible"));
    });
    await run.scenario("The Save dialog offers Presentation mode, and its tooltip tells how to get back", async () => {
      await session.step(87, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(88, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(89, "And \"Presentation mode\" input in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Presentation mode\" input in \"Save project\" dialog"), "visible"));
      await session.step(90, "When user hovers over \"Presentation mode\" input in \"Save project\" dialog", () => hoverOver(page, el("\"Presentation mode\" input in \"Save project\" dialog")));
      await session.step(91, "Then tooltip should contain text \"visualization\"", () => shouldContainText(page, el("tooltip"), "visualization"));
      await session.step(92, "And tooltip should contain text \"Design mode\"", () => shouldContainText(page, el("tooltip"), "Design mode"));
      await session.step(93, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(94, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    });
    await run.scenario("Presentation mode is turned on and saved into the project", async () => {
      await session.step(97, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(98, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(99, "When user switches on \"Presentation mode\" input in \"Save project\" dialog", () => switchOn(page, el("\"Presentation mode\" input in \"Save project\" dialog")));
      await session.step(100, "Then \"Presentation mode\" input in \"Save project\" dialog should be switched on", () => shouldBeSwitchedOn(page, el("\"Presentation mode\" input in \"Save project\" dialog")));
      await session.step(101, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(102, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(103, "And an info balloon containing 'Project \"BDDPresentProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDPresentProj{time}\" uploaded")));
      await session.step(104, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(105, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    });
    await run.scenario("The project reopens in presentation mode", async () => {
      await session.step(108, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(109, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(110, "And user enters \"BDDPresent\" into gallery search", () => enterInto(page, "BDDPresent", el("gallery search")));
      await session.step(111, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(112, "And user double-clicks on BDDPresentProj{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDPresentProj{time} gallery card"))));
      await session.step(113, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(114, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(115, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(116, "And browse tab should be hidden", () => shouldBe(page, el("browse tab"), "hidden"));
      await session.step(117, "And toolbox tab should be hidden", () => shouldBe(page, el("toolbox tab"), "hidden"));
      await session.step(118, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(119, "And an info balloon containing \"Press F7 to go back to the design mode\" should have been shown", () => infoBalloonText(page, "Press F7 to go back to the design mode"));
    });
    await run.scenario("back to design mode and F7 switch between the modes", async () => {
      await session.step(122, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")));
      await session.step(123, "Then browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
      await session.step(124, "And status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(125, "And \"back to design mode\" link should be absent", () => shouldBe(page, el("\"back to design mode\" link"), "absent"));
      await session.step(126, "And an info balloon containing \"Press F7 to go back to the presentation mode\" should have been shown", () => infoBalloonText(page, "Press F7 to go back to the presentation mode"));
      await session.step(127, "When user presses F7", () => pressKey(page, "F7"));
      await session.step(128, "Then status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(129, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(130, "When user presses F7", () => pressKey(page, "F7"));
      await session.step(131, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(132, "And \"back to design mode\" link should be absent", () => shouldBe(page, el("\"back to design mode\" link"), "absent"));
    });
    await run.scenario("A new project saved from the Dashboards panel with the switch on opens in presentation mode", async () => {
      await session.step(135, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(136, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(137, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(138, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(139, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(140, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(141, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(142, "When user clicks on Dashboards tab", () => clickOn(page, el("Dashboards tab")));
      await session.step(143, "Then \"New Dashboard > demog\" tree node inside browse tree should be visible", () => shouldBe(page, el("\"New Dashboard > demog\" tree node inside browse tree"), "visible"));
      await session.step(144, "When user clicks on Save button in \"New Dashboard\" tree node inside browse tree", () => clickOn(page, el("Save button in \"New Dashboard\" tree node inside browse tree")));
      await session.step(145, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(146, "When user enters \"BDDPresentNew{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDPresentNew{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(147, "And user switches on \"Presentation mode\" input in \"Save project\" dialog", () => switchOn(page, el("\"Presentation mode\" input in \"Save project\" dialog")));
      await session.step(148, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(149, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(150, "And an info balloon containing 'Project \"BDDPresentNew{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDPresentNew{time}\" uploaded")));
      await session.step(151, "When user clicks on OK button in \"Share BDDPresentNew{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDPresentNew{time}\" dialog"))));
      await session.step(152, "Then the \"Share BDDPresentNew{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDPresentNew{time}")));
      await session.step(154, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(155, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(156, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")));
      await session.step(157, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(158, "And BDDPresentNew{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDPresentNew{time} tree node inside browse tree")), "visible"));
      await session.step(159, "And \"New Dashboard > demog\" tree node inside browse tree should be absent", () => shouldBe(page, el("\"New Dashboard > demog\" tree node inside browse tree"), "absent"));
      await session.step(161, "When user clicks on Dashboards tab", () => clickOn(page, el("Dashboards tab")));
      await session.step(162, "Then \"New Dashboard\" tree node inside browse tree should be hidden", () => shouldBe(page, el("\"New Dashboard\" tree node inside browse tree"), "hidden"));
      await session.step(163, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(164, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(165, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(166, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(167, "And user enters \"BDDPresent\" into gallery search", () => enterInto(page, "BDDPresent", el("gallery search")));
      await session.step(168, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(169, "Then BDDPresentNew{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDPresentNew{time} gallery card")), "visible"));
      await session.step(170, "When user double-clicks on BDDPresentNew{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDPresentNew{time} gallery card"))));
      await session.step(171, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(172, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(173, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(174, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")));
      await session.step(175, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
    });
    await run.scenario("The project's address opens it in presentation mode, and ?mode=presentation any project", async () => {
      await session.step(178, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(179, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(180, "When user opens the address \"/p/Admin.BDDPresentProj{time}\"", () => openAddress(page, session.text("/p/Admin.BDDPresentProj{time}")));
      await session.step(181, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(182, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(183, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(184, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")));
      await session.step(185, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(186, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(187, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(188, "When user opens the address \"/p/Admin.BDDPresentPlain{time}\"", () => openAddress(page, session.text("/p/Admin.BDDPresentPlain{time}")));
      await session.step(189, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(190, "And grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(191, "And scatter plot viewer should be absent", () => shouldBe(page, el("scatter plot viewer"), "absent"));
      await session.step(192, "And status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(193, "And \"back to design mode\" link should be absent", () => shouldBe(page, el("\"back to design mode\" link"), "absent"));
      await session.step(194, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(195, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(196, "When user opens the address \"/p/Admin.BDDPresentPlain{time}?mode=presentation\"", () => openAddress(page, session.text("/p/Admin.BDDPresentPlain{time}?mode=presentation")));
      await session.step(197, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(198, "And scatter plot viewer should be absent", () => shouldBe(page, el("scatter plot viewer"), "absent"));
      await session.step(199, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
      await session.step(200, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
      await session.step(201, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")));
      await session.step(202, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
    });
    run.finish();
  });
});
