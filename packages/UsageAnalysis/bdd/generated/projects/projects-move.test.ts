/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-move.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, views.space, sharing.share-dialog]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/nx.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, dragTo, enterInto, isExpanded, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, dialogCloses, noProjectOnServer, noSpaceOnServer, pickSharingUser, projectsOnServer, reloadedByDataSync, spacesOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project moved between spaces, and access through a shared space", () => {
  const session = feature(test, "features/projects/projects-move.feature", import.meta.url);
  test("A project moved between spaces, and access through a shared space", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:views.space", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no project named \"BDDMoveProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDMoveProj{time}")));
    await session.step(26, "And no space named \"BDDMoveA{time}, BDDMoveAChild{time}, BDDMoveB{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDMoveA{time}, BDDMoveAChild{time}, BDDMoveB{time}")));
    await run.scenario("Two spaces and a child of the first are created", async () => {
      await session.step(29, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(30, "Then Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(31, "When user enters \"BDDMoveA{time}\" into Name input in Create Space dialog", () => enterInto(page, session.text("BDDMoveA{time}"), el("Name input in Create Space dialog")));
      await session.step(32, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(33, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(34, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(35, "Then Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(36, "When user enters \"BDDMoveB{time}\" into Name input in Create Space dialog", () => enterInto(page, session.text("BDDMoveB{time}"), el("Name input in Create Space dialog")));
      await session.step(37, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(38, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(39, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(40, "Then BDDMoveA{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDMoveA{time} tree node inside browse tree")), "visible"));
      await session.step(41, "And BDDMoveB{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDMoveB{time} tree node inside browse tree")), "visible"));
      await session.step(42, "When user picks \"Create Child Space...\" from the context menu of BDDMoveA{time} tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el(session.text("BDDMoveA{time} tree node inside browse tree"))));
      await session.step(43, "Then Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(44, "When user enters \"BDDMoveAChild{time}\" into Name input in Create Space dialog", () => enterInto(page, session.text("BDDMoveAChild{time}"), el("Name input in Create Space dialog")));
      await session.step(45, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(46, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(47, "And BDDMoveAChild{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDMoveAChild{time} tree node inside browse tree")), "visible"));
      await session.step(48, "And 1 space named \"BDDMoveA{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDMoveA{time}")));
      await session.step(49, "And 1 space named \"BDDMoveB{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDMoveB{time}")));
    });
    await run.scenario("A file is saved as a project with Data sync", async () => {
      await session.step(52, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(53, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(54, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(55, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(56, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(57, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(58, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(59, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(60, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(61, "When user enters \"BDDMoveProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDMoveProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(62, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(63, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(64, "And an info balloon containing 'Project \"BDDMoveProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDMoveProj{time}\" uploaded")));
      await session.step(65, "And \"Share BDDMoveProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDMoveProj{time}\" dialog")), "visible"));
      await session.step(66, "When user clicks on CANCEL button in \"Share BDDMoveProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDMoveProj{time}\" dialog"))));
      await session.step(67, "Then the \"Share BDDMoveProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDMoveProj{time}")));
      await session.step(68, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(69, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The project is moved into the first space from its card and opens there", async () => {
      await session.step(73, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(74, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(75, "And user enters \"BDDMoveProj{time}\" into gallery search", () => enterInto(page, session.text("BDDMoveProj{time}"), el("gallery search")));
      await session.step(76, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(77, "And user picks \"Move to Space...\" from the context menu of BDDMoveProj{time} gallery card", () => pickFromContextMenu(page, "Move to Space...", el(session.text("BDDMoveProj{time} gallery card"))));
      await session.step(78, "Then \"Move to space\" dialog should be visible", () => shouldBe(page, el("\"Move to space\" dialog"), "visible"));
      await session.step(79, "When user selects \"BDDMoveA{time}\" in Space input in \"Move to space\" dialog", () => selectIn(page, session.text("BDDMoveA{time}"), el("Space input in \"Move to space\" dialog")));
      await session.step(80, "And user clicks on OK button in \"Move to space\" dialog", () => clickOn(page, el("OK button in \"Move to space\" dialog")));
      await session.step(81, "Then the \"Move to space\" dialog should close", () => dialogCloses(page, "Move to space"));
      await session.step(82, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(83, "When user double-clicks on BDDMoveA{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveA{time} tree node inside browse tree"))));
      await session.step(84, "Then the \"BDDMoveA{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveA{time}")));
      await session.step(85, "And \"BDDMoveProj{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDMoveProj{time}\" link in gallery")), "visible"));
      await session.step(86, "When user double-clicks on \"BDDMoveProj{time}\" link in gallery", () => doubleClickOn(page, el(session.text("\"BDDMoveProj{time}\" link in gallery"))));
      await session.step(87, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(88, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(89, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(90, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(91, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(92, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    });
    await run.scenario("The space is shared, and the project's access is inherited from it", async () => {
      await session.step(95, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(96, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(97, "When user picks \"Share...\" from the context menu of BDDMoveA{time} tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el(session.text("BDDMoveA{time} tree node inside browse tree"))));
      await session.step(98, "Then \"Share BDDMoveA{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDMoveA{time}\" dialog")), "visible"));
      await session.step(99, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(100, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDMoveA{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDMoveA{time}\" dialog"))));
      await session.step(101, "And user clicks on OK button in \"Share BDDMoveA{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDMoveA{time}\" dialog"))));
      await session.step(102, "Then the \"Share BDDMoveA{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDMoveA{time}")));
      await session.step(103, "When user double-clicks on BDDMoveA{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveA{time} tree node inside browse tree"))));
      await session.step(104, "Then the \"BDDMoveA{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveA{time}")));
      await session.step(105, "When user picks \"Share...\" from the context menu of \"BDDMoveProj{time}\" link in gallery", () => pickFromContextMenu(page, "Share...", el(session.text("\"BDDMoveProj{time}\" link in gallery"))));
      await session.step(106, "Then \"Share BDDMoveProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDMoveProj{time}\" dialog")), "visible"));
      await session.step(107, "And \"Share BDDMoveProj{time}\" dialog should contain text \"Inherited from\"", () => shouldContainText(page, el(session.text("\"Share BDDMoveProj{time}\" dialog")), "Inherited from"));
      await session.step(108, "And \"Share BDDMoveProj{time}\" dialog should contain text \"BDDMoveA{time}\"", () => shouldContainText(page, el(session.text("\"Share BDDMoveProj{time}\" dialog")), session.text("BDDMoveA{time}")));
      await session.step(109, "When user clicks on CANCEL button in \"Share BDDMoveProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDMoveProj{time}\" dialog"))));
      await session.step(110, "Then the \"Share BDDMoveProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDMoveProj{time}")));
    });
    await run.scenario("The project is dragged into the child space and opens there", async () => {
      await session.step(113, "When user drags \"BDDMoveProj{time}\" link in gallery to BDDMoveAChild{time} tree node inside browse tree", () => dragTo(page, el(session.text("\"BDDMoveProj{time}\" link in gallery")), el(session.text("BDDMoveAChild{time} tree node inside browse tree"))));
      await session.step(114, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(115, "When user selects \"Move\" in Move entity dialog", () => selectIn(page, "Move", el("Move entity dialog")));
      await session.step(116, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(117, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(118, "When user double-clicks on BDDMoveAChild{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveAChild{time} tree node inside browse tree"))));
      await session.step(119, "Then the \"BDDMoveAChild{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveAChild{time}")));
      await session.step(120, "And \"BDDMoveProj{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDMoveProj{time}\" link in gallery")), "visible"));
      await session.step(121, "When user double-clicks on \"BDDMoveProj{time}\" link in gallery", () => doubleClickOn(page, el(session.text("\"BDDMoveProj{time}\" link in gallery"))));
      await session.step(122, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(123, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(124, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(125, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(126, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    });
    await run.scenario("The project is moved on to the second space and opens there", async () => {
      await session.step(129, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(130, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(131, "When user double-clicks on BDDMoveAChild{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveAChild{time} tree node inside browse tree"))));
      await session.step(132, "Then the \"BDDMoveAChild{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveAChild{time}")));
      await session.step(133, "When user picks \"Move to Space...\" from the context menu of \"BDDMoveProj{time}\" link in gallery", () => pickFromContextMenu(page, "Move to Space...", el(session.text("\"BDDMoveProj{time}\" link in gallery"))));
      await session.step(134, "Then \"Move to space\" dialog should be visible", () => shouldBe(page, el("\"Move to space\" dialog"), "visible"));
      await session.step(135, "When user selects \"BDDMoveB{time}\" in Space input in \"Move to space\" dialog", () => selectIn(page, session.text("BDDMoveB{time}"), el("Space input in \"Move to space\" dialog")));
      await session.step(136, "And user clicks on OK button in \"Move to space\" dialog", () => clickOn(page, el("OK button in \"Move to space\" dialog")));
      await session.step(137, "Then the \"Move to space\" dialog should close", () => dialogCloses(page, "Move to space"));
      await session.step(138, "When user double-clicks on BDDMoveB{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveB{time} tree node inside browse tree"))));
      await session.step(139, "Then the \"BDDMoveB{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveB{time}")));
      await session.step(140, "And \"BDDMoveProj{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDMoveProj{time}\" link in gallery")), "visible"));
      await session.step(141, "When user double-clicks on BDDMoveAChild{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveAChild{time} tree node inside browse tree"))));
      await session.step(142, "Then the \"BDDMoveAChild{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveAChild{time}")));
      await session.step(143, "And \"BDDMoveProj{time}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BDDMoveProj{time}\" link in gallery")), "absent"));
      await session.step(144, "When user double-clicks on BDDMoveB{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveB{time} tree node inside browse tree"))));
      await session.step(145, "Then the \"BDDMoveB{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveB{time}")));
      await session.step(146, "When user double-clicks on \"BDDMoveProj{time}\" link in gallery", () => doubleClickOn(page, el(session.text("\"BDDMoveProj{time}\" link in gallery"))));
      await session.step(147, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(148, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(149, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(150, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(151, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    });
    await run.scenario("The project and the spaces are deleted", async () => {
      await session.step(154, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(155, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(156, "When user double-clicks on BDDMoveB{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveB{time} tree node inside browse tree"))));
      await session.step(157, "Then the \"BDDMoveB{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveB{time}")));
      await session.step(158, "When user picks \"Delete Project\" from the context menu of \"BDDMoveProj{time}\" link in gallery", () => pickFromContextMenu(page, "Delete Project", el(session.text("\"BDDMoveProj{time}\" link in gallery"))));
      await session.step(159, "Then \"Are you sure?\" dialog should contain text 'Delete project \"BDDMoveProj{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDMoveProj{time}\"?")));
      await session.step(160, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(161, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(162, "And 0 projects named \"BDDMoveProj{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDMoveProj{time}")));
      await session.step(163, "When user picks \"Delete Space\" from the context menu of BDDMoveA{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el(session.text("BDDMoveA{time} tree node inside browse tree"))));
      await session.step(164, "Then \"Are you sure?\" dialog should contain text 'Delete space \"BDDMoveA{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete space \"BDDMoveA{time}\"?")));
      await session.step(165, "And \"Are you sure?\" dialog should contain text \"This will delete space and its related data\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "This will delete space and its related data"));
      await session.step(166, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(167, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(168, "When user picks \"Delete Space\" from the context menu of BDDMoveB{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el(session.text("BDDMoveB{time} tree node inside browse tree"))));
      await session.step(169, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(170, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(171, "And 0 spaces named \"BDDMoveA{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDMoveA{time}")));
      await session.step(172, "And 0 spaces named \"BDDMoveB{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDMoveB{time}")));
      await session.step(173, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
