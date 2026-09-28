/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-move.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, views.space]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/minimized-viewers.js';
import '../../bindings/projects-copies.js';
import '../../bindings/projects-derived.js';
import '../../bindings/projects-regressions.js';
import '../../bindings/projects-sources.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, dialogCloses, noProjectOnServer, noSpaceOnServer, projectsOnServer, spaceHoldsEntity, spaceHoldsNotEntity, spacesOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noErrorBalloon, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project moved into a space and on to another", () => {
  const session = feature(test, "features/projects/projects-move.feature", import.meta.url);
  test("A project moved into a space and on to another", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "And no project named \"BDDMove{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDMove{time}")));
    await session.step(37, "And no space named \"BDDMoveA{time}, BDDMoveB{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDMoveA{time}, BDDMoveB{time}")));
    await run.scenario("Two spaces are created", async () => {
      await session.step(40, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(41, "Then Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(42, "When user enters \"BDDMoveA{time}\" into Name input in Create Space dialog", () => enterInto(page, session.text("BDDMoveA{time}"), el("Name input in Create Space dialog")));
      await session.step(43, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(44, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(45, "And 1 space named \"BDDMoveA{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDMoveA{time}")));
      await session.step(46, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(47, "And user enters \"BDDMoveB{time}\" into Name input in Create Space dialog", () => enterInto(page, session.text("BDDMoveB{time}"), el("Name input in Create Space dialog")));
      await session.step(48, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(49, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(50, "And 1 space named \"BDDMoveB{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDMoveB{time}")));
      await session.step(51, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(52, "Then BDDMoveA{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDMoveA{time} tree node inside browse tree")), "visible"));
      await session.step(53, "And BDDMoveB{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDMoveB{time} tree node inside browse tree")), "visible"));
    });
    await run.scenario("The file is saved as a project", async () => {
      await session.step(56, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(57, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(58, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(59, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(60, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(61, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(62, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(63, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(64, "When user enters \"BDDMove{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDMove{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(65, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(66, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(67, "And an info balloon containing 'Project \"BDDMove{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDMove{time}\" uploaded")));
      await session.step(68, "And 1 project named \"BDDMove{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDMove{time}")));
      await session.step(69, "And \"Share BDDMove{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDMove{time}\" dialog")), "visible"));
      await session.step(70, "When user clicks on CANCEL button in \"Share BDDMove{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDMove{time}\" dialog"))));
      await session.step(71, "Then the \"Share BDDMove{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDMove{time}")));
      await session.step(72, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(73, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The project is moved into the first space from its card", async () => {
      await session.step(77, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(78, "And the \"BDDMoveA{time}\" space should not hold the project \"BDDMove{time}\" on the server", () => spaceHoldsNotEntity(page, session.text("BDDMoveA{time}"), "project", session.text("BDDMove{time}")));
      await session.step(79, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(80, "And user enters \"BDDMove{time}\" into gallery search", () => enterInto(page, session.text("BDDMove{time}"), el("gallery search")));
      await session.step(81, "And user picks \"Move to Space...\" from the context menu of BDDMove{time} gallery card", () => pickFromContextMenu(page, "Move to Space...", el(session.text("BDDMove{time} gallery card"))));
      await session.step(82, "Then Move to space dialog should be visible", () => shouldBe(page, el("Move to space dialog"), "visible"));
      await session.step(83, "When user selects \"BDDMoveA{time}\" in Space input in Move to space dialog", () => selectIn(page, session.text("BDDMoveA{time}"), el("Space input in Move to space dialog")));
      await session.step(84, "And user clicks on OK button in Move to space dialog", () => clickOn(page, el("OK button in Move to space dialog")));
      await session.step(85, "Then the \"Move to space\" dialog should close", () => dialogCloses(page, "Move to space"));
      await session.step(86, "And an info balloon containing \"Moved BDDMove{time} to BDDMoveA{time}\" should have been shown", () => infoBalloonText(page, session.text("Moved BDDMove{time} to BDDMoveA{time}")));
      await session.step(87, "And the \"BDDMoveA{time}\" space should hold the project \"BDDMove{time}\" on the server", () => spaceHoldsEntity(page, session.text("BDDMoveA{time}"), "project", session.text("BDDMove{time}")));
      await session.step(88, "When user double-clicks on BDDMoveA{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveA{time} tree node inside browse tree"))));
      await session.step(89, "Then the \"BDDMoveA{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveA{time}")));
      await session.step(90, "And BDDMove{time} link in space gallery should be visible", () => shouldBe(page, el(session.text("BDDMove{time} link in space gallery")), "visible"));
    });
    await run.scenario("The project opens from the first space", async () => {
      await session.step(93, "When user double-clicks on BDDMove{time} link in space gallery", () => doubleClickOn(page, el(session.text("BDDMove{time} link in space gallery"))));
      await session.step(94, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(95, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(96, "And no error balloon should have been shown", () => noErrorBalloon(page));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
      await session.step(98, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(99, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    });
    await run.scenario("The project is moved on to the second space from the first", async () => {
      await session.step(102, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(103, "When user double-clicks on BDDMoveA{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveA{time} tree node inside browse tree"))));
      await session.step(104, "Then the \"BDDMoveA{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveA{time}")));
      await session.step(105, "When user picks \"Move to Space...\" from the context menu of BDDMove{time} link in space gallery", () => pickFromContextMenu(page, "Move to Space...", el(session.text("BDDMove{time} link in space gallery"))));
      await session.step(106, "Then Move to space dialog should be visible", () => shouldBe(page, el("Move to space dialog"), "visible"));
      await session.step(107, "When user selects \"BDDMoveB{time}\" in Space input in Move to space dialog", () => selectIn(page, session.text("BDDMoveB{time}"), el("Space input in Move to space dialog")));
      await session.step(108, "And user clicks on OK button in Move to space dialog", () => clickOn(page, el("OK button in Move to space dialog")));
      await session.step(109, "Then the \"Move to space\" dialog should close", () => dialogCloses(page, "Move to space"));
      await session.step(110, "And an info balloon containing \"Moved BDDMove{time} to BDDMoveB{time}\" should have been shown", () => infoBalloonText(page, session.text("Moved BDDMove{time} to BDDMoveB{time}")));
      await session.step(111, "And the \"BDDMoveB{time}\" space should hold the project \"BDDMove{time}\" on the server", () => spaceHoldsEntity(page, session.text("BDDMoveB{time}"), "project", session.text("BDDMove{time}")));
      await session.step(112, "And the \"BDDMoveA{time}\" space should not hold the project \"BDDMove{time}\" on the server", () => spaceHoldsNotEntity(page, session.text("BDDMoveA{time}"), "project", session.text("BDDMove{time}")));
      await session.step(113, "When user double-clicks on BDDMoveB{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveB{time} tree node inside browse tree"))));
      await session.step(114, "Then the \"BDDMoveB{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveB{time}")));
      await session.step(115, "And BDDMove{time} link in space gallery should be visible", () => shouldBe(page, el(session.text("BDDMove{time} link in space gallery")), "visible"));
      await session.step(116, "When user double-clicks on BDDMoveA{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveA{time} tree node inside browse tree"))));
      await session.step(117, "Then the \"BDDMoveA{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveA{time}")));
      await session.step(118, "And BDDMove{time} link in space gallery should be absent", () => shouldBe(page, el(session.text("BDDMove{time} link in space gallery")), "absent"));
    });
    await run.scenario("The project opens from the second space", async () => {
      await session.step(121, "When user double-clicks on BDDMoveB{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveB{time} tree node inside browse tree"))));
      await session.step(122, "Then the \"BDDMoveB{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveB{time}")));
      await session.step(123, "When user double-clicks on BDDMove{time} link in space gallery", () => doubleClickOn(page, el(session.text("BDDMove{time} link in space gallery"))));
      await session.step(124, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(125, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(126, "And no error balloon should have been shown", () => noErrorBalloon(page));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
      await session.step(128, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(129, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    });
    await run.scenario("The project and both spaces are deleted", async () => {
      await session.step(132, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(133, "When user double-clicks on BDDMoveB{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDMoveB{time} tree node inside browse tree"))));
      await session.step(134, "Then the \"BDDMoveB{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDMoveB{time}")));
      await session.step(135, "When user picks \"Delete Project\" from the context menu of BDDMove{time} link in space gallery", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDMove{time} link in space gallery"))));
      await session.step(136, "Then \"Are you sure?\" dialog should contain text 'Delete project \"BDDMove{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDMove{time}\"?")));
      await session.step(137, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(138, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(139, "And 0 projects named \"BDDMove{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDMove{time}")));
      await session.step(140, "When user picks \"Delete Space\" from the context menu of BDDMoveA{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el(session.text("BDDMoveA{time} tree node inside browse tree"))));
      await session.step(141, "Then \"Are you sure?\" dialog should contain text 'Delete space \"BDDMoveA{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete space \"BDDMoveA{time}\"?")));
      await session.step(142, "And \"Are you sure?\" dialog should contain text \"This will delete space and its related data\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "This will delete space and its related data"));
      await session.step(143, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(144, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(145, "And 0 spaces named \"BDDMoveA{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDMoveA{time}")));
      await session.step(146, "When user picks \"Delete Space\" from the context menu of BDDMoveB{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el(session.text("BDDMoveB{time} tree node inside browse tree"))));
      await session.step(147, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(148, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(149, "And 0 spaces named \"BDDMoveB{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDMoveB{time}")));
    });
    run.finish();
  });
});
