/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-create.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space]
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
import {clearField, clickOn, dragTo, enterInto, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, childSpaceOnServer, dialogCloses, noSpaceOnServer, spaceOnServer, spacesOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, errorBalloonText, menuDoesNotList, menuLists, noBalloons, noErrors, openContextMenu, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Creating a space", () => {
  const session = feature(test, "features/spaces/spaces-create.feature", import.meta.url);
  test("Creating a space", {tag: ["@journey", "@serial", "@spaces", "@realizes:views.space", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And no space named \"BDD-Root, BDD-Dup, BDD-Parent, BDD-Child\" is on the server", () => noSpaceOnServer(page, "BDD-Root, BDD-Dup, BDD-Parent, BDD-Child"));
    await run.scenario("A root space is created from the Spaces node", async () => {
      await session.step(32, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(33, "Then Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(34, "When user enters \"BDD-Root\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Root", el("Name input in Create Space dialog")));
      await session.step(35, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(36, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(37, "And 1 space named \"BDD-Root\" should be on the server", () => spacesOnServer(page, 1, "BDD-Root"));
      await session.step(38, "And Spaces---BDD-Root tree node inside browse tree should be visible", () => shouldBe(page, el("Spaces---BDD-Root tree node inside browse tree"), "visible"));
    });
    await run.scenario("The space offers its actions", async () => {
      await session.step(41, "When user opens the context menu of Spaces---BDD-Root tree node inside browse tree", () => openContextMenu(page, el("Spaces---BDD-Root tree node inside browse tree")));
      await session.step(42, "Then the open menu should list \"Share...\"", () => menuLists(page, "Share..."));
      await session.step(43, "And the open menu should list \"Rename...\"", () => menuLists(page, "Rename..."));
      await session.step(44, "And the open menu should list \"Delete Space\"", () => menuLists(page, "Delete Space"));
      await session.step(45, "And the open menu should list \"Create Child Space...\"", () => menuLists(page, "Create Child Space..."));
      await session.step(46, "And the open menu should list \"Add To Favorites\"", () => menuLists(page, "Add To Favorites"));
      await session.step(47, "And the open menu should not list \"Duplicate\"", () => menuDoesNotList(page, "Duplicate"));
      await session.step(48, "When user closes the context menu", () => closeContextMenu(page));
    });
    await run.scenario("An empty name disables OK, and typing one enables it again", async () => {
      await session.step(51, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(52, "And user clears Name input in Create Space dialog", () => clearField(page, el("Name input in Create Space dialog")));
      await session.step(53, "Then OK button in Create Space dialog should be disabled", () => shouldBe(page, el("OK button in Create Space dialog"), "disabled"));
      await session.step(54, "When user enters \"BDD-Dup\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Dup", el("Name input in Create Space dialog")));
      await session.step(55, "Then OK button in Create Space dialog should be enabled", () => shouldBe(page, el("OK button in Create Space dialog"), "enabled"));
      await session.step(56, "When user clicks on CANCEL button in Create Space dialog", () => clickOn(page, el("CANCEL button in Create Space dialog")));
      await session.step(57, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(58, "And 0 spaces named \"BDD-Dup\" should be on the server", () => spacesOnServer(page, 0, "BDD-Dup"));
    });
    await run.scenario("A second root space of the same name is refused", async () => {
      await session.step(61, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(62, "And user enters \"BDD-Dup\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Dup", el("Name input in Create Space dialog")));
      await session.step(63, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(64, "Then 1 space named \"BDD-Dup\" should be on the server", () => spacesOnServer(page, 1, "BDD-Dup"));
      await session.step(65, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(66, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(67, "And user enters \"BDD-Dup\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Dup", el("Name input in Create Space dialog")));
      await session.step(68, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(69, "Then an error balloon containing \"Root project with same name already exists\" should have been shown", () => errorBalloonText(page, "Root project with same name already exists"));
      await session.step(70, "And 1 space named \"BDD-Dup\" should be on the server", () => spacesOnServer(page, 1, "BDD-Dup"));
      await session.step(71, "And Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(72, "When user clicks on CANCEL button in Create Space dialog", () => clickOn(page, el("CANCEL button in Create Space dialog")));
      await session.step(73, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
    });
    await run.scenario("A child space is created under a root space", async () => {
      await session.step(76, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(77, "And user enters \"BDD-Parent\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Parent", el("Name input in Create Space dialog")));
      await session.step(78, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(79, "Then 1 space named \"BDD-Parent\" should be on the server", () => spacesOnServer(page, 1, "BDD-Parent"));
      await session.step(80, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(81, "When user picks \"Create Child Space...\" from the context menu of Spaces---BDD-Parent tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("Spaces---BDD-Parent tree node inside browse tree")));
      await session.step(82, "And user enters \"BDD-Child\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Child", el("Name input in Create Space dialog")));
      await session.step(83, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(84, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(85, "And BDD-Child tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Child tree node inside browse tree"), "visible"));
    });
    await run.scenario("A child space dropped on the Spaces node is moved out to the root", async () => {
      await session.step(89, "Given no space named \"BDDEdgeOut{time}, BDDEdgeOutParent{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeOut{time}, BDDEdgeOutParent{time}")));
      await session.step(90, "And a space named \"BDDEdgeOutParent{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeOutParent{time}")));
      await session.step(91, "And a space named \"BDDEdgeOut{time}\" under \"BDDEdgeOutParent{time}\" is on the server", () => childSpaceOnServer(page, session.text("BDDEdgeOut{time}"), session.text("BDDEdgeOutParent{time}")));
      await session.step(92, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(93, "And Spaces---BDDEdgeOutParent{time} tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("Spaces---BDDEdgeOutParent{time} tree node inside browse tree"))));
      await session.step(94, "Then Spaces---BDDEdgeOutParent{time}---BDDEdgeOut{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Spaces---BDDEdgeOutParent{time}---BDDEdgeOut{time} tree node inside browse tree")), "visible"));
      await session.step(95, "When user drags Spaces---BDDEdgeOutParent{time}---BDDEdgeOut{time} tree node inside browse tree to Spaces tree node inside browse tree", () => dragTo(page, el(session.text("Spaces---BDDEdgeOutParent{time}---BDDEdgeOut{time} tree node inside browse tree")), el("Spaces tree node inside browse tree")));
      await session.step(96, "Then \"Move space\" dialog should be visible", () => shouldBe(page, el("\"Move space\" dialog"), "visible"));
      await session.step(97, "And \"Move space\" dialog should contain text \"Move BDDEdgeOut{time} out?\"", () => shouldContainText(page, el("\"Move space\" dialog"), session.text("Move BDDEdgeOut{time} out?")));
      await session.step(98, "When user clicks on MOVE button in \"Move space\" dialog", () => clickOn(page, el("MOVE button in \"Move space\" dialog")));
      await session.step(99, "Then the \"Move space\" dialog should close", () => dialogCloses(page, "Move space"));
      await session.step(100, "And Spaces---BDDEdgeOut{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Spaces---BDDEdgeOut{time} tree node inside browse tree")), "visible"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
      await session.step(102, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A child named like a root space is refused at the root and keeps its parent", async () => {
      await session.step(106, "Given no space named \"BDDEdgeTaken{time}, BDDEdgeParent{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeTaken{time}, BDDEdgeParent{time}")));
      await session.step(107, "And a space named \"BDDEdgeTaken{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeTaken{time}")));
      await session.step(108, "And a space named \"BDDEdgeParent{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeParent{time}")));
      await session.step(109, "And a space named \"BDDEdgeTaken{time}\" under \"BDDEdgeParent{time}\" is on the server", () => childSpaceOnServer(page, session.text("BDDEdgeTaken{time}"), session.text("BDDEdgeParent{time}")));
      await session.step(110, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(111, "And Spaces---BDDEdgeParent{time} tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("Spaces---BDDEdgeParent{time} tree node inside browse tree"))));
      await session.step(112, "Then Spaces---BDDEdgeTaken{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Spaces---BDDEdgeTaken{time} tree node inside browse tree")), "visible"));
      await session.step(113, "And Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree")), "visible"));
      await session.step(114, "When user drags Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree to Spaces tree node inside browse tree", () => dragTo(page, el(session.text("Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree")), el("Spaces tree node inside browse tree")));
      await session.step(115, "Then \"Move space\" dialog should be visible", () => shouldBe(page, el("\"Move space\" dialog"), "visible"));
      await session.step(116, "When user clicks on MOVE button in \"Move space\" dialog", () => clickOn(page, el("MOVE button in \"Move space\" dialog")));
      await session.step(117, "Then an error balloon containing \"already exists\" should have been shown", () => errorBalloonText(page, "already exists"));
      await session.step(118, "And \"Move space\" dialog should be visible", () => shouldBe(page, el("\"Move space\" dialog"), "visible"));
      await session.step(119, "When user clicks on CANCEL button in \"Move space\" dialog", () => clickOn(page, el("CANCEL button in \"Move space\" dialog")));
      await session.step(120, "Then the \"Move space\" dialog should close", () => dialogCloses(page, "Move space"));
      await session.step(121, "And Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree")), "visible"));
      await session.step(122, "And Spaces---BDDEdgeTaken{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Spaces---BDDEdgeTaken{time} tree node inside browse tree")), "visible"));
    });
    run.finish();
  });
});
