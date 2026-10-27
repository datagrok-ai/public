/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/space-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space, core#976]
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
import {clearField, clickOn, doubleClickOn, dragTo, enterInto, isExpanded, shouldBe, shouldContainText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, childSpaceOnServer, closeAllViews, dialogCloses, noProjectOnServer, noSpaceOnServer, openDataset, projectsOnServer, saveAsProject, spaceOnServer, spacesOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Space edge cases: an emptied name, a child space moved to the root, a delete that takes a linked project with it", () => {
  const session = feature(test, "features/spaces/space-edge-cases.feature", import.meta.url);
  test("An emptied name disables OK in the Rename dialog", {tag: ["@serial", "@spaces", "@edge-cases", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no space named \"BDDEdgeKeep{time}, BDDEdgeLink{time}, BDDEdgeOut{time}, BDDEdgeTaken{time}, BDDEdgeParent{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeKeep{time}, BDDEdgeLink{time}, BDDEdgeOut{time}, BDDEdgeTaken{time}, BDDEdgeParent{time}")));
    await session.step(26, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(29, "Given a space named \"BDDEdgeKeep{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeKeep{time}")));
    await session.step(30, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(31, "When user picks \"Rename...\" from the context menu of BDDEdgeKeep{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDEdgeKeep{time} tree node inside browse tree"))));
    await session.step(32, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
    await session.step(33, "And OK button in Rename project dialog should be enabled", () => shouldBe(page, el("OK button in Rename project dialog"), "enabled"));
    await session.step(34, "When user clears Name input in Rename project dialog", () => clearField(page, el("Name input in Rename project dialog")));
    await session.step(35, "Then Name input in Rename project dialog should be invalid", () => shouldBe(page, el("Name input in Rename project dialog"), "invalid"));
    await session.step(36, "And OK button in Rename project dialog should be disabled", () => shouldBe(page, el("OK button in Rename project dialog"), "disabled"));
    await session.step(37, "When user types \"x\" into Name input in Rename project dialog", () => typeInto(page, "x", el("Name input in Rename project dialog")));
    await session.step(38, "Then OK button in Rename project dialog should be enabled", () => shouldBe(page, el("OK button in Rename project dialog"), "enabled"));
    await session.step(39, "When user clicks on CANCEL button in Rename project dialog", () => clickOn(page, el("CANCEL button in Rename project dialog")));
    await session.step(40, "Then Rename project dialog should be hidden", () => shouldBe(page, el("Rename project dialog"), "hidden"));
    await session.step(41, "And 1 space named \"BDDEdgeKeep{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDEdgeKeep{time}")));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
  test("A child space dragged onto the Spaces node becomes a root space", {tag: ["@serial", "@spaces", "@edge-cases", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no space named \"BDDEdgeKeep{time}, BDDEdgeLink{time}, BDDEdgeOut{time}, BDDEdgeTaken{time}, BDDEdgeParent{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeKeep{time}, BDDEdgeLink{time}, BDDEdgeOut{time}, BDDEdgeTaken{time}, BDDEdgeParent{time}")));
    await session.step(26, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(45, "Given a space named \"BDDEdgeParent{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeParent{time}")));
    await session.step(46, "And a space named \"BDDEdgeOut{time}\" under \"BDDEdgeParent{time}\" is on the server", () => childSpaceOnServer(page, session.text("BDDEdgeOut{time}"), session.text("BDDEdgeParent{time}")));
    await session.step(47, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(48, "And Spaces---BDDEdgeParent{time} tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("Spaces---BDDEdgeParent{time} tree node inside browse tree"))));
    await session.step(49, "When user drags Spaces---BDDEdgeParent{time}---BDDEdgeOut{time} tree node inside browse tree to Spaces tree node inside browse tree", () => dragTo(page, el(session.text("Spaces---BDDEdgeParent{time}---BDDEdgeOut{time} tree node inside browse tree")), el("Spaces tree node inside browse tree")));
    await session.step(50, "Then \"Move space\" dialog should be visible", () => shouldBe(page, el("\"Move space\" dialog"), "visible"));
    await session.step(51, "And \"Move space\" dialog should contain text \"Move BDDEdgeOut{time} out?\"", () => shouldContainText(page, el("\"Move space\" dialog"), session.text("Move BDDEdgeOut{time} out?")));
    await session.step(52, "When user clicks on MOVE button in \"Move space\" dialog", () => clickOn(page, el("MOVE button in \"Move space\" dialog")));
    await session.step(53, "Then the \"Move space\" dialog should close", () => dialogCloses(page, "Move space"));
    await session.step(54, "And Spaces---BDDEdgeOut{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Spaces---BDDEdgeOut{time} tree node inside browse tree")), "visible"));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
    await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A child space named like a root space is refused at the root and keeps its parent", {tag: ["@serial", "@spaces", "@edge-cases", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no space named \"BDDEdgeKeep{time}, BDDEdgeLink{time}, BDDEdgeOut{time}, BDDEdgeTaken{time}, BDDEdgeParent{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeKeep{time}, BDDEdgeLink{time}, BDDEdgeOut{time}, BDDEdgeTaken{time}, BDDEdgeParent{time}")));
    await session.step(26, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(59, "Given a space named \"BDDEdgeTaken{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeTaken{time}")));
    await session.step(60, "And a space named \"BDDEdgeParent{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeParent{time}")));
    await session.step(61, "And a space named \"BDDEdgeTaken{time}\" under \"BDDEdgeParent{time}\" is on the server", () => childSpaceOnServer(page, session.text("BDDEdgeTaken{time}"), session.text("BDDEdgeParent{time}")));
    await session.step(62, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(63, "And Spaces---BDDEdgeParent{time} tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("Spaces---BDDEdgeParent{time} tree node inside browse tree"))));
    await session.step(64, "Then Spaces---BDDEdgeTaken{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Spaces---BDDEdgeTaken{time} tree node inside browse tree")), "visible"));
    await session.step(65, "And Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree")), "visible"));
    await session.step(66, "When user drags Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree to Spaces tree node inside browse tree", () => dragTo(page, el(session.text("Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree")), el("Spaces tree node inside browse tree")));
    await session.step(67, "Then \"Move space\" dialog should be visible", () => shouldBe(page, el("\"Move space\" dialog"), "visible"));
    await session.step(68, "When user clicks on MOVE button in \"Move space\" dialog", () => clickOn(page, el("MOVE button in \"Move space\" dialog")));
    await session.step(69, "Then an error balloon containing \"already exists\" should have been shown", () => errorBalloonText(page, "already exists"));
    await session.step(70, "And \"Move space\" dialog should be visible", () => shouldBe(page, el("\"Move space\" dialog"), "visible"));
    await session.step(71, "When user clicks on CANCEL button in \"Move space\" dialog", () => clickOn(page, el("CANCEL button in \"Move space\" dialog")));
    await session.step(72, "Then the \"Move space\" dialog should close", () => dialogCloses(page, "Move space"));
    await session.step(73, "And Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Spaces---BDDEdgeParent{time}---BDDEdgeTaken{time} tree node inside browse tree")), "visible"));
    await session.step(74, "And Spaces---BDDEdgeTaken{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Spaces---BDDEdgeTaken{time} tree node inside browse tree")), "visible"));
  });
  test("Deleting a space keeps a project that was only linked into it", {tag: ["@serial", "@spaces", "@edge-cases", "@realizes:views.space", "@known-failure:104", "@realizes:core#976"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no space named \"BDDEdgeKeep{time}, BDDEdgeLink{time}, BDDEdgeOut{time}, BDDEdgeTaken{time}, BDDEdgeParent{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeKeep{time}, BDDEdgeLink{time}, BDDEdgeOut{time}, BDDEdgeTaken{time}, BDDEdgeParent{time}")));
    await session.step(26, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await knownFailure(async () => {
      await session.step(80, "Given no project named \"bddedgelinkproj{time}\" is on the server", () => noProjectOnServer(page, session.text("bddedgelinkproj{time}")));
      await session.step(81, "And a space named \"BDDEdgeLink{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeLink{time}")));
      await session.step(82, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(83, "When user saves the current view as project \"bddedgelinkproj{time}\"", () => saveAsProject(page, session.text("bddedgelinkproj{time}")));
      await session.step(84, "And user closes all views", () => closeAllViews(page));
      await session.step(85, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(86, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(87, "And user enters \"bddedgelinkproj{time}\" into gallery search", () => enterInto(page, session.text("bddedgelinkproj{time}"), el("gallery search")));
      await session.step(88, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(89, "Then bddedgelinkproj{time} gallery card should be visible", () => shouldBe(page, el(session.text("bddedgelinkproj{time} gallery card")), "visible"));
      await session.step(90, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(91, "When user drags bddedgelinkproj{time} gallery card to BDDEdgeLink{time} tree node inside browse tree", () => dragTo(page, el(session.text("bddedgelinkproj{time} gallery card")), el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(92, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(93, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(94, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(95, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(96, "When user double-clicks on BDDEdgeLink{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(97, "Then the \"BDDEdgeLink{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeLink{time}")));
      await session.step(98, "And \"bddedgelinkproj{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"bddedgelinkproj{time}\" link in gallery")), "visible"));
      await session.step(99, "And 1 project named \"bddedgelinkproj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("bddedgelinkproj{time}")));
      await session.step(100, "When user picks \"Delete Space\" from the context menu of BDDEdgeLink{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(101, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(102, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(103, "And 0 spaces named \"BDDEdgeLink{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDEdgeLink{time}")));
      await session.step(104, "And 1 project named \"bddedgelinkproj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("bddedgelinkproj{time}")));
    }, 104);
  });
});
