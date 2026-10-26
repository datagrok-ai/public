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
import {clickOn, doubleClickOn, dragTo, enterInto, isExpanded, pressKey, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, dialogCloses, noProjectOnServer, noSpaceOnServer, openDataset, projectsOnServer, saveAsProject, spaceOnServer, spacesOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Space edge cases: Escape on its dialogs, and a delete that takes a linked project with it", () => {
  const session = feature(test, "features/spaces/space-edge-cases.feature", import.meta.url);
  test("Escape on the Create Space dialog creates nothing", {tag: ["@serial", "@spaces", "@edge-cases", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no space named \"BDDEdgeEsc{time}, BDDEdgeMade{time}, BDDEdgeKeep{time}, BDDEdgeKeepNew{time}, BDDEdgeLink{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeEsc{time}, BDDEdgeMade{time}, BDDEdgeKeep{time}, BDDEdgeKeepNew{time}, BDDEdgeLink{time}")));
    await session.step(26, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(29, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
    await session.step(30, "Then Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
    await session.step(31, "When user enters \"BDDEdgeEsc{time}\" into Name input in Create Space dialog", () => enterInto(page, session.text("BDDEdgeEsc{time}"), el("Name input in Create Space dialog")));
    await session.step(32, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(33, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
    await session.step(34, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
    await session.step(35, "And user enters \"BDDEdgeMade{time}\" into Name input in Create Space dialog", () => enterInto(page, session.text("BDDEdgeMade{time}"), el("Name input in Create Space dialog")));
    await session.step(36, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
    await session.step(37, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
    await session.step(38, "And 1 space named \"BDDEdgeMade{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDEdgeMade{time}")));
    await session.step(39, "And 0 spaces named \"BDDEdgeEsc{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDEdgeEsc{time}")));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
  });
  test("Escape on the Rename dialog keeps the old name", {tag: ["@serial", "@spaces", "@edge-cases", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no space named \"BDDEdgeEsc{time}, BDDEdgeMade{time}, BDDEdgeKeep{time}, BDDEdgeKeepNew{time}, BDDEdgeLink{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeEsc{time}, BDDEdgeMade{time}, BDDEdgeKeep{time}, BDDEdgeKeepNew{time}, BDDEdgeLink{time}")));
    await session.step(26, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(43, "Given a space named \"BDDEdgeKeep{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeKeep{time}")));
    await session.step(44, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(45, "When user picks \"Rename...\" from the context menu of BDDEdgeKeep{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDEdgeKeep{time} tree node inside browse tree"))));
    await session.step(46, "Then Name input in Rename project dialog should have value \"BDDEdgeKeep{time}\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), session.text("BDDEdgeKeep{time}")));
    await session.step(47, "When user enters \"BDDEdgeKeepNew{time}\" into Name input in Rename project dialog", () => enterInto(page, session.text("BDDEdgeKeepNew{time}"), el("Name input in Rename project dialog")));
    await session.step(48, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(49, "Then Rename project dialog should be hidden", () => shouldBe(page, el("Rename project dialog"), "hidden"));
    await session.step(50, "And 1 space named \"BDDEdgeKeep{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDEdgeKeep{time}")));
    await session.step(51, "And 0 spaces named \"BDDEdgeKeepNew{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDEdgeKeepNew{time}")));
    await session.step(52, "And BDDEdgeKeep{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDEdgeKeep{time} tree node inside browse tree")), "visible"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
  test("Escape on the delete confirmation keeps the space", {tag: ["@serial", "@spaces", "@edge-cases", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no space named \"BDDEdgeEsc{time}, BDDEdgeMade{time}, BDDEdgeKeep{time}, BDDEdgeKeepNew{time}, BDDEdgeLink{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeEsc{time}, BDDEdgeMade{time}, BDDEdgeKeep{time}, BDDEdgeKeepNew{time}, BDDEdgeLink{time}")));
    await session.step(26, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(56, "Given a space named \"BDDEdgeKeep{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeKeep{time}")));
    await session.step(57, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(58, "When user picks \"Delete Space\" from the context menu of BDDEdgeKeep{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el(session.text("BDDEdgeKeep{time} tree node inside browse tree"))));
    await session.step(59, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
    await session.step(60, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(61, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
    await session.step(62, "And 1 space named \"BDDEdgeKeep{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDEdgeKeep{time}")));
    await session.step(63, "And BDDEdgeKeep{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDEdgeKeep{time} tree node inside browse tree")), "visible"));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
  test("Deleting a space keeps a project that was only linked into it", {tag: ["@serial", "@spaces", "@edge-cases", "@realizes:views.space", "@known-failure:94", "@realizes:core#976"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no space named \"BDDEdgeEsc{time}, BDDEdgeMade{time}, BDDEdgeKeep{time}, BDDEdgeKeepNew{time}, BDDEdgeLink{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeEsc{time}, BDDEdgeMade{time}, BDDEdgeKeep{time}, BDDEdgeKeepNew{time}, BDDEdgeLink{time}")));
    await session.step(26, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await knownFailure(async () => {
      await session.step(70, "Given no project named \"bddedgelinkproj{time}\" is on the server", () => noProjectOnServer(page, session.text("bddedgelinkproj{time}")));
      await session.step(71, "And a space named \"BDDEdgeLink{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeLink{time}")));
      await session.step(72, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(73, "When user saves the current view as project \"bddedgelinkproj{time}\"", () => saveAsProject(page, session.text("bddedgelinkproj{time}")));
      await session.step(74, "And user closes all views", () => closeAllViews(page));
      await session.step(75, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(76, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(77, "And user enters \"bddedgelinkproj{time}\" into gallery search", () => enterInto(page, session.text("bddedgelinkproj{time}"), el("gallery search")));
      await session.step(78, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(79, "Then bddedgelinkproj{time} gallery card should be visible", () => shouldBe(page, el(session.text("bddedgelinkproj{time} gallery card")), "visible"));
      await session.step(80, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(81, "When user drags bddedgelinkproj{time} gallery card to BDDEdgeLink{time} tree node inside browse tree", () => dragTo(page, el(session.text("bddedgelinkproj{time} gallery card")), el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(82, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(83, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(84, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(85, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(86, "When user double-clicks on BDDEdgeLink{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(87, "Then the \"BDDEdgeLink{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeLink{time}")));
      await session.step(88, "And \"bddedgelinkproj{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"bddedgelinkproj{time}\" link in gallery")), "visible"));
      await session.step(89, "And 1 project named \"bddedgelinkproj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("bddedgelinkproj{time}")));
      await session.step(90, "When user picks \"Delete Space\" from the context menu of BDDEdgeLink{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(91, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(92, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(93, "And 0 spaces named \"BDDEdgeLink{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDEdgeLink{time}")));
      await session.step(94, "And 1 project named \"bddedgelinkproj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("bddedgelinkproj{time}")));
    }, 94);
  });
});
