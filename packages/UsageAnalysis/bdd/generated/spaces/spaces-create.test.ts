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
import {clearField, clickOn, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses, noSpaceOnServer, spacesOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, errorBalloonText, menuDoesNotList, menuLists, openContextMenu, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Creating a space", () => {
  const session = feature(test, "features/spaces/spaces-create.feature", import.meta.url);
  test("Creating a space", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And no space named \"BDD-Root, BDD-Dup, BDD-Parent, BDD-Child\" is on the server", () => noSpaceOnServer(page, "BDD-Root, BDD-Dup, BDD-Parent, BDD-Child"));
    await run.scenario("A root space is created from the Spaces node", async () => {
      await session.step(22, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(23, "Then Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(24, "When user enters \"BDD-Root\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Root", el("Name input in Create Space dialog")));
      await session.step(25, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(26, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(27, "And 1 space named \"BDD-Root\" should be on the server", () => spacesOnServer(page, 1, "BDD-Root"));
      await session.step(28, "And Spaces---BDD-Root tree node inside browse tree should be visible", () => shouldBe(page, el("Spaces---BDD-Root tree node inside browse tree"), "visible"));
    });
    await run.scenario("The space offers its actions", async () => {
      await session.step(31, "When user opens the context menu of Spaces---BDD-Root tree node inside browse tree", () => openContextMenu(page, el("Spaces---BDD-Root tree node inside browse tree")));
      await session.step(32, "Then the open menu should list \"Share...\"", () => menuLists(page, "Share..."));
      await session.step(33, "And the open menu should list \"Rename...\"", () => menuLists(page, "Rename..."));
      await session.step(34, "And the open menu should list \"Delete Space\"", () => menuLists(page, "Delete Space"));
      await session.step(35, "And the open menu should list \"Create Child Space...\"", () => menuLists(page, "Create Child Space..."));
      await session.step(36, "And the open menu should list \"Add To Favorites\"", () => menuLists(page, "Add To Favorites"));
      await session.step(37, "And the open menu should not list \"Duplicate\"", () => menuDoesNotList(page, "Duplicate"));
      await session.step(38, "When user closes the context menu", () => closeContextMenu(page));
    });
    await run.scenario("An empty name disables OK, and typing one enables it again", async () => {
      await session.step(41, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(42, "And user clears Name input in Create Space dialog", () => clearField(page, el("Name input in Create Space dialog")));
      await session.step(43, "Then OK button in Create Space dialog should be disabled", () => shouldBe(page, el("OK button in Create Space dialog"), "disabled"));
      await session.step(44, "When user enters \"BDD-Dup\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Dup", el("Name input in Create Space dialog")));
      await session.step(45, "Then OK button in Create Space dialog should be enabled", () => shouldBe(page, el("OK button in Create Space dialog"), "enabled"));
      await session.step(46, "When user clicks on CANCEL button in Create Space dialog", () => clickOn(page, el("CANCEL button in Create Space dialog")));
      await session.step(47, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(48, "And 0 spaces named \"BDD-Dup\" should be on the server", () => spacesOnServer(page, 0, "BDD-Dup"));
    });
    await run.scenario("A second root space of the same name is refused", async () => {
      await session.step(51, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(52, "And user enters \"BDD-Dup\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Dup", el("Name input in Create Space dialog")));
      await session.step(53, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(54, "Then 1 space named \"BDD-Dup\" should be on the server", () => spacesOnServer(page, 1, "BDD-Dup"));
      await session.step(55, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(56, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(57, "And user enters \"BDD-Dup\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Dup", el("Name input in Create Space dialog")));
      await session.step(58, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(59, "Then an error balloon containing \"Root project with same name already exists\" should have been shown", () => errorBalloonText(page, "Root project with same name already exists"));
      await session.step(60, "And 1 space named \"BDD-Dup\" should be on the server", () => spacesOnServer(page, 1, "BDD-Dup"));
      await session.step(61, "And Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(62, "When user clicks on CANCEL button in Create Space dialog", () => clickOn(page, el("CANCEL button in Create Space dialog")));
      await session.step(63, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
    });
    await run.scenario("A child space is created under a root space", async () => {
      await session.step(66, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(67, "And user enters \"BDD-Parent\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Parent", el("Name input in Create Space dialog")));
      await session.step(68, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(69, "Then 1 space named \"BDD-Parent\" should be on the server", () => spacesOnServer(page, 1, "BDD-Parent"));
      await session.step(70, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(71, "When user picks \"Create Child Space...\" from the context menu of Spaces---BDD-Parent tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("Spaces---BDD-Parent tree node inside browse tree")));
      await session.step(72, "And user enters \"BDD-Child\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Child", el("Name input in Create Space dialog")));
      await session.step(73, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(74, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(75, "And BDD-Child tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Child tree node inside browse tree"), "visible"));
    });
    run.finish();
  });
});
