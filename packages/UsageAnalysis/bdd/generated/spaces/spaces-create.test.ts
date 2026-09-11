/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-create.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space]
--- */
import {test} from '@playwright/test';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noSpaceOnServer, spacesOnServer} from '../../bindings/spaces.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, errorBalloonText, menuDoesNotList, menuLists, openContextMenu, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Creating a space", () => {
  const session = feature(test, "features/spaces/spaces-create.feature", import.meta.url);
  test("Creating a space", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And no space named \"BDD-Root, BDD-Dup, BDD-Parent, BDD-Child, BDD Name With Spaces\" is on the server", () => noSpaceOnServer(page, "BDD-Root, BDD-Dup, BDD-Parent, BDD-Child, BDD Name With Spaces"));
    await run.scenario("A root space is created from the Spaces node", async () => {
      await session.step(19, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(20, "Then Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(21, "When user enters \"BDD-Root\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Root", el("Name input in Create Space dialog")));
      await session.step(22, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(23, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(24, "And 1 space named \"BDD-Root\" should be on the server", () => spacesOnServer(page, 1, "BDD-Root"));
      await session.step(25, "And BDD-Root tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Root tree node inside browse tree"), "visible"));
    });
    await run.scenario("The space offers its actions", async () => {
      await session.step(28, "When user opens the context menu of BDD-Root tree node inside browse tree", () => openContextMenu(page, el("BDD-Root tree node inside browse tree")));
      await session.step(29, "Then the open menu should list \"Share...\"", () => menuLists(page, "Share..."));
      await session.step(30, "And the open menu should list \"Rename...\"", () => menuLists(page, "Rename..."));
      await session.step(31, "And the open menu should list \"Delete Space\"", () => menuLists(page, "Delete Space"));
      await session.step(32, "And the open menu should list \"Create Child Space...\"", () => menuLists(page, "Create Child Space..."));
      await session.step(33, "And the open menu should list \"Add to favorites\"", () => menuLists(page, "Add to favorites"));
      await session.step(34, "And the open menu should not list \"Duplicate\"", () => menuDoesNotList(page, "Duplicate"));
      await session.step(35, "When user closes the context menu", () => closeContextMenu(page));
    });
    await run.scenario("An empty name disables OK, and typing one enables it again", async () => {
      await session.step(38, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(39, "And user clears Name input in Create Space dialog", () => clearField(page, el("Name input in Create Space dialog")));
      await session.step(40, "Then OK button in Create Space dialog should be disabled", () => shouldBe(page, el("OK button in Create Space dialog"), "disabled"));
      await session.step(41, "When user enters \"BDD-Dup\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Dup", el("Name input in Create Space dialog")));
      await session.step(42, "Then OK button in Create Space dialog should be enabled", () => shouldBe(page, el("OK button in Create Space dialog"), "enabled"));
      await session.step(43, "When user clicks on CANCEL button in Create Space dialog", () => clickOn(page, el("CANCEL button in Create Space dialog")));
      await session.step(44, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(45, "And 0 spaces named \"BDD-Dup\" should be on the server", () => spacesOnServer(page, 0, "BDD-Dup"));
    });
    await run.scenario("A second root space of the same name is refused", async () => {
      await session.step(48, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(49, "And user enters \"BDD-Dup\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Dup", el("Name input in Create Space dialog")));
      await session.step(50, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(51, "Then 1 space named \"BDD-Dup\" should be on the server", () => spacesOnServer(page, 1, "BDD-Dup"));
      await session.step(52, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(53, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(54, "And user enters \"BDD-Dup\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Dup", el("Name input in Create Space dialog")));
      await session.step(55, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(56, "Then an error balloon containing \"Root project with same name already exists\" should have been shown", () => errorBalloonText(page, "Root project with same name already exists"));
      await session.step(57, "And 1 space named \"BDD-Dup\" should be on the server", () => spacesOnServer(page, 1, "BDD-Dup"));
      await session.step(58, "And Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(59, "When user clicks on CANCEL button in Create Space dialog", () => clickOn(page, el("CANCEL button in Create Space dialog")));
      await session.step(60, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
    });
    await run.scenario("A child space is created under a root space", async () => {
      await session.step(63, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(64, "And user enters \"BDD-Parent\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Parent", el("Name input in Create Space dialog")));
      await session.step(65, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(66, "Then 1 space named \"BDD-Parent\" should be on the server", () => spacesOnServer(page, 1, "BDD-Parent"));
      await session.step(67, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(68, "When user picks \"Create Child Space...\" from the context menu of BDD-Parent tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Parent tree node inside browse tree")));
      await session.step(69, "And user enters \"BDD-Child\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Child", el("Name input in Create Space dialog")));
      await session.step(70, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(71, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(72, "And BDD-Child tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Child tree node inside browse tree"), "visible"));
    });
    await run.scenario("A second child of the same name is refused", async () => {
      await session.step(75, "When user picks \"Create Child Space...\" from the context menu of BDD-Parent tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Parent tree node inside browse tree")));
      await session.step(76, "And user enters \"BDD-Child\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Child", el("Name input in Create Space dialog")));
      await session.step(77, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(78, "Then an error balloon containing \"already exists\" should have been shown", () => errorBalloonText(page, "already exists"));
      await session.step(79, "And Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(80, "When user clicks on CANCEL button in Create Space dialog", () => clickOn(page, el("CANCEL button in Create Space dialog")));
      await session.step(81, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
    });
    await run.scenario("A name with spaces is kept as typed", async () => {
      await session.step(84, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(85, "And user enters \"BDD Name With Spaces\" into Name input in Create Space dialog", () => enterInto(page, "BDD Name With Spaces", el("Name input in Create Space dialog")));
      await session.step(86, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(87, "Then 1 space named \"BDD Name With Spaces\" should be on the server", () => spacesOnServer(page, 1, "BDD Name With Spaces"));
      await session.step(88, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(89, "And BDD Name With Spaces tree node inside browse tree should be visible", () => shouldBe(page, el("BDD Name With Spaces tree node inside browse tree"), "visible"));
    });
    run.finish();
  });
});
