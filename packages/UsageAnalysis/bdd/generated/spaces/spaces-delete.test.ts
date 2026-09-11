/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-delete.feature
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
import {clickOn, enterInto, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Deleting a space", () => {
  const session = feature(test, "features/spaces/spaces-delete.feature", import.meta.url);
  test("Deleting a space", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(12, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(13, "And no space named \"BDD-Del, BDD-Del-Parent, BDD-Del-Child1, BDD-Del-Child2\" is on the server", () => noSpaceOnServer(page, "BDD-Del, BDD-Del-Parent, BDD-Del-Child1, BDD-Del-Child2"));
    await run.scenario("Deleting asks first", async () => {
      await session.step(16, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(17, "And user enters \"BDD-Del\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Del", el("Name input in Create Space dialog")));
      await session.step(18, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(19, "Then 1 space named \"BDD-Del\" should be on the server", () => spacesOnServer(page, 1, "BDD-Del"));
      await session.step(20, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(21, "When user picks \"Delete Space\" from the context menu of BDD-Del tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("BDD-Del tree node inside browse tree")));
      await session.step(22, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(23, "And \"Are you sure?\" dialog should contain text \"BDD-Del\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "BDD-Del"));
    });
    await run.scenario("Cancelling the confirmation keeps the space", async () => {
      await session.step(26, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
      await session.step(27, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(28, "And 1 space named \"BDD-Del\" should be on the server", () => spacesOnServer(page, 1, "BDD-Del"));
      await session.step(29, "And BDD-Del tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Del tree node inside browse tree"), "visible"));
    });
    await run.scenario("Confirming removes it from the server and the tree", async () => {
      await session.step(32, "When user picks \"Delete Space\" from the context menu of BDD-Del tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("BDD-Del tree node inside browse tree")));
      await session.step(33, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(34, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(35, "And 0 spaces named \"BDD-Del\" should be on the server", () => spacesOnServer(page, 0, "BDD-Del"));
      await session.step(36, "And BDD-Del tree node inside browse tree should be absent", () => shouldBe(page, el("BDD-Del tree node inside browse tree"), "absent"));
    });
    await run.scenario("Deleting one child leaves its sibling", async () => {
      await session.step(39, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(40, "And user enters \"BDD-Del-Parent\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Del-Parent", el("Name input in Create Space dialog")));
      await session.step(41, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(42, "Then 1 space named \"BDD-Del-Parent\" should be on the server", () => spacesOnServer(page, 1, "BDD-Del-Parent"));
      await session.step(43, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(44, "When user picks \"Create Child Space...\" from the context menu of BDD-Del-Parent tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Del-Parent tree node inside browse tree")));
      await session.step(45, "And user enters \"BDD-Del-Child1\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Del-Child1", el("Name input in Create Space dialog")));
      await session.step(46, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(47, "Then BDD-Del-Child1 tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Del-Child1 tree node inside browse tree"), "visible"));
      await session.step(48, "When user picks \"Create Child Space...\" from the context menu of BDD-Del-Parent tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Del-Parent tree node inside browse tree")));
      await session.step(49, "And user enters \"BDD-Del-Child2\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Del-Child2", el("Name input in Create Space dialog")));
      await session.step(50, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(51, "Then BDD-Del-Child2 tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Del-Child2 tree node inside browse tree"), "visible"));
      await session.step(52, "When user picks \"Delete Space\" from the context menu of BDD-Del-Child1 tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("BDD-Del-Child1 tree node inside browse tree")));
      await session.step(53, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(54, "Then BDD-Del-Child1 tree node inside browse tree should be absent", () => shouldBe(page, el("BDD-Del-Child1 tree node inside browse tree"), "absent"));
      await session.step(55, "And BDD-Del-Child2 tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Del-Child2 tree node inside browse tree"), "visible"));
    });
    await run.scenario("Deleting the parent takes the remaining child with it", async () => {
      await session.step(58, "When user picks \"Delete Space\" from the context menu of BDD-Del-Parent tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("BDD-Del-Parent tree node inside browse tree")));
      await session.step(59, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(60, "Then 0 spaces named \"BDD-Del-Parent\" should be on the server", () => spacesOnServer(page, 0, "BDD-Del-Parent"));
      await session.step(61, "And BDD-Del-Parent tree node inside browse tree should be absent", () => shouldBe(page, el("BDD-Del-Parent tree node inside browse tree"), "absent"));
      await session.step(62, "And BDD-Del-Child2 tree node inside browse tree should be absent", () => shouldBe(page, el("BDD-Del-Child2 tree node inside browse tree"), "absent"));
    });
    run.finish();
  });
});
