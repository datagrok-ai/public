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
import {noSpaceOnServer, spacesOnServer, treeHidesSpace, treeShowsSpace} from '../../bindings/spaces.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
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
    await session.step(12, "And no space named \"BDD-Del, BDD-Del-Parent, BDD-Del-Child1, BDD-Del-Child2\" is on the server", () => noSpaceOnServer(page, "BDD-Del, BDD-Del-Parent, BDD-Del-Child1, BDD-Del-Child2"));
    await run.scenario("Deleting asks first", async () => {
      await session.step(15, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(16, "And user enters \"BDD-Del\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Del", el("Name input in Create Space dialog")));
      await session.step(17, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(18, "Then 1 space named \"BDD-Del\" should be on the server", () => spacesOnServer(page, 1, "BDD-Del"));
      await session.step(19, "When user picks \"Delete Space\" from the context menu of BDD-Del tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("BDD-Del tree node inside browse tree")));
      await session.step(20, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(21, "And \"Are you sure?\" dialog should contain text \"BDD-Del\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "BDD-Del"));
    });
    await run.scenario("Cancelling the confirmation keeps the space", async () => {
      await session.step(24, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
      await session.step(25, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(26, "And 1 space named \"BDD-Del\" should be on the server", () => spacesOnServer(page, 1, "BDD-Del"));
      await session.step(27, "And the browse tree should show the \"BDD-Del\" space", () => treeShowsSpace(page, "BDD-Del"));
    });
    await run.scenario("Confirming removes it from the server and the tree", async () => {
      await session.step(30, "When user picks \"Delete Space\" from the context menu of BDD-Del tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("BDD-Del tree node inside browse tree")));
      await session.step(31, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(32, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(33, "And 0 spaces named \"BDD-Del\" should be on the server", () => spacesOnServer(page, 0, "BDD-Del"));
      await session.step(34, "And the browse tree should not show the \"BDD-Del\" space", () => treeHidesSpace(page, "BDD-Del"));
    });
    await run.scenario("Deleting one child leaves its sibling", async () => {
      await session.step(37, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(38, "And user enters \"BDD-Del-Parent\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Del-Parent", el("Name input in Create Space dialog")));
      await session.step(39, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(40, "Then 1 space named \"BDD-Del-Parent\" should be on the server", () => spacesOnServer(page, 1, "BDD-Del-Parent"));
      await session.step(41, "When user picks \"Create Child Space...\" from the context menu of BDD-Del-Parent tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Del-Parent tree node inside browse tree")));
      await session.step(42, "And user enters \"BDD-Del-Child1\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Del-Child1", el("Name input in Create Space dialog")));
      await session.step(43, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(44, "Then the browse tree should show the \"BDD-Del-Child1\" space", () => treeShowsSpace(page, "BDD-Del-Child1"));
      await session.step(45, "When user picks \"Create Child Space...\" from the context menu of BDD-Del-Parent tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Del-Parent tree node inside browse tree")));
      await session.step(46, "And user enters \"BDD-Del-Child2\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Del-Child2", el("Name input in Create Space dialog")));
      await session.step(47, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(48, "Then the browse tree should show the \"BDD-Del-Child2\" space", () => treeShowsSpace(page, "BDD-Del-Child2"));
      await session.step(49, "When user picks \"Delete Space\" from the context menu of BDD-Del-Child1 tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("BDD-Del-Child1 tree node inside browse tree")));
      await session.step(50, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(51, "Then the browse tree should not show the \"BDD-Del-Child1\" space", () => treeHidesSpace(page, "BDD-Del-Child1"));
      await session.step(52, "And the browse tree should show the \"BDD-Del-Child2\" space", () => treeShowsSpace(page, "BDD-Del-Child2"));
    });
    await run.scenario("Deleting the parent takes the remaining child with it", async () => {
      await session.step(55, "When user picks \"Delete Space\" from the context menu of BDD-Del-Parent tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("BDD-Del-Parent tree node inside browse tree")));
      await session.step(56, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(57, "Then 0 spaces named \"BDD-Del-Parent\" should be on the server", () => spacesOnServer(page, 0, "BDD-Del-Parent"));
      await session.step(58, "And the browse tree should not show the \"BDD-Del-Parent\" space", () => treeHidesSpace(page, "BDD-Del-Parent"));
      await session.step(59, "And the browse tree should not show the \"BDD-Del-Child2\" space", () => treeHidesSpace(page, "BDD-Del-Child2"));
    });
    run.finish();
  });
});
