/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-share-and-edges.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog]
--- */
import {test} from '@playwright/test';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noSpaceOnServer, sharingPaneLists, sharingPaneListsNot, spacesOnServer} from '../../bindings/spaces.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, dragTo, enterInto, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses, pickSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sharing a space, and what it refuses", () => {
  const session = feature(test, "features/spaces/spaces-share-and-edges.feature", import.meta.url);
  test("Sharing a space, and what it refuses", {tag: ["@journey", "@spaces", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(33, "Given user is logged in", () => loggedIn(page));
    await session.step(34, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(35, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(36, "And no space named \"BDD-Share, BDD-Share-Child\" is on the server", () => noSpaceOnServer(page, "BDD-Share, BDD-Share-Child"));
    await run.scenario("A space and a child to share", async () => {
      await session.step(39, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(40, "And user enters \"BDD-Share\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Share", el("Name input in Create Space dialog")));
      await session.step(41, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(42, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(43, "And 1 space named \"BDD-Share\" should be on the server", () => spacesOnServer(page, 1, "BDD-Share"));
      await session.step(44, "When user picks \"Create Child Space...\" from the context menu of BDD-Share tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Share tree node inside browse tree")));
      await session.step(45, "And user enters \"BDD-Share-Child\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Share-Child", el("Name input in Create Space dialog")));
      await session.step(46, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(47, "When user double-clicks on BDD-Share tree node inside browse tree", () => doubleClickOn(page, el("BDD-Share tree node inside browse tree")));
      await session.step(48, "Then the \"BDD-Share\" view should be current", () => viewIsCurrent(page, "BDD-Share"));
      await session.step(49, "And BDD-Share-Child link in space gallery should be visible", () => shouldBe(page, el("BDD-Share-Child link in space gallery"), "visible"));
    });
    await run.scenario("The Share dialog asks who and how much", async () => {
      await session.step(52, "When user picks \"Share...\" from the context menu of BDD-Share tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("BDD-Share tree node inside browse tree")));
      await session.step(53, "Then \"Share BDD-Share\" dialog should be visible", () => shouldBe(page, el("\"Share BDD-Share\" dialog"), "visible"));
      await session.step(54, "And \"User, group, or email\" input in \"Share BDD-Share\" dialog should be visible", () => shouldBe(page, el("\"User, group, or email\" input in \"Share BDD-Share\" dialog"), "visible"));
      await session.step(55, "And share access selector should be visible", () => shouldBe(page, el("share access selector"), "visible"));
      await session.step(56, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(57, "When user clicks on CANCEL button in \"Share BDD-Share\" dialog", () => clickOn(page, el("CANCEL button in \"Share BDD-Share\" dialog")));
      await session.step(58, "Then \"Share BDD-Share\" dialog should be hidden", () => shouldBe(page, el("\"Share BDD-Share\" dialog"), "hidden"));
    });
    await run.scenario("The space is shared with the second account", async () => {
      await session.step(61, "When user clicks on BDD-Share tree node inside browse tree", () => clickOn(page, el("BDD-Share tree node inside browse tree")));
      await session.step(62, "Then the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(63, "When user picks \"Share...\" from the context menu of BDD-Share tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("BDD-Share tree node inside browse tree")));
      await session.step(64, "Then share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(65, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDD-Share\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share BDD-Share\" dialog")));
      await session.step(66, "And user clicks on OK button in \"Share BDD-Share\" dialog", () => clickOn(page, el("OK button in \"Share BDD-Share\" dialog")));
      await session.step(67, "Then \"Share BDD-Share\" dialog should be hidden", () => shouldBe(page, el("\"Share BDD-Share\" dialog"), "hidden"));
      await session.step(68, "When user clicks on BDD-Share tree node inside browse tree", () => clickOn(page, el("BDD-Share tree node inside browse tree")));
      await session.step(69, "Then the sharing pane should list the sharing user", () => sharingPaneLists(page));
    });
    await run.scenario("A child space can be shared on its own", async () => {
      await session.step(72, "When user double-clicks on BDD-Share tree node inside browse tree", () => doubleClickOn(page, el("BDD-Share tree node inside browse tree")));
      await session.step(73, "Then the \"BDD-Share\" view should be current", () => viewIsCurrent(page, "BDD-Share"));
      await session.step(74, "When user picks \"Share...\" from the context menu of BDD-Share-Child link in space gallery", () => pickFromContextMenu(page, "Share...", el("BDD-Share-Child link in space gallery")));
      await session.step(75, "Then \"Share BDD-Share-Child\" dialog should be visible", () => shouldBe(page, el("\"Share BDD-Share-Child\" dialog"), "visible"));
      await session.step(76, "When user clicks on CANCEL button in \"Share BDD-Share-Child\" dialog", () => clickOn(page, el("CANCEL button in \"Share BDD-Share-Child\" dialog")));
      await session.step(77, "Then \"Share BDD-Share-Child\" dialog should be hidden", () => shouldBe(page, el("\"Share BDD-Share-Child\" dialog"), "hidden"));
    });
    await run.scenario("Dragging a parent onto its own child changes nothing", async () => {
      await session.step(80, "When user drags BDD-Share tree node inside browse tree to BDD-Share-Child tree node inside browse tree", () => dragTo(page, el("BDD-Share tree node inside browse tree"), el("BDD-Share-Child tree node inside browse tree")));
      await session.step(81, "Then BDD-Share tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Share tree node inside browse tree"), "visible"));
      await session.step(82, "And 1 space named \"BDD-Share\" should be on the server", () => spacesOnServer(page, 1, "BDD-Share"));
      await session.step(83, "And BDD-Share-Child tree node inside browse tree should be present", () => shouldBe(page, el("BDD-Share-Child tree node inside browse tree"), "present"));
    });
    await run.scenario("Deleting a shared space removes it whole", async () => {
      await session.step(86, "When user picks \"Delete Space\" from the context menu of BDD-Share tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("BDD-Share tree node inside browse tree")));
      await session.step(87, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(88, "Then 0 spaces named \"BDD-Share\" should be on the server", () => spacesOnServer(page, 0, "BDD-Share"));
      await session.step(89, "And BDD-Share tree node inside browse tree should be absent", () => shouldBe(page, el("BDD-Share tree node inside browse tree"), "absent"));
    });
    run.finish();
  });
});
