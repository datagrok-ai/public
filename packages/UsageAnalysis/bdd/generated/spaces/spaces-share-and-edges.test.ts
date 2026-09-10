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
import {browsePanelOpen, noSpaceOnServer, spacesOnServer} from '../../bindings/spaces.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, dragTo, enterInto, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sharing a space, and what it refuses", () => {
  const session = feature(test, "features/spaces/spaces-share-and-edges.feature", import.meta.url);
  test("Sharing a space, and what it refuses", {tag: ["@journey", "@spaces", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And no space named \"BDD-Share, BDD-Share-Child\" is on the server", () => noSpaceOnServer(page, "BDD-Share, BDD-Share-Child"));
    await run.scenario("A space and a child to share", async () => {
      await session.step(19, "When user picks \"Create Space...\" from the context menu of Spaces tree node", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node")));
      await session.step(20, "And user enters \"BDD-Share\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Share", el("Name input in Create Space dialog")));
      await session.step(21, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(22, "Then 1 space named \"BDD-Share\" should be on the server", () => spacesOnServer(page, 1, "BDD-Share"));
      await session.step(23, "When user picks \"Create Child Space...\" from the context menu of BDD-Share tree node", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Share tree node")));
      await session.step(24, "And user enters \"BDD-Share-Child\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Share-Child", el("Name input in Create Space dialog")));
      await session.step(25, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(26, "Then BDD-Share-Child tree node should be visible", () => shouldBe(page, el("BDD-Share-Child tree node"), "visible"));
    });
    await run.scenario("The Share dialog asks who and how much", async () => {
      await session.step(29, "When user picks \"Share...\" from the context menu of BDD-Share tree node", () => pickFromContextMenu(page, "Share...", el("BDD-Share tree node")));
      await session.step(30, "Then \"Share BDD-Share\" dialog should be visible", () => shouldBe(page, el("\"Share BDD-Share\" dialog"), "visible"));
      await session.step(31, "And \"User, group, or email\" input in \"Share BDD-Share\" dialog should be visible", () => shouldBe(page, el("\"User, group, or email\" input in \"Share BDD-Share\" dialog"), "visible"));
      await session.step(32, "And share access selector should be visible", () => shouldBe(page, el("share access selector"), "visible"));
      await session.step(33, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(34, "When user clicks on CANCEL button in \"Share BDD-Share\" dialog", () => clickOn(page, el("CANCEL button in \"Share BDD-Share\" dialog")));
      await session.step(35, "Then \"Share BDD-Share\" dialog should be hidden", () => shouldBe(page, el("\"Share BDD-Share\" dialog"), "hidden"));
    });
    await run.scenario("A child space can be shared on its own", async () => {
      await session.step(38, "When user double-clicks on BDD-Share tree node", () => doubleClickOn(page, el("BDD-Share tree node")));
      await session.step(39, "Then the \"BDD-Share\" view should be current", () => viewIsCurrent(page, "BDD-Share"));
      await session.step(40, "When user picks \"Share...\" from the context menu of BDD-Share-Child link in space gallery", () => pickFromContextMenu(page, "Share...", el("BDD-Share-Child link in space gallery")));
      await session.step(41, "Then \"Share BDD-Share-Child\" dialog should be visible", () => shouldBe(page, el("\"Share BDD-Share-Child\" dialog"), "visible"));
      await session.step(42, "When user clicks on CANCEL button in \"Share BDD-Share-Child\" dialog", () => clickOn(page, el("CANCEL button in \"Share BDD-Share-Child\" dialog")));
      await session.step(43, "Then \"Share BDD-Share-Child\" dialog should be hidden", () => shouldBe(page, el("\"Share BDD-Share-Child\" dialog"), "hidden"));
    });
    await run.scenario("Dragging a parent onto its own child changes nothing", async () => {
      await session.step(46, "When user drags BDD-Share tree node to BDD-Share-Child tree node", () => dragTo(page, el("BDD-Share tree node"), el("BDD-Share-Child tree node")));
      await session.step(47, "Then BDD-Share tree node should be visible", () => shouldBe(page, el("BDD-Share tree node"), "visible"));
      await session.step(48, "And 1 space named \"BDD-Share\" should be on the server", () => spacesOnServer(page, 1, "BDD-Share"));
      await session.step(49, "And BDD-Share-Child tree node should be present", () => shouldBe(page, el("BDD-Share-Child tree node"), "present"));
    });
    run.finish();
  });
});
