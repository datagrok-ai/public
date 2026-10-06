/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-share-and-edges.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, isExpanded, shouldBe, shouldContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, childSpaceOnServer, contextPanelOpen, contextPanelShows, noSpaceOnServer, pickSharingUser, sharingPaneLists, sharingPaneListsNot, spaceOnServer, spacesOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sharing a space", () => {
  const session = feature(test, "features/spaces/spaces-share-and-edges.feature", import.meta.url);
  test("Sharing a space", {tag: ["@journey", "@spaces", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(37, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(39, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(40, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(41, "And no space named \"BDD-Share, BDD-Share-Child\" is on the server", () => noSpaceOnServer(page, "BDD-Share, BDD-Share-Child"));
    await run.scenario("A space and a child to share", async () => {
      await session.step(44, "Given a space named \"BDD-Share\" is on the server", () => spaceOnServer(page, "BDD-Share"));
      await session.step(45, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(46, "And 1 space named \"BDD-Share\" should be on the server", () => spacesOnServer(page, 1, "BDD-Share"));
      await session.step(47, "Given a space named \"BDD-Share-Child\" under \"BDD-Share\" is on the server", () => childSpaceOnServer(page, "BDD-Share-Child", "BDD-Share"));
      await session.step(48, "When user double-clicks on Spaces---BDD-Share tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-Share tree node inside browse tree")));
      await session.step(49, "Then the \"BDD-Share\" view should be current", () => viewIsCurrent(page, "BDD-Share"));
      await session.step(50, "And BDD-Share-Child link in gallery should be visible", () => shouldBe(page, el("BDD-Share-Child link in gallery"), "visible"));
    });
    await run.scenario("The Share dialog asks who and how much", async () => {
      await session.step(53, "When user picks \"Share...\" from the context menu of Spaces---BDD-Share tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("Spaces---BDD-Share tree node inside browse tree")));
      await session.step(54, "Then \"Share BDD-Share\" dialog should be visible", () => shouldBe(page, el("\"Share BDD-Share\" dialog"), "visible"));
      await session.step(55, "And \"User, group, or email\" input in \"Share BDD-Share\" dialog should be visible", () => shouldBe(page, el("\"User, group, or email\" input in \"Share BDD-Share\" dialog"), "visible"));
      await session.step(56, "And share access selector should be visible", () => shouldBe(page, el("share access selector"), "visible"));
      await session.step(57, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(58, "When user clicks on CANCEL button in \"Share BDD-Share\" dialog", () => clickOn(page, el("CANCEL button in \"Share BDD-Share\" dialog")));
      await session.step(59, "Then \"Share BDD-Share\" dialog should be hidden", () => shouldBe(page, el("\"Share BDD-Share\" dialog"), "hidden"));
    });
    await run.scenario("The space is shared with the second account", async () => {
      await session.step(62, "When user clicks on Spaces---BDD-Share tree node inside browse tree", () => clickOn(page, el("Spaces---BDD-Share tree node inside browse tree")));
      await session.step(63, "Then the context panel should show \"BDD-Share\"", () => contextPanelShows(page, "BDD-Share"));
      await session.step(64, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(65, "When user picks \"Share...\" from the context menu of Spaces---BDD-Share tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("Spaces---BDD-Share tree node inside browse tree")));
      await session.step(66, "Then share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(67, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDD-Share\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share BDD-Share\" dialog")));
      await session.step(68, "And user unchecks \"Send notifications\" input in \"Share BDD-Share\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share BDD-Share\" dialog")));
      await session.step(69, "And user clicks on OK button in \"Share BDD-Share\" dialog", () => clickOn(page, el("OK button in \"Share BDD-Share\" dialog")));
      await session.step(70, "Then \"Share BDD-Share\" dialog should be hidden", () => shouldBe(page, el("\"Share BDD-Share\" dialog"), "hidden"));
      await session.step(71, "When user clicks on Spaces---BDD-Share tree node inside browse tree", () => clickOn(page, el("Spaces---BDD-Share tree node inside browse tree")));
      await session.step(72, "Then the context panel should show \"BDD-Share\"", () => contextPanelShows(page, "BDD-Share"));
      await session.step(73, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
    });
    await run.scenario("A child space can be shared on its own", async () => {
      await session.step(76, "When user double-clicks on Spaces---BDD-Share tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-Share tree node inside browse tree")));
      await session.step(77, "Then the \"BDD-Share\" view should be current", () => viewIsCurrent(page, "BDD-Share"));
      await session.step(78, "When user picks \"Share...\" from the context menu of BDD-Share-Child link in gallery", () => pickFromContextMenu(page, "Share...", el("BDD-Share-Child link in gallery")));
      await session.step(79, "Then \"Share BDD-Share-Child\" dialog should be visible", () => shouldBe(page, el("\"Share BDD-Share-Child\" dialog"), "visible"));
      await session.step(80, "When user clicks on CANCEL button in \"Share BDD-Share-Child\" dialog", () => clickOn(page, el("CANCEL button in \"Share BDD-Share-Child\" dialog")));
      await session.step(81, "Then \"Share BDD-Share-Child\" dialog should be hidden", () => shouldBe(page, el("\"Share BDD-Share-Child\" dialog"), "hidden"));
    });
    run.finish();
  });
});
