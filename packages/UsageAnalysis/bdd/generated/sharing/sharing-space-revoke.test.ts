/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/sharing-space-revoke.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog, sharing.browse-shared-with-me, views.space]
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
import {clickOn, isExpanded, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, expandShared, pickSharingUser, refreshBrowse, removeSharingUser, runningAccountSignedIn, sharedShouldBe, sharingPaneLists, sharingPaneListsNot, sharingUserLevel, sharingUserNoNotifications, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, spaceOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A space, shared and taken back", () => {
  const session = feature(test, "features/sharing/sharing-space-revoke.feature", import.meta.url);
  test("A space, shared and taken back", {tag: ["@journey", "@serial", "@realizes:sharing.share-dialog", "@realizes:sharing.browse-shared-with-me", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(24, "And a space named \"ShareRevokeSpace1\" is on the server", () => spaceOnServer(page, "ShareRevokeSpace1"));
    await session.step(25, "And a space named \"ShareRevokeSpace2\" is on the server", () => spaceOnServer(page, "ShareRevokeSpace2"));
    await session.step(26, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(28, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(29, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await run.scenario("The owner shares both spaces", async () => {
      await session.step(32, "When user picks \"Share...\" from the context menu of Spaces---ShareRevokeSpace2 tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("Spaces---ShareRevokeSpace2 tree node inside browse tree")));
      await session.step(33, "And user picks the sharing user in \"User, group, or email\" input in \"Share ShareRevokeSpace2\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareRevokeSpace2\" dialog")));
      await session.step(34, "And user unchecks \"Send notifications\" input in \"Share ShareRevokeSpace2\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareRevokeSpace2\" dialog")));
      await session.step(35, "And user clicks on OK button in \"Share ShareRevokeSpace2\" dialog", () => clickOn(page, el("OK button in \"Share ShareRevokeSpace2\" dialog")));
      await session.step(36, "Then the \"Share ShareRevokeSpace2\" dialog should close", () => dialogCloses(page, "Share ShareRevokeSpace2"));
      await session.step(37, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(38, "When user picks \"Share...\" from the context menu of Spaces---ShareRevokeSpace1 tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("Spaces---ShareRevokeSpace1 tree node inside browse tree")));
      await session.step(39, "And user picks the sharing user in \"User, group, or email\" input in \"Share ShareRevokeSpace1\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareRevokeSpace1\" dialog")));
      await session.step(40, "And user unchecks \"Send notifications\" input in \"Share ShareRevokeSpace1\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareRevokeSpace1\" dialog")));
      await session.step(41, "And user clicks on OK button in \"Share ShareRevokeSpace1\" dialog", () => clickOn(page, el("OK button in \"Share ShareRevokeSpace1\" dialog")));
      await session.step(42, "Then the \"Share ShareRevokeSpace1\" dialog should close", () => dialogCloses(page, "Share ShareRevokeSpace1"));
      await session.step(43, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
      await session.step(45, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The second account sees both spaces under Shared with me", async () => {
      await session.step(48, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(49, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(50, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(51, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(52, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(53, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(54, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(55, "Then \"ShareRevokeSpace1\" shared by the running account should be visible", () => sharedShouldBe(page, "ShareRevokeSpace1", "visible"));
      await session.step(56, "And \"ShareRevokeSpace2\" shared by the running account should be visible", () => sharedShouldBe(page, "ShareRevokeSpace2", "visible"));
    });
    await run.scenario("The owner takes the share of one space back", async () => {
      await session.step(59, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(60, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(61, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(62, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(63, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(64, "When user clicks on Spaces---ShareRevokeSpace1 tree node inside browse tree", () => clickOn(page, el("Spaces---ShareRevokeSpace1 tree node inside browse tree")));
      await session.step(65, "Then the context panel should show \"ShareRevokeSpace1\"", () => contextPanelShows(page, "ShareRevokeSpace1"));
      await session.step(66, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(67, "When user picks \"Share...\" from the context menu of Spaces---ShareRevokeSpace1 tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("Spaces---ShareRevokeSpace1 tree node inside browse tree")));
      await session.step(68, "Then \"Share ShareRevokeSpace1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareRevokeSpace1\" dialog"), "visible"));
      await session.step(69, "And the access level of the sharing user in \"Share ShareRevokeSpace1\" dialog should be \"View and use\"", () => sharingUserLevel(page, el("\"Share ShareRevokeSpace1\" dialog"), "View and use"));
      await session.step(70, "When user removes the sharing user from \"Share ShareRevokeSpace1\" dialog", () => removeSharingUser(page, el("\"Share ShareRevokeSpace1\" dialog")));
      await session.step(71, "And user clicks on OK button in \"Share ShareRevokeSpace1\" dialog", () => clickOn(page, el("OK button in \"Share ShareRevokeSpace1\" dialog")));
      await session.step(72, "Then the \"Share ShareRevokeSpace1\" dialog should close", () => dialogCloses(page, "Share ShareRevokeSpace1"));
      await session.step(73, "When user clicks on Spaces---ShareRevokeSpace1 tree node inside browse tree", () => clickOn(page, el("Spaces---ShareRevokeSpace1 tree node inside browse tree")));
      await session.step(74, "Then the context panel should show \"ShareRevokeSpace1\"", () => contextPanelShows(page, "ShareRevokeSpace1"));
      await session.step(75, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The second account no longer has the space", async () => {
      await session.step(79, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(80, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(81, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(82, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(83, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(84, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(85, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(86, "Then \"ShareRevokeSpace2\" shared by the running account should be visible", () => sharedShouldBe(page, "ShareRevokeSpace2", "visible"));
      await session.step(87, "And \"ShareRevokeSpace1\" shared by the running account should be absent", () => sharedShouldBe(page, "ShareRevokeSpace1", "absent"));
      await session.step(88, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(89, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    run.finish();
  });
});
