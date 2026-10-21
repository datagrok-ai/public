/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/sharing-connection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog, sharing.browse-shared-with-me, core#977]
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
import {clickOn, isExpanded, shouldBe, shouldContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, connectionOnServer, contextPanelOpen, contextPanelShows, dialogCloses, expandShared, pickSharingUser, refreshBrowse, removeSharingUser, runningAccountSignedIn, sharedShouldBe, sharingPaneLists, sharingPaneListsNot, sharingUserLevel, sharingUserNoNotifications, sharingUserSignedIn, signInAsSelf, signInAsSharingUser} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A connection, shared and taken back", () => {
  const session = feature(test, "features/sharing/sharing-connection.feature", import.meta.url);
  test("A connection, shared and taken back", {tag: ["@journey", "@serial", "@realizes:sharing.share-dialog", "@realizes:sharing.browse-shared-with-me", "@known-failure:66", "@realizes:core#977"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(32, "And a \"Postgres\" connection named \"ShareConn1\" is on the server", () => connectionOnServer(page, "Postgres", "ShareConn1"));
    await session.step(33, "And a \"Postgres\" connection named \"ShareConn2\" is on the server", () => connectionOnServer(page, "Postgres", "ShareConn2"));
    await session.step(34, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("The owner shares both connections from their Browse nodes", async () => {
      await session.step(38, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(39, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(40, "When user picks \"Share...\" from the context menu of Databases---Postgres---ShareConn2 tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("Databases---Postgres---ShareConn2 tree node inside browse tree")));
      await session.step(41, "And user picks the sharing user in \"User, group, or email\" input in \"Share ShareConn2\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareConn2\" dialog")));
      await session.step(42, "And user unchecks \"Send notifications\" input in \"Share ShareConn2\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareConn2\" dialog")));
      await session.step(43, "And user clicks on OK button in \"Share ShareConn2\" dialog", () => clickOn(page, el("OK button in \"Share ShareConn2\" dialog")));
      await session.step(44, "Then the \"Share ShareConn2\" dialog should close", () => dialogCloses(page, "Share ShareConn2"));
      await session.step(45, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(46, "When user picks \"Share...\" from the context menu of Databases---Postgres---ShareConn1 tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("Databases---Postgres---ShareConn1 tree node inside browse tree")));
      await session.step(47, "Then \"Share ShareConn1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareConn1\" dialog"), "visible"));
      await session.step(48, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(49, "When user picks the sharing user in \"User, group, or email\" input in \"Share ShareConn1\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareConn1\" dialog")));
      await session.step(50, "And user unchecks \"Send notifications\" input in \"Share ShareConn1\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareConn1\" dialog")));
      await session.step(51, "And user clicks on OK button in \"Share ShareConn1\" dialog", () => clickOn(page, el("OK button in \"Share ShareConn1\" dialog")));
      await session.step(52, "Then the \"Share ShareConn1\" dialog should close", () => dialogCloses(page, "Share ShareConn1"));
      await session.step(53, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(54, "When user clicks on Databases---Postgres---ShareConn1 tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---ShareConn1 tree node inside browse tree")));
      await session.step(55, "Then the context panel should show \"ShareConn1\"", () => contextPanelShows(page, "ShareConn1"));
      await session.step(56, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(57, "And Sharing section in context panel should contain text \"You are the owner\"", () => shouldContainText(page, el("Sharing section in context panel"), "You are the owner"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
      await session.step(59, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Sharing pane names the level of the connection's grant (core#977)", async () => {
      await session.step(63, "When user clicks on Databases---Postgres---ShareConn1 tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---ShareConn1 tree node inside browse tree")));
      await session.step(64, "Then the context panel should show \"ShareConn1\"", () => contextPanelShows(page, "ShareConn1"));
      await session.step(65, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(66, "And Sharing section in context panel should contain text \"can view and use\"", () => shouldContainText(page, el("Sharing section in context panel"), "can view and use"));
    }, {knownFailure: true, line: 66});
    await run.scenario("The second account finds the connection under Shared with me", async () => {
      await session.step(69, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(70, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(71, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(72, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(73, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(74, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(75, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(76, "Then \"ShareConn1\" shared by the running account should be visible", () => sharedShouldBe(page, "ShareConn1", "visible"));
      await session.step(77, "And \"ShareConn2\" shared by the running account should be visible", () => sharedShouldBe(page, "ShareConn2", "visible"));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The owner takes the share back", async () => {
      await session.step(81, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(82, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(83, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(84, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(85, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(86, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(87, "When user picks \"Share...\" from the context menu of Databases---Postgres---ShareConn1 tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("Databases---Postgres---ShareConn1 tree node inside browse tree")));
      await session.step(88, "Then \"Share ShareConn1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareConn1\" dialog"), "visible"));
      await session.step(89, "And the access level of the sharing user in \"Share ShareConn1\" dialog should be \"View and use\"", () => sharingUserLevel(page, el("\"Share ShareConn1\" dialog"), "View and use"));
      await session.step(90, "When user removes the sharing user from \"Share ShareConn1\" dialog", () => removeSharingUser(page, el("\"Share ShareConn1\" dialog")));
      await session.step(91, "And user clicks on OK button in \"Share ShareConn1\" dialog", () => clickOn(page, el("OK button in \"Share ShareConn1\" dialog")));
      await session.step(92, "Then the \"Share ShareConn1\" dialog should close", () => dialogCloses(page, "Share ShareConn1"));
      await session.step(93, "When user clicks on Databases---Postgres---ShareConn1 tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---ShareConn1 tree node inside browse tree")));
      await session.step(94, "Then the context panel should show \"ShareConn1\"", () => contextPanelShows(page, "ShareConn1"));
      await session.step(95, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The second account no longer has the connection", async () => {
      await session.step(99, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(100, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(101, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(102, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(103, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(104, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(105, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(106, "Then \"ShareConn2\" shared by the running account should be visible", () => sharedShouldBe(page, "ShareConn2", "visible"));
      await session.step(107, "And \"ShareConn1\" shared by the running account should be absent", () => sharedShouldBe(page, "ShareConn1", "absent"));
      await session.step(108, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(109, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    run.finish();
  });
});
