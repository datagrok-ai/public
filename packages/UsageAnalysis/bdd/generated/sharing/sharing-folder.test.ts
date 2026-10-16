/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/sharing-folder.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog, sharing.browse-shared-with-me, public#4268]
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
import {browsePanelOpen, connectionsOnServer, dialogCloses, expandShared, openShared, ownConnectionGone, pickSharingUser, refreshBrowse, runningAccountSignedIn, sharingUserLevel, sharingUserNoNotifications, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fixtureInUsersFiles} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A folder in your files, shared", () => {
  const session = feature(test, "features/sharing/sharing-folder.feature", import.meta.url);
  test("A folder in your files, shared", {tag: ["@journey", "@serial", "@realizes:sharing.share-dialog", "@realizes:sharing.browse-shared-with-me", "@known-failure:92", "@realizes:public#4268"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(30, "And the user's own connection \"MyFilesShareFolder1\" is removed now and at feature end", () => ownConnectionGone(page, "MyFilesShareFolder1"));
    await session.step(31, "And the \"fixtures/cars-small.csv\" file of the project is in the user's files as \"ShareFolder1/data1.csv\"", () => fixtureInUsersFiles(page, "fixtures/cars-small.csv", "ShareFolder1/data1.csv"));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(34, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(35, "And Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
    await run.scenario("Share folder opens and cancels cleanly (GROK-19264)", async () => {
      await session.step(38, "When user picks \"Share folder...\" from the context menu of Files---My-files---ShareFolder1 tree node inside browse tree", () => pickFromContextMenu(page, "Share folder...", el("Files---My-files---ShareFolder1 tree node inside browse tree")));
      await session.step(39, "Then \"Share folder\" dialog should be visible", () => shouldBe(page, el("\"Share folder\" dialog"), "visible"));
      await session.step(40, "And link icon in \"Share folder\" dialog should be hidden", () => shouldBe(page, el("link icon in \"Share folder\" dialog"), "hidden"));
      await session.step(41, "And \"User, group, or email\" input in \"Share folder\" dialog should be visible", () => shouldBe(page, el("\"User, group, or email\" input in \"Share folder\" dialog"), "visible"));
      await session.step(42, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(43, "And OK button in \"Share folder\" dialog should be visible", () => shouldBe(page, el("OK button in \"Share folder\" dialog"), "visible"));
      await session.step(44, "And CANCEL button in \"Share folder\" dialog should be visible", () => shouldBe(page, el("CANCEL button in \"Share folder\" dialog"), "visible"));
      await session.step(45, "When user clicks on CANCEL button in \"Share folder\" dialog", () => clickOn(page, el("CANCEL button in \"Share folder\" dialog")));
      await session.step(46, "Then the \"Share folder\" dialog should close", () => dialogCloses(page, "Share folder"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
      await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(49, "And 0 connections named \"MyFilesShareFolder1\" should be on the server", () => connectionsOnServer(page, 0, "MyFilesShareFolder1"));
      await session.step(50, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(51, "Given Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
      await session.step(52, "Then Files---My-files---ShareFolder1 tree node inside browse tree should be visible", () => shouldBe(page, el("Files---My-files---ShareFolder1 tree node inside browse tree"), "visible"));
      await session.step(53, "And \"My files/ShareFolder1\" tree node inside browse tree should be absent", () => shouldBe(page, el("\"My files/ShareFolder1\" tree node inside browse tree"), "absent"));
    });
    await run.scenario("The folder is shared with the second account (GROK-12795, GROK-11266)", async () => {
      await session.step(56, "When user picks \"Share folder...\" from the context menu of Files---My-files---ShareFolder1 tree node inside browse tree", () => pickFromContextMenu(page, "Share folder...", el("Files---My-files---ShareFolder1 tree node inside browse tree")));
      await session.step(57, "Then \"Share folder\" dialog should be visible", () => shouldBe(page, el("\"Share folder\" dialog"), "visible"));
      await session.step(58, "When user picks the sharing user in \"User, group, or email\" input in \"Share folder\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share folder\" dialog")));
      await session.step(59, "And user unchecks \"Send notifications\" input in \"Share folder\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share folder\" dialog")));
      await session.step(60, "And user clicks on OK button in \"Share folder\" dialog", () => clickOn(page, el("OK button in \"Share folder\" dialog")));
      await session.step(61, "Then the \"Share folder\" dialog should close", () => dialogCloses(page, "Share folder"));
      await session.step(62, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(64, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(65, "Then \"My files/ShareFolder1\" tree node inside browse tree should be visible", () => shouldBe(page, el("\"My files/ShareFolder1\" tree node inside browse tree"), "visible"));
      await session.step(66, "When user clicks on \"My files/ShareFolder1\" tree node inside browse tree", () => clickOn(page, el("\"My files/ShareFolder1\" tree node inside browse tree")));
      await session.step(67, "Then data1.csv link in gallery should be visible", () => shouldBe(page, el("data1.csv link in gallery"), "visible"));
      await session.step(68, "Given Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
      await session.step(69, "When user picks \"Share folder...\" from the context menu of Files---My-files---ShareFolder1 tree node inside browse tree", () => pickFromContextMenu(page, "Share folder...", el("Files---My-files---ShareFolder1 tree node inside browse tree")));
      await session.step(70, "Then \"Share folder\" dialog should be visible", () => shouldBe(page, el("\"Share folder\" dialog"), "visible"));
      await session.step(71, "And link icon in \"Share folder\" dialog should be visible", () => shouldBe(page, el("link icon in \"Share folder\" dialog"), "visible"));
      await session.step(72, "And the access level of the sharing user in \"Share folder\" dialog should be \"View and use\"", () => sharingUserLevel(page, el("\"Share folder\" dialog"), "View and use"));
      await session.step(73, "When user clicks on CANCEL button in \"Share folder\" dialog", () => clickOn(page, el("CANCEL button in \"Share folder\" dialog")));
      await session.step(74, "Then the \"Share folder\" dialog should close", () => dialogCloses(page, "Share folder"));
    });
    await run.scenario("The second account finds the folder under Shared with me", async () => {
      await session.step(77, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(78, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(79, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(80, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(81, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(82, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(83, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(84, "And user double-clicks on \"My files/ShareFolder1\" shared by the running account", () => openShared(page, "My files/ShareFolder1"));
      await session.step(85, "Then the \"My files/ShareFolder1\" view should be current", () => viewIsCurrent(page, "My files/ShareFolder1"));
      await session.step(86, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(87, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    });
    await run.scenario("The second account sees the folder's files, not a credentials form (public#4268)", async () => {
      await session.step(91, "When user double-clicks on \"My files/ShareFolder1\" shared by the running account", () => openShared(page, "My files/ShareFolder1"));
      await session.step(92, "Then data1.csv link in gallery should be visible", () => shouldBe(page, el("data1.csv link in gallery"), "visible"));
      await session.step(93, "And \"Credentials required\" text should be absent", () => shouldBe(page, el("\"Credentials required\" text"), "absent"));
    }, {knownFailure: true, line: 92});
    await run.scenario("The owner deletes the folder's connection and the folder", async () => {
      await session.step(96, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(97, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(98, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(99, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(100, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(101, "When user picks \"Delete...\" from the context menu of \"My files/ShareFolder1\" tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el("\"My files/ShareFolder1\" tree node inside browse tree")));
      await session.step(102, "Then \"Are you sure?\" dialog should contain text 'Delete connection \"My files/ShareFolder1\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "Delete connection \"My files/ShareFolder1\"?"));
      await session.step(103, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(104, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(105, "And 0 connections named \"MyFilesShareFolder1\" should be on the server", () => connectionsOnServer(page, 0, "MyFilesShareFolder1"));
      await session.step(106, "Given Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
      await session.step(107, "When user picks \"Delete...\" from the context menu of Files---My-files---ShareFolder1 tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el("Files---My-files---ShareFolder1 tree node inside browse tree")));
      await session.step(108, "Then \"Are you sure?\" dialog should contain text \"ShareFolder1\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "ShareFolder1"));
      await session.step(109, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(110, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(111, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(112, "Given Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
      await session.step(113, "Then \"My files/ShareFolder1\" tree node inside browse tree should be absent", () => shouldBe(page, el("\"My files/ShareFolder1\" tree node inside browse tree"), "absent"));
      await session.step(114, "And Files---My-files---ShareFolder1 tree node inside browse tree should be absent", () => shouldBe(page, el("Files---My-files---ShareFolder1 tree node inside browse tree"), "absent"));
    });
    run.finish();
  });
});
