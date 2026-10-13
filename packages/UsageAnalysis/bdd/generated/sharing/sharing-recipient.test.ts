/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/sharing-recipient.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog, sharing.browse-shared-with-me, views.queries, views.space, core#693]
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
import {clearField, clickOn, enterInto, isExpanded, pressKey, selectIn, shouldBe, shouldContainText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, dragShared, expandShared, openShared, pickFromSharedMenu, pickRunningAccount, pickSharingUser, queryOnServer, refreshBrowse, removeSharingUser, runningAccountSignedIn, sharedShouldBe, sharingPaneLists, sharingPaneListsNot, sharingUserLevel, sharingUserNoNotifications, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, spaceOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fixtureInSpace} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {errorBalloonText, infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("What the other user gets, and losing it", () => {
  const session = feature(test, "features/sharing/sharing-recipient.feature", import.meta.url);
  test("What the other user gets, and losing it", {tag: ["@journey", "@serial", "@realizes:sharing.share-dialog", "@realizes:sharing.browse-shared-with-me", "@realizes:views.queries", "@realizes:views.space", "@realizes:core#693"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(36, "Given user is logged in", () => loggedIn(page));
    await session.step(37, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(38, "And a query \"ShareRecipient1\" on \"System:Datagrok\" reads \"select 1 as x\"", () => queryOnServer(page, "ShareRecipient1", "System:Datagrok", "select 1 as x"));
    await session.step(39, "And a query \"ShareRecipient2\" on \"System:Datagrok\" reads \"select 2 as y\"", () => queryOnServer(page, "ShareRecipient2", "System:Datagrok", "select 2 as y"));
    await session.step(40, "And a space named \"ShareSpace1\" is on the server", () => spaceOnServer(page, "ShareSpace1"));
    await session.step(41, "And the \"fixtures/cars-small.csv\" file of the project is in the space \"ShareSpace1\" as \"ShareSpaceStay1.csv\"", () => fixtureInSpace(page, "fixtures/cars-small.csv", "ShareSpace1", "ShareSpaceStay1.csv"));
    await session.step(42, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(43, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("The second account's starting state: ShareRecipient2 and ShareSpace1 shared with it", async () => {
      await session.step(46, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(47, "When user picks \"Share...\" from the context menu of Spaces---ShareSpace1 tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("Spaces---ShareSpace1 tree node inside browse tree")));
      await session.step(48, "And user picks the sharing user in \"User, group, or email\" input in \"Share ShareSpace1\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareSpace1\" dialog")));
      await session.step(49, "And user unchecks \"Send notifications\" input in \"Share ShareSpace1\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareSpace1\" dialog")));
      await session.step(50, "And user clicks on OK button in \"Share ShareSpace1\" dialog", () => clickOn(page, el("OK button in \"Share ShareSpace1\" dialog")));
      await session.step(51, "Then the \"Share ShareSpace1\" dialog should close", () => dialogCloses(page, "Share ShareSpace1"));
      await session.step(52, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(53, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
      await session.step(54, "When user clicks on Platform---Functions---Queries tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Queries tree node inside browse tree")));
      await session.step(55, "Then the \"Queries\" view should be current", () => viewIsCurrent(page, "Queries"));
      await session.step(56, "When user enters \"ShareRecipient2\" into gallery search", () => enterInto(page, "ShareRecipient2", el("gallery search")));
      await session.step(57, "And user picks \"Share...\" from the context menu of \"ShareRecipient2\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareRecipient2\" gallery card")));
      await session.step(58, "And user picks the sharing user in \"User, group, or email\" input in \"Share ShareRecipient2\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareRecipient2\" dialog")));
      await session.step(59, "And user unchecks \"Send notifications\" input in \"Share ShareRecipient2\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareRecipient2\" dialog")));
      await session.step(60, "And user clicks on OK button in \"Share ShareRecipient2\" dialog", () => clickOn(page, el("OK button in \"Share ShareRecipient2\" dialog")));
      await session.step(61, "Then the \"Share ShareRecipient2\" dialog should close", () => dialogCloses(page, "Share ShareRecipient2"));
      await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The owner shares the query to view and use (GROK-19235)", async () => {
      await session.step(65, "When user enters \"ShareRecipient1\" into gallery search", () => enterInto(page, "ShareRecipient1", el("gallery search")));
      await session.step(66, "And user picks \"Share...\" from the context menu of \"ShareRecipient1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareRecipient1\" gallery card")));
      await session.step(67, "Then \"Share ShareRecipient1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareRecipient1\" dialog"), "visible"));
      await session.step(68, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(69, "When user picks the sharing user in \"User, group, or email\" input in \"Share ShareRecipient1\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareRecipient1\" dialog")));
      await session.step(70, "And user unchecks \"Send notifications\" input in \"Share ShareRecipient1\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareRecipient1\" dialog")));
      await session.step(71, "And user clicks on OK button in \"Share ShareRecipient1\" dialog", () => clickOn(page, el("OK button in \"Share ShareRecipient1\" dialog")));
      await session.step(72, "Then the \"Share ShareRecipient1\" dialog should close", () => dialogCloses(page, "Share ShareRecipient1"));
      await session.step(73, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(74, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(75, "When user clicks on \"ShareRecipient1\" gallery card", () => clickOn(page, el("\"ShareRecipient1\" gallery card")));
      await session.step(76, "Then the context panel should show \"ShareRecipient1\"", () => contextPanelShows(page, "ShareRecipient1"));
      await session.step(77, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(78, "And Sharing section in context panel should contain text \"You are the owner\"", () => shouldContainText(page, el("Sharing section in context panel"), "You are the owner"));
      await session.step(79, "And Sharing section in context panel should contain text \"On ShareRecipient1\"", () => shouldContainText(page, el("Sharing section in context panel"), "On ShareRecipient1"));
      await session.step(80, "And Sharing section in context panel should contain text \"can view and use\"", () => shouldContainText(page, el("Sharing section in context panel"), "can view and use"));
    });
    await run.scenario("The second account finds the query under Shared with me and runs it (GROK-19401, GROK-18525, GROK-18210)", async () => {
      await session.step(83, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(84, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(85, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(86, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(87, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(88, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(89, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(90, "And user double-clicks on \"ShareRecipient1\" shared by the running account", () => openShared(page, "ShareRecipient1"));
      await session.step(91, "Then the \"ShareRecipient1\" view should be current", () => viewIsCurrent(page, "ShareRecipient1"));
      await session.step(92, "And the table should have 1 row", () => rowCount(page, 1));
      await session.step(93, "And the table should have the columns \"x\"", () => columnsExactly(page, "x"));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
      await session.step(95, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The second account cannot share the query further", async () => {
      await session.step(98, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(99, "And \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(100, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(101, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(102, "And user picks \"Share...\" from the context menu of \"ShareRecipient1\" shared by the running account", () => pickFromSharedMenu(page, "Share...", "ShareRecipient1"));
      await session.step(103, "Then \"Share ShareRecipient1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareRecipient1\" dialog"), "visible"));
      await session.step(104, "When user picks the running account in \"User, group, or email\" input in \"Share ShareRecipient1\" dialog", () => pickRunningAccount(page, el("\"User, group, or email\" input in \"Share ShareRecipient1\" dialog")));
      await session.step(105, "And user unchecks \"Send notifications\" input in \"Share ShareRecipient1\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareRecipient1\" dialog")));
      await session.step(106, "And user clicks on OK button in \"Share ShareRecipient1\" dialog", () => clickOn(page, el("OK button in \"Share ShareRecipient1\" dialog")));
      await session.step(107, "Then an error balloon containing \"You don't have a permission to share this object\" should have been shown", () => errorBalloonText(page, "You don't have a permission to share this object"));
      await session.step(108, "And \"Share ShareRecipient1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareRecipient1\" dialog"), "visible"));
      await session.step(109, "When user clicks on CANCEL button in \"Share ShareRecipient1\" dialog", () => clickOn(page, el("CANCEL button in \"Share ShareRecipient1\" dialog")));
      await session.step(110, "Then the \"Share ShareRecipient1\" dialog should close", () => dialogCloses(page, "Share ShareRecipient1"));
    });
    await run.scenario("The second account cannot Move a query into a space shared for viewing", async () => {
      await session.step(113, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(114, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(115, "And \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(116, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(117, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(118, "And user drags \"ShareRecipient2\" shared by the running account to \"ShareSpace1\" shared by the running account", () => dragShared(page, "ShareRecipient2", "ShareSpace1"));
      await session.step(119, "Then \"Move entity\" dialog should be visible", () => shouldBe(page, el("\"Move entity\" dialog"), "visible"));
      await session.step(120, "And \"Move entity\" dialog should contain text \"to ShareSpace1 project\"", () => shouldContainText(page, el("\"Move entity\" dialog"), "to ShareSpace1 project"));
      await session.step(121, "And \"Move entity\" dialog should contain text \"ShareRecipient2\"", () => shouldContainText(page, el("\"Move entity\" dialog"), "ShareRecipient2"));
      await session.step(122, "And choice input in \"Move entity\" dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in \"Move entity\" dialog"), "Link"));
      await session.step(123, "When user selects \"Move\" in choice input in \"Move entity\" dialog", () => selectIn(page, "Move", el("choice input in \"Move entity\" dialog")));
      await session.step(124, "Then choice input in \"Move entity\" dialog should have value \"Move\"", () => shouldHaveValue(page, el("choice input in \"Move entity\" dialog"), "Move"));
      await session.step(125, "When user clicks on YES button in \"Move entity\" dialog", () => clickOn(page, el("YES button in \"Move entity\" dialog")));
      await session.step(126, "Then an error balloon containing \"insufficient privileges to edit the project\" should have been shown", () => errorBalloonText(page, "insufficient privileges to edit the project"));
      await session.step(127, "And \"Move entity\" dialog should be visible", () => shouldBe(page, el("\"Move entity\" dialog"), "visible"));
      await session.step(128, "When user clicks on CANCEL button in \"Move entity\" dialog", () => clickOn(page, el("CANCEL button in \"Move entity\" dialog")));
      await session.step(129, "Then the \"Move entity\" dialog should close", () => dialogCloses(page, "Move entity"));
    });
    await run.scenario("The second account asks for a Link into the same space", async () => {
      await session.step(132, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(133, "And user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(134, "And user drags \"ShareRecipient2\" shared by the running account to \"ShareSpace1\" shared by the running account", () => dragShared(page, "ShareRecipient2", "ShareSpace1"));
      await session.step(135, "Then \"Move entity\" dialog should be visible", () => shouldBe(page, el("\"Move entity\" dialog"), "visible"));
      await session.step(136, "And choice input in \"Move entity\" dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in \"Move entity\" dialog"), "Link"));
    });
    await run.scenario("The Link is refused, and the space does not get the query (core#693)", async () => {
      await session.step(140, "When user clicks on YES button in \"Move entity\" dialog", () => clickOn(page, el("YES button in \"Move entity\" dialog")));
      await session.step(141, "Then an error balloon containing \"insufficient privileges to edit the project\" should have been shown", () => errorBalloonText(page, "insufficient privileges to edit the project"));
      await session.step(142, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(143, "And user double-clicks on \"ShareSpace1\" shared by the running account", () => openShared(page, "ShareSpace1"));
      await session.step(144, "Then the \"ShareSpace1\" view should be current", () => viewIsCurrent(page, "ShareSpace1"));
      await session.step(145, "And \"ShareSpaceStay1.csv\" link in gallery should be visible", () => shouldBe(page, el("\"ShareSpaceStay1.csv\" link in gallery"), "visible"));
      await session.step(146, "And \"ShareRecipient2\" link in gallery should be absent", () => shouldBe(page, el("\"ShareRecipient2\" link in gallery"), "absent"));
    });
    await run.scenario("The owner takes the share back", async () => {
      await session.step(149, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(150, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(151, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(152, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(153, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(154, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
      await session.step(155, "When user clicks on Platform---Functions---Queries tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Queries tree node inside browse tree")));
      await session.step(156, "Then the \"Queries\" view should be current", () => viewIsCurrent(page, "Queries"));
      await session.step(157, "When user enters \"ShareRecipient1\" into gallery search", () => enterInto(page, "ShareRecipient1", el("gallery search")));
      await session.step(158, "And user picks \"Share...\" from the context menu of \"ShareRecipient1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareRecipient1\" gallery card")));
      await session.step(159, "Then \"Share ShareRecipient1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareRecipient1\" dialog"), "visible"));
      await session.step(160, "And the access level of the sharing user in \"Share ShareRecipient1\" dialog should be \"View and use\"", () => sharingUserLevel(page, el("\"Share ShareRecipient1\" dialog"), "View and use"));
      await session.step(161, "When user removes the sharing user from \"Share ShareRecipient1\" dialog", () => removeSharingUser(page, el("\"Share ShareRecipient1\" dialog")));
      await session.step(162, "And user clicks on OK button in \"Share ShareRecipient1\" dialog", () => clickOn(page, el("OK button in \"Share ShareRecipient1\" dialog")));
      await session.step(163, "Then the \"Share ShareRecipient1\" dialog should close", () => dialogCloses(page, "Share ShareRecipient1"));
      await session.step(164, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(165, "When user clicks on \"ShareRecipient1\" gallery card", () => clickOn(page, el("\"ShareRecipient1\" gallery card")));
      await session.step(166, "Then the context panel should show \"ShareRecipient1\"", () => contextPanelShows(page, "ShareRecipient1"));
      await session.step(167, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(168, "And Sharing section in context panel should contain text \"You are the owner\"", () => shouldContainText(page, el("Sharing section in context panel"), "You are the owner"));
      await session.step(169, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(170, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
    });
    await run.scenario("The second account no longer has the query", async () => {
      await session.step(173, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(174, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(175, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(176, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(177, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(178, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(179, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(180, "Then \"ShareRecipient2\" shared by the running account should be visible", () => sharedShouldBe(page, "ShareRecipient2", "visible"));
      await session.step(181, "And \"ShareRecipient1\" shared by the running account should be absent", () => sharedShouldBe(page, "ShareRecipient1", "absent"));
      await session.step(182, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(183, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    run.finish();
  });
});
