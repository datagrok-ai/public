/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/sharing-script.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog, sharing.browse-shared-with-me]
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
import {scriptResult} from '../../bindings/scripts.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, isExpanded, shouldBe, shouldContainText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, expandShared, noScriptOnServer, pickFromSharedMenu, pickSharingUser, refreshBrowse, removeSharingUser, runningAccountSignedIn, scriptOnServer, sharedShouldBe, sharingPaneLists, sharingPaneListsNot, sharingUserLevel, sharingUserNoNotifications, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A script, shared and taken back", () => {
  const session = feature(test, "features/sharing/sharing-script.feature", import.meta.url);
  test("A script, shared and taken back", {tag: ["@journey", "@serial", "@realizes:sharing.share-dialog", "@realizes:sharing.browse-shared-with-me"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(25, "And no script named \"ShareScript1\" is on the server", () => noScriptOnServer(page, "ShareScript1"));
    await session.step(26, "And no script named \"ShareScript2\" is on the server", () => noScriptOnServer(page, "ShareScript2"));
    await session.step(27, "And a script \"ShareScript1\" is on the server:", () => scriptOnServer(page, "ShareScript1", "//language: javascript\n//output: string res\nres = 'shared-ok';"));
    await session.step(33, "And a script \"ShareScript2\" is on the server:", () => scriptOnServer(page, "ShareScript2", "//language: javascript\n//output: string res\nres = 'anchor';"));
    await session.step(39, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(40, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("The owner shares both scripts to view and use", async () => {
      await session.step(43, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(44, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
      await session.step(45, "When user clicks on Platform---Functions---Scripts tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Scripts tree node inside browse tree")));
      await session.step(46, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(47, "When user enters \"ShareScript2\" into gallery search", () => enterInto(page, "ShareScript2", el("gallery search")));
      await session.step(48, "And user picks \"Share...\" from the context menu of \"ShareScript2\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareScript2\" gallery card")));
      await session.step(49, "And user picks the sharing user in \"User, group, or email\" input in \"Share ShareScript2\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareScript2\" dialog")));
      await session.step(50, "And user unchecks \"Send notifications\" input in \"Share ShareScript2\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareScript2\" dialog")));
      await session.step(51, "And user clicks on OK button in \"Share ShareScript2\" dialog", () => clickOn(page, el("OK button in \"Share ShareScript2\" dialog")));
      await session.step(52, "Then the \"Share ShareScript2\" dialog should close", () => dialogCloses(page, "Share ShareScript2"));
      await session.step(53, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(54, "When user enters \"ShareScript1\" into gallery search", () => enterInto(page, "ShareScript1", el("gallery search")));
      await session.step(55, "And user picks \"Share...\" from the context menu of \"ShareScript1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareScript1\" gallery card")));
      await session.step(56, "Then \"Share ShareScript1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareScript1\" dialog"), "visible"));
      await session.step(57, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(58, "When user picks the sharing user in \"User, group, or email\" input in \"Share ShareScript1\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareScript1\" dialog")));
      await session.step(59, "And user unchecks \"Send notifications\" input in \"Share ShareScript1\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareScript1\" dialog")));
      await session.step(60, "And user clicks on OK button in \"Share ShareScript1\" dialog", () => clickOn(page, el("OK button in \"Share ShareScript1\" dialog")));
      await session.step(61, "Then the \"Share ShareScript1\" dialog should close", () => dialogCloses(page, "Share ShareScript1"));
      await session.step(62, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(63, "When user clicks on \"ShareScript1\" gallery card", () => clickOn(page, el("\"ShareScript1\" gallery card")));
      await session.step(64, "Then the context panel should show \"ShareScript1\"", () => contextPanelShows(page, "ShareScript1"));
      await session.step(65, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(66, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(67, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
      await session.step(69, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The second account finds the script under Shared with me and runs it", async () => {
      await session.step(72, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(73, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(74, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(75, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(76, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(77, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(78, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(79, "Then \"ShareScript1\" shared by the running account should be visible", () => sharedShouldBe(page, "ShareScript1", "visible"));
      await session.step(80, "When user picks \"Run\" from the context menu of \"ShareScript1\" shared by the running account", () => pickFromSharedMenu(page, "Run", "ShareScript1"));
      await session.step(81, "Then the script results should show \"res\" as \"shared-ok\"", () => scriptResult(page, "res", "shared-ok"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
      await session.step(83, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The owner takes the share back", async () => {
      await session.step(86, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(87, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(88, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(89, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(90, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(91, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
      await session.step(92, "When user clicks on Platform---Functions---Scripts tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Scripts tree node inside browse tree")));
      await session.step(93, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(94, "When user enters \"ShareScript1\" into gallery search", () => enterInto(page, "ShareScript1", el("gallery search")));
      await session.step(95, "And user picks \"Share...\" from the context menu of \"ShareScript1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareScript1\" gallery card")));
      await session.step(96, "Then \"Share ShareScript1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareScript1\" dialog"), "visible"));
      await session.step(97, "And the access level of the sharing user in \"Share ShareScript1\" dialog should be \"View and use\"", () => sharingUserLevel(page, el("\"Share ShareScript1\" dialog"), "View and use"));
      await session.step(98, "When user removes the sharing user from \"Share ShareScript1\" dialog", () => removeSharingUser(page, el("\"Share ShareScript1\" dialog")));
      await session.step(99, "And user clicks on OK button in \"Share ShareScript1\" dialog", () => clickOn(page, el("OK button in \"Share ShareScript1\" dialog")));
      await session.step(100, "Then the \"Share ShareScript1\" dialog should close", () => dialogCloses(page, "Share ShareScript1"));
      await session.step(101, "When user clicks on \"ShareScript1\" gallery card", () => clickOn(page, el("\"ShareScript1\" gallery card")));
      await session.step(102, "Then the context panel should show \"ShareScript1\"", () => contextPanelShows(page, "ShareScript1"));
      await session.step(103, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(104, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(105, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
    });
    await run.scenario("The second account no longer has the script", async () => {
      await session.step(108, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(109, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(110, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(111, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(112, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(113, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(114, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(115, "Then \"ShareScript2\" shared by the running account should be visible", () => sharedShouldBe(page, "ShareScript2", "visible"));
      await session.step(116, "And \"ShareScript1\" shared by the running account should be absent", () => sharedShouldBe(page, "ShareScript1", "absent"));
      await session.step(117, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(118, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    run.finish();
  });
});
