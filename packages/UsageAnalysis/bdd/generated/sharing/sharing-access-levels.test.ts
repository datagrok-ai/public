/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/sharing-access-levels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog]
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
import {check, clearField, clickOn, enterInto, isExpanded, pressKey, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, pickSharingUser, queryOnServer, sharingPaneLists, sharingPaneListsNot, sharingUserLevel, sharingUserNoNotifications, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sharing at Full access, and with All users", () => {
  const session = feature(test, "features/sharing/sharing-access-levels.feature", import.meta.url);
  test("Sharing at Full access, and with All users", {tag: ["@journey", "@serial", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(26, "And a query \"ShareLevel1\" on \"System:Datagrok\" reads \"select 1 as x\"", () => queryOnServer(page, "ShareLevel1", "System:Datagrok", "select 1 as x"));
    await session.step(27, "And a query \"ShareLevel2\" on \"System:Datagrok\" reads \"select 2 as y\"", () => queryOnServer(page, "ShareLevel2", "System:Datagrok", "select 2 as y"));
    await session.step(28, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(31, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(32, "When user clicks on Platform---Functions---Queries tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Queries tree node inside browse tree")));
    await session.step(33, "Then the \"Queries\" view should be current", () => viewIsCurrent(page, "Queries"));
    await session.step(34, "When user enters \"ShareLevel\" into gallery search", () => enterInto(page, "ShareLevel", el("gallery search")));
    await session.step(35, "Then \"ShareLevel1\" gallery card should be visible", () => shouldBe(page, el("\"ShareLevel1\" gallery card"), "visible"));
    await session.step(36, "And \"ShareLevel2\" gallery card should be visible", () => shouldBe(page, el("\"ShareLevel2\" gallery card"), "visible"));
    await run.scenario("A level checked in the list is the level the grant gets", async () => {
      await session.step(39, "When user picks \"Share...\" from the context menu of \"ShareLevel1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareLevel1\" gallery card")));
      await session.step(40, "Then \"Share ShareLevel1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareLevel1\" dialog"), "visible"));
      await session.step(41, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(42, "When user picks the sharing user in \"User, group, or email\" input in \"Share ShareLevel1\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareLevel1\" dialog")));
      await session.step(43, "And user clicks on share access selector", () => clickOn(page, el("share access selector")));
      await session.step(44, "And user checks \"Full access\" tree node inside privilege tree", () => check(page, el("\"Full access\" tree node inside privilege tree")));
      await session.step(45, "And user presses Escape", () => pressKey(page, "Escape"));
      await session.step(46, "Then privilege tree should be hidden", () => shouldBe(page, el("privilege tree"), "hidden"));
      await session.step(47, "And share access selector should contain text \"Full access\"", () => shouldContainText(page, el("share access selector"), "Full access"));
      await session.step(48, "When user unchecks \"Send notifications\" input in \"Share ShareLevel1\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareLevel1\" dialog")));
      await session.step(49, "And user clicks on OK button in \"Share ShareLevel1\" dialog", () => clickOn(page, el("OK button in \"Share ShareLevel1\" dialog")));
      await session.step(50, "Then the \"Share ShareLevel1\" dialog should close", () => dialogCloses(page, "Share ShareLevel1"));
      await session.step(51, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(52, "When user clicks on \"ShareLevel2\" gallery card", () => clickOn(page, el("\"ShareLevel2\" gallery card")));
      await session.step(53, "Then the context panel should show \"ShareLevel2\"", () => contextPanelShows(page, "ShareLevel2"));
      await session.step(54, "When user clicks on \"ShareLevel1\" gallery card", () => clickOn(page, el("\"ShareLevel1\" gallery card")));
      await session.step(55, "Then the context panel should show \"ShareLevel1\"", () => contextPanelShows(page, "ShareLevel1"));
      await session.step(56, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(57, "And Sharing section in context panel should contain text \"has full access\"", () => shouldContainText(page, el("Sharing section in context panel"), "has full access"));
      await session.step(58, "When user picks \"Share...\" from the context menu of \"ShareLevel1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareLevel1\" gallery card")));
      await session.step(59, "Then \"Share ShareLevel1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareLevel1\" dialog"), "visible"));
      await session.step(60, "And the access level of the sharing user in \"Share ShareLevel1\" dialog should be \"Full access\"", () => sharingUserLevel(page, el("\"Share ShareLevel1\" dialog"), "Full access"));
      await session.step(61, "When user clicks on CANCEL button in \"Share ShareLevel1\" dialog", () => clickOn(page, el("CANCEL button in \"Share ShareLevel1\" dialog")));
      await session.step(62, "Then the \"Share ShareLevel1\" dialog should close", () => dialogCloses(page, "Share ShareLevel1"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
      await session.step(64, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("All users is shared with as a user is", async () => {
      await session.step(67, "When user clicks on \"ShareLevel1\" gallery card", () => clickOn(page, el("\"ShareLevel1\" gallery card")));
      await session.step(68, "Then the context panel should show \"ShareLevel1\"", () => contextPanelShows(page, "ShareLevel1"));
      await session.step(69, "When user clicks on \"ShareLevel2\" gallery card", () => clickOn(page, el("\"ShareLevel2\" gallery card")));
      await session.step(70, "Then the context panel should show \"ShareLevel2\"", () => contextPanelShows(page, "ShareLevel2"));
      await session.step(71, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(72, "And Sharing section in context panel should not contain text \"All users\"", () => shouldNotContainText(page, el("Sharing section in context panel"), "All users"));
      await session.step(73, "When user picks \"Share...\" from the context menu of \"ShareLevel2\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareLevel2\" gallery card")));
      await session.step(74, "Then \"Share ShareLevel2\" dialog should be visible", () => shouldBe(page, el("\"Share ShareLevel2\" dialog"), "visible"));
      await session.step(75, "When user types \"All\" into \"User, group, or email\" input in \"Share ShareLevel2\" dialog", () => typeInto(page, "All", el("\"User, group, or email\" input in \"Share ShareLevel2\" dialog")));
      await session.step(76, "And user clicks on \"All users\" table row in \"Share ShareLevel2\" dialog", () => clickOn(page, el("\"All users\" table row in \"Share ShareLevel2\" dialog")));
      await session.step(77, "And user unchecks \"Send notifications\" input in \"Share ShareLevel2\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareLevel2\" dialog")));
      await session.step(78, "And user clicks on OK button in \"Share ShareLevel2\" dialog", () => clickOn(page, el("OK button in \"Share ShareLevel2\" dialog")));
      await session.step(79, "Then the \"Share ShareLevel2\" dialog should close", () => dialogCloses(page, "Share ShareLevel2"));
      await session.step(80, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(81, "When user clicks on \"ShareLevel1\" gallery card", () => clickOn(page, el("\"ShareLevel1\" gallery card")));
      await session.step(82, "Then the context panel should show \"ShareLevel1\"", () => contextPanelShows(page, "ShareLevel1"));
      await session.step(83, "When user clicks on \"ShareLevel2\" gallery card", () => clickOn(page, el("\"ShareLevel2\" gallery card")));
      await session.step(84, "Then the context panel should show \"ShareLevel2\"", () => contextPanelShows(page, "ShareLevel2"));
      await session.step(85, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(86, "And Sharing section in context panel should contain text \"All users\"", () => shouldContainText(page, el("Sharing section in context panel"), "All users"));
      await session.step(87, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(88, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
      await session.step(90, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
