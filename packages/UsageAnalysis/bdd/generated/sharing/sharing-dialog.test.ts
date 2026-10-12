/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/sharing-dialog.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, followingShouldBe, isExpanded, pressKey, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, pickSharingUser, queryOnServer, sharingPaneListsNot, shellSettingIs, suggestionsNotOfferOwn, typeOwnLogin, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Share dialog and the Sharing pane", () => {
  const session = feature(test, "features/sharing/sharing-dialog.feature", import.meta.url);
  test("The Share dialog and the Sharing pane", {tag: ["@journey", "@serial", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And a query \"ShareDialog1\" on \"System:Datagrok\" reads \"select 1 as x\"", () => queryOnServer(page, "ShareDialog1", "System:Datagrok", "select 1 as x"));
    await session.step(31, "And a query \"ShareDialog2\" on \"System:Datagrok\" reads \"select 2 as y\"", () => queryOnServer(page, "ShareDialog2", "System:Datagrok", "select 2 as y"));
    await session.step(32, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(34, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(35, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(36, "When user clicks on Platform---Functions---Queries tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Queries tree node inside browse tree")));
    await session.step(37, "Then the \"Queries\" view should be current", () => viewIsCurrent(page, "Queries"));
    await session.step(38, "When user enters \"ShareDialog\" into gallery search", () => enterInto(page, "ShareDialog", el("gallery search")));
    await session.step(39, "Then \"ShareDialog1\" gallery card should be visible", () => shouldBe(page, el("\"ShareDialog1\" gallery card"), "visible"));
    await session.step(40, "And \"ShareDialog2\" gallery card should be visible", () => shouldBe(page, el("\"ShareDialog2\" gallery card"), "visible"));
    await run.scenario("The Sharing pane names the owner and opens the dialog", async () => {
      await session.step(43, "When user clicks on \"ShareDialog1\" gallery card", () => clickOn(page, el("\"ShareDialog1\" gallery card")));
      await session.step(44, "Then the context panel should show \"ShareDialog1\"", () => contextPanelShows(page, "ShareDialog1"));
      await session.step(45, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(46, "And Sharing section in context panel should contain text \"You are the owner\"", () => shouldContainText(page, el("Sharing section in context panel"), "You are the owner"));
      await session.step(47, "And \"Share...\" button in Sharing section in context panel should be visible", () => shouldBe(page, el("\"Share...\" button in Sharing section in context panel"), "visible"));
      await session.step(48, "And \"Calculate resulting permissions for this entity\" button in Sharing section in context panel should be visible", () => shouldBe(page, el("\"Calculate resulting permissions for this entity\" button in Sharing section in context panel"), "visible"));
      await session.step(49, "When user clicks on \"Share...\" button in Sharing section in context panel", () => clickOn(page, el("\"Share...\" button in Sharing section in context panel")));
      await session.step(50, "Then \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareDialog1\" dialog"), "visible"));
      await session.step(51, "When user clicks on CANCEL button in \"Share ShareDialog1\" dialog", () => clickOn(page, el("CANCEL button in \"Share ShareDialog1\" dialog")));
      await session.step(52, "Then the \"Share ShareDialog1\" dialog should close", () => dialogCloses(page, "Share ShareDialog1"));
      await session.step(53, "When user clicks on \"ShareDialog2\" gallery card", () => clickOn(page, el("\"ShareDialog2\" gallery card")));
      await session.step(54, "Then the context panel should show \"ShareDialog2\"", () => contextPanelShows(page, "ShareDialog2"));
      await session.step(55, "When user clicks on \"ShareDialog1\" gallery card", () => clickOn(page, el("\"ShareDialog1\" gallery card")));
      await session.step(56, "Then the context panel should show \"ShareDialog1\"", () => contextPanelShows(page, "ShareDialog1"));
      await session.step(57, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(58, "And Sharing section in context panel should contain text \"You are the owner\"", () => shouldContainText(page, el("Sharing section in context panel"), "You are the owner"));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The dialog asks who and how much, and offers a message once someone is picked", async () => {
      await session.step(62, "When user picks \"Share...\" from the context menu of \"ShareDialog1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareDialog1\" gallery card")));
      await session.step(63, "Then \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareDialog1\" dialog"), "visible"));
      await session.step(64, "And \"User, group, or email\" input in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"User, group, or email\" input in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(65, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(66, "And OK button in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("OK button in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(67, "And CANCEL button in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("CANCEL button in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(68, "And \"Send notifications\" input in \"Share ShareDialog1\" dialog should be hidden", () => shouldBe(page, el("\"Send notifications\" input in \"Share ShareDialog1\" dialog"), "hidden"));
      await session.step(69, "And share message should be hidden", () => shouldBe(page, el("share message"), "hidden"));
      await session.step(70, "And \"Remove\" icon in \"Share ShareDialog1\" dialog should be hidden", () => shouldBe(page, el("\"Remove\" icon in \"Share ShareDialog1\" dialog"), "hidden"));
      await session.step(71, "When user picks the sharing user in \"User, group, or email\" input in \"Share ShareDialog1\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareDialog1\" dialog")));
      await session.step(72, "Then \"Remove\" icon in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Remove\" icon in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(73, "And \"Send notifications\" input in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Send notifications\" input in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(74, "And share message should be visible", () => shouldBe(page, el("share message"), "visible"));
      await session.step(75, "And the \"sendNotificationOnShareCheckbox\" shell setting should be true", () => shellSettingIs(page, "sendNotificationOnShareCheckbox", "true"));
      await session.step(76, "And \"Send notifications\" input in \"Share ShareDialog1\" dialog should be checked", () => shouldBe(page, el("\"Send notifications\" input in \"Share ShareDialog1\" dialog"), "checked"));
      await session.step(77, "And \"Share ShareDialog1\" dialog should contain text \"(Will be sent to 1 user)\"", () => shouldContainText(page, el("\"Share ShareDialog1\" dialog"), "(Will be sent to 1 user)"));
      await session.step(78, "When user clicks on CANCEL button in \"Share ShareDialog1\" dialog", () => clickOn(page, el("CANCEL button in \"Share ShareDialog1\" dialog")));
      await session.step(79, "Then the \"Share ShareDialog1\" dialog should close", () => dialogCloses(page, "Share ShareDialog1"));
      await session.step(80, "When user clicks on \"ShareDialog2\" gallery card", () => clickOn(page, el("\"ShareDialog2\" gallery card")));
      await session.step(81, "Then the context panel should show \"ShareDialog2\"", () => contextPanelShows(page, "ShareDialog2"));
      await session.step(82, "When user clicks on \"ShareDialog1\" gallery card", () => clickOn(page, el("\"ShareDialog1\" gallery card")));
      await session.step(83, "Then the context panel should show \"ShareDialog1\"", () => contextPanelShows(page, "ShareDialog1"));
      await session.step(84, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Escape closes the list of levels and leaves the dialog open, and CANCEL closes it while the suggestions are open (GROK-20758)", async () => {
      await session.step(88, "When user picks \"Share...\" from the context menu of \"ShareDialog1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareDialog1\" gallery card")));
      await session.step(89, "Then \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareDialog1\" dialog"), "visible"));
      await session.step(90, "When user clicks on share access selector", () => clickOn(page, el("share access selector")));
      await session.step(91, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Full access\" tree node inside privilege tree"],["\"View and use\" tree node inside privilege tree"],["\"Edit\" tree node inside privilege tree"],["\"Delete\" tree node inside privilege tree"],["\"Share\" tree node inside privilege tree"]]), [["\"Full access\" tree node inside privilege tree"],["\"View and use\" tree node inside privilege tree"],["\"Edit\" tree node inside privilege tree"],["\"Delete\" tree node inside privilege tree"],["\"Share\" tree node inside privilege tree"]]);
      await session.step(97, "And \"View and use\" tree node inside privilege tree should be checked", () => shouldBe(page, el("\"View and use\" tree node inside privilege tree"), "checked"));
      await session.step(98, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(99, "Then privilege tree should be hidden", () => shouldBe(page, el("privilege tree"), "hidden"));
      await session.step(100, "And \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareDialog1\" dialog"), "visible"));
      await session.step(101, "And \"All users\" table row in \"Share ShareDialog1\" dialog should be hidden", () => shouldBe(page, el("\"All users\" table row in \"Share ShareDialog1\" dialog"), "hidden"));
      await session.step(102, "When user types \"All\" into \"User, group, or email\" input in \"Share ShareDialog1\" dialog", () => typeInto(page, "All", el("\"User, group, or email\" input in \"Share ShareDialog1\" dialog")));
      await session.step(103, "Then \"All users\" table row in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"All users\" table row in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(104, "When user clicks on CANCEL button in \"Share ShareDialog1\" dialog", () => clickOn(page, el("CANCEL button in \"Share ShareDialog1\" dialog")));
      await session.step(105, "Then the \"Share ShareDialog1\" dialog should close", () => dialogCloses(page, "Share ShareDialog1"));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Your own account is not offered (GROK-20755)", async () => {
      await session.step(109, "When user picks \"Share...\" from the context menu of \"ShareDialog1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareDialog1\" gallery card")));
      await session.step(110, "Then \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareDialog1\" dialog"), "visible"));
      await session.step(111, "When user types the current user's login into \"User, group, or email\" input in \"Share ShareDialog1\" dialog", () => typeOwnLogin(page, el("\"User, group, or email\" input in \"Share ShareDialog1\" dialog")));
      await session.step(112, "Then the user suggestions should not offer the current user", () => suggestionsNotOfferOwn(page));
      await session.step(113, "When user clicks on CANCEL button in \"Share ShareDialog1\" dialog", () => clickOn(page, el("CANCEL button in \"Share ShareDialog1\" dialog")));
      await session.step(114, "Then the \"Share ShareDialog1\" dialog should close", () => dialogCloses(page, "Share ShareDialog1"));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Resulting permissions shows what a user gets", async () => {
      await session.step(118, "When user clicks on \"ShareDialog1\" gallery card", () => clickOn(page, el("\"ShareDialog1\" gallery card")));
      await session.step(119, "Then the context panel should show \"ShareDialog1\"", () => contextPanelShows(page, "ShareDialog1"));
      await session.step(120, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(121, "When user clicks on \"Calculate resulting permissions for this entity\" button in Sharing section in context panel", () => clickOn(page, el("\"Calculate resulting permissions for this entity\" button in Sharing section in context panel")));
      await session.step(122, "Then \"Resulting permissions\" dialog should be visible", () => shouldBe(page, el("\"Resulting permissions\" dialog"), "visible"));
      await session.step(123, "And \"Resulting permissions\" dialog should contain text \"ShareDialog1\"", () => shouldContainText(page, el("\"Resulting permissions\" dialog"), "ShareDialog1"));
      await session.step(124, "And \"User, group or role\" input in \"Resulting permissions\" dialog should be visible", () => shouldBe(page, el("\"User, group or role\" input in \"Resulting permissions\" dialog"), "visible"));
      await session.step(125, "And CANCEL button in \"Resulting permissions\" dialog should be visible", () => shouldBe(page, el("CANCEL button in \"Resulting permissions\" dialog"), "visible"));
      await session.step(126, "And \"Resulting permissions\" dialog should not contain text \"No access\"", () => shouldNotContainText(page, el("\"Resulting permissions\" dialog"), "No access"));
      await session.step(127, "When user picks the sharing user in \"User, group or role\" input in \"Resulting permissions\" dialog", () => pickSharingUser(page, el("\"User, group or role\" input in \"Resulting permissions\" dialog")));
      await session.step(128, "Then \"Resulting permissions\" dialog should contain text \"No access\"", () => shouldContainText(page, el("\"Resulting permissions\" dialog"), "No access"));
      await session.step(129, "When user clicks on CANCEL button in \"Resulting permissions\" dialog", () => clickOn(page, el("CANCEL button in \"Resulting permissions\" dialog")));
      await session.step(130, "Then the \"Resulting permissions\" dialog should close", () => dialogCloses(page, "Resulting permissions"));
      await session.step(131, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(132, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
      await session.step(133, "And no errors should have been logged", () => noErrors(page));
      await session.step(134, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
