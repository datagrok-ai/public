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
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, pickSharingUser, queryOnServer, sharingPaneListsNot, shellSettingIs, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Share dialog and the Sharing pane", () => {
  const session = feature(test, "features/sharing/sharing-dialog.feature", import.meta.url);
  test("The Share dialog and the Sharing pane", {tag: ["@journey", "@serial", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And a query \"ShareDialog1\" on \"System:Datagrok\" reads \"select 1 as x\"", () => queryOnServer(page, "ShareDialog1", "System:Datagrok", "select 1 as x"));
    await session.step(34, "And a query \"ShareDialog2\" on \"System:Datagrok\" reads \"select 2 as y\"", () => queryOnServer(page, "ShareDialog2", "System:Datagrok", "select 2 as y"));
    await session.step(35, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(36, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(37, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(38, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(39, "When user clicks on Platform---Functions---Queries tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Queries tree node inside browse tree")));
    await session.step(40, "Then the \"Queries\" view should be current", () => viewIsCurrent(page, "Queries"));
    await session.step(41, "When user enters \"ShareDialog\" into gallery search", () => enterInto(page, "ShareDialog", el("gallery search")));
    await session.step(42, "Then \"ShareDialog1\" gallery card should be visible", () => shouldBe(page, el("\"ShareDialog1\" gallery card"), "visible"));
    await session.step(43, "And \"ShareDialog2\" gallery card should be visible", () => shouldBe(page, el("\"ShareDialog2\" gallery card"), "visible"));
    await run.scenario("The Sharing pane names the owner and opens the dialog", async () => {
      await session.step(46, "When user clicks on \"ShareDialog1\" gallery card", () => clickOn(page, el("\"ShareDialog1\" gallery card")));
      await session.step(47, "Then the context panel should show \"ShareDialog1\"", () => contextPanelShows(page, "ShareDialog1"));
      await session.step(48, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(49, "And Sharing section in context panel should contain text \"You are the owner\"", () => shouldContainText(page, el("Sharing section in context panel"), "You are the owner"));
      await session.step(50, "And \"Share...\" button in Sharing section in context panel should be visible", () => shouldBe(page, el("\"Share...\" button in Sharing section in context panel"), "visible"));
      await session.step(51, "And \"Calculate resulting permissions for this entity\" button in Sharing section in context panel should be visible", () => shouldBe(page, el("\"Calculate resulting permissions for this entity\" button in Sharing section in context panel"), "visible"));
      await session.step(52, "When user clicks on \"Share...\" button in Sharing section in context panel", () => clickOn(page, el("\"Share...\" button in Sharing section in context panel")));
      await session.step(53, "Then \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareDialog1\" dialog"), "visible"));
      await session.step(54, "When user clicks on CANCEL button in \"Share ShareDialog1\" dialog", () => clickOn(page, el("CANCEL button in \"Share ShareDialog1\" dialog")));
      await session.step(55, "Then the \"Share ShareDialog1\" dialog should close", () => dialogCloses(page, "Share ShareDialog1"));
      await session.step(56, "When user clicks on \"ShareDialog2\" gallery card", () => clickOn(page, el("\"ShareDialog2\" gallery card")));
      await session.step(57, "Then the context panel should show \"ShareDialog2\"", () => contextPanelShows(page, "ShareDialog2"));
      await session.step(58, "When user clicks on \"ShareDialog1\" gallery card", () => clickOn(page, el("\"ShareDialog1\" gallery card")));
      await session.step(59, "Then the context panel should show \"ShareDialog1\"", () => contextPanelShows(page, "ShareDialog1"));
      await session.step(60, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(61, "And Sharing section in context panel should contain text \"You are the owner\"", () => shouldContainText(page, el("Sharing section in context panel"), "You are the owner"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The dialog asks who and how much, and offers a message once someone is picked", async () => {
      await session.step(65, "When user picks \"Share...\" from the context menu of \"ShareDialog1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareDialog1\" gallery card")));
      await session.step(66, "Then \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareDialog1\" dialog"), "visible"));
      await session.step(67, "And \"User, group, or email\" input in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"User, group, or email\" input in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(68, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(69, "And OK button in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("OK button in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(70, "And CANCEL button in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("CANCEL button in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(71, "And \"Send notifications\" input in \"Share ShareDialog1\" dialog should be hidden", () => shouldBe(page, el("\"Send notifications\" input in \"Share ShareDialog1\" dialog"), "hidden"));
      await session.step(72, "And \"Remove\" icon in \"Share ShareDialog1\" dialog should be hidden", () => shouldBe(page, el("\"Remove\" icon in \"Share ShareDialog1\" dialog"), "hidden"));
      await session.step(73, "When user picks the sharing user in \"User, group, or email\" input in \"Share ShareDialog1\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareDialog1\" dialog")));
      await session.step(74, "Then \"Remove\" icon in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Remove\" icon in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(75, "And \"Send notifications\" input in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Send notifications\" input in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(76, "And the \"sendNotificationOnShareCheckbox\" shell setting should be true", () => shellSettingIs(page, "sendNotificationOnShareCheckbox", "true"));
      await session.step(77, "And \"Send notifications\" input in \"Share ShareDialog1\" dialog should be checked", () => shouldBe(page, el("\"Send notifications\" input in \"Share ShareDialog1\" dialog"), "checked"));
      await session.step(78, "And \"Share ShareDialog1\" dialog should contain text \"(Will be sent to 1 user)\"", () => shouldContainText(page, el("\"Share ShareDialog1\" dialog"), "(Will be sent to 1 user)"));
      await session.step(79, "When user clicks on CANCEL button in \"Share ShareDialog1\" dialog", () => clickOn(page, el("CANCEL button in \"Share ShareDialog1\" dialog")));
      await session.step(80, "Then the \"Share ShareDialog1\" dialog should close", () => dialogCloses(page, "Share ShareDialog1"));
      await session.step(81, "When user clicks on \"ShareDialog2\" gallery card", () => clickOn(page, el("\"ShareDialog2\" gallery card")));
      await session.step(82, "Then the context panel should show \"ShareDialog2\"", () => contextPanelShows(page, "ShareDialog2"));
      await session.step(83, "When user clicks on \"ShareDialog1\" gallery card", () => clickOn(page, el("\"ShareDialog1\" gallery card")));
      await session.step(84, "Then the context panel should show \"ShareDialog1\"", () => contextPanelShows(page, "ShareDialog1"));
      await session.step(85, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Escape closes the list of levels and leaves the dialog open, and CANCEL closes it while the suggestions are open (GROK-20758)", async () => {
      await session.step(89, "When user picks \"Share...\" from the context menu of \"ShareDialog1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareDialog1\" gallery card")));
      await session.step(90, "Then \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareDialog1\" dialog"), "visible"));
      await session.step(91, "When user clicks on share access selector", () => clickOn(page, el("share access selector")));
      await session.step(92, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Full access\" tree node"],["\"View and use\" tree node"],["\"Edit\" tree node"],["\"Delete\" tree node"],["\"Share\" tree node"]]), [["\"Full access\" tree node"],["\"View and use\" tree node"],["\"Edit\" tree node"],["\"Delete\" tree node"],["\"Share\" tree node"]]);
      await session.step(98, "And \"View and use\" tree node should be checked", () => shouldBe(page, el("\"View and use\" tree node"), "checked"));
      await session.step(99, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(100, "Then \"Full access\" tree node should be hidden", () => shouldBe(page, el("\"Full access\" tree node"), "hidden"));
      await session.step(101, "And \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareDialog1\" dialog"), "visible"));
      await session.step(102, "And \"All users\" table row in \"Share ShareDialog1\" dialog should be hidden", () => shouldBe(page, el("\"All users\" table row in \"Share ShareDialog1\" dialog"), "hidden"));
      await session.step(103, "When user types \"All\" into \"User, group, or email\" input in \"Share ShareDialog1\" dialog", () => typeInto(page, "All", el("\"User, group, or email\" input in \"Share ShareDialog1\" dialog")));
      await session.step(104, "Then \"All users\" table row in \"Share ShareDialog1\" dialog should be visible", () => shouldBe(page, el("\"All users\" table row in \"Share ShareDialog1\" dialog"), "visible"));
      await session.step(105, "When user clicks on CANCEL button in \"Share ShareDialog1\" dialog", () => clickOn(page, el("CANCEL button in \"Share ShareDialog1\" dialog")));
      await session.step(106, "Then the \"Share ShareDialog1\" dialog should close", () => dialogCloses(page, "Share ShareDialog1"));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Resulting permissions shows what a user gets", async () => {
      await session.step(110, "When user clicks on \"ShareDialog1\" gallery card", () => clickOn(page, el("\"ShareDialog1\" gallery card")));
      await session.step(111, "Then the context panel should show \"ShareDialog1\"", () => contextPanelShows(page, "ShareDialog1"));
      await session.step(112, "And the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(113, "When user clicks on \"Calculate resulting permissions for this entity\" button in Sharing section in context panel", () => clickOn(page, el("\"Calculate resulting permissions for this entity\" button in Sharing section in context panel")));
      await session.step(114, "Then \"Resulting permissions\" dialog should be visible", () => shouldBe(page, el("\"Resulting permissions\" dialog"), "visible"));
      await session.step(115, "And \"Resulting permissions\" dialog should contain text \"ShareDialog1\"", () => shouldContainText(page, el("\"Resulting permissions\" dialog"), "ShareDialog1"));
      await session.step(116, "And \"User, group or role\" input in \"Resulting permissions\" dialog should be visible", () => shouldBe(page, el("\"User, group or role\" input in \"Resulting permissions\" dialog"), "visible"));
      await session.step(117, "And CANCEL button in \"Resulting permissions\" dialog should be visible", () => shouldBe(page, el("CANCEL button in \"Resulting permissions\" dialog"), "visible"));
      await session.step(118, "And \"Resulting permissions\" dialog should not contain text \"No access\"", () => shouldNotContainText(page, el("\"Resulting permissions\" dialog"), "No access"));
      await session.step(119, "When user picks the sharing user in \"User, group or role\" input in \"Resulting permissions\" dialog", () => pickSharingUser(page, el("\"User, group or role\" input in \"Resulting permissions\" dialog")));
      await session.step(120, "Then \"Resulting permissions\" dialog should contain text \"No access\"", () => shouldContainText(page, el("\"Resulting permissions\" dialog"), "No access"));
      await session.step(121, "When user clicks on CANCEL button in \"Resulting permissions\" dialog", () => clickOn(page, el("CANCEL button in \"Resulting permissions\" dialog")));
      await session.step(122, "Then the \"Resulting permissions\" dialog should close", () => dialogCloses(page, "Resulting permissions"));
      await session.step(123, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(124, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
      await session.step(125, "And no errors should have been logged", () => noErrors(page));
      await session.step(126, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
