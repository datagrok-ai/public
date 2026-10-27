/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/share-dialog-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog, core#980]
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
import {clearField, clickOn, enterInto, isExpanded, pressKeyIn, shouldBe, shouldHaveValue, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, pickSharingUser, queryOnServer, removeSharingUser, sharingPaneLists, sharingUserLevel, sharingUserNoNotifications, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Share dialog edge cases: a cancelled removal, OK with no recipient accepted", () => {
  const session = feature(test, "features/sharing/share-dialog-edge-cases.feature", import.meta.url);
  test("A grant removed in the dialog and then cancelled stays", {tag: ["@serial", "@edge-cases", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(26, "And a query \"ShareEdgeA{time}\" on \"System:Datagrok\" reads \"select 1 as x\"", () => queryOnServer(page, session.text("ShareEdgeA{time}"), "System:Datagrok", "select 1 as x"));
    await session.step(27, "And a query \"ShareEdgeB{time}\" on \"System:Datagrok\" reads \"select 2 as y\"", () => queryOnServer(page, session.text("ShareEdgeB{time}"), "System:Datagrok", "select 2 as y"));
    await session.step(28, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(31, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(32, "When user clicks on Platform---Functions---Queries tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Queries tree node inside browse tree")));
    await session.step(33, "Then the \"Queries\" view should be current", () => viewIsCurrent(page, "Queries"));
    await session.step(34, "When user enters \"ShareEdge\" into gallery search", () => enterInto(page, "ShareEdge", el("gallery search")));
    await session.step(35, "Then \"ShareEdgeA{time}\" gallery card should be visible", () => shouldBe(page, el(session.text("\"ShareEdgeA{time}\" gallery card")), "visible"));
    await session.step(36, "And \"ShareEdgeB{time}\" gallery card should be visible", () => shouldBe(page, el(session.text("\"ShareEdgeB{time}\" gallery card")), "visible"));
    await session.step(39, "When user picks \"Share...\" from the context menu of \"ShareEdgeB{time}\" gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("\"ShareEdgeB{time}\" gallery card"))));
    await session.step(40, "And user picks the sharing user in \"User, group, or email\" input in \"Share ShareEdgeB{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share ShareEdgeB{time}\" dialog"))));
    await session.step(41, "And user unchecks \"Send notifications\" input in \"Share ShareEdgeB{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share ShareEdgeB{time}\" dialog"))));
    await session.step(42, "And user clicks on OK button in \"Share ShareEdgeB{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share ShareEdgeB{time}\" dialog"))));
    await session.step(43, "Then the \"Share ShareEdgeB{time}\" dialog should close", () => dialogCloses(page, session.text("Share ShareEdgeB{time}")));
    await session.step(44, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
    await session.step(45, "When user picks \"Share...\" from the context menu of \"ShareEdgeB{time}\" gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("\"ShareEdgeB{time}\" gallery card"))));
    await session.step(46, "Then \"Share ShareEdgeB{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share ShareEdgeB{time}\" dialog")), "visible"));
    await session.step(47, "And the access level of the sharing user in \"Share ShareEdgeB{time}\" dialog should be \"View and use\"", () => sharingUserLevel(page, el(session.text("\"Share ShareEdgeB{time}\" dialog")), "View and use"));
    await session.step(48, "When user removes the sharing user from \"Share ShareEdgeB{time}\" dialog", () => removeSharingUser(page, el(session.text("\"Share ShareEdgeB{time}\" dialog"))));
    await session.step(49, "And user clicks on CANCEL button in \"Share ShareEdgeB{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share ShareEdgeB{time}\" dialog"))));
    await session.step(50, "Then the \"Share ShareEdgeB{time}\" dialog should close", () => dialogCloses(page, session.text("Share ShareEdgeB{time}")));
    await session.step(51, "When user clicks on \"ShareEdgeA{time}\" gallery card", () => clickOn(page, el(session.text("\"ShareEdgeA{time}\" gallery card"))));
    await session.step(52, "Then the context panel should show \"ShareEdgeA{time}\"", () => contextPanelShows(page, session.text("ShareEdgeA{time}")));
    await session.step(53, "When user clicks on \"ShareEdgeB{time}\" gallery card", () => clickOn(page, el(session.text("\"ShareEdgeB{time}\" gallery card"))));
    await session.step(54, "Then the context panel should show \"ShareEdgeB{time}\"", () => contextPanelShows(page, session.text("ShareEdgeB{time}")));
    await session.step(55, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
    await session.step(56, "When user picks \"Share...\" from the context menu of \"ShareEdgeB{time}\" gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("\"ShareEdgeB{time}\" gallery card"))));
    await session.step(57, "Then the access level of the sharing user in \"Share ShareEdgeB{time}\" dialog should be \"View and use\"", () => sharingUserLevel(page, el(session.text("\"Share ShareEdgeB{time}\" dialog")), "View and use"));
    await session.step(58, "When user clicks on CANCEL button in \"Share ShareEdgeB{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share ShareEdgeB{time}\" dialog"))));
    await session.step(59, "Then the \"Share ShareEdgeB{time}\" dialog should close", () => dialogCloses(page, session.text("Share ShareEdgeB{time}")));
    await session.step(60, "When user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(61, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
    await session.step(62, "And no errors should have been logged", () => noErrors(page));
  });
  test("OK stays disabled until a recipient is accepted", {tag: ["@serial", "@edge-cases", "@realizes:sharing.share-dialog", "@known-failure:70", "@realizes:core#980"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(26, "And a query \"ShareEdgeA{time}\" on \"System:Datagrok\" reads \"select 1 as x\"", () => queryOnServer(page, session.text("ShareEdgeA{time}"), "System:Datagrok", "select 1 as x"));
    await session.step(27, "And a query \"ShareEdgeB{time}\" on \"System:Datagrok\" reads \"select 2 as y\"", () => queryOnServer(page, session.text("ShareEdgeB{time}"), "System:Datagrok", "select 2 as y"));
    await session.step(28, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(31, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(32, "When user clicks on Platform---Functions---Queries tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Queries tree node inside browse tree")));
    await session.step(33, "Then the \"Queries\" view should be current", () => viewIsCurrent(page, "Queries"));
    await session.step(34, "When user enters \"ShareEdge\" into gallery search", () => enterInto(page, "ShareEdge", el("gallery search")));
    await session.step(35, "Then \"ShareEdgeA{time}\" gallery card should be visible", () => shouldBe(page, el(session.text("\"ShareEdgeA{time}\" gallery card")), "visible"));
    await session.step(36, "And \"ShareEdgeB{time}\" gallery card should be visible", () => shouldBe(page, el(session.text("\"ShareEdgeB{time}\" gallery card")), "visible"));
    await knownFailure(async () => {
      await session.step(68, "When user picks \"Share...\" from the context menu of \"ShareEdgeA{time}\" gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("\"ShareEdgeA{time}\" gallery card"))));
      await session.step(69, "Then \"Share ShareEdgeA{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share ShareEdgeA{time}\" dialog")), "visible"));
      await session.step(70, "And OK button in \"Share ShareEdgeA{time}\" dialog should be disabled", () => shouldBe(page, el(session.text("OK button in \"Share ShareEdgeA{time}\" dialog")), "disabled"));
      await session.step(71, "When user types \"bdd-share-edge-{time}@example.com\" into \"User, group, or email\" input in \"Share ShareEdgeA{time}\" dialog", () => typeInto(page, session.text("bdd-share-edge-{time}@example.com"), el(session.text("\"User, group, or email\" input in \"Share ShareEdgeA{time}\" dialog"))));
      await session.step(72, "Then OK button in \"Share ShareEdgeA{time}\" dialog should be disabled", () => shouldBe(page, el(session.text("OK button in \"Share ShareEdgeA{time}\" dialog")), "disabled"));
      await session.step(73, "When user presses Enter in \"User, group, or email\" input in \"Share ShareEdgeA{time}\" dialog", () => pressKeyIn(page, "Enter", el(session.text("\"User, group, or email\" input in \"Share ShareEdgeA{time}\" dialog"))));
      await session.step(74, "Then \"Remove\" icon in \"Share ShareEdgeA{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Remove\" icon in \"Share ShareEdgeA{time}\" dialog")), "visible"));
      await session.step(75, "And OK button in \"Share ShareEdgeA{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"Share ShareEdgeA{time}\" dialog")), "enabled"));
      await session.step(76, "When user clicks on CANCEL button in \"Share ShareEdgeA{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share ShareEdgeA{time}\" dialog"))));
      await session.step(77, "Then the \"Share ShareEdgeA{time}\" dialog should close", () => dialogCloses(page, session.text("Share ShareEdgeA{time}")));
      await session.step(78, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(79, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
    }, 70);
  });
});
