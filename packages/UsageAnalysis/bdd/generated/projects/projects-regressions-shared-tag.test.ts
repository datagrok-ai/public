/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-shared-tag.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, sharing.share-dialog, core#296, core#1162]
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
import {clickOn, doubleClickOn, enterInto, hoverOver, isExpanded, pressKey, shouldBe, shouldContainText, typeAtCaret, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, contextPanelOpen, contextPanelShows, dialogCloses, noProjectOnServer, pickSharingUser, projectNotTaggedOnServer, runningAccountSignedIn, sharingUserNoNotifications, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, menuListsDisabled, openContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a tag added to a project shared for viewing", () => {
  const session = feature(test, "features/projects/projects-regressions-shared-tag.feature", import.meta.url);
  test("The second account adding a tag to a project shared for viewing is told it cannot", {tag: ["@realizes:views.projects", "@realizes:sharing.share-dialog", "@known-failure:59", "@realizes:core#296"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(22, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And no project named \"BDDSharedTag{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSharedTag{time}")));
    await knownFailure(async () => {
      await session.step(29, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(30, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(31, "When user double-clicks Files---Demo---demog-1000.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog-1000.csv tree node inside browse tree")));
      await session.step(32, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
      await session.step(33, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(34, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(35, "When user enters \"BDDSharedTag{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDSharedTag{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(36, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(37, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(38, "And \"Share BDDSharedTag{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDSharedTag{time}\" dialog")), "visible"));
      await session.step(39, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(40, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDSharedTag{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDSharedTag{time}\" dialog"))));
      await session.step(41, "And user unchecks \"Send notifications\" input in \"Share BDDSharedTag{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDSharedTag{time}\" dialog"))));
      await session.step(42, "And user clicks on OK button in \"Share BDDSharedTag{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDSharedTag{time}\" dialog"))));
      await session.step(43, "Then the \"Share BDDSharedTag{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDSharedTag{time}")));
      await session.step(44, "When user closes all views", () => closeAllViews(page));
      await session.step(45, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(46, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(47, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(48, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(49, "And user enters \"BDDSharedTag{time}\" into gallery search", () => enterInto(page, session.text("BDDSharedTag{time}"), el("gallery search")));
      await session.step(50, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(51, "And user clicks on BDDSharedTag{time} gallery card", () => clickOn(page, el(session.text("BDDSharedTag{time} gallery card"))));
      await session.step(52, "Then the context panel should show \"BDDSharedTag{time}\"", () => contextPanelShows(page, session.text("BDDSharedTag{time}")));
      await session.step(53, "When user hovers over Tags details row in context panel", () => hoverOver(page, el("Tags details row in context panel")));
      await session.step(54, "And user clicks on edit icon of Tags details row in context panel", () => clickOn(page, el("edit icon of Tags details row in context panel")));
      await session.step(55, "And user clicks on value of Tags details row in context panel", () => clickOn(page, el("value of Tags details row in context panel")));
      await session.step(56, "And user types \"bddsharedtag{time}\" at the caret", () => typeAtCaret(page, session.text("bddsharedtag{time}")));
      await session.step(57, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(58, "And user presses Escape", () => pressKey(page, "Escape"));
      await session.step(59, "Then an error balloon containing \"privileges\" should have been shown", () => errorBalloonText(page, "privileges"));
      await session.step(60, "And the project named \"BDDSharedTag{time}\" should not be tagged \"bddsharedtag{time}\" on the server", () => projectNotTaggedOnServer(page, session.text("BDDSharedTag{time}"), session.text("bddsharedtag{time}")));
      await session.step(61, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(62, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    }, 59);
  });
  test("A View and use recipient sees Rename disabled in the card menu", {tag: ["@realizes:views.projects", "@realizes:sharing.share-dialog", "@edge-cases", "@known-failure:93", "@realizes:core#1162"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(22, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And no project named \"BDDSharedTag{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSharedTag{time}")));
    await knownFailure(async () => {
      await session.step(68, "Given no project named \"BDDSharedMenu{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSharedMenu{time}")));
      await session.step(69, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(70, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(71, "When user double-clicks Files---Demo---cars.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---cars.csv tree node inside browse tree")));
      await session.step(72, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
      await session.step(73, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(74, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(75, "When user enters \"BDDSharedMenu{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDSharedMenu{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(76, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(77, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(78, "And \"Share BDDSharedMenu{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDSharedMenu{time}\" dialog")), "visible"));
      await session.step(79, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(80, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDSharedMenu{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDSharedMenu{time}\" dialog"))));
      await session.step(81, "And user unchecks \"Send notifications\" input in \"Share BDDSharedMenu{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDSharedMenu{time}\" dialog"))));
      await session.step(82, "And user clicks on OK button in \"Share BDDSharedMenu{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDSharedMenu{time}\" dialog"))));
      await session.step(83, "Then the \"Share BDDSharedMenu{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDSharedMenu{time}")));
      await session.step(84, "When user closes all views", () => closeAllViews(page));
      await session.step(85, "And user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(86, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(87, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(88, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(89, "And user enters \"BDDSharedMenu{time}\" into gallery search", () => enterInto(page, session.text("BDDSharedMenu{time}"), el("gallery search")));
      await session.step(90, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(91, "Then BDDSharedMenu{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDSharedMenu{time} gallery card")), "visible"));
      await session.step(92, "When user opens the context menu of BDDSharedMenu{time} gallery card", () => openContextMenu(page, el(session.text("BDDSharedMenu{time} gallery card"))));
      await session.step(93, "Then the open menu should list \"Rename...\" as disabled", () => menuListsDisabled(page, "Rename..."));
    }, 93);
  });
});
