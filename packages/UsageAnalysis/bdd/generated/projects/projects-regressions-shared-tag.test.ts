/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-shared-tag.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, sharing.share-dialog, core#296]
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
import {errorBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a tag added to a project shared for viewing", () => {
  const session = feature(test, "features/projects/projects-regressions-shared-tag.feature", import.meta.url);
  test("The second account adding a tag to a project shared for viewing is told it cannot", {tag: ["@realizes:views.projects", "@realizes:sharing.share-dialog", "@known-failure:54", "@realizes:core#296"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(17, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And no project named \"BDDSharedTag{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSharedTag{time}")));
    await knownFailure(async () => {
      await session.step(24, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(25, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(26, "When user double-clicks Files---Demo---demog-1000.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog-1000.csv tree node inside browse tree")));
      await session.step(27, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
      await session.step(28, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(29, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(30, "When user enters \"BDDSharedTag{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDSharedTag{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(31, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(32, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(33, "And \"Share BDDSharedTag{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDSharedTag{time}\" dialog")), "visible"));
      await session.step(34, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(35, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDSharedTag{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDSharedTag{time}\" dialog"))));
      await session.step(36, "And user unchecks \"Send notifications\" input in \"Share BDDSharedTag{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDSharedTag{time}\" dialog"))));
      await session.step(37, "And user clicks on OK button in \"Share BDDSharedTag{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDSharedTag{time}\" dialog"))));
      await session.step(38, "Then the \"Share BDDSharedTag{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDSharedTag{time}")));
      await session.step(39, "When user closes all views", () => closeAllViews(page));
      await session.step(40, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(41, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(42, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(43, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(44, "And user enters \"BDDSharedTag{time}\" into gallery search", () => enterInto(page, session.text("BDDSharedTag{time}"), el("gallery search")));
      await session.step(45, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(46, "And user clicks on BDDSharedTag{time} gallery card", () => clickOn(page, el(session.text("BDDSharedTag{time} gallery card"))));
      await session.step(47, "Then the context panel should show \"BDDSharedTag{time}\"", () => contextPanelShows(page, session.text("BDDSharedTag{time}")));
      await session.step(48, "When user hovers over Tags details row in context panel", () => hoverOver(page, el("Tags details row in context panel")));
      await session.step(49, "And user clicks on edit icon of Tags details row in context panel", () => clickOn(page, el("edit icon of Tags details row in context panel")));
      await session.step(50, "And user clicks on value of Tags details row in context panel", () => clickOn(page, el("value of Tags details row in context panel")));
      await session.step(51, "And user types \"bddsharedtag{time}\" at the caret", () => typeAtCaret(page, session.text("bddsharedtag{time}")));
      await session.step(52, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(53, "And user presses Escape", () => pressKey(page, "Escape"));
      await session.step(54, "Then an error balloon containing \"privileges\" should have been shown", () => errorBalloonText(page, "privileges"));
      await session.step(55, "And the project named \"BDDSharedTag{time}\" should not be tagged \"bddsharedtag{time}\" on the server", () => projectNotTaggedOnServer(page, session.text("BDDSharedTag{time}"), session.text("bddsharedtag{time}")));
      await session.step(56, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(57, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    }, 54);
  });
});
