/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/home/home-as-sharing-user.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.view.welcome, powerpack.dashboard.spotlight]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {everyWidgetHasContent, homeWidgetsAre, widgetsAllShown} from '../../bindings/home.js';
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, shouldBe, shouldContainText, shouldHaveText, tabShowing} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses, noProjectOnServer, openDataset, pickSharingUser, refreshBrowse, saveAsProject, sharingUserNoNotifications, signInAsSelf, signInAsSharingUser, signedInNotMember, unreadNotifications} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Home page of a user who is neither a developer nor an administrator", () => {
  const session = feature(test, "features/home/home-as-sharing-user.feature", import.meta.url);
  test("The Home page of a user who is neither a developer nor an administrator", {tag: ["@journey", "@serial", "@realizes:powerpack.view.welcome", "@realizes:powerpack.dashboard.spotlight"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And every widget of the Home page is stored as shown, now and when the feature ends", () => widgetsAllShown(page));
    await run.scenario("A project shared with the sharing user, with notifications on", async () => {
      await session.step(34, "Given the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
      await session.step(35, "And user opens demog dataset", () => openDataset(page, ds("demog")));
      await session.step(36, "And no project named \"bdd-home-shared-{time}\" is on the server", () => noProjectOnServer(page, session.text("bdd-home-shared-{time}")));
      await session.step(37, "And user saves the current view as project \"bdd-home-shared-{time}\"", () => saveAsProject(page, session.text("bdd-home-shared-{time}")));
      await session.step(38, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(39, "And user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(40, "And \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(41, "And \"My stuff > My dashboards\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > My dashboards\" tree node inside browse tree")));
      await session.step(42, "When user picks \"Share...\" from the context menu of \"My stuff > My dashboards > bdd-home-shared-{time}\" tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el(session.text("\"My stuff > My dashboards > bdd-home-shared-{time}\" tree node inside browse tree"))));
      await session.step(43, "And user picks the sharing user in \"User, group, or email\" input in \"Share bdd-home-shared-{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share bdd-home-shared-{time}\" dialog"))));
      await session.step(44, "Then \"Send notifications\" input in \"Share bdd-home-shared-{time}\" dialog should be checked", () => shouldBe(page, el(session.text("\"Send notifications\" input in \"Share bdd-home-shared-{time}\" dialog")), "checked"));
      await session.step(45, "When user clicks on OK button in \"Share bdd-home-shared-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share bdd-home-shared-{time}\" dialog"))));
      await session.step(46, "Then the \"Share bdd-home-shared-{time}\" dialog should close", () => dialogCloses(page, session.text("Share bdd-home-shared-{time}")));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The sharing user's Home page has Spotlight and Community only", async () => {
      await session.step(50, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(51, "Then the signed-in user should not be a member of \"Developers\"", () => signedInNotMember(page, "Developers"));
      await session.step(52, "And the signed-in user should not be a member of \"Administrators\"", () => signedInNotMember(page, "Administrators"));
      await session.step(53, "And the Home page should show the widgets \"Spotlight, Community\"", () => homeWidgetsAre(page, "Spotlight, Community"));
      await session.step(54, "And every widget of the Home page should show content", () => everyWidgetHasContent(page));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
      await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The share arrives as the one unread notification, on the server and in Spotlight", async () => {
      await session.step(59, "Then the signed-in user should have 1 unread notification on the server", () => unreadNotifications(page, 1));
      await session.step(60, "When user reloads the page", () => reloadPage(page));
      await session.step(61, "Then badge of Spotlight home widget should have text \"1\"", () => shouldHaveText(page, el("badge of Spotlight home widget"), "1"));
      await session.step(62, "When user clicks on Notifications tab in Spotlight home widget", () => clickOn(page, el("Notifications tab in Spotlight home widget")));
      await session.step(63, "Then the \"Notifications\" tab of Spotlight home widget should be showing", () => tabShowing(page, "Notifications", el("Spotlight home widget")));
      await session.step(64, "And notifications page of Spotlight home widget should contain text \"bdd-home-shared-{time}\"", () => shouldContainText(page, el("notifications page of Spotlight home widget"), session.text("bdd-home-shared-{time}")));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Back in the running account, the Home page has all four widgets", async () => {
      await session.step(68, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(69, "Then the Home page should show the widgets \"Spotlight, Reports, Usage, Community\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage, Community"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
