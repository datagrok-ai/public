/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/06-dashboards-and-sharing/04-share-a-dashboard.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, hoverOver, hoverText, pressKey, shouldBe, shouldContainText, shouldHaveValue, typeInto, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {balloonsClosed, browsePanelOpen, closeAllViews, contextPanelClosed, contextPanelOpen, contextPanelShows, dialogCloses, doubleClickOpensProject, noProjectOnServer, openDataset, refreshBrowse, sharingPaneLists, sharingUserNoNotifications, sharingUserSignedIn, signInAsSharingUser, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, infoBalloonText, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Dashboards", () => {
  const session = feature(test, "features/howto/06-dashboards-and-sharing/04-share-a-dashboard.feature", import.meta.url);
  test("Share a dashboard", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And simple mode is off", () => simpleModeOff(page));
    await session.step(26, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(27, "And no project named \"Demographics overview\" is on the server", () => noProjectOnServer(page, "Demographics overview"));
    await session.step(28, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(29, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]);
    await session.step(33, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"]]), [["splitColumnName","RACE"]]);
    await session.step(35, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(36, "And user enters \"Demographics overview\" into Name text input in \"Save project\" dialog", () => enterInto(page, "Demographics overview", el("Name text input in \"Save project\" dialog")));
    await session.step(37, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(38, "And the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(39, "And \"Share Demographics overview\" dialog should be visible", () => shouldBe(page, el("\"Share Demographics overview\" dialog"), "visible"));
    await session.step(40, "And user clicks on CANCEL button in \"Share Demographics overview\" dialog", () => clickOn(page, el("CANCEL button in \"Share Demographics overview\" dialog")));
    await session.step(41, "And the \"Share Demographics overview\" dialog should close", () => dialogCloses(page, "Share Demographics overview"));
    await session.step(42, "And user closes all views", () => closeAllViews(page));
    await session.step(43, "And the notification balloons are closed", () => balloonsClosed(page));
    await session.step(44, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(45, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(46, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(47, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(48, "And user enters \"Demographics overview\" into gallery search", () => enterInto(page, "Demographics overview", el("gallery search")));
    await session.step(49, "And \"Demographics overview\" gallery card should be visible", () => shouldBe(page, el("\"Demographics overview\" gallery card"), "visible"));
    await session.step(50, "When the video says \"A new dashboard is yours alone. Share it to let others open it\"", () => videoSays(page, "A new dashboard is yours alone. Share it to let others open it"));
    await session.step(52, "And user picks \"Share...\" from the context menu of \"Demographics overview\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"Demographics overview\" gallery card")), undefined, "Right-click it, then Share...");
    await session.step(53, "Then \"Share Demographics overview\" dialog should be visible", () => shouldBe(page, el("\"Share Demographics overview\" dialog"), "visible"));
    await session.step(55, "When user types \"bddsecond\" into \"User, group, or email\" input in \"Share Demographics overview\" dialog", () => typeInto(page, "bddsecond", el("\"User, group, or email\" input in \"Share Demographics overview\" dialog")), undefined, "Type a user, a group or an email");
    await session.step(56, "And user clicks on \"bddsecond\" table row in \"Share Demographics overview\" dialog", () => clickOn(page, el("\"bddsecond\" table row in \"Share Demographics overview\" dialog")));
    await session.step(58, "And user clicks on share access selector", () => clickOn(page, el("share access selector")), undefined, "Choose what they may do: view, edit, or full access");
    await session.step(59, "Then privilege tree should be visible", () => shouldBe(page, el("privilege tree"), "visible"));
    await session.step(60, "When user hovers over \"Edit\" tree node inside privilege tree", () => hoverOver(page, el("\"Edit\" tree node inside privilege tree")));
    await session.step(61, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(62, "Then privilege tree should be hidden", () => shouldBe(page, el("privilege tree"), "hidden"));
    await session.step(63, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
    await session.step(65, "When user unchecks \"Send notifications\" input in \"Share Demographics overview\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share Demographics overview\" dialog")), undefined, "Uncheck Send notifications to share without a message");
    await session.step(66, "And user clicks on OK button in \"Share Demographics overview\" dialog", () => clickOn(page, el("OK button in \"Share Demographics overview\" dialog")));
    await session.step(67, "Then the \"Share Demographics overview\" dialog should close", () => dialogCloses(page, "Share Demographics overview"));
    await session.step(68, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
    await session.step(69, "When user clicks on \"Demographics overview\" gallery card", () => clickOn(page, el("\"Demographics overview\" gallery card")));
    await session.step(70, "Then the context panel should show \"Demographics overview\"", () => contextPanelShows(page, "Demographics overview"));
    await session.step(72, "When user hovers over \"Links...\" link in Details pane in context panel", () => hoverOver(page, el("\"Links...\" link in Details pane in context panel")), undefined, "Links... gives its address, for the people who have access");
    await session.step(74, "And user expands Sharing section in context panel", () => expand(page, el("Sharing section in context panel")), undefined, "The Sharing section shows who can open it");
    await session.step(75, "Then the sharing pane should list the sharing user", () => sharingPaneLists(page));
    await session.step(76, "When user hovers over the text \"can view and use\" in context panel", () => hoverText(page, "can view and use", el("context panel")));
    await session.step(77, "Given user closes all views", () => closeAllViews(page));
    await session.step(78, "And user signs in as the sharing user", () => signInAsSharingUser(page));
    await session.step(79, "And the sharing user should be signed in", () => sharingUserSignedIn(page));
    await session.step(80, "And simple mode is off", () => simpleModeOff(page));
    await session.step(81, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(82, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(83, "And user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(84, "And the notification balloons are closed", () => balloonsClosed(page));
    await session.step(85, "When the video says \"Now signed in as the colleague it was shared with\"", () => videoSays(page, "Now signed in as the colleague it was shared with"));
    await session.step(87, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")), undefined, "They find it in their own Browse > Dashboards");
    await session.step(88, "And user enters \"Demographics overview\" into gallery search", () => enterInto(page, "Demographics overview", el("gallery search")));
    await session.step(89, "Then \"Demographics overview\" gallery card should be visible", () => shouldBe(page, el("\"Demographics overview\" gallery card"), "visible"));
    await session.step(91, "When user double-clicks on \"Demographics overview\" gallery card and the project opens", () => doubleClickOpensProject(page, el("\"Demographics overview\" gallery card")), undefined, "Double-click to open it");
    await session.step(92, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(93, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(94, "And bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
    await session.step(96, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")), undefined, "With View and use they cannot save over it, only save a copy");
    await session.step(97, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(98, "And Name text input in \"Save project\" dialog should have value \"Copy of Demographics overview\"", () => shouldHaveValue(page, el("Name text input in \"Save project\" dialog"), "Copy of Demographics overview"));
    await session.step(99, "When user hovers over \"Save original project\" text in \"Save project\" dialog", () => hoverOver(page, el("\"Save original project\" text in \"Save project\" dialog")));
    await session.step(100, "And user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
    await session.step(101, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(102, "And no errors should have been logged", () => noErrors(page));
  });
});
