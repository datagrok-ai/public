/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/06-dashboards-and-sharing/07-share-a-space.feature
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
import {clickOn, doubleClickOn, enterInto, expand, hoverOver, selectIn, shouldBe, shouldContainText, typeInto, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {balloonsClosed, browsePanelOpen, closeAllViews, contextPanelClosed, dialogCloses, noProjectOnServer, noSpaceOnServer, openDataset, sharingUserNoNotifications, simpleModeOff, spaceOnServer, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Dashboards", () => {
  const session = feature(test, "features/howto/06-dashboards-and-sharing/07-share-a-space.feature", import.meta.url);
  test("Share a space with your team", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(21, "And no project named \"Demographics overview\" is on the server", () => noProjectOnServer(page, "Demographics overview"));
    await session.step(22, "And no space named \"Clinical team\" is on the server", () => noSpaceOnServer(page, "Clinical team"));
    await session.step(23, "And a space named \"Clinical team\" is on the server", () => spaceOnServer(page, "Clinical team"));
    await session.step(24, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(25, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]);
    await session.step(29, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"]]), [["splitColumnName","RACE"]]);
    await session.step(31, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(32, "And user enters \"Demographics overview\" into Name text input in \"Save project\" dialog", () => enterInto(page, "Demographics overview", el("Name text input in \"Save project\" dialog")));
    await session.step(33, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(34, "And the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(35, "And \"Share Demographics overview\" dialog should be visible", () => shouldBe(page, el("\"Share Demographics overview\" dialog"), "visible"));
    await session.step(36, "And user clicks on CANCEL button in \"Share Demographics overview\" dialog", () => clickOn(page, el("CANCEL button in \"Share Demographics overview\" dialog")));
    await session.step(37, "And the \"Share Demographics overview\" dialog should close", () => dialogCloses(page, "Share Demographics overview"));
    await session.step(38, "And user closes all views", () => closeAllViews(page));
    await session.step(39, "And the notification balloons are closed", () => balloonsClosed(page));
    await session.step(40, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(41, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(42, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(44, "And user enters \"Demographics overview\" into gallery search", () => enterInto(page, "Demographics overview", el("gallery search")));
    await session.step(45, "And \"Demographics overview\" gallery card should be visible", () => shouldBe(page, el("\"Demographics overview\" gallery card"), "visible"));
    await session.step(46, "When the video says \"To share many dashboards at once, put them in a space\"", () => videoSays(page, "To share many dashboards at once, put them in a space"));
    await session.step(48, "And user picks \"Move to Space...\" from the context menu of \"Demographics overview\" gallery card", () => pickFromContextMenu(page, "Move to Space...", el("\"Demographics overview\" gallery card")), undefined, "Move the dashboard into your team's space");
    await session.step(49, "Then \"Move to space\" dialog should be visible", () => shouldBe(page, el("\"Move to space\" dialog"), "visible"));
    await session.step(50, "When user selects \"Clinical team\" in Space input in \"Move to space\" dialog", () => selectIn(page, "Clinical team", el("Space input in \"Move to space\" dialog")));
    await session.step(51, "And user clicks on OK button in \"Move to space\" dialog", () => clickOn(page, el("OK button in \"Move to space\" dialog")));
    await session.step(52, "Then the \"Move to space\" dialog should close", () => dialogCloses(page, "Move to space"));
    await session.step(54, "When user expands Spaces tree node inside browse tree", () => expand(page, el("Spaces tree node inside browse tree")), undefined, "A space keeps a team's dashboards, files and queries together");
    await session.step(55, "And user double-clicks on \"Spaces > Clinical team\" tree node inside browse tree", () => doubleClickOn(page, el("\"Spaces > Clinical team\" tree node inside browse tree")));
    await session.step(56, "Then the \"Clinical team\" view should be current", () => viewIsCurrent(page, "Clinical team"));
    await session.step(57, "And \"Demographics overview\" link in gallery should be visible", () => shouldBe(page, el("\"Demographics overview\" link in gallery"), "visible"));
    await session.step(59, "When user picks \"Share...\" from the context menu of \"Spaces > Clinical team\" tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("\"Spaces > Clinical team\" tree node inside browse tree")), undefined, "Share the space once");
    await session.step(60, "Then \"Share Clinical team\" dialog should be visible", () => shouldBe(page, el("\"Share Clinical team\" dialog"), "visible"));
    await session.step(61, "When user types \"bddsecond\" into \"User, group, or email\" input in \"Share Clinical team\" dialog", () => typeInto(page, "bddsecond", el("\"User, group, or email\" input in \"Share Clinical team\" dialog")));
    await session.step(62, "And user clicks on \"bddsecond\" table row in \"Share Clinical team\" dialog", () => clickOn(page, el("\"bddsecond\" table row in \"Share Clinical team\" dialog")));
    await session.step(63, "And user unchecks \"Send notifications\" input in \"Share Clinical team\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share Clinical team\" dialog")));
    await session.step(64, "And user clicks on OK button in \"Share Clinical team\" dialog", () => clickOn(page, el("OK button in \"Share Clinical team\" dialog")));
    await session.step(65, "Then the \"Share Clinical team\" dialog should close", () => dialogCloses(page, "Share Clinical team"));
    await session.step(67, "When user picks \"Share...\" from the context menu of \"Demographics overview\" link in gallery", () => pickFromContextMenu(page, "Share...", el("\"Demographics overview\" link in gallery")), undefined, "Now look at the dashboard inside");
    await session.step(68, "Then \"Share Demographics overview\" dialog should be visible", () => shouldBe(page, el("\"Share Demographics overview\" dialog"), "visible"));
    await session.step(69, "And \"Share Demographics overview\" dialog should contain text \"Inherited from\"", () => shouldContainText(page, el("\"Share Demographics overview\" dialog"), "Inherited from"));
    await session.step(70, "And \"Share Demographics overview\" dialog should contain text \"Clinical team\"", () => shouldContainText(page, el("\"Share Demographics overview\" dialog"), "Clinical team"));
    await session.step(72, "When user hovers over \"Share Demographics overview\" dialog", () => hoverOver(page, el("\"Share Demographics overview\" dialog")), undefined, "It has the same access, inherited from the space");
    await session.step(73, "Then no errors should have been logged", () => noErrors(page));
  });
});
