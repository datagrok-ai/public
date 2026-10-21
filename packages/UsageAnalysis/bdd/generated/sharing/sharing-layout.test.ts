/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/sharing-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog, views.layouts]
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
import {clickOn, isExpanded, pressKey, shouldBe, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelOpen, contextPanelShows, dialogCloses, layoutsDeleted, openDatasetRowsAs, pickSharingUser, removeSharingUser, runningAccountSignedIn, sharingUserLevel, sharingUserNoNotifications, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, infoBalloonText, noBalloons, noErrors, pickFromContextMenu, propertyShouldBe, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A layout, shared, applied by the other user, and taken back", () => {
  const session = feature(test, "features/sharing/sharing-layout.feature", import.meta.url);
  test("A layout, shared, applied by the other user, and taken back", {tag: ["@journey", "@serial", "@realizes:sharing.share-dialog", "@realizes:views.layouts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(27, "And the layouts named \"ShareLayout-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("ShareLayout-{time}")));
    await session.step(28, "And the layouts named \"ShareLayoutKeep-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("ShareLayoutKeep-{time}")));
    await session.step(29, "And the layouts named \"share-layout-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("share-layout-{time}")));
    await session.step(30, "And user opens demog dataset keeping the first 100 rows as \"share-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("share-layout-{time}")));
    await session.step(31, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]);
    await session.step(34, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(35, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(36, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await run.scenario("The owner saves two layouts from the Layouts pane and shares them from their cards", async () => {
      await session.step(39, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
      await session.step(40, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
      await session.step(41, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
      await session.step(42, "Then the context panel should show \"share-layout-{time}\"", () => contextPanelShows(page, session.text("share-layout-{time}")));
      await session.step(43, "When user types \"ShareLayout-{time}\" into Name details field in context panel", () => typeInto(page, session.text("ShareLayout-{time}"), el("Name details field in context panel")));
      await session.step(44, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(45, "Then \"ShareLayout-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"ShareLayout-{time}\" layout card")), "visible"));
      await session.step(46, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
      await session.step(47, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
      await session.step(48, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
      await session.step(49, "Then the context panel should show \"share-layout-{time}\"", () => contextPanelShows(page, session.text("share-layout-{time}")));
      await session.step(50, "When user types \"ShareLayoutKeep-{time}\" into Name details field in context panel", () => typeInto(page, session.text("ShareLayoutKeep-{time}"), el("Name details field in context panel")));
      await session.step(51, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(52, "Then \"ShareLayout-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"ShareLayout-{time}\" layout card")), "visible"));
      await session.step(53, "And \"ShareLayoutKeep-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"ShareLayoutKeep-{time}\" layout card")), "visible"));
      await session.step(54, "When user picks \"Share...\" from the context menu of \"ShareLayoutKeep-{time}\" layout card", () => pickFromContextMenu(page, "Share...", el(session.text("\"ShareLayoutKeep-{time}\" layout card"))));
      await session.step(55, "And user picks the sharing user in \"User, group, or email\" input in \"Share ShareLayoutKeep-{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share ShareLayoutKeep-{time}\" dialog"))));
      await session.step(56, "And user unchecks \"Send notifications\" input in \"Share ShareLayoutKeep-{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share ShareLayoutKeep-{time}\" dialog"))));
      await session.step(57, "And user clicks on OK button in \"Share ShareLayoutKeep-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share ShareLayoutKeep-{time}\" dialog"))));
      await session.step(58, "Then the \"Share ShareLayoutKeep-{time}\" dialog should close", () => dialogCloses(page, session.text("Share ShareLayoutKeep-{time}")));
      await session.step(59, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(60, "When user picks \"Share...\" from the context menu of \"ShareLayout-{time}\" layout card", () => pickFromContextMenu(page, "Share...", el(session.text("\"ShareLayout-{time}\" layout card"))));
      await session.step(61, "And user picks the sharing user in \"User, group, or email\" input in \"Share ShareLayout-{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share ShareLayout-{time}\" dialog"))));
      await session.step(62, "And user unchecks \"Send notifications\" input in \"Share ShareLayout-{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share ShareLayout-{time}\" dialog"))));
      await session.step(63, "And user clicks on OK button in \"Share ShareLayout-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share ShareLayout-{time}\" dialog"))));
      await session.step(64, "Then the \"Share ShareLayout-{time}\" dialog should close", () => dialogCloses(page, session.text("Share ShareLayout-{time}")));
      await session.step(65, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The second account gets the layout's card for a table of the same columns and applies it", async () => {
      await session.step(70, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(71, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(72, "Given user opens demog dataset keeping the first 100 rows as \"share-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("share-layout-{time}")));
      await session.step(73, "And the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(74, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
      await session.step(75, "Then \"ShareLayout-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"ShareLayout-{time}\" layout card")), "visible"));
      await session.step(76, "And \"ShareLayoutKeep-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"ShareLayoutKeep-{time}\" layout card")), "visible"));
      await session.step(77, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(78, "When user clicks on \"ShareLayout-{time}\" layout card", () => clickOn(page, el(session.text("\"ShareLayout-{time}\" layout card"))));
      await session.step(79, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(80, "And \"xColumnName\" property of scatter plot viewer should be \"HEIGHT\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "HEIGHT"));
      await session.step(81, "And \"yColumnName\" property of scatter plot viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer"), "WEIGHT"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
      await session.step(83, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The owner takes the share back", async () => {
      await session.step(86, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(87, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(88, "Given user opens demog dataset keeping the first 100 rows as \"share-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("share-layout-{time}")));
      await session.step(89, "And the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(90, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
      await session.step(91, "When user picks \"Share...\" from the context menu of \"ShareLayout-{time}\" layout card", () => pickFromContextMenu(page, "Share...", el(session.text("\"ShareLayout-{time}\" layout card"))));
      await session.step(92, "Then \"Share ShareLayout-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share ShareLayout-{time}\" dialog")), "visible"));
      await session.step(93, "And the access level of the sharing user in \"Share ShareLayout-{time}\" dialog should be \"View and use\"", () => sharingUserLevel(page, el(session.text("\"Share ShareLayout-{time}\" dialog")), "View and use"));
      await session.step(94, "When user removes the sharing user from \"Share ShareLayout-{time}\" dialog", () => removeSharingUser(page, el(session.text("\"Share ShareLayout-{time}\" dialog"))));
      await session.step(95, "And user clicks on OK button in \"Share ShareLayout-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share ShareLayout-{time}\" dialog"))));
      await session.step(96, "Then the \"Share ShareLayout-{time}\" dialog should close", () => dialogCloses(page, session.text("Share ShareLayout-{time}")));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The second account no longer gets the layout", async () => {
      await session.step(100, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(101, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(102, "Given user opens demog dataset keeping the first 100 rows as \"share-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("share-layout-{time}")));
      await session.step(103, "And the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(104, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
      await session.step(105, "Then \"ShareLayoutKeep-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"ShareLayoutKeep-{time}\" layout card")), "visible"));
      await session.step(106, "And \"ShareLayout-{time}\" layout card should be absent", () => shouldBe(page, el(session.text("\"ShareLayout-{time}\" layout card")), "absent"));
      await session.step(107, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(108, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    run.finish();
  });
});
