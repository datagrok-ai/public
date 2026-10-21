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
import {contextPanelOpen, contextPanelShows, dialogCloses, layoutsDeleted, openDatasetRowsAs, pickSharingUser, runningAccountSignedIn, sharingUserNoNotifications, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, infoBalloonText, noBalloons, noErrors, pickFromContextMenu, propertyShouldBe, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A layout, shared and applied by the other user", () => {
  const session = feature(test, "features/sharing/sharing-layout.feature", import.meta.url);
  test("A layout, shared and applied by the other user", {tag: ["@journey", "@serial", "@realizes:sharing.share-dialog", "@realizes:views.layouts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(24, "And the layouts named \"ShareLayout-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("ShareLayout-{time}")));
    await session.step(25, "And the layouts named \"share-layout-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("share-layout-{time}")));
    await session.step(26, "And user opens demog dataset keeping the first 100 rows as \"share-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("share-layout-{time}")));
    await session.step(27, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]);
    await session.step(30, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(31, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(32, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await run.scenario("The owner saves a layout from the Layouts pane and shares it from its card", async () => {
      await session.step(35, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
      await session.step(36, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
      await session.step(37, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
      await session.step(38, "Then the context panel should show \"share-layout-{time}\"", () => contextPanelShows(page, session.text("share-layout-{time}")));
      await session.step(39, "When user types \"ShareLayout-{time}\" into Name details field in context panel", () => typeInto(page, session.text("ShareLayout-{time}"), el("Name details field in context panel")));
      await session.step(40, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(41, "Then \"ShareLayout-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"ShareLayout-{time}\" layout card")), "visible"));
      await session.step(42, "When user picks \"Share...\" from the context menu of \"ShareLayout-{time}\" layout card", () => pickFromContextMenu(page, "Share...", el(session.text("\"ShareLayout-{time}\" layout card"))));
      await session.step(43, "And user picks the sharing user in \"User, group, or email\" input in \"Share ShareLayout-{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share ShareLayout-{time}\" dialog"))));
      await session.step(44, "And user unchecks \"Send notifications\" input in \"Share ShareLayout-{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share ShareLayout-{time}\" dialog"))));
      await session.step(45, "And user clicks on OK button in \"Share ShareLayout-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share ShareLayout-{time}\" dialog"))));
      await session.step(46, "Then the \"Share ShareLayout-{time}\" dialog should close", () => dialogCloses(page, session.text("Share ShareLayout-{time}")));
      await session.step(47, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
      await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The second account gets the layout's card for a table of the same columns and applies it", async () => {
      await session.step(52, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(53, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(54, "Given user opens demog dataset keeping the first 100 rows as \"share-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 100, session.text("share-layout-{time}")));
      await session.step(55, "And the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(56, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
      await session.step(57, "Then \"ShareLayout-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"ShareLayout-{time}\" layout card")), "visible"));
      await session.step(58, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(59, "When user clicks on \"ShareLayout-{time}\" layout card", () => clickOn(page, el(session.text("\"ShareLayout-{time}\" layout card"))));
      await session.step(60, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(61, "And \"xColumnName\" property of scatter plot viewer should be \"HEIGHT\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "HEIGHT"));
      await session.step(62, "And \"yColumnName\" property of scatter plot viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer"), "WEIGHT"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
      await session.step(64, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(65, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(66, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    run.finish();
  });
});
