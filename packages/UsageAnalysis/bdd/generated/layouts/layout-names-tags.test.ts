/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/layouts/layout-names-tags.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.layouts, sharing.share-dialog]
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
import {clearField, clickOn, hoverOver, isExpanded, pressKey, shouldBe, shouldContainText, shouldHaveText, typeAtCaret, typeInto, uncheck, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, contextPanelOpen, contextPanelShows, contextPanelTitle, dialogCloses, layoutsDeleted, openDatasetRowsAs, pickSharingUser, sharingUserLevel, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, infoBalloonText, noBalloons, noErrors, pickFromContextMenu, propertyShouldBe, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A layout saved from the toolbox gets a name and tags", () => {
  const session = feature(test, "features/layouts/layout-names-tags.feature", import.meta.url);
  test("A layout saved from the toolbox gets a name and tags", {tag: ["@realizes:views.layouts", "@journey", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the layouts named \"bdd-layout-{time}, BDD height vs weight-{time}, BDD tags-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-layout-{time}, BDD height vs weight-{time}, BDD tags-{time}")));
    await session.step(20, "And user opens demog-1000 dataset keeping the first 100 rows as \"bdd-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog-1000"), 100, session.text("bdd-layout-{time}")));
    await session.step(21, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]);
    await session.step(24, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(25, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(26, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(27, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
    await session.step(28, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
    await run.scenario("Save in the Layouts pane names the layout after the table and opens it in the context panel", async () => {
      await session.step(31, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
      await session.step(32, "Then \"bdd-layout-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-layout-{time}\" layout card")), "visible"));
      await session.step(33, "And the context panel should show \"bdd-layout-{time}\"", () => contextPanelShows(page, session.text("bdd-layout-{time}")));
      await session.step(34, "And the title of context panel should be \"bdd-layout-{time}\"", () => contextPanelTitle(page, session.text("bdd-layout-{time}")));
      await session.step(35, "And context panel should contain text \"Layout saved\"", () => shouldContainText(page, el("context panel"), "Layout saved"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The name typed into the context panel renames the layout on its card and in the panel", async () => {
      await session.step(39, "When user types \"BDD height vs weight-{time}\" into Name details field in context panel", () => typeInto(page, session.text("BDD height vs weight-{time}"), el("Name details field in context panel")));
      await session.step(40, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(41, "Then \"BDD height vs weight-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"BDD height vs weight-{time}\" layout card")), "visible"));
      await session.step(42, "And \"bdd-layout-{time}\" layout card should be absent", () => shouldBe(page, el(session.text("\"bdd-layout-{time}\" layout card")), "absent"));
      await session.step(43, "And value of Name details field in context panel should have text \"BDD height vs weight-{time}\"", () => shouldHaveText(page, el("value of Name details field in context panel"), session.text("BDD height vs weight-{time}")));
      await session.step(44, "And the context panel should show \"BDD height vs weight-{time}\"", () => contextPanelShows(page, session.text("BDD height vs weight-{time}")));
      await session.step(45, "And the title of context panel should be \"BDD height vs weight-{time}\"", () => contextPanelTitle(page, session.text("BDD height vs weight-{time}")));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The pencil next to Tags adds a tag that the card shows", async () => {
      await session.step(49, "Given the context panel should show \"BDD height vs weight-{time}\"", () => contextPanelShows(page, session.text("BDD height vs weight-{time}")));
      await session.step(50, "When user hovers over Tags details row in context panel", () => hoverOver(page, el("Tags details row in context panel")));
      await session.step(51, "And user clicks on edit icon of Tags details row in context panel", () => clickOn(page, el("edit icon of Tags details row in context panel")));
      await session.step(52, "And user types \"anthropometry\" at the caret", () => typeAtCaret(page, "anthropometry"));
      await session.step(53, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(54, "Then \"BDD height vs weight-{time}\" layout card should contain text \"#anthropometry\"", () => shouldContainText(page, el(session.text("\"BDD height vs weight-{time}\" layout card")), "#anthropometry"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The filter of the Layouts pane finds the layout by part of its name and by the tag", async () => {
      await session.step(58, "When user types \"bdd-layout-{time}\" into layouts filter", () => typeInto(page, session.text("bdd-layout-{time}"), el("layouts filter")));
      await session.step(59, "Then there should be 0 visible layout cards", () => visibleCount(page, 0, el("layout cards")));
      await session.step(60, "When user types \"height vs weight-{time}\" into layouts filter", () => typeInto(page, session.text("height vs weight-{time}"), el("layouts filter")));
      await session.step(61, "Then \"BDD height vs weight-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"BDD height vs weight-{time}\" layout card")), "visible"));
      await session.step(62, "And there should be 1 visible layout card", () => visibleCount(page, 1, el("layout card")));
      await session.step(63, "When user types \"bdd-layout-{time}\" into layouts filter", () => typeInto(page, session.text("bdd-layout-{time}"), el("layouts filter")));
      await session.step(64, "Then there should be 0 visible layout cards", () => visibleCount(page, 0, el("layout cards")));
      await session.step(65, "When user types \"#anthropometry\" into layouts filter", () => typeInto(page, "#anthropometry", el("layouts filter")));
      await session.step(66, "Then \"BDD height vs weight-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"BDD height vs weight-{time}\" layout card")), "visible"));
      await session.step(67, "When user types \"bdd-layout-{time}\" into layouts filter", () => typeInto(page, session.text("bdd-layout-{time}"), el("layouts filter")));
      await session.step(68, "Then there should be 0 visible layout cards", () => visibleCount(page, 0, el("layout cards")));
      await session.step(69, "When user clears layouts filter", () => clearField(page, el("layouts filter")));
      await session.step(70, "Then \"BDD height vs weight-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"BDD height vs weight-{time}\" layout card")), "visible"));
    });
    await run.scenario("The table opened again lists the renamed, tagged layout, and a click on its card applies it", async () => {
      await session.step(73, "When user closes all views", () => closeAllViews(page));
      await session.step(74, "And user opens demog-1000 dataset keeping the first 100 rows as \"bdd-layout-{time}\"", () => openDatasetRowsAs(page, ds("demog-1000"), 100, session.text("bdd-layout-{time}")));
      await session.step(75, "And the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(76, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
      await session.step(77, "Then \"BDD height vs weight-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"BDD height vs weight-{time}\" layout card")), "visible"));
      await session.step(78, "And \"bdd-layout-{time}\" layout card should be absent", () => shouldBe(page, el(session.text("\"bdd-layout-{time}\" layout card")), "absent"));
      await session.step(79, "And \"BDD height vs weight-{time}\" layout card should contain text \"#anthropometry\"", () => shouldContainText(page, el(session.text("\"BDD height vs weight-{time}\" layout card")), "#anthropometry"));
      await session.step(80, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(81, "When user clicks on \"BDD height vs weight-{time}\" layout card", () => clickOn(page, el(session.text("\"BDD height vs weight-{time}\" layout card"))));
      await session.step(82, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(83, "And \"xColumnName\" property of scatter plot viewer should be \"HEIGHT\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "HEIGHT"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second layout saved and renamed takes a tag through Edit tags... on its card", async () => {
      await session.step(87, "When user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
      await session.step(88, "Then the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
      await session.step(89, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
      await session.step(90, "Then the context panel should show \"bdd-layout-{time}\"", () => contextPanelShows(page, session.text("bdd-layout-{time}")));
      await session.step(91, "When user types \"BDD tags-{time}\" into Name details field in context panel", () => typeInto(page, session.text("BDD tags-{time}"), el("Name details field in context panel")));
      await session.step(92, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(93, "Then \"BDD tags-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"BDD tags-{time}\" layout card")), "visible"));
      await session.step(94, "When user picks \"Edit tags...\" from the context menu of \"BDD tags-{time}\" layout card", () => pickFromContextMenu(page, "Edit tags...", el(session.text("\"BDD tags-{time}\" layout card"))));
      await session.step(95, "Then \"Edit tags\" dialog should be visible", () => shouldBe(page, el("\"Edit tags\" dialog"), "visible"));
      await session.step(96, "When user types \"anthropometry2\" into Tags input in \"Edit tags\" dialog", () => typeInto(page, "anthropometry2", el("Tags input in \"Edit tags\" dialog")));
      await session.step(97, "And user clicks on OK button in \"Edit tags\" dialog", () => clickOn(page, el("OK button in \"Edit tags\" dialog")));
      await session.step(98, "Then \"Edit tags\" dialog should be absent", () => shouldBe(page, el("\"Edit tags\" dialog"), "absent"));
      await session.step(99, "And \"BDD tags-{time}\" layout card should contain text \"#anthropometry2\"", () => shouldContainText(page, el(session.text("\"BDD tags-{time}\" layout card")), "#anthropometry2"));
      await session.step(100, "And \"bdd-layout-{time}\" layout card should be absent", () => shouldBe(page, el(session.text("\"bdd-layout-{time}\" layout card")), "absent"));
      await session.step(101, "When user types \"#anthropometry2\" into layouts filter", () => typeInto(page, "#anthropometry2", el("layouts filter")));
      await session.step(102, "Then there should be 1 visible layout card", () => visibleCount(page, 1, el("layout card")));
      await session.step(103, "And \"BDD tags-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"BDD tags-{time}\" layout card")), "visible"));
      await session.step(104, "When user clears layouts filter", () => clearField(page, el("layouts filter")));
    });
    await run.scenario("Share... on a layout card shares the layout with another user", async () => {
      await session.step(108, "When user picks \"Share...\" from the context menu of \"BDD tags-{time}\" layout card", () => pickFromContextMenu(page, "Share...", el(session.text("\"BDD tags-{time}\" layout card"))));
      await session.step(109, "Then \"Share BDD tags-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD tags-{time}\" dialog")), "visible"));
      await session.step(110, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(111, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDD tags-{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDD tags-{time}\" dialog"))));
      await session.step(112, "And user unchecks \"Send notifications\" input in \"Share BDD tags-{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDD tags-{time}\" dialog"))));
      await session.step(113, "And user clicks on OK button in \"Share BDD tags-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDD tags-{time}\" dialog"))));
      await session.step(114, "Then the \"Share BDD tags-{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDD tags-{time}")));
      await session.step(115, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(116, "When user picks \"Share...\" from the context menu of \"BDD tags-{time}\" layout card", () => pickFromContextMenu(page, "Share...", el(session.text("\"BDD tags-{time}\" layout card"))));
      await session.step(117, "Then \"Share BDD tags-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD tags-{time}\" dialog")), "visible"));
      await session.step(118, "And the access level of the sharing user in \"Share BDD tags-{time}\" dialog should be \"View and use\"", () => sharingUserLevel(page, el(session.text("\"Share BDD tags-{time}\" dialog")), "View and use"));
      await session.step(119, "When user clicks on CANCEL button in \"Share BDD tags-{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDD tags-{time}\" dialog"))));
      await session.step(120, "Then the \"Share BDD tags-{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDD tags-{time}")));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
      await session.step(122, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
