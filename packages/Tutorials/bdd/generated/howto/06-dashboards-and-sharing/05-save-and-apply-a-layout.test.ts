/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/06-dashboards-and-sharing/05-save-and-apply-a-layout.feature
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
import {ownLayoutsGone} from '../../../bindings/howto-dashboards.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, expand, hoverOver, isExpanded, pressKey, shouldBe, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelOpen, contextPanelShows, layoutsDeleted, openDatasetRowsAs, simpleModeOff, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, noErrors, propertyShouldBe, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Dashboards", () => {
  const session = feature(test, "features/howto/06-dashboards-and-sharing/05-save-and-apply-a-layout.feature", import.meta.url);
  test("Save a layout and apply it to another table", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the layouts named \"Demographics, Height and weight\" are deleted when the feature ends", () => layoutsDeleted(page, "Demographics, Height and weight"));
    await session.step(22, "And the user's own layouts named \"Height and weight\" are removed now and at feature end", () => ownLayoutsGone(page, "Height and weight"));
    await session.step(23, "And user opens demog dataset keeping the first 5850 rows as \"Demographics\"", () => openDatasetRowsAs(page, ds("demog"), 5850, "Demographics"));
    await session.step(24, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]);
    await session.step(28, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"]]), [["valueColumnName","AGE"]]);
    await session.step(30, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(31, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(32, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(33, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(34, "And user clicks on the \"header AGE\" area of grid", () => clickArea(page, "header AGE", el("grid")));
    await session.step(35, "And the context panel should show \"AGE\"", () => contextPanelShows(page, "AGE"));
    await session.step(36, "When the video says \"A layout is the arrangement of viewers, without the data\"", () => videoSays(page, "A layout is the arrangement of viewers, without the data"));
    await session.step(38, "And user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")), undefined, "Toolbox > Layouts > Save keeps this arrangement");
    await session.step(39, "Then \"Demographics\" layout card should be visible", () => shouldBe(page, el("\"Demographics\" layout card"), "visible"));
    await session.step(40, "And the context panel should show \"Demographics\"", () => contextPanelShows(page, "Demographics"));
    await session.step(42, "When user types \"Height and weight\" into Name details field in context panel", () => typeInto(page, "Height and weight", el("Name details field in context panel")), undefined, "It is named after the table. Give it a better name");
    await session.step(43, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(44, "Then \"Height and weight\" layout card should be visible", () => shouldBe(page, el("\"Height and weight\" layout card"), "visible"));
    await session.step(45, "When the video says \"Now another table with the same columns\"", () => videoSays(page, "Now another table with the same columns"));
    await session.step(46, "And user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(47, "Then browse panel should be visible", () => shouldBe(page, el("browse panel"), "visible"));
    await session.step(49, "When user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")), undefined, "Browse > Files > Demo: open demog-1000.csv");
    await session.step(50, "And user expands Files---Demo tree node inside browse tree", () => expand(page, el("Files---Demo tree node inside browse tree")));
    await session.step(51, "And user double-clicks on Files---Demo---demog-1000.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog-1000.csv tree node inside browse tree")));
    await session.step(52, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(53, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
    await session.step(54, "And toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
    await session.step(56, "When user expands Layouts accordion header in toolbox", () => expand(page, el("Layouts accordion header in toolbox")), undefined, "Its Layouts section offers the saved layout");
    await session.step(57, "And user hovers over \"Height and weight\" layout card", () => hoverOver(page, el("\"Height and weight\" layout card")));
    await session.step(59, "And user clicks on \"Height and weight\" layout card", () => clickOn(page, el("\"Height and weight\" layout card")), undefined, "Click it to arrange this table the same way");
    await session.step(60, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(61, "And the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
    await session.step(62, "And \"xColumnName\" property of scatter plot viewer should be \"HEIGHT\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "HEIGHT"));
    await session.step(63, "When the video says \"Same viewers, now showing this table's 1,000 rows\"", () => videoSays(page, "Same viewers, now showing this table's 1,000 rows"));
    await session.step(64, "And user hovers over scatter plot viewer", () => hoverOver(page, el("scatter plot viewer")));
    await session.step(65, "Then no errors should have been logged", () => noErrors(page));
  });
});
