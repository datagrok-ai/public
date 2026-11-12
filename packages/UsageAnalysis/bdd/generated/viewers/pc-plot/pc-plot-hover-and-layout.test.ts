/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot-hover-and-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pc-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {axesShouldBe, dragAxisLabel} from '../../../bindings/pc-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaShorter, areaTaller, dragSelectionOverArea, hasArea, hasNoArea, hoverArea, legendDocked, legendLists, legendMiniIcon, lessInk, moreHighlight, moreInk, noErrors, notRepainted, oneTooltip, pointerAway, readingIs, repaintedBy, resizeTo, restoreSize, setProperty, showsRows, takeSnapshot, tooltipValue, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot — the hovered line, the layout at a small size, a selection after a reorder, and Undo", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-hover-and-layout.feature", import.meta.url);
  test("PC plot — the hovered line, the layout at a small size, a selection after a reorder, and Undo", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]), [["Column Names","AGE, HEIGHT, WEIGHT"]]);
    await session.step(19, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(20, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(21, "And pc plot viewer should have a \"line of row 1000\" area", () => hasArea(page, el("pc plot viewer"), "line of row 1000"));
    await session.step(22, "And pc plot viewer should have an \"axis label \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis label \"AGE\""));
    await run.scenario("Hovering a polyline raises the tooltip of its row", async () => {
      await session.step(25, "When user hovers over the \"line of row 1000\" area of pc plot viewer", () => hoverArea(page, "line of row 1000", el("pc plot viewer")));
      await session.step(26, "Then the \"hovered row\" reading of pc plot viewer should be 1000", () => readingIs(page, "hovered row", el("pc plot viewer"), 1000));
      await session.step(27, "And exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(28, "And the tooltip should show \"AGE\" as \"63\"", () => tooltipValue(page, "AGE", "63"));
      await session.step(29, "And the tooltip should show \"WEIGHT\" as \"84.00\"", () => tooltipValue(page, "WEIGHT", "84.00"));
      await session.step(30, "When user hovers over the \"empty space\" area of pc plot viewer", () => hoverArea(page, "empty space", el("pc plot viewer")));
      await session.step(31, "Then the \"hovered row\" reading of pc plot viewer should be 0", () => readingIs(page, "hovered row", el("pc plot viewer"), 0));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Mouse Over Line Width thickens the hovered polyline", async () => {
      await session.step(35, "When user hovers over the \"line of row 1000\" area of pc plot viewer", () => hoverArea(page, "line of row 1000", el("pc plot viewer")));
      await session.step(36, "Then the \"hovered row\" reading of pc plot viewer should be 1000", () => readingIs(page, "hovered row", el("pc plot viewer"), 1000));
      await session.step(37, "When user sets \"Mouse Over Line Width\" property of pc plot viewer to \"8\"", () => setProperty(page, "Mouse Over Line Width", el("pc plot viewer"), "8"));
      await session.step(38, "Then the \"hovered row\" reading of pc plot viewer should be 1000", () => readingIs(page, "hovered row", el("pc plot viewer"), 1000));
      await session.step(39, "And pc plot viewer should have more ink than before", () => moreInk(page, el("pc plot viewer")));
      await session.step(40, "When user sets \"Mouse Over Line Width\" property of pc plot viewer to \"2\"", () => setProperty(page, "Mouse Over Line Width", el("pc plot viewer"), "2"));
      await session.step(41, "Then pc plot viewer should have less ink than before", () => lessInk(page, el("pc plot viewer")));
      await session.step(42, "When user hovers over the \"empty space\" area of pc plot viewer", () => hoverArea(page, "empty space", el("pc plot viewer")));
      await session.step(43, "Then the \"hovered row\" reading of pc plot viewer should be 0", () => readingIs(page, "hovered row", el("pc plot viewer"), 0));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A category hovered in another viewer fades the plot behind its rows, and the flag switches it off", async () => {
      await session.step(47, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
      await session.step(49, "And user takes a snapshot of pc plot viewer", () => takeSnapshot(page, el("pc plot viewer")));
      await session.step(50, "And user hovers over the \"bar Asian\" area of bar chart viewer", () => hoverArea(page, "bar Asian", el("bar chart viewer")));
      await session.step(51, "Then pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(52, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
      await session.step(53, "And user sets \"Show Mouse Over Row Group\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Mouse Over Row Group", el("pc plot viewer"), "false"));
      await session.step(54, "And user takes a snapshot of pc plot viewer", () => takeSnapshot(page, el("pc plot viewer")));
      await session.step(55, "And user hovers over the \"bar Asian\" area of bar chart viewer", () => hoverArea(page, "bar Asian", el("bar chart viewer")));
      await session.step(56, "Then pc plot viewer should not have repainted", () => notRepainted(page, el("pc plot viewer")));
      await session.step(57, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
      await session.step(58, "And user sets \"Show Mouse Over Row Group\" property of pc plot viewer to \"true\"", () => setProperty(page, "Show Mouse Over Row Group", el("pc plot viewer"), "true"));
      await session.step(59, "And user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(60, "Then bar chart viewer should be absent", () => shouldBe(page, el("bar chart viewer"), "absent"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Min Max Orientation Vert turns the min and max captions upright", async () => {
      await session.step(64, "Then pc plot viewer should have an \"axis min \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis min \"AGE\""));
      await session.step(65, "When user sets \"Min Max Orientation\" property of pc plot viewer to \"Vert\"", () => setProperty(page, "Min Max Orientation", el("pc plot viewer"), "Vert"));
      await session.step(66, "Then the \"axis min \\\"AGE\\\"\" area of pc plot viewer should be taller than before", () => areaTaller(page, "axis min \"AGE\"", el("pc plot viewer")));
      await session.step(67, "And the \"view\" area of pc plot viewer should be shorter than before", () => areaShorter(page, "view", el("pc plot viewer")));
      await session.step(68, "When user sets \"Min Max Orientation\" property of pc plot viewer to \"Auto\"", () => setProperty(page, "Min Max Orientation", el("pc plot viewer"), "Auto"));
      await session.step(69, "Then the \"axis min \\\"AGE\\\"\" area of pc plot viewer should be shorter than before", () => areaShorter(page, "axis min \"AGE\"", el("pc plot viewer")));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shrunk below the layout gate the plot drops its labels and min/max, and Auto Layout off keeps the labels", async () => {
      await session.step(73, "When user resizes pc plot viewer to 200 by 200", () => resizeTo(page, el("pc plot viewer"), 200, 200));
      await session.step(74, "Then pc plot viewer should not have an \"axis label \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis label \"AGE\""));
      await session.step(75, "And pc plot viewer should not have an \"axis min \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis min \"AGE\""));
      await session.step(76, "And the \"axes\" reading of pc plot viewer should be 3", () => readingIs(page, "axes", el("pc plot viewer"), 3));
      await session.step(77, "When user sets \"Auto Layout\" property of pc plot viewer to \"false\"", () => setProperty(page, "Auto Layout", el("pc plot viewer"), "false"));
      await session.step(78, "Then pc plot viewer should have an \"axis label \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis label \"AGE\""));
      await session.step(79, "When user sets \"Auto Layout\" property of pc plot viewer to \"true\"", () => setProperty(page, "Auto Layout", el("pc plot viewer"), "true"));
      await session.step(80, "Then pc plot viewer should not have an \"axis label \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis label \"AGE\""));
      await session.step(81, "When user restores the size of pc plot viewer", () => restoreSize(page, el("pc plot viewer")));
      await session.step(82, "Then pc plot viewer should have an \"axis label \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis label \"AGE\""));
      await session.step(83, "And pc plot viewer should have an \"axis min \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis min \"AGE\""));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shrunk small, a coloured plot folds its legend into the mini icon", async () => {
      await session.step(87, "When user sets \"Color\" property of pc plot viewer to \"RACE\"", () => setProperty(page, "Color", el("pc plot viewer"), "RACE"));
      await session.step(88, "Then the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(89, "And the legend of pc plot viewer should be docked", () => legendDocked(page, el("pc plot viewer")));
      await session.step(90, "When user resizes pc plot viewer to 200 by 200", () => resizeTo(page, el("pc plot viewer"), 200, 200));
      await session.step(91, "Then the legend of pc plot viewer should be collapsed to the mini icon", () => legendMiniIcon(page, el("pc plot viewer")));
      await session.step(92, "When user restores the size of pc plot viewer", () => restoreSize(page, el("pc plot viewer")));
      await session.step(93, "Then the legend of pc plot viewer should be docked", () => legendDocked(page, el("pc plot viewer")));
      await session.step(94, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(95, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Shift-drag selects over the band a label reorder made", async () => {
      await session.step(98, "When user drags the \"WEIGHT\" axis label of pc plot viewer onto the \"AGE\" axis label", () => dragAxisLabel(page, "WEIGHT", "AGE"));
      await session.step(99, "Then the axes of pc plot viewer should be \"WEIGHT, AGE, HEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "WEIGHT, AGE, HEIGHT"));
      await session.step(100, "And pc plot viewer should have a \"band \\\"WEIGHT\\\" - \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "band \"WEIGHT\" - \"AGE\""));
      await session.step(103, "When user drags a selection box over the \"band \\\"WEIGHT\\\" - \\\"AGE\\\"\" area of pc plot viewer", () => dragSelectionOverArea(page, "band \"WEIGHT\" - \"AGE\"", el("pc plot viewer")));
      await session.step(104, "Then 986 rows should be selected", () => selectedRowCount(page, 986));
      await session.step(105, "And the \"rows selected\" reading of pc plot viewer should be 986", () => readingIs(page, "rows selected", el("pc plot viewer"), 986));
      await session.step(106, "And pc plot viewer should show more selection highlight than before", () => moreHighlight(page, el("pc plot viewer")));
      await session.step(107, "When user clears the row selection", () => clearSelection(page));
      await session.step(108, "And user sets \"Column Names\" property of pc plot viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "Column Names", el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(109, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A plot closed from its title bar comes back with Ctrl+Z", async () => {
      await session.step(113, "When user clicks on close icon of pc plot viewer", () => clickOn(page, el("close icon of pc plot viewer")));
      await session.step(114, "Then the open tableview should have 0 pc plot viewers", () => viewerCount(page, 0, "pc plot"));
      await session.step(115, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(116, "Then the open tableview should have 1 pc plot viewer", () => viewerCount(page, 1, "pc plot"));
      await session.step(117, "And the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(118, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(119, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
