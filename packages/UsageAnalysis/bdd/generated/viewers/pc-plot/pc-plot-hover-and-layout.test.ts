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
import {clickOn, pressKey} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaShorter, areaTaller, dragSelectionOverArea, hasArea, hasNoArea, hoverArea, legendDocked, legendLists, legendMiniIcon, lessInk, moreHighlight, moreInk, noErrors, oneTooltip, readingIs, resizeTo, restoreSize, setProperty, showsRows, tooltipValue, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot — the hovered line, the layout at a small size, a selection after a reorder, and Undo", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-hover-and-layout.feature", import.meta.url);
  test("PC plot — the hovered line, the layout at a small size, a selection after a reorder, and Undo", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]), [["Column Names","AGE, HEIGHT, WEIGHT"]]);
    await session.step(18, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(19, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(20, "And pc plot viewer should have a \"line of row 1000\" area", () => hasArea(page, el("pc plot viewer"), "line of row 1000"));
    await session.step(21, "And pc plot viewer should have an \"axis label \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis label \"AGE\""));
    await run.scenario("Hovering a polyline raises the tooltip of its row", async () => {
      await session.step(24, "When user hovers over the \"line of row 1000\" area of pc plot viewer", () => hoverArea(page, "line of row 1000", el("pc plot viewer")));
      await session.step(25, "Then the \"hovered row\" reading of pc plot viewer should be 1000", () => readingIs(page, "hovered row", el("pc plot viewer"), 1000));
      await session.step(26, "And exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(27, "And the tooltip should show \"AGE\" as \"63\"", () => tooltipValue(page, "AGE", "63"));
      await session.step(28, "And the tooltip should show \"WEIGHT\" as \"84.00\"", () => tooltipValue(page, "WEIGHT", "84.00"));
      await session.step(29, "When user hovers over the \"empty space\" area of pc plot viewer", () => hoverArea(page, "empty space", el("pc plot viewer")));
      await session.step(30, "Then the \"hovered row\" reading of pc plot viewer should be 0", () => readingIs(page, "hovered row", el("pc plot viewer"), 0));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Mouse Over Line Width thickens the hovered polyline", async () => {
      await session.step(34, "When user hovers over the \"line of row 1000\" area of pc plot viewer", () => hoverArea(page, "line of row 1000", el("pc plot viewer")));
      await session.step(35, "Then the \"hovered row\" reading of pc plot viewer should be 1000", () => readingIs(page, "hovered row", el("pc plot viewer"), 1000));
      await session.step(36, "When user sets \"Mouse Over Line Width\" property of pc plot viewer to \"8\"", () => setProperty(page, "Mouse Over Line Width", el("pc plot viewer"), "8"));
      await session.step(37, "Then the \"hovered row\" reading of pc plot viewer should be 1000", () => readingIs(page, "hovered row", el("pc plot viewer"), 1000));
      await session.step(38, "And pc plot viewer should have more ink than before", () => moreInk(page, el("pc plot viewer")));
      await session.step(39, "When user sets \"Mouse Over Line Width\" property of pc plot viewer to \"2\"", () => setProperty(page, "Mouse Over Line Width", el("pc plot viewer"), "2"));
      await session.step(40, "Then pc plot viewer should have less ink than before", () => lessInk(page, el("pc plot viewer")));
      await session.step(41, "When user hovers over the \"empty space\" area of pc plot viewer", () => hoverArea(page, "empty space", el("pc plot viewer")));
      await session.step(42, "Then the \"hovered row\" reading of pc plot viewer should be 0", () => readingIs(page, "hovered row", el("pc plot viewer"), 0));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Min Max Orientation Vert turns the min and max captions upright", async () => {
      await session.step(46, "Then pc plot viewer should have an \"axis min \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis min \"AGE\""));
      await session.step(47, "When user sets \"Min Max Orientation\" property of pc plot viewer to \"Vert\"", () => setProperty(page, "Min Max Orientation", el("pc plot viewer"), "Vert"));
      await session.step(48, "Then the \"axis min \\\"AGE\\\"\" area of pc plot viewer should be taller than before", () => areaTaller(page, "axis min \"AGE\"", el("pc plot viewer")));
      await session.step(49, "And the \"view\" area of pc plot viewer should be shorter than before", () => areaShorter(page, "view", el("pc plot viewer")));
      await session.step(50, "When user sets \"Min Max Orientation\" property of pc plot viewer to \"Auto\"", () => setProperty(page, "Min Max Orientation", el("pc plot viewer"), "Auto"));
      await session.step(51, "Then the \"axis min \\\"AGE\\\"\" area of pc plot viewer should be shorter than before", () => areaShorter(page, "axis min \"AGE\"", el("pc plot viewer")));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shrunk below the layout gate the plot drops its labels and min/max, and Auto Layout off keeps the labels", async () => {
      await session.step(55, "When user resizes pc plot viewer to 200 by 200", () => resizeTo(page, el("pc plot viewer"), 200, 200));
      await session.step(56, "Then pc plot viewer should not have an \"axis label \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis label \"AGE\""));
      await session.step(57, "And pc plot viewer should not have an \"axis min \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis min \"AGE\""));
      await session.step(58, "And the \"axes\" reading of pc plot viewer should be 3", () => readingIs(page, "axes", el("pc plot viewer"), 3));
      await session.step(59, "When user sets \"Auto Layout\" property of pc plot viewer to \"false\"", () => setProperty(page, "Auto Layout", el("pc plot viewer"), "false"));
      await session.step(60, "Then pc plot viewer should have an \"axis label \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis label \"AGE\""));
      await session.step(61, "When user sets \"Auto Layout\" property of pc plot viewer to \"true\"", () => setProperty(page, "Auto Layout", el("pc plot viewer"), "true"));
      await session.step(62, "Then pc plot viewer should not have an \"axis label \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis label \"AGE\""));
      await session.step(63, "When user restores the size of pc plot viewer", () => restoreSize(page, el("pc plot viewer")));
      await session.step(64, "Then pc plot viewer should have an \"axis label \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis label \"AGE\""));
      await session.step(65, "And pc plot viewer should have an \"axis min \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis min \"AGE\""));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shrunk small, a coloured plot folds its legend into the mini icon", async () => {
      await session.step(69, "When user sets \"Color\" property of pc plot viewer to \"RACE\"", () => setProperty(page, "Color", el("pc plot viewer"), "RACE"));
      await session.step(70, "Then the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(71, "And the legend of pc plot viewer should be docked", () => legendDocked(page, el("pc plot viewer")));
      await session.step(72, "When user resizes pc plot viewer to 200 by 200", () => resizeTo(page, el("pc plot viewer"), 200, 200));
      await session.step(73, "Then the legend of pc plot viewer should be collapsed to the mini icon", () => legendMiniIcon(page, el("pc plot viewer")));
      await session.step(74, "When user restores the size of pc plot viewer", () => restoreSize(page, el("pc plot viewer")));
      await session.step(75, "Then the legend of pc plot viewer should be docked", () => legendDocked(page, el("pc plot viewer")));
      await session.step(76, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(77, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Shift-drag selects over the band a label reorder made", async () => {
      await session.step(80, "When user drags the \"WEIGHT\" axis label of pc plot viewer onto the \"AGE\" axis label", () => dragAxisLabel(page, "WEIGHT", "AGE"));
      await session.step(81, "Then the axes of pc plot viewer should be \"WEIGHT, AGE, HEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "WEIGHT, AGE, HEIGHT"));
      await session.step(82, "And pc plot viewer should have a \"band \\\"WEIGHT\\\" - \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "band \"WEIGHT\" - \"AGE\""));
      await session.step(85, "When user drags a selection box over the \"band \\\"WEIGHT\\\" - \\\"AGE\\\"\" area of pc plot viewer", () => dragSelectionOverArea(page, "band \"WEIGHT\" - \"AGE\"", el("pc plot viewer")));
      await session.step(86, "Then 986 rows should be selected", () => selectedRowCount(page, 986));
      await session.step(87, "And the \"rows selected\" reading of pc plot viewer should be 986", () => readingIs(page, "rows selected", el("pc plot viewer"), 986));
      await session.step(88, "And pc plot viewer should show more selection highlight than before", () => moreHighlight(page, el("pc plot viewer")));
      await session.step(89, "When user clears the row selection", () => clearSelection(page));
      await session.step(90, "And user sets \"Column Names\" property of pc plot viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "Column Names", el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(91, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A plot closed from its title bar comes back with Ctrl+Z", async () => {
      await session.step(95, "When user clicks on close icon of pc plot viewer", () => clickOn(page, el("close icon of pc plot viewer")));
      await session.step(96, "Then the open tableview should have 0 pc plot viewers", () => viewerCount(page, 0, "pc plot"));
      await session.step(97, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(98, "Then the open tableview should have 1 pc plot viewer", () => viewerCount(page, 1, "pc plot"));
      await session.step(99, "And the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(100, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
