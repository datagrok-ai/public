/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-selection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {currentRowIs, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, filterPasses, selectNoRows, selectWhereIs, selectWhereOneOf, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, dragDeselectionOverArea, dragSelectionOverArea, lessHighlight, lessInk, moreHighlight, moreInk, noErrors, noHighlight, openContextMenu, pickFromContextMenu, propertyShouldBe, readingIs, readingReads, repainted, reportsNoError, setProperty, someHighlight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragLasso} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart selection and the row markers", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-selection.feature", import.meta.url);
  test("Line chart selection and the row markers", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(26, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"],["showSelectedRows","true"]]));
    await session.step(30, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(31, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(32, "And the \"rows selected\" reading of line chart viewer should be 0", () => readingIs(page, "rows selected", el("line chart viewer"), 0));
    await session.step(33, "And line chart viewer should show no selection highlight", () => noHighlight(page, el("line chart viewer")));
    await session.step(34, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("A Shift-drag selects every row whose X falls in the band it spans", async () => {
      await session.step(37, "Given user clears the row selection", () => clearSelection(page));
      await session.step(38, "When user drags a selection box over the \"plot\" area of line chart viewer", () => dragSelectionOverArea(page, "plot", el("line chart viewer")));
      await session.step(39, "Then the \"rows selected\" reading of line chart viewer should be 82", () => readingIs(page, "rows selected", el("line chart viewer"), 82));
      await session.step(40, "And 82 rows should be selected", () => selectedRowCount(page, 82));
      await session.step(41, "And line chart viewer should show more selection highlight than before", () => moreHighlight(page, el("line chart viewer")));
      await session.step(42, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(43, "When user drags a deselection box over the \"plot\" area of line chart viewer", () => dragDeselectionOverArea(page, "plot", el("line chart viewer")));
      await session.step(44, "Then the \"rows selected\" reading of line chart viewer should be 0", () => readingIs(page, "rows selected", el("line chart viewer"), 0));
      await session.step(45, "And line chart viewer should show no selection highlight", () => noHighlight(page, el("line chart viewer")));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Lasso Tool property does not change what a Shift-drag selects", async () => {
      await session.step(49, "Given user clears the row selection", () => clearSelection(page));
      await session.step(50, "When user sets \"lassoTool\" property of line chart viewer to \"true\"", () => setProperty(page, "lassoTool", el("line chart viewer"), "true"));
      await session.step(51, "And user drags a selection box over the \"plot\" area of line chart viewer", () => dragSelectionOverArea(page, "plot", el("line chart viewer")));
      await session.step(52, "Then the \"rows selected\" reading of line chart viewer should be 82", () => readingIs(page, "rows selected", el("line chart viewer"), 82));
      await session.step(53, "When user clears the row selection", () => clearSelection(page));
      await session.step(54, "And user sets \"lassoTool\" property of line chart viewer to \"false\"", () => setProperty(page, "lassoTool", el("line chart viewer"), "false"));
      await session.step(55, "And user drags a selection box over the \"plot\" area of line chart viewer", () => dragSelectionOverArea(page, "plot", el("line chart viewer")));
      await session.step(56, "Then the \"rows selected\" reading of line chart viewer should be 82", () => readingIs(page, "rows selected", el("line chart viewer"), 82));
      await session.step(57, "When user clears the row selection", () => clearSelection(page));
      await session.step(58, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Lasso Tool is offered under Tools, and picking it flips the property", async () => {
      await session.step(61, "Then \"lassoTool\" property of line chart viewer should be \"false\"", () => propertyShouldBe(page, "lassoTool", el("line chart viewer"), "false"));
      await session.step(62, "When user picks \"Tools > Lasso Tool\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Tools > Lasso Tool", el("line chart viewer")));
      await session.step(63, "Then \"lassoTool\" property of line chart viewer should be \"true\"", () => propertyShouldBe(page, "lassoTool", el("line chart viewer"), "true"));
      await session.step(64, "When user picks \"Tools > Lasso Tool\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Tools > Lasso Tool", el("line chart viewer")));
      await session.step(65, "Then \"lassoTool\" property of line chart viewer should be \"false\"", () => propertyShouldBe(page, "lassoTool", el("line chart viewer"), "false"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A drag that ends where it started selects an empty band", async () => {
      await session.step(69, "Given user clears the row selection", () => clearSelection(page));
      await session.step(70, "When user drags a lasso over the \"plot\" area of line chart viewer", () => dragLasso(page, "plot", el("line chart viewer")));
      await session.step(71, "Then the \"rows selected\" reading of line chart viewer should be 0", () => readingIs(page, "rows selected", el("line chart viewer"), 0));
      await session.step(72, "And line chart viewer should show no selection highlight", () => noHighlight(page, el("line chart viewer")));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Selected Rows decides whether the selection is painted at all", async () => {
      await session.step(76, "When user selects rows where \"Stereo Category\" is \"R_ONE\"", () => selectWhereIs(page, "Stereo Category", "R_ONE"));
      await session.step(77, "Then the \"rows selected\" reading of line chart viewer should be 36", () => readingIs(page, "rows selected", el("line chart viewer"), 36));
      await session.step(78, "And line chart viewer should show a selection highlight", () => someHighlight(page, el("line chart viewer")));
      await session.step(79, "When user sets \"showSelectedRows\" property of line chart viewer to \"false\"", () => setProperty(page, "showSelectedRows", el("line chart viewer"), "false"));
      await session.step(80, "Then line chart viewer should show less selection highlight than before", () => lessHighlight(page, el("line chart viewer")));
      await session.step(81, "And the \"rows selected\" reading of line chart viewer should be 36", () => readingIs(page, "rows selected", el("line chart viewer"), 36));
      await session.step(82, "When user sets \"showSelectedRows\" property of line chart viewer to \"true\"", () => setProperty(page, "showSelectedRows", el("line chart viewer"), "true"));
      await session.step(83, "Then line chart viewer should show more selection highlight than before", () => moreHighlight(page, el("line chart viewer")));
      await session.step(84, "When user clears the row selection", () => clearSelection(page));
      await session.step(85, "Then the \"rows selected\" reading of line chart viewer should be 0", () => readingIs(page, "rows selected", el("line chart viewer"), 0));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Opening and dismissing the context menu leaves the selection and the current row alone", async () => {
      await session.step(89, "When user selects rows where \"Stereo Category\" is \"S_ABS\"", () => selectWhereIs(page, "Stereo Category", "S_ABS"));
      await session.step(90, "And user makes row 7 current", () => makeRowCurrent(page, 7));
      await session.step(91, "Then the \"rows selected\" reading of line chart viewer should be 2", () => readingIs(page, "rows selected", el("line chart viewer"), 2));
      await session.step(92, "When user opens the context menu of line chart viewer", () => openContextMenu(page, el("line chart viewer")));
      await session.step(93, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(94, "Then the \"rows selected\" reading of line chart viewer should be 2", () => readingIs(page, "rows selected", el("line chart viewer"), 2));
      await session.step(95, "And row 7 should be current", () => currentRowIs(page, 7));
      await session.step(96, "And 2 rows should be selected", () => selectedRowCount(page, 2));
      await session.step(97, "When user clears the row selection", () => clearSelection(page));
      await session.step(98, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The current-row line is drawn only when it is asked for", async () => {
      await session.step(101, "When user makes row 40 current", () => makeRowCurrent(page, 40));
      await session.step(102, "And user sets \"showCurrentRowLine\" property of line chart viewer to \"true\"", () => setProperty(page, "showCurrentRowLine", el("line chart viewer"), "true"));
      await session.step(103, "Then line chart viewer should have more ink than before", () => moreInk(page, el("line chart viewer")));
      await session.step(104, "When user makes row 80 current", () => makeRowCurrent(page, 80));
      await session.step(105, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(106, "When user sets \"showCurrentRowLine\" property of line chart viewer to \"false\"", () => setProperty(page, "showCurrentRowLine", el("line chart viewer"), "false"));
      await session.step(107, "Then line chart viewer should have less ink than before", () => lessInk(page, el("line chart viewer")));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Selecting through the data and reading it off the chart agree", async () => {
      await session.step(111, "When user selects rows where \"Stereo Category\" is one of \"S_ABS, S_PART\"", () => selectWhereOneOf(page, "Stereo Category", "S_ABS, S_PART"));
      await session.step(112, "Then 12 rows should be selected", () => selectedRowCount(page, 12));
      await session.step(113, "And the \"rows selected\" reading of line chart viewer should be 12", () => readingIs(page, "rows selected", el("line chart viewer"), 12));
      await session.step(114, "And line chart viewer should show a selection highlight", () => someHighlight(page, el("line chart viewer")));
      await session.step(115, "When user selects no rows", () => selectNoRows(page));
      await session.step(116, "Then the \"rows selected\" reading of line chart viewer should be 0", () => readingIs(page, "rows selected", el("line chart viewer"), 0));
      await session.step(117, "And line chart viewer should show no selection highlight", () => noHighlight(page, el("line chart viewer")));
      await session.step(118, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("After Tools > Lasso Tool the next Shift-drag selects nothing", async () => {
      await session.step(129, "Given user clears the row selection", () => clearSelection(page));
      await session.step(130, "When user picks \"Tools > Lasso Tool\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Tools > Lasso Tool", el("line chart viewer")));
      await session.step(131, "Then \"lassoTool\" property of line chart viewer should be \"true\"", () => propertyShouldBe(page, "lassoTool", el("line chart viewer"), "true"));
      await session.step(132, "And the \"region drawing mode\" reading of line chart viewer should be \"false\"", () => readingReads(page, "region drawing mode", el("line chart viewer"), "false"));
      await session.step(133, "When user drags a selection box over the \"plot\" area of line chart viewer", () => dragSelectionOverArea(page, "plot", el("line chart viewer")));
      await session.step(134, "Then the \"rows selected\" reading of line chart viewer should be 82", () => readingIs(page, "rows selected", el("line chart viewer"), 82));
      await session.step(135, "When user clears the row selection", () => clearSelection(page));
      await session.step(136, "And user sets \"lassoTool\" property of line chart viewer to \"false\"", () => setProperty(page, "lassoTool", el("line chart viewer"), "false"));
      await session.step(137, "Then no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true});
    run.finish();
  });
});
