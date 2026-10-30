/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart, core#1100]
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {deleteSelected, filterBetween, filterPasses, noneSelected, removeColumn, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, dragSelectionOverArea, hoverArea, noErrors, painted, pointerAway, propertyShouldBe, readingBetween, readingIs, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart edge cases — a log Y axis below zero, an empty row set, a removed column, SPC under a filter", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-edge-cases.feature", import.meta.url);
  test("A logarithmic Y axis over a column that starts below zero drops the non-positive points", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
    await session.step(19, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(20, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(21, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(24, "When user sets \"yAxisType\" property of line chart viewer to \"logarithmic\"", () => setProperty(page, "yAxisType", el("line chart viewer"), "logarithmic"));
    await session.step(25, "Then the \"markers drawn\" reading of line chart viewer should be 59", () => readingIs(page, "markers drawn", el("line chart viewer"), 59));
    await session.step(26, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(27, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a Shift-drag over the plot selects nothing and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
    await session.step(19, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(20, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(21, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(31, "When user filters rows where \"CAST Idea ID\" is between 1 and 2", () => filterBetween(page, "CAST Idea ID", 1, 2));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of line chart viewer should be 0", () => readingIs(page, "rows shown", el("line chart viewer"), 0));
    await session.step(34, "And the \"markers drawn\" reading of line chart viewer should be 0", () => readingIs(page, "markers drawn", el("line chart viewer"), 0));
    await session.step(35, "When user drags a selection box over the \"plot\" area of line chart viewer", () => dragSelectionOverArea(page, "plot", el("line chart viewer")));
    await session.step(36, "Then no rows should be selected", () => noneSelected(page));
    await session.step(37, "And the \"rows selected\" reading of line chart viewer should be 0", () => readingIs(page, "rows selected", el("line chart viewer"), 0));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("Once every row of the table is deleted the chart draws no point and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
    await session.step(19, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(20, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(21, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(41, "When user selects all rows", () => selectAllRows(page));
    await session.step(42, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(43, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(44, "And the \"rows shown\" reading of line chart viewer should be 0", () => readingIs(page, "rows shown", el("line chart viewer"), 0));
    await session.step(45, "And the \"markers drawn\" reading of line chart viewer should be 0", () => readingIs(page, "markers drawn", el("line chart viewer"), 0));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
  test("The Y column removed from the table while the chart shows it raises no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
    await session.step(19, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(20, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(21, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(49, "When user removes \"Chemical Space X\" column", () => removeColumn(page, "Chemical Space X"));
    await session.step(50, "Then the table should not have a column \"Chemical Space X\"", () => hasNoColumn(page, "Chemical Space X"));
    await session.step(51, "When user hovers over the \"view\" area of line chart viewer", () => hoverArea(page, "view", el("line chart viewer")));
    await session.step(52, "And user moves the pointer away from line chart viewer", () => pointerAway(page, el("line chart viewer")));
    await session.step(53, "Then no errors should have been logged", () => noErrors(page));
  });
  test("SPC computes its process average over the rows that pass the filter (core#1100)", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.line-chart", "@realizes:core#1100"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
    await session.step(19, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(20, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(21, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(60, "When user sets \"showStatisticalProcessControl\" property of line chart viewer to \"true\"", () => setProperty(page, "showStatisticalProcessControl", el("line chart viewer"), "true"));
    await session.step(61, "Then the \"spc average\" reading of line chart viewer should be between 3.4 and 3.6", () => readingBetween(page, "spc average", el("line chart viewer"), 3.4, 3.6));
    await session.step(62, "When user filters rows where \"CAST Idea ID\" is between 634783 and 634834", () => filterBetween(page, "CAST Idea ID", 634783, 634834));
    await session.step(63, "Then 50 rows should pass the filter", () => filterPasses(page, 50));
    await session.step(64, "And the \"spc average\" reading of line chart viewer should be between 2.4 and 2.6", () => readingBetween(page, "spc average", el("line chart viewer"), 2.4, 2.6));
  });
  test("X Min and X Max outside every X value are kept and the data stays in view", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
    await session.step(19, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(20, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(21, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(67, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["X Min","700000"],["X Max","800000"]]), [["X Min","700000"],["X Max","800000"]]);
    await session.step(70, "Then \"X Min\" property of line chart viewer should be \"700000\"", () => propertyShouldBe(page, "X Min", el("line chart viewer"), "700000"));
    await session.step(71, "And \"X Max\" property of line chart viewer should be \"800000\"", () => propertyShouldBe(page, "X Max", el("line chart viewer"), "800000"));
    await session.step(72, "And the \"markers drawn\" reading of line chart viewer should be 100", () => readingIs(page, "markers drawn", el("line chart viewer"), 100));
    await session.step(73, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(74, "And no errors should have been logged", () => noErrors(page));
  });
});
