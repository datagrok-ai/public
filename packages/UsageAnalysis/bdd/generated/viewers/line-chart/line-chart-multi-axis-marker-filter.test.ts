/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-multi-axis-marker-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickLegendItem, clickLegendItemHolding, noErrors, readingAsRemembered, readingIs, rememberReading, reportsNoError} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart in Multi Axis mode filtered from its marker legend", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-multi-axis-marker-filter.feature", import.meta.url);
  test("One marker clicked in the legend keeps a marker for every Y column at the heights the unfiltered chart drew it", {tag: ["@viewers", "@realizes:viewers.line-chart", "@core-299"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","STARTED"],["xMap","quarter"],["yColumnNames","AGE, HEIGHT, WEIGHT"],["markersColumnName","STARTED"],["markersMap","quarter"],["multiAxis","true"]]), [["xColumnName","STARTED"],["xMap","quarter"],["yColumnNames","AGE, HEIGHT, WEIGHT"],["markersColumnName","STARTED"],["markersMap","quarter"],["multiAxis","true"]]);
    await session.step(22, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
    await session.step(23, "And the \"rows shown\" reading of line chart viewer should be 1000", () => readingIs(page, "rows shown", el("line chart viewer"), 1000));
    await session.step(24, "And the \"x categories\" reading of line chart viewer should be 4", () => readingIs(page, "x categories", el("line chart viewer"), 4));
    await session.step(25, "And the \"lines\" reading of line chart viewer should be 3", () => readingIs(page, "lines", el("line chart viewer"), 3));
    await session.step(26, "And the \"markers drawn\" reading of line chart viewer should be 12", () => readingIs(page, "markers drawn", el("line chart viewer"), 12));
    await session.step(27, "When user remembers the 'y axis min of \"AGE\"' reading of line chart viewer", () => rememberReading(page, "y axis min of \"AGE\"", el("line chart viewer")));
    await session.step(28, "And user remembers the 'y axis max of \"AGE\"' reading of line chart viewer", () => rememberReading(page, "y axis max of \"AGE\"", el("line chart viewer")));
    await session.step(29, "And user remembers the 'y axis min of \"HEIGHT\"' reading of line chart viewer", () => rememberReading(page, "y axis min of \"HEIGHT\"", el("line chart viewer")));
    await session.step(30, "And user remembers the 'y axis max of \"HEIGHT\"' reading of line chart viewer", () => rememberReading(page, "y axis max of \"HEIGHT\"", el("line chart viewer")));
    await session.step(31, "And user remembers the 'y axis min of \"WEIGHT\"' reading of line chart viewer", () => rememberReading(page, "y axis min of \"WEIGHT\"", el("line chart viewer")));
    await session.step(32, "And user remembers the 'y axis max of \"WEIGHT\"' reading of line chart viewer", () => rememberReading(page, "y axis max of \"WEIGHT\"", el("line chart viewer")));
    await session.step(35, "When user clicks on \"Q3\" item in the legend of line chart viewer", () => clickLegendItem(page, "Q3", el("line chart viewer")));
    await session.step(36, "Then the \"rows shown\" reading of line chart viewer should be 259", () => readingIs(page, "rows shown", el("line chart viewer"), 259));
    await session.step(37, "And the \"x categories\" reading of line chart viewer should be 1", () => readingIs(page, "x categories", el("line chart viewer"), 1));
    await session.step(38, "And the \"lines\" reading of line chart viewer should be 3", () => readingIs(page, "lines", el("line chart viewer"), 3));
    await session.step(39, "And the \"markers drawn\" reading of line chart viewer should be 3", () => readingIs(page, "markers drawn", el("line chart viewer"), 3));
    await session.step(40, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(41, "And the 'y axis min of \"AGE\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis min of \"AGE\"", el("line chart viewer")));
    await session.step(42, "And the 'y axis max of \"AGE\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis max of \"AGE\"", el("line chart viewer")));
    await session.step(43, "And the 'y axis min of \"HEIGHT\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis min of \"HEIGHT\"", el("line chart viewer")));
    await session.step(44, "And the 'y axis max of \"HEIGHT\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis max of \"HEIGHT\"", el("line chart viewer")));
    await session.step(45, "And the 'y axis min of \"WEIGHT\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis min of \"WEIGHT\"", el("line chart viewer")));
    await session.step(46, "And the 'y axis max of \"WEIGHT\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis max of \"WEIGHT\"", el("line chart viewer")));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
  test("Two markers clicked in the legend keep every Y column's line where it was", {tag: ["@viewers", "@realizes:viewers.line-chart", "@core-299"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","STARTED"],["xMap","quarter"],["yColumnNames","AGE, HEIGHT, WEIGHT"],["markersColumnName","STARTED"],["markersMap","quarter"],["multiAxis","true"]]), [["xColumnName","STARTED"],["xMap","quarter"],["yColumnNames","AGE, HEIGHT, WEIGHT"],["markersColumnName","STARTED"],["markersMap","quarter"],["multiAxis","true"]]);
    await session.step(22, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
    await session.step(23, "And the \"rows shown\" reading of line chart viewer should be 1000", () => readingIs(page, "rows shown", el("line chart viewer"), 1000));
    await session.step(24, "And the \"x categories\" reading of line chart viewer should be 4", () => readingIs(page, "x categories", el("line chart viewer"), 4));
    await session.step(25, "And the \"lines\" reading of line chart viewer should be 3", () => readingIs(page, "lines", el("line chart viewer"), 3));
    await session.step(26, "And the \"markers drawn\" reading of line chart viewer should be 12", () => readingIs(page, "markers drawn", el("line chart viewer"), 12));
    await session.step(27, "When user remembers the 'y axis min of \"AGE\"' reading of line chart viewer", () => rememberReading(page, "y axis min of \"AGE\"", el("line chart viewer")));
    await session.step(28, "And user remembers the 'y axis max of \"AGE\"' reading of line chart viewer", () => rememberReading(page, "y axis max of \"AGE\"", el("line chart viewer")));
    await session.step(29, "And user remembers the 'y axis min of \"HEIGHT\"' reading of line chart viewer", () => rememberReading(page, "y axis min of \"HEIGHT\"", el("line chart viewer")));
    await session.step(30, "And user remembers the 'y axis max of \"HEIGHT\"' reading of line chart viewer", () => rememberReading(page, "y axis max of \"HEIGHT\"", el("line chart viewer")));
    await session.step(31, "And user remembers the 'y axis min of \"WEIGHT\"' reading of line chart viewer", () => rememberReading(page, "y axis min of \"WEIGHT\"", el("line chart viewer")));
    await session.step(32, "And user remembers the 'y axis max of \"WEIGHT\"' reading of line chart viewer", () => rememberReading(page, "y axis max of \"WEIGHT\"", el("line chart viewer")));
    await session.step(50, "When user clicks on \"Q2\" item in the legend of line chart viewer", () => clickLegendItem(page, "Q2", el("line chart viewer")));
    await session.step(51, "Then the \"rows shown\" reading of line chart viewer should be 254", () => readingIs(page, "rows shown", el("line chart viewer"), 254));
    await session.step(52, "When user clicks on \"Q3\" item in the legend of line chart viewer holding Control", () => clickLegendItemHolding(page, "Q3", el("line chart viewer"), "Control"));
    await session.step(53, "Then the \"rows shown\" reading of line chart viewer should be 513", () => readingIs(page, "rows shown", el("line chart viewer"), 513));
    await session.step(54, "And the \"x categories\" reading of line chart viewer should be 2", () => readingIs(page, "x categories", el("line chart viewer"), 2));
    await session.step(55, "And the \"lines\" reading of line chart viewer should be 3", () => readingIs(page, "lines", el("line chart viewer"), 3));
    await session.step(56, "And the \"markers drawn\" reading of line chart viewer should be 6", () => readingIs(page, "markers drawn", el("line chart viewer"), 6));
    await session.step(57, "And the 'y axis min of \"AGE\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis min of \"AGE\"", el("line chart viewer")));
    await session.step(58, "And the 'y axis max of \"AGE\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis max of \"AGE\"", el("line chart viewer")));
    await session.step(59, "And the 'y axis min of \"HEIGHT\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis min of \"HEIGHT\"", el("line chart viewer")));
    await session.step(60, "And the 'y axis max of \"HEIGHT\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis max of \"HEIGHT\"", el("line chart viewer")));
    await session.step(61, "And the 'y axis min of \"WEIGHT\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis min of \"WEIGHT\"", el("line chart viewer")));
    await session.step(62, "And the 'y axis max of \"WEIGHT\"' reading of line chart viewer should be as remembered", () => readingAsRemembered(page, "y axis max of \"WEIGHT\"", el("line chart viewer")));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
  });
});
