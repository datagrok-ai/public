/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/row-source/row-source-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, viewers.histogram, viewers.pie-chart]
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
import {filterPasses, filterTo, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, errorBalloonText, noErrors, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Row Source edge cases — a Filter formula that is not a condition, a row set nothing is left in", () => {
  const session = feature(test, "features/viewers/row-source/row-source-edge-cases.feature", import.meta.url);
  test("A scatter plot Filter formula that is not a condition is refused with an error balloon and ignored [viewer=scatter plot, column property=xColumnName]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot", "@realizes:viewers.histogram", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"]]), [["xColumnName","AGE"]]);
    await session.step(18, "Then the \"rows shown\" reading of scatter plot viewer should be 1000", () => readingIs(page, "rows shown", el("scatter plot viewer"), 1000));
    await session.step(19, "When user sets \"filter\" property of scatter plot viewer to \"${AGE} + 1\"", () => setProperty(page, "filter", el("scatter plot viewer"), "${AGE} + 1"));
    await session.step(20, "Then an error balloon containing \"Invalid formula\" should have been shown", () => errorBalloonText(page, "Invalid formula"));
    await session.step(21, "And the \"rows shown\" reading of scatter plot viewer should be 1000", () => readingIs(page, "rows shown", el("scatter plot viewer"), 1000));
    await session.step(22, "When user sets \"filter\" property of scatter plot viewer to \"${AGE} > 44\"", () => setProperty(page, "filter", el("scatter plot viewer"), "${AGE} > 44"));
    await session.step(23, "Then the \"rows shown\" reading of scatter plot viewer should be 519", () => readingIs(page, "rows shown", el("scatter plot viewer"), 519));
    await session.step(24, "And no errors should have been logged", () => noErrors(page));
  });
  test("A histogram Filter formula that is not a condition is refused with an error balloon and ignored [viewer=histogram, column property=valueColumnName]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot", "@realizes:viewers.histogram", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "When user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"]]), [["valueColumnName","AGE"]]);
    await session.step(18, "Then the \"rows shown\" reading of histogram viewer should be 1000", () => readingIs(page, "rows shown", el("histogram viewer"), 1000));
    await session.step(19, "When user sets \"filter\" property of histogram viewer to \"${AGE} + 1\"", () => setProperty(page, "filter", el("histogram viewer"), "${AGE} + 1"));
    await session.step(20, "Then an error balloon containing \"Invalid formula\" should have been shown", () => errorBalloonText(page, "Invalid formula"));
    await session.step(21, "And the \"rows shown\" reading of histogram viewer should be 1000", () => readingIs(page, "rows shown", el("histogram viewer"), 1000));
    await session.step(22, "When user sets \"filter\" property of histogram viewer to \"${AGE} > 44\"", () => setProperty(page, "filter", el("histogram viewer"), "${AGE} > 44"));
    await session.step(23, "Then the \"rows shown\" reading of histogram viewer should be 519", () => readingIs(page, "rows shown", el("histogram viewer"), 519));
    await session.step(24, "And no errors should have been logged", () => noErrors(page));
  });
  test("A pie chart Filter formula that is not a condition is refused with an error balloon and ignored [viewer=pie chart, column property=categoryColumnName]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot", "@realizes:viewers.histogram", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","AGE"]]), [["categoryColumnName","AGE"]]);
    await session.step(18, "Then the \"rows shown\" reading of pie chart viewer should be 1000", () => readingIs(page, "rows shown", el("pie chart viewer"), 1000));
    await session.step(19, "When user sets \"filter\" property of pie chart viewer to \"${AGE} + 1\"", () => setProperty(page, "filter", el("pie chart viewer"), "${AGE} + 1"));
    await session.step(20, "Then an error balloon containing \"Invalid formula\" should have been shown", () => errorBalloonText(page, "Invalid formula"));
    await session.step(21, "And the \"rows shown\" reading of pie chart viewer should be 1000", () => readingIs(page, "rows shown", el("pie chart viewer"), 1000));
    await session.step(22, "When user sets \"filter\" property of pie chart viewer to \"${AGE} > 44\"", () => setProperty(page, "filter", el("pie chart viewer"), "${AGE} > 44"));
    await session.step(23, "Then the \"rows shown\" reading of pie chart viewer should be 519", () => readingIs(page, "rows shown", el("pie chart viewer"), 519));
    await session.step(24, "And no errors should have been logged", () => noErrors(page));
  });
  test("Under FilteredSelected with every selected row filtered out the viewer draws no row and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot", "@realizes:viewers.histogram", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["rowSource","FilteredSelected"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["rowSource","FilteredSelected"]]);
    await session.step(37, "And user selects rows where \"SEX\" is \"M\"", () => selectWhereIs(page, "SEX", "M"));
    await session.step(38, "Then 447 rows should be selected", () => selectedRowCount(page, 447));
    await session.step(39, "And the \"rows shown\" reading of scatter plot viewer should be 447", () => readingIs(page, "rows shown", el("scatter plot viewer"), 447));
    await session.step(40, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
    await session.step(41, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
    await session.step(42, "And the \"rows shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "rows shown", el("scatter plot viewer"), 0));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
});
