/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart, core#974]
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
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {deleteSelected, filterBetween, filterPasses, filterPassesAll, filterPassesFewer, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, clickArea, hoverArea, noErrors, readingIs, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickEmptySpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart edge cases", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-edge-cases.feature", import.meta.url);
  test("A vertical bar chart over a table with every row deleted draws no bars and logs no error", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Orientation","vertical"]]), [["Split","RACE"],["Orientation","vertical"]]);
    await session.step(21, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(22, "When user selects all rows", () => selectAllRows(page));
    await session.step(23, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(24, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(25, "And the \"bars\" reading of bar chart viewer should be 0", () => readingIs(page, "bars", el("bar chart viewer"), 0));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
  });
  test("A negative Bar Corner Radius draws the bars and logs no error", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(31, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Stack","SEX"]]), [["Split","RACE"],["Stack","SEX"]]);
    await session.step(34, "When user sets \"Bar Corner Radius\" property of bar chart viewer to \"-5\"", () => setProperty(page, "Bar Corner Radius", el("bar chart viewer"), "-5"));
    await session.step(35, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(36, "And the \"bar Caucasian\" area of bar chart viewer should be painted", () => areaPainted(page, "bar Caucasian", el("bar chart viewer")));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("The tooltip shows a category with markup characters as typed", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "Given user opens a table \"markup-categories\" with:", () => openTableOf(page, "markup-categories", [["Category","Value"],["a<b","1"],["a<b","2"],["x & y","3"]]), [["Category","Value"],["a<b","1"],["a<b","2"],["x & y","3"]]);
    await session.step(46, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Category"]]), [["Split","Category"]]);
    await session.step(48, "Then the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
    await session.step(49, "When user hovers over the \"bar a<b\" area of bar chart viewer", () => hoverArea(page, "bar a<b", el("bar chart viewer")));
    await session.step(50, "Then tooltip should contain text \"a<b\"", () => shouldContainText(page, el("tooltip"), "a<b"));
    await session.step(51, "And tooltip should contain text \"2 rows\"", () => shouldContainText(page, el("tooltip"), "2 rows"));
    await session.step(52, "When user hovers over the \"bar x & y\" area of bar chart viewer", () => hoverArea(page, "bar x & y", el("bar chart viewer")));
    await session.step(53, "Then tooltip should contain text \"x & y\"", () => shouldContainText(page, el("tooltip"), "x & y"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
  });
  test("A Filter click on a chart stacked by a date column is released after an unrelated setting changes", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(58, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(59, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Stack","STARTED"],["On Click","Filter"]]), [["Split","RACE"],["Stack","STARTED"],["On Click","Filter"]]);
    await session.step(63, "When user clicks on the \"bar Asian\" area of bar chart viewer", () => clickArea(page, "bar Asian", el("bar chart viewer")));
    await session.step(64, "Then 6 rows should pass the filter", () => filterPasses(page, 6));
    await session.step(65, "When user sets \"Bar Corner Radius\" property of bar chart viewer to \"4\"", () => setProperty(page, "Bar Corner Radius", el("bar chart viewer"), "4"));
    await session.step(66, "And user clicks on empty plot space of bar chart viewer", () => clickEmptySpace(page, el("bar chart viewer")));
    await session.step(67, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
  });
  test("With Row Source All and Show Filtered Rows every bar stays when every row is filtered out", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(72, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(73, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Row Source","All"],["Show Filtered Rows","true"]]), [["Split","RACE"],["Row Source","All"],["Show Filtered Rows","true"]]);
    await session.step(77, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(78, "When user filters rows where \"AGE\" is between 200 and 300", () => filterBetween(page, "AGE", 200, 300));
    await session.step(79, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(80, "And the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(81, "And no errors should have been logged", () => noErrors(page));
  });
  test("Closing the bar chart gives back the rows its clicked bar filtered out", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(85, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(86, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["On Click","Filter"]]), [["Split","RACE"],["On Click","Filter"]]);
    await session.step(89, "When user clicks on the \"bar Asian\" area of bar chart viewer", () => clickArea(page, "bar Asian", el("bar chart viewer")));
    await session.step(90, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
    await session.step(91, "When user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
    await session.step(92, "Then the open tableview should have 0 bar chart viewers", () => viewerCount(page, 0, "bar chart"));
    await session.step(93, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(94, "And no errors should have been logged", () => noErrors(page));
  });
});
