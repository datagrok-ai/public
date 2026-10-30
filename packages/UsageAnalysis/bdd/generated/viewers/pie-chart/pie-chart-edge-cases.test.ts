/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart, core#982]
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
import {filterIsExactlyCategory, filterPasses, filterPassesAll, rowCount, setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hasArea, noErrors, painted, readingIs, readingReads, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart edge cases — empty and one-row tables, a click under a Filter formula", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-edge-cases.feature", import.meta.url);
  test("On a one-row table a single wedge holds the whole circle", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "Given user opens demog-1000 dataset keeping the first 1 rows", () => openDatasetRows(page, ds("demog-1000"), 1));
    await session.step(17, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"]]), [["Category","RACE"]]);
    await session.step(19, "Then the table should have 1 row", () => rowCount(page, 1));
    await session.step(20, "And the \"slices\" reading of pie chart viewer should be 1", () => readingIs(page, "slices", el("pie chart viewer"), 1));
    await session.step(21, "And the \"share of Caucasian\" reading of pie chart viewer should be 100", () => readingIs(page, "share of Caucasian", el("pie chart viewer"), 100));
    await session.step(22, "And pie chart viewer should be painted", () => painted(page, el("pie chart viewer")));
    await session.step(23, "And no errors should have been logged", () => noErrors(page));
  });
  test("A Filter formula no row passes empties the pie, and a click on its empty space leaves the table filter alone", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(27, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["On Click","Filter"]]), [["Category","RACE"],["On Click","Filter"]]);
    await session.step(30, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(31, "When user sets \"Filter\" property of pie chart viewer to \"${AGE} > 1000\"", () => setProperty(page, "Filter", el("pie chart viewer"), "${AGE} > 1000"));
    await session.step(32, "Then pie chart viewer should show 0 rows", () => showsRows(page, el("pie chart viewer"), 0));
    await session.step(33, "And the \"slices\" reading of pie chart viewer should be 0", () => readingIs(page, "slices", el("pie chart viewer"), 0));
    await session.step(34, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(35, "When user clicks on the \"view\" area of pie chart viewer", () => clickArea(page, "view", el("pie chart viewer")));
    await session.step(36, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("Under a Filter formula a wedge click in Filter mode keeps exactly the clicked category (core#982)", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart", "@realizes:core#982"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(42, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["On Click","Filter"],["Filter","${AGE} > 40"]]), [["Category","RACE"],["On Click","Filter"],["Filter","${AGE} > 40"]]);
    await session.step(46, "Then pie chart viewer should show 635 rows", () => showsRows(page, el("pie chart viewer"), 635));
    await session.step(47, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 562", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 562));
    await session.step(48, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
    await session.step(49, "Then the filter should pass exactly the rows where \"RACE\" is \"Caucasian\"", () => filterIsExactlyCategory(page, "RACE", "Caucasian"));
    await session.step(50, "And 896 rows should pass the filter", () => filterPasses(page, 896));
    await session.step(51, "When user clicks on the \"empty space\" area of pie chart viewer", () => clickArea(page, "empty space", el("pie chart viewer")));
    await session.step(52, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
  test("An angle column holding NaN leaves that value out, as an empty one [value=NaN]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(56, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(57, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]), [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]);
    await session.step(61, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(62, "And the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(63, "When user sets \"WEIGHT\" column in row 1 to \"NaN\"", () => setCell(page, "WEIGHT", 1, "NaN"));
    await session.step(64, "Then the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(65, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(66, "And pie chart viewer should have a \"pie\" area", () => hasArea(page, el("pie chart viewer"), "pie"));
    await session.step(67, "And no errors should have been logged", () => noErrors(page));
  });
  test("An angle column holding Infinity leaves that value out, as an empty one [value=Infinity]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(56, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(57, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]), [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]);
    await session.step(61, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(62, "And the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(63, "When user sets \"WEIGHT\" column in row 1 to \"Infinity\"", () => setCell(page, "WEIGHT", 1, "Infinity"));
    await session.step(64, "Then the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(65, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(66, "And pie chart viewer should have a \"pie\" area", () => hasArea(page, el("pie chart viewer"), "pie"));
    await session.step(67, "And no errors should have been logged", () => noErrors(page));
  });
  test("An angle column holding -Infinity leaves that value out, as an empty one [value=-Infinity]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(56, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(57, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]), [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]);
    await session.step(61, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(62, "And the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(63, "When user sets \"WEIGHT\" column in row 1 to \"-Infinity\"", () => setCell(page, "WEIGHT", 1, "-Infinity"));
    await session.step(64, "Then the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(65, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(66, "And pie chart viewer should have a \"pie\" area", () => hasArea(page, el("pie chart viewer"), "pie"));
    await session.step(67, "And no errors should have been logged", () => noErrors(page));
  });
});
