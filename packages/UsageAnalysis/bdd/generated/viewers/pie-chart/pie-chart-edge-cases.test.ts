/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart, core#982, core#1200]
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
import {clickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {deleteSelected, filterIsExactlyCategory, filterPasses, filterPassesAll, rowCount, selectAllRows, setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, clickArea, hasArea, noErrors, painted, propertyShouldBe, readingIs, readingReads, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart edge cases — empty and one-row tables, a click under a Filter formula", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-edge-cases.feature", import.meta.url);
  test("Once every row of the table is deleted the pie draws no wedge and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"]]), [["Category","RACE"]]);
    await session.step(20, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(21, "When user selects all rows", () => selectAllRows(page));
    await session.step(22, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(23, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(24, "And pie chart viewer should show 0 rows", () => showsRows(page, el("pie chart viewer"), 0));
    await session.step(25, "And the \"slices\" reading of pie chart viewer should be 0", () => readingIs(page, "slices", el("pie chart viewer"), 0));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
  });
  test("On a one-row table a single wedge holds the whole circle", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "Given user opens demog-1000 dataset keeping the first 1 rows", () => openDatasetRows(page, ds("demog-1000"), 1));
    await session.step(30, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"]]), [["Category","RACE"]]);
    await session.step(32, "Then the table should have 1 row", () => rowCount(page, 1));
    await session.step(33, "And the \"slices\" reading of pie chart viewer should be 1", () => readingIs(page, "slices", el("pie chart viewer"), 1));
    await session.step(34, "And the \"share of Caucasian\" reading of pie chart viewer should be 100", () => readingIs(page, "share of Caucasian", el("pie chart viewer"), 100));
    await session.step(35, "And pie chart viewer should be painted", () => painted(page, el("pie chart viewer")));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
  });
  test("A Filter formula no row passes empties the pie, and a click on its empty space leaves the table filter alone", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(39, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(40, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["On Click","Filter"]]), [["Category","RACE"],["On Click","Filter"]]);
    await session.step(43, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(44, "When user sets \"Filter\" property of pie chart viewer to \"${AGE} > 1000\"", () => setProperty(page, "Filter", el("pie chart viewer"), "${AGE} > 1000"));
    await session.step(45, "Then pie chart viewer should show 0 rows", () => showsRows(page, el("pie chart viewer"), 0));
    await session.step(46, "And the \"slices\" reading of pie chart viewer should be 0", () => readingIs(page, "slices", el("pie chart viewer"), 0));
    await session.step(47, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(48, "When user clicks on the \"view\" area of pie chart viewer", () => clickArea(page, "view", el("pie chart viewer")));
    await session.step(49, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
  });
  test("Under a Filter formula a wedge click in Filter mode keeps exactly the clicked category (core#982)", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart", "@realizes:core#982"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(54, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(55, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["On Click","Filter"],["Filter","${AGE} > 40"]]), [["Category","RACE"],["On Click","Filter"],["Filter","${AGE} > 40"]]);
    await session.step(59, "Then pie chart viewer should show 635 rows", () => showsRows(page, el("pie chart viewer"), 635));
    await session.step(60, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 562", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 562));
    await session.step(61, "When user clicks on the \"slice Caucasian\" area of pie chart viewer", () => clickArea(page, "slice Caucasian", el("pie chart viewer")));
    await session.step(62, "Then the filter should pass exactly the rows where \"RACE\" is \"Caucasian\"", () => filterIsExactlyCategory(page, "RACE", "Caucasian"));
    await session.step(63, "And 896 rows should pass the filter", () => filterPasses(page, 896));
    await session.step(64, "When user clicks on the \"empty space\" area of pie chart viewer", () => clickArea(page, "empty space", el("pie chart viewer")));
    await session.step(65, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
  test("An angle column holding NaN leaves that value out, as an empty one [value=NaN]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(69, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(70, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]), [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]);
    await session.step(74, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(75, "And the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(76, "When user sets \"WEIGHT\" column in row 1 to \"NaN\"", () => setCell(page, "WEIGHT", 1, "NaN"));
    await session.step(77, "Then the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(78, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(79, "And pie chart viewer should have a \"pie\" area", () => hasArea(page, el("pie chart viewer"), "pie"));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
  test("An angle column holding Infinity leaves that value out, as an empty one [value=Infinity]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(69, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(70, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]), [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]);
    await session.step(74, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(75, "And the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(76, "When user sets \"WEIGHT\" column in row 1 to \"Infinity\"", () => setCell(page, "WEIGHT", 1, "Infinity"));
    await session.step(77, "Then the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(78, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(79, "And pie chart viewer should have a \"pie\" area", () => hasArea(page, el("pie chart viewer"), "pie"));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
  test("An angle column holding -Infinity leaves that value out, as an empty one [value=-Infinity]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(69, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(70, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]), [["Category","RACE"],["Segment Angle Column","WEIGHT"],["Segment Angle Aggr Type","sum"]]);
    await session.step(74, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(75, "And the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(76, "When user sets \"WEIGHT\" column in row 1 to \"-Infinity\"", () => setCell(page, "WEIGHT", 1, "-Infinity"));
    await session.step(77, "Then the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
    await session.step(78, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(79, "And pie chart viewer should have a \"pie\" area", () => hasArea(page, el("pie chart viewer"), "pie"));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
  test("Start Angle typed as -30 in the context panel is outside its bounds of 0 to 360 and is not taken [typed=-30]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart", "@known-failure:95", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(89, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(90, "When user adds a pie chart viewer", () => addViewer(page, "pie chart"));
      await session.step(91, "And user clicks on settings icon of pie chart viewer", () => clickOn(page, el("settings icon of pie chart viewer")));
      await session.step(92, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(93, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(94, "When user enters \"-30\" into \"Start Angle\" property in context panel", () => enterInto(page, "-30", el("\"Start Angle\" property in context panel")));
      await session.step(95, "Then \"Start Angle\" property of pie chart viewer should be \"0\"", () => propertyShouldBe(page, "Start Angle", el("pie chart viewer"), "0"));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    }, 95);
  });
  test("Start Angle typed as 400 in the context panel is outside its bounds of 0 to 360 and is not taken [typed=400]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pie-chart", "@known-failure:95", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(89, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(90, "When user adds a pie chart viewer", () => addViewer(page, "pie chart"));
      await session.step(91, "And user clicks on settings icon of pie chart viewer", () => clickOn(page, el("settings icon of pie chart viewer")));
      await session.step(92, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(93, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(94, "When user enters \"400\" into \"Start Angle\" property in context panel", () => enterInto(page, "400", el("\"Start Angle\" property in context panel")));
      await session.step(95, "Then \"Start Angle\" property of pie chart viewer should be \"0\"", () => propertyShouldBe(page, "Start Angle", el("pie chart viewer"), "0"));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    }, 95);
  });
});
