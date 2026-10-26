/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot, core#1057]
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
import {clickOn, enterInto, isExpanded, pressKeyIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, filterBetween, filterPasses, filterPassesAll, renameColumn, resetFilter, selectedPassFilter, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, boundTable, clickArea, hasArea, noErrors, propertyShouldBe, readingIs, readingReads, reportsNoError, setProperties, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot edge cases: no rows, a renamed split column, Tiles Per Row 0, keys after a table swap, a label click under a filter", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-edge-cases.feature", import.meta.url);
  test("A filter that removes every row leaves no row and no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(22, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await session.step(23, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await session.step(26, "When user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(27, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(28, "And the \"rows shown\" reading of trellis plot viewer should be 0", () => readingIs(page, "rows shown", el("trellis plot viewer"), 0));
    await session.step(29, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
    await session.step(30, "When user resets the filter", () => resetFilter(page));
    await session.step(31, "Then the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await session.step(32, "And no errors should have been logged", () => noErrors(page));
  });
  test("A split column renamed is followed by the split", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(22, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await session.step(23, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await session.step(35, "When user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(36, "Then \"Y Column Names\" property of trellis plot viewer should be \"ORIGIN\"", () => propertyShouldBe(page, "Y Column Names", el("trellis plot viewer"), "ORIGIN"));
    await session.step(37, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await session.step(38, "And trellis plot viewer should have a \"cell F | Caucasian\" area", () => hasArea(page, el("trellis plot viewer"), "cell F | Caucasian"));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("Tiles Per Row 0 typed in the context panel on a one-column trellis logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.trellis-plot", "@known-failure:52", "@realizes:core#1057"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(22, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await session.step(23, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await knownFailure(async () => {
      await session.step(45, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Y Column Names",""]]), [["Y Column Names",""]]);
      await session.step(47, "Then the \"y categories\" reading of trellis plot viewer should be 1", () => readingIs(page, "y categories", el("trellis plot viewer"), 1));
      await session.step(48, "When user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
      await session.step(49, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(50, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(51, "When user enters \"0\" into \"Tiles Per Row\" property in context panel", () => enterInto(page, "0", el("\"Tiles Per Row\" property in context panel")));
      await session.step(52, "Then no errors should have been logged", () => noErrors(page));
      await session.step(53, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
    }, 52);
  });
  test("The arrow keys still move the current cell after the trellis is set to another table and back", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.trellis-plot", "@known-failure:73", "@realizes:core#1057"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(22, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await session.step(23, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await knownFailure(async () => {
      await session.step(60, "Given user opens cars dataset", () => openDataset(page, ds("cars")));
      await session.step(61, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(62, "When user sets \"table\" property of trellis plot viewer to \"cars\"", () => setProperty(page, "table", el("trellis plot viewer"), "cars"));
      await session.step(63, "Then trellis plot viewer should be bound to table \"cars\"", () => boundTable(page, el("trellis plot viewer"), "cars"));
      await session.step(64, "When user sets \"table\" property of trellis plot viewer to \"demog-1000\"", () => setProperty(page, "table", el("trellis plot viewer"), "demog-1000"));
      await session.step(65, "And user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
      await session.step(69, "Then trellis plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("trellis plot viewer"), "demog-1000"));
      await session.step(70, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(71, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(72, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(73, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Caucasian"));
    }, 73);
  });
  test("A category label click under a filter selects only the rows that pass the filter", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.trellis-plot", "@known-failure:83", "@realizes:core#1057"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(22, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await session.step(23, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await knownFailure(async () => {
      await session.step(79, "Given user clears the row selection", () => clearSelection(page));
      await session.step(80, "When user filters rows where \"AGE\" is between 20 and 40", () => filterBetween(page, "AGE", 20, 40));
      await session.step(81, "And user clicks on the \"x label F\" area of trellis plot viewer", () => clickArea(page, "x label F", el("trellis plot viewer")));
      await session.step(82, "Then some rows should be selected", () => someSelected(page));
      await session.step(83, "And every selected row should pass the filter", () => selectedPassFilter(page));
    }, 83);
  });
  test("Closing the trellis gives back the rows its clicked cell filtered out", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(22, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await session.step(23, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await session.step(86, "When user sets \"On Click\" property of trellis plot viewer to \"Filter\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Filter"));
    await session.step(87, "And user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
    await session.step(88, "Then 480 rows should pass the filter", () => filterPasses(page, 480));
    await session.step(89, "When user clicks on close icon of trellis plot viewer", () => clickOn(page, el("close icon of trellis plot viewer")));
    await session.step(90, "Then the open tableview should have 0 trellis plot viewers", () => viewerCount(page, 0, "trellis plot"));
    await session.step(91, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(92, "And no errors should have been logged", () => noErrors(page));
  });
});
