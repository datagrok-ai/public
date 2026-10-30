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
import {clearSelection, filterBetween, filterPasses, filterPassesAll, selectedPassFilter, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, boundTable, clickArea, noErrors, readingIs, readingReads, reportsNoError, setProperties, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot edge cases: Tiles Per Row 0, keys after a table swap, a label click under a filter", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-edge-cases.feature", import.meta.url);
  test("Tiles Per Row 0 typed in the context panel on a one-column trellis logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.trellis-plot", "@realizes:core#1057"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(19, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await session.step(20, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await session.step(26, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Y Column Names",""]]), [["Y Column Names",""]]);
    await session.step(28, "Then the \"y categories\" reading of trellis plot viewer should be 1", () => readingIs(page, "y categories", el("trellis plot viewer"), 1));
    await session.step(29, "When user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
    await session.step(30, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(31, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
    await session.step(32, "When user enters \"0\" into \"Tiles Per Row\" property in context panel", () => enterInto(page, "0", el("\"Tiles Per Row\" property in context panel")));
    await session.step(33, "Then no errors should have been logged", () => noErrors(page));
    await session.step(34, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
  });
  test("The arrow keys still move the current cell after the trellis is set to another table and back", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.trellis-plot", "@realizes:core#1057"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(19, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await session.step(20, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await session.step(41, "Given user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(42, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
    await session.step(43, "When user sets \"table\" property of trellis plot viewer to \"cars\"", () => setProperty(page, "table", el("trellis plot viewer"), "cars"));
    await session.step(44, "Then trellis plot viewer should be bound to table \"cars\"", () => boundTable(page, el("trellis plot viewer"), "cars"));
    await session.step(45, "When user sets \"table\" property of trellis plot viewer to \"demog-1000\"", () => setProperty(page, "table", el("trellis plot viewer"), "demog-1000"));
    await session.step(46, "And user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(50, "Then trellis plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("trellis plot viewer"), "demog-1000"));
    await session.step(51, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
    await session.step(52, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
    await session.step(53, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
    await session.step(54, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Caucasian"));
  });
  test("A category label click under a filter selects only the rows that pass the filter", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.trellis-plot", "@realizes:core#1057"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(19, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await session.step(20, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await session.step(60, "Given user clears the row selection", () => clearSelection(page));
    await session.step(61, "When user filters rows where \"AGE\" is between 20 and 40", () => filterBetween(page, "AGE", 20, 40));
    await session.step(62, "And user clicks on the \"x label F\" area of trellis plot viewer", () => clickArea(page, "x label F", el("trellis plot viewer")));
    await session.step(63, "Then some rows should be selected", () => someSelected(page));
    await session.step(64, "And every selected row should pass the filter", () => selectedPassFilter(page));
  });
  test("Closing the trellis gives back the rows its clicked cell filtered out", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(19, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await session.step(20, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await session.step(67, "When user sets \"On Click\" property of trellis plot viewer to \"Filter\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Filter"));
    await session.step(68, "And user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
    await session.step(69, "Then 480 rows should pass the filter", () => filterPasses(page, 480));
    await session.step(70, "When user clicks on close icon of trellis plot viewer", () => clickOn(page, el("close icon of trellis plot viewer")));
    await session.step(71, "Then the open tableview should have 0 trellis plot viewers", () => viewerCount(page, 0, "trellis plot"));
    await session.step(72, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(73, "And no errors should have been logged", () => noErrors(page));
  });
});
