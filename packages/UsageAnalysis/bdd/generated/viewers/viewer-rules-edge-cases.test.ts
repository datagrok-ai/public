/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/viewer-rules-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [core#1213, public#4452, core#233, core#1174, core#1200, core#1217]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterBetween, filterPasses, filterPassesAll, filterPassesFewer, removeColumn, renameColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, clickArea, loadLayout, noErrors, propertyShouldBe, propertyShouldNotBe, readingAtLeast, readingIs, reportsNoError, saveLayout, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Viewer edge cases shared by every viewer: no rows, a column removed or renamed, a layout with a missing column, a filter left by a closed viewer, values out of bounds, empty categories", () => {
  const session = feature(test, "features/viewers/viewer-rules-edge-cases.feature", import.meta.url);
  test("With every row filtered out a bar chart shows no row and no message, and the rows come back [viewer=bar chart, property=splitColumnName, value=RACE, drawn=rows shown]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"]]), [["splitColumnName","RACE"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of bar chart viewer should be 0", () => readingIs(page, "rows shown", el("bar chart viewer"), 0));
    await session.step(34, "And the \"rows shown\" reading of bar chart viewer should be 0", () => readingIs(page, "rows shown", el("bar chart viewer"), 0));
    await session.step(35, "And bar chart viewer should report no error", () => reportsNoError(page, el("bar chart viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of bar chart viewer should be 1000", () => readingIs(page, "rows shown", el("bar chart viewer"), 1000));
    await session.step(38, "And the \"rows shown\" reading of bar chart viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("bar chart viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a box plot shows no row and no message, and the rows come back [viewer=box plot, property=valueColumnName, value=AGE, drawn=rows shown]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["valueColumnName","AGE"]]), [["valueColumnName","AGE"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of box plot viewer should be 0", () => readingIs(page, "rows shown", el("box plot viewer"), 0));
    await session.step(34, "And the \"rows shown\" reading of box plot viewer should be 0", () => readingIs(page, "rows shown", el("box plot viewer"), 0));
    await session.step(35, "And box plot viewer should report no error", () => reportsNoError(page, el("box plot viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of box plot viewer should be 1000", () => readingIs(page, "rows shown", el("box plot viewer"), 1000));
    await session.step(38, "And the \"rows shown\" reading of box plot viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("box plot viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a heat map shows no row and no message, and the rows come back [viewer=heat map, property=rowSource, value=Filtered, drawn=rows shown]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a heat map viewer with:", () => addViewerWith(page, "heat map", [["rowSource","Filtered"]]), [["rowSource","Filtered"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of heat map viewer should be 0", () => readingIs(page, "rows shown", el("heat map viewer"), 0));
    await session.step(34, "And the \"rows shown\" reading of heat map viewer should be 0", () => readingIs(page, "rows shown", el("heat map viewer"), 0));
    await session.step(35, "And heat map viewer should report no error", () => reportsNoError(page, el("heat map viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
    await session.step(38, "And the \"rows shown\" reading of heat map viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("heat map viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a histogram shows no row and no message, and the rows come back [viewer=histogram, property=valueColumnName, value=AGE, drawn=rows shown]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"]]), [["valueColumnName","AGE"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of histogram viewer should be 0", () => readingIs(page, "rows shown", el("histogram viewer"), 0));
    await session.step(34, "And the \"rows shown\" reading of histogram viewer should be 0", () => readingIs(page, "rows shown", el("histogram viewer"), 0));
    await session.step(35, "And histogram viewer should report no error", () => reportsNoError(page, el("histogram viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of histogram viewer should be 1000", () => readingIs(page, "rows shown", el("histogram viewer"), 1000));
    await session.step(38, "And the \"rows shown\" reading of histogram viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("histogram viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a matrix plot shows no row and no message, and the rows come back [viewer=matrix plot, property=rowSource, value=Filtered, drawn=rows shown]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a matrix plot viewer with:", () => addViewerWith(page, "matrix plot", [["rowSource","Filtered"]]), [["rowSource","Filtered"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of matrix plot viewer should be 0", () => readingIs(page, "rows shown", el("matrix plot viewer"), 0));
    await session.step(34, "And the \"rows shown\" reading of matrix plot viewer should be 0", () => readingIs(page, "rows shown", el("matrix plot viewer"), 0));
    await session.step(35, "And matrix plot viewer should report no error", () => reportsNoError(page, el("matrix plot viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of matrix plot viewer should be 1000", () => readingIs(page, "rows shown", el("matrix plot viewer"), 1000));
    await session.step(38, "And the \"rows shown\" reading of matrix plot viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("matrix plot viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a network diagram shows no row and no message, and the rows come back [viewer=network diagram, property=node1ColumnName, value=RACE, drawn=rows shown]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a network diagram viewer with:", () => addViewerWith(page, "network diagram", [["node1ColumnName","RACE"]]), [["node1ColumnName","RACE"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of network diagram viewer should be 0", () => readingIs(page, "rows shown", el("network diagram viewer"), 0));
    await session.step(34, "And the \"rows shown\" reading of network diagram viewer should be 0", () => readingIs(page, "rows shown", el("network diagram viewer"), 0));
    await session.step(35, "And network diagram viewer should report no error", () => reportsNoError(page, el("network diagram viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of network diagram viewer should be 1000", () => readingIs(page, "rows shown", el("network diagram viewer"), 1000));
    await session.step(38, "And the \"rows shown\" reading of network diagram viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("network diagram viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a pc plot shows no row and no message, and the rows come back [viewer=pc plot, property=columnNames, value=AGE, HEIGHT, WEIGHT, drawn=lines drawn]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["columnNames","AGE, HEIGHT, WEIGHT"]]), [["columnNames","AGE, HEIGHT, WEIGHT"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of pc plot viewer should be 0", () => readingIs(page, "rows shown", el("pc plot viewer"), 0));
    await session.step(34, "And the \"lines drawn\" reading of pc plot viewer should be 0", () => readingIs(page, "lines drawn", el("pc plot viewer"), 0));
    await session.step(35, "And pc plot viewer should report no error", () => reportsNoError(page, el("pc plot viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of pc plot viewer should be 1000", () => readingIs(page, "rows shown", el("pc plot viewer"), 1000));
    await session.step(38, "And the \"lines drawn\" reading of pc plot viewer should be at least 1", () => readingAtLeast(page, "lines drawn", el("pc plot viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a pie chart shows no row and no message, and the rows come back [viewer=pie chart, property=categoryColumnName, value=RACE, drawn=slices]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","RACE"]]), [["categoryColumnName","RACE"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of pie chart viewer should be 0", () => readingIs(page, "rows shown", el("pie chart viewer"), 0));
    await session.step(34, "And the \"slices\" reading of pie chart viewer should be 0", () => readingIs(page, "slices", el("pie chart viewer"), 0));
    await session.step(35, "And pie chart viewer should report no error", () => reportsNoError(page, el("pie chart viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of pie chart viewer should be 1000", () => readingIs(page, "rows shown", el("pie chart viewer"), 1000));
    await session.step(38, "And the \"slices\" reading of pie chart viewer should be at least 1", () => readingAtLeast(page, "slices", el("pie chart viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a pivot table shows no row and no message, and the rows come back [viewer=pivot table, property=rowSource, value=Filtered, drawn=aggregated rows]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a pivot table viewer with:", () => addViewerWith(page, "pivot table", [["rowSource","Filtered"]]), [["rowSource","Filtered"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of pivot table viewer should be 0", () => readingIs(page, "rows shown", el("pivot table viewer"), 0));
    await session.step(34, "And the \"aggregated rows\" reading of pivot table viewer should be 0", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 0));
    await session.step(35, "And pivot table viewer should report no error", () => reportsNoError(page, el("pivot table viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
    await session.step(38, "And the \"aggregated rows\" reading of pivot table viewer should be at least 1", () => readingAtLeast(page, "aggregated rows", el("pivot table viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a scatter plot shows no row and no message, and the rows come back [viewer=scatter plot, property=xColumnName, value=AGE, drawn=rows shown]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"]]), [["xColumnName","AGE"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "rows shown", el("scatter plot viewer"), 0));
    await session.step(34, "And the \"rows shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "rows shown", el("scatter plot viewer"), 0));
    await session.step(35, "And scatter plot viewer should report no error", () => reportsNoError(page, el("scatter plot viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of scatter plot viewer should be 1000", () => readingIs(page, "rows shown", el("scatter plot viewer"), 1000));
    await session.step(38, "And the \"rows shown\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("scatter plot viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a tile viewer shows no row and no message, and the rows come back [viewer=tile viewer, property=lanesColumnName, value=RACE, drawn=tiles]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a tile viewer viewer with:", () => addViewerWith(page, "tile viewer", [["lanesColumnName","RACE"]]), [["lanesColumnName","RACE"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of tile viewer viewer should be 0", () => readingIs(page, "rows shown", el("tile viewer viewer"), 0));
    await session.step(34, "And the \"tiles\" reading of tile viewer viewer should be 0", () => readingIs(page, "tiles", el("tile viewer viewer"), 0));
    await session.step(35, "And tile viewer viewer should report no error", () => reportsNoError(page, el("tile viewer viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of tile viewer viewer should be 1000", () => readingIs(page, "rows shown", el("tile viewer viewer"), 1000));
    await session.step(38, "And the \"tiles\" reading of tile viewer viewer should be at least 1", () => readingAtLeast(page, "tiles", el("tile viewer viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a tree map shows no row and no message, and the rows come back [viewer=tree map, property=splitByColumnNames, value=RACE, SEX, drawn=leaves]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE, SEX"]]), [["splitByColumnNames","RACE, SEX"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of tree map viewer should be 0", () => readingIs(page, "rows shown", el("tree map viewer"), 0));
    await session.step(34, "And the \"leaves\" reading of tree map viewer should be 0", () => readingIs(page, "leaves", el("tree map viewer"), 0));
    await session.step(35, "And tree map viewer should report no error", () => reportsNoError(page, el("tree map viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await session.step(38, "And the \"leaves\" reading of tree map viewer should be at least 1", () => readingAtLeast(page, "leaves", el("tree map viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out a trellis plot shows no row and no message, and the rows come back [viewer=trellis plot, property=yColumnNames, value=RACE, drawn=rows shown]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "When user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["yColumnNames","RACE"]]), [["yColumnNames","RACE"]]);
    await session.step(31, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(32, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(33, "And the \"rows shown\" reading of trellis plot viewer should be 0", () => readingIs(page, "rows shown", el("trellis plot viewer"), 0));
    await session.step(34, "And the \"rows shown\" reading of trellis plot viewer should be 0", () => readingIs(page, "rows shown", el("trellis plot viewer"), 0));
    await session.step(35, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
    await session.step(36, "When user resets the filter", () => resetFilter(page));
    await session.step(37, "Then the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await session.step(38, "And the \"rows shown\" reading of trellis plot viewer should be at least 1", () => readingAtLeast(page, "rows shown", el("trellis plot viewer"), 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("The colorColumnName column removed from the table leaves the settings of the scatter plot [viewer=scatter plot, property=colorColumnName, column=RACE]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(58, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["colorColumnName","RACE"]]), [["colorColumnName","RACE"]]);
    await session.step(60, "And user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(61, "Then no errors should have been logged", () => noErrors(page));
    await session.step(62, "And \"colorColumnName\" property of scatter plot viewer should not be \"RACE\"", () => propertyShouldNotBe(page, "colorColumnName", el("scatter plot viewer"), "RACE"));
    await session.step(63, "And scatter plot viewer should report no error", () => reportsNoError(page, el("scatter plot viewer")));
  });
  test("The sizeColumnName column removed from the table leaves the settings of the scatter plot [viewer=scatter plot, property=sizeColumnName, column=HEIGHT]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(58, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["sizeColumnName","HEIGHT"]]), [["sizeColumnName","HEIGHT"]]);
    await session.step(60, "And user removes \"HEIGHT\" column", () => removeColumn(page, "HEIGHT"));
    await session.step(61, "Then no errors should have been logged", () => noErrors(page));
    await session.step(62, "And \"sizeColumnName\" property of scatter plot viewer should not be \"HEIGHT\"", () => propertyShouldNotBe(page, "sizeColumnName", el("scatter plot viewer"), "HEIGHT"));
    await session.step(63, "And scatter plot viewer should report no error", () => reportsNoError(page, el("scatter plot viewer")));
  });
  test("The lanesColumnName column removed from the table leaves the settings of the tile viewer [viewer=tile viewer, property=lanesColumnName, column=RACE]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(58, "When user adds a tile viewer viewer with:", () => addViewerWith(page, "tile viewer", [["lanesColumnName","RACE"]]), [["lanesColumnName","RACE"]]);
    await session.step(60, "And user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(61, "Then no errors should have been logged", () => noErrors(page));
    await session.step(62, "And \"lanesColumnName\" property of tile viewer viewer should not be \"RACE\"", () => propertyShouldNotBe(page, "lanesColumnName", el("tile viewer viewer"), "RACE"));
    await session.step(63, "And tile viewer viewer should report no error", () => reportsNoError(page, el("tile viewer viewer")));
  });
  test("The colorColumnName column removed from the table leaves the settings of the tree map [viewer=tree map, property=colorColumnName, column=AGE]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(58, "When user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["colorColumnName","AGE"]]), [["colorColumnName","AGE"]]);
    await session.step(60, "And user removes \"AGE\" column", () => removeColumn(page, "AGE"));
    await session.step(61, "Then no errors should have been logged", () => noErrors(page));
    await session.step(62, "And \"colorColumnName\" property of tree map viewer should not be \"AGE\"", () => propertyShouldNotBe(page, "colorColumnName", el("tree map viewer"), "AGE"));
    await session.step(63, "And tree map viewer should report no error", () => reportsNoError(page, el("tree map viewer")));
  });
  test("The sizeColumnName column removed from the table leaves the settings of the tree map [viewer=tree map, property=sizeColumnName, column=WEIGHT]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(58, "When user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["sizeColumnName","WEIGHT"]]), [["sizeColumnName","WEIGHT"]]);
    await session.step(60, "And user removes \"WEIGHT\" column", () => removeColumn(page, "WEIGHT"));
    await session.step(61, "Then no errors should have been logged", () => noErrors(page));
    await session.step(62, "And \"sizeColumnName\" property of tree map viewer should not be \"WEIGHT\"", () => propertyShouldNotBe(page, "sizeColumnName", el("tree map viewer"), "WEIGHT"));
    await session.step(63, "And tree map viewer should report no error", () => reportsNoError(page, el("tree map viewer")));
  });
  test("The category1ColumnName column removed from the table leaves the settings of the box plot [viewer=box plot, property=category1ColumnName, column=SEX]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(58, "When user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["category1ColumnName","SEX"]]), [["category1ColumnName","SEX"]]);
    await session.step(60, "And user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(61, "Then no errors should have been logged", () => noErrors(page));
    await session.step(62, "And \"category1ColumnName\" property of box plot viewer should not be \"SEX\"", () => propertyShouldNotBe(page, "category1ColumnName", el("box plot viewer"), "SEX"));
    await session.step(63, "And box plot viewer should report no error", () => reportsNoError(page, el("box plot viewer")));
  });
  test("The AGE column renamed is followed by the xColumnName setting of the scatter plot [viewer=scatter plot, property=xColumnName, value=AGE, column=AGE, new name=AGE_YEARS, renamed=AGE_YEARS]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(75, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"]]), [["xColumnName","AGE"]]);
    await session.step(77, "And user renames \"AGE\" column to \"AGE_YEARS\"", () => renameColumn(page, "AGE", "AGE_YEARS"));
    await session.step(78, "Then \"xColumnName\" property of scatter plot viewer should be \"AGE_YEARS\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "AGE_YEARS"));
    await session.step(79, "And scatter plot viewer should report no error", () => reportsNoError(page, el("scatter plot viewer")));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
  test("The RACE column renamed is followed by the colorColumnName setting of the scatter plot [viewer=scatter plot, property=colorColumnName, value=RACE, column=RACE, new name=ORIGIN, renamed=ORIGIN]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(75, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["colorColumnName","RACE"]]), [["colorColumnName","RACE"]]);
    await session.step(77, "And user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(78, "Then \"colorColumnName\" property of scatter plot viewer should be \"ORIGIN\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "ORIGIN"));
    await session.step(79, "And scatter plot viewer should report no error", () => reportsNoError(page, el("scatter plot viewer")));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
  test("The RACE column renamed is followed by the categoryColumnName setting of the pie chart [viewer=pie chart, property=categoryColumnName, value=RACE, column=RACE, new name=ORIGIN, renamed=ORIGIN]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(75, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","RACE"]]), [["categoryColumnName","RACE"]]);
    await session.step(77, "And user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(78, "Then \"categoryColumnName\" property of pie chart viewer should be \"ORIGIN\"", () => propertyShouldBe(page, "categoryColumnName", el("pie chart viewer"), "ORIGIN"));
    await session.step(79, "And pie chart viewer should report no error", () => reportsNoError(page, el("pie chart viewer")));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
  test("The HEIGHT column renamed is followed by the columnNames setting of the pc plot [viewer=pc plot, property=columnNames, value=AGE, HEIGHT, WEIGHT, column=HEIGHT, new name=TALL, renamed=AGE, TALL, WEIGHT]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(75, "When user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["columnNames","AGE, HEIGHT, WEIGHT"]]), [["columnNames","AGE, HEIGHT, WEIGHT"]]);
    await session.step(77, "And user renames \"HEIGHT\" column to \"TALL\"", () => renameColumn(page, "HEIGHT", "TALL"));
    await session.step(78, "Then \"columnNames\" property of pc plot viewer should be \"AGE, TALL, WEIGHT\"", () => propertyShouldBe(page, "columnNames", el("pc plot viewer"), "AGE, TALL, WEIGHT"));
    await session.step(79, "And pc plot viewer should report no error", () => reportsNoError(page, el("pc plot viewer")));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
  test("The RACE column renamed is followed by the lanesColumnName setting of the tile viewer [viewer=tile viewer, property=lanesColumnName, value=RACE, column=RACE, new name=ORIGIN, renamed=ORIGIN]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(75, "When user adds a tile viewer viewer with:", () => addViewerWith(page, "tile viewer", [["lanesColumnName","RACE"]]), [["lanesColumnName","RACE"]]);
    await session.step(77, "And user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(78, "Then \"lanesColumnName\" property of tile viewer viewer should be \"ORIGIN\"", () => propertyShouldBe(page, "lanesColumnName", el("tile viewer viewer"), "ORIGIN"));
    await session.step(79, "And tile viewer viewer should report no error", () => reportsNoError(page, el("tile viewer viewer")));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
  test("The SEX column renamed is followed by the splitByColumnNames setting of the tree map [viewer=tree map, property=splitByColumnNames, value=RACE, SEX, column=SEX, new name=GENDER, renamed=RACE, GENDER]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(75, "When user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE, SEX"]]), [["splitByColumnNames","RACE, SEX"]]);
    await session.step(77, "And user renames \"SEX\" column to \"GENDER\"", () => renameColumn(page, "SEX", "GENDER"));
    await session.step(78, "Then \"splitByColumnNames\" property of tree map viewer should be \"RACE, GENDER\"", () => propertyShouldBe(page, "splitByColumnNames", el("tree map viewer"), "RACE, GENDER"));
    await session.step(79, "And tree map viewer should report no error", () => reportsNoError(page, el("tree map viewer")));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
  test("The RACE column renamed is followed by the yColumnNames setting of the trellis plot [viewer=trellis plot, property=yColumnNames, value=RACE, column=RACE, new name=ORIGIN, renamed=ORIGIN]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(75, "When user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["yColumnNames","RACE"]]), [["yColumnNames","RACE"]]);
    await session.step(77, "And user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(78, "Then \"yColumnNames\" property of trellis plot viewer should be \"ORIGIN\"", () => propertyShouldBe(page, "yColumnNames", el("trellis plot viewer"), "ORIGIN"));
    await session.step(79, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
  test("The AGE column renamed is followed by the filter setting of the scatter plot [viewer=scatter plot, property=filter, value=${AGE} > 44, column=AGE, new name=YEARS, renamed=${YEARS} > 44]", {tag: ["@viewers", "@edge-cases", "@known-failure:78", "@realizes:core#1213"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(75, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["filter","${AGE} > 44"]]), [["filter","${AGE} > 44"]]);
      await session.step(77, "And user renames \"AGE\" column to \"YEARS\"", () => renameColumn(page, "AGE", "YEARS"));
      await session.step(78, "Then \"filter\" property of scatter plot viewer should be \"${YEARS} > 44\"", () => propertyShouldBe(page, "filter", el("scatter plot viewer"), "${YEARS} > 44"));
      await session.step(79, "And scatter plot viewer should report no error", () => reportsNoError(page, el("scatter plot viewer")));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    }, 78);
  });
  test("A layout naming a RACE column the table no longer has is applied to the scatter plot, which draws without it [viewer=scatter plot, property=colorColumnName, value=RACE, column=RACE, left=RACE]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(98, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["colorColumnName","RACE"]]), [["colorColumnName","RACE"]]);
    await session.step(100, "And user saves the layout of the current table view", () => saveLayout(page));
    await session.step(101, "And user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(102, "And user loads the saved layout", () => loadLayout(page));
    await session.step(103, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(104, "And \"colorColumnName\" property of scatter plot viewer should be \"RACE\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "RACE"));
    await session.step(105, "And scatter plot viewer should report no error", () => reportsNoError(page, el("scatter plot viewer")));
    await session.step(106, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout naming a HEIGHT column the table no longer has is applied to the pc plot, which draws without it [viewer=pc plot, property=columnNames, value=AGE, HEIGHT, WEIGHT, column=HEIGHT, left=AGE, HEIGHT, WEIGHT]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(98, "When user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["columnNames","AGE, HEIGHT, WEIGHT"]]), [["columnNames","AGE, HEIGHT, WEIGHT"]]);
    await session.step(100, "And user saves the layout of the current table view", () => saveLayout(page));
    await session.step(101, "And user removes \"HEIGHT\" column", () => removeColumn(page, "HEIGHT"));
    await session.step(102, "And user loads the saved layout", () => loadLayout(page));
    await session.step(103, "Then pc plot viewer should be visible", () => shouldBe(page, el("pc plot viewer"), "visible"));
    await session.step(104, "And \"columnNames\" property of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => propertyShouldBe(page, "columnNames", el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(105, "And pc plot viewer should report no error", () => reportsNoError(page, el("pc plot viewer")));
    await session.step(106, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout naming a SEX column the table no longer has is applied to the tree map, which draws without it [viewer=tree map, property=splitByColumnNames, value=RACE, SEX, column=SEX, left=RACE, SEX]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(98, "When user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE, SEX"]]), [["splitByColumnNames","RACE, SEX"]]);
    await session.step(100, "And user saves the layout of the current table view", () => saveLayout(page));
    await session.step(101, "And user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(102, "And user loads the saved layout", () => loadLayout(page));
    await session.step(103, "Then tree map viewer should be visible", () => shouldBe(page, el("tree map viewer"), "visible"));
    await session.step(104, "And \"splitByColumnNames\" property of tree map viewer should be \"RACE, SEX\"", () => propertyShouldBe(page, "splitByColumnNames", el("tree map viewer"), "RACE, SEX"));
    await session.step(105, "And tree map viewer should report no error", () => reportsNoError(page, el("tree map viewer")));
    await session.step(106, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout naming a HEIGHT column the table no longer has is applied to the forms, which draws without it [viewer=forms, property=fields, value=AGE, HEIGHT, SEX, column=HEIGHT, left=AGE, HEIGHT, SEX]", {tag: ["@viewers", "@edge-cases", "@known-failure:106", "@realizes:public#4452"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(98, "When user adds a forms viewer with:", () => addViewerWith(page, "forms", [["fields","AGE, HEIGHT, SEX"]]), [["fields","AGE, HEIGHT, SEX"]]);
      await session.step(100, "And user saves the layout of the current table view", () => saveLayout(page));
      await session.step(101, "And user removes \"HEIGHT\" column", () => removeColumn(page, "HEIGHT"));
      await session.step(102, "And user loads the saved layout", () => loadLayout(page));
      await session.step(103, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
      await session.step(104, "And \"fields\" property of forms viewer should be \"AGE, HEIGHT, SEX\"", () => propertyShouldBe(page, "fields", el("forms viewer"), "AGE, HEIGHT, SEX"));
      await session.step(105, "And forms viewer should report no error", () => reportsNoError(page, el("forms viewer")));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    }, 106);
  });
  test("Closing a bar chart gives back the rows its click filtered out [viewer=bar chart, property=splitColumnName, value=RACE, area=bar Asian]", {tag: ["@viewers", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(120, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"]]), [["splitColumnName","RACE"]]);
    await session.step(122, "And user sets \"onClick\" property of bar chart viewer to \"Filter\"", () => setProperty(page, "onClick", el("bar chart viewer"), "Filter"));
    await session.step(123, "And user clicks on the \"bar Asian\" area of bar chart viewer", () => clickArea(page, "bar Asian", el("bar chart viewer")));
    await session.step(124, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
    await session.step(125, "When user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
    await session.step(126, "Then the open tableview should have 0 bar chart viewers", () => viewerCount(page, 0, "bar chart"));
    await session.step(127, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(128, "And no errors should have been logged", () => noErrors(page));
  });
  test("Closing a calendar gives back the rows its click filtered out [viewer=calendar, property=onClick, value=Filter, area=month 1990-01]", {tag: ["@viewers", "@edge-cases", "@known-failure:127", "@realizes:core#233"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(120, "When user adds a calendar viewer with:", () => addViewerWith(page, "calendar", [["onClick","Filter"]]), [["onClick","Filter"]]);
      await session.step(122, "And user sets \"onClick\" property of calendar viewer to \"Filter\"", () => setProperty(page, "onClick", el("calendar viewer"), "Filter"));
      await session.step(123, "And user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(124, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(125, "When user clicks on close icon of calendar viewer", () => clickOn(page, el("close icon of calendar viewer")));
      await session.step(126, "Then the open tableview should have 0 calendar viewers", () => viewerCount(page, 0, "calendar"));
      await session.step(127, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(128, "And no errors should have been logged", () => noErrors(page));
    }, 127);
  });
  test("Removing the column of the category a bar chart click filtered by gives the rows back [viewer=bar chart, property=splitColumnName, value=RACE, area=bar Asian, column=RACE]", {tag: ["@viewers", "@edge-cases", "@known-failure:146", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(140, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"]]), [["splitColumnName","RACE"]]);
      await session.step(142, "And user sets \"onClick\" property of bar chart viewer to \"Filter\"", () => setProperty(page, "onClick", el("bar chart viewer"), "Filter"));
      await session.step(143, "And user clicks on the \"bar Asian\" area of bar chart viewer", () => clickArea(page, "bar Asian", el("bar chart viewer")));
      await session.step(144, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(145, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(146, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(147, "And no errors should have been logged", () => noErrors(page));
    }, 146);
  });
  test("Bins typed in the context panel outside its declared bounds is not taken by the histogram [viewer=histogram, category=X Axis, caption=Bins, typed=201, kept=20]", {tag: ["@viewers", "@edge-cases", "@known-failure:160", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(155, "When user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(156, "And user clicks on settings icon of histogram viewer", () => clickOn(page, el("settings icon of histogram viewer")));
      await session.step(157, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(158, "Given \"X Axis\" category in context panel is expanded", () => isExpanded(page, el("\"X Axis\" category in context panel")));
      await session.step(159, "When user enters \"201\" into \"Bins\" property in context panel", () => enterInto(page, "201", el("\"Bins\" property in context panel")));
      await session.step(160, "Then \"Bins\" property of histogram viewer should be \"20\"", () => propertyShouldBe(page, "Bins", el("histogram viewer"), "20"));
      await session.step(161, "And no errors should have been logged", () => noErrors(page));
    }, 160);
  });
  test("Max Radius typed in the context panel outside its declared bounds is not taken by the pie chart [viewer=pie chart, category=Misc, caption=Max Radius, typed=600, kept=150]", {tag: ["@viewers", "@edge-cases", "@known-failure:160", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(155, "When user adds a pie chart viewer", () => addViewer(page, "pie chart"));
      await session.step(156, "And user clicks on settings icon of pie chart viewer", () => clickOn(page, el("settings icon of pie chart viewer")));
      await session.step(157, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(158, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(159, "When user enters \"600\" into \"Max Radius\" property in context panel", () => enterInto(page, "600", el("\"Max Radius\" property in context panel")));
      await session.step(160, "Then \"Max Radius\" property of pie chart viewer should be \"150\"", () => propertyShouldBe(page, "Max Radius", el("pie chart viewer"), "150"));
      await session.step(161, "And no errors should have been logged", () => noErrors(page));
    }, 160);
  });
  test("Tiles Per Row typed in the context panel outside its declared bounds is not taken by the trellis plot [viewer=trellis plot, category=Misc, caption=Tiles Per Row, typed=11, kept=4]", {tag: ["@viewers", "@edge-cases", "@known-failure:160", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(155, "When user adds a trellis plot viewer", () => addViewer(page, "trellis plot"));
      await session.step(156, "And user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
      await session.step(157, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(158, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(159, "When user enters \"11\" into \"Tiles Per Row\" property in context panel", () => enterInto(page, "11", el("\"Tiles Per Row\" property in context panel")));
      await session.step(160, "Then \"Tiles Per Row\" property of trellis plot viewer should be \"4\"", () => propertyShouldBe(page, "Tiles Per Row", el("trellis plot viewer"), "4"));
      await session.step(161, "And no errors should have been logged", () => noErrors(page));
    }, 160);
  });
  test("bins set through the API outside its declared bounds is not taken by the histogram [viewer=histogram, property=bins, set=201, kept=20]", {tag: ["@viewers", "@edge-cases", "@known-failure:173", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(171, "When user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(172, "And user sets \"bins\" property of histogram viewer to \"201\"", () => setProperty(page, "bins", el("histogram viewer"), "201"));
      await session.step(173, "Then \"bins\" property of histogram viewer should be \"20\"", () => propertyShouldBe(page, "bins", el("histogram viewer"), "20"));
      await session.step(174, "And no errors should have been logged", () => noErrors(page));
    }, 173);
  });
  test("bins set through the API outside its declared bounds is not taken by the histogram [viewer=histogram, property=bins, set=0, kept=20]", {tag: ["@viewers", "@edge-cases", "@known-failure:173", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(171, "When user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(172, "And user sets \"bins\" property of histogram viewer to \"0\"", () => setProperty(page, "bins", el("histogram viewer"), "0"));
      await session.step(173, "Then \"bins\" property of histogram viewer should be \"20\"", () => propertyShouldBe(page, "bins", el("histogram viewer"), "20"));
      await session.step(174, "And no errors should have been logged", () => noErrors(page));
    }, 173);
  });
  test("maxRadius set through the API outside its declared bounds is not taken by the pie chart [viewer=pie chart, property=maxRadius, set=600, kept=150]", {tag: ["@viewers", "@edge-cases", "@known-failure:173", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(171, "When user adds a pie chart viewer", () => addViewer(page, "pie chart"));
      await session.step(172, "And user sets \"maxRadius\" property of pie chart viewer to \"600\"", () => setProperty(page, "maxRadius", el("pie chart viewer"), "600"));
      await session.step(173, "Then \"maxRadius\" property of pie chart viewer should be \"150\"", () => propertyShouldBe(page, "maxRadius", el("pie chart viewer"), "150"));
      await session.step(174, "And no errors should have been logged", () => noErrors(page));
    }, 173);
  });
  test("tilesPerRow set through the API outside its declared bounds is not taken by the trellis plot [viewer=trellis plot, property=tilesPerRow, set=0, kept=4]", {tag: ["@viewers", "@edge-cases", "@known-failure:173", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(171, "When user adds a trellis plot viewer", () => addViewer(page, "trellis plot"));
      await session.step(172, "And user sets \"tilesPerRow\" property of trellis plot viewer to \"0\"", () => setProperty(page, "tilesPerRow", el("trellis plot viewer"), "0"));
      await session.step(173, "Then \"tilesPerRow\" property of trellis plot viewer should be \"4\"", () => propertyShouldBe(page, "tilesPerRow", el("trellis plot viewer"), "4"));
      await session.step(174, "And no errors should have been logged", () => noErrors(page));
    }, 173);
  });
  test("An empty value is the (empty) category of the tile viewer [viewer=tile viewer, property=lanesColumnName, shown=tiles, reading=lane names]", {tag: ["@viewers", "@edge-cases", "@known-failure:195", "@realizes:core#1217"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(185, "Given user opens a table \"with-empty\" with:", () => openTableOf(page, "with-empty", [["site","arm"],["A","x"],["A",""],["B","y"]]), [["site","arm"],["A","x"],["A",""],["B","y"]]);
      await session.step(190, "When user adds a tile viewer viewer with:", () => addViewerWith(page, "tile viewer", [["lanesColumnName","arm"]]), [["lanesColumnName","arm"]]);
      await session.step(192, "Then no errors should have been logged", () => noErrors(page));
      await session.step(193, "And tile viewer viewer should report no error", () => reportsNoError(page, el("tile viewer viewer")));
      await session.step(194, "And the \"tiles\" reading of tile viewer viewer should be 3", () => readingIs(page, "tiles", el("tile viewer viewer"), 3));
      await session.step(195, "And the \"lane names\" reading of tile viewer viewer should contain \"(empty)\"", () => readingContains(page, "lane names", el("tile viewer viewer"), "(empty)"));
    }, 195);
  });
});
