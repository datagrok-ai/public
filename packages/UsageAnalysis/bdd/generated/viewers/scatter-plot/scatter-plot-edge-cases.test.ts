/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
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
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterBetween, filterPasses, filterPassesAll, filterPassesFewer, removeColumn, renameColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, legendLists, lessInk, loadLayout, moreInk, noErrors, propertyShouldBe, propertyShouldNotBe, reportsNoError, saveLayout, setProperty, showsRows, viewerCount, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot edge cases: empty sizes, no rows, one row, columns removed or renamed, a layout with a missing column", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-edge-cases.feature", import.meta.url);
  test("Rows with an empty size are hidden until Show Markers With Empty Size is on", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(22, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(25, "When user sets \"Size\" property of scatter plot viewer to \"HEIGHT\"", () => setProperty(page, "Size", el("scatter plot viewer"), "HEIGHT"));
    await session.step(26, "And user sets \"Show Markers With Empty Size\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Markers With Empty Size", el("scatter plot viewer"), "true"));
    await session.step(27, "Then scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
    await session.step(28, "When user sets \"Show Markers With Empty Size\" property of scatter plot viewer to \"false\"", () => setProperty(page, "Show Markers With Empty Size", el("scatter plot viewer"), "false"));
    await session.step(29, "Then scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
  });
  test("A filter that removes every row leaves the plot empty without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(22, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(33, "When user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(34, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(35, "And scatter plot viewer should show 0 rows", () => showsRows(page, el("scatter plot viewer"), 0));
    await session.step(36, "And scatter plot viewer should report no error", () => reportsNoError(page, el("scatter plot viewer")));
    await session.step(37, "When user resets the filter", () => resetFilter(page));
    await session.step(38, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("A one-row table draws its one marker", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(22, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(42, "Given user opens a table \"one-row\" with:", () => openTableOf(page, "one-row", [["name","x","y"],["a","1","2"]]), [["name","x","y"],["a","1","2"]]);
    await session.step(45, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","x"],["Y","y"]]), [["X","x"],["Y","y"]]);
    await session.step(48, "Then scatter plot viewer should show 1 rows", () => showsRows(page, el("scatter plot viewer"), 1));
    await session.step(49, "And scatter plot viewer should have a \"marker of row 1\" area", () => hasArea(page, el("scatter plot viewer"), "marker of row 1"));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the X column drops it from the plot's settings without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(22, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(53, "When user removes \"AGE\" column", () => removeColumn(page, "AGE"));
    await session.step(54, "Then the table should not have a column \"AGE\"", () => hasNoColumn(page, "AGE"));
    await session.step(55, "And \"X\" property of scatter plot viewer should not be \"AGE\"", () => propertyShouldNotBe(page, "X", el("scatter plot viewer"), "AGE"));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the Color column drops it from the plot's settings without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(22, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(59, "When user sets \"Color\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE"));
    await session.step(60, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(61, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(62, "Then \"Color\" property of scatter plot viewer should not be \"RACE\"", () => propertyShouldNotBe(page, "Color", el("scatter plot viewer"), "RACE"));
    await session.step(63, "And scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
  test("Renaming the X and Color columns is followed by the plot's settings", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(22, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(67, "When user sets \"Color\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE"));
    await session.step(68, "And user renames \"AGE\" column to \"AGE_YEARS\"", () => renameColumn(page, "AGE", "AGE_YEARS"));
    await session.step(69, "And user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(70, "Then \"X\" property of scatter plot viewer should be \"AGE_YEARS\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "AGE_YEARS"));
    await session.step(71, "And \"Color\" property of scatter plot viewer should be \"ORIGIN\"", () => propertyShouldBe(page, "Color", el("scatter plot viewer"), "ORIGIN"));
    await session.step(72, "And scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(73, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(74, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout naming a column the table no longer has is applied without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(22, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(77, "When user sets \"Color\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE"));
    await session.step(78, "And user saves the layout of the current table view", () => saveLayout(page));
    await session.step(79, "And user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(80, "And user loads the saved layout", () => loadLayout(page));
    await session.step(81, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(82, "And scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(83, "And no errors should have been logged", () => noErrors(page));
  });
  test("Closing the plot gives back the rows its zoom filtered out", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(22, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(86, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
    await session.step(87, "And user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
    await session.step(88, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
    await session.step(89, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
    await session.step(90, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
    await session.step(91, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(92, "And no errors should have been logged", () => noErrors(page));
  });
});
