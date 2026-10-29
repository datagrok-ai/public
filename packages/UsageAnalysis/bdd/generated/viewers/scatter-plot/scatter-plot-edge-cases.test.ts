/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, core#1200]
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
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterPassesAll, filterPassesFewer, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, lessInk, moreInk, noErrors, propertyShouldBe, propertyShouldNotBe, setProperty, showsRows, viewerCount, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot edge cases: empty sizes, one row, the X column removed, a filter by zoom left on close", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-edge-cases.feature", import.meta.url);
  test("Rows with an empty size are hidden until Show Markers With Empty Size is on", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(21, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(24, "When user sets \"Size\" property of scatter plot viewer to \"HEIGHT\"", () => setProperty(page, "Size", el("scatter plot viewer"), "HEIGHT"));
    await session.step(25, "And user sets \"Show Markers With Empty Size\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Markers With Empty Size", el("scatter plot viewer"), "true"));
    await session.step(26, "Then scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
    await session.step(27, "When user sets \"Show Markers With Empty Size\" property of scatter plot viewer to \"false\"", () => setProperty(page, "Show Markers With Empty Size", el("scatter plot viewer"), "false"));
    await session.step(28, "Then scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
  });
  test("A one-row table draws its one marker", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(21, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(32, "Given user opens a table \"one-row\" with:", () => openTableOf(page, "one-row", [["name","x","y"],["a","1","2"]]), [["name","x","y"],["a","1","2"]]);
    await session.step(35, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","x"],["Y","y"]]), [["X","x"],["Y","y"]]);
    await session.step(38, "Then scatter plot viewer should show 1 rows", () => showsRows(page, el("scatter plot viewer"), 1));
    await session.step(39, "And scatter plot viewer should have a \"marker of row 1\" area", () => hasArea(page, el("scatter plot viewer"), "marker of row 1"));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the X column drops it from the plot's settings without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(21, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(43, "When user removes \"AGE\" column", () => removeColumn(page, "AGE"));
    await session.step(44, "Then the table should not have a column \"AGE\"", () => hasNoColumn(page, "AGE"));
    await session.step(45, "And \"X\" property of scatter plot viewer should not be \"AGE\"", () => propertyShouldNotBe(page, "X", el("scatter plot viewer"), "AGE"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
  test("Closing the plot gives back the rows its zoom filtered out", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(21, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(49, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
    await session.step(50, "And user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
    await session.step(51, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
    await session.step(52, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
    await session.step(53, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
    await session.step(54, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
  });
  test("Moving Average Window typed in the context panel as 1 is not taken [typed=1]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot", "@known-failure:63", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(21, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await knownFailure(async () => {
      await session.step(58, "When user sets \"Show Moving Average Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Line", el("scatter plot viewer"), "true"));
      await session.step(59, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(60, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(61, "Given \"Statistics\" category in context panel is expanded", () => isExpanded(page, el("\"Statistics\" category in context panel")));
      await session.step(62, "When user enters \"1\" into \"Moving Average Window\" property in context panel", () => enterInto(page, "1", el("\"Moving Average Window\" property in context panel")));
      await session.step(63, "Then \"Moving Average Window\" property of scatter plot viewer should be \"10\"", () => propertyShouldBe(page, "Moving Average Window", el("scatter plot viewer"), "10"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    }, 63);
  });
  test("Moving Average Window typed in the context panel as 101 is not taken [typed=101]", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot", "@known-failure:63", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(21, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await knownFailure(async () => {
      await session.step(58, "When user sets \"Show Moving Average Line\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Show Moving Average Line", el("scatter plot viewer"), "true"));
      await session.step(59, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(60, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(61, "Given \"Statistics\" category in context panel is expanded", () => isExpanded(page, el("\"Statistics\" category in context panel")));
      await session.step(62, "When user enters \"101\" into \"Moving Average Window\" property in context panel", () => enterInto(page, "101", el("\"Moving Average Window\" property in context panel")));
      await session.step(63, "Then \"Moving Average Window\" property of scatter plot viewer should be \"10\"", () => propertyShouldBe(page, "Moving Average Window", el("scatter plot viewer"), "10"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    }, 63);
  });
});
