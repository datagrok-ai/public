/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/legend/legend-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.legend]
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
import {shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, filterBetween, filterPasses, removeColumn, renameColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickLegendItem, legendLists, noErrors, propertyShouldBe, propertyShouldNotBe, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Legend edge cases: the color column removed under a legend selection, every row filtered out, a constant column", () => {
  const session = feature(test, "features/viewers/legend/legend-edge-cases.feature", import.meta.url);
  test("Removing the color column under a legend selection gives every row back", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.legend"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Position","Right"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Position","Right"]]);
    await session.step(20, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(21, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(24, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
    await session.step(25, "Then scatter plot viewer should show 25 rows", () => showsRows(page, el("scatter plot viewer"), 25));
    await session.step(26, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(27, "Then \"Color\" property of scatter plot viewer should not be \"RACE\"", () => propertyShouldNotBe(page, "Color", el("scatter plot viewer"), "RACE"));
    await session.step(28, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
  });
  test("A legend selection outlives a filter that removes every row", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.legend"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Position","Right"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Position","Right"]]);
    await session.step(20, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(21, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(32, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
    await session.step(33, "Then scatter plot viewer should show 25 rows", () => showsRows(page, el("scatter plot viewer"), 25));
    await session.step(34, "When user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(35, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(36, "And scatter plot viewer should show 0 rows", () => showsRows(page, el("scatter plot viewer"), 0));
    await session.step(37, "When user resets the filter", () => resetFilter(page));
    await session.step(38, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(39, "And \"Black\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
    await session.step(40, "And scatter plot viewer should show 25 rows", () => showsRows(page, el("scatter plot viewer"), 25));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
  test("A color column holding one value gives a legend of one item", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.legend"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Position","Right"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Position","Right"]]);
    await session.step(20, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(21, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(44, "When user adds a calculated column \"ONE\" with formula \"'same'\"", () => addCalculated(page, "ONE", "'same'"));
    await session.step(45, "And user sets \"Color\" property of scatter plot viewer to \"ONE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "ONE"));
    await session.step(46, "Then the legend of scatter plot viewer should list 1 item", () => legendLists(page, el("scatter plot viewer"), 1));
    await session.step(47, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
  test("Renaming the color column releases the legend selection on the legend and on the plot", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.legend"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Position","Right"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Position","Right"]]);
    await session.step(20, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(21, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(51, "When user clicks on \"Black\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Black", el("scatter plot viewer")));
    await session.step(52, "Then scatter plot viewer should show 25 rows", () => showsRows(page, el("scatter plot viewer"), 25));
    await session.step(53, "When user renames \"RACE\" column to \"ETHNICITY\"", () => renameColumn(page, "RACE", "ETHNICITY"));
    await session.step(54, "Then \"Color\" property of scatter plot viewer should be \"ETHNICITY\"", () => propertyShouldBe(page, "Color", el("scatter plot viewer"), "ETHNICITY"));
    await session.step(55, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(56, "And \"Black\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
    await session.step(57, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(58, "And no errors should have been logged", () => noErrors(page));
  });
});
