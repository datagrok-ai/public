/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-legend-no-value.feature
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
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, addCategoricalFilter, filterPasses, filterPassesAll, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickLegendItem, legendLists, noErrors, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot — the \"(no value)\" legend item as a click target", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-legend-no-value.feature", import.meta.url);
  test("Scatter plot — the \"(no value)\" legend item as a click target", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(19, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await run.scenario("A click on \"(no value)\" draws only the rows with no colour value, and a panel card narrows them further", async () => {
      await session.step(22, "When user adds a calculated column \"DIS_POP but Other\" with formula \"if(${RACE}=='Other', null, ${DIS_POP})\"", () => addCalculated(page, "DIS_POP but Other", "if(${RACE}=='Other', null, ${DIS_POP})"));
      await session.step(23, "And user sets \"Color\" property of scatter plot viewer to \"DIS_POP but Other\"", () => setProperty(page, "Color", el("scatter plot viewer"), "DIS_POP but Other"));
      await session.step(24, "Then the legend of scatter plot viewer should list 6 items", () => legendLists(page, el("scatter plot viewer"), 6));
      await session.step(25, "And \"(no value)\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"(no value)\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(26, "When user clicks on \"(no value)\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "(no value)", el("scatter plot viewer")));
      await session.step(27, "Then scatter plot viewer should show 59 rows", () => showsRows(page, el("scatter plot viewer"), 59));
      await session.step(28, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(29, "And \"(no value)\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"(no value)\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(30, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(31, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(32, "And scatter plot viewer should show 12 rows", () => showsRows(page, el("scatter plot viewer"), 12));
      await session.step(33, "When user clicks on \"(no value)\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "(no value)", el("scatter plot viewer")));
      await session.step(34, "Then scatter plot viewer should show 373 rows", () => showsRows(page, el("scatter plot viewer"), 373));
      await session.step(35, "And 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(36, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(37, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(38, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(39, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(40, "When user sets \"Color\" property of scatter plot viewer to \"\"", () => setProperty(page, "Color", el("scatter plot viewer"), ""));
      await session.step(41, "And user removes \"DIS_POP but Other\" column", () => removeColumn(page, "DIS_POP but Other"));
      await session.step(42, "Then legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
