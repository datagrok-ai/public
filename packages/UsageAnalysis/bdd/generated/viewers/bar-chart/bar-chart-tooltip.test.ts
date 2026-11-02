/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-tooltip.feature
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
import {clickOn, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hoverArea, noErrors, pointerAway, readingIs, setProperties, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart bar tooltip", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-tooltip.feature", import.meta.url);
  test("Bar chart bar tooltip", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Value Aggr Type","count"]]), [["Split","RACE"],["Value Aggr Type","count"]]);
    await run.scenario("The bar tooltip reports the category, its count and its rows", async () => {
      await session.step(19, "When user hovers over the \"bar Caucasian\" area of bar chart viewer", () => hoverArea(page, "bar Caucasian", el("bar chart viewer")));
      await session.step(20, "Then tooltip should contain text \"Caucasian\"", () => shouldContainText(page, el("tooltip"), "Caucasian"));
      await session.step(21, "And tooltip should contain text \"count: 896\"", () => shouldContainText(page, el("tooltip"), "count: 896"));
      await session.step(22, "And tooltip should contain text \"896 rows\"", () => shouldContainText(page, el("tooltip"), "896 rows"));
      await session.step(23, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(24, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Switched to the average of WEIGHT, the bar tooltip reports the average instead of the count", async () => {
      await session.step(27, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Value","WEIGHT"],["Value Aggr Type","avg"]]), [["Value","WEIGHT"],["Value Aggr Type","avg"]]);
      await session.step(30, "And user hovers over the \"bar Caucasian\" area of bar chart viewer", () => hoverArea(page, "bar Caucasian", el("bar chart viewer")));
      await session.step(31, "Then tooltip should contain text \"avg: 79.49\"", () => shouldContainText(page, el("tooltip"), "avg: 79.49"));
      await session.step(32, "And tooltip should not contain text \"count: 896\"", () => shouldNotContainText(page, el("tooltip"), "count: 896"));
      await session.step(33, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(34, "And user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Value Aggr Type","count"]]), [["Value Aggr Type","count"]]);
      await session.step(36, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the bar chart gives back the rows its Filter click filtered out", async () => {
      await session.step(39, "When user sets \"onClick\" property of bar chart viewer to \"Filter\"", () => setProperty(page, "onClick", el("bar chart viewer"), "Filter"));
      await session.step(40, "And user clicks on the \"bar Asian\" area of bar chart viewer", () => clickArea(page, "bar Asian", el("bar chart viewer")));
      await session.step(41, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(42, "When user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(43, "Then the open tableview should have 0 bar chart viewers", () => viewerCount(page, 0, "bar chart"));
      await session.step(44, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A category with markup characters is shown in the bar tooltip as typed (core#974)", async () => {
      await session.step(49, "Given user opens a table \"markup-categories\" with:", () => openTableOf(page, "markup-categories", [["Category","Value"],["a<b","1"],["a<b","2"],["x &amp; y","3"]]), [["Category","Value"],["a<b","1"],["a<b","2"],["x &amp; y","3"]]);
      await session.step(54, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Category"]]), [["Split","Category"]]);
      await session.step(56, "Then the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
      await session.step(57, "When user hovers over the \"bar a<b\" area of bar chart viewer", () => hoverArea(page, "bar a<b", el("bar chart viewer")));
      await session.step(58, "Then tooltip should contain text \"a<b\"", () => shouldContainText(page, el("tooltip"), "a<b"));
      await session.step(59, "And tooltip should contain text \"2 rows\"", () => shouldContainText(page, el("tooltip"), "2 rows"));
      await session.step(60, "When user hovers over the \"bar x &amp; y\" area of bar chart viewer", () => hoverArea(page, "bar x &amp; y", el("bar chart viewer")));
      await session.step(61, "Then tooltip should contain text \"x &amp; y\"", () => shouldContainText(page, el("tooltip"), "x &amp; y"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
