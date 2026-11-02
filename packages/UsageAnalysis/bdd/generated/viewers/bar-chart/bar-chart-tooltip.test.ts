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
import {shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeCurrentView, openDataset, openTableOf, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hoverArea, noErrors, pointerAway, readingIs, setProperties} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart bar tooltip", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-tooltip.feature", import.meta.url);
  test("Bar chart bar tooltip", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Value Aggr Type","count"]]), [["Split","RACE"],["Value Aggr Type","count"]]);
    await run.scenario("The bar tooltip reports the category, its count and its rows", async () => {
      await session.step(20, "When user hovers over the \"bar Caucasian\" area of bar chart viewer", () => hoverArea(page, "bar Caucasian", el("bar chart viewer")));
      await session.step(21, "Then tooltip should contain text \"Caucasian\"", () => shouldContainText(page, el("tooltip"), "Caucasian"));
      await session.step(22, "And tooltip should contain text \"count: 896\"", () => shouldContainText(page, el("tooltip"), "count: 896"));
      await session.step(23, "And tooltip should contain text \"896 rows\"", () => shouldContainText(page, el("tooltip"), "896 rows"));
      await session.step(24, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(25, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Switched to the average of WEIGHT, the bar tooltip reports the average instead of the count", async () => {
      await session.step(28, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Value","WEIGHT"],["Value Aggr Type","avg"]]), [["Value","WEIGHT"],["Value Aggr Type","avg"]]);
      await session.step(31, "And user hovers over the \"bar Caucasian\" area of bar chart viewer", () => hoverArea(page, "bar Caucasian", el("bar chart viewer")));
      await session.step(32, "Then tooltip should contain text \"avg: 79.49\"", () => shouldContainText(page, el("tooltip"), "avg: 79.49"));
      await session.step(33, "And tooltip should not contain text \"count: 896\"", () => shouldNotContainText(page, el("tooltip"), "count: 896"));
      await session.step(34, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(35, "And user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Value Aggr Type","count"]]), [["Value Aggr Type","count"]]);
      await session.step(37, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A category with markup characters is shown in the bar tooltip as typed (core#974)", async () => {
      await session.step(41, "Given user opens a table \"markup-categories\" with:", () => openTableOf(page, "markup-categories", [["Category","Value"],["a<b","1"],["a<b","2"],["x &amp; y","3"]]), [["Category","Value"],["a<b","1"],["a<b","2"],["x &amp; y","3"]]);
      await session.step(46, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Category"]]), [["Split","Category"]]);
      await session.step(48, "Then the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
      await session.step(49, "When user hovers over the \"bar a<b\" area of bar chart viewer", () => hoverArea(page, "bar a<b", el("bar chart viewer")));
      await session.step(50, "Then tooltip should contain text \"a<b\"", () => shouldContainText(page, el("tooltip"), "a<b"));
      await session.step(51, "And tooltip should contain text \"2 rows\"", () => shouldContainText(page, el("tooltip"), "2 rows"));
      await session.step(52, "When user hovers over the \"bar x &amp; y\" area of bar chart viewer", () => hoverArea(page, "bar x &amp; y", el("bar chart viewer")));
      await session.step(53, "Then tooltip should contain text \"x &amp; y\"", () => shouldContainText(page, el("tooltip"), "x &amp; y"));
      await session.step(54, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
      await session.step(55, "And user closes the current view", () => closeCurrentView(page));
      await session.step(56, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(57, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
