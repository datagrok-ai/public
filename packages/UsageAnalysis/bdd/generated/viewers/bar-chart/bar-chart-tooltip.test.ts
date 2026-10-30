/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-tooltip.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart]
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hoverArea, noErrors, pointerAway, setProperties} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart bar tooltip", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-tooltip.feature", import.meta.url);
  test("Bar chart bar tooltip", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Value Aggr Type","count"]]), [["Split","RACE"],["Value Aggr Type","count"]]);
    await run.scenario("The bar tooltip reports the category, its count and its rows", async () => {
      await session.step(16, "When user hovers over the \"bar Caucasian\" area of bar chart viewer", () => hoverArea(page, "bar Caucasian", el("bar chart viewer")));
      await session.step(17, "Then tooltip should contain text \"Caucasian\"", () => shouldContainText(page, el("tooltip"), "Caucasian"));
      await session.step(18, "And tooltip should contain text \"count: 896\"", () => shouldContainText(page, el("tooltip"), "count: 896"));
      await session.step(19, "And tooltip should contain text \"896 rows\"", () => shouldContainText(page, el("tooltip"), "896 rows"));
      await session.step(20, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(21, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Switched to the average of WEIGHT, the bar tooltip reports the average instead of the count", async () => {
      await session.step(24, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Value","WEIGHT"],["Value Aggr Type","avg"]]), [["Value","WEIGHT"],["Value Aggr Type","avg"]]);
      await session.step(27, "And user hovers over the \"bar Caucasian\" area of bar chart viewer", () => hoverArea(page, "bar Caucasian", el("bar chart viewer")));
      await session.step(28, "Then tooltip should contain text \"avg: 79.49\"", () => shouldContainText(page, el("tooltip"), "avg: 79.49"));
      await session.step(29, "And tooltip should not contain text \"count: 896\"", () => shouldNotContainText(page, el("tooltip"), "count: 896"));
      await session.step(30, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(31, "And user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Value Aggr Type","count"]]), [["Value Aggr Type","count"]]);
      await session.step(33, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
