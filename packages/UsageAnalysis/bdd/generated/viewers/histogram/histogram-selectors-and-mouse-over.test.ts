/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/histogram/histogram-selectors-and-mouse-over.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.histogram]
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
import {hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {mouseOverRowIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaRepainted, hasArea, hasNoArea, hoverArea, noErrors, notRepainted, pointerAway, readingIs, resizeTo, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Histogram split selector and mouse-over row", () => {
  const session = feature(test, "features/viewers/histogram/histogram-selectors-and-mouse-over.feature", import.meta.url);
  test("Histogram split selector and mouse-over row", {tag: ["@journey", "@viewers", "@realizes:viewers.histogram"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"]]), [["Value","AGE"]]);
    await session.step(14, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(15, "Then the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
    await run.scenario("Show Split Selector takes the split selector off the plot and brings it back", async () => {
      await session.step(18, "When user hovers over histogram viewer", () => hoverOver(page, el("histogram viewer")));
      await session.step(19, "Then Split column input in histogram viewer should be visible", () => shouldBe(page, el("Split column input in histogram viewer"), "visible"));
      await session.step(20, "When user sets \"Show Split Selector\" property of histogram viewer to \"false\"", () => setProperty(page, "Show Split Selector", el("histogram viewer"), "false"));
      await session.step(21, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(22, "And user hovers over histogram viewer", () => hoverOver(page, el("histogram viewer")));
      await session.step(23, "Then Split column input in histogram viewer should be hidden", () => shouldBe(page, el("Split column input in histogram viewer"), "hidden"));
      await session.step(24, "And Value column input in histogram viewer should be visible", () => shouldBe(page, el("Value column input in histogram viewer"), "visible"));
      await session.step(25, "When user sets \"Show Split Selector\" property of histogram viewer to \"true\"", () => setProperty(page, "Show Split Selector", el("histogram viewer"), "true"));
      await session.step(26, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(27, "And user hovers over histogram viewer", () => hoverOver(page, el("histogram viewer")));
      await session.step(28, "Then Split column input in histogram viewer should be visible", () => shouldBe(page, el("Split column input in histogram viewer"), "visible"));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Mouse Over Row decides whether a grid hover puts the marker on the baseline", async () => {
      await session.step(32, "Then histogram viewer should not have a \"mouse over row marker\" area", () => hasNoArea(page, el("histogram viewer"), "mouse over row marker"));
      await session.step(33, "When user takes a snapshot of histogram viewer", () => takeSnapshot(page, el("histogram viewer")));
      await session.step(34, "And user hovers over the \"cell 8 of AGE\" area of grid", () => hoverArea(page, "cell 8 of AGE", el("grid")));
      await session.step(35, "Then histogram viewer should have a \"mouse over row marker\" area", () => hasArea(page, el("histogram viewer"), "mouse over row marker"));
      await session.step(36, "And the \"mouse over row marker\" area of histogram viewer should have repainted", () => areaRepainted(page, "mouse over row marker", el("histogram viewer")));
      await session.step(37, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(38, "Then histogram viewer should not have a \"mouse over row marker\" area", () => hasNoArea(page, el("histogram viewer"), "mouse over row marker"));
      await session.step(39, "When user sets \"Show Mouse Over Row\" property of histogram viewer to \"false\"", () => setProperty(page, "Show Mouse Over Row", el("histogram viewer"), "false"));
      await session.step(40, "And user takes a snapshot of histogram viewer", () => takeSnapshot(page, el("histogram viewer")));
      await session.step(41, "And user hovers over the \"cell 8 of AGE\" area of grid", () => hoverArea(page, "cell 8 of AGE", el("grid")));
      await session.step(42, "Then the mouse-over row of the table should be 8", () => mouseOverRowIs(page, 8));
      await session.step(43, "And histogram viewer should not have a \"mouse over row marker\" area", () => hasNoArea(page, el("histogram viewer"), "mouse over row marker"));
      await session.step(44, "And histogram viewer should not have repainted", () => notRepainted(page, el("histogram viewer")));
      await session.step(45, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(46, "And user sets \"Show Mouse Over Row\" property of histogram viewer to \"true\"", () => setProperty(page, "Show Mouse Over Row", el("histogram viewer"), "true"));
      await session.step(47, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
