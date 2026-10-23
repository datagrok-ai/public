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
import {addViewerWith, areaAtLeastTall, areaRepainted, areaShorter, areaTaller, dragAreaBy, hasArea, hasNoArea, hoverArea, noErrors, notRepainted, pointerAway, propertyShouldBe, readingIs, resizeTo, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Histogram split selector, bins slider, x axis height and mouse-over row", () => {
  const session = feature(test, "features/viewers/histogram/histogram-selectors-and-mouse-over.feature", import.meta.url);
  test("Histogram split selector, bins slider, x axis height and mouse-over row", {tag: ["@journey", "@viewers", "@realizes:viewers.histogram"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"]]), [["Value","AGE"]]);
    await session.step(16, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(17, "Then the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
    await run.scenario("Show Split Selector takes the split selector off the plot and brings it back", async () => {
      await session.step(20, "When user hovers over histogram viewer", () => hoverOver(page, el("histogram viewer")));
      await session.step(21, "Then Split column input in histogram viewer should be visible", () => shouldBe(page, el("Split column input in histogram viewer"), "visible"));
      await session.step(22, "When user sets \"Show Split Selector\" property of histogram viewer to \"false\"", () => setProperty(page, "Show Split Selector", el("histogram viewer"), "false"));
      await session.step(23, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(24, "And user hovers over histogram viewer", () => hoverOver(page, el("histogram viewer")));
      await session.step(25, "Then Split column input in histogram viewer should be hidden", () => shouldBe(page, el("Split column input in histogram viewer"), "hidden"));
      await session.step(26, "And Value column input in histogram viewer should be visible", () => shouldBe(page, el("Value column input in histogram viewer"), "visible"));
      await session.step(27, "When user sets \"Show Split Selector\" property of histogram viewer to \"true\"", () => setProperty(page, "Show Split Selector", el("histogram viewer"), "true"));
      await session.step(28, "And user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(29, "And user hovers over histogram viewer", () => hoverOver(page, el("histogram viewer")));
      await session.step(30, "Then Split column input in histogram viewer should be visible", () => shouldBe(page, el("Split column input in histogram viewer"), "visible"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The bins slider and the Bins property agree on the bin count (github-2296)", async () => {
      await session.step(34, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(35, "And user hovers over histogram viewer", () => hoverOver(page, el("histogram viewer")));
      await session.step(36, "And user drags the \"bins slider\" area of histogram viewer by 40 pixels to the right", () => dragAreaBy(page, "bins slider", el("histogram viewer"), 40, "right"));
      await session.step(37, "Then \"Bins\" property of histogram viewer should be \"169\"", () => propertyShouldBe(page, "Bins", el("histogram viewer"), "169"));
      await session.step(38, "And the \"bins shown\" reading of histogram viewer should be 169", () => readingIs(page, "bins shown", el("histogram viewer"), 169));
      await session.step(39, "When user sets \"Bins\" property of histogram viewer to \"20\"", () => setProperty(page, "Bins", el("histogram viewer"), "20"));
      await session.step(40, "Then the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Auto Layout off, X Axis Height makes the x axis box taller and gives the height back", async () => {
      await session.step(44, "When user sets \"Auto Layout\" property of histogram viewer to \"false\"", () => setProperty(page, "Auto Layout", el("histogram viewer"), "false"));
      await session.step(45, "And user sets \"X Axis Height\" property of histogram viewer to \"60\"", () => setProperty(page, "X Axis Height", el("histogram viewer"), "60"));
      await session.step(46, "Then the \"x axis\" area of histogram viewer should be taller than before", () => areaTaller(page, "x axis", el("histogram viewer")));
      await session.step(47, "And the \"x axis\" area of histogram viewer should be at least 60 pixels tall", () => areaAtLeastTall(page, "x axis", el("histogram viewer"), 60));
      await session.step(48, "When user sets \"X Axis Height\" property of histogram viewer to \"30\"", () => setProperty(page, "X Axis Height", el("histogram viewer"), "30"));
      await session.step(49, "Then the \"x axis\" area of histogram viewer should be shorter than before", () => areaShorter(page, "x axis", el("histogram viewer")));
      await session.step(50, "When user sets \"Auto Layout\" property of histogram viewer to \"true\"", () => setProperty(page, "Auto Layout", el("histogram viewer"), "true"));
      await session.step(51, "Then histogram viewer should have an \"x axis\" area", () => hasArea(page, el("histogram viewer"), "x axis"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Mouse Over Row decides whether a grid hover puts the marker on the baseline", async () => {
      await session.step(55, "Given histogram viewer should not have a \"mouse over row marker\" area", () => hasNoArea(page, el("histogram viewer"), "mouse over row marker"));
      await session.step(56, "When user takes a snapshot of histogram viewer", () => takeSnapshot(page, el("histogram viewer")));
      await session.step(57, "And user hovers over the \"cell 8 of AGE\" area of grid", () => hoverArea(page, "cell 8 of AGE", el("grid")));
      await session.step(58, "Then histogram viewer should have a \"mouse over row marker\" area", () => hasArea(page, el("histogram viewer"), "mouse over row marker"));
      await session.step(59, "And the \"mouse over row marker\" area of histogram viewer should have repainted", () => areaRepainted(page, "mouse over row marker", el("histogram viewer")));
      await session.step(60, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(61, "Then histogram viewer should not have a \"mouse over row marker\" area", () => hasNoArea(page, el("histogram viewer"), "mouse over row marker"));
      await session.step(62, "When user sets \"Show Mouse Over Row\" property of histogram viewer to \"false\"", () => setProperty(page, "Show Mouse Over Row", el("histogram viewer"), "false"));
      await session.step(63, "And user takes a snapshot of histogram viewer", () => takeSnapshot(page, el("histogram viewer")));
      await session.step(64, "And user hovers over the \"cell 8 of AGE\" area of grid", () => hoverArea(page, "cell 8 of AGE", el("grid")));
      await session.step(65, "Then the mouse-over row of the table should be 8", () => mouseOverRowIs(page, 8));
      await session.step(66, "And histogram viewer should not have a \"mouse over row marker\" area", () => hasNoArea(page, el("histogram viewer"), "mouse over row marker"));
      await session.step(67, "And histogram viewer should not have repainted", () => notRepainted(page, el("histogram viewer")));
      await session.step(68, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(69, "And user sets \"Show Mouse Over Row\" property of histogram viewer to \"true\"", () => setProperty(page, "Show Mouse Over Row", el("histogram viewer"), "true"));
      await session.step(70, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
