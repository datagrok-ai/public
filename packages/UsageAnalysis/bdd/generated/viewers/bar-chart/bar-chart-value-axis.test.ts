/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-value-axis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/nx.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, painted, pointerAway, propertyShouldBe, readingIs, repainted, repaintedBy, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart value axis range and scale", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-value-axis.feature", import.meta.url);
  test("Bar chart value axis range and scale", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(13, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Primary Series Name"],["Value","CAST Idea ID"],["Value Aggr Type","count"],["Show Clipped Bar Indicators","true"]]), [["Split","Primary Series Name"],["Value","CAST Idea ID"],["Value Aggr Type","count"],["Show Clipped Bar Indicators","true"]]);
    await session.step(18, "Then the \"bars\" reading of bar chart viewer should be 5", () => readingIs(page, "bars", el("bar chart viewer"), 5));
    await session.step(19, "And the \"clipped bars\" reading of bar chart viewer should be 0", () => readingIs(page, "clipped bars", el("bar chart viewer"), 0));
    await run.scenario("Min above the shortest bars clips them and the indicators show it", async () => {
      await session.step(22, "When user sets \"Min\" property of bar chart viewer to \"10\"", () => setProperty(page, "Min", el("bar chart viewer"), "10"));
      await session.step(23, "Then the \"clipped bars\" reading of bar chart viewer should be 3", () => readingIs(page, "clipped bars", el("bar chart viewer"), 3));
      await session.step(24, "When user sets \"Show Clipped Bar Indicators\" property of bar chart viewer to \"false\"", () => setProperty(page, "Show Clipped Bar Indicators", el("bar chart viewer"), "false"));
      await session.step(25, "Then bar chart viewer should have repainted by at least 150 pixels", () => repaintedBy(page, el("bar chart viewer"), 150));
      await session.step(26, "When user sets \"Show Clipped Bar Indicators\" property of bar chart viewer to \"true\"", () => setProperty(page, "Show Clipped Bar Indicators", el("bar chart viewer"), "true"));
      await session.step(27, "Then bar chart viewer should have repainted by at least 150 pixels", () => repaintedBy(page, el("bar chart viewer"), 150));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Max below the tallest bar clips it too", async () => {
      await session.step(31, "When user sets \"Max\" property of bar chart viewer to \"60\"", () => setProperty(page, "Max", el("bar chart viewer"), "60"));
      await session.step(32, "Then the \"clipped bars\" reading of bar chart viewer should be 4", () => readingIs(page, "clipped bars", el("bar chart viewer"), 4));
      await session.step(33, "When user sets \"Show Clipped Bar Indicators\" property of bar chart viewer to \"false\"", () => setProperty(page, "Show Clipped Bar Indicators", el("bar chart viewer"), "false"));
      await session.step(34, "Then bar chart viewer should have repainted by at least 150 pixels", () => repaintedBy(page, el("bar chart viewer"), 150));
      await session.step(35, "When user sets \"Show Clipped Bar Indicators\" property of bar chart viewer to \"true\"", () => setProperty(page, "Show Clipped Bar Indicators", el("bar chart viewer"), "true"));
      await session.step(36, "Then bar chart viewer should have repainted by at least 150 pixels", () => repaintedBy(page, el("bar chart viewer"), 150));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The value-axis slider is there under the pointer on the constrained range", async () => {
      await session.step(40, "When user hovers over bar chart viewer", () => hoverOver(page, el("bar chart viewer")));
      await session.step(41, "Then x-slider range slider in bar chart viewer should be visible", () => shouldBe(page, el("x-slider range slider in bar chart viewer"), "visible"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
      await session.step(43, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
    });
    await run.scenario("A logarithmic axis keeps the clipping and raises no error", async () => {
      await session.step(46, "Then \"Axis Type\" property of bar chart viewer should be \"linear\"", () => propertyShouldBe(page, "Axis Type", el("bar chart viewer"), "linear"));
      await session.step(47, "When user sets \"Axis Type\" property of bar chart viewer to \"logarithmic\"", () => setProperty(page, "Axis Type", el("bar chart viewer"), "logarithmic"));
      await session.step(48, "Then \"Axis Type\" property of bar chart viewer should be \"logarithmic\"", () => propertyShouldBe(page, "Axis Type", el("bar chart viewer"), "logarithmic"));
      await session.step(49, "And bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
      await session.step(50, "And \"Min\" property of bar chart viewer should be \"10\"", () => propertyShouldBe(page, "Min", el("bar chart viewer"), "10"));
      await session.step(51, "And \"Show Clipped Bar Indicators\" property of bar chart viewer should be \"true\"", () => propertyShouldBe(page, "Show Clipped Bar Indicators", el("bar chart viewer"), "true"));
      await session.step(52, "And the \"clipped bars\" reading of bar chart viewer should be 4", () => readingIs(page, "clipped bars", el("bar chart viewer"), 4));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
      await session.step(54, "When user sets \"Axis Type\" property of bar chart viewer to \"linear\"", () => setProperty(page, "Axis Type", el("bar chart viewer"), "linear"));
      await session.step(55, "Then bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clearing Min and Max unclips every bar", async () => {
      await session.step(59, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Min",""],["Max",""]]), [["Min",""],["Max",""]]);
      await session.step(62, "Then \"Axis Type\" property of bar chart viewer should be \"linear\"", () => propertyShouldBe(page, "Axis Type", el("bar chart viewer"), "linear"));
      await session.step(63, "And the \"clipped bars\" reading of bar chart viewer should be 0", () => readingIs(page, "clipped bars", el("bar chart viewer"), 0));
      await session.step(64, "And bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
      await session.step(65, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
      await session.step(66, "Then x-slider range slider in bar chart viewer should be hidden", () => shouldBe(page, el("x-slider range slider in bar chart viewer"), "hidden"));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
