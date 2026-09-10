/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-value-axis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart]
--- */
import {test} from '@playwright/test';
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
import {addViewerWith, noErrors, painted, pointerAway, propertiesShouldBe, propertyShouldBe, readingIs, repainted, repaintedBy, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart value axis range and scale", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-value-axis.feature", import.meta.url);
  test("Bar chart value axis range and scale", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(12, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Primary Series Name"],["Value","CAST Idea ID"],["Value Aggr Type","count"],["Show Clipped Bar Indicators","true"]]));
    await session.step(17, "Then the \"bars\" reading of bar chart viewer should be 5", () => readingIs(page, "bars", el("bar chart viewer"), 5));
    await session.step(18, "And the \"clipped bars\" reading of bar chart viewer should be 0", () => readingIs(page, "clipped bars", el("bar chart viewer"), 0));
    await run.scenario("Min above the shortest bars clips them and the indicators show it", async () => {
      await session.step(21, "When user sets \"Min\" property of bar chart viewer to \"10\"", () => setProperty(page, "Min", el("bar chart viewer"), "10"));
      await session.step(22, "Then \"Min\" property of bar chart viewer should be \"10\"", () => propertyShouldBe(page, "Min", el("bar chart viewer"), "10"));
      await session.step(23, "And the \"clipped bars\" reading of bar chart viewer should be 3", () => readingIs(page, "clipped bars", el("bar chart viewer"), 3));
      await session.step(24, "When user sets \"Show Clipped Bar Indicators\" property of bar chart viewer to \"false\"", () => setProperty(page, "Show Clipped Bar Indicators", el("bar chart viewer"), "false"));
      await session.step(25, "Then bar chart viewer should have repainted by at least 150 pixels", () => repaintedBy(page, el("bar chart viewer"), 150));
      await session.step(26, "When user sets \"Show Clipped Bar Indicators\" property of bar chart viewer to \"true\"", () => setProperty(page, "Show Clipped Bar Indicators", el("bar chart viewer"), "true"));
      await session.step(27, "Then bar chart viewer should have repainted by at least 150 pixels", () => repaintedBy(page, el("bar chart viewer"), 150));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Max below the tallest bar clips it too", async () => {
      await session.step(31, "When user sets \"Max\" property of bar chart viewer to \"60\"", () => setProperty(page, "Max", el("bar chart viewer"), "60"));
      await session.step(32, "Then \"Max\" property of bar chart viewer should be \"60\"", () => propertyShouldBe(page, "Max", el("bar chart viewer"), "60"));
      await session.step(33, "And the \"clipped bars\" reading of bar chart viewer should be 4", () => readingIs(page, "clipped bars", el("bar chart viewer"), 4));
      await session.step(34, "When user sets \"Show Clipped Bar Indicators\" property of bar chart viewer to \"false\"", () => setProperty(page, "Show Clipped Bar Indicators", el("bar chart viewer"), "false"));
      await session.step(35, "Then bar chart viewer should have repainted by at least 150 pixels", () => repaintedBy(page, el("bar chart viewer"), 150));
      await session.step(36, "When user sets \"Show Clipped Bar Indicators\" property of bar chart viewer to \"true\"", () => setProperty(page, "Show Clipped Bar Indicators", el("bar chart viewer"), "true"));
      await session.step(37, "Then bar chart viewer should have repainted by at least 150 pixels", () => repaintedBy(page, el("bar chart viewer"), 150));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The value-axis scroll bar shows on the constrained range", async () => {
      await session.step(41, "When user hovers over bar chart viewer", () => hoverOver(page, el("bar chart viewer")));
      await session.step(42, "Then x-slider range slider in bar chart viewer should be visible", () => shouldBe(page, el("x-slider range slider in bar chart viewer"), "visible"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
      await session.step(44, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
    });
    await run.scenario("A logarithmic axis keeps the clipping and raises no error", async () => {
      await session.step(47, "Then \"Axis Type\" property of bar chart viewer should be \"linear\"", () => propertyShouldBe(page, "Axis Type", el("bar chart viewer"), "linear"));
      await session.step(48, "When user sets \"Axis Type\" property of bar chart viewer to \"logarithmic\"", () => setProperty(page, "Axis Type", el("bar chart viewer"), "logarithmic"));
      await session.step(49, "Then \"Axis Type\" property of bar chart viewer should be \"logarithmic\"", () => propertyShouldBe(page, "Axis Type", el("bar chart viewer"), "logarithmic"));
      await session.step(50, "And bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
      await session.step(51, "And \"Min\" property of bar chart viewer should be \"10\"", () => propertyShouldBe(page, "Min", el("bar chart viewer"), "10"));
      await session.step(52, "And \"Show Clipped Bar Indicators\" property of bar chart viewer should be \"true\"", () => propertyShouldBe(page, "Show Clipped Bar Indicators", el("bar chart viewer"), "true"));
      await session.step(53, "And the \"clipped bars\" reading of bar chart viewer should be 4", () => readingIs(page, "clipped bars", el("bar chart viewer"), 4));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
      await session.step(55, "When user sets \"Axis Type\" property of bar chart viewer to \"linear\"", () => setProperty(page, "Axis Type", el("bar chart viewer"), "linear"));
      await session.step(56, "Then bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clearing Min and Max restores the full range", async () => {
      await session.step(60, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Min",""],["Max",""]]));
      await session.step(63, "Then properties of bar chart viewer should be:", () => propertiesShouldBe(page, el("bar chart viewer"), [["Min",""],["Max",""],["Axis Type","linear"]]));
      await session.step(67, "And the \"clipped bars\" reading of bar chart viewer should be 0", () => readingIs(page, "clipped bars", el("bar chart viewer"), 0));
      await session.step(68, "And bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
      await session.step(69, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
      await session.step(70, "Then x-slider range slider in bar chart viewer should be hidden", () => shouldBe(page, el("x-slider range slider in bar chart viewer"), "hidden"));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
