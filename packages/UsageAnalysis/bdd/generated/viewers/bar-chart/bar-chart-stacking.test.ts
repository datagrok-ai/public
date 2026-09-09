/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-stacking.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {barsDiffer, barsEqual} from '../../../bindings/bar-chart.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, shouldHaveItems} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaMoreInk, areaPainted, hasArea, noBalloons, noErrors, notRepainted, painted, propertyShouldBe, readingHigher, readingIs, repaintedBy, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart stacking and relative values", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-stacking.feature", import.meta.url);
  test("Bar chart stacking and relative values", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(12, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Primary Series Name"],["Value","Chemical Space X"],["Value Aggr Type","sum"]]));
    await session.step(16, "Then the \"bars\" reading of bar chart viewer should be 5", () => readingIs(page, "bars", el("bar chart viewer"), 5));
    await session.step(17, "And the \"stack segments\" reading of bar chart viewer should be 0", () => readingIs(page, "stack segments", el("bar chart viewer"), 0));
    await session.step(18, "And \"Relative Values\" property of bar chart viewer should be \"false\"", () => propertyShouldBe(page, "Relative Values", el("bar chart viewer"), "false"));
    await run.scenario("The stacked bars render on negative sums without errors", async () => {
      await session.step(21, "When user sets \"Stack\" property of bar chart viewer to \"Scaffold Names\"", () => setProperty(page, "Stack", el("bar chart viewer"), "Scaffold Names"));
      await session.step(22, "Then the \"stack segments\" reading of bar chart viewer should be higher than before", () => readingHigher(page, "stack segments", el("bar chart viewer")));
      await session.step(23, "And legend of bar chart viewer should be visible", () => shouldBe(page, el("legend of bar chart viewer"), "visible"));
      await session.step(24, "And legend of bar chart viewer should have 6 items", () => shouldHaveItems(page, el("legend of bar chart viewer"), 6));
      await session.step(25, "And bar chart viewer should have a \"bar Triazoles | TRISUBSTITUTED\" area", () => hasArea(page, el("bar chart viewer"), "bar Triazoles | TRISUBSTITUTED"));
      await session.step(26, "And the \"bar Triazoles | TRISUBSTITUTED\" area of bar chart viewer should be painted", () => areaPainted(page, "bar Triazoles | TRISUBSTITUTED", el("bar chart viewer")));
      await session.step(27, "And bar chart viewer should have a \"bar Pyrrolidines | AMINOPYRROLIDINE\" area", () => hasArea(page, el("bar chart viewer"), "bar Pyrrolidines | AMINOPYRROLIDINE"));
      await session.step(28, "And the bars of bar chart viewer should differ in length", () => barsDiffer(page, el("bar chart viewer")));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
      await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Relative Values with a Stack column normalizes the bars", async () => {
      await session.step(33, "When user sets \"Relative Values\" property of bar chart viewer to \"true\"", () => setProperty(page, "Relative Values", el("bar chart viewer"), "true"));
      await session.step(34, "Then bar chart viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("bar chart viewer"), 2000));
      await session.step(35, "And the bars of bar chart viewer should be of equal length", () => barsEqual(page, el("bar chart viewer")));
      await session.step(36, "And the \"bar Pyrrolidines | AMINOPYRROLIDINE\" area of bar chart viewer should have more ink than before", () => areaMoreInk(page, "bar Pyrrolidines | AMINOPYRROLIDINE", el("bar chart viewer")));
      await session.step(37, "And bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Relative Values off restores the absolute lengths", async () => {
      await session.step(41, "When user sets \"Relative Values\" property of bar chart viewer to \"false\"", () => setProperty(page, "Relative Values", el("bar chart viewer"), "false"));
      await session.step(42, "Then bar chart viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("bar chart viewer"), 2000));
      await session.step(43, "And the bars of bar chart viewer should differ in length", () => barsDiffer(page, el("bar chart viewer")));
      await session.step(44, "And \"Stack\" property of bar chart viewer should be \"Scaffold Names\"", () => propertyShouldBe(page, "Stack", el("bar chart viewer"), "Scaffold Names"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the Stack collapses the bars to single segments", async () => {
      await session.step(48, "When user sets \"Stack\" property of bar chart viewer to \"\"", () => setProperty(page, "Stack", el("bar chart viewer"), ""));
      await session.step(49, "Then the \"stack segments\" reading of bar chart viewer should be 0", () => readingIs(page, "stack segments", el("bar chart viewer"), 0));
      await session.step(50, "And legend of bar chart viewer should be hidden", () => shouldBe(page, el("legend of bar chart viewer"), "hidden"));
      await session.step(51, "And bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Relative Values alone is inert", async () => {
      await session.step(55, "When user sets \"Relative Values\" property of bar chart viewer to \"true\"", () => setProperty(page, "Relative Values", el("bar chart viewer"), "true"));
      await session.step(56, "Then \"Stack\" property of bar chart viewer should be \"\"", () => propertyShouldBe(page, "Stack", el("bar chart viewer"), ""));
      await session.step(57, "And bar chart viewer should not have repainted", () => notRepainted(page, el("bar chart viewer")));
      await session.step(58, "And the bars of bar chart viewer should differ in length", () => barsDiffer(page, el("bar chart viewer")));
      await session.step(59, "And legend of bar chart viewer should be hidden", () => shouldBe(page, el("legend of bar chart viewer"), "hidden"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Stack column activates it at once, and again after it is removed", async () => {
      await session.step(63, "When user sets \"Stack\" property of bar chart viewer to \"Scaffold Names\"", () => setProperty(page, "Stack", el("bar chart viewer"), "Scaffold Names"));
      await session.step(64, "Then the bars of bar chart viewer should be of equal length", () => barsEqual(page, el("bar chart viewer")));
      await session.step(65, "And legend of bar chart viewer should have 6 items", () => shouldHaveItems(page, el("legend of bar chart viewer"), 6));
      await session.step(66, "And the \"stack segments\" reading of bar chart viewer should be higher than before", () => readingHigher(page, "stack segments", el("bar chart viewer")));
      await session.step(67, "When user sets \"Stack\" property of bar chart viewer to \"\"", () => setProperty(page, "Stack", el("bar chart viewer"), ""));
      await session.step(68, "Then the bars of bar chart viewer should differ in length", () => barsDiffer(page, el("bar chart viewer")));
      await session.step(69, "And legend of bar chart viewer should be hidden", () => shouldBe(page, el("legend of bar chart viewer"), "hidden"));
      await session.step(70, "When user sets \"Stack\" property of bar chart viewer to \"Scaffold Names\"", () => setProperty(page, "Stack", el("bar chart viewer"), "Scaffold Names"));
      await session.step(71, "Then the bars of bar chart viewer should be of equal length", () => barsEqual(page, el("bar chart viewer")));
      await session.step(72, "And legend of bar chart viewer should have 6 items", () => shouldHaveItems(page, el("legend of bar chart viewer"), 6));
      await session.step(73, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Stack",""],["Relative Values","false"]]));
      await session.step(76, "Then \"Relative Values\" property of bar chart viewer should be \"false\"", () => propertyShouldBe(page, "Relative Values", el("bar chart viewer"), "false"));
      await session.step(77, "And the bars of bar chart viewer should differ in length", () => barsDiffer(page, el("bar chart viewer")));
      await session.step(78, "And legend of bar chart viewer should be hidden", () => shouldBe(page, el("legend of bar chart viewer"), "hidden"));
      await session.step(79, "And bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
