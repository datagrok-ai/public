/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-sort-orientation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {barsAscend, barsDescend, barsStacked, tallestReaches} from '../../../bindings/bar-chart.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, shouldHaveItems} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, noBalloons, noErrors, painted, propertiesShouldBe, propertyShouldBe, readingHigher, readingIs, repainted, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaHangsBelow} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart sorting and orientation", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-sort-orientation.feature", import.meta.url);
  test("Bar chart sorting and orientation", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(12, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Primary Series Name"],["Value","Chemical Space X"],["Value Aggr Type","sum"]]));
    await session.step(16, "Then the \"bars\" reading of bar chart viewer should be 5", () => readingIs(page, "bars", el("bar chart viewer"), 5));
    await session.step(17, "And the bars of bar chart viewer should lie one under another", () => barsStacked(page, el("bar chart viewer")));
    await run.scenario("Vertical and descending by value puts the tallest bar on the left", async () => {
      await session.step(20, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Orientation","vertical"],["Bar Sort Type","by value"],["Bar Sort Order","desc"]]));
      await session.step(24, "Then bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
      await session.step(25, "And the bars of bar chart viewer should descend from left to right", () => barsDescend(page, el("bar chart viewer")));
      await session.step(26, "And the \"bar Triazoles\" area of bar chart viewer should hang below the \"bar Aminopiperidines\" area", () => areaHangsBelow(page, "bar Triazoles", el("bar chart viewer"), "bar Aminopiperidines"));
      await session.step(27, "And no errors should have been logged", () => noErrors(page));
      await session.step(28, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Stacking holds on a negative-sum value", async () => {
      await session.step(31, "When user sets \"Legend Visibility\" property of bar chart viewer to \"Always\"", () => setProperty(page, "Legend Visibility", el("bar chart viewer"), "Always"));
      await session.step(32, "Then legend of bar chart viewer should be hidden", () => shouldBe(page, el("legend of bar chart viewer"), "hidden"));
      await session.step(33, "When user sets \"Stack\" property of bar chart viewer to \"Stereo Category\"", () => setProperty(page, "Stack", el("bar chart viewer"), "Stereo Category"));
      await session.step(34, "Then legend of bar chart viewer should be visible", () => shouldBe(page, el("legend of bar chart viewer"), "visible"));
      await session.step(35, "And legend of bar chart viewer should have 5 items", () => shouldHaveItems(page, el("legend of bar chart viewer"), 5));
      await session.step(36, "And the \"stack segments\" reading of bar chart viewer should be higher than before", () => readingHigher(page, "stack segments", el("bar chart viewer")));
      await session.step(37, "And bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
      await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(40, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Stack",""],["Legend Visibility","Auto"]]));
      await session.step(43, "Then legend of bar chart viewer should be hidden", () => shouldBe(page, el("legend of bar chart viewer"), "hidden"));
      await session.step(44, "And the \"stack segments\" reading of bar chart viewer should be 0", () => readingIs(page, "stack segments", el("bar chart viewer"), 0));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ascending swaps the tall side, horizontal puts the bars back under each other", async () => {
      await session.step(48, "When user sets \"Bar Sort Order\" property of bar chart viewer to \"asc\"", () => setProperty(page, "Bar Sort Order", el("bar chart viewer"), "asc"));
      await session.step(49, "Then bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
      await session.step(50, "And the bars of bar chart viewer should ascend from left to right", () => barsAscend(page, el("bar chart viewer")));
      await session.step(51, "And the \"bar Triazoles\" area of bar chart viewer should hang below the \"bar Aminopiperidines\" area", () => areaHangsBelow(page, "bar Triazoles", el("bar chart viewer"), "bar Aminopiperidines"));
      await session.step(52, "When user sets \"Orientation\" property of bar chart viewer to \"horizontal\"", () => setProperty(page, "Orientation", el("bar chart viewer"), "horizontal"));
      await session.step(53, "Then bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
      await session.step(54, "And the bars of bar chart viewer should lie one under another", () => barsStacked(page, el("bar chart viewer")));
      await session.step(55, "And \"Bar Sort Order\" property of bar chart viewer should be \"asc\"", () => propertyShouldBe(page, "Bar Sort Order", el("bar chart viewer"), "asc"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Vertical bars with small counts fill the view", async () => {
      await session.step(59, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Value","CAST Idea ID"],["Value Aggr Type","count"]]));
      await session.step(62, "And user sets \"Orientation\" property of bar chart viewer to \"vertical\"", () => setProperty(page, "Orientation", el("bar chart viewer"), "vertical"));
      await session.step(63, "Then the tallest bar of bar chart viewer should start within the top 40% of the view", () => tallestReaches(page, el("bar chart viewer"), 40));
      await session.step(64, "And the \"bar Triazoles\" area of bar chart viewer should be painted", () => areaPainted(page, "bar Triazoles", el("bar chart viewer")));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
      await session.step(66, "When user sets \"Orientation\" property of bar chart viewer to \"horizontal\"", () => setProperty(page, "Orientation", el("bar chart viewer"), "horizontal"));
      await session.step(67, "Then bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
      await session.step(68, "And the bars of bar chart viewer should lie one under another", () => barsStacked(page, el("bar chart viewer")));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
      await session.step(70, "When user sets \"Orientation\" property of bar chart viewer to \"auto\"", () => setProperty(page, "Orientation", el("bar chart viewer"), "auto"));
      await session.step(71, "Then properties of bar chart viewer should be:", () => propertiesShouldBe(page, el("bar chart viewer"), [["Orientation","auto"],["Value","CAST Idea ID"],["Value Aggr Type","count"]]));
      await session.step(75, "And the bars of bar chart viewer should lie one under another", () => barsStacked(page, el("bar chart viewer")));
      await session.step(76, "And bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
      await session.step(78, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Value","Chemical Space X"],["Value Aggr Type","sum"],["Bar Sort Type","by category"],["Bar Sort Order","asc"]]));
    });
    run.finish();
  });
});
