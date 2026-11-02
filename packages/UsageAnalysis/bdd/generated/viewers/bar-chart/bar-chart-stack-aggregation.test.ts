/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-stack-aggregation.feature
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
import {clickOn, shouldBe, shouldContainText, shouldHaveItems, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterBetween, filterPasses, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hoverArea, noErrors, propertyShouldBe, readingIs, repainted, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart stack aggregation and datetime split", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-stack-aggregation.feature", import.meta.url);
  test("Bar chart stack aggregation and datetime split", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"],["Stack","SEX"]]), [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"],["Stack","SEX"]]);
    await session.step(17, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await run.scenario("A Stack under avg draws no segments and no legend", async () => {
      await session.step(20, "Then the \"stack segments\" reading of bar chart viewer should be 0", () => readingIs(page, "stack segments", el("bar chart viewer"), 0));
      await session.step(21, "And legend of bar chart viewer should be hidden", () => shouldBe(page, el("legend of bar chart viewer"), "hidden"));
      await session.step(22, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sum builds the stack and its legend", async () => {
      await session.step(25, "When user sets \"Value Aggr Type\" property of bar chart viewer to \"sum\"", () => setProperty(page, "Value Aggr Type", el("bar chart viewer"), "sum"));
      await session.step(26, "Then the \"stack segments\" reading of bar chart viewer should be 8", () => readingIs(page, "stack segments", el("bar chart viewer"), 8));
      await session.step(27, "And legend of bar chart viewer should be visible", () => shouldBe(page, el("legend of bar chart viewer"), "visible"));
      await session.step(28, "And legend of bar chart viewer should have 2 items", () => shouldHaveItems(page, el("legend of bar chart viewer"), 2));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Count keeps the stack; removing the Stack collapses it", async () => {
      await session.step(32, "When user sets \"Value Aggr Type\" property of bar chart viewer to \"count\"", () => setProperty(page, "Value Aggr Type", el("bar chart viewer"), "count"));
      await session.step(33, "Then the \"stack segments\" reading of bar chart viewer should be 8", () => readingIs(page, "stack segments", el("bar chart viewer"), 8));
      await session.step(34, "And legend of bar chart viewer should have 2 items", () => shouldHaveItems(page, el("legend of bar chart viewer"), 2));
      await session.step(35, "When user sets \"Stack\" property of bar chart viewer to \"\"", () => setProperty(page, "Stack", el("bar chart viewer"), ""));
      await session.step(36, "Then the \"stack segments\" reading of bar chart viewer should be 0", () => readingIs(page, "stack segments", el("bar chart viewer"), 0));
      await session.step(37, "And legend of bar chart viewer should be hidden", () => shouldBe(page, el("legend of bar chart viewer"), "hidden"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A string Split has no Split Map", async () => {
      await session.step(41, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
      await session.step(42, "Then \"Split Map\" property in context panel should be disabled", () => shouldBe(page, el("\"Split Map\" property in context panel"), "disabled"));
      await session.step(43, "And the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A datetime Split enables the Split Map, which re-categorizes the bars", async () => {
      await session.step(47, "When user sets \"Split\" property of bar chart viewer to \"STARTED\"", () => setProperty(page, "Split", el("bar chart viewer"), "STARTED"));
      await session.step(48, "Then \"Split Map\" property in context panel should be enabled", () => shouldBe(page, el("\"Split Map\" property in context panel"), "enabled"));
      await session.step(49, "And Split column input in bar chart viewer should contain text \"STARTED\"", () => shouldContainText(page, el("Split column input in bar chart viewer"), "STARTED"));
      await session.step(50, "And Split column input in bar chart viewer should not contain text \"RACE\"", () => shouldNotContainText(page, el("Split column input in bar chart viewer"), "RACE"));
      await session.step(51, "When user sets \"Split Map\" property of bar chart viewer to \"month\"", () => setProperty(page, "Split Map", el("bar chart viewer"), "month"));
      await session.step(52, "Then the \"bars\" reading of bar chart viewer should be 12", () => readingIs(page, "bars", el("bar chart viewer"), 12));
      await session.step(53, "And bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
      await session.step(54, "When user sets \"Split Map\" property of bar chart viewer to \"year\"", () => setProperty(page, "Split Map", el("bar chart viewer"), "year"));
      await session.step(55, "Then the \"bars\" reading of bar chart viewer should be 3", () => readingIs(page, "bars", el("bar chart viewer"), 3));
      await session.step(56, "And bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
      await session.step(57, "And \"Split\" property of bar chart viewer should be \"STARTED\"", () => propertyShouldBe(page, "Split", el("bar chart viewer"), "STARTED"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
      await session.step(59, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Split","RACE"],["Value Aggr Type","avg"]]), [["Split","RACE"],["Value Aggr Type","avg"]]);
      await session.step(62, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    });
    await run.scenario("A vertical bar chart with every row filtered out draws no bar and logs nothing (core#974)", async () => {
      await session.step(66, "When user sets \"Orientation\" property of bar chart viewer to \"vertical\"", () => setProperty(page, "Orientation", el("bar chart viewer"), "vertical"));
      await session.step(67, "And user filters rows where \"AGE\" is between 200 and 300", () => filterBetween(page, "AGE", 200, 300));
      await session.step(68, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(69, "And bar chart viewer should show 0 rows", () => showsRows(page, el("bar chart viewer"), 0));
      await session.step(70, "And the \"bars\" reading of bar chart viewer should be 0", () => readingIs(page, "bars", el("bar chart viewer"), 0));
      await session.step(71, "When user hovers over the \"view\" area of bar chart viewer", () => hoverArea(page, "view", el("bar chart viewer")));
      await session.step(72, "Then no errors should have been logged", () => noErrors(page));
      await session.step(73, "When user resets the filter", () => resetFilter(page));
      await session.step(74, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
      await session.step(75, "And bar chart viewer should show 1000 rows", () => showsRows(page, el("bar chart viewer"), 1000));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
