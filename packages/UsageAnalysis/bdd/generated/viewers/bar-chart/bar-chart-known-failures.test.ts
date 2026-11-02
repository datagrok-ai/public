/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart, core#1200, core#1174]
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
import {clickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, pickFromAreaContextMenu, propertyShouldBe, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart known failures", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-known-failures.feature", import.meta.url);
  test("Max Bar Height 0 typed in the context panel outside its declared bounds is not taken (core#1200) [caption=Max Bar Height, typed=0, kept=50]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:26", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"]]), [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"]]);
    await session.step(17, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await knownFailure(async () => {
      await session.step(21, "Given user sets \"Max Bar Height\" property of bar chart viewer to \"50\"", () => setProperty(page, "Max Bar Height", el("bar chart viewer"), "50"));
      await session.step(22, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
      await session.step(23, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(24, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(25, "When user enters \"0\" into \"Max Bar Height\" property in context panel", () => enterInto(page, "0", el("\"Max Bar Height\" property in context panel")));
      await session.step(26, "Then \"Max Bar Height\" property of bar chart viewer should be \"50\"", () => propertyShouldBe(page, "Max Bar Height", el("bar chart viewer"), "50"));
    }, 26);
  });
  test("Max Bar Height -5 typed in the context panel outside its declared bounds is not taken (core#1200) [caption=Max Bar Height, typed=-5, kept=50]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:26", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"]]), [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"]]);
    await session.step(17, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await knownFailure(async () => {
      await session.step(21, "Given user sets \"Max Bar Height\" property of bar chart viewer to \"50\"", () => setProperty(page, "Max Bar Height", el("bar chart viewer"), "50"));
      await session.step(22, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
      await session.step(23, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(24, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(25, "When user enters \"-5\" into \"Max Bar Height\" property in context panel", () => enterInto(page, "-5", el("\"Max Bar Height\" property in context panel")));
      await session.step(26, "Then \"Max Bar Height\" property of bar chart viewer should be \"50\"", () => propertyShouldBe(page, "Max Bar Height", el("bar chart viewer"), "50"));
    }, 26);
  });
  test("Max Category Width 0 typed in the context panel outside its declared bounds is not taken (core#1200) [caption=Max Category Width, typed=0, kept=300]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:26", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"]]), [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"]]);
    await session.step(17, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await knownFailure(async () => {
      await session.step(21, "Given user sets \"Max Category Width\" property of bar chart viewer to \"300\"", () => setProperty(page, "Max Category Width", el("bar chart viewer"), "300"));
      await session.step(22, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
      await session.step(23, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(24, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(25, "When user enters \"0\" into \"Max Category Width\" property in context panel", () => enterInto(page, "0", el("\"Max Category Width\" property in context panel")));
      await session.step(26, "Then \"Max Category Width\" property of bar chart viewer should be \"300\"", () => propertyShouldBe(page, "Max Category Width", el("bar chart viewer"), "300"));
    }, 26);
  });
  test("Category Value Width 0 typed in the context panel outside its declared bounds is not taken (core#1200) [caption=Category Value Width, typed=0, kept=50]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:26", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"]]), [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"]]);
    await session.step(17, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await knownFailure(async () => {
      await session.step(21, "Given user sets \"Category Value Width\" property of bar chart viewer to \"50\"", () => setProperty(page, "Category Value Width", el("bar chart viewer"), "50"));
      await session.step(22, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
      await session.step(23, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(24, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(25, "When user enters \"0\" into \"Category Value Width\" property in context panel", () => enterInto(page, "0", el("\"Category Value Width\" property in context panel")));
      await session.step(26, "Then \"Category Value Width\" property of bar chart viewer should be \"50\"", () => propertyShouldBe(page, "Category Value Width", el("bar chart viewer"), "50"));
    }, 26);
  });
  test("Removing the column of the category a click filtered by gives the rows back (core#1174)", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:41", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"]]), [["Split","RACE"],["Value","AGE"],["Value Aggr Type","avg"]]);
    await session.step(17, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await knownFailure(async () => {
      await session.step(37, "When user sets \"onClick\" property of bar chart viewer to \"Filter\"", () => setProperty(page, "onClick", el("bar chart viewer"), "Filter"));
      await session.step(38, "And user clicks on the \"bar Asian\" area of bar chart viewer", () => clickArea(page, "bar Asian", el("bar chart viewer")));
      await session.step(39, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(40, "When user picks \"Remove\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header RACE", el("grid")));
      await session.step(41, "Then all rows should pass the filter", () => filterPassesAll(page));
    }, 41);
  });
});
