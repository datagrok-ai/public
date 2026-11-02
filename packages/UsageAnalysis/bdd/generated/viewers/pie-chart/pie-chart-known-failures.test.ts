/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart, core#1200]
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, propertyShouldBe, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart — the known failures", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-known-failures.feature", import.meta.url);
  test("Max Radius 600 typed in the context panel outside its declared bounds is not taken (core#1200) [caption=Max Radius, typed=600, kept=150]", {tag: ["@viewers", "@realizes:viewers.pie-chart", "@edge-cases", "@known-failure:26", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"]]), [["Category","RACE"]]);
    await session.step(17, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await knownFailure(async () => {
      await session.step(21, "Given user sets \"Max Radius\" property of pie chart viewer to \"150\"", () => setProperty(page, "Max Radius", el("pie chart viewer"), "150"));
      await session.step(22, "When user clicks on settings icon of pie chart viewer", () => clickOn(page, el("settings icon of pie chart viewer")));
      await session.step(23, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(24, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(25, "When user enters \"600\" into \"Max Radius\" property in context panel", () => enterInto(page, "600", el("\"Max Radius\" property in context panel")));
      await session.step(26, "Then \"Max Radius\" property of pie chart viewer should be \"150\"", () => propertyShouldBe(page, "Max Radius", el("pie chart viewer"), "150"));
    }, 26);
  });
  test("Start Angle -30 typed in the context panel outside its declared bounds is not taken (core#1200) [caption=Start Angle, typed=-30, kept=0]", {tag: ["@viewers", "@realizes:viewers.pie-chart", "@edge-cases", "@known-failure:26", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"]]), [["Category","RACE"]]);
    await session.step(17, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await knownFailure(async () => {
      await session.step(21, "Given user sets \"Start Angle\" property of pie chart viewer to \"0\"", () => setProperty(page, "Start Angle", el("pie chart viewer"), "0"));
      await session.step(22, "When user clicks on settings icon of pie chart viewer", () => clickOn(page, el("settings icon of pie chart viewer")));
      await session.step(23, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(24, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(25, "When user enters \"-30\" into \"Start Angle\" property in context panel", () => enterInto(page, "-30", el("\"Start Angle\" property in context panel")));
      await session.step(26, "Then \"Start Angle\" property of pie chart viewer should be \"0\"", () => propertyShouldBe(page, "Start Angle", el("pie chart viewer"), "0"));
    }, 26);
  });
  test("Start Angle 400 typed in the context panel outside its declared bounds is not taken (core#1200) [caption=Start Angle, typed=400, kept=0]", {tag: ["@viewers", "@realizes:viewers.pie-chart", "@edge-cases", "@known-failure:26", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"]]), [["Category","RACE"]]);
    await session.step(17, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await knownFailure(async () => {
      await session.step(21, "Given user sets \"Start Angle\" property of pie chart viewer to \"0\"", () => setProperty(page, "Start Angle", el("pie chart viewer"), "0"));
      await session.step(22, "When user clicks on settings icon of pie chart viewer", () => clickOn(page, el("settings icon of pie chart viewer")));
      await session.step(23, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(24, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(25, "When user enters \"400\" into \"Start Angle\" property in context panel", () => enterInto(page, "400", el("\"Start Angle\" property in context panel")));
      await session.step(26, "Then \"Start Angle\" property of pie chart viewer should be \"0\"", () => propertyShouldBe(page, "Start Angle", el("pie chart viewer"), "0"));
    }, 26);
  });
});
