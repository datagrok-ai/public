/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-category-order.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart]
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
import {clickOn, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, readingBetween, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart category order", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-category-order.feature", import.meta.url);
  test("Pie chart category order", {tag: ["@journey", "@viewers", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"]]), [["Category","RACE"]]);
    await session.step(16, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await run.scenario("Category Order orders the slices by row count, and Default gives the column's order", async () => {
      await session.step(19, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(20, "When user sets \"Pie Sort Type\" property of pie chart viewer to \"custom\"", () => setProperty(page, "Pie Sort Type", el("pie chart viewer"), "custom"));
      await session.step(21, "And user clicks on settings icon of pie chart viewer", () => clickOn(page, el("settings icon of pie chart viewer")));
      await session.step(22, "And \"Data\" category in context panel is expanded", () => isExpanded(page, el("\"Data\" category in context panel")));
      await session.step(23, "And user clicks on button in \"Category Order\" property in context panel", () => clickOn(page, el("button in \"Category Order\" property in context panel")));
      await session.step(24, "Then \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("\"Order categories: RACE\" dialog"), "visible"));
      await session.step(25, "And the \"start angle of Asian\" reading of pie chart viewer should be between -0.01 and 0.01", () => readingBetween(page, "start angle of Asian", el("pie chart viewer"), -0.01, 0.01));
      await session.step(26, "And the \"start angle of Black\" reading of pie chart viewer should be between 5.39 and 5.41", () => readingBetween(page, "start angle of Black", el("pie chart viewer"), 5.39, 5.41));
      await session.step(27, "And the \"start angle of Caucasian\" reading of pie chart viewer should be between 15.11 and 15.13", () => readingBetween(page, "start angle of Caucasian", el("pie chart viewer"), 15.11, 15.13));
      await session.step(28, "When user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(29, "And user selects \"Descending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Descending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(30, "Then the \"start angle of Caucasian\" reading of pie chart viewer should be between -0.01 and 0.01", () => readingBetween(page, "start angle of Caucasian", el("pie chart viewer"), -0.01, 0.01));
      await session.step(31, "And the \"start angle of Other\" reading of pie chart viewer should be between 322.55 and 322.57", () => readingBetween(page, "start angle of Other", el("pie chart viewer"), 322.55, 322.57));
      await session.step(32, "And the \"start angle of Black\" reading of pie chart viewer should be between 344.87 and 344.89", () => readingBetween(page, "start angle of Black", el("pie chart viewer"), 344.87, 344.89));
      await session.step(33, "And the \"start angle of Asian\" reading of pie chart viewer should be between 354.59 and 354.61", () => readingBetween(page, "start angle of Asian", el("pie chart viewer"), 354.59, 354.61));
      await session.step(34, "When user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(35, "Then the \"start angle of Asian\" reading of pie chart viewer should be between -0.01 and 0.01", () => readingBetween(page, "start angle of Asian", el("pie chart viewer"), -0.01, 0.01));
      await session.step(36, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(37, "And user sets \"Pie Sort Type\" property of pie chart viewer to \"by value\"", () => setProperty(page, "Pie Sort Type", el("pie chart viewer"), "by value"));
      await session.step(38, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
