/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map, core#1202]
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
import {addViewer, hasArea, hasNoArea, readingHigherThanRemembered, readingIs, readingReads, rememberReading, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map known failures", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map-known-failures.feature", import.meta.url);
  test("A column hidden under Max Heatmap Columns stays hidden when the cap is raised (core#1202)", {tag: ["@viewers", "@realizes:viewers.heat-map", "@edge-cases", "@known-failure:32", "@realizes:core#1202"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(13, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(14, "And the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
    await knownFailure(async () => {
      await session.step(18, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"5\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "5"));
      await session.step(19, "Then the \"columns shown\" reading of heat map viewer should be 5", () => readingIs(page, "columns shown", el("heat map viewer"), 5));
      await session.step(20, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
      await session.step(21, "When user clicks on settings icon of heat map viewer", () => clickOn(page, el("settings icon of heat map viewer")));
      await session.step(22, "Given \"Columns\" category in context panel is expanded", () => isExpanded(page, el("\"Columns\" category in context panel")));
      await session.step(23, "When user clicks on first button in \"Columns\" property in context panel", () => clickOn(page, el("first button in \"Columns\" property in context panel")));
      await session.step(24, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
      await session.step(25, "When user toggles the \"AGE\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "AGE", el("Order or Hide Columns dialog")));
      await session.step(26, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
      await session.step(27, "Then the \"columns shown\" reading of heat map viewer should be 4", () => readingIs(page, "columns shown", el("heat map viewer"), 4));
      await session.step(28, "And heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
      await session.step(29, "When user remembers the \"columns shown\" reading of heat map viewer", () => rememberReading(page, "columns shown", el("heat map viewer")));
      await session.step(30, "And user enters \"6\" into \"Max Heatmap Columns\" property in context panel", () => enterInto(page, "6", el("\"Max Heatmap Columns\" property in context panel")));
      await session.step(31, "Then the \"columns shown\" reading of heat map viewer should be higher than remembered", () => readingHigherThanRemembered(page, "columns shown", el("heat map viewer")));
      await session.step(32, "And heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
    }, 32);
  });
});
