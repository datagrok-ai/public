/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map-edge-cases.feature
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
import {clickOn, expand, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {noneSelected, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickAreaHolding, doubleClickArea, dragZoomOverArea, hasArea, hasNoArea, noErrors, painted, readingBetween, readingHigher, readingIs, reportsNoError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map edge cases — a column removed", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map-edge-cases.feature", import.meta.url);
  test("A column removed from the table leaves the map without its band and keeps the others", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(12, "Then the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
    await session.step(13, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(14, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(17, "Then heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(18, "When user removes \"AGE\" column", () => removeColumn(page, "AGE"));
    await session.step(19, "Then the table should not have a column \"AGE\"", () => hasNoColumn(page, "AGE"));
    await session.step(20, "And the \"columns shown\" reading of heat map viewer should be 10", () => readingIs(page, "columns shown", el("heat map viewer"), 10));
    await session.step(21, "And heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
    await session.step(22, "And heat map viewer should have a \"column WEIGHT\" area", () => hasArea(page, el("heat map viewer"), "column WEIGHT"));
    await session.step(23, "And heat map viewer should report no error", () => reportsNoError(page, el("heat map viewer")));
    await session.step(24, "And no errors should have been logged", () => noErrors(page));
  });
  test("An Alt-click without moving zooms the map to one row, an Alt-drag to the rows it spans", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(12, "Then the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
    await session.step(13, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(14, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(27, "Then the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
    await session.step(28, "When user clicks on the \"column AGE\" area of heat map viewer holding Alt", () => clickAreaHolding(page, "column AGE", el("heat map viewer"), "Alt"));
    await session.step(29, "Then the \"y scroll span\" reading of heat map viewer should be 0", () => readingIs(page, "y scroll span", el("heat map viewer"), 0));
    await session.step(30, "And the \"row height\" reading of heat map viewer should be higher than before", () => readingHigher(page, "row height", el("heat map viewer")));
    await session.step(31, "And no rows should be selected", () => noneSelected(page));
    await session.step(32, "When user double-clicks on the \"y scroll slider\" area of heat map viewer", () => doubleClickArea(page, "y scroll slider", el("heat map viewer")));
    await session.step(33, "Then the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
    await session.step(34, "When user drags a zoom box over the \"column AGE\" area of heat map viewer", () => dragZoomOverArea(page, "column AGE", el("heat map viewer")));
    await session.step(35, "Then the \"y scroll span\" reading of heat map viewer should be between 0.1 and 0.99", () => readingBetween(page, "y scroll span", el("heat map viewer"), 0.1, 0.99));
    await session.step(36, "And no rows should be selected", () => noneSelected(page));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column hidden under Max Heatmap Columns stays hidden when the limit is raised (core#1202)", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.heat-map", "@known-failure:53", "@realizes:core#1202"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(12, "Then the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
    await session.step(13, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(14, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await knownFailure(async () => {
      await session.step(41, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"5\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "5"));
      await session.step(42, "Then the \"columns shown\" reading of heat map viewer should be 5", () => readingIs(page, "columns shown", el("heat map viewer"), 5));
      await session.step(43, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
      await session.step(44, "When user clicks on settings icon of heat map viewer", () => clickOn(page, el("settings icon of heat map viewer")));
      await session.step(45, "And user expands \"Columns\" category in context panel", () => expand(page, el("\"Columns\" category in context panel")));
      await session.step(46, "And user clicks on first button in \"Columns\" property in context panel", () => clickOn(page, el("first button in \"Columns\" property in context panel")));
      await session.step(47, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
      await session.step(48, "When user toggles the \"AGE\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "AGE", el("Order or Hide Columns dialog")));
      await session.step(49, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
      await session.step(50, "Then the \"columns shown\" reading of heat map viewer should be 4", () => readingIs(page, "columns shown", el("heat map viewer"), 4));
      await session.step(51, "And heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
      await session.step(52, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"6\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "6"));
      await session.step(53, "Then heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    }, 53);
  });
});
