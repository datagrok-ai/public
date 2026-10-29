/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map]
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
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {deleteSelected, noneSelected, removeColumn, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickAreaHolding, doubleClickArea, dragZoomOverArea, hasArea, hasNoArea, noErrors, painted, readingBetween, readingHigher, readingIs, reportsNoError} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map edge cases — every row deleted, a column removed", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map-edge-cases.feature", import.meta.url);
  test("Once every row of the table is deleted the map logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(12, "Then the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
    await session.step(13, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(14, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(17, "When user selects all rows", () => selectAllRows(page));
    await session.step(18, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(19, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(20, "And the \"rows\" reading of heat map viewer should be 0", () => readingIs(page, "rows", el("heat map viewer"), 0));
    await session.step(21, "And heat map viewer should report no error", () => reportsNoError(page, el("heat map viewer")));
    await session.step(22, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column removed from the table leaves the map without its band and keeps the others", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(12, "Then the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
    await session.step(13, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(14, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(25, "Then heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(26, "When user removes \"AGE\" column", () => removeColumn(page, "AGE"));
    await session.step(27, "Then the table should not have a column \"AGE\"", () => hasNoColumn(page, "AGE"));
    await session.step(28, "And the \"columns shown\" reading of heat map viewer should be 10", () => readingIs(page, "columns shown", el("heat map viewer"), 10));
    await session.step(29, "And heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
    await session.step(30, "And heat map viewer should have a \"column WEIGHT\" area", () => hasArea(page, el("heat map viewer"), "column WEIGHT"));
    await session.step(31, "And heat map viewer should report no error", () => reportsNoError(page, el("heat map viewer")));
    await session.step(32, "And no errors should have been logged", () => noErrors(page));
  });
  test("An Alt-click without moving zooms the map to one row, an Alt-drag to the rows it spans", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(12, "Then the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
    await session.step(13, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(14, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await session.step(35, "Then the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
    await session.step(36, "When user clicks on the \"column AGE\" area of heat map viewer holding Alt", () => clickAreaHolding(page, "column AGE", el("heat map viewer"), "Alt"));
    await session.step(37, "Then the \"y scroll span\" reading of heat map viewer should be 0", () => readingIs(page, "y scroll span", el("heat map viewer"), 0));
    await session.step(38, "And the \"row height\" reading of heat map viewer should be higher than before", () => readingHigher(page, "row height", el("heat map viewer")));
    await session.step(39, "And no rows should be selected", () => noneSelected(page));
    await session.step(40, "When user double-clicks on the \"y scroll slider\" area of heat map viewer", () => doubleClickArea(page, "y scroll slider", el("heat map viewer")));
    await session.step(41, "Then the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
    await session.step(42, "When user drags a zoom box over the \"column AGE\" area of heat map viewer", () => dragZoomOverArea(page, "column AGE", el("heat map viewer")));
    await session.step(43, "Then the \"y scroll span\" reading of heat map viewer should be between 0.1 and 0.99", () => readingBetween(page, "y scroll span", el("heat map viewer"), 0.1, 0.99));
    await session.step(44, "And no rows should be selected", () => noneSelected(page));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
});
