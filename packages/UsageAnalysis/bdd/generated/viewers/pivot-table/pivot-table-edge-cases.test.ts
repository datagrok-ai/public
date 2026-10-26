/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pivot-table]
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
import {deleteSelected, filterBetween, filterPasses, removeColumn, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table edge cases — its Group by column removed, an empty row set, every row deleted", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-edge-cases.feature", import.meta.url);
  test("The Group by column removed from the table is replaced and the pivot aggregates again", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(16, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
    await session.step(17, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
    await session.step(20, "When user removes \"DIS_POP\" column", () => removeColumn(page, "DIS_POP"));
    await session.step(21, "Then the table should not have a column \"DIS_POP\"", () => hasNoColumn(page, "DIS_POP"));
    await session.step(22, "And the \"group by\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "group by", el("pivot table viewer"), "SEVERITY"));
    await session.step(23, "And the \"error\" reading of pivot table viewer should be \"\"", () => readingReads(page, "error", el("pivot table viewer"), ""));
    await session.step(24, "And the \"aggregated rows\" reading of pivot table viewer should be 5", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 5));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out the pivot holds no aggregated row and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(16, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
    await session.step(17, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
    await session.step(28, "When user filters rows where \"AGE\" is between 200 and 300", () => filterBetween(page, "AGE", 200, 300));
    await session.step(29, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(30, "And the \"rows shown\" reading of pivot table viewer should be 0", () => readingIs(page, "rows shown", el("pivot table viewer"), 0));
    await session.step(31, "And the \"aggregated rows\" reading of pivot table viewer should be 0", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 0));
    await session.step(32, "And the \"error\" reading of pivot table viewer should be \"\"", () => readingReads(page, "error", el("pivot table viewer"), ""));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
  });
  test("Once every row of the table is deleted the pivot holds no aggregated row and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(16, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
    await session.step(17, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
    await session.step(36, "When user selects all rows", () => selectAllRows(page));
    await session.step(37, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(38, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(39, "And the \"aggregated rows\" reading of pivot table viewer should be 0", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 0));
    await session.step(40, "And the \"error\" reading of pivot table viewer should be \"\"", () => readingReads(page, "error", el("pivot table viewer"), ""));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
});
