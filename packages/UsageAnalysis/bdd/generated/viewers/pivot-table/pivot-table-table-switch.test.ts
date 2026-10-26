/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-table-switch.feature
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
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, hasArea, hasNoArea, noErrors, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — the Table property switched to another table and back", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-table-switch.feature", import.meta.url);
  test("Switched to spgi-100 and back, the Data row names one table each time and the cross tab follows", {tag: ["@viewers", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(15, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
    await session.step(16, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(17, "Then the \"data\" reading of pivot table viewer should be \"demog-1000\"", () => readingReads(page, "data", el("pivot table viewer"), "demog-1000"));
    await session.step(18, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(21, "When user sets \"table\" property of pivot table viewer to \"spgi-100\"", () => setProperty(page, "table", el("pivot table viewer"), "spgi-100"));
    await session.step(22, "Then pivot table viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("pivot table viewer"), "spgi-100"));
    await session.step(23, "And the \"data\" reading of pivot table viewer should be \"spgi-100\"", () => readingReads(page, "data", el("pivot table viewer"), "spgi-100"));
    await session.step(24, "And the \"rows shown\" reading of pivot table viewer should be 100", () => readingIs(page, "rows shown", el("pivot table viewer"), 100));
    await session.step(25, "And the \"group by\" reading of pivot table viewer should be \"Cellular assay 1 Curve ID\"", () => readingReads(page, "group by", el("pivot table viewer"), "Cellular assay 1 Curve ID"));
    await session.step(26, "And the \"aggregated rows\" reading of pivot table viewer should be 87", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 87));
    await session.step(27, "And the \"error\" reading of pivot table viewer should be \"\"", () => readingReads(page, "error", el("pivot table viewer"), ""));
    await session.step(28, "When user sets \"table\" property of pivot table viewer to \"demog-1000\"", () => setProperty(page, "table", el("pivot table viewer"), "demog-1000"));
    await session.step(29, "Then pivot table viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("pivot table viewer"), "demog-1000"));
    await session.step(30, "And the \"data\" reading of pivot table viewer should be \"demog-1000\"", () => readingReads(page, "data", el("pivot table viewer"), "demog-1000"));
    await session.step(31, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
    await session.step(32, "And the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
    await session.step(33, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(34, "And pivot table viewer should have a \"data chip demog-1000\" area", () => hasArea(page, el("pivot table viewer"), "data chip demog-1000"));
    await session.step(35, "And pivot table viewer should not have a \"data chip spgi-100\" area", () => hasNoArea(page, el("pivot table viewer"), "data chip spgi-100"));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
  });
});
