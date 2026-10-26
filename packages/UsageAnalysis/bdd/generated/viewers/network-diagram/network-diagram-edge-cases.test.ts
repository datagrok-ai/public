/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/network-diagram/network-diagram-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.network-diagram, core#1041]
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
import {deleteSelected, filterBetween, filterPasses, filterTo, resetFilter, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, readingIs, readingReads, reportsNoError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Network diagram edge cases — an empty row set, every row deleted, one column for both nodes", () => {
  const session = feature(test, "features/viewers/network-diagram/network-diagram-edge-cases.feature", import.meta.url);
  test("With every row filtered out the diagram hides every node and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.network-diagram"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a network diagram viewer", () => addViewer(page, "network diagram"));
    await session.step(14, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(15, "And the \"edges\" reading of network diagram viewer should be 4", () => readingIs(page, "edges", el("network diagram viewer"), 4));
    await session.step(18, "When user filters rows where \"AGE\" is between 200 and 300", () => filterBetween(page, "AGE", 200, 300));
    await session.step(19, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(20, "And the \"rows shown\" reading of network diagram viewer should be 0", () => readingIs(page, "rows shown", el("network diagram viewer"), 0));
    await session.step(21, "And the \"nodes\" reading of network diagram viewer should be 0", () => readingIs(page, "nodes", el("network diagram viewer"), 0));
    await session.step(22, "And no errors should have been logged", () => noErrors(page));
    await session.step(23, "When user resets the filter", () => resetFilter(page));
    await session.step(24, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
  });
  test("Once every row of the table is deleted the diagram logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.network-diagram"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a network diagram viewer", () => addViewer(page, "network diagram"));
    await session.step(14, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(15, "And the \"edges\" reading of network diagram viewer should be 4", () => readingIs(page, "edges", el("network diagram viewer"), 4));
    await session.step(27, "When user selects all rows", () => selectAllRows(page));
    await session.step(28, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(29, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(30, "And the \"rows shown\" reading of network diagram viewer should be 0", () => readingIs(page, "rows shown", el("network diagram viewer"), 0));
    await session.step(31, "And network diagram viewer should report no error", () => reportsNoError(page, el("network diagram viewer")));
    await session.step(32, "And no errors should have been logged", () => noErrors(page));
  });
  test("With the same column as Node 1 and Node 2 a filter change raises no error (core#1041)", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.network-diagram", "@known-failure:41", "@realizes:core#1041"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a network diagram viewer", () => addViewer(page, "network diagram"));
    await session.step(14, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(15, "And the \"edges\" reading of network diagram viewer should be 4", () => readingIs(page, "edges", el("network diagram viewer"), 4));
    await knownFailure(async () => {
      await session.step(38, "When user sets \"node2ColumnName\" property of network diagram viewer to \"SEX\"", () => setProperty(page, "node2ColumnName", el("network diagram viewer"), "SEX"));
      await session.step(39, "Then the \"node 2 column\" reading of network diagram viewer should be \"SEX\"", () => readingReads(page, "node 2 column", el("network diagram viewer"), "SEX"));
      await session.step(40, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(41, "Then no errors should have been logged", () => noErrors(page));
    }, 41);
  });
});
