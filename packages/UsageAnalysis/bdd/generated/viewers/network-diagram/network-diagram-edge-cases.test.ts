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
import {shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueContains} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, filterTo, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, noErrors, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Network diagram edge cases — one column for both nodes, more than 10000 nodes", () => {
  const session = feature(test, "features/viewers/network-diagram/network-diagram-edge-cases.feature", import.meta.url);
  test("With the same column as Node 1 and Node 2 a filter change raises no error (core#1041)", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.network-diagram", "@realizes:core#1041"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a network diagram viewer", () => addViewer(page, "network diagram"));
    await session.step(12, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(13, "And the \"edges\" reading of network diagram viewer should be 4", () => readingIs(page, "edges", el("network diagram viewer"), 4));
    await session.step(19, "When user sets \"node2ColumnName\" property of network diagram viewer to \"SEX\"", () => setProperty(page, "node2ColumnName", el("network diagram viewer"), "SEX"));
    await session.step(20, "Then the \"node 2 column\" reading of network diagram viewer should be \"SEX\"", () => readingReads(page, "node 2 column", el("network diagram viewer"), "SEX"));
    await session.step(21, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
    await session.step(22, "Then no errors should have been logged", () => noErrors(page));
  });
  test("More than 10000 nodes are not drawn and the diagram says so", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.network-diagram"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a network diagram viewer", () => addViewer(page, "network diagram"));
    await session.step(12, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(13, "And the \"edges\" reading of network diagram viewer should be 4", () => readingIs(page, "edges", el("network diagram viewer"), 4));
    await session.step(28, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(29, "Then the table should have 5850 rows", () => rowCount(page, 5850));
    await session.step(30, "When user adds a calculated column \"PARTNER\" with formula \"${USUBJID} + '-b'\"", () => addCalculated(page, "PARTNER", "${USUBJID} + '-b'"));
    await session.step(31, "Then every value of \"PARTNER\" column should contain \"-b\"", () => everyValueContains(page, "PARTNER", "-b"));
    await session.step(32, "When user adds a network diagram viewer with:", () => addViewerWith(page, "network diagram", [["node1ColumnName","USUBJID"],["node2ColumnName","PARTNER"]]), [["node1ColumnName","USUBJID"],["node2ColumnName","PARTNER"]]);
    await session.step(35, "Then network diagram viewer should contain the text \"Too many nodes to display, max 10000\"", () => shouldContainText(page, el("network diagram viewer"), "Too many nodes to display, max 10000"));
    await session.step(36, "And the \"nodes\" reading of network diagram viewer should be 0", () => readingIs(page, "nodes", el("network diagram viewer"), 0));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
});
