/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/network-diagram/network-diagram-column-selector.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.network-diagram]
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, noErrors, propertyShouldBe, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInAreaSelector, readingContains, readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Network diagram — Node 1 picked on the viewer's own column selector", () => {
  const session = feature(test, "features/viewers/network-diagram/network-diagram-column-selector.feature", import.meta.url);
  test("Picking RACE in the Node 1 selector rebuilds the graph around the races, and SEX puts it back", {tag: ["@viewers", "@realizes:viewers.network-diagram"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a network diagram viewer", () => addViewer(page, "network diagram"));
    await session.step(12, "Then the \"node 1 column\" reading of network diagram viewer should be \"SEX\"", () => readingReads(page, "node 1 column", el("network diagram viewer"), "SEX"));
    await session.step(13, "And the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(14, "And network diagram viewer should have a \"node 1 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 1 selector"));
    await session.step(15, "When user picks \"RACE\" in the column selector at the \"node 1 selector\" area of network diagram viewer", () => pickInAreaSelector(page, "RACE", "node 1 selector", el("network diagram viewer")));
    await session.step(16, "Then the \"node 1 column\" reading of network diagram viewer should be \"RACE\"", () => readingReads(page, "node 1 column", el("network diagram viewer"), "RACE"));
    await session.step(17, "And \"node1ColumnName\" property of network diagram viewer should be \"RACE\"", () => propertyShouldBe(page, "node1ColumnName", el("network diagram viewer"), "RACE"));
    await session.step(18, "And the \"node 2 column\" reading of network diagram viewer should be \"CONTROL\"", () => readingReads(page, "node 2 column", el("network diagram viewer"), "CONTROL"));
    await session.step(19, "And the \"nodes\" reading of network diagram viewer should be 6", () => readingIs(page, "nodes", el("network diagram viewer"), 6));
    await session.step(20, "And the \"edges\" reading of network diagram viewer should be 5", () => readingIs(page, "edges", el("network diagram viewer"), 5));
    await session.step(21, "And the \"node names\" reading of network diagram viewer should contain \"Caucasian\"", () => readingContains(page, "node names", el("network diagram viewer"), "Caucasian"));
    await session.step(22, "And the \"node names\" reading of network diagram viewer should not contain \"F\"", () => readingNotContains(page, "node names", el("network diagram viewer"), "F"));
    await session.step(23, "And the \"rows of node \\\"Caucasian\\\"\" reading of network diagram viewer should be 896", () => readingIs(page, "rows of node \"Caucasian\"", el("network diagram viewer"), 896));
    await session.step(24, "And network diagram viewer should have a \"node \\\"Caucasian\\\"\" area", () => hasArea(page, el("network diagram viewer"), "node \"Caucasian\""));
    await session.step(25, "When user picks \"SEX\" in the column selector at the \"node 1 selector\" area of network diagram viewer", () => pickInAreaSelector(page, "SEX", "node 1 selector", el("network diagram viewer")));
    await session.step(26, "Then the \"node 1 column\" reading of network diagram viewer should be \"SEX\"", () => readingReads(page, "node 1 column", el("network diagram viewer"), "SEX"));
    await session.step(27, "And the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(28, "And the \"rows of node \\\"F\\\"\" reading of network diagram viewer should be 553", () => readingIs(page, "rows of node \"F\"", el("network diagram viewer"), 553));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
  });
});
