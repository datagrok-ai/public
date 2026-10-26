/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/network-diagram/network-diagram-regressions.feature
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
import {addViewer, areaPainted, hasArea, hasNoArea, noErrors, painted, readingIs, repainted, reportsNoError, resizeTo, restoreSize, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Network diagram regressions of closed issues", () => {
  const session = feature(test, "features/viewers/network-diagram/network-diagram-regressions.feature", import.meta.url);
  test("Network diagram regressions of closed issues", {tag: ["@journey", "@viewers", "@realizes:viewers.network-diagram"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a network diagram viewer", () => addViewer(page, "network diagram"));
    await session.step(14, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(15, "And network diagram viewer should be painted", () => painted(page, el("network diagram viewer")));
    await run.scenario("The icon node shape draws every node and logs nothing per node (core#869)", async () => {
      await session.step(18, "When user sets \"nodeShape\" property of network diagram viewer to \"icon\"", () => setProperty(page, "nodeShape", el("network diagram viewer"), "icon"));
      await session.step(19, "Then network diagram viewer should have repainted", () => repainted(page, el("network diagram viewer")));
      await session.step(20, "And the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
      await session.step(21, "And the \"node names\" reading of network diagram viewer should contain \"F\"", () => readingContains(page, "node names", el("network diagram viewer"), "F"));
      await session.step(22, "And the \"node \\\"F\\\"\" area of network diagram viewer should be painted", () => areaPainted(page, "node \"F\"", el("network diagram viewer")));
      await session.step(23, "And the \"node \\\"true\\\"\" area of network diagram viewer should be painted", () => areaPainted(page, "node \"true\"", el("network diagram viewer")));
      await session.step(24, "And network diagram viewer should report no error", () => reportsNoError(page, el("network diagram viewer")));
      await session.step(25, "And no errors should have been logged", () => noErrors(page));
      await session.step(26, "When user sets \"nodeShape\" property of network diagram viewer to \"dot\"", () => setProperty(page, "nodeShape", el("network diagram viewer"), "dot"));
      await session.step(27, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    });
    await run.scenario("Auto Layout drops the column selectors on a small diagram, and Auto Layout off keeps them (public#2129)", async () => {
      await session.step(30, "Then network diagram viewer should have a \"node 1 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(31, "And network diagram viewer should have a \"node 2 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 2 selector"));
      await session.step(32, "When user resizes network diagram viewer to 180 by 180", () => resizeTo(page, el("network diagram viewer"), 180, 180));
      await session.step(33, "Then network diagram viewer should not have a \"node 1 selector\" area", () => hasNoArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(34, "And network diagram viewer should not have a \"node 2 selector\" area", () => hasNoArea(page, el("network diagram viewer"), "node 2 selector"));
      await session.step(35, "When user sets \"autoLayout\" property of network diagram viewer to \"false\"", () => setProperty(page, "autoLayout", el("network diagram viewer"), "false"));
      await session.step(36, "Then network diagram viewer should have a \"node 1 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(37, "And network diagram viewer should have a \"node 2 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 2 selector"));
      await session.step(38, "When user sets \"autoLayout\" property of network diagram viewer to \"true\"", () => setProperty(page, "autoLayout", el("network diagram viewer"), "true"));
      await session.step(39, "Then network diagram viewer should not have a \"node 1 selector\" area", () => hasNoArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(40, "When user restores the size of network diagram viewer", () => restoreSize(page, el("network diagram viewer")));
      await session.step(41, "Then network diagram viewer should have a \"node 1 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
