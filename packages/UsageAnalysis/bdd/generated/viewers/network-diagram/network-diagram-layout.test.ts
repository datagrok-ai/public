/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/network-diagram/network-diagram-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.network-diagram]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearSelection, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, hasArea, noErrors, pickFromContextMenu, propertyShouldBe, readingAsRemembered, readingIs, readingNotAsRemembered, rememberReading, repainted, reportsNoError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Suspend Simulation freezes the layout, and the Layout menu says so too", () => {
  const session = feature(test, "features/viewers/network-diagram/network-diagram-layout.feature", import.meta.url);
  test("Suspend Simulation freezes the layout, and the Layout menu says so too", {tag: ["@journey", "@viewers", "@realizes:viewers.network-diagram"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a network diagram viewer", () => addViewer(page, "network diagram"));
    await session.step(19, "And user clears the row selection", () => clearSelection(page));
    await session.step(20, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(21, "And \"suspendSimulation\" property of network diagram viewer should be \"false\"", () => propertyShouldBe(page, "suspendSimulation", el("network diagram viewer"), "false"));
    await session.step(22, "And network diagram viewer should report no error", () => reportsNoError(page, el("network diagram viewer")));
    await run.scenario("With the simulation suspended the layout does not move", async () => {
      await session.step(25, "When user sets \"suspendSimulation\" property of network diagram viewer to \"true\"", () => setProperty(page, "suspendSimulation", el("network diagram viewer"), "true"));
      await session.step(26, "And user remembers the \"layout signature\" reading of network diagram viewer", () => rememberReading(page, "layout signature", el("network diagram viewer")));
      await session.step(27, "And user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(28, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(29, "And network diagram viewer should have repainted", () => repainted(page, el("network diagram viewer")));
      await session.step(30, "And the \"layout signature\" reading of network diagram viewer should be as remembered", () => readingAsRemembered(page, "layout signature", el("network diagram viewer")));
      await session.step(31, "When user clears the row selection", () => clearSelection(page));
      await session.step(32, "Then the \"layout signature\" reading of network diagram viewer should be as remembered", () => readingAsRemembered(page, "layout signature", el("network diagram viewer")));
      await session.step(33, "When user sets \"suspendSimulation\" property of network diagram viewer to \"false\"", () => setProperty(page, "suspendSimulation", el("network diagram viewer"), "false"));
      await session.step(34, "Then \"suspendSimulation\" property of network diagram viewer should be \"false\"", () => propertyShouldBe(page, "suspendSimulation", el("network diagram viewer"), "false"));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A rebuild of the graph moves the layout, so an unchanged signature means something", async () => {
      await session.step(38, "When user remembers the \"layout signature\" reading of network diagram viewer", () => rememberReading(page, "layout signature", el("network diagram viewer")));
      await session.step(39, "And user sets \"node1ColumnName\" property of network diagram viewer to \"RACE\"", () => setProperty(page, "node1ColumnName", el("network diagram viewer"), "RACE"));
      await session.step(40, "Then the \"nodes\" reading of network diagram viewer should be 6", () => readingIs(page, "nodes", el("network diagram viewer"), 6));
      await session.step(41, "And the \"layout signature\" reading of network diagram viewer should not be as remembered", () => readingNotAsRemembered(page, "layout signature", el("network diagram viewer")));
      await session.step(42, "When user sets \"node1ColumnName\" property of network diagram viewer to \"SEX\"", () => setProperty(page, "node1ColumnName", el("network diagram viewer"), "SEX"));
      await session.step(43, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Layout > Suspend simulation is the same setting, reached from the viewer's own menu", async () => {
      await session.step(47, "Then \"suspendSimulation\" property of network diagram viewer should be \"false\"", () => propertyShouldBe(page, "suspendSimulation", el("network diagram viewer"), "false"));
      await session.step(48, "When user picks \"Layout > Suspend simulation\" from the context menu of network diagram viewer", () => pickFromContextMenu(page, "Layout > Suspend simulation", el("network diagram viewer")));
      await session.step(49, "Then \"suspendSimulation\" property of network diagram viewer should be \"true\"", () => propertyShouldBe(page, "suspendSimulation", el("network diagram viewer"), "true"));
      await session.step(50, "When user remembers the \"layout signature\" reading of network diagram viewer", () => rememberReading(page, "layout signature", el("network diagram viewer")));
      await session.step(51, "And user clicks on the \"node \\\"M\\\"\" area of network diagram viewer", () => clickArea(page, "node \"M\"", el("network diagram viewer")));
      await session.step(52, "Then 447 rows should be selected", () => selectedRowCount(page, 447));
      await session.step(53, "And the \"layout signature\" reading of network diagram viewer should be as remembered", () => readingAsRemembered(page, "layout signature", el("network diagram viewer")));
      await session.step(54, "When user clears the row selection", () => clearSelection(page));
      await session.step(55, "And user picks \"Layout > Suspend simulation\" from the context menu of network diagram viewer", () => pickFromContextMenu(page, "Layout > Suspend simulation", el("network diagram viewer")));
      await session.step(56, "Then \"suspendSimulation\" property of network diagram viewer should be \"false\"", () => propertyShouldBe(page, "suspendSimulation", el("network diagram viewer"), "false"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Reset View leaves the graph and its counts where they were", async () => {
      await session.step(60, "When user remembers the \"layout signature\" reading of network diagram viewer", () => rememberReading(page, "layout signature", el("network diagram viewer")));
      await session.step(61, "And user picks \"Reset View\" from the context menu of network diagram viewer", () => pickFromContextMenu(page, "Reset View", el("network diagram viewer")));
      await session.step(62, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
      await session.step(63, "And the \"edges\" reading of network diagram viewer should be 4", () => readingIs(page, "edges", el("network diagram viewer"), 4));
      await session.step(64, "And the \"rows of node \\\"F\\\"\" reading of network diagram viewer should be 553", () => readingIs(page, "rows of node \"F\"", el("network diagram viewer"), 553));
      await session.step(65, "And the \"layout signature\" reading of network diagram viewer should be as remembered", () => readingAsRemembered(page, "layout signature", el("network diagram viewer")));
      await session.step(66, "And network diagram viewer should have a \"node \\\"F\\\"\" area", () => hasArea(page, el("network diagram viewer"), "node \"F\""));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
