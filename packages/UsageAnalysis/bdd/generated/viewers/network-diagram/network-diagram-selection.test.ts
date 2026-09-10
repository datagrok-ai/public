/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/network-diagram/network-diagram-selection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.network-diagram]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, clearSelection, filterPasses, noneSelected, onlyOfSelected, selectedPassFilter, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, eventFired, hasArea, hasNoArea, listenFor, noErrors, readingIs, repainted, reportsNoError, setProperties} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Clicking a network diagram node selects the rows behind it", () => {
  const session = feature(test, "features/viewers/network-diagram/network-diagram-selection.feature", import.meta.url);
  test("Clicking a network diagram node selects the rows behind it", {tag: ["@journey", "@viewers", "@realizes:viewers.network-diagram"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a network diagram viewer", () => addViewer(page, "network diagram"));
    await session.step(21, "And user clears the row selection", () => clearSelection(page));
    await session.step(22, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(23, "And the \"rows of node \\\"F\\\"\" reading of network diagram viewer should be 553", () => readingIs(page, "rows of node \"F\"", el("network diagram viewer"), 553));
    await session.step(24, "And no rows should be selected", () => noneSelected(page));
    await session.step(25, "And network diagram viewer should report no error", () => reportsNoError(page, el("network diagram viewer")));
    await run.scenario("Clicking a node selects exactly the rows it stands for", async () => {
      await session.step(28, "When user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(29, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(30, "And only rows where \"SEX\" is \"F\" should be selected", () => onlyOfSelected(page, "SEX", "F"));
      await session.step(31, "And network diagram viewer should have repainted", () => repainted(page, el("network diagram viewer")));
      await session.step(32, "When user clears the row selection", () => clearSelection(page));
      await session.step(33, "And user clicks on the \"node \\\"M\\\"\" area of network diagram viewer", () => clickArea(page, "node \"M\"", el("network diagram viewer")));
      await session.step(34, "Then 447 rows should be selected", () => selectedRowCount(page, 447));
      await session.step(35, "And only rows where \"SEX\" is \"M\" should be selected", () => onlyOfSelected(page, "SEX", "M"));
      await session.step(36, "When user clears the row selection", () => clearSelection(page));
      await session.step(37, "Then no rows should be selected", () => noneSelected(page));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A node of the second column selects the rows of that column's value", async () => {
      await session.step(41, "When user clicks on the \"node \\\"true\\\"\" area of network diagram viewer", () => clickArea(page, "node \"true\"", el("network diagram viewer")));
      await session.step(42, "Then 6 rows should be selected", () => selectedRowCount(page, 6));
      await session.step(43, "And only rows where \"CONTROL\" is \"true\" should be selected", () => onlyOfSelected(page, "CONTROL", "true"));
      await session.step(44, "When user clears the row selection", () => clearSelection(page));
      await session.step(45, "And user clicks on the \"node \\\"false\\\"\" area of network diagram viewer", () => clickArea(page, "node \"false\"", el("network diagram viewer")));
      await session.step(46, "Then 994 rows should be selected", () => selectedRowCount(page, 994));
      await session.step(47, "And only rows where \"CONTROL\" is \"false\" should be selected", () => onlyOfSelected(page, "CONTROL", "false"));
      await session.step(48, "When user clears the row selection", () => clearSelection(page));
      await session.step(49, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A node click is announced as an event", async () => {
      await session.step(52, "Given user listens for \"d4-network-diagram-node-click\" event on network diagram viewer", () => listenFor(page, "d4-network-diagram-node-click", el("network diagram viewer")));
      await session.step(53, "When user clicks on the \"node \\\"M\\\"\" area of network diagram viewer", () => clickArea(page, "node \"M\"", el("network diagram viewer")));
      await session.step(54, "Then \"d4-network-diagram-node-click\" event should have fired on network diagram viewer", () => eventFired(page, "d4-network-diagram-node-click", el("network diagram viewer")));
      await session.step(55, "And 447 rows should be selected", () => selectedRowCount(page, 447));
      await session.step(56, "When user clears the row selection", () => clearSelection(page));
      await session.step(57, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Select Rows On Click off the same click selects nothing", async () => {
      await session.step(60, "When user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(61, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(62, "When user clears the row selection", () => clearSelection(page));
      await session.step(63, "And user sets properties of network diagram viewer:", () => setProperties(page, el("network diagram viewer"), [["selectRowsOnClick","false"],["selectEdgesOnClick","false"]]));
      await session.step(66, "And user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(67, "Then no rows should be selected", () => noneSelected(page));
      await session.step(68, "When user sets properties of network diagram viewer:", () => setProperties(page, el("network diagram viewer"), [["selectRowsOnClick","true"],["selectEdgesOnClick","true"]]));
      await session.step(71, "And user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(72, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(73, "When user clears the row selection", () => clearSelection(page));
      await session.step(74, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filtered-away node cannot be clicked because it is not there", async () => {
      await session.step(77, "When user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
      await session.step(78, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(79, "And network diagram viewer should not have a \"node \\\"M\\\"\" area", () => hasNoArea(page, el("network diagram viewer"), "node \"M\""));
      await session.step(80, "When user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(81, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(82, "And every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(83, "When user clears the row selection", () => clearSelection(page));
      await session.step(84, "And user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(85, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(86, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(87, "And network diagram viewer should have a \"node \\\"M\\\"\" area", () => hasArea(page, el("network diagram viewer"), "node \"M\""));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
