/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/network-diagram/network-diagram-selection.feature
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
import {clickOn, hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, allOfSelected, clearSelection, filterPasses, noneOfSelected, noneSelected, onlyOfSelected, selectedPassFilter, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, clickAreaHolding, eventFired, hasArea, hasNoArea, listenFor, noErrors, readingIs, repainted, reportsNoError, setProperties} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Clicking a network diagram node selects the rows behind it", () => {
  const session = feature(test, "features/viewers/network-diagram/network-diagram-selection.feature", import.meta.url);
  test("Clicking a network diagram node selects the rows behind it", {tag: ["@journey", "@viewers", "@realizes:viewers.network-diagram"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(27, "And user adds a network diagram viewer", () => addViewer(page, "network diagram"));
    await session.step(28, "And user clears the row selection", () => clearSelection(page));
    await session.step(29, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(30, "And the \"rows of node \\\"F\\\"\" reading of network diagram viewer should be 553", () => readingIs(page, "rows of node \"F\"", el("network diagram viewer"), 553));
    await session.step(31, "And no rows should be selected", () => noneSelected(page));
    await session.step(32, "And network diagram viewer should report no error", () => reportsNoError(page, el("network diagram viewer")));
    await run.scenario("Clicking a node selects exactly the rows it stands for", async () => {
      await session.step(35, "When user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(36, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(37, "And only rows where \"SEX\" is \"F\" should be selected", () => onlyOfSelected(page, "SEX", "F"));
      await session.step(38, "And network diagram viewer should have repainted", () => repainted(page, el("network diagram viewer")));
      await session.step(39, "When user clears the row selection", () => clearSelection(page));
      await session.step(40, "And user clicks on the \"node \\\"M\\\"\" area of network diagram viewer", () => clickArea(page, "node \"M\"", el("network diagram viewer")));
      await session.step(41, "Then 447 rows should be selected", () => selectedRowCount(page, 447));
      await session.step(42, "And only rows where \"SEX\" is \"M\" should be selected", () => onlyOfSelected(page, "SEX", "M"));
      await session.step(43, "When user clears the row selection", () => clearSelection(page));
      await session.step(44, "Then no rows should be selected", () => noneSelected(page));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A node of the second column selects the rows of that column's value", async () => {
      await session.step(48, "When user clicks on the \"node \\\"true\\\"\" area of network diagram viewer", () => clickArea(page, "node \"true\"", el("network diagram viewer")));
      await session.step(49, "Then 6 rows should be selected", () => selectedRowCount(page, 6));
      await session.step(50, "And only rows where \"CONTROL\" is \"true\" should be selected", () => onlyOfSelected(page, "CONTROL", "true"));
      await session.step(51, "When user clears the row selection", () => clearSelection(page));
      await session.step(52, "And user clicks on the \"node \\\"false\\\"\" area of network diagram viewer", () => clickArea(page, "node \"false\"", el("network diagram viewer")));
      await session.step(53, "Then 994 rows should be selected", () => selectedRowCount(page, 994));
      await session.step(54, "And only rows where \"CONTROL\" is \"false\" should be selected", () => onlyOfSelected(page, "CONTROL", "false"));
      await session.step(55, "When user clears the row selection", () => clearSelection(page));
      await session.step(56, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shift-click adds a node's rows, Control-click toggles them, Control+Shift-click takes them off", async () => {
      await session.step(59, "When user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(60, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(61, "When user clicks on the \"node \\\"true\\\"\" area of network diagram viewer holding Shift", () => clickAreaHolding(page, "node \"true\"", el("network diagram viewer"), "Shift"));
      await session.step(62, "Then 554 rows should be selected", () => selectedRowCount(page, 554));
      await session.step(63, "And all rows where \"CONTROL\" is \"true\" should be selected", () => allOfSelected(page, "CONTROL", "true"));
      await session.step(64, "When user clicks on the \"node \\\"F\\\"\" area of network diagram viewer holding Control", () => clickAreaHolding(page, "node \"F\"", el("network diagram viewer"), "Control"));
      await session.step(65, "Then 1 row should be selected", () => selectedRowCount(page, 1));
      await session.step(66, "And no rows where \"SEX\" is \"F\" should be selected", () => noneOfSelected(page, "SEX", "F"));
      await session.step(67, "When user clicks on the \"node \\\"M\\\"\" area of network diagram viewer holding Shift", () => clickAreaHolding(page, "node \"M\"", el("network diagram viewer"), "Shift"));
      await session.step(68, "Then 447 rows should be selected", () => selectedRowCount(page, 447));
      await session.step(69, "When user clicks on the \"node \\\"true\\\"\" area of network diagram viewer holding Control+Shift", () => clickAreaHolding(page, "node \"true\"", el("network diagram viewer"), "Control+Shift"));
      await session.step(70, "Then 446 rows should be selected", () => selectedRowCount(page, 446));
      await session.step(71, "And no rows where \"CONTROL\" is \"true\" should be selected", () => noneOfSelected(page, "CONTROL", "true"));
      await session.step(72, "When user clears the row selection", () => clearSelection(page));
      await session.step(73, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A node click is announced as an event", async () => {
      await session.step(76, "Given user listens for \"d4-network-diagram-node-click\" event on network diagram viewer", () => listenFor(page, "d4-network-diagram-node-click", el("network diagram viewer")));
      await session.step(77, "When user clicks on the \"node \\\"M\\\"\" area of network diagram viewer", () => clickArea(page, "node \"M\"", el("network diagram viewer")));
      await session.step(78, "Then \"d4-network-diagram-node-click\" event should have fired on network diagram viewer", () => eventFired(page, "d4-network-diagram-node-click", el("network diagram viewer")));
      await session.step(79, "And 447 rows should be selected", () => selectedRowCount(page, 447));
      await session.step(80, "When user clears the row selection", () => clearSelection(page));
      await session.step(81, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Select Rows On Click off the same click selects nothing", async () => {
      await session.step(84, "When user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(85, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(86, "When user clears the row selection", () => clearSelection(page));
      await session.step(87, "And user sets properties of network diagram viewer:", () => setProperties(page, el("network diagram viewer"), [["selectRowsOnClick","false"],["selectEdgesOnClick","false"]]), [["selectRowsOnClick","false"],["selectEdgesOnClick","false"]]);
      await session.step(90, "And user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(91, "Then no rows should be selected", () => noneSelected(page));
      await session.step(92, "When user sets properties of network diagram viewer:", () => setProperties(page, el("network diagram viewer"), [["selectRowsOnClick","true"],["selectEdgesOnClick","true"]]), [["selectRowsOnClick","true"],["selectEdgesOnClick","true"]]);
      await session.step(95, "And user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(96, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(97, "When user clears the row selection", () => clearSelection(page));
      await session.step(98, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filtered-away node cannot be clicked because it is not there", async () => {
      await session.step(101, "When user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
      await session.step(102, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(103, "And network diagram viewer should not have a \"node \\\"M\\\"\" area", () => hasNoArea(page, el("network diagram viewer"), "node \"M\""));
      await session.step(104, "When user clicks on the \"node \\\"F\\\"\" area of network diagram viewer", () => clickArea(page, "node \"F\"", el("network diagram viewer")));
      await session.step(105, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(106, "And every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(107, "When user clears the row selection", () => clearSelection(page));
      await session.step(108, "And user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(109, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(110, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(111, "And network diagram viewer should have a \"node \\\"M\\\"\" area", () => hasArea(page, el("network diagram viewer"), "node \"M\""));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
