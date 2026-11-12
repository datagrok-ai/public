/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/network-diagram/network-diagram.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.network-diagram, entities.viewer.action.close-viewer, core#869, public#2129]
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
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaPainted, areaWider, areasDiffer, areasSame, hasArea, hasNoArea, noErrors, painted, readingIs, readingReads, repainted, reportsNoError, resizeTo, restoreSize, setProperties, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains, readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Network diagram graph shape, node columns, filtering and chrome", () => {
  const session = feature(test, "features/viewers/network-diagram/network-diagram.feature", import.meta.url);
  test("Network diagram graph shape, node columns, filtering and chrome", {tag: ["@journey", "@viewers", "@realizes:viewers.network-diagram", "@realizes:entities.viewer.action.close-viewer", "@realizes:core#869", "@realizes:public#2129"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(33, "Given user is logged in", () => loggedIn(page));
    await session.step(34, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(35, "And user adds a network diagram viewer", () => addViewer(page, "network diagram"));
    await session.step(36, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(37, "And the \"node 1 column\" reading of network diagram viewer should be \"SEX\"", () => readingReads(page, "node 1 column", el("network diagram viewer"), "SEX"));
    await session.step(38, "And the \"node 2 column\" reading of network diagram viewer should be \"CONTROL\"", () => readingReads(page, "node 2 column", el("network diagram viewer"), "CONTROL"));
    await session.step(39, "And the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
    await session.step(40, "And the \"edges\" reading of network diagram viewer should be 4", () => readingIs(page, "edges", el("network diagram viewer"), 4));
    await session.step(41, "And the \"rows shown\" reading of network diagram viewer should be 1000", () => readingIs(page, "rows shown", el("network diagram viewer"), 1000));
    await session.step(42, "And network diagram viewer should be painted", () => painted(page, el("network diagram viewer")));
    await session.step(43, "And network diagram viewer should report no error", () => reportsNoError(page, el("network diagram viewer")));
    await run.scenario("Each node stands for the rows of the edges hanging off it", async () => {
      await session.step(46, "Then the \"node names\" reading of network diagram viewer should contain \"F\"", () => readingContains(page, "node names", el("network diagram viewer"), "F"));
      await session.step(47, "And the \"node names\" reading of network diagram viewer should contain \"M\"", () => readingContains(page, "node names", el("network diagram viewer"), "M"));
      await session.step(48, "And the \"rows of node \\\"F\\\"\" reading of network diagram viewer should be 553", () => readingIs(page, "rows of node \"F\"", el("network diagram viewer"), 553));
      await session.step(49, "And the \"rows of node \\\"M\\\"\" reading of network diagram viewer should be 447", () => readingIs(page, "rows of node \"M\"", el("network diagram viewer"), 447));
      await session.step(50, "And the \"rows of node \\\"false\\\"\" reading of network diagram viewer should be 994", () => readingIs(page, "rows of node \"false\"", el("network diagram viewer"), 994));
      await session.step(51, "And the \"rows of node \\\"true\\\"\" reading of network diagram viewer should be 6", () => readingIs(page, "rows of node \"true\"", el("network diagram viewer"), 6));
      await session.step(52, "And network diagram viewer should have a \"node \\\"F\\\"\" area", () => hasArea(page, el("network diagram viewer"), "node \"F\""));
      await session.step(53, "And network diagram viewer should have a \"node \\\"M\\\"\" area", () => hasArea(page, el("network diagram viewer"), "node \"M\""));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Changing Node 1 rebuilds the graph around the other column's values", async () => {
      await session.step(57, "When user sets \"node1ColumnName\" property of network diagram viewer to \"RACE\"", () => setProperty(page, "node1ColumnName", el("network diagram viewer"), "RACE"));
      await session.step(58, "Then the \"node 1 column\" reading of network diagram viewer should be \"RACE\"", () => readingReads(page, "node 1 column", el("network diagram viewer"), "RACE"));
      await session.step(59, "And the \"nodes\" reading of network diagram viewer should be 6", () => readingIs(page, "nodes", el("network diagram viewer"), 6));
      await session.step(60, "And the \"edges\" reading of network diagram viewer should be 5", () => readingIs(page, "edges", el("network diagram viewer"), 5));
      await session.step(61, "And the \"node names\" reading of network diagram viewer should contain \"Caucasian\"", () => readingContains(page, "node names", el("network diagram viewer"), "Caucasian"));
      await session.step(62, "And the \"node names\" reading of network diagram viewer should not contain \"F\"", () => readingNotContains(page, "node names", el("network diagram viewer"), "F"));
      await session.step(63, "And the \"rows of node \\\"Caucasian\\\"\" reading of network diagram viewer should be 896", () => readingIs(page, "rows of node \"Caucasian\"", el("network diagram viewer"), 896));
      await session.step(64, "And the \"rows of node \\\"Other\\\"\" reading of network diagram viewer should be 62", () => readingIs(page, "rows of node \"Other\"", el("network diagram viewer"), 62));
      await session.step(65, "And the \"rows of node \\\"Black\\\"\" reading of network diagram viewer should be 27", () => readingIs(page, "rows of node \"Black\"", el("network diagram viewer"), 27));
      await session.step(66, "And the \"rows of node \\\"Asian\\\"\" reading of network diagram viewer should be 15", () => readingIs(page, "rows of node \"Asian\"", el("network diagram viewer"), 15));
      await session.step(67, "And network diagram viewer should have a \"node \\\"Caucasian\\\"\" area", () => hasArea(page, el("network diagram viewer"), "node \"Caucasian\""));
      await session.step(68, "When user sets \"node1ColumnName\" property of network diagram viewer to \"SEX\"", () => setProperty(page, "node1ColumnName", el("network diagram viewer"), "SEX"));
      await session.step(69, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filtered-away node is not drawn, and Show Filtered Out Nodes brings it back", async () => {
      await session.step(73, "When user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
      await session.step(74, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(75, "And the \"rows shown\" reading of network diagram viewer should be 553", () => readingIs(page, "rows shown", el("network diagram viewer"), 553));
      await session.step(76, "And the \"nodes\" reading of network diagram viewer should be 3", () => readingIs(page, "nodes", el("network diagram viewer"), 3));
      await session.step(77, "And the \"edges\" reading of network diagram viewer should be 2", () => readingIs(page, "edges", el("network diagram viewer"), 2));
      await session.step(78, "And the \"node names\" reading of network diagram viewer should not contain \"M\"", () => readingNotContains(page, "node names", el("network diagram viewer"), "M"));
      await session.step(79, "And network diagram viewer should not have a \"node \\\"M\\\"\" area", () => hasNoArea(page, el("network diagram viewer"), "node \"M\""));
      await session.step(80, "And network diagram viewer should have a \"node \\\"F\\\"\" area", () => hasArea(page, el("network diagram viewer"), "node \"F\""));
      await session.step(81, "When user sets \"showFilteredOutNodes\" property of network diagram viewer to \"true\"", () => setProperty(page, "showFilteredOutNodes", el("network diagram viewer"), "true"));
      await session.step(82, "Then the \"filtered out nodes shown\" reading of network diagram viewer should be \"true\"", () => readingReads(page, "filtered out nodes shown", el("network diagram viewer"), "true"));
      await session.step(83, "And the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
      await session.step(84, "And the \"edges\" reading of network diagram viewer should be 4", () => readingIs(page, "edges", el("network diagram viewer"), 4));
      await session.step(85, "And network diagram viewer should have a \"node \\\"M\\\"\" area", () => hasArea(page, el("network diagram viewer"), "node \"M\""));
      await session.step(86, "And the \"rows shown\" reading of network diagram viewer should be 553", () => readingIs(page, "rows shown", el("network diagram viewer"), 553));
      await session.step(87, "When user sets \"showFilteredOutNodes\" property of network diagram viewer to \"false\"", () => setProperty(page, "showFilteredOutNodes", el("network diagram viewer"), "false"));
      await session.step(88, "Then the \"nodes\" reading of network diagram viewer should be 3", () => readingIs(page, "nodes", el("network diagram viewer"), 3));
      await session.step(89, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(90, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(91, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(92, "And the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
      await session.step(93, "And the \"rows shown\" reading of network diagram viewer should be 1000", () => readingIs(page, "rows shown", el("network diagram viewer"), 1000));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Column Selectors takes the two on-viewer selectors away", async () => {
      await session.step(97, "Then network diagram viewer should have a \"node 1 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(98, "And network diagram viewer should have a \"node 2 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 2 selector"));
      await session.step(99, "When user sets \"showColumnSelectors\" property of network diagram viewer to \"false\"", () => setProperty(page, "showColumnSelectors", el("network diagram viewer"), "false"));
      await session.step(100, "Then network diagram viewer should not have a \"node 1 selector\" area", () => hasNoArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(101, "And network diagram viewer should not have a \"node 2 selector\" area", () => hasNoArea(page, el("network diagram viewer"), "node 2 selector"));
      await session.step(102, "And the \"node 1 column\" reading of network diagram viewer should be \"SEX\"", () => readingReads(page, "node 1 column", el("network diagram viewer"), "SEX"));
      await session.step(103, "And network diagram viewer should have a \"node \\\"F\\\"\" area", () => hasArea(page, el("network diagram viewer"), "node \"F\""));
      await session.step(104, "When user sets \"showColumnSelectors\" property of network diagram viewer to \"true\"", () => setProperty(page, "showColumnSelectors", el("network diagram viewer"), "true"));
      await session.step(105, "Then network diagram viewer should have a \"node 1 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Arrows is a setting the viewer reports, and a rebuild keeps it", async () => {
      await session.step(109, "Then the \"arrows\" reading of network diagram viewer should be \"none\"", () => readingReads(page, "arrows", el("network diagram viewer"), "none"));
      await session.step(110, "When user sets \"showArrows\" property of network diagram viewer to \"to\"", () => setProperty(page, "showArrows", el("network diagram viewer"), "to"));
      await session.step(111, "Then the \"arrows\" reading of network diagram viewer should be \"to\"", () => readingReads(page, "arrows", el("network diagram viewer"), "to"));
      await session.step(112, "And network diagram viewer should have repainted", () => repainted(page, el("network diagram viewer")));
      await session.step(113, "When user sets \"node1ColumnName\" property of network diagram viewer to \"RACE\"", () => setProperty(page, "node1ColumnName", el("network diagram viewer"), "RACE"));
      await session.step(114, "Then the \"nodes\" reading of network diagram viewer should be 6", () => readingIs(page, "nodes", el("network diagram viewer"), 6));
      await session.step(115, "And the \"arrows\" reading of network diagram viewer should be \"to\"", () => readingReads(page, "arrows", el("network diagram viewer"), "to"));
      await session.step(116, "When user sets properties of network diagram viewer:", () => setProperties(page, el("network diagram viewer"), [["node1ColumnName","SEX"],["showArrows","none"]]), [["node1ColumnName","SEX"],["showArrows","none"]]);
      await session.step(119, "Then the \"arrows\" reading of network diagram viewer should be \"none\"", () => readingReads(page, "arrows", el("network diagram viewer"), "none"));
      await session.step(120, "And the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Node 1 Color and Size code the first column's nodes, and Edge Color recolors the edges", async () => {
      await session.step(124, "Then the \"node \\\"F\\\"\" and \"node \\\"M\\\"\" areas of network diagram viewer should be painted in the same colors", () => areasSame(page, "node \"F\"", "node \"M\"", el("network diagram viewer")));
      await session.step(125, "When user sets \"node1ColorColumnName\" property of network diagram viewer to \"SEX\"", () => setProperty(page, "node1ColorColumnName", el("network diagram viewer"), "SEX"));
      await session.step(126, "Then the \"node \\\"F\\\"\" and \"node \\\"M\\\"\" areas of network diagram viewer should be painted in different colors", () => areasDiffer(page, "node \"F\"", "node \"M\"", el("network diagram viewer")));
      await session.step(127, "And the \"node \\\"false\\\"\" and \"node \\\"true\\\"\" areas of network diagram viewer should be painted in the same colors", () => areasSame(page, "node \"false\"", "node \"true\"", el("network diagram viewer")));
      await session.step(128, "When user sets \"node1SizeColumnName\" property of network diagram viewer to \"AGE\"", () => setProperty(page, "node1SizeColumnName", el("network diagram viewer"), "AGE"));
      await session.step(129, "Then the \"node \\\"F\\\"\" area of network diagram viewer should be wider than before", () => areaWider(page, "node \"F\"", el("network diagram viewer")));
      await session.step(130, "And the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
      await session.step(131, "When user sets \"edgeColorColumnName\" property of network diagram viewer to \"AGE\"", () => setProperty(page, "edgeColorColumnName", el("network diagram viewer"), "AGE"));
      await session.step(132, "Then network diagram viewer should have repainted", () => repainted(page, el("network diagram viewer")));
      await session.step(133, "And the \"edges\" reading of network diagram viewer should be 4", () => readingIs(page, "edges", el("network diagram viewer"), 4));
      await session.step(134, "When user sets properties of network diagram viewer:", () => setProperties(page, el("network diagram viewer"), [["node1ColorColumnName",""],["node1SizeColumnName",""],["edgeColorColumnName",""]]), [["node1ColorColumnName",""],["node1SizeColumnName",""],["edgeColorColumnName",""]]);
      await session.step(138, "Then the \"node \\\"F\\\"\" and \"node \\\"M\\\"\" areas of network diagram viewer should be painted in the same colors", () => areasSame(page, "node \"F\"", "node \"M\"", el("network diagram viewer")));
      await session.step(139, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The icon node shape draws every node and logs nothing per node (core#869)", async () => {
      await session.step(143, "When user sets \"nodeShape\" property of network diagram viewer to \"icon\"", () => setProperty(page, "nodeShape", el("network diagram viewer"), "icon"));
      await session.step(144, "Then network diagram viewer should have repainted", () => repainted(page, el("network diagram viewer")));
      await session.step(145, "And the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
      await session.step(146, "And the \"node names\" reading of network diagram viewer should contain \"F\"", () => readingContains(page, "node names", el("network diagram viewer"), "F"));
      await session.step(147, "And the \"node \\\"F\\\"\" area of network diagram viewer should be painted", () => areaPainted(page, "node \"F\"", el("network diagram viewer")));
      await session.step(148, "And the \"node \\\"true\\\"\" area of network diagram viewer should be painted", () => areaPainted(page, "node \"true\"", el("network diagram viewer")));
      await session.step(149, "And network diagram viewer should report no error", () => reportsNoError(page, el("network diagram viewer")));
      await session.step(150, "When user sets \"nodeShape\" property of network diagram viewer to \"dot\"", () => setProperty(page, "nodeShape", el("network diagram viewer"), "dot"));
      await session.step(151, "Then the \"nodes\" reading of network diagram viewer should be 4", () => readingIs(page, "nodes", el("network diagram viewer"), 4));
      await session.step(152, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Auto Layout drops the column selectors on a small diagram, and Auto Layout off keeps them (public#2129)", async () => {
      await session.step(156, "Then network diagram viewer should have a \"node 1 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(157, "And network diagram viewer should have a \"node 2 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 2 selector"));
      await session.step(158, "When user resizes network diagram viewer to 180 by 180", () => resizeTo(page, el("network diagram viewer"), 180, 180));
      await session.step(159, "Then network diagram viewer should not have a \"node 1 selector\" area", () => hasNoArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(160, "And network diagram viewer should not have a \"node 2 selector\" area", () => hasNoArea(page, el("network diagram viewer"), "node 2 selector"));
      await session.step(161, "When user sets \"autoLayout\" property of network diagram viewer to \"false\"", () => setProperty(page, "autoLayout", el("network diagram viewer"), "false"));
      await session.step(162, "Then network diagram viewer should have a \"node 1 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(163, "And network diagram viewer should have a \"node 2 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 2 selector"));
      await session.step(164, "When user sets \"autoLayout\" property of network diagram viewer to \"true\"", () => setProperty(page, "autoLayout", el("network diagram viewer"), "true"));
      await session.step(165, "Then network diagram viewer should not have a \"node 1 selector\" area", () => hasNoArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(166, "When user restores the size of network diagram viewer", () => restoreSize(page, el("network diagram viewer")));
      await session.step(167, "Then network diagram viewer should have a \"node 1 selector\" area", () => hasArea(page, el("network diagram viewer"), "node 1 selector"));
      await session.step(168, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title bar closes the network diagram", async () => {
      await session.step(171, "When user clicks on close icon of network diagram viewer", () => clickOn(page, el("close icon of network diagram viewer")));
      await session.step(172, "Then network diagram viewer should be absent", () => shouldBe(page, el("network diagram viewer"), "absent"));
      await session.step(173, "And the open tableview should have 0 network diagram viewers", () => viewerCount(page, 0, "network diagram"));
      await session.step(174, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
