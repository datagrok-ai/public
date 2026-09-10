/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/map/map.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.map-viewer]
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
import {clickOn, hoverOver, pressKey, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, clearSelection, filterPasses, noneSelected, selectedPassFilter, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, dragBoxHolding, hasArea, hasNoArea, hoverArea, noErrors, oneTooltip, pointerAway, propertyShouldBe, readingAsRemembered, readingHigher, readingIs, readingLower, readingNotAsRemembered, readingReads, readingsEqual, rememberReading, reportsNoError, setProperties, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {moveAcrossPlace, rememberPlace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Map viewer layers, zoom, selection and the point tooltip", () => {
  const session = feature(test, "features/viewers/map/map.feature", import.meta.url);
  test("Map viewer layers, zoom, selection and the point tooltip", {tag: ["@journey", "@viewers", "@realizes:viewers.map-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And user opens earthquakes dataset", () => openDataset(page, ds("earthquakes")));
    await session.step(32, "And user adds a map viewer", () => addViewer(page, "map"));
    await session.step(33, "Then map viewer should be visible", () => shouldBe(page, el("map viewer"), "visible"));
    await session.step(34, "And 2426 rows should pass the filter", () => filterPasses(page, 2426));
    await session.step(35, "And the \"markers\" reading of map viewer should be 2426", () => readingIs(page, "markers", el("map viewer"), 2426));
    await session.step(36, "And the \"rows shown\" reading of map viewer should be 2426", () => readingIs(page, "rows shown", el("map viewer"), 2426));
    await run.scenario("The geo columns are detected and every filtered row becomes a marker", async () => {
      await session.step(39, "Then \"latitudeColumnName\" property of map viewer should be \"Latitude\"", () => propertyShouldBe(page, "latitudeColumnName", el("map viewer"), "Latitude"));
      await session.step(40, "And \"longitudeColumnName\" property of map viewer should be \"Longitude\"", () => propertyShouldBe(page, "longitudeColumnName", el("map viewer"), "Longitude"));
      await session.step(41, "And the \"markers\" and \"rows shown\" readings of map viewer should be the same", () => readingsEqual(page, "markers", "rows shown", el("map viewer")));
      await session.step(42, "And the \"render type\" reading of map viewer should be \"markers\"", () => readingReads(page, "render type", el("map viewer"), "markers"));
      await session.step(43, "And the \"layer \\\"Markers GL\\\" visible\" reading of map viewer should be \"true\"", () => readingReads(page, "layer \"Markers GL\" visible", el("map viewer"), "true"));
      await session.step(44, "And the \"layer \\\"Heatmap\\\" visible\" reading of map viewer should be \"false\"", () => readingReads(page, "layer \"Heatmap\" visible", el("map viewer"), "false"));
      await session.step(45, "And the \"layers\" reading of map viewer should be 5", () => readingIs(page, "layers", el("map viewer"), 5));
      await session.step(46, "And map viewer should report no error", () => reportsNoError(page, el("map viewer")));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A colour column and a size column keep every feature on the map", async () => {
      await session.step(50, "When user sets properties of map viewer:", () => setProperties(page, el("map viewer"), [["colorColumnName","Magnitude"],["sizeColumnName","Depth"]]));
      await session.step(53, "Then \"colorColumnName\" property of map viewer should be \"Magnitude\"", () => propertyShouldBe(page, "colorColumnName", el("map viewer"), "Magnitude"));
      await session.step(54, "And \"sizeColumnName\" property of map viewer should be \"Depth\"", () => propertyShouldBe(page, "sizeColumnName", el("map viewer"), "Depth"));
      await session.step(55, "And the \"markers\" reading of map viewer should be 2426", () => readingIs(page, "markers", el("map viewer"), 2426));
      await session.step(56, "And the \"rows shown\" reading of map viewer should be 2426", () => readingIs(page, "rows shown", el("map viewer"), 2426));
      await session.step(57, "And the \"layer \\\"Markers GL\\\" visible\" reading of map viewer should be \"true\"", () => readingReads(page, "layer \"Markers GL\" visible", el("map viewer"), "true"));
      await session.step(58, "When user sets properties of map viewer:", () => setProperties(page, el("map viewer"), [["colorColumnName",""],["sizeColumnName",""]]));
      await session.step(61, "Then the \"markers\" reading of map viewer should be 2426", () => readingIs(page, "markers", el("map viewer"), 2426));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The layers button reveals a row per layer, and closes them again", async () => {
      await session.step(65, "Then map viewer should not have a \"layers panel\" area", () => hasNoArea(page, el("map viewer"), "layers panel"));
      await session.step(66, "And map viewer should not have a \"layer \\\"Heatmap\\\"\" area", () => hasNoArea(page, el("map viewer"), "layer \"Heatmap\""));
      await session.step(67, "When user clicks on \"Map layers\" button", () => clickOn(page, el("\"Map layers\" button")));
      await session.step(68, "Then map viewer should have a \"layers panel\" area", () => hasArea(page, el("map viewer"), "layers panel"));
      await session.step(69, "And map viewer should have a \"layer \\\"Bing sat\\\"\" area", () => hasArea(page, el("map viewer"), "layer \"Bing sat\""));
      await session.step(70, "And map viewer should have a \"layer \\\"BaseLayer\\\"\" area", () => hasArea(page, el("map viewer"), "layer \"BaseLayer\""));
      await session.step(71, "And map viewer should have a \"layer \\\"Heatmap\\\"\" area", () => hasArea(page, el("map viewer"), "layer \"Heatmap\""));
      await session.step(72, "And map viewer should have a \"layer \\\"Markers GL\\\"\" area", () => hasArea(page, el("map viewer"), "layer \"Markers GL\""));
      await session.step(73, "And map viewer should have a \"layer \\\"Heatmap\\\" visibility\" area", () => hasArea(page, el("map viewer"), "layer \"Heatmap\" visibility"));
      await session.step(74, "When user clicks on \"Map layers\" button", () => clickOn(page, el("\"Map layers\" button")));
      await session.step(75, "Then map viewer should not have a \"layers panel\" area", () => hasNoArea(page, el("map viewer"), "layers panel"));
      await session.step(76, "And map viewer should not have a \"layer \\\"Heatmap\\\"\" area", () => hasNoArea(page, el("map viewer"), "layer \"Heatmap\""));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The checkbox cell of a row toggles that layer, the row itself makes it current", async () => {
      await session.step(80, "When user clicks on \"Map layers\" button", () => clickOn(page, el("\"Map layers\" button")));
      await session.step(81, "Then the \"layer \\\"Heatmap\\\" visible\" reading of map viewer should be \"false\"", () => readingReads(page, "layer \"Heatmap\" visible", el("map viewer"), "false"));
      await session.step(82, "When user clicks on the \"layer \\\"Heatmap\\\" visibility\" area of map viewer", () => clickArea(page, "layer \"Heatmap\" visibility", el("map viewer")));
      await session.step(83, "Then the \"layer \\\"Heatmap\\\" visible\" reading of map viewer should be \"true\"", () => readingReads(page, "layer \"Heatmap\" visible", el("map viewer"), "true"));
      await session.step(84, "And the \"visible layers\" reading of map viewer should be 5", () => readingIs(page, "visible layers", el("map viewer"), 5));
      await session.step(85, "When user clicks on the \"layer \\\"Markers GL\\\" visibility\" area of map viewer", () => clickArea(page, "layer \"Markers GL\" visibility", el("map viewer")));
      await session.step(86, "Then the \"layer \\\"Markers GL\\\" visible\" reading of map viewer should be \"false\"", () => readingReads(page, "layer \"Markers GL\" visible", el("map viewer"), "false"));
      await session.step(87, "And the \"visible layers\" reading of map viewer should be 4", () => readingIs(page, "visible layers", el("map viewer"), 4));
      await session.step(88, "When user clicks on the \"layer \\\"BaseLayer\\\"\" area of map viewer", () => clickArea(page, "layer \"BaseLayer\"", el("map viewer")));
      await session.step(89, "Then \"currentLayer\" property of map viewer should be \"BaseLayer\"", () => propertyShouldBe(page, "currentLayer", el("map viewer"), "BaseLayer"));
      await session.step(90, "And the \"layer \\\"BaseLayer\\\" visible\" reading of map viewer should be \"true\"", () => readingReads(page, "layer \"BaseLayer\" visible", el("map viewer"), "true"));
      await session.step(91, "When user clicks on the \"layer \\\"Markers GL\\\" visibility\" area of map viewer", () => clickArea(page, "layer \"Markers GL\" visibility", el("map viewer")));
      await session.step(92, "And user clicks on the \"layer \\\"Heatmap\\\" visibility\" area of map viewer", () => clickArea(page, "layer \"Heatmap\" visibility", el("map viewer")));
      await session.step(93, "Then the \"layer \\\"Markers GL\\\" visible\" reading of map viewer should be \"true\"", () => readingReads(page, "layer \"Markers GL\" visible", el("map viewer"), "true"));
      await session.step(94, "And the \"layer \\\"Heatmap\\\" visible\" reading of map viewer should be \"false\"", () => readingReads(page, "layer \"Heatmap\" visible", el("map viewer"), "false"));
      await session.step(95, "When user clicks on \"Map layers\" button", () => clickOn(page, el("\"Map layers\" button")));
      await session.step(96, "Then map viewer should not have a \"layers panel\" area", () => hasNoArea(page, el("map viewer"), "layers panel"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Render Type moves the two product layers and leaves the base maps alone", async () => {
      await session.step(100, "When user sets \"renderType\" property of map viewer to \"heatmap\"", () => setProperty(page, "renderType", el("map viewer"), "heatmap"));
      await session.step(101, "Then the \"render type\" reading of map viewer should be \"heatmap\"", () => readingReads(page, "render type", el("map viewer"), "heatmap"));
      await session.step(102, "And the \"layer \\\"Heatmap\\\" visible\" reading of map viewer should be \"true\"", () => readingReads(page, "layer \"Heatmap\" visible", el("map viewer"), "true"));
      await session.step(103, "And the \"layer \\\"Markers GL\\\" visible\" reading of map viewer should be \"false\"", () => readingReads(page, "layer \"Markers GL\" visible", el("map viewer"), "false"));
      await session.step(104, "And the \"layer \\\"BaseLayer\\\" visible\" reading of map viewer should be \"true\"", () => readingReads(page, "layer \"BaseLayer\" visible", el("map viewer"), "true"));
      await session.step(105, "When user sets \"renderType\" property of map viewer to \"both\"", () => setProperty(page, "renderType", el("map viewer"), "both"));
      await session.step(106, "Then the \"render type\" reading of map viewer should be \"both\"", () => readingReads(page, "render type", el("map viewer"), "both"));
      await session.step(107, "And the \"layer \\\"Heatmap\\\" visible\" reading of map viewer should be \"true\"", () => readingReads(page, "layer \"Heatmap\" visible", el("map viewer"), "true"));
      await session.step(108, "And the \"layer \\\"Markers GL\\\" visible\" reading of map viewer should be \"true\"", () => readingReads(page, "layer \"Markers GL\" visible", el("map viewer"), "true"));
      await session.step(109, "And the \"visible layers\" reading of map viewer should be 5", () => readingIs(page, "visible layers", el("map viewer"), 5));
      await session.step(110, "When user sets \"renderType\" property of map viewer to \"markers\"", () => setProperty(page, "renderType", el("map viewer"), "markers"));
      await session.step(111, "Then the \"render type\" reading of map viewer should be \"markers\"", () => readingReads(page, "render type", el("map viewer"), "markers"));
      await session.step(112, "And the \"layer \\\"Heatmap\\\" visible\" reading of map viewer should be \"false\"", () => readingReads(page, "layer \"Heatmap\" visible", el("map viewer"), "false"));
      await session.step(113, "And the \"layer \\\"Markers GL\\\" visible\" reading of map viewer should be \"true\"", () => readingReads(page, "layer \"Markers GL\" visible", el("map viewer"), "true"));
      await session.step(114, "And the \"visible layers\" reading of map viewer should be 4", () => readingIs(page, "visible layers", el("map viewer"), 4));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The zoom buttons move the view by one level each way", async () => {
      await session.step(118, "Given user remembers the \"zoom\" reading of map viewer", () => rememberReading(page, "zoom", el("map viewer")));
      await session.step(119, "When user clicks on the \"zoom in\" area of map viewer", () => clickArea(page, "zoom in", el("map viewer")));
      await session.step(120, "Then the \"zoom\" reading of map viewer should be higher than before", () => readingHigher(page, "zoom", el("map viewer")));
      await session.step(121, "And the \"zoom\" reading of map viewer should not be as remembered", () => readingNotAsRemembered(page, "zoom", el("map viewer")));
      await session.step(122, "When user clicks on the \"zoom out\" area of map viewer", () => clickArea(page, "zoom out", el("map viewer")));
      await session.step(123, "Then the \"zoom\" reading of map viewer should be lower than before", () => readingLower(page, "zoom", el("map viewer")));
      await session.step(124, "And the \"zoom\" reading of map viewer should be as remembered", () => readingAsRemembered(page, "zoom", el("map viewer")));
      await session.step(125, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl and a drag select the points inside the rectangle, Escape clears them", async () => {
      await session.step(128, "Given user clears the row selection", () => clearSelection(page));
      await session.step(129, "When user drags a box over the \"view\" area of map viewer holding Control", () => dragBoxHolding(page, "view", el("map viewer"), "Control"));
      await session.step(130, "Then some rows should be selected", () => someSelected(page));
      await session.step(131, "And every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(132, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(133, "Then no rows should be selected", () => noneSelected(page));
      await session.step(134, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A table filter narrows what the map holds", async () => {
      await session.step(137, "When user adds a categorical filter on \"MagType\" keeping \"Mw\"", () => addCategoricalFilter(page, "MagType", "Mw"));
      await session.step(138, "Then 2312 rows should pass the filter", () => filterPasses(page, 2312));
      await session.step(139, "And the \"markers\" reading of map viewer should be 2312", () => readingIs(page, "markers", el("map viewer"), 2312));
      await session.step(140, "And the \"rows shown\" reading of map viewer should be 2312", () => readingIs(page, "rows shown", el("map viewer"), 2312));
      await session.step(141, "And the \"markers\" and \"rows shown\" readings of map viewer should be the same", () => readingsEqual(page, "markers", "rows shown", el("map viewer")));
      await session.step(142, "When user hovers over \"MagType\" filter card", () => hoverOver(page, el("\"MagType\" filter card")));
      await session.step(143, "And user clicks on close of \"MagType\" filter card", () => clickOn(page, el("close of \"MagType\" filter card")));
      await session.step(144, "Then 2426 rows should pass the filter", () => filterPasses(page, 2426));
      await session.step(145, "And the \"markers\" reading of map viewer should be 2426", () => readingIs(page, "markers", el("map viewer"), 2426));
      await session.step(146, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Tooltip decides whether a point under the pointer says anything", async () => {
      await session.step(149, "Then \"showTooltip\" property of map viewer should be \"false\"", () => propertyShouldBe(page, "showTooltip", el("map viewer"), "false"));
      await session.step(150, "When user hovers over the \"point 5\" area of map viewer", () => hoverArea(page, "point 5", el("map viewer")));
      await session.step(151, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(152, "When user sets \"showTooltip\" property of map viewer to \"true\"", () => setProperty(page, "showTooltip", el("map viewer"), "true"));
      await session.step(153, "And user moves the pointer away from map viewer", () => pointerAway(page, el("map viewer")));
      await session.step(154, "And user hovers over the \"point 5\" area of map viewer", () => hoverArea(page, "point 5", el("map viewer")));
      await session.step(155, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(156, "And tooltip should contain text \"Latitude\"", () => shouldContainText(page, el("tooltip"), "Latitude"));
      await session.step(157, "And tooltip should contain text \"Longitude\"", () => shouldContainText(page, el("tooltip"), "Longitude"));
      await session.step(158, "And tooltip should contain text \"Magnitude\"", () => shouldContainText(page, el("tooltip"), "Magnitude"));
      await session.step(159, "And tooltip should contain text \"MagType\"", () => shouldContainText(page, el("tooltip"), "MagType"));
      await session.step(160, "When user sets \"showTooltip\" property of map viewer to \"false\"", () => setProperty(page, "showTooltip", el("map viewer"), "false"));
      await session.step(161, "And user moves the pointer away from map viewer", () => pointerAway(page, el("map viewer")));
      await session.step(162, "And user hovers over the \"point 5\" area of map viewer", () => hoverArea(page, "point 5", el("map viewer")));
      await session.step(163, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(164, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the map under the pointer disposes it without an error", async () => {
      await session.step(167, "Given user remembers the place of map viewer", () => rememberPlace(page, el("map viewer")));
      await session.step(168, "When user clicks on close icon of map viewer", () => clickOn(page, el("close icon of map viewer")));
      await session.step(169, "Then map viewer should be absent", () => shouldBe(page, el("map viewer"), "absent"));
      await session.step(170, "When user moves the pointer across the remembered place", () => moveAcrossPlace(page));
      await session.step(171, "Then the open tableview should have 0 map viewers", () => viewerCount(page, 0, "map"));
      await session.step(172, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
