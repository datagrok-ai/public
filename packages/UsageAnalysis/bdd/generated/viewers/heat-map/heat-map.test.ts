/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaShorter, areaTaller, hasArea, hasNoArea, noErrors, painted, paintedInColors, readingAsRemembered, readingBetween, readingHigher, readingIs, readingLower, readingReads, rememberReading, repainted, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map layout, column labels, column cap and scrollbars", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map.feature", import.meta.url);
  test("Heat map layout, column labels, column cap and scrollbars", {tag: ["@journey", "@viewers", "@realizes:viewers.heat-map", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(18, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(19, "And the \"rows\" reading of heat map viewer should be 1000", () => readingIs(page, "rows", el("heat map viewer"), 1000));
    await session.step(20, "And the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
    await session.step(21, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await run.scenario("Every row and every data column is on screen, as bands rather than cells", async () => {
      await session.step(24, "Then the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(25, "And the \"column order\" reading of heat map viewer should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("heat map viewer"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(26, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
      await session.step(27, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
      await session.step(28, "And heat map viewer should have a \"header AGE\" area", () => hasArea(page, el("heat map viewer"), "header AGE"));
      await session.step(29, "And heat map viewer should not have a \"cell 1 of AGE\" area", () => hasNoArea(page, el("heat map viewer"), "cell 1 of AGE"));
      await session.step(30, "And heat map viewer should be painted in at least 3 colors", () => paintedInColors(page, el("heat map viewer"), 3));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("In heat map mode the Row Height property row is disabled (property_grid_lib.dart:583-612)", async () => {
      await session.step(43, "When user clicks on settings icon of heat map viewer", () => clickOn(page, el("settings icon of heat map viewer")));
      await session.step(44, "Then \"Row Height\" property in context panel should be present", () => shouldBe(page, el("\"Row Height\" property in context panel"), "present"));
      await session.step(45, "And \"Row Height\" property in context panel should be disabled", () => shouldBe(page, el("\"Row Height\" property in context panel"), "disabled"));
    }, {knownFailure: true});
    await run.scenario("Col Labels Orientation is what was asked for, and what the layout resolved", async () => {
      await session.step(48, "Then the \"col labels orientation\" reading of heat map viewer should be \"Auto\"", () => readingReads(page, "col labels orientation", el("heat map viewer"), "Auto"));
      await session.step(49, "When user remembers the \"effective col labels orientation\" reading of heat map viewer", () => rememberReading(page, "effective col labels orientation", el("heat map viewer")));
      await session.step(50, "And user sets \"colLabelsOrientation\" property of heat map viewer to \"Vert\"", () => setProperty(page, "colLabelsOrientation", el("heat map viewer"), "Vert"));
      await session.step(51, "Then the \"col labels orientation\" reading of heat map viewer should be \"Vert\"", () => readingReads(page, "col labels orientation", el("heat map viewer"), "Vert"));
      await session.step(52, "And the \"effective col labels orientation\" reading of heat map viewer should be \"Vert\"", () => readingReads(page, "effective col labels orientation", el("heat map viewer"), "Vert"));
      await session.step(53, "And the \"header AGE\" area of heat map viewer should be taller than before", () => areaTaller(page, "header AGE", el("heat map viewer")));
      await session.step(54, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(55, "When user sets \"colLabelsOrientation\" property of heat map viewer to \"Horz\"", () => setProperty(page, "colLabelsOrientation", el("heat map viewer"), "Horz"));
      await session.step(56, "Then the \"col labels orientation\" reading of heat map viewer should be \"Horz\"", () => readingReads(page, "col labels orientation", el("heat map viewer"), "Horz"));
      await session.step(57, "And the \"effective col labels orientation\" reading of heat map viewer should be \"Horz\"", () => readingReads(page, "effective col labels orientation", el("heat map viewer"), "Horz"));
      await session.step(58, "And the \"header AGE\" area of heat map viewer should be shorter than before", () => areaShorter(page, "header AGE", el("heat map viewer")));
      await session.step(59, "When user sets \"colLabelsOrientation\" property of heat map viewer to \"Auto\"", () => setProperty(page, "colLabelsOrientation", el("heat map viewer"), "Auto"));
      await session.step(60, "Then the \"col labels orientation\" reading of heat map viewer should be \"Auto\"", () => readingReads(page, "col labels orientation", el("heat map viewer"), "Auto"));
      await session.step(61, "And the \"effective col labels orientation\" reading of heat map viewer should be as remembered", () => readingAsRemembered(page, "effective col labels orientation", el("heat map viewer")));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Heatmap Scrollbars decides whether the sliders are drawn, not what they span", async () => {
      await session.step(65, "Then heat map viewer should have an \"x scroll slider\" area", () => hasArea(page, el("heat map viewer"), "x scroll slider"));
      await session.step(66, "And heat map viewer should have a \"y scroll slider\" area", () => hasArea(page, el("heat map viewer"), "y scroll slider"));
      await session.step(67, "And heat map viewer should have a \"y scroll min handle\" area", () => hasArea(page, el("heat map viewer"), "y scroll min handle"));
      await session.step(68, "And heat map viewer should have a \"y scroll max handle\" area", () => hasArea(page, el("heat map viewer"), "y scroll max handle"));
      await session.step(69, "And the \"x scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "x scroll span", el("heat map viewer"), 1));
      await session.step(70, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
      await session.step(71, "When user sets \"showHeatmapScrollbars\" property of heat map viewer to \"false\"", () => setProperty(page, "showHeatmapScrollbars", el("heat map viewer"), "false"));
      await session.step(72, "Then heat map viewer should not have an \"x scroll slider\" area", () => hasNoArea(page, el("heat map viewer"), "x scroll slider"));
      await session.step(73, "And heat map viewer should not have a \"y scroll slider\" area", () => hasNoArea(page, el("heat map viewer"), "y scroll slider"));
      await session.step(74, "And heat map viewer should not have a \"y scroll min handle\" area", () => hasNoArea(page, el("heat map viewer"), "y scroll min handle"));
      await session.step(75, "And the \"x scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "x scroll span", el("heat map viewer"), 1));
      await session.step(76, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
      await session.step(77, "When user sets \"showHeatmapScrollbars\" property of heat map viewer to \"true\"", () => setProperty(page, "showHeatmapScrollbars", el("heat map viewer"), "true"));
      await session.step(78, "Then heat map viewer should have an \"x scroll slider\" area", () => hasArea(page, el("heat map viewer"), "x scroll slider"));
      await session.step(79, "And heat map viewer should have a \"y scroll slider\" area", () => hasArea(page, el("heat map viewer"), "y scroll slider"));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter on the table leaves fewer, thicker rows on screen", async () => {
      await session.step(83, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(84, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(85, "And the \"rows shown\" reading of heat map viewer should be 447", () => readingIs(page, "rows shown", el("heat map viewer"), 447));
      await session.step(86, "And the \"row height\" reading of heat map viewer should be higher than before", () => readingHigher(page, "row height", el("heat map viewer")));
      await session.step(87, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(88, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(89, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(90, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(91, "And the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
      await session.step(92, "And the \"row height\" reading of heat map viewer should be lower than before", () => readingLower(page, "row height", el("heat map viewer")));
      await session.step(93, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Max Heatmap Columns set through the property panel takes columns off the screen", async () => {
      await session.step(102, "When user clicks on settings icon of heat map viewer", () => clickOn(page, el("settings icon of heat map viewer")));
      await session.step(103, "Then the \"max heatmap columns\" reading of heat map viewer should be 100", () => readingIs(page, "max heatmap columns", el("heat map viewer"), 100));
      await session.step(104, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(105, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"3\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "3"));
      await session.step(106, "Then the \"max heatmap columns\" reading of heat map viewer should be 3", () => readingIs(page, "max heatmap columns", el("heat map viewer"), 3));
      await session.step(107, "And the \"columns shown\" reading of heat map viewer should be 3", () => readingIs(page, "columns shown", el("heat map viewer"), 3));
      await session.step(108, "And the \"column order\" reading of heat map viewer should be \"AGE, HEIGHT, WEIGHT\"", () => readingReads(page, "column order", el("heat map viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(109, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(110, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"100\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "100"));
      await session.step(111, "Then the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title bar closes the map", async () => {
      await session.step(115, "When user clicks on close icon of heat map viewer", () => clickOn(page, el("close icon of heat map viewer")));
      await session.step(116, "Then heat map viewer should be absent", () => shouldBe(page, el("heat map viewer"), "absent"));
      await session.step(117, "And the open tableview should have 0 heat map viewers", () => viewerCount(page, 0, "heat map"));
      await session.step(118, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The same cap written through the JS API alone rebuilds nothing (grid_look.dart:331)", async () => {
      await session.step(129, "Given user adds a heat map viewer", () => addViewer(page, "heat map"));
      await session.step(130, "Then the \"max heatmap columns\" reading of heat map viewer should be 100", () => readingIs(page, "max heatmap columns", el("heat map viewer"), 100));
      await session.step(131, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(132, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"3\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "3"));
      await session.step(133, "Then the \"max heatmap columns\" reading of heat map viewer should be 3", () => readingIs(page, "max heatmap columns", el("heat map viewer"), 3));
      await session.step(134, "And the \"columns shown\" reading of heat map viewer should be 3", () => readingIs(page, "columns shown", el("heat map viewer"), 3));
    }, {knownFailure: true});
    run.finish();
  });
});
