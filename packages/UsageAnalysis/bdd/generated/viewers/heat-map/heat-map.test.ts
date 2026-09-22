/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map]
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
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaShorter, areaTaller, hasArea, hasNoArea, noErrors, painted, paintedInColors, readingAsRemembered, readingBetween, readingHigher, readingIs, readingLower, readingReads, rememberReading, repainted, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map layout, column labels, column cap and scrollbars", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map.feature", import.meta.url);
  test("Heat map layout, column labels, column cap and scrollbars", {tag: ["@journey", "@viewers", "@realizes:viewers.heat-map"]}, async ({browser}) => {
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
    await run.scenario("In heat map mode the Row Height property row is disabled", async () => {
      await session.step(39, "When user clicks on settings icon of heat map viewer", () => clickOn(page, el("settings icon of heat map viewer")));
      await session.step(40, "Then \"Row Height\" property in context panel should be present", () => shouldBe(page, el("\"Row Height\" property in context panel"), "present"));
      await session.step(41, "And \"Row Height\" property in context panel should be disabled", () => shouldBe(page, el("\"Row Height\" property in context panel"), "disabled"));
    });
    await run.scenario("Col Labels Orientation is what was asked for, and what the layout resolved", async () => {
      await session.step(44, "Then the \"col labels orientation\" reading of heat map viewer should be \"Auto\"", () => readingReads(page, "col labels orientation", el("heat map viewer"), "Auto"));
      await session.step(45, "When user remembers the \"effective col labels orientation\" reading of heat map viewer", () => rememberReading(page, "effective col labels orientation", el("heat map viewer")));
      await session.step(46, "And user sets \"colLabelsOrientation\" property of heat map viewer to \"Vert\"", () => setProperty(page, "colLabelsOrientation", el("heat map viewer"), "Vert"));
      await session.step(47, "Then the \"col labels orientation\" reading of heat map viewer should be \"Vert\"", () => readingReads(page, "col labels orientation", el("heat map viewer"), "Vert"));
      await session.step(48, "And the \"effective col labels orientation\" reading of heat map viewer should be \"Vert\"", () => readingReads(page, "effective col labels orientation", el("heat map viewer"), "Vert"));
      await session.step(49, "And the \"header AGE\" area of heat map viewer should be taller than before", () => areaTaller(page, "header AGE", el("heat map viewer")));
      await session.step(50, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(51, "When user sets \"colLabelsOrientation\" property of heat map viewer to \"Horz\"", () => setProperty(page, "colLabelsOrientation", el("heat map viewer"), "Horz"));
      await session.step(52, "Then the \"col labels orientation\" reading of heat map viewer should be \"Horz\"", () => readingReads(page, "col labels orientation", el("heat map viewer"), "Horz"));
      await session.step(53, "And the \"effective col labels orientation\" reading of heat map viewer should be \"Horz\"", () => readingReads(page, "effective col labels orientation", el("heat map viewer"), "Horz"));
      await session.step(54, "And the \"header AGE\" area of heat map viewer should be shorter than before", () => areaShorter(page, "header AGE", el("heat map viewer")));
      await session.step(55, "When user sets \"colLabelsOrientation\" property of heat map viewer to \"Auto\"", () => setProperty(page, "colLabelsOrientation", el("heat map viewer"), "Auto"));
      await session.step(56, "Then the \"col labels orientation\" reading of heat map viewer should be \"Auto\"", () => readingReads(page, "col labels orientation", el("heat map viewer"), "Auto"));
      await session.step(57, "And the \"effective col labels orientation\" reading of heat map viewer should be as remembered", () => readingAsRemembered(page, "effective col labels orientation", el("heat map viewer")));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Heatmap Scrollbars decides whether the sliders are drawn, not what they span", async () => {
      await session.step(61, "Then heat map viewer should have an \"x scroll slider\" area", () => hasArea(page, el("heat map viewer"), "x scroll slider"));
      await session.step(62, "And heat map viewer should have a \"y scroll slider\" area", () => hasArea(page, el("heat map viewer"), "y scroll slider"));
      await session.step(63, "And heat map viewer should have a \"y scroll min handle\" area", () => hasArea(page, el("heat map viewer"), "y scroll min handle"));
      await session.step(64, "And heat map viewer should have a \"y scroll max handle\" area", () => hasArea(page, el("heat map viewer"), "y scroll max handle"));
      await session.step(65, "And the \"x scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "x scroll span", el("heat map viewer"), 1));
      await session.step(66, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
      await session.step(67, "When user sets \"showHeatmapScrollbars\" property of heat map viewer to \"false\"", () => setProperty(page, "showHeatmapScrollbars", el("heat map viewer"), "false"));
      await session.step(68, "Then heat map viewer should not have an \"x scroll slider\" area", () => hasNoArea(page, el("heat map viewer"), "x scroll slider"));
      await session.step(69, "And heat map viewer should not have a \"y scroll slider\" area", () => hasNoArea(page, el("heat map viewer"), "y scroll slider"));
      await session.step(70, "And heat map viewer should not have a \"y scroll min handle\" area", () => hasNoArea(page, el("heat map viewer"), "y scroll min handle"));
      await session.step(71, "And the \"x scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "x scroll span", el("heat map viewer"), 1));
      await session.step(72, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
      await session.step(73, "When user sets \"showHeatmapScrollbars\" property of heat map viewer to \"true\"", () => setProperty(page, "showHeatmapScrollbars", el("heat map viewer"), "true"));
      await session.step(74, "Then heat map viewer should have an \"x scroll slider\" area", () => hasArea(page, el("heat map viewer"), "x scroll slider"));
      await session.step(75, "And heat map viewer should have a \"y scroll slider\" area", () => hasArea(page, el("heat map viewer"), "y scroll slider"));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter on the table leaves fewer, thicker rows on screen", async () => {
      await session.step(79, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(80, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(81, "And the \"rows shown\" reading of heat map viewer should be 447", () => readingIs(page, "rows shown", el("heat map viewer"), 447));
      await session.step(82, "And the \"row height\" reading of heat map viewer should be higher than before", () => readingHigher(page, "row height", el("heat map viewer")));
      await session.step(83, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(84, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(85, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(86, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(87, "And the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
      await session.step(88, "And the \"row height\" reading of heat map viewer should be lower than before", () => readingLower(page, "row height", el("heat map viewer")));
      await session.step(89, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Max Heatmap Columns set after opening settings takes columns off the screen", async () => {
      await session.step(96, "When user clicks on settings icon of heat map viewer", () => clickOn(page, el("settings icon of heat map viewer")));
      await session.step(97, "Then the \"max heatmap columns\" reading of heat map viewer should be 100", () => readingIs(page, "max heatmap columns", el("heat map viewer"), 100));
      await session.step(98, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(99, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"3\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "3"));
      await session.step(100, "Then the \"max heatmap columns\" reading of heat map viewer should be 3", () => readingIs(page, "max heatmap columns", el("heat map viewer"), 3));
      await session.step(101, "And the \"columns shown\" reading of heat map viewer should be 3", () => readingIs(page, "columns shown", el("heat map viewer"), 3));
      await session.step(102, "And the \"column order\" reading of heat map viewer should be \"AGE, HEIGHT, WEIGHT\"", () => readingReads(page, "column order", el("heat map viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(103, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(104, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"100\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "100"));
      await session.step(105, "Then the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title bar closes the map", async () => {
      await session.step(109, "When user clicks on close icon of heat map viewer", () => clickOn(page, el("close icon of heat map viewer")));
      await session.step(110, "Then heat map viewer should be absent", () => shouldBe(page, el("heat map viewer"), "absent"));
      await session.step(111, "And the open tableview should have 0 heat map viewers", () => viewerCount(page, 0, "heat map"));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Max Heatmap Columns applies before opening settings", async () => {
      await session.step(118, "Given user adds a heat map viewer", () => addViewer(page, "heat map"));
      await session.step(119, "Then the \"max heatmap columns\" reading of heat map viewer should be 100", () => readingIs(page, "max heatmap columns", el("heat map viewer"), 100));
      await session.step(120, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(121, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"3\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "3"));
      await session.step(122, "Then the \"max heatmap columns\" reading of heat map viewer should be 3", () => readingIs(page, "max heatmap columns", el("heat map viewer"), 3));
      await session.step(123, "And the \"columns shown\" reading of heat map viewer should be 3", () => readingIs(page, "columns shown", el("heat map viewer"), 3));
    });
    run.finish();
  });
});
