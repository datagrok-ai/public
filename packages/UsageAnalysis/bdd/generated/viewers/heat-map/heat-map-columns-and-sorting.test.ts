/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map-columns-and-sorting.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map]
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
import {clickOn, expand, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {maxInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areasSameSize, doubleClickArea, dragAreaBy, hasArea, hasNoArea, noErrors, painted, readingBetween, readingHigherThanRemembered, readingIs, readingReads, rememberReading, repainted, reportsNoError, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {checkedInColumnList, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map columns, sorting and column widths", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map-columns-and-sorting.feature", import.meta.url);
  test("Heat map columns, sorting and column widths", {tag: ["@journey", "@viewers", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(21, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(22, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
    await session.step(23, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await run.scenario("The Columns property counts the visible columns and its dialog hides one from the map (core#376)", async () => {
      await session.step(26, "When user clicks on settings icon of heat map viewer", () => clickOn(page, el("settings icon of heat map viewer")));
      await session.step(27, "And user expands \"Columns\" category in context panel", () => expand(page, el("\"Columns\" category in context panel")));
      await session.step(28, "Then \"Columns\" property in context panel should have value \"11 visible\"", () => shouldHaveValue(page, el("\"Columns\" property in context panel"), "11 visible"));
      await session.step(29, "When user clicks on first button in \"Columns\" property in context panel", () => clickOn(page, el("first button in \"Columns\" property in context panel")));
      await session.step(30, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
      await session.step(31, "And the \"AGE\" column should be checked in the column list of Order or Hide Columns dialog", () => checkedInColumnList(page, "AGE", el("Order or Hide Columns dialog")));
      await session.step(32, "When user takes a snapshot of heat map viewer", () => takeSnapshot(page, el("heat map viewer")));
      await session.step(33, "And user toggles the \"AGE\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "AGE", el("Order or Hide Columns dialog")));
      await session.step(34, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
      await session.step(35, "Then the \"columns shown\" reading of heat map viewer should be 10", () => readingIs(page, "columns shown", el("heat map viewer"), 10));
      await session.step(36, "And the \"column order\" reading of heat map viewer should be \"USUBJID, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("heat map viewer"), "USUBJID, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(37, "And heat map viewer should not have a \"header AGE\" area", () => hasNoArea(page, el("heat map viewer"), "header AGE"));
      await session.step(38, "And heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
      await session.step(39, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(40, "And \"Columns\" property in context panel should have value \"10 visible\"", () => shouldHaveValue(page, el("\"Columns\" property in context panel"), "10 visible"));
      await session.step(41, "When user clicks on first button in \"Columns\" property in context panel", () => clickOn(page, el("first button in \"Columns\" property in context panel")));
      await session.step(42, "And user toggles the \"AGE\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "AGE", el("Order or Hide Columns dialog")));
      await session.step(43, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
      await session.step(44, "Then the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(45, "And the \"column order\" reading of heat map viewer should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("heat map viewer"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(46, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
      await session.step(47, "And \"Columns\" property in context panel should have value \"11 visible\"", () => shouldHaveValue(page, el("\"Columns\" property in context panel"), "11 visible"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A header double-click sorts descending, then ascending, then puts the table's order back (public#182)", async () => {
      await session.step(51, "Then the \"sort column\" reading of heat map viewer should be \"\"", () => readingReads(page, "sort column", el("heat map viewer"), ""));
      await session.step(52, "And the \"row order\" reading of heat map viewer should be \"1, 2, 3, 4, 5, 6, 7, 8, 9, 10\"", () => readingReads(page, "row order", el("heat map viewer"), "1, 2, 3, 4, 5, 6, 7, 8, 9, 10"));
      await session.step(53, "When user double-clicks on the \"header AGE\" area of heat map viewer", () => doubleClickArea(page, "header AGE", el("heat map viewer")));
      await session.step(54, "Then the \"sort column\" reading of heat map viewer should be \"AGE\"", () => readingReads(page, "sort column", el("heat map viewer"), "AGE"));
      await session.step(55, "And the \"sort direction\" reading of heat map viewer should be \"descending\"", () => readingReads(page, "sort direction", el("heat map viewer"), "descending"));
      await session.step(56, "And the \"row order\" reading of heat map viewer should be \"600, 720, 597, 469, 577, 639, 645, 848, 620, 440\"", () => readingReads(page, "row order", el("heat map viewer"), "600, 720, 597, 469, 577, 639, 645, 848, 620, 440"));
      await session.step(57, "And \"AGE\" column should have its maximum in row 600", () => maxInRow(page, "AGE", 600));
      await session.step(58, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(59, "When user double-clicks on the \"header AGE\" area of heat map viewer", () => doubleClickArea(page, "header AGE", el("heat map viewer")));
      await session.step(60, "Then the \"sort direction\" reading of heat map viewer should be \"ascending\"", () => readingReads(page, "sort direction", el("heat map viewer"), "ascending"));
      await session.step(61, "And the \"row order\" reading of heat map viewer should be \"692, 816, 89, 460, 479, 988, 62, 672, 795, 893\"", () => readingReads(page, "row order", el("heat map viewer"), "692, 816, 89, 460, 479, 988, 62, 672, 795, 893"));
      await session.step(62, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(63, "When user double-clicks on the \"header AGE\" area of heat map viewer", () => doubleClickArea(page, "header AGE", el("heat map viewer")));
      await session.step(64, "Then the \"sort column\" reading of heat map viewer should be \"\"", () => readingReads(page, "sort column", el("heat map viewer"), ""));
      await session.step(65, "And the \"sort direction\" reading of heat map viewer should be \"\"", () => readingReads(page, "sort direction", el("heat map viewer"), ""));
      await session.step(66, "And the \"row order\" reading of heat map viewer should be \"1, 2, 3, 4, 5, 6, 7, 8, 9, 10\"", () => readingReads(page, "row order", el("heat map viewer"), "1, 2, 3, 4, 5, 6, 7, 8, 9, 10"));
      await session.step(67, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column border drag does not resize a heat map column, and does once Is Heatmap is off (GROK-18672)", async () => {
      await session.step(73, "Then the \"column width of AGE\" reading of heat map viewer should be 30", () => readingIs(page, "column width of AGE", el("heat map viewer"), 30));
      await session.step(74, "When user drags the \"column resizer AGE\" area of heat map viewer by 40 pixels to the right", () => dragAreaBy(page, "column resizer AGE", el("heat map viewer"), 40, "right"));
      await session.step(75, "Then the \"column width of AGE\" reading of heat map viewer should be 30", () => readingIs(page, "column width of AGE", el("heat map viewer"), 30));
      await session.step(76, "And the \"column width of SEX\" reading of heat map viewer should be 30", () => readingIs(page, "column width of SEX", el("heat map viewer"), 30));
      await session.step(77, "And the \"column AGE\" and \"column SEX\" areas of heat map viewer should be the same width", () => areasSameSize(page, "column AGE", "column SEX", el("heat map viewer"), "width"));
      await session.step(78, "When user sets \"isHeatmap\" property of heat map viewer to \"false\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "false"));
      await session.step(79, "Then the \"is heatmap\" reading of heat map viewer should be \"false\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "false"));
      await session.step(80, "When user remembers the \"column width of AGE\" reading of heat map viewer", () => rememberReading(page, "column width of AGE", el("heat map viewer")));
      await session.step(81, "And user drags the \"column resizer AGE\" area of heat map viewer by 40 pixels to the right", () => dragAreaBy(page, "column resizer AGE", el("heat map viewer"), 40, "right"));
      await session.step(82, "Then the \"column width of AGE\" reading of heat map viewer should be higher than remembered", () => readingHigherThanRemembered(page, "column width of AGE", el("heat map viewer")));
      await session.step(83, "When user sets \"isHeatmap\" property of heat map viewer to \"true\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "true"));
      await session.step(84, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
      await session.step(85, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
      await session.step(86, "And the \"column width of AGE\" reading of heat map viewer should be 30", () => readingIs(page, "column width of AGE", el("heat map viewer"), 30));
      await session.step(87, "And the \"column AGE\" and \"column SEX\" areas of heat map viewer should be the same width", () => areasSameSize(page, "column AGE", "column SEX", el("heat map viewer"), "width"));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Max Heatmap Columns 0 leaves no column on the map and no error (GROK-17776)", async () => {
      await session.step(91, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"0\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "0"));
      await session.step(92, "Then the \"max heatmap columns\" reading of heat map viewer should be 0", () => readingIs(page, "max heatmap columns", el("heat map viewer"), 0));
      await session.step(93, "And the \"columns shown\" reading of heat map viewer should be 0", () => readingIs(page, "columns shown", el("heat map viewer"), 0));
      await session.step(94, "And the \"column order\" reading of heat map viewer should be \"\"", () => readingReads(page, "column order", el("heat map viewer"), ""));
      await session.step(95, "And heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
      await session.step(96, "And heat map viewer should report no error", () => reportsNoError(page, el("heat map viewer")));
      await session.step(97, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(98, "When user sets \"maxHeatmapColumns\" property of heat map viewer to \"100\"", () => setProperty(page, "maxHeatmapColumns", el("heat map viewer"), "100"));
      await session.step(99, "Then the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(100, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
      await session.step(101, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
