/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-selection-overlays.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart, core#1174]
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
import {clearSelection, filterPasses, filterPassesAll, filterTo, noneSelected, onlyOfAnySelected, onlyOfSelected, removeColumn, resetFilter, rowCount, selectWhereIs, selectWhereOneOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaLessInk, clickArea, hasArea, hasNoArea, noErrors, noHighlight, repaintedBy, setProperties, setProperty, someHighlight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart selected and filtered rows overlays", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-selection-overlays.feature", import.meta.url);
  test("Bar chart selected and filtered rows overlays", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:112", "@realizes:core#1174"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Value","AGE"],["Value Aggr Type","count"],["Row Source","All"],["Show Selected Rows","true"],["Show Filtered Rows","true"]]), [["Split","RACE"],["Value","AGE"],["Value Aggr Type","count"],["Row Source","All"],["Show Selected Rows","true"],["Show Filtered Rows","true"]]);
    await session.step(24, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(25, "And bar chart viewer should show no selection highlight", () => noHighlight(page, el("bar chart viewer")));
    await session.step(26, "And bar chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("bar chart viewer"), "selected Asian"));
    await run.scenario("Under count a selection paints the selected-rows overlay on its bar", async () => {
      await session.step(29, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(30, "Then only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(31, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(32, "And bar chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("bar chart viewer"), "selected Asian"));
      await session.step(33, "And bar chart viewer should not have a \"selected Caucasian\" area", () => hasNoArea(page, el("bar chart viewer"), "selected Caucasian"));
      await session.step(34, "And the \"selected Asian\" area of bar chart viewer should contain the color \"#FF8C00\"", () => areaColor(page, "selected Asian", el("bar chart viewer"), "#FF8C00"));
      await session.step(35, "And bar chart viewer should show a selection highlight", () => someHighlight(page, el("bar chart viewer")));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter shrinks the filtered share and leaves the selection alone", async () => {
      await session.step(39, "Then bar chart viewer should have a \"filtered Caucasian\" area", () => hasArea(page, el("bar chart viewer"), "filtered Caucasian"));
      await session.step(40, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(41, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(42, "And only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(43, "And bar chart viewer should have repainted by at least 100 pixels", () => repaintedBy(page, el("bar chart viewer"), 100));
      await session.step(44, "And the \"filtered Caucasian\" area of bar chart viewer should have less ink than before", () => areaLessInk(page, "filtered Caucasian", el("bar chart viewer")));
      await session.step(45, "And the \"filtered Caucasian\" area of bar chart viewer should contain the color \"#0000A0\"", () => areaColor(page, "filtered Caucasian", el("bar chart viewer"), "#0000A0"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Min drops both overlays and keeps the state", async () => {
      await session.step(49, "When user sets \"Value Aggr Type\" property of bar chart viewer to \"min\"", () => setProperty(page, "Value Aggr Type", el("bar chart viewer"), "min"));
      await session.step(50, "Then bar chart viewer should show no selection highlight", () => noHighlight(page, el("bar chart viewer")));
      await session.step(51, "And bar chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("bar chart viewer"), "selected Asian"));
      await session.step(52, "And bar chart viewer should not have a \"filtered Caucasian\" area", () => hasNoArea(page, el("bar chart viewer"), "filtered Caucasian"));
      await session.step(53, "And only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(54, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sum restores the overlays for a fresh selection and filter", async () => {
      await session.step(58, "When user clears the row selection", () => clearSelection(page));
      await session.step(59, "And user resets the filter", () => resetFilter(page));
      await session.step(60, "Then no rows should be selected", () => noneSelected(page));
      await session.step(61, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(62, "When user sets \"Value Aggr Type\" property of bar chart viewer to \"sum\"", () => setProperty(page, "Value Aggr Type", el("bar chart viewer"), "sum"));
      await session.step(63, "Then bar chart viewer should show no selection highlight", () => noHighlight(page, el("bar chart viewer")));
      await session.step(64, "And bar chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("bar chart viewer"), "selected Asian"));
      await session.step(65, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(66, "And user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(67, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(68, "And only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(69, "And bar chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("bar chart viewer"), "selected Asian"));
      await session.step(70, "And the \"selected Asian\" area of bar chart viewer should contain the color \"#FF8C00\"", () => areaColor(page, "selected Asian", el("bar chart viewer"), "#FF8C00"));
      await session.step(71, "And the \"filtered Caucasian\" area of bar chart viewer should have less ink than before", () => areaLessInk(page, "filtered Caucasian", el("bar chart viewer")));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Value count is cumulative too", async () => {
      await session.step(75, "When user sets \"Value Aggr Type\" property of bar chart viewer to \"values\"", () => setProperty(page, "Value Aggr Type", el("bar chart viewer"), "values"));
      await session.step(76, "And user selects rows where \"RACE\" is one of \"Asian, Caucasian\"", () => selectWhereOneOf(page, "RACE", "Asian, Caucasian"));
      await session.step(77, "Then only rows where \"RACE\" is one of \"Asian, Caucasian\" should be selected", () => onlyOfAnySelected(page, "RACE", "Asian, Caucasian"));
      await session.step(78, "And bar chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("bar chart viewer"), "selected Asian"));
      await session.step(79, "And bar chart viewer should have a \"selected Caucasian\" area", () => hasArea(page, el("bar chart viewer"), "selected Caucasian"));
      await session.step(80, "And bar chart viewer should show a selection highlight", () => someHighlight(page, el("bar chart viewer")));
      await session.step(81, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
      await session.step(82, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(83, "And bar chart viewer should have repainted by at least 100 pixels", () => repaintedBy(page, el("bar chart viewer"), 100));
      await session.step(84, "And only rows where \"RACE\" is one of \"Asian, Caucasian\" should be selected", () => onlyOfAnySelected(page, "RACE", "Asian, Caucasian"));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
      await session.step(86, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
    });
    await run.scenario("Average is not", async () => {
      await session.step(89, "When user sets \"Value Aggr Type\" property of bar chart viewer to \"avg\"", () => setProperty(page, "Value Aggr Type", el("bar chart viewer"), "avg"));
      await session.step(90, "Then bar chart viewer should show no selection highlight", () => noHighlight(page, el("bar chart viewer")));
      await session.step(91, "And bar chart viewer should not have a \"selected Caucasian\" area", () => hasNoArea(page, el("bar chart viewer"), "selected Caucasian"));
      await session.step(92, "And bar chart viewer should not have a \"filtered Caucasian\" area", () => hasNoArea(page, el("bar chart viewer"), "filtered Caucasian"));
      await session.step(93, "And only rows where \"RACE\" is one of \"Asian, Caucasian\" should be selected", () => onlyOfAnySelected(page, "RACE", "Asian, Caucasian"));
      await session.step(94, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
      await session.step(96, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Show Selected Rows","false"],["Show Filtered Rows","false"],["Value Aggr Type","count"]]), [["Show Selected Rows","false"],["Show Filtered Rows","false"],["Value Aggr Type","count"]]);
      await session.step(100, "And user clears the row selection", () => clearSelection(page));
      await session.step(101, "And user resets the filter", () => resetFilter(page));
      await session.step(102, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(103, "And no rows should be selected", () => noneSelected(page));
      await session.step(104, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the column of the category a click filtered by gives the rows back (core#1174)", async () => {
      await session.step(108, "When user sets \"onClick\" property of bar chart viewer to \"Filter\"", () => setProperty(page, "onClick", el("bar chart viewer"), "Filter"));
      await session.step(109, "And user clicks on the \"bar Asian\" area of bar chart viewer", () => clickArea(page, "bar Asian", el("bar chart viewer")));
      await session.step(110, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(111, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(112, "Then all rows should pass the filter", () => filterPassesAll(page));
    }, {knownFailure: true, line: 112});
    run.finish();
  });
});
