/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-selection-overlays.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearSelection, filterPasses, filterPassesAll, filterTo, noneSelected, onlyOfAnySelected, onlyOfSelected, resetFilter, rowCount, selectWhereIs, selectWhereOneOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaLessInk, hasArea, hasNoArea, noErrors, noHighlight, repaintedBy, setProperties, setProperty, someHighlight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart selected and filtered rows overlays", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-selection-overlays.feature", import.meta.url);
  test("Bar chart selected and filtered rows overlays", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Value","AGE"],["Value Aggr Type","count"],["Row Source","All"],["Show Selected Rows","true"],["Show Filtered Rows","true"]]));
    await session.step(20, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(21, "And bar chart viewer should show no selection highlight", () => noHighlight(page, el("bar chart viewer")));
    await session.step(22, "And bar chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("bar chart viewer"), "selected Asian"));
    await run.scenario("Under count a selection paints the selected-rows overlay on its bar", async () => {
      await session.step(25, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(26, "Then only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(27, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(28, "And bar chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("bar chart viewer"), "selected Asian"));
      await session.step(29, "And bar chart viewer should not have a \"selected Caucasian\" area", () => hasNoArea(page, el("bar chart viewer"), "selected Caucasian"));
      await session.step(30, "And the \"selected Asian\" area of bar chart viewer should contain the color \"#FF8C00\"", () => areaColor(page, "selected Asian", el("bar chart viewer"), "#FF8C00"));
      await session.step(31, "And bar chart viewer should show a selection highlight", () => someHighlight(page, el("bar chart viewer")));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter shrinks the filtered share and leaves the selection alone", async () => {
      await session.step(35, "Then bar chart viewer should have a \"filtered Caucasian\" area", () => hasArea(page, el("bar chart viewer"), "filtered Caucasian"));
      await session.step(36, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(37, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(38, "And only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(39, "And bar chart viewer should have repainted by at least 100 pixels", () => repaintedBy(page, el("bar chart viewer"), 100));
      await session.step(40, "And the \"filtered Caucasian\" area of bar chart viewer should have less ink than before", () => areaLessInk(page, "filtered Caucasian", el("bar chart viewer")));
      await session.step(41, "And the \"filtered Caucasian\" area of bar chart viewer should contain the color \"#0000A0\"", () => areaColor(page, "filtered Caucasian", el("bar chart viewer"), "#0000A0"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Min drops both overlays and keeps the state", async () => {
      await session.step(45, "When user sets \"Value Aggr Type\" property of bar chart viewer to \"min\"", () => setProperty(page, "Value Aggr Type", el("bar chart viewer"), "min"));
      await session.step(46, "Then bar chart viewer should show no selection highlight", () => noHighlight(page, el("bar chart viewer")));
      await session.step(47, "And bar chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("bar chart viewer"), "selected Asian"));
      await session.step(48, "And bar chart viewer should not have a \"filtered Caucasian\" area", () => hasNoArea(page, el("bar chart viewer"), "filtered Caucasian"));
      await session.step(49, "And only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(50, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sum restores the overlays for a fresh selection and filter", async () => {
      await session.step(54, "When user clears the row selection", () => clearSelection(page));
      await session.step(55, "And user resets the filter", () => resetFilter(page));
      await session.step(56, "Then no rows should be selected", () => noneSelected(page));
      await session.step(57, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(58, "When user sets \"Value Aggr Type\" property of bar chart viewer to \"sum\"", () => setProperty(page, "Value Aggr Type", el("bar chart viewer"), "sum"));
      await session.step(59, "Then bar chart viewer should show no selection highlight", () => noHighlight(page, el("bar chart viewer")));
      await session.step(60, "And bar chart viewer should not have a \"selected Asian\" area", () => hasNoArea(page, el("bar chart viewer"), "selected Asian"));
      await session.step(61, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(62, "And user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(63, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(64, "And only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(65, "And bar chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("bar chart viewer"), "selected Asian"));
      await session.step(66, "And the \"selected Asian\" area of bar chart viewer should contain the color \"#FF8C00\"", () => areaColor(page, "selected Asian", el("bar chart viewer"), "#FF8C00"));
      await session.step(67, "And the \"filtered Caucasian\" area of bar chart viewer should have less ink than before", () => areaLessInk(page, "filtered Caucasian", el("bar chart viewer")));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Value count is cumulative too", async () => {
      await session.step(71, "When user sets \"Value Aggr Type\" property of bar chart viewer to \"values\"", () => setProperty(page, "Value Aggr Type", el("bar chart viewer"), "values"));
      await session.step(72, "And user selects rows where \"RACE\" is one of \"Asian, Caucasian\"", () => selectWhereOneOf(page, "RACE", "Asian, Caucasian"));
      await session.step(73, "Then only rows where \"RACE\" is one of \"Asian, Caucasian\" should be selected", () => onlyOfAnySelected(page, "RACE", "Asian, Caucasian"));
      await session.step(74, "And bar chart viewer should have a \"selected Asian\" area", () => hasArea(page, el("bar chart viewer"), "selected Asian"));
      await session.step(75, "And bar chart viewer should have a \"selected Caucasian\" area", () => hasArea(page, el("bar chart viewer"), "selected Caucasian"));
      await session.step(76, "And bar chart viewer should show a selection highlight", () => someHighlight(page, el("bar chart viewer")));
      await session.step(77, "When user filters rows where \"SEX\" is \"M\"", () => filterTo(page, "SEX", "M"));
      await session.step(78, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(79, "And bar chart viewer should have repainted by at least 100 pixels", () => repaintedBy(page, el("bar chart viewer"), 100));
      await session.step(80, "And only rows where \"RACE\" is one of \"Asian, Caucasian\" should be selected", () => onlyOfAnySelected(page, "RACE", "Asian, Caucasian"));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
      await session.step(82, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
    });
    await run.scenario("Average is not", async () => {
      await session.step(85, "When user sets \"Value Aggr Type\" property of bar chart viewer to \"avg\"", () => setProperty(page, "Value Aggr Type", el("bar chart viewer"), "avg"));
      await session.step(86, "Then bar chart viewer should show no selection highlight", () => noHighlight(page, el("bar chart viewer")));
      await session.step(87, "And bar chart viewer should not have a \"selected Caucasian\" area", () => hasNoArea(page, el("bar chart viewer"), "selected Caucasian"));
      await session.step(88, "And bar chart viewer should not have a \"filtered Caucasian\" area", () => hasNoArea(page, el("bar chart viewer"), "filtered Caucasian"));
      await session.step(89, "And only rows where \"RACE\" is one of \"Asian, Caucasian\" should be selected", () => onlyOfAnySelected(page, "RACE", "Asian, Caucasian"));
      await session.step(90, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
      await session.step(92, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Show Selected Rows","false"],["Show Filtered Rows","false"],["Value Aggr Type","count"]]));
      await session.step(96, "And user clears the row selection", () => clearSelection(page));
      await session.step(97, "And user resets the filter", () => resetFilter(page));
      await session.step(98, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(99, "And no rows should be selected", () => noneSelected(page));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
