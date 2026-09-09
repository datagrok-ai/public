/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-nulls-and-highlight.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, noneOfSelected, noneSelected, removeColumn, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hasArea, hasNoArea, hoverArea, legendLists, moreHighlight, noErrors, painted, pointerAway, readingBetween, readingIs, repainted, repaintedBy, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart missing values and the mouse-over row group", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-nulls-and-highlight.feature", import.meta.url);
  test("Pie chart missing values and the mouse-over row group", {tag: ["@journey", "@viewers", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a calculated column \"RACE_GAPS\" with formula \"if(Mod(${AGE}, 10) == 0, null, ${RACE})\"", () => addCalculated(page, "RACE_GAPS", "if(Mod(${AGE}, 10) == 0, null, ${RACE})"));
    await session.step(18, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE_GAPS"],["Legend Visibility","Always"]]));
    await session.step(21, "Then the \"slices\" reading of pie chart viewer should be 5", () => readingIs(page, "slices", el("pie chart viewer"), 5));
    await session.step(22, "And pie chart viewer should show 1000 rows", () => showsRows(page, el("pie chart viewer"), 1000));
    await run.scenario("The blank rows get a wedge of their own", async () => {
      await session.step(25, "Then pie chart viewer should have a \"slice (no value)\" area", () => hasArea(page, el("pie chart viewer"), "slice (no value)"));
      await session.step(26, "And the \"angle value of (no value)\" reading of pie chart viewer should be 112", () => readingIs(page, "angle value of (no value)", el("pie chart viewer"), 112));
      await session.step(27, "And the \"share of (no value)\" reading of pie chart viewer should be 11.2", () => readingIs(page, "share of (no value)", el("pie chart viewer"), 11.2));
      await session.step(28, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 793", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 793));
      await session.step(29, "And pie chart viewer should be painted", () => painted(page, el("pie chart viewer")));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking the blank wedge selects exactly the blank rows", async () => {
      await session.step(33, "When user clicks on the \"slice (no value)\" area of pie chart viewer", () => clickArea(page, "slice (no value)", el("pie chart viewer")));
      await session.step(34, "Then 112 rows should be selected", () => selectedRowCount(page, 112));
      await session.step(35, "And no rows where \"RACE_GAPS\" is \"Caucasian\" should be selected", () => noneOfSelected(page, "RACE_GAPS", "Caucasian"));
      await session.step(36, "And no rows where \"RACE_GAPS\" is \"Asian\" should be selected", () => noneOfSelected(page, "RACE_GAPS", "Asian"));
      await session.step(37, "And pie chart viewer should have a \"selected (no value)\" area", () => hasArea(page, el("pie chart viewer"), "selected (no value)"));
      await session.step(38, "And pie chart viewer should show more selection highlight than before", () => moreHighlight(page, el("pie chart viewer")));
      await session.step(39, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(40, "Then no rows should be selected", () => noneSelected(page));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Include Nulls off drops the wedge and the legend item with it", async () => {
      await session.step(44, "Then pie chart viewer should have a \"slice (no value)\" area", () => hasArea(page, el("pie chart viewer"), "slice (no value)"));
      await session.step(45, "And the legend of pie chart viewer should list 5 items", () => legendLists(page, el("pie chart viewer"), 5));
      await session.step(46, "When user sets \"Include Nulls\" property of pie chart viewer to \"false\"", () => setProperty(page, "Include Nulls", el("pie chart viewer"), "false"));
      await session.step(47, "Then pie chart viewer should not have a \"slice (no value)\" area", () => hasNoArea(page, el("pie chart viewer"), "slice (no value)"));
      await session.step(48, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(49, "And the legend of pie chart viewer should list 4 items", () => legendLists(page, el("pie chart viewer"), 4));
      await session.step(50, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 793", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 793));
      await session.step(51, "And the \"share of Caucasian\" reading of pie chart viewer should be between 89.3 and 89.31", () => readingBetween(page, "share of Caucasian", el("pie chart viewer"), 89.3, 89.31));
      await session.step(52, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Include Nulls on brings both back", async () => {
      await session.step(56, "When user sets \"Include Nulls\" property of pie chart viewer to \"true\"", () => setProperty(page, "Include Nulls", el("pie chart viewer"), "true"));
      await session.step(57, "Then pie chart viewer should have a \"slice (no value)\" area", () => hasArea(page, el("pie chart viewer"), "slice (no value)"));
      await session.step(58, "And the \"slices\" reading of pie chart viewer should be 5", () => readingIs(page, "slices", el("pie chart viewer"), 5));
      await session.step(59, "And the legend of pie chart viewer should list 5 items", () => legendLists(page, el("pie chart viewer"), 5));
      await session.step(60, "And the \"share of (no value)\" reading of pie chart viewer should be 11.2", () => readingIs(page, "share of (no value)", el("pie chart viewer"), 11.2));
      await session.step(61, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Show Mouse Over Row Group off a hovered wedge paints no group overlay", async () => {
      await session.step(65, "When user sets \"Show Mouse Over Row Group\" property of pie chart viewer to \"false\"", () => setProperty(page, "Show Mouse Over Row Group", el("pie chart viewer"), "false"));
      await session.step(66, "And user moves the pointer away from pie chart viewer", () => pointerAway(page, el("pie chart viewer")));
      await session.step(67, "Then pie chart viewer should not have a \"mouse over Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "mouse over Caucasian"));
      await session.step(68, "When user hovers over the \"slice Caucasian\" area of pie chart viewer", () => hoverArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(69, "Then pie chart viewer should not have a \"mouse over Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "mouse over Caucasian"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With it on the hovered wedge lights up and goes out again", async () => {
      await session.step(73, "When user moves the pointer away from pie chart viewer", () => pointerAway(page, el("pie chart viewer")));
      await session.step(74, "And user sets \"Show Mouse Over Row Group\" property of pie chart viewer to \"true\"", () => setProperty(page, "Show Mouse Over Row Group", el("pie chart viewer"), "true"));
      await session.step(75, "And user hovers over the \"slice Caucasian\" area of pie chart viewer", () => hoverArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(76, "Then pie chart viewer should have a \"mouse over Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "mouse over Caucasian"));
      await session.step(77, "And pie chart viewer should not have a \"mouse over Asian\" area", () => hasNoArea(page, el("pie chart viewer"), "mouse over Asian"));
      await session.step(78, "And pie chart viewer should have repainted", () => repainted(page, el("pie chart viewer")));
      await session.step(79, "When user hovers over the \"slice Asian\" area of pie chart viewer", () => hoverArea(page, "slice Asian", el("pie chart viewer")));
      await session.step(80, "Then pie chart viewer should have a \"mouse over Asian\" area", () => hasArea(page, el("pie chart viewer"), "mouse over Asian"));
      await session.step(81, "And pie chart viewer should not have a \"mouse over Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "mouse over Caucasian"));
      await session.step(82, "When user moves the pointer away from pie chart viewer", () => pointerAway(page, el("pie chart viewer")));
      await session.step(83, "Then pie chart viewer should not have a \"mouse over Asian\" area", () => hasNoArea(page, el("pie chart viewer"), "mouse over Asian"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Back on the gapless column the blank wedge is gone for good", async () => {
      await session.step(87, "When user sets \"Category\" property of pie chart viewer to \"RACE\"", () => setProperty(page, "Category", el("pie chart viewer"), "RACE"));
      await session.step(88, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(89, "And pie chart viewer should not have a \"slice (no value)\" area", () => hasNoArea(page, el("pie chart viewer"), "slice (no value)"));
      await session.step(90, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 896", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 896));
      await session.step(91, "When user removes \"RACE_GAPS\" column", () => removeColumn(page, "RACE_GAPS"));
      await session.step(92, "Then the table should not have a column \"RACE_GAPS\"", () => hasNoColumn(page, "RACE_GAPS"));
      await session.step(93, "And pie chart viewer should be painted", () => painted(page, el("pie chart viewer")));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
