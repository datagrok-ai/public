/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/grid-visible-columns-tooltip.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, isExpanded, selectIn, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {mouseOverRowIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dragAreaBy, hasArea, hasNoArea, hoverArea, noErrors, pointerAway, propertyShouldBe, readingHigherThanRemembered, readingReads, rememberReading, tooltipColumns, tooltipValue} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The grid's tooltip and the columns it cannot show", () => {
  const session = feature(test, "features/viewers/tooltips/grid-visible-columns-tooltip.feature", import.meta.url);
  test("The grid's tooltip and the columns it cannot show", {tag: ["@journey", "@viewers", "@realizes:viewers.tooltips"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens a table \"energy_uk\" with:", () => openTableOf(page, "energy_uk", [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]), [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]);
    await run.scenario("By default the grid shows no row tooltip, even for the columns pushed off it", async () => {
      await session.step(31, "Then \"Show Tooltip\" property of grid should be \"show custom tooltip\"", () => propertyShouldBe(page, "Show Tooltip", el("grid"), "show custom tooltip"));
      await session.step(32, "And \"Row Tooltip\" property of grid should be \"\"", () => propertyShouldBe(page, "Row Tooltip", el("grid"), ""));
      await session.step(33, "And \"Show Visible Columns In Tooltip\" property of grid should be \"false\"", () => propertyShouldBe(page, "Show Visible Columns In Tooltip", el("grid"), "false"));
      await session.step(34, "When user drags the \"column resizer value\" area of grid by 1000 pixels to the right", () => dragAreaBy(page, "column resizer value", el("grid"), 1000, "right"));
      await session.step(35, "Then grid should not have a \"cell 1 of target\" area", () => hasNoArea(page, el("grid"), "cell 1 of target"));
      await session.step(36, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(37, "Then the mouse-over row of the table should be 1", () => mouseOverRowIs(page, 1));
      await session.step(38, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(39, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(40, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The grid takes the table's tooltip, and Show Visible Columns In Tooltip is off", async () => {
      await session.step(43, "When user closes all views", () => closeAllViews(page));
      await session.step(44, "And user opens a table \"energy_uk\" with:", () => openTableOf(page, "energy_uk", [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]), [["value","source","target"],["124.72899627685547","Agricultural 'waste'","Bio-conversion"],["0.597000002861023","Bio-conversion","Liquid"],["26.86199951171875","Bio-conversion","Losses"],["280.3219909667969","Bio-conversion","Solid"],["81.14399719238281","Bio-conversion","Gas"]]);
      await session.step(51, "Then \"Show Tooltip\" property of grid should be \"show custom tooltip\"", () => propertyShouldBe(page, "Show Tooltip", el("grid"), "show custom tooltip"));
      await session.step(52, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(53, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(54, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(55, "Given \"Tooltip\" category in context panel is expanded", () => isExpanded(page, el("\"Tooltip\" category in context panel")));
      await session.step(56, "When user selects \"inherit from table\" in \"Show Tooltip\" property in context panel", () => selectIn(page, "inherit from table", el("\"Show Tooltip\" property in context panel")));
      await session.step(57, "Then \"Show Tooltip\" property of grid should be \"inherit from table\"", () => propertyShouldBe(page, "Show Tooltip", el("grid"), "inherit from table"));
      await session.step(58, "When user selects \"Always\" in \"Show Column Names\" property in context panel", () => selectIn(page, "Always", el("\"Show Column Names\" property in context panel")));
      await session.step(59, "Then \"Show Column Names\" property of grid should be \"Always\"", () => propertyShouldBe(page, "Show Column Names", el("grid"), "Always"));
      await session.step(60, "And \"Show Visible Columns In Tooltip\" property in context panel should be unchecked", () => shouldBe(page, el("\"Show Visible Columns In Tooltip\" property in context panel"), "unchecked"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Off, no tooltip while every column fits", async () => {
      await session.step(64, "Then the \"column order\" reading of grid should be \"value, source, target\"", () => readingReads(page, "column order", el("grid"), "value, source, target"));
      await session.step(65, "And grid should have a \"cell 1 of target\" area", () => hasArea(page, el("grid"), "cell 1 of target"));
      await session.step(66, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(67, "Then the mouse-over row of the table should be 1", () => mouseOverRowIs(page, 1));
      await session.step(68, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(69, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(70, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Switched on, the tooltip lists every column while every column fits", async () => {
      await session.step(73, "When user checks \"Show Visible Columns In Tooltip\" property in context panel", () => check(page, el("\"Show Visible Columns In Tooltip\" property in context panel")));
      await session.step(74, "Then \"Show Visible Columns In Tooltip\" property of grid should be \"true\"", () => propertyShouldBe(page, "Show Visible Columns In Tooltip", el("grid"), "true"));
      await session.step(75, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(76, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(77, "And the tooltip should show columns \"value, source, target\"", () => tooltipColumns(page, "value, source, target"));
      await session.step(78, "And the tooltip should show \"source\" as \"Agricultural 'waste'\"", () => tooltipValue(page, "source", "Agricultural 'waste'"));
      await session.step(79, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(80, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Switched on, the tooltip stays the same with columns pushed off the grid", async () => {
      await session.step(83, "When user remembers the \"column width of value\" reading of grid", () => rememberReading(page, "column width of value", el("grid")));
      await session.step(84, "And user drags the \"column resizer value\" area of grid by 1000 pixels to the right", () => dragAreaBy(page, "column resizer value", el("grid"), 1000, "right"));
      await session.step(85, "Then the \"column width of value\" reading of grid should be higher than remembered", () => readingHigherThanRemembered(page, "column width of value", el("grid")));
      await session.step(86, "And grid should not have a \"cell 1 of source\" area", () => hasNoArea(page, el("grid"), "cell 1 of source"));
      await session.step(87, "And grid should not have a \"cell 1 of target\" area", () => hasNoArea(page, el("grid"), "cell 1 of target"));
      await session.step(88, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(89, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(90, "And the tooltip should show columns \"value, source, target\"", () => tooltipColumns(page, "value, source, target"));
      await session.step(91, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(92, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Off again, the tooltip lists only the columns pushed off the grid", async () => {
      await session.step(95, "When user unchecks \"Show Visible Columns In Tooltip\" property in context panel", () => uncheck(page, el("\"Show Visible Columns In Tooltip\" property in context panel")));
      await session.step(96, "Then \"Show Visible Columns In Tooltip\" property of grid should be \"false\"", () => propertyShouldBe(page, "Show Visible Columns In Tooltip", el("grid"), "false"));
      await session.step(97, "When user hovers over the \"cell 1 of value\" area of grid", () => hoverArea(page, "cell 1 of value", el("grid")));
      await session.step(98, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(99, "And the tooltip should show columns \"source, target\"", () => tooltipColumns(page, "source, target"));
      await session.step(100, "And the tooltip should show \"target\" as \"Bio-conversion\"", () => tooltipValue(page, "target", "Bio-conversion"));
      await session.step(101, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(102, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
