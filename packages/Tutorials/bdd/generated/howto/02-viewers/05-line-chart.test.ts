/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/05-line-chart.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {pickAggregation, pickLineChartAggregation, viewerShare} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, expand, selectIn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hoverArea, noErrors, propertyShouldBe, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {addLineChartSplit, pickInColumnSelector, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/05-line-chart.feature", import.meta.url);
  test("Configure a line chart", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(24, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"]]);
    await session.step(27, "And line chart viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("line chart viewer"), 65));
    await session.step(28, "When the video says \"A line chart shows how a value changes along X\"", () => videoSays(page, "A line chart shows how a value changes along X"));
    await session.step(30, "And user picks \"STARTED\" in the \"x\" column selector of line chart viewer", () => pickInColumnSelector(page, "STARTED", "x", el("line chart viewer")), undefined, "X can be a number or a date");
    await session.step(31, "Then the \"x column\" reading of line chart viewer should be \"STARTED\"", () => readingReads(page, "x column", el("line chart viewer"), "STARTED"));
    await session.step(33, "When user picks the \"year month\" aggregation in the \"x\" column selector of line chart viewer", () => pickAggregation(page, "year month", "x", el("line chart viewer")), undefined, "Group the dates, here by month");
    await session.step(34, "Then the \"aggregated\" reading of line chart viewer should be \"true\"", () => readingReads(page, "aggregated", el("line chart viewer"), "true"));
    await session.step(35, "And the \"aggregation\" reading of line chart viewer should be \"avg\"", () => readingReads(page, "aggregation", el("line chart viewer"), "avg"));
    await session.step(37, "When user adds \"SEX\" to the splits of line chart viewer", () => addLineChartSplit(page, "SEX", el("line chart viewer")), undefined, "Split: one line for each category");
    await session.step(38, "Then the \"lines\" reading of line chart viewer should be 2", () => readingIs(page, "lines", el("line chart viewer"), 2));
    await session.step(40, "When user hovers over the \"x axis\" area of line chart viewer", () => hoverArea(page, "x axis", el("line chart viewer")), undefined, "The gear opens every setting, grouped in sections");
    await session.step(41, "And user clicks on settings icon of line chart viewer", () => clickOn(page, el("settings icon of line chart viewer")));
    await session.step(42, "Then the context panel should show \"Line chart\"", () => contextPanelShows(page, "Line chart"));
    await session.step(44, "When user expands \"Y Axis\" category in context panel", () => expand(page, el("\"Y Axis\" category in context panel")), undefined, "Y Axis: add more columns, each gets its own chart");
    await session.step(45, "And user clicks on \"...\" button in \"Y\" property in context panel", () => clickOn(page, el("\"...\" button in \"Y\" property in context panel")));
    await session.step(46, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(47, "When user toggles the \"HEIGHT\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "HEIGHT", el("\"Select columns...\" dialog")));
    await session.step(48, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(49, "Then the \"charts\" reading of line chart viewer should be 2", () => readingIs(page, "charts", el("line chart viewer"), 2));
    await session.step(51, "When user checks \"Multi Axis\" property in context panel", () => check(page, el("\"Multi Axis\" property in context panel")), undefined, "Data: Multi Axis draws them in one chart");
    await session.step(52, "Then the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
    await session.step(53, "And the \"y axes\" reading of line chart viewer should be 2", () => readingIs(page, "y axes", el("line chart viewer"), 2));
    await session.step(55, "When user picks \"med\" in the Y Aggr Type selector of line chart viewer", () => pickLineChartAggregation(page, "med", el("line chart viewer")), undefined, "Y Aggr Type: how each month's rows are combined");
    await session.step(56, "Then the \"aggregation\" reading of line chart viewer should be \"med, med\"", () => readingReads(page, "aggregation", el("line chart viewer"), "med, med"));
    await session.step(57, "And \"Aggr Type\" property of line chart viewer should be \"med\"", () => propertyShouldBe(page, "Aggr Type", el("line chart viewer"), "med"));
    await session.step(59, "When user expands \"Lines\" category in context panel", () => expand(page, el("\"Lines\" category in context panel")), undefined, "Lines: smooth the lines");
    await session.step(60, "And user selects \"Spline\" in \"Interpolation\" property in context panel", () => selectIn(page, "Spline", el("\"Interpolation\" property in context panel")));
    await session.step(61, "Then \"Interpolation\" property of line chart viewer should be \"Spline\"", () => propertyShouldBe(page, "Interpolation", el("line chart viewer"), "Spline"));
    await session.step(62, "And no errors should have been logged", () => noErrors(page));
  });
});
