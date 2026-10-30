/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart, core#974]
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
import {barsStacked, zoomCategories} from '../../../bindings/bar-chart.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterBetween, filterPasses, filterPassesAll, filterTo, removeColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, clickArea, hasArea, hasNoArea, hoverArea, noErrors, propertyShouldBe, readingAtLeast, readingIs, reportsNoError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies, clickEmptySpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart edge cases", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-edge-cases.feature", import.meta.url);
  test("A negative Bar Corner Radius draws the bars and logs no error", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Stack","SEX"]]), [["Split","RACE"],["Stack","SEX"]]);
    await session.step(21, "When user sets \"Bar Corner Radius\" property of bar chart viewer to \"-5\"", () => setProperty(page, "Bar Corner Radius", el("bar chart viewer"), "-5"));
    await session.step(22, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(23, "And the \"bar Caucasian\" area of bar chart viewer should be painted", () => areaPainted(page, "bar Caucasian", el("bar chart viewer")));
    await session.step(24, "And no errors should have been logged", () => noErrors(page));
  });
  test("The tooltip shows a category with markup characters as typed", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "Given user opens a table \"markup-categories\" with:", () => openTableOf(page, "markup-categories", [["Category","Value"],["a<b","1"],["a<b","2"],["x & y","3"]]), [["Category","Value"],["a<b","1"],["a<b","2"],["x & y","3"]]);
    await session.step(33, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Category"]]), [["Split","Category"]]);
    await session.step(35, "Then the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
    await session.step(36, "When user hovers over the \"bar a<b\" area of bar chart viewer", () => hoverArea(page, "bar a<b", el("bar chart viewer")));
    await session.step(37, "Then tooltip should contain text \"a<b\"", () => shouldContainText(page, el("tooltip"), "a<b"));
    await session.step(38, "And tooltip should contain text \"2 rows\"", () => shouldContainText(page, el("tooltip"), "2 rows"));
    await session.step(39, "When user hovers over the \"bar x & y\" area of bar chart viewer", () => hoverArea(page, "bar x & y", el("bar chart viewer")));
    await session.step(40, "Then tooltip should contain text \"x & y\"", () => shouldContainText(page, el("tooltip"), "x & y"));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
  test("A Filter click on a chart stacked by a date column is released after an unrelated setting changes", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(45, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(46, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Stack","STARTED"],["On Click","Filter"]]), [["Split","RACE"],["Stack","STARTED"],["On Click","Filter"]]);
    await session.step(50, "When user clicks on the \"bar Asian\" area of bar chart viewer", () => clickArea(page, "bar Asian", el("bar chart viewer")));
    await session.step(51, "Then 6 rows should pass the filter", () => filterPasses(page, 6));
    await session.step(52, "When user sets \"Bar Corner Radius\" property of bar chart viewer to \"4\"", () => setProperty(page, "Bar Corner Radius", el("bar chart viewer"), "4"));
    await session.step(53, "And user clicks on empty plot space of bar chart viewer", () => clickEmptySpace(page, el("bar chart viewer")));
    await session.step(54, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
  });
  test("With Row Source All and Show Filtered Rows every bar stays when every row is filtered out", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(59, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(60, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Row Source","All"],["Show Filtered Rows","true"]]), [["Split","RACE"],["Row Source","All"],["Show Filtered Rows","true"]]);
    await session.step(64, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(65, "When user filters rows where \"AGE\" is between 200 and 300", () => filterBetween(page, "AGE", 200, 300));
    await session.step(66, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(67, "And the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
  });
  test("\"abc\" typed into the numeric Min property is not taken and the old value stays [category=Y Axis, caption=Min, kept=]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(72, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(73, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
    await session.step(75, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
    await session.step(76, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(77, "Given \"Y Axis\" category in context panel is expanded", () => isExpanded(page, el("\"Y Axis\" category in context panel")));
    await session.step(78, "When user enters \"abc\" into \"Min\" property in context panel", () => enterInto(page, "abc", el("\"Min\" property in context panel")));
    await session.step(79, "Then \"Min\" property of bar chart viewer should be \"\"", () => propertyShouldBe(page, "Min", el("bar chart viewer"), ""));
    await session.step(80, "And the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(81, "And no errors should have been logged", () => noErrors(page));
  });
  test("\"abc\" typed into the numeric Max Bar Height property is not taken and the old value stays [category=Style, caption=Max Bar Height, kept=50]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(72, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(73, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
    await session.step(75, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
    await session.step(76, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(77, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
    await session.step(78, "When user enters \"abc\" into \"Max Bar Height\" property in context panel", () => enterInto(page, "abc", el("\"Max Bar Height\" property in context panel")));
    await session.step(79, "Then \"Max Bar Height\" property of bar chart viewer should be \"50\"", () => propertyShouldBe(page, "Max Bar Height", el("bar chart viewer"), "50"));
    await session.step(80, "And the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(81, "And no errors should have been logged", () => noErrors(page));
  });
  test("\"abc\" typed into the numeric Max Category Width property is not taken and the old value stays [category=Style, caption=Max Category Width, kept=300]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(72, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(73, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
    await session.step(75, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
    await session.step(76, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(77, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
    await session.step(78, "When user enters \"abc\" into \"Max Category Width\" property in context panel", () => enterInto(page, "abc", el("\"Max Category Width\" property in context panel")));
    await session.step(79, "Then \"Max Category Width\" property of bar chart viewer should be \"300\"", () => propertyShouldBe(page, "Max Category Width", el("bar chart viewer"), "300"));
    await session.step(80, "And the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(81, "And no errors should have been logged", () => noErrors(page));
  });
  test("Zoomed to two bars, then filtered to one category, the bar chart resets the zoom and shows the remaining bar", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(91, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(92, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
    await session.step(94, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(95, "When user zooms into the categories from the \"bar Asian\" area to the \"bar Black\" area of bar chart viewer", () => zoomCategories(page, "bar Asian", "bar Black", el("bar chart viewer")));
    await session.step(96, "Then the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
    await session.step(97, "When user filters rows where \"RACE\" is \"Other\"", () => filterTo(page, "RACE", "Other"));
    await session.step(98, "Then the \"bars\" reading of bar chart viewer should be 1", () => readingIs(page, "bars", el("bar chart viewer"), 1));
    await session.step(99, "And bar chart viewer should have a \"bar Other\" area", () => hasArea(page, el("bar chart viewer"), "bar Other"));
    await session.step(100, "When user resets the filter", () => resetFilter(page));
    await session.step(101, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(102, "And no errors should have been logged", () => noErrors(page));
  });
  test("Sorted by value in asc order, the category whose aggregate is empty has a bar only with Show Empty Bars, as the lowest value [order=asc, where=above]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(106, "Given user opens a table \"empty-aggregate\" with:", () => openTableOf(page, "empty-aggregate", [["Category","Value"],["A","1"],["A","3"],["B","5"],["B","7"],["C",""],["C",""]]), [["Category","Value"],["A","1"],["A","3"],["B","5"],["B","7"],["C",""],["C",""]]);
    await session.step(114, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Category"],["Value","Value"],["Value Aggr Type","avg"],["Orientation","horizontal"],["Bar Sort Type","by value"],["Bar Sort Order","asc"]]), [["Split","Category"],["Value","Value"],["Value Aggr Type","avg"],["Orientation","horizontal"],["Bar Sort Type","by value"],["Bar Sort Order","asc"]]);
    await session.step(121, "Then the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
    await session.step(122, "And bar chart viewer should not have a \"bar C\" area", () => hasNoArea(page, el("bar chart viewer"), "bar C"));
    await session.step(123, "When user sets \"showEmptyBars\" property of bar chart viewer to \"true\"", () => setProperty(page, "showEmptyBars", el("bar chart viewer"), "true"));
    await session.step(124, "Then the \"bars\" reading of bar chart viewer should be 3", () => readingIs(page, "bars", el("bar chart viewer"), 3));
    await session.step(125, "And bar chart viewer should have a \"bar C\" area", () => hasArea(page, el("bar chart viewer"), "bar C"));
    await session.step(126, "And the bars of bar chart viewer should lie one under another", () => barsStacked(page, el("bar chart viewer")));
    await session.step(127, "And the \"bar C\" area of bar chart viewer should lie above the \"bar A\" area", () => areaLies(page, "bar C", el("bar chart viewer"), "above", "bar A"));
    await session.step(128, "And the \"bar C\" area of bar chart viewer should lie above the \"bar B\" area", () => areaLies(page, "bar C", el("bar chart viewer"), "above", "bar B"));
    await session.step(129, "And no errors should have been logged", () => noErrors(page));
  });
  test("Sorted by value in desc order, the category whose aggregate is empty has a bar only with Show Empty Bars, as the lowest value [order=desc, where=below]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(106, "Given user opens a table \"empty-aggregate\" with:", () => openTableOf(page, "empty-aggregate", [["Category","Value"],["A","1"],["A","3"],["B","5"],["B","7"],["C",""],["C",""]]), [["Category","Value"],["A","1"],["A","3"],["B","5"],["B","7"],["C",""],["C",""]]);
    await session.step(114, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Category"],["Value","Value"],["Value Aggr Type","avg"],["Orientation","horizontal"],["Bar Sort Type","by value"],["Bar Sort Order","desc"]]), [["Split","Category"],["Value","Value"],["Value Aggr Type","avg"],["Orientation","horizontal"],["Bar Sort Type","by value"],["Bar Sort Order","desc"]]);
    await session.step(121, "Then the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
    await session.step(122, "And bar chart viewer should not have a \"bar C\" area", () => hasNoArea(page, el("bar chart viewer"), "bar C"));
    await session.step(123, "When user sets \"showEmptyBars\" property of bar chart viewer to \"true\"", () => setProperty(page, "showEmptyBars", el("bar chart viewer"), "true"));
    await session.step(124, "Then the \"bars\" reading of bar chart viewer should be 3", () => readingIs(page, "bars", el("bar chart viewer"), 3));
    await session.step(125, "And bar chart viewer should have a \"bar C\" area", () => hasArea(page, el("bar chart viewer"), "bar C"));
    await session.step(126, "And the bars of bar chart viewer should lie one under another", () => barsStacked(page, el("bar chart viewer")));
    await session.step(127, "And the \"bar C\" area of bar chart viewer should lie below the \"bar A\" area", () => areaLies(page, "bar C", el("bar chart viewer"), "below", "bar A"));
    await session.step(128, "And the \"bar C\" area of bar chart viewer should lie below the \"bar B\" area", () => areaLies(page, "bar C", el("bar chart viewer"), "below", "bar B"));
    await session.step(129, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the Stack column stops the stacking and logs no error", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(138, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(139, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Stack","SEX"]]), [["Split","RACE"],["Stack","SEX"]]);
    await session.step(142, "Then the \"stack segments\" reading of bar chart viewer should be at least 1", () => readingAtLeast(page, "stack segments", el("bar chart viewer"), 1));
    await session.step(143, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(144, "Then the \"stack segments\" reading of bar chart viewer should be 0", () => readingIs(page, "stack segments", el("bar chart viewer"), 0));
    await session.step(145, "And the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(146, "And bar chart viewer should report no error", () => reportsNoError(page, el("bar chart viewer")));
    await session.step(147, "And no errors should have been logged", () => noErrors(page));
  });
});
