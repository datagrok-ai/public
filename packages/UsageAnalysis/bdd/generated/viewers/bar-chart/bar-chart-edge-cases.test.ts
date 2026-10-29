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
import {deleteSelected, filterBetween, filterPasses, filterPassesAll, filterTo, resetFilter, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, clickArea, hasArea, hoverArea, noErrors, propertyShouldBe, readingIs, reportsNoError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies, clickEmptySpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart edge cases", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-edge-cases.feature", import.meta.url);
  test("A vertical bar chart over a table with every row deleted draws no bars and logs no error", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Orientation","vertical"]]), [["Split","RACE"],["Orientation","vertical"]]);
    await session.step(21, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(22, "When user selects all rows", () => selectAllRows(page));
    await session.step(23, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(24, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(25, "And the \"bars\" reading of bar chart viewer should be 0", () => readingIs(page, "bars", el("bar chart viewer"), 0));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
  });
  test("A negative Bar Corner Radius draws the bars and logs no error", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(31, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Stack","SEX"]]), [["Split","RACE"],["Stack","SEX"]]);
    await session.step(34, "When user sets \"Bar Corner Radius\" property of bar chart viewer to \"-5\"", () => setProperty(page, "Bar Corner Radius", el("bar chart viewer"), "-5"));
    await session.step(35, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(36, "And the \"bar Caucasian\" area of bar chart viewer should be painted", () => areaPainted(page, "bar Caucasian", el("bar chart viewer")));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("The tooltip shows a category with markup characters as typed", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "Given user opens a table \"markup-categories\" with:", () => openTableOf(page, "markup-categories", [["Category","Value"],["a<b","1"],["a<b","2"],["x & y","3"]]), [["Category","Value"],["a<b","1"],["a<b","2"],["x & y","3"]]);
    await session.step(46, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Category"]]), [["Split","Category"]]);
    await session.step(48, "Then the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
    await session.step(49, "When user hovers over the \"bar a<b\" area of bar chart viewer", () => hoverArea(page, "bar a<b", el("bar chart viewer")));
    await session.step(50, "Then tooltip should contain text \"a<b\"", () => shouldContainText(page, el("tooltip"), "a<b"));
    await session.step(51, "And tooltip should contain text \"2 rows\"", () => shouldContainText(page, el("tooltip"), "2 rows"));
    await session.step(52, "When user hovers over the \"bar x & y\" area of bar chart viewer", () => hoverArea(page, "bar x & y", el("bar chart viewer")));
    await session.step(53, "Then tooltip should contain text \"x & y\"", () => shouldContainText(page, el("tooltip"), "x & y"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
  });
  test("A Filter click on a chart stacked by a date column is released after an unrelated setting changes", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(58, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(59, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Stack","STARTED"],["On Click","Filter"]]), [["Split","RACE"],["Stack","STARTED"],["On Click","Filter"]]);
    await session.step(63, "When user clicks on the \"bar Asian\" area of bar chart viewer", () => clickArea(page, "bar Asian", el("bar chart viewer")));
    await session.step(64, "Then 6 rows should pass the filter", () => filterPasses(page, 6));
    await session.step(65, "When user sets \"Bar Corner Radius\" property of bar chart viewer to \"4\"", () => setProperty(page, "Bar Corner Radius", el("bar chart viewer"), "4"));
    await session.step(66, "And user clicks on empty plot space of bar chart viewer", () => clickEmptySpace(page, el("bar chart viewer")));
    await session.step(67, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
  });
  test("With Row Source All and Show Filtered Rows every bar stays when every row is filtered out", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@realizes:core#974"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(72, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(73, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"],["Row Source","All"],["Show Filtered Rows","true"]]), [["Split","RACE"],["Row Source","All"],["Show Filtered Rows","true"]]);
    await session.step(77, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(78, "When user filters rows where \"AGE\" is between 200 and 300", () => filterBetween(page, "AGE", 200, 300));
    await session.step(79, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(80, "And the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(81, "And no errors should have been logged", () => noErrors(page));
  });
  test("Max Bar Height typed as 0 in the context panel is not taken by the bar chart [caption=Max Bar Height, typed=0, kept=50]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:92", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(85, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(86, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
      await session.step(88, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
      await session.step(89, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(90, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(91, "When user enters \"0\" into \"Max Bar Height\" property in context panel", () => enterInto(page, "0", el("\"Max Bar Height\" property in context panel")));
      await session.step(92, "Then \"Max Bar Height\" property of bar chart viewer should be \"50\"", () => propertyShouldBe(page, "Max Bar Height", el("bar chart viewer"), "50"));
      await session.step(93, "And bar chart viewer should report no error", () => reportsNoError(page, el("bar chart viewer")));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    }, 92);
  });
  test("Max Bar Height typed as -5 in the context panel is not taken by the bar chart [caption=Max Bar Height, typed=-5, kept=50]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:92", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(85, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(86, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
      await session.step(88, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
      await session.step(89, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(90, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(91, "When user enters \"-5\" into \"Max Bar Height\" property in context panel", () => enterInto(page, "-5", el("\"Max Bar Height\" property in context panel")));
      await session.step(92, "Then \"Max Bar Height\" property of bar chart viewer should be \"50\"", () => propertyShouldBe(page, "Max Bar Height", el("bar chart viewer"), "50"));
      await session.step(93, "And bar chart viewer should report no error", () => reportsNoError(page, el("bar chart viewer")));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    }, 92);
  });
  test("Max Category Width typed as 0 in the context panel is not taken by the bar chart [caption=Max Category Width, typed=0, kept=300]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:92", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(85, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(86, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
      await session.step(88, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
      await session.step(89, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(90, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(91, "When user enters \"0\" into \"Max Category Width\" property in context panel", () => enterInto(page, "0", el("\"Max Category Width\" property in context panel")));
      await session.step(92, "Then \"Max Category Width\" property of bar chart viewer should be \"300\"", () => propertyShouldBe(page, "Max Category Width", el("bar chart viewer"), "300"));
      await session.step(93, "And bar chart viewer should report no error", () => reportsNoError(page, el("bar chart viewer")));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    }, 92);
  });
  test("Category Value Width typed as 0 in the context panel is not taken by the bar chart [caption=Category Value Width, typed=0, kept=50]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:92", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(85, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(86, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
      await session.step(88, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
      await session.step(89, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(90, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(91, "When user enters \"0\" into \"Category Value Width\" property in context panel", () => enterInto(page, "0", el("\"Category Value Width\" property in context panel")));
      await session.step(92, "Then \"Category Value Width\" property of bar chart viewer should be \"50\"", () => propertyShouldBe(page, "Category Value Width", el("bar chart viewer"), "50"));
      await session.step(93, "And bar chart viewer should report no error", () => reportsNoError(page, el("bar chart viewer")));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    }, 92);
  });
  test("\"abc\" typed into the numeric Min property is not taken and the old value stays [category=Y Axis, caption=Min, kept=]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(106, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(107, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
    await session.step(109, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
    await session.step(110, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(111, "Given \"Y Axis\" category in context panel is expanded", () => isExpanded(page, el("\"Y Axis\" category in context panel")));
    await session.step(112, "When user enters \"abc\" into \"Min\" property in context panel", () => enterInto(page, "abc", el("\"Min\" property in context panel")));
    await session.step(113, "Then \"Min\" property of bar chart viewer should be \"\"", () => propertyShouldBe(page, "Min", el("bar chart viewer"), ""));
    await session.step(114, "And the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(115, "And no errors should have been logged", () => noErrors(page));
  });
  test("\"abc\" typed into the numeric Max Bar Height property is not taken and the old value stays [category=Style, caption=Max Bar Height, kept=50]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(106, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(107, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
    await session.step(109, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
    await session.step(110, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(111, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
    await session.step(112, "When user enters \"abc\" into \"Max Bar Height\" property in context panel", () => enterInto(page, "abc", el("\"Max Bar Height\" property in context panel")));
    await session.step(113, "Then \"Max Bar Height\" property of bar chart viewer should be \"50\"", () => propertyShouldBe(page, "Max Bar Height", el("bar chart viewer"), "50"));
    await session.step(114, "And the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(115, "And no errors should have been logged", () => noErrors(page));
  });
  test("\"abc\" typed into the numeric Max Category Width property is not taken and the old value stays [category=Style, caption=Max Category Width, kept=300]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(106, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(107, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
    await session.step(109, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
    await session.step(110, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(111, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
    await session.step(112, "When user enters \"abc\" into \"Max Category Width\" property in context panel", () => enterInto(page, "abc", el("\"Max Category Width\" property in context panel")));
    await session.step(113, "Then \"Max Category Width\" property of bar chart viewer should be \"300\"", () => propertyShouldBe(page, "Max Category Width", el("bar chart viewer"), "300"));
    await session.step(114, "And the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(115, "And no errors should have been logged", () => noErrors(page));
  });
  test("Zoomed to two bars, then filtered to one category, the bar chart resets the zoom and shows the remaining bar", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(125, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(126, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
    await session.step(128, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(129, "When user zooms into the categories from the \"bar Asian\" area to the \"bar Black\" area of bar chart viewer", () => zoomCategories(page, "bar Asian", "bar Black", el("bar chart viewer")));
    await session.step(130, "Then the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
    await session.step(131, "When user filters rows where \"RACE\" is \"Other\"", () => filterTo(page, "RACE", "Other"));
    await session.step(132, "Then the \"bars\" reading of bar chart viewer should be 1", () => readingIs(page, "bars", el("bar chart viewer"), 1));
    await session.step(133, "And bar chart viewer should have a \"bar Other\" area", () => hasArea(page, el("bar chart viewer"), "bar Other"));
    await session.step(134, "When user resets the filter", () => resetFilter(page));
    await session.step(135, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(136, "And no errors should have been logged", () => noErrors(page));
  });
  test("Sorted by value in asc order, the bar whose aggregate is empty comes last [order=asc]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:156", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(140, "Given user opens a table \"empty-aggregate\" with:", () => openTableOf(page, "empty-aggregate", [["Category","Value"],["A","1"],["A","3"],["B","5"],["B","7"],["C",""],["C",""]]), [["Category","Value"],["A","1"],["A","3"],["B","5"],["B","7"],["C",""],["C",""]]);
      await session.step(148, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Category"],["Value","Value"],["Value Aggr Type","avg"],["Orientation","horizontal"],["Bar Sort Type","by value"],["Bar Sort Order","asc"]]), [["Split","Category"],["Value","Value"],["Value Aggr Type","avg"],["Orientation","horizontal"],["Bar Sort Type","by value"],["Bar Sort Order","asc"]]);
      await session.step(155, "Then the bars of bar chart viewer should lie one under another", () => barsStacked(page, el("bar chart viewer")));
      await session.step(156, "And the \"bar C\" area of bar chart viewer should lie below the \"bar A\" area", () => areaLies(page, "bar C", el("bar chart viewer"), "below", "bar A"));
      await session.step(157, "And the \"bar C\" area of bar chart viewer should lie below the \"bar B\" area", () => areaLies(page, "bar C", el("bar chart viewer"), "below", "bar B"));
      await session.step(158, "And no errors should have been logged", () => noErrors(page));
    }, 156);
  });
  test("Sorted by value in desc order, the bar whose aggregate is empty comes last [order=desc]", {tag: ["@viewers", "@realizes:viewers.bar-chart", "@edge-cases", "@known-failure:156", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(140, "Given user opens a table \"empty-aggregate\" with:", () => openTableOf(page, "empty-aggregate", [["Category","Value"],["A","1"],["A","3"],["B","5"],["B","7"],["C",""],["C",""]]), [["Category","Value"],["A","1"],["A","3"],["B","5"],["B","7"],["C",""],["C",""]]);
      await session.step(148, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Category"],["Value","Value"],["Value Aggr Type","avg"],["Orientation","horizontal"],["Bar Sort Type","by value"],["Bar Sort Order","desc"]]), [["Split","Category"],["Value","Value"],["Value Aggr Type","avg"],["Orientation","horizontal"],["Bar Sort Type","by value"],["Bar Sort Order","desc"]]);
      await session.step(155, "Then the bars of bar chart viewer should lie one under another", () => barsStacked(page, el("bar chart viewer")));
      await session.step(156, "And the \"bar C\" area of bar chart viewer should lie below the \"bar A\" area", () => areaLies(page, "bar C", el("bar chart viewer"), "below", "bar A"));
      await session.step(157, "And the \"bar C\" area of bar chart viewer should lie below the \"bar B\" area", () => areaLies(page, "bar C", el("bar chart viewer"), "below", "bar B"));
      await session.step(158, "And no errors should have been logged", () => noErrors(page));
    }, 156);
  });
});
