/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-category-order.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart]
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
import {barsStacked} from '../../../bindings/bar-chart.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, menuDoesNotList, menuLists, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingIs, rightClickArea, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart category order", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-category-order.feature", import.meta.url);
  test("Bar chart category order", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
    await session.step(17, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
    await session.step(18, "And the bars of bar chart viewer should lie one under another", () => barsStacked(page, el("bar chart viewer")));
    await run.scenario("Custom Order... orders the bars by row count in both directions, and Default gives the column's order", async () => {
      await session.step(21, "When user picks \"Order > Custom Order...\" from the context menu of the \"bar Caucasian\" area of bar chart viewer", () => pickFromAreaContextMenu(page, "Order > Custom Order...", "bar Caucasian", el("bar chart viewer")));
      await session.step(22, "Then \"Order categories: RACE\" dialog should be visible", () => shouldBe(page, el("\"Order categories: RACE\" dialog"), "visible"));
      await session.step(23, "And \"Bar Sort Type\" property of bar chart viewer should be \"custom\"", () => propertyShouldBe(page, "Bar Sort Type", el("bar chart viewer"), "custom"));
      await session.step(24, "And the \"bar Other\" area of bar chart viewer should lie above the \"bar Caucasian\" area", () => areaLies(page, "bar Other", el("bar chart viewer"), "above", "bar Caucasian"));
      await session.step(25, "And the \"bar Black\" area of bar chart viewer should lie above the \"bar Asian\" area", () => areaLies(page, "bar Black", el("bar chart viewer"), "above", "bar Asian"));
      await session.step(26, "When user selects \"Row count\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Row count", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(27, "And user selects \"Descending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Descending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(28, "Then the \"bar Caucasian\" area of bar chart viewer should lie above the \"bar Other\" area", () => areaLies(page, "bar Caucasian", el("bar chart viewer"), "above", "bar Other"));
      await session.step(29, "And the \"bar Other\" area of bar chart viewer should lie above the \"bar Black\" area", () => areaLies(page, "bar Other", el("bar chart viewer"), "above", "bar Black"));
      await session.step(30, "And the \"bar Black\" area of bar chart viewer should lie above the \"bar Asian\" area", () => areaLies(page, "bar Black", el("bar chart viewer"), "above", "bar Asian"));
      await session.step(31, "When user selects \"Ascending\" in \"Direction\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Ascending", el("\"Direction\" input in \"Order categories: RACE\" dialog")));
      await session.step(32, "Then the \"bar Asian\" area of bar chart viewer should lie above the \"bar Black\" area", () => areaLies(page, "bar Asian", el("bar chart viewer"), "above", "bar Black"));
      await session.step(33, "And the \"bar Other\" area of bar chart viewer should lie above the \"bar Caucasian\" area", () => areaLies(page, "bar Other", el("bar chart viewer"), "above", "bar Caucasian"));
      await session.step(34, "When user selects \"Default\" in \"Order by\" input in \"Order categories: RACE\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: RACE\" dialog")));
      await session.step(35, "Then the \"bar Other\" area of bar chart viewer should lie above the \"bar Caucasian\" area", () => areaLies(page, "bar Other", el("bar chart viewer"), "above", "bar Caucasian"));
      await session.step(36, "And the \"bar Caucasian\" area of bar chart viewer should lie above the \"bar Black\" area", () => areaLies(page, "bar Caucasian", el("bar chart viewer"), "above", "bar Black"));
      await session.step(37, "When user clicks on \"CLOSE\" button in \"Order categories: RACE\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: RACE\" dialog")));
      await session.step(38, "And user sets \"Bar Sort Type\" property of bar chart viewer to \"by category\"", () => setProperty(page, "Bar Sort Type", el("bar chart viewer"), "by category"));
      await session.step(39, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bar Sort Order reverses the column's order, but not an order of the viewer", async () => {
      await session.step(42, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Bar Sort Type","custom"],["Bar Sort Order","asc"]]), [["Bar Sort Type","custom"],["Bar Sort Order","asc"]]);
      await session.step(45, "Then the \"bar Asian\" area of bar chart viewer should lie above the \"bar Black\" area", () => areaLies(page, "bar Asian", el("bar chart viewer"), "above", "bar Black"));
      await session.step(46, "And the \"bar Caucasian\" area of bar chart viewer should lie above the \"bar Other\" area", () => areaLies(page, "bar Caucasian", el("bar chart viewer"), "above", "bar Other"));
      await session.step(47, "When user sets \"categoryOrder\" property of bar chart viewer to '{\"column\":\"AGE\",\"aggr\":\"avg\",\"order\":\"desc\"}'", () => setProperty(page, "categoryOrder", el("bar chart viewer"), "{\"column\":\"AGE\",\"aggr\":\"avg\",\"order\":\"desc\"}"));
      await session.step(48, "Then the \"bar Other\" area of bar chart viewer should lie above the \"bar Black\" area", () => areaLies(page, "bar Other", el("bar chart viewer"), "above", "bar Black"));
      await session.step(49, "And the \"bar Black\" area of bar chart viewer should lie above the \"bar Caucasian\" area", () => areaLies(page, "bar Black", el("bar chart viewer"), "above", "bar Caucasian"));
      await session.step(50, "And the \"bar Caucasian\" area of bar chart viewer should lie above the \"bar Asian\" area", () => areaLies(page, "bar Caucasian", el("bar chart viewer"), "above", "bar Asian"));
      await session.step(51, "When user sets \"Bar Sort Order\" property of bar chart viewer to \"desc\"", () => setProperty(page, "Bar Sort Order", el("bar chart viewer"), "desc"));
      await session.step(52, "Then the \"bar Other\" area of bar chart viewer should lie above the \"bar Black\" area", () => areaLies(page, "bar Other", el("bar chart viewer"), "above", "bar Black"));
      await session.step(53, "And the \"bar Caucasian\" area of bar chart viewer should lie above the \"bar Asian\" area", () => areaLies(page, "bar Caucasian", el("bar chart viewer"), "above", "bar Asian"));
      await session.step(54, "When user sets \"categoryOrder\" property of bar chart viewer to \"\"", () => setProperty(page, "categoryOrder", el("bar chart viewer"), ""));
      await session.step(55, "And user sets \"Bar Sort Type\" property of bar chart viewer to \"by category\"", () => setProperty(page, "Bar Sort Type", el("bar chart viewer"), "by category"));
      await session.step(56, "Then the \"bar Other\" area of bar chart viewer should lie above the \"bar Caucasian\" area", () => areaLies(page, "bar Other", el("bar chart viewer"), "above", "bar Caucasian"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A numerical split offers no Custom Order...", async () => {
      await session.step(60, "When user sets \"Split\" property of bar chart viewer to \"AGE\"", () => setProperty(page, "Split", el("bar chart viewer"), "AGE"));
      await session.step(61, "And user right-clicks on the \"view\" area of bar chart viewer", () => rightClickArea(page, "view", el("bar chart viewer")));
      await session.step(62, "Then the open menu should list \"Order > by value\"", () => menuLists(page, "Order > by value"));
      await session.step(63, "And the open menu should not list \"Order > Custom Order...\"", () => menuDoesNotList(page, "Order > Custom Order..."));
      await session.step(64, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(65, "And user sets \"Split\" property of bar chart viewer to \"RACE\"", () => setProperty(page, "Split", el("bar chart viewer"), "RACE"));
      await session.step(66, "Then the \"bars\" reading of bar chart viewer should be 4", () => readingIs(page, "bars", el("bar chart viewer"), 4));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
