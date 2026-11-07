/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-interaction.feature
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
import {barsDiffer, zoomCategories} from '../../../bindings/bar-chart.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {allOfSelected, clearSelection, colorCategorical, colorCodedCategorically, colorOff, filterIsExactlyCategory, filterPasses, filterPassesAll, noColorCoding, noneSelected, onlyOfAnySelected, onlyOfSelected, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaPainted, clickArea, clickAreaHolding, dragSelectionBetweenAreas, eventFired, hasArea, hasNoArea, listenFor, loadLayout, moreHighlight, noErrors, propertyShouldBe, readingIs, readingLower, repaintedBy, saveLayoutToServer, setProperties, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickEmptySpace, doubleClickEmptySpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart setup and interaction", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-interaction.feature", import.meta.url);
  test("Bar chart setup and interaction", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(17, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Primary Series Name"],["Value","CAST Idea ID"],["Value Aggr Type","count"]]), [["Split","Primary Series Name"],["Value","CAST Idea ID"],["Value Aggr Type","count"]]);
    await session.step(21, "Then the table should have 100 rows", () => rowCount(page, 100));
    await session.step(22, "And the \"bars\" reading of bar chart viewer should be 5", () => readingIs(page, "bars", el("bar chart viewer"), 5));
    await run.scenario("A bar per category", async () => {
      await session.step(25, "Then bar chart viewer should have a \"bar Triazoles\" area", () => hasArea(page, el("bar chart viewer"), "bar Triazoles"));
      await session.step(26, "And bar chart viewer should have a \"bar Pyrrolidines\" area", () => hasArea(page, el("bar chart viewer"), "bar Pyrrolidines"));
      await session.step(27, "And the \"bar Triazoles\" area of bar chart viewer should be painted", () => areaPainted(page, "bar Triazoles", el("bar chart viewer")));
      await session.step(28, "And the \"bar Pyrrolidines\" area of bar chart viewer should be painted", () => areaPainted(page, "bar Pyrrolidines", el("bar chart viewer")));
      await session.step(29, "And the bars of bar chart viewer should differ in length", () => barsDiffer(page, el("bar chart viewer")));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On Click Filter filters the table to the bar's category", async () => {
      await session.step(33, "When user sets \"On Click\" property of bar chart viewer to \"Filter\"", () => setProperty(page, "On Click", el("bar chart viewer"), "Filter"));
      await session.step(34, "Then \"Row Source\" property of bar chart viewer should be \"Filtered\"", () => propertyShouldBe(page, "Row Source", el("bar chart viewer"), "Filtered"));
      await session.step(35, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(36, "When user clicks on the \"bar Triazoles\" area of bar chart viewer", () => clickArea(page, "bar Triazoles", el("bar chart viewer")));
      await session.step(37, "Then the filter should pass exactly the rows where \"Primary Series Name\" is \"Triazoles\"", () => filterIsExactlyCategory(page, "Primary Series Name", "Triazoles"));
      await session.step(38, "And 64 rows should pass the filter", () => filterPasses(page, 64));
      await session.step(39, "And the \"bars\" reading of bar chart viewer should be 1", () => readingIs(page, "bars", el("bar chart viewer"), 1));
      await session.step(40, "When user clicks on empty plot space of bar chart viewer", () => clickEmptySpace(page, el("bar chart viewer")));
      await session.step(41, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(42, "And the \"bars\" reading of bar chart viewer should be 5", () => readingIs(page, "bars", el("bar chart viewer"), 5));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Row Source All the other bars stay and a click on one of them switches the filter", async () => {
      await session.step(46, "When user sets \"Row Source\" property of bar chart viewer to \"All\"", () => setProperty(page, "Row Source", el("bar chart viewer"), "All"));
      await session.step(47, "And user clicks on the \"bar Triazoles\" area of bar chart viewer", () => clickArea(page, "bar Triazoles", el("bar chart viewer")));
      await session.step(48, "Then the filter should pass exactly the rows where \"Primary Series Name\" is \"Triazoles\"", () => filterIsExactlyCategory(page, "Primary Series Name", "Triazoles"));
      await session.step(49, "And the \"bars\" reading of bar chart viewer should be 5", () => readingIs(page, "bars", el("bar chart viewer"), 5));
      await session.step(50, "When user clicks on the \"bar Pyrrolidines\" area of bar chart viewer", () => clickArea(page, "bar Pyrrolidines", el("bar chart viewer")));
      await session.step(51, "Then the filter should pass exactly the rows where \"Primary Series Name\" is \"Pyrrolidines\"", () => filterIsExactlyCategory(page, "Primary Series Name", "Pyrrolidines"));
      await session.step(52, "And 21 rows should pass the filter", () => filterPasses(page, 21));
      await session.step(53, "When user clicks on empty plot space of bar chart viewer", () => clickEmptySpace(page, el("bar chart viewer")));
      await session.step(54, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
      await session.step(56, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Row Source","Filtered"],["On Click","Select"]]), [["Row Source","Filtered"],["On Click","Select"]]);
    });
    await run.scenario("An Alt-drag zooms the categories and a double-click resets the view", async () => {
      await session.step(61, "Given user listens for \"d4-bar-chart-reset-view\" event on bar chart viewer", () => listenFor(page, "d4-bar-chart-reset-view", el("bar chart viewer")));
      await session.step(62, "When user zooms into the categories from the \"bar Pyrrolidines\" area to the \"bar Triazoles\" area of bar chart viewer", () => zoomCategories(page, "bar Pyrrolidines", "bar Triazoles", el("bar chart viewer")));
      await session.step(63, "Then the \"bars\" reading of bar chart viewer should be lower than before", () => readingLower(page, "bars", el("bar chart viewer")));
      await session.step(64, "And bar chart viewer should have a \"bar Triazoles\" area", () => hasArea(page, el("bar chart viewer"), "bar Triazoles"));
      await session.step(65, "And bar chart viewer should not have a \"bar Aminopiperidines\" area", () => hasNoArea(page, el("bar chart viewer"), "bar Aminopiperidines"));
      await session.step(66, "When user double-clicks on empty plot space of bar chart viewer", () => doubleClickEmptySpace(page, el("bar chart viewer")));
      await session.step(67, "Then \"d4-bar-chart-reset-view\" event should have fired on bar chart viewer", () => eventFired(page, "d4-bar-chart-reset-view", el("bar chart viewer")));
      await session.step(68, "And the \"bars\" reading of bar chart viewer should be 5", () => readingIs(page, "bars", el("bar chart viewer"), 5));
      await session.step(69, "And bar chart viewer should have a \"bar Aminopiperidines\" area", () => hasArea(page, el("bar chart viewer"), "bar Aminopiperidines"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On Click Select selects the category without filtering", async () => {
      await session.step(73, "Then \"On Click\" property of bar chart viewer should be \"Select\"", () => propertyShouldBe(page, "On Click", el("bar chart viewer"), "Select"));
      await session.step(74, "When user clears the row selection", () => clearSelection(page));
      await session.step(75, "And user clicks on the \"bar Triazoles\" area of bar chart viewer", () => clickArea(page, "bar Triazoles", el("bar chart viewer")));
      await session.step(76, "Then only rows where \"Primary Series Name\" is \"Triazoles\" should be selected", () => onlyOfSelected(page, "Primary Series Name", "Triazoles"));
      await session.step(77, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(78, "And bar chart viewer should show more selection highlight than before", () => moreHighlight(page, el("bar chart viewer")));
      await session.step(79, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(80, "Then no rows should be selected", () => noneSelected(page));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control adds a category to the selection", async () => {
      await session.step(84, "When user clicks on the \"bar Triazoles\" area of bar chart viewer holding Control", () => clickAreaHolding(page, "bar Triazoles", el("bar chart viewer"), "Control"));
      await session.step(85, "Then only rows where \"Primary Series Name\" is \"Triazoles\" should be selected", () => onlyOfSelected(page, "Primary Series Name", "Triazoles"));
      await session.step(86, "When user clicks on the \"bar Pyrrolidines\" area of bar chart viewer holding Control", () => clickAreaHolding(page, "bar Pyrrolidines", el("bar chart viewer"), "Control"));
      await session.step(87, "Then only rows where \"Primary Series Name\" is one of \"Triazoles, Pyrrolidines\" should be selected", () => onlyOfAnySelected(page, "Primary Series Name", "Triazoles, Pyrrolidines"));
      await session.step(88, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(89, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(90, "Then no rows should be selected", () => noneSelected(page));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Shift-drag selects the bars it covers", async () => {
      await session.step(94, "When user drags a selection box from the \"bar Triazoles\" area to the \"bar Pyrrolidines\" area of bar chart viewer", () => dragSelectionBetweenAreas(page, "bar Triazoles", "bar Pyrrolidines", el("bar chart viewer")));
      await session.step(95, "Then all rows where \"Primary Series Name\" is \"Triazoles\" should be selected", () => allOfSelected(page, "Primary Series Name", "Triazoles"));
      await session.step(96, "And all rows where \"Primary Series Name\" is \"Pyrrolidines\" should be selected", () => allOfSelected(page, "Primary Series Name", "Pyrrolidines"));
      await session.step(97, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(98, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(99, "Then no rows should be selected", () => noneSelected(page));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Split column's color coding drives the bar colors and survives a layout round-trip", async () => {
      await session.step(103, "Then \"Primary Series Name\" column should have no color coding", () => noColorCoding(page, "Primary Series Name"));
      await session.step(104, "When user colors \"Primary Series Name\" column categorically:", () => colorCategorical(page, "Primary Series Name", [["Triazoles","#FF0000"],["Pyrrolidines","#0000FF"]]), [["Triazoles","#FF0000"],["Pyrrolidines","#0000FF"]]);
      await session.step(107, "Then \"Primary Series Name\" column should be color-coded categorically", () => colorCodedCategorically(page, "Primary Series Name"));
      await session.step(108, "And bar chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("bar chart viewer"), 500));
      await session.step(109, "And the \"bar Triazoles\" area of bar chart viewer should contain the color \"#FF0000\"", () => areaColor(page, "bar Triazoles", el("bar chart viewer"), "#FF0000"));
      await session.step(110, "And the \"bar Pyrrolidines\" area of bar chart viewer should contain the color \"#0000FF\"", () => areaColor(page, "bar Pyrrolidines", el("bar chart viewer"), "#0000FF"));
      await session.step(111, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(112, "And user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(113, "Then bar chart viewer should be absent", () => shouldBe(page, el("bar chart viewer"), "absent"));
      await session.step(114, "When user removes the coloring of \"Primary Series Name\" column", () => colorOff(page, "Primary Series Name"));
      await session.step(115, "Then \"Primary Series Name\" column should have no color coding", () => noColorCoding(page, "Primary Series Name"));
      await session.step(116, "When user loads the saved layout", () => loadLayout(page));
      await session.step(117, "Then bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
      await session.step(118, "And \"Primary Series Name\" column should be color-coded categorically", () => colorCodedCategorically(page, "Primary Series Name"));
      await session.step(119, "And the \"bar Triazoles\" area of bar chart viewer should contain the color \"#FF0000\"", () => areaColor(page, "bar Triazoles", el("bar chart viewer"), "#FF0000"));
      await session.step(120, "And the \"bar Pyrrolidines\" area of bar chart viewer should contain the color \"#0000FF\"", () => areaColor(page, "bar Pyrrolidines", el("bar chart viewer"), "#0000FF"));
      await session.step(121, "When user removes the coloring of \"Primary Series Name\" column", () => colorOff(page, "Primary Series Name"));
      await session.step(122, "Then \"Primary Series Name\" column should have no color coding", () => noColorCoding(page, "Primary Series Name"));
      await session.step(123, "And bar chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("bar chart viewer"), 500));
      await session.step(124, "And the \"bar Triazoles\" area of bar chart viewer should contain the color \"#96D794\"", () => areaColor(page, "bar Triazoles", el("bar chart viewer"), "#96D794"));
      await session.step(125, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the bar chart gives back the rows its Filter click filtered out", async () => {
      await session.step(129, "When user sets \"On Click\" property of bar chart viewer to \"Filter\"", () => setProperty(page, "On Click", el("bar chart viewer"), "Filter"));
      await session.step(130, "And user clicks on the \"bar Triazoles\" area of bar chart viewer", () => clickArea(page, "bar Triazoles", el("bar chart viewer")));
      await session.step(131, "Then 64 rows should pass the filter", () => filterPasses(page, 64));
      await session.step(132, "When user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(133, "Then the open tableview should have 0 bar chart viewers", () => viewerCount(page, 0, "bar chart"));
      await session.step(134, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
