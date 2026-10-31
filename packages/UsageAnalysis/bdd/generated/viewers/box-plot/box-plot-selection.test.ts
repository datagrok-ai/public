/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-selection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot, public#2764]
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
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, clearSelection, currentRowValue, deleteSelected, filterOut, filterPasses, filterPassesAll, hasCurrentRow, noRowsWhere, noneOfSelected, noneSelected, onlyOfAnySelected, onlyOfSelected, resetFilter, selectedRowCount, someOfSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaPainted, clickArea, clickAreaHolding, dragSelectionOverArea, eventFired, eventNotFired, hasNoArea, hoverArea, lessHighlight, listenFor, moreHighlight, noBalloons, noErrors, noHighlight, notRepainted, painted, pointerAway, propertyShouldBe, repainted, repaintedBy, setProperties, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickEmptySpace, doubleClickEmptySpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot selection and highlight", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-selection.feature", import.meta.url);
  test("Box plot selection and highlight", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot", "@realizes:public#2764"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"],["Marker Size","10"]]), [["Value","AGE"],["Category 1","RACE"],["Marker Size","10"]]);
    await run.scenario("A marker click sets the current row", async () => {
      await session.step(21, "Given user listens for \"d4-boxplot-point-click\" event on box plot viewer", () => listenFor(page, "d4-boxplot-point-click", el("box plot viewer")));
      await session.step(22, "When user clicks on the \"marker\" area of box plot viewer", () => clickArea(page, "marker", el("box plot viewer")));
      await session.step(23, "Then \"d4-boxplot-point-click\" event should have fired on box plot viewer", () => eventFired(page, "d4-boxplot-point-click", el("box plot viewer")));
      await session.step(24, "And the table should have a current row", () => hasCurrentRow(page));
      await session.step(25, "And \"RACE\" of the current row should be \"Asian\"", () => currentRowValue(page, "RACE", "Asian"));
    });
    await run.scenario("A Shift-drag selects a band and highlights it", async () => {
      await session.step(28, "When user clears the row selection", () => clearSelection(page));
      await session.step(29, "And user drags a selection box over the \"Caucasian values\" area of box plot viewer", () => dragSelectionOverArea(page, "Caucasian values", el("box plot viewer")));
      await session.step(30, "Then some rows where \"RACE\" is \"Caucasian\" should be selected", () => someOfSelected(page, "RACE", "Caucasian"));
      await session.step(31, "And box plot viewer should show more selection highlight than before", () => moreHighlight(page, el("box plot viewer")));
    });
    await run.scenario("The highlight survives a categorical coloring", async () => {
      await session.step(34, "When user clears the row selection", () => clearSelection(page));
      await session.step(35, "And user sets \"Marker Color Column\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "SEX"));
      await session.step(36, "And user drags a selection box over the \"Caucasian values\" area of box plot viewer", () => dragSelectionOverArea(page, "Caucasian values", el("box plot viewer")));
      await session.step(37, "Then some rows where \"RACE\" is \"Caucasian\" should be selected", () => someOfSelected(page, "RACE", "Caucasian"));
      await session.step(38, "And box plot viewer should show more selection highlight than before", () => moreHighlight(page, el("box plot viewer")));
      await session.step(39, "When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("Category labels select categories", async () => {
      await session.step(42, "When user clears the row selection", () => clearSelection(page));
      await session.step(43, "And user clicks on the \"category Caucasian\" area of box plot viewer", () => clickArea(page, "category Caucasian", el("box plot viewer")));
      await session.step(44, "Then only rows where \"RACE\" is \"Caucasian\" should be selected", () => onlyOfSelected(page, "RACE", "Caucasian"));
      await session.step(45, "When user clicks on the \"category Black\" area of box plot viewer holding Control", () => clickAreaHolding(page, "category Black", el("box plot viewer"), "Control"));
      await session.step(46, "Then only rows where \"RACE\" is one of \"Black, Caucasian\" should be selected", () => onlyOfAnySelected(page, "RACE", "Black, Caucasian"));
      await session.step(47, "When user clicks on the \"category Asian\" area of box plot viewer", () => clickArea(page, "category Asian", el("box plot viewer")));
      await session.step(48, "Then only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(49, "Given user listens for \"d4-boxplot-reset-view\" event on box plot viewer", () => listenFor(page, "d4-boxplot-reset-view", el("box plot viewer")));
      await session.step(50, "When user clicks on empty plot space of box plot viewer", () => clickEmptySpace(page, el("box plot viewer")));
      await session.step(51, "Then no rows should be selected", () => noneSelected(page));
      await session.step(52, "And \"d4-boxplot-reset-view\" event should not have fired on box plot viewer", () => eventNotFired(page, "d4-boxplot-reset-view", el("box plot viewer")));
      await session.step(53, "When user double-clicks on empty plot space of box plot viewer", () => doubleClickEmptySpace(page, el("box plot viewer")));
      await session.step(54, "Then \"d4-boxplot-reset-view\" event should have fired on box plot viewer", () => eventFired(page, "d4-boxplot-reset-view", el("box plot viewer")));
    });
    await run.scenario("No selection leaks into a filtered-out category", async () => {
      await session.step(57, "When user clears the row selection", () => clearSelection(page));
      await session.step(58, "And user filters out rows where \"RACE\" is \"Asian\"", () => filterOut(page, "RACE", "Asian"));
      await session.step(59, "And user drags a selection box over the \"view\" area of box plot viewer", () => dragSelectionOverArea(page, "view", el("box plot viewer")));
      await session.step(60, "Then some rows should be selected", () => someSelected(page));
      await session.step(61, "And no rows where \"RACE\" is \"Asian\" should be selected", () => noneOfSelected(page, "RACE", "Asian"));
      await session.step(62, "When user resets the filter", () => resetFilter(page));
      await session.step(63, "And user clears the row selection", () => clearSelection(page));
    });
    await run.scenario("A Shift-drag picks its category when the categories to its left are filtered out (public#2764)", async () => {
      await session.step(67, "When user clears the row selection", () => clearSelection(page));
      await session.step(68, "And user adds a categorical filter on \"RACE\" keeping \"Caucasian, Other\"", () => addCategoricalFilter(page, "RACE", "Caucasian, Other"));
      await session.step(69, "Then 958 rows should pass the filter", () => filterPasses(page, 958));
      await session.step(70, "And box plot viewer should not have a \"category Asian\" area", () => hasNoArea(page, el("box plot viewer"), "category Asian"));
      await session.step(71, "When user drags a selection box over the \"Other values\" area of box plot viewer", () => dragSelectionOverArea(page, "Other values", el("box plot viewer")));
      await session.step(72, "Then some rows where \"RACE\" is \"Other\" should be selected", () => someOfSelected(page, "RACE", "Other"));
      await session.step(73, "And no rows where \"RACE\" is \"Caucasian\" should be selected", () => noneOfSelected(page, "RACE", "Caucasian"));
      await session.step(74, "And 51 rows should be selected", () => selectedRowCount(page, 51));
      await session.step(75, "When user clears the row selection", () => clearSelection(page));
      await session.step(76, "And user drags a selection box over the \"Caucasian values\" area of box plot viewer", () => dragSelectionOverArea(page, "Caucasian values", el("box plot viewer")));
      await session.step(77, "Then some rows where \"RACE\" is \"Caucasian\" should be selected", () => someOfSelected(page, "RACE", "Caucasian"));
      await session.step(78, "And no rows where \"RACE\" is \"Other\" should be selected", () => noneOfSelected(page, "RACE", "Other"));
      await session.step(79, "And 768 rows should be selected", () => selectedRowCount(page, 768));
      await session.step(80, "When user adds a categorical filter on \"RACE\" keeping \"Asian, Black, Caucasian, Other\"", () => addCategoricalFilter(page, "RACE", "Asian, Black, Caucasian, Other"));
      await session.step(81, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(82, "When user clears the row selection", () => clearSelection(page));
    });
    await run.scenario("Show Selected Rows and Row Source gate the highlight", async () => {
      await session.step(85, "When user sets \"Marker Color Column\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(86, "And user clicks on the \"category Caucasian\" area of box plot viewer", () => clickArea(page, "category Caucasian", el("box plot viewer")));
      await session.step(87, "Then only rows where \"RACE\" is \"Caucasian\" should be selected", () => onlyOfSelected(page, "RACE", "Caucasian"));
      await session.step(88, "And box plot viewer should show more selection highlight than before", () => moreHighlight(page, el("box plot viewer")));
      await session.step(89, "When user sets \"Show Selected Rows\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Selected Rows", el("box plot viewer"), "false"));
      await session.step(90, "Then box plot viewer should show less selection highlight than before", () => lessHighlight(page, el("box plot viewer")));
      await session.step(91, "And \"Marker Color Column\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(92, "When user sets \"Show Selected Rows\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Selected Rows", el("box plot viewer"), "true"));
      await session.step(93, "Then box plot viewer should show more selection highlight than before", () => moreHighlight(page, el("box plot viewer")));
      await session.step(94, "When user sets \"Row Source\" property of box plot viewer to \"Selected\"", () => setProperty(page, "Row Source", el("box plot viewer"), "Selected"));
      await session.step(95, "Then the \"Caucasian values\" area of box plot viewer should be painted", () => areaPainted(page, "Caucasian values", el("box plot viewer")));
      await session.step(96, "And box plot viewer should show no selection highlight", () => noHighlight(page, el("box plot viewer")));
      await session.step(97, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Row Source","All"],["Marker Color Column",""]]), [["Row Source","All"],["Marker Color Column",""]]);
      await session.step(100, "And user clears the row selection", () => clearSelection(page));
    });
    await run.scenario("The hover tooltip and Show Mouse Over Point", async () => {
      await session.step(103, "When user takes a snapshot of box plot viewer", () => takeSnapshot(page, el("box plot viewer")));
      await session.step(104, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(105, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(106, "And box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(107, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(108, "And user sets \"Show Mouse Over Point\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Mouse Over Point", el("box plot viewer"), "false"));
      await session.step(109, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(110, "Then box plot viewer should not have repainted", () => notRepainted(page, el("box plot viewer")));
      await session.step(111, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(112, "And user sets \"Show Mouse Over Point\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Mouse Over Point", el("box plot viewer"), "true"));
    });
    await run.scenario("A hover on a bar chart highlights its rows in the box plot only while the row group is shown (github-3065)", async () => {
      await session.step(115, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
      await session.step(117, "And user sets \"Marker Color Column\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(118, "Then the \"Caucasian values\" area of box plot viewer should contain the color \"#2CA02C\"", () => areaColor(page, "Caucasian values", el("box plot viewer"), "#2CA02C"));
      await session.step(119, "When user takes a snapshot of box plot viewer", () => takeSnapshot(page, el("box plot viewer")));
      await session.step(120, "And user hovers over the \"bar Caucasian\" area of bar chart viewer", () => hoverArea(page, "bar Caucasian", el("bar chart viewer")));
      await session.step(121, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(122, "And the \"Caucasian values\" area of box plot viewer should contain the color \"#2CA02C\"", () => areaColor(page, "Caucasian values", el("box plot viewer"), "#2CA02C"));
      await session.step(123, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
      await session.step(124, "And user sets \"Show Mouse Over Row Group\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Mouse Over Row Group", el("box plot viewer"), "false"));
      await session.step(125, "And user takes a snapshot of box plot viewer", () => takeSnapshot(page, el("box plot viewer")));
      await session.step(126, "And user hovers over the \"bar Caucasian\" area of bar chart viewer", () => hoverArea(page, "bar Caucasian", el("bar chart viewer")));
      await session.step(127, "Then box plot viewer should not have repainted", () => notRepainted(page, el("box plot viewer")));
      await session.step(128, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
      await session.step(129, "And user sets \"Show Mouse Over Row Group\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Mouse Over Row Group", el("box plot viewer"), "true"));
      await session.step(130, "And user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(131, "Then bar chart viewer should be absent", () => shouldBe(page, el("bar chart viewer"), "absent"));
      await session.step(132, "When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("A categorical coloring survives deleting the selected rows", async () => {
      await session.step(135, "When user sets \"Marker Color Column\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(136, "And user clicks on the \"category Other\" area of box plot viewer", () => clickArea(page, "category Other", el("box plot viewer")));
      await session.step(137, "Then only rows where \"RACE\" is \"Other\" should be selected", () => onlyOfSelected(page, "RACE", "Other"));
      await session.step(138, "When user deletes the selected rows", () => deleteSelected(page));
      await session.step(139, "Then the table should have no rows where \"RACE\" is \"Other\"", () => noRowsWhere(page, "RACE", "Other"));
      await session.step(140, "And box plot viewer should be painted", () => painted(page, el("box plot viewer")));
      await session.step(141, "And box plot viewer should not have a \"category Other\" area", () => hasNoArea(page, el("box plot viewer"), "category Other"));
      await session.step(142, "And no errors should have been logged", () => noErrors(page));
      await session.step(143, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(144, "And \"Marker Color Column\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "RACE"));
    });
    run.finish();
  });
});
