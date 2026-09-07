/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot-selection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clickEmptySpace, doubleClickEmptySpace} from '../../bindings/box-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {allOfSelected, clearSelection, deleteSelected, filterOut, hasCurrentRow, noRowsWhere, noneOfSelected, noneSelected, onlyOfSelected, resetFilter, someOfSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, clickAreaHolding, dragSelectionOverArea, eventFired, eventNotFired, hasNoArea, hoverArea, lessHighlight, listenFor, moreHighlight, noBalloons, noErrors, noHighlight, notRepainted, painted, pointerAway, propertyShouldBe, setProperties, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot selection and highlight", () => {
  const session = feature(test, "features/viewers/box-plot-selection.feature", import.meta.url);
  test("Box plot selection and highlight", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"],["Marker Size","10"]]));
    await run.scenario("A marker click sets the current row", async () => {
      await session.step(20, "Given user listens for \"d4-boxplot-point-click\" event on box plot viewer", () => listenFor(page, "d4-boxplot-point-click", el("box plot viewer")));
      await session.step(21, "When user clicks on the \"marker\" area of box plot viewer", () => clickArea(page, "marker", el("box plot viewer")));
      await session.step(22, "Then \"d4-boxplot-point-click\" event should have fired on box plot viewer", () => eventFired(page, "d4-boxplot-point-click", el("box plot viewer")));
      await session.step(23, "And the table should have a current row", () => hasCurrentRow(page));
    });
    await run.scenario("A Shift-drag selects a band and highlights it", async () => {
      await session.step(26, "When user clears the row selection", () => clearSelection(page));
      await session.step(27, "And user drags a selection box over the \"Caucasian values\" area of box plot viewer", () => dragSelectionOverArea(page, "Caucasian values", el("box plot viewer")));
      await session.step(28, "Then some rows where \"RACE\" is \"Caucasian\" should be selected", () => someOfSelected(page, "RACE", "Caucasian"));
      await session.step(29, "And box plot viewer should show more selection highlight than before", () => moreHighlight(page, el("box plot viewer")));
    });
    await run.scenario("The highlight survives a categorical coloring", async () => {
      await session.step(32, "When user clears the row selection", () => clearSelection(page));
      await session.step(33, "And user sets \"Marker Color Column\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "SEX"));
      await session.step(34, "And user drags a selection box over the \"Caucasian values\" area of box plot viewer", () => dragSelectionOverArea(page, "Caucasian values", el("box plot viewer")));
      await session.step(35, "Then some rows where \"RACE\" is \"Caucasian\" should be selected", () => someOfSelected(page, "RACE", "Caucasian"));
      await session.step(36, "And box plot viewer should show more selection highlight than before", () => moreHighlight(page, el("box plot viewer")));
      await session.step(37, "When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("Category labels select categories", async () => {
      await session.step(40, "When user clears the row selection", () => clearSelection(page));
      await session.step(41, "And user drags a selection box over the \"Caucasian values\" area of box plot viewer", () => dragSelectionOverArea(page, "Caucasian values", el("box plot viewer")));
      await session.step(42, "Then some rows where \"RACE\" is \"Caucasian\" should be selected", () => someOfSelected(page, "RACE", "Caucasian"));
      await session.step(43, "When user clicks on the \"category Black\" area of box plot viewer holding Control", () => clickAreaHolding(page, "category Black", el("box plot viewer"), "Control"));
      await session.step(44, "Then all rows where \"RACE\" is \"Black\" should be selected", () => allOfSelected(page, "RACE", "Black"));
      await session.step(45, "And some rows where \"RACE\" is \"Caucasian\" should be selected", () => someOfSelected(page, "RACE", "Caucasian"));
      await session.step(46, "When user clicks on the \"category Asian\" area of box plot viewer", () => clickArea(page, "category Asian", el("box plot viewer")));
      await session.step(47, "Then only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(48, "Given user listens for \"d4-boxplot-reset-view\" event on box plot viewer", () => listenFor(page, "d4-boxplot-reset-view", el("box plot viewer")));
      await session.step(49, "When user clicks on empty plot space of box plot viewer", () => clickEmptySpace(page, el("box plot viewer")));
      await session.step(50, "Then no rows should be selected", () => noneSelected(page));
      await session.step(51, "And \"d4-boxplot-reset-view\" event should not have fired on box plot viewer", () => eventNotFired(page, "d4-boxplot-reset-view", el("box plot viewer")));
      await session.step(52, "When user double-clicks on empty plot space of box plot viewer", () => doubleClickEmptySpace(page, el("box plot viewer")));
      await session.step(53, "Then \"d4-boxplot-reset-view\" event should have fired on box plot viewer", () => eventFired(page, "d4-boxplot-reset-view", el("box plot viewer")));
    });
    await run.scenario("No selection leaks into a filtered-out category", async () => {
      await session.step(56, "When user clears the row selection", () => clearSelection(page));
      await session.step(57, "And user filters out rows where \"RACE\" is \"Asian\"", () => filterOut(page, "RACE", "Asian"));
      await session.step(58, "And user drags a selection box over the \"view\" area of box plot viewer", () => dragSelectionOverArea(page, "view", el("box plot viewer")));
      await session.step(59, "Then some rows should be selected", () => someSelected(page));
      await session.step(60, "And no rows where \"RACE\" is \"Asian\" should be selected", () => noneOfSelected(page, "RACE", "Asian"));
      await session.step(61, "When user resets the filter", () => resetFilter(page));
      await session.step(62, "And user clears the row selection", () => clearSelection(page));
    });
    await run.scenario("Show Selected Rows and Row Source gate the highlight", async () => {
      await session.step(65, "When user sets \"Marker Color Column\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(66, "And user clicks on the \"category Caucasian\" area of box plot viewer", () => clickArea(page, "category Caucasian", el("box plot viewer")));
      await session.step(67, "Then only rows where \"RACE\" is \"Caucasian\" should be selected", () => onlyOfSelected(page, "RACE", "Caucasian"));
      await session.step(68, "And box plot viewer should show more selection highlight than before", () => moreHighlight(page, el("box plot viewer")));
      await session.step(69, "When user sets \"Show Selected Rows\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Selected Rows", el("box plot viewer"), "false"));
      await session.step(70, "Then box plot viewer should show less selection highlight than before", () => lessHighlight(page, el("box plot viewer")));
      await session.step(71, "When user sets \"Show Selected Rows\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Selected Rows", el("box plot viewer"), "true"));
      await session.step(72, "And user sets \"Row Source\" property of box plot viewer to \"Selected\"", () => setProperty(page, "Row Source", el("box plot viewer"), "Selected"));
      await session.step(73, "Then box plot viewer should show no selection highlight", () => noHighlight(page, el("box plot viewer")));
      await session.step(74, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Row Source","All"],["Marker Color Column",""]]));
      await session.step(77, "And user clears the row selection", () => clearSelection(page));
    });
    await run.scenario("The hover tooltip and Show Mouse Over Point", async () => {
      await session.step(80, "When user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(81, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(82, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(83, "And user sets \"Show Mouse Over Point\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Mouse Over Point", el("box plot viewer"), "false"));
      await session.step(84, "And user hovers over the \"marker\" area of box plot viewer", () => hoverArea(page, "marker", el("box plot viewer")));
      await session.step(85, "Then box plot viewer should not have repainted", () => notRepainted(page, el("box plot viewer")));
      await session.step(86, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(87, "And user sets \"Show Mouse Over Point\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Mouse Over Point", el("box plot viewer"), "true"));
    });
    await run.scenario("A hover on a bar chart leaves the box plot alone", async () => {
      await session.step(90, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]));
      await session.step(92, "And user sets \"Marker Color Column\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(93, "And user takes a snapshot of box plot viewer", () => takeSnapshot(page, el("box plot viewer")));
      await session.step(94, "And user hovers over the \"bar Caucasian\" area of bar chart viewer", () => hoverArea(page, "bar Caucasian", el("bar chart viewer")));
      await session.step(95, "Then box plot viewer should not have repainted", () => notRepainted(page, el("box plot viewer")));
      await session.step(96, "When user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(97, "Then bar chart viewer should be absent", () => shouldBe(page, el("bar chart viewer"), "absent"));
      await session.step(98, "When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("A categorical coloring survives deleting the selected rows", async () => {
      await session.step(101, "When user sets \"Marker Color Column\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(102, "And user clicks on the \"category Other\" area of box plot viewer", () => clickArea(page, "category Other", el("box plot viewer")));
      await session.step(103, "Then only rows where \"RACE\" is \"Other\" should be selected", () => onlyOfSelected(page, "RACE", "Other"));
      await session.step(104, "When user deletes the selected rows", () => deleteSelected(page));
      await session.step(105, "Then the table should have no rows where \"RACE\" is \"Other\"", () => noRowsWhere(page, "RACE", "Other"));
      await session.step(106, "And box plot viewer should be painted", () => painted(page, el("box plot viewer")));
      await session.step(107, "And box plot viewer should not have a \"category Other\" area", () => hasNoArea(page, el("box plot viewer"), "category Other"));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
      await session.step(109, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(110, "And \"Marker Color Column\" property of box plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "RACE"));
    });
    run.finish();
  });
});
