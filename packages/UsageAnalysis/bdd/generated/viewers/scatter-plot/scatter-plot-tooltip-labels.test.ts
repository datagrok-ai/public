/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-tooltip-labels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, public#3189, public#3176]
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
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, clearSelection, filterPasses, filterPassesAll, noneSelected, selectFirstRows, selectWhereIs, selectedPassFilter, selectedPassingCount, selectedRowCount, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, dragSelectionOverArea, hoverArea, narrowerRange, noErrors, painted, pickFromContextMenu, pointerAway, propertiesShouldBe, readingAtLeast, readingHigher, readingIs, rememberRange, setProperties, setProperty, showsRows, tooltipColumns, tooltipSomeColumns, tooltipValue, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot tooltip and marker labels", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-tooltip-labels.feature", import.meta.url);
  test("Scatter plot tooltip and marker labels", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot", "@realizes:public#3189", "@realizes:public#3176"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(24, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(25, "And properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["Show Tooltip","inherit from table"],["Label Columns",""]]), [["Show Tooltip","inherit from table"],["Label Columns",""]]);
    await run.scenario("Hovering a marker makes its row hovered and shows the table's tooltip", async () => {
      await session.step(30, "Then the \"hovered row\" reading of scatter plot viewer should be 0", () => readingIs(page, "hovered row", el("scatter plot viewer"), 0));
      await session.step(31, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(32, "Then the \"hovered row\" reading of scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("scatter plot viewer"), 11));
      await session.step(33, "And tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(34, "And the tooltip should show some columns", () => tooltipSomeColumns(page));
      await session.step(35, "And the tooltip should show \"HEIGHT\" as \"168.974\"", () => tooltipValue(page, "HEIGHT", "168.974"));
      await session.step(36, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(37, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A custom column list replaces what the tooltip shows", async () => {
      await session.step(41, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Tooltip","show custom tooltip"],["Row Tooltip","AGE\nSEX"],["Data Values","Do not add"],["showLabels","Always"]]), [["Show Tooltip","show custom tooltip"],["Row Tooltip","AGE\nSEX"],["Data Values","Do not add"],["showLabels","Always"]]);
      await session.step(46, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(47, "Then the \"hovered row\" reading of scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("scatter plot viewer"), 11));
      await session.step(48, "And the tooltip should show columns \"AGE, SEX\"", () => tooltipColumns(page, "AGE, SEX"));
      await session.step(49, "And the tooltip should show \"AGE\" as \"46\"", () => tooltipValue(page, "AGE", "46"));
      await session.step(50, "And the tooltip should show \"SEX\" as \"F\"", () => tooltipValue(page, "SEX", "F"));
      await session.step(51, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(52, "And user sets \"Show Tooltip\" property of scatter plot viewer to \"do not show\"", () => setProperty(page, "Show Tooltip", el("scatter plot viewer"), "do not show"));
      await session.step(53, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(54, "Then the \"hovered row\" reading of scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("scatter plot viewer"), 11));
      await session.step(55, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(56, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(57, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Tooltip","inherit from table"],["Row Tooltip",""],["Data Values","Merge"],["showLabels","Auto"]]), [["Show Tooltip","inherit from table"],["Row Tooltip",""],["Data Values","Merge"],["showLabels","Auto"]]);
      await session.step(62, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(63, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(64, "And the tooltip should show some columns", () => tooltipSomeColumns(page));
      await session.step(65, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(66, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Labels are drawn for every row, then only for the selected ones", async () => {
      await session.step(69, "Then no rows should be selected", () => noneSelected(page));
      await session.step(70, "And the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(71, "When user sets \"Label Columns\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Label Columns", el("scatter plot viewer"), "AGE"));
      await session.step(72, "Then the \"labels shown\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "labels shown", el("scatter plot viewer"), 1));
      await session.step(73, "When user sets \"Show Labels For\" property of scatter plot viewer to \"Selected\"", () => setProperty(page, "Show Labels For", el("scatter plot viewer"), "Selected"));
      await session.step(74, "Then the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(75, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(76, "Then some rows should be selected", () => someSelected(page));
      await session.step(77, "And the \"labels shown\" reading of scatter plot viewer should be higher than before", () => readingHigher(page, "labels shown", el("scatter plot viewer")));
      await session.step(78, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(79, "Then no rows should be selected", () => noneSelected(page));
      await session.step(80, "And the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(81, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Label Columns",""],["Show Labels For","All"]]), [["Label Columns",""],["Show Labels For","All"]]);
      await session.step(84, "Then the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Labels on a datetime axis survive a filter that keeps one category and a selection", async () => {
      await session.step(88, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","STARTED"],["Y","RACE"],["Label Columns","USUBJID"],["Show Labels For","Selected"]]), [["X","STARTED"],["Y","RACE"],["Label Columns","USUBJID"],["Show Labels For","Selected"]]);
      await session.step(93, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
      await session.step(94, "When user adds a categorical filter on \"RACE\" keeping \"Caucasian\"", () => addCategoricalFilter(page, "RACE", "Caucasian"));
      await session.step(95, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(96, "And scatter plot viewer should show 896 rows", () => showsRows(page, el("scatter plot viewer"), 896));
      await session.step(97, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(98, "Then some rows should be selected", () => someSelected(page));
      await session.step(99, "And every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(100, "And the \"labels shown\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "labels shown", el("scatter plot viewer"), 1));
      await session.step(101, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
      await session.step(103, "When user adds a categorical filter on \"RACE\" keeping \"Asian, Black, Caucasian, Other\"", () => addCategoricalFilter(page, "RACE", "Asian, Black, Caucasian, Other"));
      await session.step(104, "And user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(105, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(106, "And no rows should be selected", () => noneSelected(page));
      await session.step(107, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","WEIGHT"],["Y","HEIGHT"],["Label Columns",""],["Show Labels For","All"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Label Columns",""],["Show Labels For","All"]]);
      await session.step(112, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Labels for selected rows are not drawn for the rows the filter took out (public#3189)", async () => {
      await session.step(117, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Label Columns","AGE"],["Show Labels For","Selected"]]), [["Label Columns","AGE"],["Show Labels For","Selected"]]);
      await session.step(120, "And user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(121, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(122, "When user selects rows where \"SEX\" is \"M\"", () => selectWhereIs(page, "SEX", "M"));
      await session.step(123, "Then the \"labels shown\" reading of scatter plot viewer should be 373", () => readingIs(page, "labels shown", el("scatter plot viewer"), 373));
      await session.step(124, "When user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
      await session.step(125, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(126, "And 0 selected rows should pass the filter", () => selectedPassingCount(page, 0));
      await session.step(127, "And the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(128, "When user clears the row selection", () => clearSelection(page));
      await session.step(129, "And user adds a categorical filter on \"SEX\" keeping \"F, M\"", () => addCategoricalFilter(page, "SEX", "F, M"));
      await session.step(130, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Label Columns",""],["Show Labels For","All"]]), [["Label Columns",""],["Show Labels For","All"]]);
      await session.step(133, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(134, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every label of the selected rows comes back after a zoom in and out (public#3176)", async () => {
      await session.step(139, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Zoom and Filter","no action"],["Label Columns","AGE"],["Show Labels For","Selected"]]), [["Zoom and Filter","no action"],["Label Columns","AGE"],["Show Labels For","Selected"]]);
      await session.step(143, "And user selects the first 5 rows", () => selectFirstRows(page, 5));
      await session.step(144, "Then the \"labels shown\" reading of scatter plot viewer should be 5", () => readingIs(page, "labels shown", el("scatter plot viewer"), 5));
      await session.step(145, "When user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(146, "And user scrolls the mouse wheel up 4 times over the \"view\" area of scatter plot viewer", () => wheelOverAreaTimes(page, "up", 4, "view", el("scatter plot viewer")));
      await session.step(147, "Then scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(148, "When user scrolls the mouse wheel down 4 times over the \"view\" area of scatter plot viewer", () => wheelOverAreaTimes(page, "down", 4, "view", el("scatter plot viewer")));
      await session.step(149, "Then the \"labels shown\" reading of scatter plot viewer should be 5", () => readingIs(page, "labels shown", el("scatter plot viewer"), 5));
      await session.step(150, "When user clears the row selection", () => clearSelection(page));
      await session.step(151, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(152, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Label Columns",""],["Show Labels For","All"],["Zoom and Filter","filter by zoom"]]), [["Label Columns",""],["Show Labels For","All"],["Zoom and Filter","filter by zoom"]]);
      await session.step(156, "Then the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(157, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(158, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
