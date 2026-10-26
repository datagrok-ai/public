/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
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
import {clickOn, hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnIncomplete, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, addCategoricalFilter, addRangeFilter, clearSelection, filterPasses, filterPassesAll, filterPassesFewer, noneSelected, openFilterPanel, removeColumn, selectFirstRows, selectNoRows, selectWhereIs, selectedPassingCount, selectedRowCount, someOfSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickFirstArea, dragAcrossArea, dragSelectionOverArea, narrowerRange, noErrors, noHighlight, pickFromContextMenu, propertyShouldBe, readingAsRemembered, readingAtLeast, readingDoesNotRead, readingFinite, readingHigherThanRemembered, readingIs, readingNotAsRemembered, rememberRange, rememberReading, rememberedRange, reportsNoError, resizeTo, restoreSize, setProperties, setProperty, showsFewerRows, showsRows, someHighlight, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {currentRowIsHovered, dragRangeHandle} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot regressions — filtering, zoom, selection and labels", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-regressions.feature", import.meta.url);
  test("Scatter plot regressions — filtering, zoom, selection and labels", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(22, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(23, "And no rows should be selected", () => noneSelected(page));
    await run.scenario("Labels for selected rows are not drawn for the rows the filter took out (public#3189)", async () => {
      await session.step(27, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Label Columns","AGE"],["Show Labels For","Selected"]]), [["Label Columns","AGE"],["Show Labels For","Selected"]]);
      await session.step(30, "And user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(31, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(32, "When user selects rows where \"SEX\" is \"M\"", () => selectWhereIs(page, "SEX", "M"));
      await session.step(33, "Then the \"labels shown\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "labels shown", el("scatter plot viewer"), 1));
      await session.step(34, "When user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
      await session.step(35, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(36, "And 0 selected rows should pass the filter", () => selectedPassingCount(page, 0));
      await session.step(37, "And the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(38, "When user clears the row selection", () => clearSelection(page));
      await session.step(39, "And user adds a categorical filter on \"SEX\" keeping \"F, M\"", () => addCategoricalFilter(page, "SEX", "F, M"));
      await session.step(40, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Label Columns",""],["Show Labels For","All"]]), [["Label Columns",""],["Show Labels For","All"]]);
      await session.step(43, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(44, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A viewer Filter keeps its rows when a calculated column is added (public#1156)", async () => {
      await session.step(49, "When user sets \"Filter\" property of scatter plot viewer to '${RACE} in [\"Asian\", \"Black\"]'", () => setProperty(page, "Filter", el("scatter plot viewer"), "${RACE} in [\"Asian\", \"Black\"]"));
      await session.step(50, "Then scatter plot viewer should show 38 rows", () => showsRows(page, el("scatter plot viewer"), 38));
      await session.step(51, "When user adds a calculated column \"AGE2\" with formula \"${AGE} * 2\"", () => addCalculated(page, "AGE2", "${AGE} * 2"));
      await session.step(52, "Then the table should have a column \"AGE2\"", () => hasColumn(page, "AGE2"));
      await session.step(53, "And scatter plot viewer should show 38 rows", () => showsRows(page, el("scatter plot viewer"), 38));
      await session.step(54, "And scatter plot viewer should report no error", () => reportsNoError(page, el("scatter plot viewer")));
      await session.step(55, "When user removes \"AGE2\" column", () => removeColumn(page, "AGE2"));
      await session.step(56, "And user sets \"Filter\" property of scatter plot viewer to \"\"", () => setProperty(page, "Filter", el("scatter plot viewer"), ""));
      await session.step(57, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A zoom under Zoom and Filter = no action survives filtering in the filter panel (public#552)", async () => {
      await session.step(62, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Zoom and Filter","no action"],["axesFollowFilter","false"]]), [["Zoom and Filter","no action"],["axesFollowFilter","false"]]);
      await session.step(65, "And user drags the min handle of the \"x\" range slider of scatter plot viewer by 30 pixels", () => dragRangeHandle(page, "min", "x", el("scatter plot viewer"), 30));
      await session.step(66, "Then scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(67, "When user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(68, "And user opens the filter panel", () => openFilterPanel(page));
      await session.step(69, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(70, "When user adds a range filter on \"AGE\" from 20 to 40", () => addRangeFilter(page, "AGE", 20, 40));
      await session.step(71, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(72, "And scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(73, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(74, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(75, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(76, "And scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(77, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Zoom and Filter","filter by zoom"],["axesFollowFilter","true"]]), [["Zoom and Filter","filter by zoom"],["axesFollowFilter","true"]]);
      await session.step(80, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(81, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(82, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A plot shrunk almost to nothing filters no row, and the panel's reset gives every row back (public#1348)", async () => {
      await session.step(87, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(88, "And user resizes scatter plot viewer to 40 by 40", () => resizeTo(page, el("scatter plot viewer"), 40, 40));
      await session.step(89, "And user adds a range filter on \"WEIGHT\" from 100 to 170", () => addRangeFilter(page, "WEIGHT", 100, 170));
      await session.step(90, "Then 143 rows should pass the filter", () => filterPasses(page, 143));
      await session.step(91, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(92, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(93, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(94, "When user restores the size of scatter plot viewer", () => restoreSize(page, el("scatter plot viewer")));
      await session.step(95, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(96, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Mouse Drag = Select makes a plain drag select, and Pan makes it pan again (public#2089)", async () => {
      await session.step(101, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"no action\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "no action"));
      await session.step(102, "And user remembers the \"x axis min\" reading of scatter plot viewer", () => rememberReading(page, "x axis min", el("scatter plot viewer")));
      await session.step(103, "And user sets \"mouseDrag\" property of scatter plot viewer to \"Select\"", () => setProperty(page, "mouseDrag", el("scatter plot viewer"), "Select"));
      await session.step(104, "And user drags across the \"view\" area of scatter plot viewer", () => dragAcrossArea(page, "view", el("scatter plot viewer")));
      await session.step(105, "Then some rows should be selected", () => someSelected(page));
      await session.step(106, "And the \"x axis min\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis min", el("scatter plot viewer")));
      await session.step(107, "When user selects no rows", () => selectNoRows(page));
      await session.step(108, "And user sets \"mouseDrag\" property of scatter plot viewer to \"Pan\"", () => setProperty(page, "mouseDrag", el("scatter plot viewer"), "Pan"));
      await session.step(109, "And user drags across the \"view\" area of scatter plot viewer", () => dragAcrossArea(page, "view", el("scatter plot viewer")));
      await session.step(110, "Then the \"x axis min\" reading of scatter plot viewer should not be as remembered", () => readingNotAsRemembered(page, "x axis min", el("scatter plot viewer")));
      await session.step(111, "And no rows should be selected", () => noneSelected(page));
      await session.step(112, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(113, "And user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(114, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(115, "And the \"x axis min\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis min", el("scatter plot viewer")));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Selected Rows = false draws the selected markers without the highlight (public#3613)", async () => {
      await session.step(120, "When user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
      await session.step(121, "Then scatter plot viewer should show a selection highlight", () => someHighlight(page, el("scatter plot viewer")));
      await session.step(122, "When user sets \"showSelectedRows\" property of scatter plot viewer to \"false\"", () => setProperty(page, "showSelectedRows", el("scatter plot viewer"), "false"));
      await session.step(123, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(124, "And scatter plot viewer should show no selection highlight", () => noHighlight(page, el("scatter plot viewer")));
      await session.step(125, "When user sets \"showSelectedRows\" property of scatter plot viewer to \"true\"", () => setProperty(page, "showSelectedRows", el("scatter plot viewer"), "true"));
      await session.step(126, "Then scatter plot viewer should show a selection highlight", () => someHighlight(page, el("scatter plot viewer")));
      await session.step(127, "When user clears the row selection", () => clearSelection(page));
      await session.step(128, "Then no rows should be selected", () => noneSelected(page));
      await session.step(129, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every label of the selected rows comes back after a zoom in and out (public#3176)", async () => {
      await session.step(133, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Zoom and Filter","no action"],["Label Columns","AGE"],["Show Labels For","Selected"]]), [["Zoom and Filter","no action"],["Label Columns","AGE"],["Show Labels For","Selected"]]);
      await session.step(137, "And user selects the first 5 rows", () => selectFirstRows(page, 5));
      await session.step(138, "Then the \"labels shown\" reading of scatter plot viewer should be 5", () => readingIs(page, "labels shown", el("scatter plot viewer"), 5));
      await session.step(139, "When user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(140, "And user scrolls the mouse wheel up 4 times over the \"view\" area of scatter plot viewer", () => wheelOverAreaTimes(page, "up", 4, "view", el("scatter plot viewer")));
      await session.step(141, "Then scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(142, "When user scrolls the mouse wheel down 4 times over the \"view\" area of scatter plot viewer", () => wheelOverAreaTimes(page, "down", 4, "view", el("scatter plot viewer")));
      await session.step(143, "Then the \"labels shown\" reading of scatter plot viewer should be 5", () => readingIs(page, "labels shown", el("scatter plot viewer"), 5));
      await session.step(144, "When user clears the row selection", () => clearSelection(page));
      await session.step(145, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(146, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Label Columns",""],["Show Labels For","All"],["Zoom and Filter","filter by zoom"]]), [["Label Columns",""],["Show Labels For","All"],["Zoom and Filter","filter by zoom"]]);
      await session.step(150, "Then the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(151, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(152, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("After a zoom a click on a marker makes its row current (public#801)", async () => {
      await session.step(156, "When user scrolls the mouse wheel up 2 times over the \"view\" area of scatter plot viewer", () => wheelOverAreaTimes(page, "up", 2, "view", el("scatter plot viewer")));
      await session.step(157, "Then scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(158, "When user clicks on the first \"marker of row\" area of scatter plot viewer", () => clickFirstArea(page, "marker of row", el("scatter plot viewer")));
      await session.step(159, "Then the \"hovered row\" reading of scatter plot viewer should be the current row", () => currentRowIsHovered(page, "hovered row", el("scatter plot viewer")));
      await session.step(160, "And the \"hovered row\" reading of scatter plot viewer should not be \"1\"", () => readingDoesNotRead(page, "hovered row", el("scatter plot viewer"), "1"));
      await session.step(161, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(162, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(163, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(164, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A selection box takes the markers whose Color value is empty (public#1174)", async () => {
      await session.step(168, "When user adds a calculated column \"RACE_F\" with formula \"if(${SEX} == \\\"F\\\", ${RACE}, null)\"", () => addCalculated(page, "RACE_F", "if(${SEX} == \"F\", ${RACE}, null)"));
      await session.step(169, "And user sets \"Color\" property of scatter plot viewer to \"RACE_F\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE_F"));
      await session.step(170, "And user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(171, "Then some rows where \"SEX\" is \"M\" should be selected", () => someOfSelected(page, "SEX", "M"));
      await session.step(172, "And some rows where \"SEX\" is \"F\" should be selected", () => someOfSelected(page, "SEX", "F"));
      await session.step(173, "When user clears the row selection", () => clearSelection(page));
      await session.step(174, "And user sets \"Color\" property of scatter plot viewer to \"\"", () => setProperty(page, "Color", el("scatter plot viewer"), ""));
      await session.step(175, "And user removes \"RACE_F\" column", () => removeColumn(page, "RACE_F"));
      await session.step(176, "Then no rows should be selected", () => noneSelected(page));
      await session.step(177, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The range slider moves a logarithmic axis (public#376)", async () => {
      await session.step(181, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Zoom and Filter","no action"],["X Axis Type","logarithmic"]]), [["Zoom and Filter","no action"],["X Axis Type","logarithmic"]]);
      await session.step(184, "And user remembers the \"x axis min\" reading of scatter plot viewer", () => rememberReading(page, "x axis min", el("scatter plot viewer")));
      await session.step(185, "And user drags the min handle of the \"x\" range slider of scatter plot viewer by 30 pixels", () => dragRangeHandle(page, "min", "x", el("scatter plot viewer"), 30));
      await session.step(186, "Then the \"x axis min\" reading of scatter plot viewer should be higher than remembered", () => readingHigherThanRemembered(page, "x axis min", el("scatter plot viewer")));
      await session.step(187, "And the \"x axis min\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "x axis min", el("scatter plot viewer")));
      await session.step(188, "And the \"x axis max\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "x axis max", el("scatter plot viewer")));
      await session.step(189, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(190, "Then the \"x axis min\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis min", el("scatter plot viewer")));
      await session.step(191, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X Axis Type","linear"],["Zoom and Filter","filter by zoom"]]), [["X Axis Type","linear"],["Zoom and Filter","filter by zoom"]]);
      await session.step(194, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(195, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(196, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic axis on an all-empty column under filter by zoom filters no row (public#1744)", async () => {
      await session.step(200, "When user adds a calculated column \"NO_VALUES\" with formula \"if(${AGE} > 1000, ${WEIGHT}, null)\"", () => addCalculated(page, "NO_VALUES", "if(${AGE} > 1000, ${WEIGHT}, null)"));
      await session.step(201, "Then \"NO_VALUES\" column should have missing values", () => columnIncomplete(page, "NO_VALUES"));
      await session.step(202, "And \"Zoom and Filter\" property of scatter plot viewer should be \"filter by zoom\"", () => propertyShouldBe(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(203, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y Axis Type","logarithmic"],["Y","NO_VALUES"]]), [["Y Axis Type","logarithmic"],["Y","NO_VALUES"]]);
      await session.step(206, "Then scatter plot viewer should show 0 rows", () => showsRows(page, el("scatter plot viewer"), 0));
      await session.step(207, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(208, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y","HEIGHT"],["Y Axis Type","linear"]]), [["Y","HEIGHT"],["Y Axis Type","linear"]]);
      await session.step(211, "And user removes \"NO_VALUES\" column", () => removeColumn(page, "NO_VALUES"));
      await session.step(212, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(213, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(214, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(215, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
