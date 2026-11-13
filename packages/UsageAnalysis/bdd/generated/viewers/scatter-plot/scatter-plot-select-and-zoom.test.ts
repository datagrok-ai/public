/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-select-and-zoom.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, public#3613, public#801, public#1174, public#376, public#2089]
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
import {pressKey} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, clearSelection, filterPasses, filterPassesAll, filterTo, noneSelected, onlyOfSelected, removeColumn, resetFilter, selectFirstRows, selectWhereIs, selectedPassFilter, selectedRowCount, someOfSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, clickFirstArea, doubleClickArea, dragAcrossArea, dragDeselectionOverArea, dragSelectionBetweenAreas, dragSelectionOverArea, dragZoomOverArea, eventFired, hoverArea, listenFor, narrowerRange, noErrors, noHighlight, pickFromContextMenu, propertyShouldBe, readingAsRemembered, readingDiffers, readingDoesNotRead, readingFinite, readingHigher, readingHigherThanRemembered, readingIs, readingLower, readingSame, rememberRange, rememberReading, rememberedRange, repainted, setProperties, setProperty, showsFewerRows, showsRows, someHighlight, takeSnapshot, wheelOverArea, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {currentRowIsHovered, dragLasso, dragRangeHandle} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot selection and viewport navigation", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-select-and-zoom.feature", import.meta.url);
  test("Scatter plot selection and viewport navigation", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot", "@realizes:public#3613", "@realizes:public#801", "@realizes:public#1174", "@realizes:public#376", "@realizes:public#2089"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 12, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(23, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(24, "And no rows should be selected", () => noneSelected(page));
    await run.scenario("A drag selects the markers inside it, a second drag adds, Control+Shift takes them out", async () => {
      await session.step(27, "Then \"Reset Selection On Background Click\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Reset Selection On Background Click", el("scatter plot viewer"), "true"));
      await session.step(28, "When user drags a selection box from the \"marker of row 11\" area to the \"marker of row 3\" area of scatter plot viewer", () => dragSelectionBetweenAreas(page, "marker of row 11", "marker of row 3", el("scatter plot viewer")));
      await session.step(29, "Then some rows should be selected", () => someSelected(page));
      await session.step(30, "And the \"rows selected\" reading of scatter plot viewer should be higher than before", () => readingHigher(page, "rows selected", el("scatter plot viewer")));
      await session.step(31, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(32, "Then the \"rows selected\" reading of scatter plot viewer should be higher than before", () => readingHigher(page, "rows selected", el("scatter plot viewer")));
      await session.step(33, "And every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(34, "When user drags a deselection box over the \"view\" area of scatter plot viewer", () => dragDeselectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(35, "Then the \"rows selected\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "rows selected", el("scatter plot viewer")));
      await session.step(36, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(37, "Then some rows should be selected", () => someSelected(page));
      await session.step(38, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(39, "Then no rows should be selected", () => noneSelected(page));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The selection survives a jitter change", async () => {
      await session.step(43, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Jitter Size","20"],["Jitter Size Y","15"]]), [["Jitter Size","20"],["Jitter Size Y","15"]]);
      await session.step(46, "And user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(47, "Then some rows should be selected", () => someSelected(page));
      await session.step(48, "When user remembers the \"rows selected\" reading of scatter plot viewer", () => rememberReading(page, "rows selected", el("scatter plot viewer")));
      await session.step(49, "And user sets \"Jitter Size\" property of scatter plot viewer to \"30\"", () => setProperty(page, "Jitter Size", el("scatter plot viewer"), "30"));
      await session.step(50, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(51, "And the \"rows selected\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "rows selected", el("scatter plot viewer")));
      await session.step(52, "When user drags a deselection box over the \"view\" area of scatter plot viewer", () => dragDeselectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(53, "Then the \"rows selected\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "rows selected", el("scatter plot viewer")));
      await session.step(54, "When user clears the row selection", () => clearSelection(page));
      await session.step(55, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Jitter Size","0"],["Jitter Size Y","0"]]), [["Jitter Size","0"],["Jitter Size Y","0"]]);
      await session.step(58, "Then no rows should be selected", () => noneSelected(page));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An Alt-drag zooms, a plain drag pans, the wheel zooms and Reset View comes home", async () => {
      await session.step(62, "When user selects the first 20 rows", () => selectFirstRows(page, 20));
      await session.step(63, "Then 20 rows should be selected", () => selectedRowCount(page, 20));
      await session.step(64, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"no action\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "no action"));
      await session.step(65, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(66, "And user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(67, "And user remembers the \"x axis span\" reading of scatter plot viewer", () => rememberReading(page, "x axis span", el("scatter plot viewer")));
      await session.step(68, "And user drags a zoom box over the \"view\" area of scatter plot viewer", () => dragZoomOverArea(page, "view", el("scatter plot viewer")));
      await session.step(69, "Then the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(70, "And scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(71, "And the \"x axis span\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "x axis span", el("scatter plot viewer")));
      await session.step(72, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(73, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(74, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(75, "And the \"x axis span\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis span", el("scatter plot viewer")));
      await session.step(76, "When user drags across the \"view\" area of scatter plot viewer", () => dragAcrossArea(page, "view", el("scatter plot viewer")));
      await session.step(77, "Then the \"x axis min\" reading of scatter plot viewer should differ from before", () => readingDiffers(page, "x axis min", el("scatter plot viewer")));
      await session.step(78, "And the \"x axis span\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "x axis span", el("scatter plot viewer")));
      await session.step(79, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(80, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(81, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(82, "Then the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(83, "And scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(84, "And the \"x axis span\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "x axis span", el("scatter plot viewer")));
      await session.step(85, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(86, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(87, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(88, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(89, "And user clears the row selection", () => clearSelection(page));
      await session.step(90, "Then no rows should be selected", () => noneSelected(page));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The X range slider narrows the axis and a double-click on empty space resets the view", async () => {
      await session.step(94, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"no action\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "no action"));
      await session.step(95, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(96, "And user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(97, "And user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(98, "And user drags the min handle of the \"x\" range slider of scatter plot viewer by 30 pixels", () => dragRangeHandle(page, "min", "x", el("scatter plot viewer"), 30));
      await session.step(99, "Then the \"x axis min\" reading of scatter plot viewer should be higher than before", () => readingHigher(page, "x axis min", el("scatter plot viewer")));
      await session.step(100, "And the \"x axis span\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "x axis span", el("scatter plot viewer")));
      await session.step(101, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(102, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(103, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(104, "Then the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(105, "And scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(106, "When user listens for \"d4-scatterplot-reset-view\" event on scatter plot viewer", () => listenFor(page, "d4-scatterplot-reset-view", el("scatter plot viewer")));
      await session.step(107, "And user hovers over the \"empty space\" area of scatter plot viewer", () => hoverArea(page, "empty space", el("scatter plot viewer")));
      await session.step(108, "Then the \"hovered row\" reading of scatter plot viewer should be 0", () => readingIs(page, "hovered row", el("scatter plot viewer"), 0));
      await session.step(109, "When user double-clicks on the \"empty space\" area of scatter plot viewer", () => doubleClickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(110, "Then \"d4-scatterplot-reset-view\" event should have fired on scatter plot viewer", () => eventFired(page, "d4-scatterplot-reset-view", el("scatter plot viewer")));
      await session.step(111, "And scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(112, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(113, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+A selects the rows the filter passes, Ctrl+Shift+A and Escape clear them", async () => {
      await session.step(116, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(117, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(118, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(119, "And user presses Control+a", () => pressKey(page, "Control+a"));
      await session.step(120, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(121, "And only rows where \"SEX\" is \"F\" should be selected", () => onlyOfSelected(page, "SEX", "F"));
      await session.step(122, "When user presses Control+Shift+a", () => pressKey(page, "Control+Shift+a"));
      await session.step(123, "Then no rows should be selected", () => noneSelected(page));
      await session.step(124, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(125, "Then some rows should be selected", () => someSelected(page));
      await session.step(126, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(127, "Then no rows should be selected", () => noneSelected(page));
      await session.step(128, "When user resets the filter", () => resetFilter(page));
      await session.step(129, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(130, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(131, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("H restores the viewport after a wheel zoom", async () => {
      await session.step(134, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(135, "And user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(136, "And user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(137, "Then the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(138, "And scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(139, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(140, "And user presses h", () => pressKey(page, "h"));
      await session.step(141, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(142, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(143, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("L turns on the Lasso Tool and a lasso selects the markers it encloses", async () => {
      await session.step(146, "Then \"Lasso Tool\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "Lasso Tool", el("scatter plot viewer"), "false"));
      await session.step(147, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(148, "And user presses l", () => pressKey(page, "l"));
      await session.step(149, "Then \"Lasso Tool\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Lasso Tool", el("scatter plot viewer"), "true"));
      await session.step(150, "When user drags a lasso over the \"view\" area of scatter plot viewer", () => dragLasso(page, "view", el("scatter plot viewer")));
      await session.step(151, "Then some rows should be selected", () => someSelected(page));
      await session.step(152, "And every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(153, "When user presses l", () => pressKey(page, "l"));
      await session.step(154, "Then \"Lasso Tool\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "Lasso Tool", el("scatter plot viewer"), "false"));
      await session.step(155, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(156, "Then no rows should be selected", () => noneSelected(page));
      await session.step(157, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Selected Rows = false draws the selected markers without the highlight (public#3613)", async () => {
      await session.step(161, "When user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
      await session.step(162, "Then scatter plot viewer should show a selection highlight", () => someHighlight(page, el("scatter plot viewer")));
      await session.step(163, "When user sets \"showSelectedRows\" property of scatter plot viewer to \"false\"", () => setProperty(page, "showSelectedRows", el("scatter plot viewer"), "false"));
      await session.step(164, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(165, "And scatter plot viewer should show no selection highlight", () => noHighlight(page, el("scatter plot viewer")));
      await session.step(166, "When user sets \"showSelectedRows\" property of scatter plot viewer to \"true\"", () => setProperty(page, "showSelectedRows", el("scatter plot viewer"), "true"));
      await session.step(167, "Then scatter plot viewer should show a selection highlight", () => someHighlight(page, el("scatter plot viewer")));
      await session.step(168, "When user clears the row selection", () => clearSelection(page));
      await session.step(169, "Then no rows should be selected", () => noneSelected(page));
      await session.step(170, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("After a zoom a click on a marker makes its row current (public#801)", async () => {
      await session.step(174, "When user scrolls the mouse wheel up 2 times over the \"view\" area of scatter plot viewer", () => wheelOverAreaTimes(page, "up", 2, "view", el("scatter plot viewer")));
      await session.step(175, "Then scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(176, "When user clicks on the first \"marker of row\" area of scatter plot viewer", () => clickFirstArea(page, "marker of row", el("scatter plot viewer")));
      await session.step(177, "Then the \"hovered row\" reading of scatter plot viewer should be the current row", () => currentRowIsHovered(page, "hovered row", el("scatter plot viewer")));
      await session.step(178, "And the \"hovered row\" reading of scatter plot viewer should not be \"1\"", () => readingDoesNotRead(page, "hovered row", el("scatter plot viewer"), "1"));
      await session.step(179, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(180, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(181, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(182, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A selection box takes the markers whose Color value is empty (public#1174)", async () => {
      await session.step(186, "When user adds a calculated column \"RACE_F\" with formula \"if(${SEX} == \\\"F\\\", ${RACE}, null)\"", () => addCalculated(page, "RACE_F", "if(${SEX} == \"F\", ${RACE}, null)"));
      await session.step(187, "And user sets \"Color\" property of scatter plot viewer to \"RACE_F\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE_F"));
      await session.step(188, "And user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(189, "Then some rows where \"SEX\" is \"M\" should be selected", () => someOfSelected(page, "SEX", "M"));
      await session.step(190, "And some rows where \"SEX\" is \"F\" should be selected", () => someOfSelected(page, "SEX", "F"));
      await session.step(191, "When user clears the row selection", () => clearSelection(page));
      await session.step(192, "And user sets \"Color\" property of scatter plot viewer to \"\"", () => setProperty(page, "Color", el("scatter plot viewer"), ""));
      await session.step(193, "And user removes \"RACE_F\" column", () => removeColumn(page, "RACE_F"));
      await session.step(194, "Then no rows should be selected", () => noneSelected(page));
      await session.step(195, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The range slider moves a logarithmic axis (public#376)", async () => {
      await session.step(199, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Zoom and Filter","no action"],["X Axis Type","logarithmic"]]), [["Zoom and Filter","no action"],["X Axis Type","logarithmic"]]);
      await session.step(202, "And user remembers the \"x axis min\" reading of scatter plot viewer", () => rememberReading(page, "x axis min", el("scatter plot viewer")));
      await session.step(203, "And user drags the min handle of the \"x\" range slider of scatter plot viewer by 30 pixels", () => dragRangeHandle(page, "min", "x", el("scatter plot viewer"), 30));
      await session.step(204, "Then the \"x axis min\" reading of scatter plot viewer should be higher than remembered", () => readingHigherThanRemembered(page, "x axis min", el("scatter plot viewer")));
      await session.step(205, "And the \"x axis min\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "x axis min", el("scatter plot viewer")));
      await session.step(206, "And the \"x axis max\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "x axis max", el("scatter plot viewer")));
      await session.step(207, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(208, "Then the \"x axis min\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis min", el("scatter plot viewer")));
      await session.step(209, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X Axis Type","linear"],["Zoom and Filter","filter by zoom"]]), [["X Axis Type","linear"],["Zoom and Filter","filter by zoom"]]);
      await session.step(212, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(213, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(214, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Mouse Drag = Select makes a plain drag select instead of pan (public#2089)", async () => {
      await session.step(218, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"no action\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "no action"));
      await session.step(219, "And user remembers the \"x axis min\" reading of scatter plot viewer", () => rememberReading(page, "x axis min", el("scatter plot viewer")));
      await session.step(220, "And user sets \"mouseDrag\" property of scatter plot viewer to \"Select\"", () => setProperty(page, "mouseDrag", el("scatter plot viewer"), "Select"));
      await session.step(221, "And user drags across the \"view\" area of scatter plot viewer", () => dragAcrossArea(page, "view", el("scatter plot viewer")));
      await session.step(222, "Then some rows should be selected", () => someSelected(page));
      await session.step(223, "And the \"x axis min\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis min", el("scatter plot viewer")));
      await session.step(224, "When user clears the row selection", () => clearSelection(page));
      await session.step(225, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["mouseDrag","Pan"],["Zoom and Filter","filter by zoom"]]), [["mouseDrag","Pan"],["Zoom and Filter","filter by zoom"]]);
      await session.step(228, "Then no rows should be selected", () => noneSelected(page));
      await session.step(229, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
