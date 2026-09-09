/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-select-and-zoom.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {dragLasso, dragRangeHandle} from '../../../bindings/scatter-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, filterPasses, filterPassesAll, filterTo, noneSelected, onlyOfSelected, resetFilter, selectFirstRows, selectedPassFilter, selectedRowCount, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, doubleClickArea, dragAcrossArea, dragDeselectionOverArea, dragSelectionBetweenAreas, dragSelectionOverArea, dragZoomOverArea, eventFired, hoverArea, listenFor, narrowerRange, noErrors, pickFromContextMenu, propertyShouldBe, readingAsRemembered, readingDiffers, readingHigher, readingIs, readingLower, readingSame, rememberRange, rememberReading, rememberedRange, repainted, setProperties, setProperty, showsRows, takeSnapshot, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot selection and viewport navigation", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-select-and-zoom.feature", import.meta.url);
  test("Scatter plot selection and viewport navigation", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]));
    await session.step(20, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(21, "And no rows should be selected", () => noneSelected(page));
    await run.scenario("A drag selects the markers inside it, a second drag adds, Control+Shift takes them out", async () => {
      await session.step(24, "Then \"Reset Selection On Background Click\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Reset Selection On Background Click", el("scatter plot viewer"), "true"));
      await session.step(25, "When user drags a selection box from the \"marker of row 11\" area to the \"marker of row 3\" area of scatter plot viewer", () => dragSelectionBetweenAreas(page, "marker of row 11", "marker of row 3", el("scatter plot viewer")));
      await session.step(26, "Then some rows should be selected", () => someSelected(page));
      await session.step(27, "And the \"rows selected\" reading of scatter plot viewer should be higher than before", () => readingHigher(page, "rows selected", el("scatter plot viewer")));
      await session.step(28, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(29, "Then the \"rows selected\" reading of scatter plot viewer should be higher than before", () => readingHigher(page, "rows selected", el("scatter plot viewer")));
      await session.step(30, "And every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(31, "When user drags a deselection box over the \"view\" area of scatter plot viewer", () => dragDeselectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(32, "Then the \"rows selected\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "rows selected", el("scatter plot viewer")));
      await session.step(33, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(34, "Then some rows should be selected", () => someSelected(page));
      await session.step(35, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(36, "Then no rows should be selected", () => noneSelected(page));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The selection survives a jitter change", async () => {
      await session.step(40, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Jitter Size","20"],["Jitter Size Y","15"]]));
      await session.step(43, "And user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(44, "Then some rows should be selected", () => someSelected(page));
      await session.step(45, "When user remembers the \"rows selected\" reading of scatter plot viewer", () => rememberReading(page, "rows selected", el("scatter plot viewer")));
      await session.step(46, "And user sets \"Jitter Size\" property of scatter plot viewer to \"30\"", () => setProperty(page, "Jitter Size", el("scatter plot viewer"), "30"));
      await session.step(47, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(48, "And the \"rows selected\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "rows selected", el("scatter plot viewer")));
      await session.step(49, "When user drags a deselection box over the \"view\" area of scatter plot viewer", () => dragDeselectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(50, "Then the \"rows selected\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "rows selected", el("scatter plot viewer")));
      await session.step(51, "When user clears the row selection", () => clearSelection(page));
      await session.step(52, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Jitter Size","0"],["Jitter Size Y","0"]]));
      await session.step(55, "Then no rows should be selected", () => noneSelected(page));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An Alt-drag zooms, a plain drag pans, the wheel zooms and Reset View comes home", async () => {
      await session.step(59, "When user selects the first 20 rows", () => selectFirstRows(page, 20));
      await session.step(60, "Then 20 rows should be selected", () => selectedRowCount(page, 20));
      await session.step(61, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"no action\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "no action"));
      await session.step(62, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(63, "And user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(64, "And user remembers the \"x axis span\" reading of scatter plot viewer", () => rememberReading(page, "x axis span", el("scatter plot viewer")));
      await session.step(65, "And user drags a zoom box over the \"view\" area of scatter plot viewer", () => dragZoomOverArea(page, "view", el("scatter plot viewer")));
      await session.step(66, "Then the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(67, "And scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(68, "And the \"x axis span\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "x axis span", el("scatter plot viewer")));
      await session.step(69, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(70, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(71, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(72, "And the \"x axis span\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis span", el("scatter plot viewer")));
      await session.step(73, "When user drags across the \"view\" area of scatter plot viewer", () => dragAcrossArea(page, "view", el("scatter plot viewer")));
      await session.step(74, "Then the \"x axis min\" reading of scatter plot viewer should differ from before", () => readingDiffers(page, "x axis min", el("scatter plot viewer")));
      await session.step(75, "And the \"x axis span\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "x axis span", el("scatter plot viewer")));
      await session.step(76, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(77, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(78, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(79, "Then the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(80, "And scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(81, "And the \"x axis span\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "x axis span", el("scatter plot viewer")));
      await session.step(82, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(83, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(84, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(85, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(86, "And user clears the row selection", () => clearSelection(page));
      await session.step(87, "Then no rows should be selected", () => noneSelected(page));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The X range slider narrows the axis and a double-click on empty space resets the view", async () => {
      await session.step(91, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"no action\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "no action"));
      await session.step(92, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(93, "And user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(94, "And user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(95, "And user drags the min handle of the x range slider of scatter plot viewer by 30 pixels", () => dragRangeHandle(page, "min", "x", 30));
      await session.step(96, "Then the \"x axis min\" reading of scatter plot viewer should be higher than before", () => readingHigher(page, "x axis min", el("scatter plot viewer")));
      await session.step(97, "And the \"x axis span\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "x axis span", el("scatter plot viewer")));
      await session.step(98, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(99, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(100, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(101, "Then the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(102, "And scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(103, "When user listens for \"d4-scatterplot-reset-view\" event on scatter plot viewer", () => listenFor(page, "d4-scatterplot-reset-view", el("scatter plot viewer")));
      await session.step(104, "And user hovers over the \"empty space\" area of scatter plot viewer", () => hoverArea(page, "empty space", el("scatter plot viewer")));
      await session.step(105, "Then the \"hovered row\" reading of scatter plot viewer should be 0", () => readingIs(page, "hovered row", el("scatter plot viewer"), 0));
      await session.step(106, "When user double-clicks on the \"empty space\" area of scatter plot viewer", () => doubleClickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(107, "Then \"d4-scatterplot-reset-view\" event should have fired on scatter plot viewer", () => eventFired(page, "d4-scatterplot-reset-view", el("scatter plot viewer")));
      await session.step(108, "And scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(109, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(110, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+A selects the rows the filter passes, Ctrl+Shift+A and Escape clear them", async () => {
      await session.step(113, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(114, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(115, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(116, "And user presses Control+a", () => pressKey(page, "Control+a"));
      await session.step(117, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(118, "And only rows where \"SEX\" is \"F\" should be selected", () => onlyOfSelected(page, "SEX", "F"));
      await session.step(119, "When user presses Control+Shift+a", () => pressKey(page, "Control+Shift+a"));
      await session.step(120, "Then no rows should be selected", () => noneSelected(page));
      await session.step(121, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(122, "Then some rows should be selected", () => someSelected(page));
      await session.step(123, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(124, "Then no rows should be selected", () => noneSelected(page));
      await session.step(125, "When user resets the filter", () => resetFilter(page));
      await session.step(126, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(127, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(128, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("H restores the viewport after a wheel zoom", async () => {
      await session.step(131, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(132, "And user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(133, "And user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(134, "Then the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(135, "And scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(136, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(137, "And user presses h", () => pressKey(page, "h"));
      await session.step(138, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(139, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(140, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("L turns on the Lasso Tool and a lasso selects the markers it encloses", async () => {
      await session.step(143, "Then \"Lasso Tool\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "Lasso Tool", el("scatter plot viewer"), "false"));
      await session.step(144, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(145, "And user presses l", () => pressKey(page, "l"));
      await session.step(146, "Then \"Lasso Tool\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Lasso Tool", el("scatter plot viewer"), "true"));
      await session.step(147, "When user drags a lasso over the \"view\" area of scatter plot viewer", () => dragLasso(page, "view"));
      await session.step(148, "Then some rows should be selected", () => someSelected(page));
      await session.step(149, "And every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(150, "When user presses l", () => pressKey(page, "l"));
      await session.step(151, "Then \"Lasso Tool\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "Lasso Tool", el("scatter plot viewer"), "false"));
      await session.step(152, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(153, "Then no rows should be selected", () => noneSelected(page));
      await session.step(154, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
