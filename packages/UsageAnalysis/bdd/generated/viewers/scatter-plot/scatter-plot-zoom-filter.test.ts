/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-zoom-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, public#552, public#1348, public#1744]
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
import {columnIncomplete} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, addCategoricalFilter, addRangeFilter, filterBetween, filterIsExactly, filterIsExactlyCategory, filterPasses, filterPassesAll, filterPassesFewer, openFilterPanel, removeColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, narrowerRange, noErrors, pickFromContextMenu, propertyShouldBe, readingAsRemembered, readingDiffers, readingLower, readingSame, rememberRange, rememberReading, rememberedRange, repainted, resizeTo, restoreSize, setProperties, setProperty, showsFewerRows, showsRows, viewerCount, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragRangeHandle, notShowsText, pickInColumnSelector, showsText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot zoom and filter synchronization", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-zoom-filter.feature", import.meta.url);
  test("Scatter plot zoom and filter synchronization", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot", "@realizes:public#552", "@realizes:public#1348", "@realizes:public#1744"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(27, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(28, "And all rows should pass the filter", () => filterPassesAll(page));
    await run.scenario("A zoom filters the table and Reset View gives it back", async () => {
      await session.step(31, "Then \"Zoom and Filter\" property of scatter plot viewer should be \"filter by zoom\"", () => propertyShouldBe(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(32, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(33, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(34, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(35, "And the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(36, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(37, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(38, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pack and zoom by filter does not leave the table filtered", async () => {
      await session.step(42, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(43, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(44, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"pack and zoom by filter\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "pack and zoom by filter"));
      await session.step(45, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(46, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(47, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(48, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(49, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(50, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Under zoom by filter an external filter narrows the viewport", async () => {
      await session.step(54, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"zoom by filter\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "zoom by filter"));
      await session.step(55, "And user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(56, "And user remembers the \"x axis span\" reading of scatter plot viewer", () => rememberReading(page, "x axis span", el("scatter plot viewer")));
      await session.step(57, "And user filters rows where \"WEIGHT\" is between 90.96 and 115.64", () => filterBetween(page, "WEIGHT", 90.96, 115.64));
      await session.step(58, "Then 212 rows should pass the filter", () => filterPasses(page, 212));
      await session.step(59, "And the \"x axis span\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "x axis span", el("scatter plot viewer")));
      await session.step(60, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(61, "And the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(62, "When user resets the filter", () => resetFilter(page));
      await session.step(63, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(64, "And the \"x axis span\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis span", el("scatter plot viewer")));
      await session.step(65, "And scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(66, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(67, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The filter panel's reset undoes what the zoom filtered", async () => {
      await session.step(70, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(71, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(72, "And status bar should not show the text \"Filtered:\"", () => notShowsText(page, el("status bar"), "Filtered:"));
      await session.step(73, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(74, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(75, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(76, "And the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(77, "And status bar should show the text \"Filtered:\"", () => showsText(page, el("status bar"), "Filtered:"));
      await session.step(78, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(79, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(80, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(81, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(82, "And status bar should not show the text \"Filtered:\"", () => notShowsText(page, el("status bar"), "Filtered:"));
      await session.step(83, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(84, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic axis alone filters nothing; Filter Out Invalid drops the non-positive rows", async () => {
      await session.step(88, "When user adds a calculated column \"AGE_SHIFT\" with formula \"${AGE} - 40\"", () => addCalculated(page, "AGE_SHIFT", "${AGE} - 40"));
      await session.step(89, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","WEIGHT"],["Y","AGE_SHIFT"]]), [["X","WEIGHT"],["Y","AGE_SHIFT"]]);
      await session.step(92, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
      await session.step(93, "And \"Filter Out Invalid\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "Filter Out Invalid", el("scatter plot viewer"), "false"));
      await session.step(94, "When user sets \"Y Axis Type\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "Y Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(95, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(96, "When user sets \"Filter Out Invalid\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Filter Out Invalid", el("scatter plot viewer"), "true"));
      await session.step(97, "Then 635 rows should pass the filter", () => filterPasses(page, 635));
      await session.step(98, "And the filter should pass exactly the rows where \"AGE\" is between 40.5 and 100", () => filterIsExactly(page, "AGE", 40.5, 100));
      await session.step(99, "When user sets \"Filter Out Invalid\" property of scatter plot viewer to \"false\"", () => setProperty(page, "Filter Out Invalid", el("scatter plot viewer"), "false"));
      await session.step(100, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(101, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y Axis Type","linear"],["Filter Out Invalid","true"]]), [["Y Axis Type","linear"],["Filter Out Invalid","true"]]);
      await session.step(104, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(105, "And scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
      await session.step(106, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Filter Out Invalid","false"],["X","WEIGHT"],["Y","HEIGHT"]]), [["Filter Out Invalid","false"],["X","WEIGHT"],["Y","HEIGHT"]]);
      await session.step(110, "And user removes \"AGE_SHIFT\" column", () => removeColumn(page, "AGE_SHIFT"));
      await session.step(111, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(112, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(113, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The axis type switched from the menu on a datetime axis leaves the filter as it was", async () => {
      await session.step(117, "When user picks \"STARTED\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "STARTED", "x", el("scatter plot viewer")));
      await session.step(118, "And user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
      await session.step(119, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(120, "And \"X Axis Type\" property of scatter plot viewer should be \"linear\"", () => propertyShouldBe(page, "X Axis Type", el("scatter plot viewer"), "linear"));
      await session.step(121, "When user picks \"Properties... > X Axis > X Axis Type > Logarithmic\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Properties... > X Axis > X Axis Type > Logarithmic", el("scatter plot viewer")));
      await session.step(122, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(123, "And \"X Axis Type\" property of scatter plot viewer should be \"logarithmic\"", () => propertyShouldBe(page, "X Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(124, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(125, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
      await session.step(126, "When user picks \"Properties... > X Axis > X Axis Type > Linear\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Properties... > X Axis > X Axis Type > Linear", el("scatter plot viewer")));
      await session.step(127, "Then \"X Axis Type\" property of scatter plot viewer should be \"linear\"", () => propertyShouldBe(page, "X Axis Type", el("scatter plot viewer"), "linear"));
      await session.step(128, "When user adds a categorical filter on \"SEX\" keeping \"F, M\"", () => addCategoricalFilter(page, "SEX", "F, M"));
      await session.step(129, "And user sets \"X\" property of scatter plot viewer to \"WEIGHT\"", () => setProperty(page, "X", el("scatter plot viewer"), "WEIGHT"));
      await session.step(130, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(131, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A large jitter on a logarithmic axis filters no row", async () => {
      await session.step(135, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["jitterSize","16"],["jitterSizeY","17"]]), [["jitterSize","16"],["jitterSizeY","17"]]);
      await session.step(138, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(139, "When user sets \"Y Axis Type\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "Y Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(140, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(141, "And the \"y axis min\" reading of scatter plot viewer should differ from before", () => readingDiffers(page, "y axis min", el("scatter plot viewer")));
      await session.step(142, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(143, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(144, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y Axis Type","linear"],["jitterSize","0"],["jitterSizeY","0"]]), [["Y Axis Type","linear"],["jitterSize","0"],["jitterSizeY","0"]]);
      await session.step(148, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(149, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(150, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A zoom under Zoom and Filter = no action survives filtering in the filter panel (public#552)", async () => {
      await session.step(154, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Zoom and Filter","no action"],["axesFollowFilter","false"]]), [["Zoom and Filter","no action"],["axesFollowFilter","false"]]);
      await session.step(157, "And user drags the min handle of the \"x\" range slider of scatter plot viewer by 30 pixels", () => dragRangeHandle(page, "min", "x", el("scatter plot viewer"), 30));
      await session.step(158, "Then scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(159, "When user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(160, "And user opens the filter panel", () => openFilterPanel(page));
      await session.step(161, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(162, "When user adds a range filter on \"AGE\" from 20 to 40", () => addRangeFilter(page, "AGE", 20, 40));
      await session.step(163, "Then 359 rows should pass the filter", () => filterPasses(page, 359));
      await session.step(164, "And scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(165, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(166, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(167, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(168, "And scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(169, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Zoom and Filter","filter by zoom"],["axesFollowFilter","true"]]), [["Zoom and Filter","filter by zoom"],["axesFollowFilter","true"]]);
      await session.step(172, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(173, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(174, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(175, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A plot shrunk almost to nothing filters no row, and the panel's reset gives every row back (public#1348)", async () => {
      await session.step(179, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(180, "And user resizes scatter plot viewer to 40 by 40", () => resizeTo(page, el("scatter plot viewer"), 40, 40));
      await session.step(181, "And user adds a range filter on \"WEIGHT\" from 100 to 170", () => addRangeFilter(page, "WEIGHT", 100, 170));
      await session.step(182, "Then 143 rows should pass the filter", () => filterPasses(page, 143));
      await session.step(183, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(184, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(185, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(186, "When user restores the size of scatter plot viewer", () => restoreSize(page, el("scatter plot viewer")));
      await session.step(187, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(188, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(189, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic axis on an all-empty column under filter by zoom filters no row (public#1744)", async () => {
      await session.step(193, "When user adds a calculated column \"NO_VALUES\" with formula \"if(${AGE} > 1000, ${WEIGHT}, null)\"", () => addCalculated(page, "NO_VALUES", "if(${AGE} > 1000, ${WEIGHT}, null)"));
      await session.step(194, "Then \"NO_VALUES\" column should have missing values", () => columnIncomplete(page, "NO_VALUES"));
      await session.step(195, "And \"Zoom and Filter\" property of scatter plot viewer should be \"filter by zoom\"", () => propertyShouldBe(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(196, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y Axis Type","logarithmic"],["Y","NO_VALUES"]]), [["Y Axis Type","logarithmic"],["Y","NO_VALUES"]]);
      await session.step(199, "Then scatter plot viewer should show 0 rows", () => showsRows(page, el("scatter plot viewer"), 0));
      await session.step(200, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(201, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y","HEIGHT"],["Y Axis Type","linear"]]), [["Y","HEIGHT"],["Y Axis Type","linear"]]);
      await session.step(204, "And user removes \"NO_VALUES\" column", () => removeColumn(page, "NO_VALUES"));
      await session.step(205, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(206, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(207, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(208, "And no errors should have been logged", () => noErrors(page));
      await session.step(209, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(210, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(211, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
      await session.step(212, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(213, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(214, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
