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
import {addViewerWith, narrowerRange, noErrors, pickFromContextMenu, propertyShouldBe, readingAsRemembered, readingDiffers, readingLower, readingSame, rememberRange, rememberReading, rememberedRange, repainted, resizeTo, restoreSize, setProperties, setProperty, showsFewerRows, showsRows, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragRangeHandle, notShowsText, pickInColumnSelector, showsText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot zoom and filter synchronization", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-zoom-filter.feature", import.meta.url);
  test("Scatter plot zoom and filter synchronization", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot", "@realizes:public#552", "@realizes:public#1348", "@realizes:public#1744"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(25, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(26, "And all rows should pass the filter", () => filterPassesAll(page));
    await run.scenario("A zoom filters the table and Reset View gives it back", async () => {
      await session.step(29, "Then \"Zoom and Filter\" property of scatter plot viewer should be \"filter by zoom\"", () => propertyShouldBe(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(30, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(31, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(32, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(33, "And the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(34, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(35, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(36, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pack and zoom by filter does not leave the table filtered", async () => {
      await session.step(40, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(41, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(42, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"pack and zoom by filter\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "pack and zoom by filter"));
      await session.step(43, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(44, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(45, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(46, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(47, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(48, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Under zoom by filter an external filter narrows the viewport", async () => {
      await session.step(52, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"zoom by filter\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "zoom by filter"));
      await session.step(53, "And user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(54, "And user remembers the \"x axis span\" reading of scatter plot viewer", () => rememberReading(page, "x axis span", el("scatter plot viewer")));
      await session.step(55, "And user filters rows where \"WEIGHT\" is between 90.96 and 115.64", () => filterBetween(page, "WEIGHT", 90.96, 115.64));
      await session.step(56, "Then 212 rows should pass the filter", () => filterPasses(page, 212));
      await session.step(57, "And the \"x axis span\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "x axis span", el("scatter plot viewer")));
      await session.step(58, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(59, "And the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(60, "When user resets the filter", () => resetFilter(page));
      await session.step(61, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(62, "And the \"x axis span\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis span", el("scatter plot viewer")));
      await session.step(63, "And scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(64, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(65, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The filter panel's reset undoes what the zoom filtered", async () => {
      await session.step(68, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(69, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(70, "And status bar should not show the text \"Filtered:\"", () => notShowsText(page, el("status bar"), "Filtered:"));
      await session.step(71, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(72, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(73, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(74, "And the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(75, "And status bar should show the text \"Filtered:\"", () => showsText(page, el("status bar"), "Filtered:"));
      await session.step(76, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(77, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(78, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(79, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(80, "And status bar should not show the text \"Filtered:\"", () => notShowsText(page, el("status bar"), "Filtered:"));
      await session.step(81, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(82, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic axis alone filters nothing; Filter Out Invalid drops the non-positive rows", async () => {
      await session.step(86, "When user adds a calculated column \"AGE_SHIFT\" with formula \"${AGE} - 40\"", () => addCalculated(page, "AGE_SHIFT", "${AGE} - 40"));
      await session.step(87, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","WEIGHT"],["Y","AGE_SHIFT"]]), [["X","WEIGHT"],["Y","AGE_SHIFT"]]);
      await session.step(90, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
      await session.step(91, "And \"Filter Out Invalid\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "Filter Out Invalid", el("scatter plot viewer"), "false"));
      await session.step(92, "When user sets \"Y Axis Type\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "Y Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(93, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(94, "When user sets \"Filter Out Invalid\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Filter Out Invalid", el("scatter plot viewer"), "true"));
      await session.step(95, "Then 635 rows should pass the filter", () => filterPasses(page, 635));
      await session.step(96, "And the filter should pass exactly the rows where \"AGE\" is between 40.5 and 100", () => filterIsExactly(page, "AGE", 40.5, 100));
      await session.step(97, "When user sets \"Filter Out Invalid\" property of scatter plot viewer to \"false\"", () => setProperty(page, "Filter Out Invalid", el("scatter plot viewer"), "false"));
      await session.step(98, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(99, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y Axis Type","linear"],["Filter Out Invalid","true"]]), [["Y Axis Type","linear"],["Filter Out Invalid","true"]]);
      await session.step(102, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(103, "And scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
      await session.step(104, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Filter Out Invalid","false"],["X","WEIGHT"],["Y","HEIGHT"]]), [["Filter Out Invalid","false"],["X","WEIGHT"],["Y","HEIGHT"]]);
      await session.step(108, "And user removes \"AGE_SHIFT\" column", () => removeColumn(page, "AGE_SHIFT"));
      await session.step(109, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(110, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(111, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The axis type switched from the menu on a datetime axis leaves the filter as it was", async () => {
      await session.step(115, "When user picks \"STARTED\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "STARTED", "x", el("scatter plot viewer")));
      await session.step(116, "And user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
      await session.step(117, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(118, "And \"X Axis Type\" property of scatter plot viewer should be \"linear\"", () => propertyShouldBe(page, "X Axis Type", el("scatter plot viewer"), "linear"));
      await session.step(119, "When user picks \"Properties... > X Axis > X Axis Type > Logarithmic\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Properties... > X Axis > X Axis Type > Logarithmic", el("scatter plot viewer")));
      await session.step(120, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(121, "And \"X Axis Type\" property of scatter plot viewer should be \"logarithmic\"", () => propertyShouldBe(page, "X Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(122, "And 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(123, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
      await session.step(124, "When user picks \"Properties... > X Axis > X Axis Type > Linear\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Properties... > X Axis > X Axis Type > Linear", el("scatter plot viewer")));
      await session.step(125, "Then \"X Axis Type\" property of scatter plot viewer should be \"linear\"", () => propertyShouldBe(page, "X Axis Type", el("scatter plot viewer"), "linear"));
      await session.step(126, "When user adds a categorical filter on \"SEX\" keeping \"F, M\"", () => addCategoricalFilter(page, "SEX", "F, M"));
      await session.step(127, "And user sets \"X\" property of scatter plot viewer to \"WEIGHT\"", () => setProperty(page, "X", el("scatter plot viewer"), "WEIGHT"));
      await session.step(128, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(129, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(130, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A large jitter on a logarithmic axis filters no row", async () => {
      await session.step(133, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["jitterSize","16"],["jitterSizeY","17"]]), [["jitterSize","16"],["jitterSizeY","17"]]);
      await session.step(136, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(137, "When user sets \"Y Axis Type\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "Y Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(138, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(139, "And the \"y axis min\" reading of scatter plot viewer should differ from before", () => readingDiffers(page, "y axis min", el("scatter plot viewer")));
      await session.step(140, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(141, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(142, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y Axis Type","linear"],["jitterSize","0"],["jitterSizeY","0"]]), [["Y Axis Type","linear"],["jitterSize","0"],["jitterSizeY","0"]]);
      await session.step(146, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(147, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(148, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A zoom under Zoom and Filter = no action survives filtering in the filter panel (public#552)", async () => {
      await session.step(152, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Zoom and Filter","no action"],["axesFollowFilter","false"]]), [["Zoom and Filter","no action"],["axesFollowFilter","false"]]);
      await session.step(155, "And user drags the min handle of the \"x\" range slider of scatter plot viewer by 30 pixels", () => dragRangeHandle(page, "min", "x", el("scatter plot viewer"), 30));
      await session.step(156, "Then scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(157, "When user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(158, "And user opens the filter panel", () => openFilterPanel(page));
      await session.step(159, "Then scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(160, "When user adds a range filter on \"AGE\" from 20 to 40", () => addRangeFilter(page, "AGE", 20, 40));
      await session.step(161, "Then 359 rows should pass the filter", () => filterPasses(page, 359));
      await session.step(162, "And scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(163, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(164, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(165, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(166, "And scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(167, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Zoom and Filter","filter by zoom"],["axesFollowFilter","true"]]), [["Zoom and Filter","filter by zoom"],["axesFollowFilter","true"]]);
      await session.step(170, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(171, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(172, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(173, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A plot shrunk almost to nothing filters no row, and the panel's reset gives every row back (public#1348)", async () => {
      await session.step(177, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(178, "And user resizes scatter plot viewer to 40 by 40", () => resizeTo(page, el("scatter plot viewer"), 40, 40));
      await session.step(179, "And user adds a range filter on \"WEIGHT\" from 100 to 170", () => addRangeFilter(page, "WEIGHT", 100, 170));
      await session.step(180, "Then 143 rows should pass the filter", () => filterPasses(page, 143));
      await session.step(181, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(182, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(183, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(184, "When user restores the size of scatter plot viewer", () => restoreSize(page, el("scatter plot viewer")));
      await session.step(185, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(186, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(187, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic axis on an all-empty column under filter by zoom filters no row (public#1744)", async () => {
      await session.step(191, "When user adds a calculated column \"NO_VALUES\" with formula \"if(${AGE} > 1000, ${WEIGHT}, null)\"", () => addCalculated(page, "NO_VALUES", "if(${AGE} > 1000, ${WEIGHT}, null)"));
      await session.step(192, "Then \"NO_VALUES\" column should have missing values", () => columnIncomplete(page, "NO_VALUES"));
      await session.step(193, "And \"Zoom and Filter\" property of scatter plot viewer should be \"filter by zoom\"", () => propertyShouldBe(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(194, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y Axis Type","logarithmic"],["Y","NO_VALUES"]]), [["Y Axis Type","logarithmic"],["Y","NO_VALUES"]]);
      await session.step(197, "Then scatter plot viewer should show 0 rows", () => showsRows(page, el("scatter plot viewer"), 0));
      await session.step(198, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(199, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y","HEIGHT"],["Y Axis Type","linear"]]), [["Y","HEIGHT"],["Y Axis Type","linear"]]);
      await session.step(202, "And user removes \"NO_VALUES\" column", () => removeColumn(page, "NO_VALUES"));
      await session.step(203, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(204, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(205, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(206, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
