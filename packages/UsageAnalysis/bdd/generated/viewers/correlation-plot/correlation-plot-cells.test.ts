/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/correlation-plot/correlation-plot-cells.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.correlation-plot]
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
import {clearField, clickOn, shouldBe, shouldContainText, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueBetween, hasColumn, maxInRow, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, closeContextMenu, doubleClickArea, hasArea, hasNoArea, hoverArea, menuLists, noErrors, oneTooltip, painted, pickFromAreaContextMenu, pointerAway, propertyShouldBe, readingAsRemembered, readingBetween, readingIs, readingLower, readingReads, rememberReading, resizeTo, restoreSize, rightClickArea, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Correlation plot — a cell as a click, hover and menu target", () => {
  const session = feature(test, "features/viewers/correlation-plot/correlation-plot-cells.feature", import.meta.url);
  test("Correlation plot — a cell as a click, hover and menu target", {tag: ["@journey", "@viewers", "@realizes:viewers.correlation-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 12, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a correlation plot viewer", () => addViewer(page, "correlation plot"));
    await session.step(21, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
    await session.step(22, "And the \"show pearson r\" reading of correlation plot viewer should be \"true\"", () => readingReads(page, "show pearson r", el("correlation plot viewer"), "true"));
    await session.step(23, "And the \"last clicked cell\" reading of correlation plot viewer should be \"\"", () => readingReads(page, "last clicked cell", el("correlation plot viewer"), ""));
    await session.step(24, "And correlation plot viewer should be painted", () => painted(page, el("correlation plot viewer")));
    await run.scenario("A click records which cell it hit and the value that cell holds", async () => {
      await session.step(27, "When user clicks on the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => clickArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(28, "Then the \"last clicked cell\" reading of correlation plot viewer should be \"HEIGHT x AGE\"", () => readingReads(page, "last clicked cell", el("correlation plot viewer"), "HEIGHT x AGE"));
      await session.step(29, "And the \"last clicked value\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "last clicked value", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(30, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(31, "And the \"current column\" reading of correlation plot viewer should be \"HEIGHT\"", () => readingReads(page, "current column", el("correlation plot viewer"), "HEIGHT"));
      await session.step(32, "When user clicks on the \"cell WEIGHT x AGE\" area of correlation plot viewer", () => clickArea(page, "cell WEIGHT x AGE", el("correlation plot viewer")));
      await session.step(33, "Then the \"last clicked cell\" reading of correlation plot viewer should be \"WEIGHT x AGE\"", () => readingReads(page, "last clicked cell", el("correlation plot viewer"), "WEIGHT x AGE"));
      await session.step(34, "And the \"last clicked value\" reading of correlation plot viewer should be between 0.0647 and 0.0649", () => readingBetween(page, "last clicked value", el("correlation plot viewer"), 0.0647, 0.0649));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The click puts that pair in the context panel with a scatter plot of its own", async () => {
      await session.step(38, "When user clicks on the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => clickArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(39, "Then context panel should contain the text \"HEIGHT vs AGE\"", () => shouldContainText(page, el("context panel"), "HEIGHT vs AGE"));
      await session.step(40, "And \"Scatter plot\" section in context panel should be present", () => shouldBe(page, el("\"Scatter plot\" section in context panel"), "present"));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A double-click opens the pair as a real scatter plot, and closing it puts the view back", async () => {
      await session.step(44, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(45, "When user double-clicks on the \"cell WEIGHT x AGE\" area of correlation plot viewer", () => doubleClickArea(page, "cell WEIGHT x AGE", el("correlation plot viewer")));
      await session.step(46, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(47, "And \"xColumnName\" property of scatter plot viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "WEIGHT"));
      await session.step(48, "And \"yColumnName\" property of scatter plot viewer should be \"AGE\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer"), "AGE"));
      await session.step(49, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
      await session.step(50, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ignore Double Click suppresses the same gesture on the same cell", async () => {
      await session.step(54, "When user sets \"ignoreDoubleClick\" property of correlation plot viewer to \"true\"", () => setProperty(page, "ignoreDoubleClick", el("correlation plot viewer"), "true"));
      await session.step(55, "And user double-clicks on the \"cell WEIGHT x AGE\" area of correlation plot viewer", () => doubleClickArea(page, "cell WEIGHT x AGE", el("correlation plot viewer")));
      await session.step(56, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(57, "When user sets \"ignoreDoubleClick\" property of correlation plot viewer to \"false\"", () => setProperty(page, "ignoreDoubleClick", el("correlation plot viewer"), "false"));
      await session.step(58, "And user double-clicks on the \"cell WEIGHT x AGE\" area of correlation plot viewer", () => doubleClickArea(page, "cell WEIGHT x AGE", el("correlation plot viewer")));
      await session.step(59, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(60, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
      await session.step(61, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The cell menu is the grid's, with the plot's own items on it", async () => {
      await session.step(71, "When user right-clicks on the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => rightClickArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(72, "Then the open menu should list \"Show Coefficients\"", () => menuLists(page, "Show Coefficients"));
      await session.step(73, "And the open menu should list \"Open as table\"", () => menuLists(page, "Open as table"));
      await session.step(74, "And the open menu should list \"Columns\"", () => menuLists(page, "Columns"));
      await session.step(75, "And the open menu should list \"Tooltip\"", () => menuLists(page, "Tooltip"));
      await session.step(76, "And the open menu should list \"Grid\"", () => menuLists(page, "Grid"));
      await session.step(77, "And the open menu should list \"'AGE' column\"", () => menuLists(page, "'AGE' column"));
      await session.step(78, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(79, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Grid > Order or Hide Columns from the cell menu searches the matrix's columns and closes, the name column still pinned (GROK-9310)", async () => {
      await session.step(82, "When user picks \"Grid > Order or Hide Columns...\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Grid > Order or Hide Columns...", "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(83, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
      await session.step(84, "When user remembers the \"rows shown\" reading of Grid viewer in Order or Hide Columns dialog", () => rememberReading(page, "rows shown", el("Grid viewer in Order or Hide Columns dialog")));
      await session.step(85, "And user types \"WEIGHT\" into \"Search\" input in Order or Hide Columns dialog", () => typeInto(page, "WEIGHT", el("\"Search\" input in Order or Hide Columns dialog")));
      await session.step(86, "Then the \"rows shown\" reading of Grid viewer in Order or Hide Columns dialog should be 1", () => readingIs(page, "rows shown", el("Grid viewer in Order or Hide Columns dialog"), 1));
      await session.step(87, "When user clears \"Search\" input in Order or Hide Columns dialog", () => clearField(page, el("\"Search\" input in Order or Hide Columns dialog")));
      await session.step(88, "Then the \"rows shown\" reading of Grid viewer in Order or Hide Columns dialog should be as remembered", () => readingAsRemembered(page, "rows shown", el("Grid viewer in Order or Hide Columns dialog")));
      await session.step(89, "When user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
      await session.step(90, "Then Order or Hide Columns dialog should be hidden", () => shouldBe(page, el("Order or Hide Columns dialog"), "hidden"));
      await session.step(91, "And correlation plot viewer should have a \"pinned band\" area", () => hasArea(page, el("correlation plot viewer"), "pinned band"));
      await session.step(92, "And correlation plot viewer should have a \"row header HEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "row header HEIGHT"));
      await session.step(93, "And the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Coefficients from the cell menu does what the property does", async () => {
      await session.step(97, "When user picks \"Show Coefficients\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Show Coefficients", "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(98, "Then the \"show pearson r\" reading of correlation plot viewer should be \"false\"", () => readingReads(page, "show pearson r", el("correlation plot viewer"), "false"));
      await session.step(99, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), ""));
      await session.step(100, "And the \"cell width\" reading of correlation plot viewer should be 20", () => readingIs(page, "cell width", el("correlation plot viewer"), 20));
      await session.step(101, "When user picks \"Show Coefficients\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Show Coefficients", "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(102, "Then the \"show pearson r\" reading of correlation plot viewer should be \"true\"", () => readingReads(page, "show pearson r", el("correlation plot viewer"), "true"));
      await session.step(103, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
      await session.step(104, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Open as table hands over the matrix, coefficients and all (GROK-19053)", async () => {
      await session.step(107, "When user picks \"Open as table\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Open as table", "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(108, "Then table \"corr\" should be open", () => tableOpen(page, "corr"));
      await session.step(109, "And table \"corr\" should have 4 rows", () => tableRows(page, "corr", 4));
      await session.step(110, "When user switches to the \"corr\" table view", () => switchTableView(page, "corr"));
      await session.step(111, "Then the table should have a column \"__name\"", () => hasColumn(page, "__name"));
      await session.step(112, "And the table should have a column \"AGE\"", () => hasColumn(page, "AGE"));
      await session.step(113, "And the value of \"__name\" column in row 2 should be \"HEIGHT\"", () => valueInRow(page, "__name", 2, "HEIGHT"));
      await session.step(114, "And every value of \"AGE\" column should lie between -1 and 1", () => everyValueBetween(page, "AGE", -1, 1));
      await session.step(115, "And \"AGE\" column should have its maximum in row 3", () => maxInRow(page, "AGE", 3));
      await session.step(116, "When user closes the current view", () => closeCurrentView(page));
      await session.step(117, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(118, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(119, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Narrowed to a fraction of its width, the matrix scrolls and keeps the name column", async () => {
      await session.step(122, "Then correlation plot viewer should not have a \"horz scroll\" area", () => hasNoArea(page, el("correlation plot viewer"), "horz scroll"));
      await session.step(123, "And the \"columns shown\" reading of correlation plot viewer should be 6", () => readingIs(page, "columns shown", el("correlation plot viewer"), 6));
      await session.step(124, "When user resizes correlation plot viewer to 190 by 400", () => resizeTo(page, el("correlation plot viewer"), 190, 400));
      await session.step(125, "Then correlation plot viewer should have a \"horz scroll\" area", () => hasArea(page, el("correlation plot viewer"), "horz scroll"));
      await session.step(126, "And the \"columns shown\" reading of correlation plot viewer should be lower than before", () => readingLower(page, "columns shown", el("correlation plot viewer")));
      await session.step(127, "And correlation plot viewer should have a \"row header HEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "row header HEIGHT"));
      await session.step(128, "And correlation plot viewer should have a \"pinned band\" area", () => hasArea(page, el("correlation plot viewer"), "pinned band"));
      await session.step(129, "When user restores the size of correlation plot viewer", () => restoreSize(page, el("correlation plot viewer")));
      await session.step(130, "Then correlation plot viewer should not have a \"horz scroll\" area", () => hasNoArea(page, el("correlation plot viewer"), "horz scroll"));
      await session.step(131, "And the \"columns shown\" reading of correlation plot viewer should be 6", () => readingIs(page, "columns shown", el("correlation plot viewer"), 6));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The cell tooltip names the pair, its coefficient and the plot drawn inside it (GROK-20125)", async () => {
      await session.step(135, "When user hovers over the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => hoverArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(136, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(137, "And tooltip should contain the text \"Pearson R: -0.235\"", () => shouldContainText(page, el("tooltip"), "Pearson R: -0.235"));
      await session.step(138, "And the \"tooltip plot x column\" reading of correlation plot viewer should be \"HEIGHT\"", () => readingReads(page, "tooltip plot x column", el("correlation plot viewer"), "HEIGHT"));
      await session.step(139, "And the \"tooltip plot y column\" reading of correlation plot viewer should be \"AGE\"", () => readingReads(page, "tooltip plot y column", el("correlation plot viewer"), "AGE"));
      await session.step(140, "When user hovers over the \"cell WEIGHT x HEIGHT\" area of correlation plot viewer", () => hoverArea(page, "cell WEIGHT x HEIGHT", el("correlation plot viewer")));
      await session.step(141, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(142, "And tooltip should contain the text \"Pearson R: 0.412\"", () => shouldContainText(page, el("tooltip"), "Pearson R: 0.412"));
      await session.step(143, "And the \"tooltip plot x column\" reading of correlation plot viewer should be \"WEIGHT\"", () => readingReads(page, "tooltip plot x column", el("correlation plot viewer"), "WEIGHT"));
      await session.step(144, "And the \"tooltip plot y column\" reading of correlation plot viewer should be \"HEIGHT\"", () => readingReads(page, "tooltip plot y column", el("correlation plot viewer"), "HEIGHT"));
      await session.step(145, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(146, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The pinned name column hovers the column's statistics, not a coefficient", async () => {
      await session.step(149, "When user hovers over the \"row header HEIGHT\" area of correlation plot viewer", () => hoverArea(page, "row header HEIGHT", el("correlation plot viewer")));
      await session.step(150, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(151, "And tooltip should contain the text \"min:\"", () => shouldContainText(page, el("tooltip"), "min:"));
      await session.step(152, "And tooltip should contain the text \"nulls: 128\"", () => shouldContainText(page, el("tooltip"), "nulls: 128"));
      await session.step(153, "And tooltip should not contain the text \"Pearson R\"", () => shouldNotContainText(page, el("tooltip"), "Pearson R"));
      await session.step(154, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(155, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Tooltip off hides the cell tooltip", async () => {
      await session.step(160, "When user sets \"showTooltip\" property of correlation plot viewer to \"false\"", () => setProperty(page, "showTooltip", el("correlation plot viewer"), "false"));
      await session.step(161, "And user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(162, "And user hovers over the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => hoverArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(163, "Then tooltip should not contain the text \"Pearson R\"", () => shouldNotContainText(page, el("tooltip"), "Pearson R"));
      await session.step(164, "When user sets \"showTooltip\" property of correlation plot viewer to \"true\"", () => setProperty(page, "showTooltip", el("correlation plot viewer"), "true"));
      await session.step(165, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
