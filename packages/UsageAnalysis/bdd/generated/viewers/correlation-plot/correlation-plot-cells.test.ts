/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/correlation-plot/correlation-plot-cells.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.correlation-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueBetween, hasColumn, maxInRow, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, closeContextMenu, doubleClickArea, hasArea, hasNoArea, hoverArea, menuLists, noErrors, oneTooltip, painted, pickFromAreaContextMenu, pointerAway, propertyShouldBe, readingBetween, readingIs, readingLower, readingReads, resizeTo, restoreSize, rightClickArea, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Correlation plot — a cell as a click, hover and menu target", () => {
  const session = feature(test, "features/viewers/correlation-plot/correlation-plot-cells.feature", import.meta.url);
  test("Correlation plot — a cell as a click, hover and menu target", {tag: ["@journey", "@viewers", "@realizes:viewers.correlation-plot", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a correlation plot viewer", () => addViewer(page, "correlation plot"));
    await session.step(20, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
    await session.step(21, "And the \"show pearson r\" reading of correlation plot viewer should be \"true\"", () => readingReads(page, "show pearson r", el("correlation plot viewer"), "true"));
    await session.step(22, "And the \"last clicked cell\" reading of correlation plot viewer should be \"\"", () => readingReads(page, "last clicked cell", el("correlation plot viewer"), ""));
    await session.step(23, "And correlation plot viewer should be painted", () => painted(page, el("correlation plot viewer")));
    await run.scenario("A click records which cell it hit and the value that cell holds", async () => {
      await session.step(26, "When user clicks on the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => clickArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(27, "Then the \"last clicked cell\" reading of correlation plot viewer should be \"HEIGHT x AGE\"", () => readingReads(page, "last clicked cell", el("correlation plot viewer"), "HEIGHT x AGE"));
      await session.step(28, "And the \"last clicked value\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "last clicked value", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(29, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(30, "And the \"current column\" reading of correlation plot viewer should be \"HEIGHT\"", () => readingReads(page, "current column", el("correlation plot viewer"), "HEIGHT"));
      await session.step(31, "When user clicks on the \"cell WEIGHT x AGE\" area of correlation plot viewer", () => clickArea(page, "cell WEIGHT x AGE", el("correlation plot viewer")));
      await session.step(32, "Then the \"last clicked cell\" reading of correlation plot viewer should be \"WEIGHT x AGE\"", () => readingReads(page, "last clicked cell", el("correlation plot viewer"), "WEIGHT x AGE"));
      await session.step(33, "And the \"last clicked value\" reading of correlation plot viewer should be between 0.0647 and 0.0649", () => readingBetween(page, "last clicked value", el("correlation plot viewer"), 0.0647, 0.0649));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The click puts that pair in the context panel with a scatter plot of its own", async () => {
      await session.step(37, "When user clicks on the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => clickArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(38, "Then context panel should contain the text \"HEIGHT vs AGE\"", () => shouldContainText(page, el("context panel"), "HEIGHT vs AGE"));
      await session.step(39, "And \"Scatter plot\" section in context panel should be present", () => shouldBe(page, el("\"Scatter plot\" section in context panel"), "present"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A double-click opens the pair as a real scatter plot, and closing it puts the view back", async () => {
      await session.step(43, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(44, "When user double-clicks on the \"cell WEIGHT x AGE\" area of correlation plot viewer", () => doubleClickArea(page, "cell WEIGHT x AGE", el("correlation plot viewer")));
      await session.step(45, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(46, "And \"xColumnName\" property of scatter plot viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "xColumnName", el("scatter plot viewer"), "WEIGHT"));
      await session.step(47, "And \"yColumnName\" property of scatter plot viewer should be \"AGE\"", () => propertyShouldBe(page, "yColumnName", el("scatter plot viewer"), "AGE"));
      await session.step(48, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
      await session.step(49, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ignore Double Click suppresses the same gesture on the same cell", async () => {
      await session.step(53, "When user sets \"ignoreDoubleClick\" property of correlation plot viewer to \"true\"", () => setProperty(page, "ignoreDoubleClick", el("correlation plot viewer"), "true"));
      await session.step(54, "And user double-clicks on the \"cell WEIGHT x AGE\" area of correlation plot viewer", () => doubleClickArea(page, "cell WEIGHT x AGE", el("correlation plot viewer")));
      await session.step(55, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(56, "When user sets \"ignoreDoubleClick\" property of correlation plot viewer to \"false\"", () => setProperty(page, "ignoreDoubleClick", el("correlation plot viewer"), "false"));
      await session.step(57, "And user double-clicks on the \"cell WEIGHT x AGE\" area of correlation plot viewer", () => doubleClickArea(page, "cell WEIGHT x AGE", el("correlation plot viewer")));
      await session.step(58, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(59, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
      await session.step(60, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The cell menu is the grid's, with the plot's own items on it", async () => {
      await session.step(70, "When user right-clicks on the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => rightClickArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(71, "Then the open menu should list \"Show Pearson R\"", () => menuLists(page, "Show Pearson R"));
      await session.step(72, "And the open menu should list \"Open as table\"", () => menuLists(page, "Open as table"));
      await session.step(73, "And the open menu should list \"Columns\"", () => menuLists(page, "Columns"));
      await session.step(74, "And the open menu should list \"Tooltip\"", () => menuLists(page, "Tooltip"));
      await session.step(75, "And the open menu should list \"Grid\"", () => menuLists(page, "Grid"));
      await session.step(76, "And the open menu should list \"'AGE' column\"", () => menuLists(page, "'AGE' column"));
      await session.step(77, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(78, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Pearson R from the cell menu does what the property does", async () => {
      await session.step(81, "When user picks \"Show Pearson R\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Show Pearson R", "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(82, "Then the \"show pearson r\" reading of correlation plot viewer should be \"false\"", () => readingReads(page, "show pearson r", el("correlation plot viewer"), "false"));
      await session.step(83, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), ""));
      await session.step(84, "And the \"cell width\" reading of correlation plot viewer should be 20", () => readingIs(page, "cell width", el("correlation plot viewer"), 20));
      await session.step(85, "When user picks \"Show Pearson R\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Show Pearson R", "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(86, "Then the \"show pearson r\" reading of correlation plot viewer should be \"true\"", () => readingReads(page, "show pearson r", el("correlation plot viewer"), "true"));
      await session.step(87, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Open as table hands over the matrix, coefficients and all (GROK-19053)", async () => {
      await session.step(91, "When user picks \"Open as table\" from the context menu of the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => pickFromAreaContextMenu(page, "Open as table", "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(92, "Then table \"corr\" should be open", () => tableOpen(page, "corr"));
      await session.step(93, "And table \"corr\" should have 4 rows", () => tableRows(page, "corr", 4));
      await session.step(94, "When user switches to the \"corr\" table view", () => switchTableView(page, "corr"));
      await session.step(95, "Then the table should have a column \"__name\"", () => hasColumn(page, "__name"));
      await session.step(96, "And the table should have a column \"AGE\"", () => hasColumn(page, "AGE"));
      await session.step(97, "And the value of \"__name\" column in row 2 should be \"HEIGHT\"", () => valueInRow(page, "__name", 2, "HEIGHT"));
      await session.step(98, "And every value of \"AGE\" column should lie between -1 and 1", () => everyValueBetween(page, "AGE", -1, 1));
      await session.step(99, "And \"AGE\" column should have its maximum in row 3", () => maxInRow(page, "AGE", 3));
      await session.step(100, "When user closes the current view", () => closeCurrentView(page));
      await session.step(101, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(102, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Narrowed to a fraction of its width, the matrix scrolls and keeps the name column", async () => {
      await session.step(106, "Then correlation plot viewer should not have a \"horz scroll\" area", () => hasNoArea(page, el("correlation plot viewer"), "horz scroll"));
      await session.step(107, "And the \"columns shown\" reading of correlation plot viewer should be 6", () => readingIs(page, "columns shown", el("correlation plot viewer"), 6));
      await session.step(108, "When user resizes correlation plot viewer to 190 by 400", () => resizeTo(page, el("correlation plot viewer"), 190, 400));
      await session.step(109, "Then correlation plot viewer should have a \"horz scroll\" area", () => hasArea(page, el("correlation plot viewer"), "horz scroll"));
      await session.step(110, "And the \"columns shown\" reading of correlation plot viewer should be lower than before", () => readingLower(page, "columns shown", el("correlation plot viewer")));
      await session.step(111, "And correlation plot viewer should have a \"row header HEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "row header HEIGHT"));
      await session.step(112, "And correlation plot viewer should have a \"pinned band\" area", () => hasArea(page, el("correlation plot viewer"), "pinned band"));
      await session.step(113, "When user restores the size of correlation plot viewer", () => restoreSize(page, el("correlation plot viewer")));
      await session.step(114, "Then correlation plot viewer should not have a \"horz scroll\" area", () => hasNoArea(page, el("correlation plot viewer"), "horz scroll"));
      await session.step(115, "And the \"columns shown\" reading of correlation plot viewer should be 6", () => readingIs(page, "columns shown", el("correlation plot viewer"), 6));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The cell tooltip names the pair, its coefficient and the plot drawn inside it (GROK-20125)", async () => {
      await session.step(119, "When user hovers over the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => hoverArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(120, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(121, "And tooltip should contain the text \"Pearson R: -0.235\"", () => shouldContainText(page, el("tooltip"), "Pearson R: -0.235"));
      await session.step(122, "And the \"tooltip plot x column\" reading of correlation plot viewer should be \"HEIGHT\"", () => readingReads(page, "tooltip plot x column", el("correlation plot viewer"), "HEIGHT"));
      await session.step(123, "And the \"tooltip plot y column\" reading of correlation plot viewer should be \"AGE\"", () => readingReads(page, "tooltip plot y column", el("correlation plot viewer"), "AGE"));
      await session.step(124, "When user hovers over the \"cell WEIGHT x HEIGHT\" area of correlation plot viewer", () => hoverArea(page, "cell WEIGHT x HEIGHT", el("correlation plot viewer")));
      await session.step(125, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(126, "And tooltip should contain the text \"Pearson R: 0.412\"", () => shouldContainText(page, el("tooltip"), "Pearson R: 0.412"));
      await session.step(127, "And the \"tooltip plot x column\" reading of correlation plot viewer should be \"WEIGHT\"", () => readingReads(page, "tooltip plot x column", el("correlation plot viewer"), "WEIGHT"));
      await session.step(128, "And the \"tooltip plot y column\" reading of correlation plot viewer should be \"HEIGHT\"", () => readingReads(page, "tooltip plot y column", el("correlation plot viewer"), "HEIGHT"));
      await session.step(129, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(130, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The pinned name column hovers the column's statistics, not a coefficient", async () => {
      await session.step(133, "When user hovers over the \"row header HEIGHT\" area of correlation plot viewer", () => hoverArea(page, "row header HEIGHT", el("correlation plot viewer")));
      await session.step(134, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(135, "And tooltip should contain the text \"min:\"", () => shouldContainText(page, el("tooltip"), "min:"));
      await session.step(136, "And tooltip should contain the text \"nulls: 128\"", () => shouldContainText(page, el("tooltip"), "nulls: 128"));
      await session.step(137, "And tooltip should not contain the text \"Pearson R\"", () => shouldNotContainText(page, el("tooltip"), "Pearson R"));
      await session.step(138, "When user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(139, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Tooltip off still shows the cell tooltip (correlation_plot_core.dart:57)", async () => {
      await session.step(149, "When user sets \"showTooltip\" property of correlation plot viewer to \"false\"", () => setProperty(page, "showTooltip", el("correlation plot viewer"), "false"));
      await session.step(150, "And user moves the pointer away from correlation plot viewer", () => pointerAway(page, el("correlation plot viewer")));
      await session.step(151, "And user hovers over the \"cell HEIGHT x AGE\" area of correlation plot viewer", () => hoverArea(page, "cell HEIGHT x AGE", el("correlation plot viewer")));
      await session.step(152, "Then tooltip should not contain the text \"Pearson R\"", () => shouldNotContainText(page, el("tooltip"), "Pearson R"));
      await session.step(153, "When user sets \"showTooltip\" property of correlation plot viewer to \"true\"", () => setProperty(page, "showTooltip", el("correlation plot viewer"), "true"));
      await session.step(154, "Then no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true});
    run.finish();
  });
});
