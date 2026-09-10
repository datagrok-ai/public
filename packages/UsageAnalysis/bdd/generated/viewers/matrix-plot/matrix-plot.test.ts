/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/matrix-plot/matrix-plot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.matrix-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, filterPasses, selectFirstRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, hasArea, hasNoArea, hoverArea, noErrors, oneTooltip, pointerAway, readingAsRemembered, readingDiffers, readingIs, readingNotAsRemembered, readingReads, readingSame, rememberReading, setProperties, setProperty, viewerCount, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Matrix plot — the cells, the plots inside them and the rows they keep", () => {
  const session = feature(test, "features/viewers/matrix-plot/matrix-plot.feature", import.meta.url);
  test("Matrix plot — the cells, the plots inside them and the rows they keep", {tag: ["@journey", "@viewers", "@realizes:viewers.matrix-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a matrix plot viewer", () => addViewer(page, "matrix plot"));
    await session.step(20, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(21, "And the \"rows shown\" reading of matrix plot viewer should be 1000", () => readingIs(page, "rows shown", el("matrix plot viewer"), 1000));
    await session.step(22, "And the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
    await session.step(23, "And the cells of matrix plot viewer should be 4 wide and 4 tall", () => cellsWideTall(page, el("matrix plot viewer"), 4, 4));
    await session.step(24, "And the \"cell viewer type\" reading of matrix plot viewer should be \"Density plot\"", () => readingReads(page, "cell viewer type", el("matrix plot viewer"), "Density plot"));
    await session.step(25, "And the \"cells drawn\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 16));
    await run.scenario("One cell per pair of numerical columns, each with a picture of its own", async () => {
      await session.step(28, "Then the \"x columns\" reading of matrix plot viewer should be 4", () => readingIs(page, "x columns", el("matrix plot viewer"), 4));
      await session.step(29, "And the \"y columns\" reading of matrix plot viewer should be 4", () => readingIs(page, "y columns", el("matrix plot viewer"), 4));
      await session.step(30, "And matrix plot viewer should have a \"cell HEIGHT x AGE\" area", () => hasArea(page, el("matrix plot viewer"), "cell HEIGHT x AGE"));
      await session.step(31, "And matrix plot viewer should have a \"cell body HEIGHT x AGE\" area", () => hasArea(page, el("matrix plot viewer"), "cell body HEIGHT x AGE"));
      await session.step(32, "And matrix plot viewer should have a \"cell 1,1\" area", () => hasArea(page, el("matrix plot viewer"), "cell 1,1"));
      await session.step(33, "And matrix plot viewer should not have a \"cell SEX x AGE\" area", () => hasNoArea(page, el("matrix plot viewer"), "cell SEX x AGE"));
      await session.step(34, "And the \"blank cells\" reading of matrix plot viewer should be 0", () => readingIs(page, "blank cells", el("matrix plot viewer"), 0));
      await session.step(35, "And the \"distinct cell signatures\" reading of matrix plot viewer should be 16", () => readingIs(page, "distinct cell signatures", el("matrix plot viewer"), 16));
      await session.step(36, "And the \"error\" reading of matrix plot viewer should be \"\"", () => readingReads(page, "error", el("matrix plot viewer"), ""));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The diagonal hosts a histogram whatever Cell Plot Type says", async () => {
      await session.step(40, "Then the \"cell viewer type of AGE x AGE\" reading of matrix plot viewer should be \"Histogram\"", () => readingReads(page, "cell viewer type of AGE x AGE", el("matrix plot viewer"), "Histogram"));
      await session.step(41, "And the \"cell viewer type of WEIGHT x WEIGHT\" reading of matrix plot viewer should be \"Histogram\"", () => readingReads(page, "cell viewer type of WEIGHT x WEIGHT", el("matrix plot viewer"), "Histogram"));
      await session.step(42, "And the \"cell viewer type of HEIGHT x AGE\" reading of matrix plot viewer should be \"Density plot\"", () => readingReads(page, "cell viewer type of HEIGHT x AGE", el("matrix plot viewer"), "Density plot"));
      await session.step(43, "When user sets \"cellPlotType\" property of matrix plot viewer to \"Scatter plot\"", () => setProperty(page, "cellPlotType", el("matrix plot viewer"), "Scatter plot"));
      await session.step(44, "Then the \"cell viewer type of HEIGHT x AGE\" reading of matrix plot viewer should be \"Scatter plot\"", () => readingReads(page, "cell viewer type of HEIGHT x AGE", el("matrix plot viewer"), "Scatter plot"));
      await session.step(45, "And the \"cell viewer type of AGE x AGE\" reading of matrix plot viewer should be \"Histogram\"", () => readingReads(page, "cell viewer type of AGE x AGE", el("matrix plot viewer"), "Histogram"));
      await session.step(46, "And the \"cell viewer type of WEIGHT x WEIGHT\" reading of matrix plot viewer should be \"Histogram\"", () => readingReads(page, "cell viewer type of WEIGHT x WEIGHT", el("matrix plot viewer"), "Histogram"));
      await session.step(47, "When user sets \"cellPlotType\" property of matrix plot viewer to \"Density plot\"", () => setProperty(page, "cellPlotType", el("matrix plot viewer"), "Density plot"));
      await session.step(48, "Then the \"cell viewer type of HEIGHT x AGE\" reading of matrix plot viewer should be \"Density plot\"", () => readingReads(page, "cell viewer type of HEIGHT x AGE", el("matrix plot viewer"), "Density plot"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cell Plot Type redraws every off-diagonal cell and the round trip returns to the picture it started from", async () => {
      await session.step(52, "Then the \"distinct cell signatures\" reading of matrix plot viewer should be 16", () => readingIs(page, "distinct cell signatures", el("matrix plot viewer"), 16));
      await session.step(53, "When user remembers the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer", () => rememberReading(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(54, "And user sets \"cellPlotType\" property of matrix plot viewer to \"Scatter plot\"", () => setProperty(page, "cellPlotType", el("matrix plot viewer"), "Scatter plot"));
      await session.step(55, "Then the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(56, "And the \"cell signature WEIGHT x HEIGHT\" reading of matrix plot viewer should differ from before", () => readingDiffers(page, "cell signature WEIGHT x HEIGHT", el("matrix plot viewer")));
      await session.step(57, "And the \"cell signature AGE x AGE\" reading of matrix plot viewer should be the same as before", () => readingSame(page, "cell signature AGE x AGE", el("matrix plot viewer")));
      await session.step(58, "And the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(59, "And the \"blank cells\" reading of matrix plot viewer should be 0", () => readingIs(page, "blank cells", el("matrix plot viewer"), 0));
      await session.step(60, "When user sets \"cellPlotType\" property of matrix plot viewer to \"Density plot\"", () => setProperty(page, "cellPlotType", el("matrix plot viewer"), "Density plot"));
      await session.step(61, "Then the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Each cell keeps its own rows, so a blank in either column drops the row from that cell only", async () => {
      await session.step(65, "Then the \"rows shown\" reading of matrix plot viewer should be 1000", () => readingIs(page, "rows shown", el("matrix plot viewer"), 1000));
      await session.step(66, "And the \"cell rows shown of HEIGHT x AGE\" reading of matrix plot viewer should be 872", () => readingIs(page, "cell rows shown of HEIGHT x AGE", el("matrix plot viewer"), 872));
      await session.step(67, "And the \"cell rows shown of WEIGHT x AGE\" reading of matrix plot viewer should be 1000", () => readingIs(page, "cell rows shown of WEIGHT x AGE", el("matrix plot viewer"), 1000));
      await session.step(68, "And the \"cell rows shown of AGE x AGE\" reading of matrix plot viewer should be 1000", () => readingIs(page, "cell rows shown of AGE x AGE", el("matrix plot viewer"), 1000));
      await session.step(69, "And the \"cell rows shown of STARTED x HEIGHT\" reading of matrix plot viewer should be 872", () => readingIs(page, "cell rows shown of STARTED x HEIGHT", el("matrix plot viewer"), 872));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Narrowing X re-tiles the grid and takes the labels with it", async () => {
      await session.step(73, "When user sets \"xColumnNames\" property of matrix plot viewer to \"AGE, HEIGHT\"", () => setProperty(page, "xColumnNames", el("matrix plot viewer"), "AGE, HEIGHT"));
      await session.step(74, "Then the \"cells\" reading of matrix plot viewer should be 8", () => readingIs(page, "cells", el("matrix plot viewer"), 8));
      await session.step(75, "And the cells of matrix plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("matrix plot viewer"), 2, 4));
      await session.step(76, "And the \"cells drawn\" reading of matrix plot viewer should be 8", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 8));
      await session.step(77, "And matrix plot viewer should not have a \"cell WEIGHT x AGE\" area", () => hasNoArea(page, el("matrix plot viewer"), "cell WEIGHT x AGE"));
      await session.step(78, "And matrix plot viewer should have a \"x label HEIGHT\" area", () => hasArea(page, el("matrix plot viewer"), "x label HEIGHT"));
      await session.step(79, "And matrix plot viewer should not have a \"x label WEIGHT\" area", () => hasNoArea(page, el("matrix plot viewer"), "x label WEIGHT"));
      await session.step(80, "And matrix plot viewer should have a \"y label WEIGHT\" area", () => hasArea(page, el("matrix plot viewer"), "y label WEIGHT"));
      await session.step(81, "When user sets \"xColumnNames\" property of matrix plot viewer to \"AGE, HEIGHT, WEIGHT, STARTED\"", () => setProperty(page, "xColumnNames", el("matrix plot viewer"), "AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(82, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(83, "And matrix plot viewer should have a \"cell WEIGHT x AGE\" area", () => hasArea(page, el("matrix plot viewer"), "cell WEIGHT x AGE"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cycling the column sets ends where it started, with no error (GROK-16473)", async () => {
      await session.step(87, "When user sets \"xColumnNames\" property of matrix plot viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "xColumnNames", el("matrix plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(88, "Then the \"cells\" reading of matrix plot viewer should be 12", () => readingIs(page, "cells", el("matrix plot viewer"), 12));
      await session.step(89, "When user sets \"yColumnNames\" property of matrix plot viewer to \"AGE, HEIGHT\"", () => setProperty(page, "yColumnNames", el("matrix plot viewer"), "AGE, HEIGHT"));
      await session.step(90, "Then the \"cells\" reading of matrix plot viewer should be 6", () => readingIs(page, "cells", el("matrix plot viewer"), 6));
      await session.step(91, "And the cells of matrix plot viewer should be 3 wide and 2 tall", () => cellsWideTall(page, el("matrix plot viewer"), 3, 2));
      await session.step(92, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]));
      await session.step(95, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(96, "And the \"cells drawn\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 16));
      await session.step(97, "And the \"blank cells\" reading of matrix plot viewer should be 0", () => readingIs(page, "blank cells", el("matrix plot viewer"), 0));
      await session.step(98, "And the \"error\" reading of matrix plot viewer should be \"\"", () => readingReads(page, "error", el("matrix plot viewer"), ""));
      await session.step(99, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer's own filter narrows every cell and leaves the table alone", async () => {
      await session.step(102, "When user sets \"filter\" property of matrix plot viewer to \"${AGE} > 30\"", () => setProperty(page, "filter", el("matrix plot viewer"), "${AGE} > 30"));
      await session.step(103, "Then the \"rows shown\" reading of matrix plot viewer should be 844", () => readingIs(page, "rows shown", el("matrix plot viewer"), 844));
      await session.step(104, "And 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(105, "And the \"cell rows shown of HEIGHT x AGE\" reading of matrix plot viewer should be 745", () => readingIs(page, "cell rows shown of HEIGHT x AGE", el("matrix plot viewer"), 745));
      await session.step(106, "And the \"cell rows shown of AGE x AGE\" reading of matrix plot viewer should be 844", () => readingIs(page, "cell rows shown of AGE x AGE", el("matrix plot viewer"), 844));
      await session.step(107, "And the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer should differ from before", () => readingDiffers(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(108, "When user sets \"filter\" property of matrix plot viewer to \"\"", () => setProperty(page, "filter", el("matrix plot viewer"), ""));
      await session.step(109, "Then the \"rows shown\" reading of matrix plot viewer should be 1000", () => readingIs(page, "rows shown", el("matrix plot viewer"), 1000));
      await session.step(110, "And the \"cell rows shown of HEIGHT x AGE\" reading of matrix plot viewer should be 872", () => readingIs(page, "cell rows shown of HEIGHT x AGE", el("matrix plot viewer"), 872));
      await session.step(111, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source Selected redraws every cell over the selection alone", async () => {
      await session.step(114, "Then the \"distinct cell signatures\" reading of matrix plot viewer should be 16", () => readingIs(page, "distinct cell signatures", el("matrix plot viewer"), 16));
      await session.step(115, "When user remembers the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer", () => rememberReading(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(116, "And user selects the first 50 rows", () => selectFirstRows(page, 50));
      await session.step(117, "Then 50 rows should be selected", () => selectedRowCount(page, 50));
      await session.step(118, "When user sets \"rowSource\" property of matrix plot viewer to \"Selected\"", () => setProperty(page, "rowSource", el("matrix plot viewer"), "Selected"));
      await session.step(119, "Then the \"rows shown\" reading of matrix plot viewer should be 50", () => readingIs(page, "rows shown", el("matrix plot viewer"), 50));
      await session.step(120, "And the \"cell rows shown of HEIGHT x AGE\" reading of matrix plot viewer should be 50", () => readingIs(page, "cell rows shown of HEIGHT x AGE", el("matrix plot viewer"), 50));
      await session.step(121, "And the \"cell rows shown of AGE x AGE\" reading of matrix plot viewer should be 50", () => readingIs(page, "cell rows shown of AGE x AGE", el("matrix plot viewer"), 50));
      await session.step(122, "And the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(123, "When user sets \"rowSource\" property of matrix plot viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("matrix plot viewer"), "Filtered"));
      await session.step(124, "And user clears the row selection", () => clearSelection(page));
      await session.step(125, "Then the \"rows shown\" reading of matrix plot viewer should be 1000", () => readingIs(page, "rows shown", el("matrix plot viewer"), 1000));
      await session.step(126, "And the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The cell tooltip names the pair the cell stands for", async () => {
      await session.step(130, "When user hovers over the \"cell HEIGHT x AGE\" area of matrix plot viewer", () => hoverArea(page, "cell HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(131, "Then exactly one tooltip should be shown", () => oneTooltip(page));
      await session.step(132, "And tooltip should contain the text \"X: HEIGHT\"", () => shouldContainText(page, el("tooltip"), "X: HEIGHT"));
      await session.step(133, "And tooltip should contain the text \"Y: AGE\"", () => shouldContainText(page, el("tooltip"), "Y: AGE"));
      await session.step(134, "When user moves the pointer away from matrix plot viewer", () => pointerAway(page, el("matrix plot viewer")));
      await session.step(135, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The expand icon shows on hover and opens the cell as the viewer it hosts", async () => {
      await session.step(138, "When user moves the pointer away from matrix plot viewer", () => pointerAway(page, el("matrix plot viewer")));
      await session.step(139, "Then matrix plot viewer should not have an \"expand icon of HEIGHT x AGE\" area", () => hasNoArea(page, el("matrix plot viewer"), "expand icon of HEIGHT x AGE"));
      await session.step(140, "And the open tableview should have 0 density plot viewers", () => viewerCount(page, 0, "density plot"));
      await session.step(141, "When user hovers over the \"cell HEIGHT x AGE\" area of matrix plot viewer", () => hoverArea(page, "cell HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(142, "Then matrix plot viewer should have an \"expand icon of HEIGHT x AGE\" area", () => hasArea(page, el("matrix plot viewer"), "expand icon of HEIGHT x AGE"));
      await session.step(143, "When user clicks on the \"expand icon of HEIGHT x AGE\" area of matrix plot viewer", () => clickArea(page, "expand icon of HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(144, "Then the open tableview should have 1 density plot viewer", () => viewerCount(page, 1, "density plot"));
      await session.step(145, "When user clicks on close icon of density plot viewer", () => clickOn(page, el("close icon of density plot viewer")));
      await session.step(146, "Then the open tableview should have 0 density plot viewers", () => viewerCount(page, 0, "density plot"));
      await session.step(147, "When user hovers over the \"cell AGE x AGE\" area of matrix plot viewer", () => hoverArea(page, "cell AGE x AGE", el("matrix plot viewer")));
      await session.step(148, "And user clicks on the \"expand icon of AGE x AGE\" area of matrix plot viewer", () => clickArea(page, "expand icon of AGE x AGE", el("matrix plot viewer")));
      await session.step(149, "Then the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
      await session.step(150, "When user clicks on close icon of histogram viewer", () => clickOn(page, el("close icon of histogram viewer")));
      await session.step(151, "Then the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
      await session.step(152, "And user moves the pointer away from matrix plot viewer", () => pointerAway(page, el("matrix plot viewer")));
      await session.step(153, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The wheel zooms the cell under the pointer and leaves its neighbour alone", async () => {
      await session.step(161, "Then the \"distinct cell signatures\" reading of matrix plot viewer should be 16", () => readingIs(page, "distinct cell signatures", el("matrix plot viewer"), 16));
      await session.step(162, "When user remembers the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer", () => rememberReading(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(163, "And user scrolls the mouse wheel up over the \"cell body HEIGHT x AGE\" area of matrix plot viewer", () => wheelOverArea(page, "up", "cell body HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(164, "Then the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(165, "And the \"cell signature WEIGHT x AGE\" reading of matrix plot viewer should be the same as before", () => readingSame(page, "cell signature WEIGHT x AGE", el("matrix plot viewer")));
      await session.step(166, "When user scrolls the mouse wheel down over the \"cell body HEIGHT x AGE\" area of matrix plot viewer", () => wheelOverArea(page, "down", "cell body HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(167, "Then the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(168, "When user moves the pointer away from matrix plot viewer", () => pointerAway(page, el("matrix plot viewer")));
      await session.step(169, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
