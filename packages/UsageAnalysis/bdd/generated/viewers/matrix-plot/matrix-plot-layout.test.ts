/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/matrix-plot/matrix-plot-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.matrix-plot]
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
import {computedStyleContains, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, loadLayout, noErrors, propertyShouldBe, readingIs, readingReads, resizeTo, restoreSize, saveLayout, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall, descriptionAbove, descriptionBelow} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Matrix plot — the chrome around the grid and the layout that decides it", () => {
  const session = feature(test, "features/viewers/matrix-plot/matrix-plot-layout.feature", import.meta.url);
  test("Matrix plot — the chrome around the grid and the layout that decides it", {tag: ["@journey", "@viewers", "@realizes:viewers.matrix-plot", "@GROK-18736"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a matrix plot viewer", () => addViewer(page, "matrix plot"));
    await session.step(18, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
    await session.step(19, "And the \"auto layout\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "auto layout", el("matrix plot viewer"), "true"));
    await session.step(20, "And the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
    await session.step(21, "And the \"y axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "y axis shown", el("matrix plot viewer"), "true"));
    await run.scenario("Show X Axes and Show Y Axes take their strips away and bring them back", async () => {
      await session.step(24, "Then matrix plot viewer should have an \"x axis\" area", () => hasArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(25, "And matrix plot viewer should have a \"y axis\" area", () => hasArea(page, el("matrix plot viewer"), "y axis"));
      await session.step(26, "And matrix plot viewer should have an \"x axis cell 1\" area", () => hasArea(page, el("matrix plot viewer"), "x axis cell 1"));
      await session.step(27, "When user sets \"showXAxes\" property of matrix plot viewer to \"false\"", () => setProperty(page, "showXAxes", el("matrix plot viewer"), "false"));
      await session.step(28, "Then the \"x axis shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "false"));
      await session.step(29, "And matrix plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(30, "And matrix plot viewer should not have an \"x axis cell 1\" area", () => hasNoArea(page, el("matrix plot viewer"), "x axis cell 1"));
      await session.step(31, "And the \"y axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "y axis shown", el("matrix plot viewer"), "true"));
      await session.step(32, "And matrix plot viewer should have a \"y axis\" area", () => hasArea(page, el("matrix plot viewer"), "y axis"));
      await session.step(33, "When user sets \"showYAxes\" property of matrix plot viewer to \"false\"", () => setProperty(page, "showYAxes", el("matrix plot viewer"), "false"));
      await session.step(34, "Then the \"y axis shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "y axis shown", el("matrix plot viewer"), "false"));
      await session.step(35, "And matrix plot viewer should not have a \"y axis\" area", () => hasNoArea(page, el("matrix plot viewer"), "y axis"));
      await session.step(36, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["showXAxes","true"],["showYAxes","true"]]), [["showXAxes","true"],["showYAxes","true"]]);
      await session.step(39, "Then the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
      await session.step(40, "And matrix plot viewer should have an \"x axis\" area", () => hasArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(41, "And matrix plot viewer should have a \"y axis\" area", () => hasArea(page, el("matrix plot viewer"), "y axis"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Auto Layout drops the strips that no longer fit, and turning it off puts them back (GROK-19106)", async () => {
      await session.step(45, "Then the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
      await session.step(46, "And the \"x labels shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x labels shown", el("matrix plot viewer"), "true"));
      await session.step(47, "When user resizes matrix plot viewer to 260 by 240", () => resizeTo(page, el("matrix plot viewer"), 260, 240));
      await session.step(48, "Then \"showXAxes\" property of matrix plot viewer should be \"true\"", () => propertyShouldBe(page, "showXAxes", el("matrix plot viewer"), "true"));
      await session.step(49, "And \"showYAxes\" property of matrix plot viewer should be \"true\"", () => propertyShouldBe(page, "showYAxes", el("matrix plot viewer"), "true"));
      await session.step(50, "And the \"x axis shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "false"));
      await session.step(51, "And the \"y axis shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "y axis shown", el("matrix plot viewer"), "false"));
      await session.step(52, "And the \"x labels shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "x labels shown", el("matrix plot viewer"), "false"));
      await session.step(53, "And the \"y labels shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "y labels shown", el("matrix plot viewer"), "false"));
      await session.step(54, "And matrix plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(55, "And matrix plot viewer should not have an \"x label AGE\" area", () => hasNoArea(page, el("matrix plot viewer"), "x label AGE"));
      await session.step(56, "When user sets \"autoLayout\" property of matrix plot viewer to \"false\"", () => setProperty(page, "autoLayout", el("matrix plot viewer"), "false"));
      await session.step(57, "Then the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
      await session.step(58, "And the \"y axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "y axis shown", el("matrix plot viewer"), "true"));
      await session.step(59, "And the \"x labels shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x labels shown", el("matrix plot viewer"), "true"));
      await session.step(60, "And matrix plot viewer should have an \"x axis\" area", () => hasArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(61, "And matrix plot viewer should have an \"x label AGE\" area", () => hasArea(page, el("matrix plot viewer"), "x label AGE"));
      await session.step(62, "When user sets \"autoLayout\" property of matrix plot viewer to \"true\"", () => setProperty(page, "autoLayout", el("matrix plot viewer"), "true"));
      await session.step(63, "And user restores the size of matrix plot viewer", () => restoreSize(page, el("matrix plot viewer")));
      await session.step(64, "Then the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
      await session.step(65, "And the \"x labels shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x labels shown", el("matrix plot viewer"), "true"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The label strips carry one box per column of the viewport", async () => {
      await session.step(69, "Then the \"x labels shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x labels shown", el("matrix plot viewer"), "true"));
      await session.step(70, "And matrix plot viewer should have an \"x labels\" area", () => hasArea(page, el("matrix plot viewer"), "x labels"));
      await session.step(71, "And matrix plot viewer should have an \"x label AGE\" area", () => hasArea(page, el("matrix plot viewer"), "x label AGE"));
      await session.step(72, "And matrix plot viewer should have an \"x label STARTED\" area", () => hasArea(page, el("matrix plot viewer"), "x label STARTED"));
      await session.step(73, "And matrix plot viewer should have a \"y labels\" area", () => hasArea(page, el("matrix plot viewer"), "y labels"));
      await session.step(74, "And matrix plot viewer should have a \"y label WEIGHT\" area", () => hasArea(page, el("matrix plot viewer"), "y label WEIGHT"));
      await session.step(75, "And matrix plot viewer should not have an \"x label SEX\" area", () => hasNoArea(page, el("matrix plot viewer"), "x label SEX"));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Font property reaches the column labels", async () => {
      await session.step(80, "Then the \"font-size\" style of \"AGE\" text in matrix plot viewer should contain \"10px\"", () => computedStyleContains(page, "font-size", el("\"AGE\" text in matrix plot viewer"), "10px"));
      await session.step(81, "When user sets \"font\" property of matrix plot viewer to 'normal normal 24px \"Roboto\"'", () => setProperty(page, "font", el("matrix plot viewer"), "normal normal 24px \"Roboto\""));
      await session.step(82, "Then the \"font-size\" style of \"AGE\" text in matrix plot viewer should contain \"24px\"", () => computedStyleContains(page, "font-size", el("\"AGE\" text in matrix plot viewer"), "24px"));
      await session.step(83, "And the \"font-size\" style of \"WEIGHT\" text in matrix plot viewer should contain \"24px\"", () => computedStyleContains(page, "font-size", el("\"WEIGHT\" text in matrix plot viewer"), "24px"));
      await session.step(84, "When user sets \"font\" property of matrix plot viewer to 'normal normal 10px \"Roboto\"'", () => setProperty(page, "font", el("matrix plot viewer"), "normal normal 10px \"Roboto\""));
      await session.step(85, "Then the \"font-size\" style of \"AGE\" text in matrix plot viewer should contain \"10px\"", () => computedStyleContains(page, "font-size", el("\"AGE\" text in matrix plot viewer"), "10px"));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title shows and the description sits above the grid, then below it", async () => {
      await session.step(89, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["showTitle","true"],["title","Pairwise"],["description","Every numerical pair"]]), [["showTitle","true"],["title","Pairwise"],["description","Every numerical pair"]]);
      await session.step(93, "Then title of matrix plot viewer should have text \"Pairwise\"", () => shouldHaveText(page, el("title of matrix plot viewer"), "Pairwise"));
      await session.step(94, "And description of matrix plot viewer should be visible", () => shouldBe(page, el("description of matrix plot viewer"), "visible"));
      await session.step(95, "And description of matrix plot viewer should have text \"Every numerical pair\"", () => shouldHaveText(page, el("description of matrix plot viewer"), "Every numerical pair"));
      await session.step(96, "And the description of matrix plot viewer should be above its content", () => descriptionAbove(page, el("matrix plot viewer")));
      await session.step(97, "When user sets \"descriptionPosition\" property of matrix plot viewer to \"Bottom\"", () => setProperty(page, "descriptionPosition", el("matrix plot viewer"), "Bottom"));
      await session.step(98, "Then the description of matrix plot viewer should be below its content", () => descriptionBelow(page, el("matrix plot viewer")));
      await session.step(99, "When user sets \"descriptionVisibilityMode\" property of matrix plot viewer to \"Never\"", () => setProperty(page, "descriptionVisibilityMode", el("matrix plot viewer"), "Never"));
      await session.step(100, "Then description of matrix plot viewer should be hidden", () => shouldBe(page, el("description of matrix plot viewer"), "hidden"));
      await session.step(101, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["descriptionVisibilityMode","Auto"],["descriptionPosition","Top"],["description",""],["title",""],["showTitle","false"]]), [["descriptionVisibilityMode","Auto"],["descriptionPosition","Top"],["description",""],["title",""],["showTitle","false"]]);
      await session.step(107, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the configured grid back", async () => {
      await session.step(111, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["cellPlotType","Scatter plot"],["xColumnNames","AGE, HEIGHT"],["yColumnNames","AGE, HEIGHT, WEIGHT"]]), [["cellPlotType","Scatter plot"],["xColumnNames","AGE, HEIGHT"],["yColumnNames","AGE, HEIGHT, WEIGHT"]]);
      await session.step(115, "Then the \"cells\" reading of matrix plot viewer should be 6", () => readingIs(page, "cells", el("matrix plot viewer"), 6));
      await session.step(116, "And the cells of matrix plot viewer should be 2 wide and 3 tall", () => cellsWideTall(page, el("matrix plot viewer"), 2, 3));
      await session.step(117, "And user saves the layout of the current table view", () => saveLayout(page));
      await session.step(118, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["cellPlotType","Density plot"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]), [["cellPlotType","Density plot"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]);
      await session.step(122, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(123, "When user loads the saved layout", () => loadLayout(page));
      await session.step(124, "Then matrix plot viewer should be visible", () => shouldBe(page, el("matrix plot viewer"), "visible"));
      await session.step(125, "And the \"cells\" reading of matrix plot viewer should be 6", () => readingIs(page, "cells", el("matrix plot viewer"), 6));
      await session.step(126, "And the cells of matrix plot viewer should be 2 wide and 3 tall", () => cellsWideTall(page, el("matrix plot viewer"), 2, 3));
      await session.step(127, "And the \"cell viewer type\" reading of matrix plot viewer should be \"Scatter plot\"", () => readingReads(page, "cell viewer type", el("matrix plot viewer"), "Scatter plot"));
      await session.step(128, "And the \"cell viewer type of HEIGHT x AGE\" reading of matrix plot viewer should be \"Scatter plot\"", () => readingReads(page, "cell viewer type of HEIGHT x AGE", el("matrix plot viewer"), "Scatter plot"));
      await session.step(129, "And the \"cell viewer type of AGE x AGE\" reading of matrix plot viewer should be \"Histogram\"", () => readingReads(page, "cell viewer type of AGE x AGE", el("matrix plot viewer"), "Histogram"));
      await session.step(130, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["cellPlotType","Density plot"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]), [["cellPlotType","Density plot"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]);
      await session.step(134, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project round-trip brings it back too (GROK-10925)", async () => {
      await session.step(138, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["cellPlotType","Scatter plot"],["xColumnNames","AGE, HEIGHT"],["yColumnNames","AGE, HEIGHT, WEIGHT"]]), [["cellPlotType","Scatter plot"],["xColumnNames","AGE, HEIGHT"],["yColumnNames","AGE, HEIGHT, WEIGHT"]]);
      await session.step(142, "Then the \"cells\" reading of matrix plot viewer should be 6", () => readingIs(page, "cells", el("matrix plot viewer"), 6));
      await session.step(143, "When user saves the current view as project \"bdd matrix plot grid\"", () => saveAsProject(page, "bdd matrix plot grid"));
      await session.step(144, "And user closes all views", () => closeAllViews(page));
      await session.step(145, "And user opens the \"bdd matrix plot grid\" project", () => openProject(page, "bdd matrix plot grid"));
      await session.step(146, "Then matrix plot viewer should be visible", () => shouldBe(page, el("matrix plot viewer"), "visible"));
      await session.step(147, "And the \"cells\" reading of matrix plot viewer should be 6", () => readingIs(page, "cells", el("matrix plot viewer"), 6));
      await session.step(148, "And the cells of matrix plot viewer should be 2 wide and 3 tall", () => cellsWideTall(page, el("matrix plot viewer"), 2, 3));
      await session.step(149, "And the \"cell viewer type\" reading of matrix plot viewer should be \"Scatter plot\"", () => readingReads(page, "cell viewer type", el("matrix plot viewer"), "Scatter plot"));
      await session.step(150, "And the \"rows shown\" reading of matrix plot viewer should be 1000", () => readingIs(page, "rows shown", el("matrix plot viewer"), 1000));
      await session.step(151, "And the \"cell rows shown of HEIGHT x AGE\" reading of matrix plot viewer should be 872", () => readingIs(page, "cell rows shown of HEIGHT x AGE", el("matrix plot viewer"), 872));
      await session.step(152, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["cellPlotType","Density plot"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]), [["cellPlotType","Density plot"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]);
      await session.step(156, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(157, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
