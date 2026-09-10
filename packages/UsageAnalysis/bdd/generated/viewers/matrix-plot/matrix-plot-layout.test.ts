/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/matrix-plot/matrix-plot-layout.feature
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
import {shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, loadLayout, noErrors, propertyShouldBe, readingIs, readingReads, resizeTo, restoreSize, saveLayout, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall, descriptionAbove, descriptionBelow} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Matrix plot — the chrome around the grid and the layout that decides it", () => {
  const session = feature(test, "features/viewers/matrix-plot/matrix-plot-layout.feature", import.meta.url);
  test("Matrix plot — the chrome around the grid and the layout that decides it", {tag: ["@journey", "@viewers", "@realizes:viewers.matrix-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a matrix plot viewer", () => addViewer(page, "matrix plot"));
    await session.step(20, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
    await session.step(21, "And the \"auto layout\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "auto layout", el("matrix plot viewer"), "true"));
    await session.step(22, "And the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
    await session.step(23, "And the \"y axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "y axis shown", el("matrix plot viewer"), "true"));
    await run.scenario("Show X Axes and Show Y Axes take their strips away and bring them back", async () => {
      await session.step(26, "Then matrix plot viewer should have an \"x axis\" area", () => hasArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(27, "And matrix plot viewer should have a \"y axis\" area", () => hasArea(page, el("matrix plot viewer"), "y axis"));
      await session.step(28, "And matrix plot viewer should have an \"x axis cell 1\" area", () => hasArea(page, el("matrix plot viewer"), "x axis cell 1"));
      await session.step(29, "When user sets \"showXAxes\" property of matrix plot viewer to \"false\"", () => setProperty(page, "showXAxes", el("matrix plot viewer"), "false"));
      await session.step(30, "Then the \"x axis shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "false"));
      await session.step(31, "And matrix plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(32, "And matrix plot viewer should not have an \"x axis cell 1\" area", () => hasNoArea(page, el("matrix plot viewer"), "x axis cell 1"));
      await session.step(33, "And the \"y axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "y axis shown", el("matrix plot viewer"), "true"));
      await session.step(34, "And matrix plot viewer should have a \"y axis\" area", () => hasArea(page, el("matrix plot viewer"), "y axis"));
      await session.step(35, "When user sets \"showYAxes\" property of matrix plot viewer to \"false\"", () => setProperty(page, "showYAxes", el("matrix plot viewer"), "false"));
      await session.step(36, "Then the \"y axis shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "y axis shown", el("matrix plot viewer"), "false"));
      await session.step(37, "And matrix plot viewer should not have a \"y axis\" area", () => hasNoArea(page, el("matrix plot viewer"), "y axis"));
      await session.step(38, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["showXAxes","true"],["showYAxes","true"]]));
      await session.step(41, "Then the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
      await session.step(42, "And matrix plot viewer should have an \"x axis\" area", () => hasArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(43, "And matrix plot viewer should have a \"y axis\" area", () => hasArea(page, el("matrix plot viewer"), "y axis"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Auto Layout drops the strips that no longer fit, and turning it off puts them back (GROK-19106)", async () => {
      await session.step(47, "Then the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
      await session.step(48, "And the \"x labels shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x labels shown", el("matrix plot viewer"), "true"));
      await session.step(49, "When user resizes matrix plot viewer to 260 by 240", () => resizeTo(page, el("matrix plot viewer"), 260, 240));
      await session.step(50, "Then \"showXAxes\" property of matrix plot viewer should be \"true\"", () => propertyShouldBe(page, "showXAxes", el("matrix plot viewer"), "true"));
      await session.step(51, "And \"showYAxes\" property of matrix plot viewer should be \"true\"", () => propertyShouldBe(page, "showYAxes", el("matrix plot viewer"), "true"));
      await session.step(52, "And the \"x axis shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "false"));
      await session.step(53, "And the \"y axis shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "y axis shown", el("matrix plot viewer"), "false"));
      await session.step(54, "And the \"x labels shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "x labels shown", el("matrix plot viewer"), "false"));
      await session.step(55, "And the \"y labels shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "y labels shown", el("matrix plot viewer"), "false"));
      await session.step(56, "And matrix plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(57, "And matrix plot viewer should not have an \"x label AGE\" area", () => hasNoArea(page, el("matrix plot viewer"), "x label AGE"));
      await session.step(58, "When user sets \"autoLayout\" property of matrix plot viewer to \"false\"", () => setProperty(page, "autoLayout", el("matrix plot viewer"), "false"));
      await session.step(59, "Then the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
      await session.step(60, "And the \"y axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "y axis shown", el("matrix plot viewer"), "true"));
      await session.step(61, "And the \"x labels shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x labels shown", el("matrix plot viewer"), "true"));
      await session.step(62, "And matrix plot viewer should have an \"x axis\" area", () => hasArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(63, "And matrix plot viewer should have an \"x label AGE\" area", () => hasArea(page, el("matrix plot viewer"), "x label AGE"));
      await session.step(64, "When user sets \"autoLayout\" property of matrix plot viewer to \"true\"", () => setProperty(page, "autoLayout", el("matrix plot viewer"), "true"));
      await session.step(65, "And user restores the size of matrix plot viewer", () => restoreSize(page, el("matrix plot viewer")));
      await session.step(66, "Then the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
      await session.step(67, "And the \"x labels shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x labels shown", el("matrix plot viewer"), "true"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The label strips carry one box per column of the viewport", async () => {
      await session.step(71, "Then the \"x labels shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x labels shown", el("matrix plot viewer"), "true"));
      await session.step(72, "And matrix plot viewer should have an \"x labels\" area", () => hasArea(page, el("matrix plot viewer"), "x labels"));
      await session.step(73, "And matrix plot viewer should have an \"x label AGE\" area", () => hasArea(page, el("matrix plot viewer"), "x label AGE"));
      await session.step(74, "And matrix plot viewer should have an \"x label STARTED\" area", () => hasArea(page, el("matrix plot viewer"), "x label STARTED"));
      await session.step(75, "And matrix plot viewer should have a \"y labels\" area", () => hasArea(page, el("matrix plot viewer"), "y labels"));
      await session.step(76, "And matrix plot viewer should have a \"y label WEIGHT\" area", () => hasArea(page, el("matrix plot viewer"), "y label WEIGHT"));
      await session.step(77, "And matrix plot viewer should not have an \"x label SEX\" area", () => hasNoArea(page, el("matrix plot viewer"), "x label SEX"));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title shows and the description sits above the grid, then below it", async () => {
      await session.step(81, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["showTitle","true"],["title","Pairwise"],["description","Every numerical pair"]]));
      await session.step(85, "Then title of matrix plot viewer should have text \"Pairwise\"", () => shouldHaveText(page, el("title of matrix plot viewer"), "Pairwise"));
      await session.step(86, "And description of matrix plot viewer should be visible", () => shouldBe(page, el("description of matrix plot viewer"), "visible"));
      await session.step(87, "And description of matrix plot viewer should have text \"Every numerical pair\"", () => shouldHaveText(page, el("description of matrix plot viewer"), "Every numerical pair"));
      await session.step(88, "And the description of matrix plot viewer should be above its content", () => descriptionAbove(page, el("matrix plot viewer")));
      await session.step(89, "When user sets \"descriptionPosition\" property of matrix plot viewer to \"Bottom\"", () => setProperty(page, "descriptionPosition", el("matrix plot viewer"), "Bottom"));
      await session.step(90, "Then the description of matrix plot viewer should be below its content", () => descriptionBelow(page, el("matrix plot viewer")));
      await session.step(91, "When user sets \"descriptionVisibilityMode\" property of matrix plot viewer to \"Never\"", () => setProperty(page, "descriptionVisibilityMode", el("matrix plot viewer"), "Never"));
      await session.step(92, "Then description of matrix plot viewer should be hidden", () => shouldBe(page, el("description of matrix plot viewer"), "hidden"));
      await session.step(93, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["descriptionVisibilityMode","Auto"],["descriptionPosition","Top"],["description",""],["title",""],["showTitle","false"]]));
      await session.step(99, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the configured grid back", async () => {
      await session.step(103, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["cellPlotType","Scatter plot"],["xColumnNames","AGE, HEIGHT"],["yColumnNames","AGE, HEIGHT, WEIGHT"]]));
      await session.step(107, "Then the \"cells\" reading of matrix plot viewer should be 6", () => readingIs(page, "cells", el("matrix plot viewer"), 6));
      await session.step(108, "And the cells of matrix plot viewer should be 2 wide and 3 tall", () => cellsWideTall(page, el("matrix plot viewer"), 2, 3));
      await session.step(109, "And user saves the layout of the current table view", () => saveLayout(page));
      await session.step(110, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["cellPlotType","Density plot"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]));
      await session.step(114, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(115, "When user loads the saved layout", () => loadLayout(page));
      await session.step(116, "Then matrix plot viewer should be visible", () => shouldBe(page, el("matrix plot viewer"), "visible"));
      await session.step(117, "And the \"cells\" reading of matrix plot viewer should be 6", () => readingIs(page, "cells", el("matrix plot viewer"), 6));
      await session.step(118, "And the cells of matrix plot viewer should be 2 wide and 3 tall", () => cellsWideTall(page, el("matrix plot viewer"), 2, 3));
      await session.step(119, "And the \"cell viewer type\" reading of matrix plot viewer should be \"Scatter plot\"", () => readingReads(page, "cell viewer type", el("matrix plot viewer"), "Scatter plot"));
      await session.step(120, "And the \"cell viewer type of HEIGHT x AGE\" reading of matrix plot viewer should be \"Scatter plot\"", () => readingReads(page, "cell viewer type of HEIGHT x AGE", el("matrix plot viewer"), "Scatter plot"));
      await session.step(121, "And the \"cell viewer type of AGE x AGE\" reading of matrix plot viewer should be \"Histogram\"", () => readingReads(page, "cell viewer type of AGE x AGE", el("matrix plot viewer"), "Histogram"));
      await session.step(122, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["cellPlotType","Density plot"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]));
      await session.step(126, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project round-trip brings it back too (GROK-10925)", async () => {
      await session.step(130, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["cellPlotType","Scatter plot"],["xColumnNames","AGE, HEIGHT"],["yColumnNames","AGE, HEIGHT, WEIGHT"]]));
      await session.step(134, "Then the \"cells\" reading of matrix plot viewer should be 6", () => readingIs(page, "cells", el("matrix plot viewer"), 6));
      await session.step(135, "When user saves the current view as project \"bdd matrix plot grid\"", () => saveAsProject(page, "bdd matrix plot grid"));
      await session.step(136, "And user closes all views", () => closeAllViews(page));
      await session.step(137, "And user opens the \"bdd matrix plot grid\" project", () => openProject(page, "bdd matrix plot grid"));
      await session.step(138, "Then matrix plot viewer should be visible", () => shouldBe(page, el("matrix plot viewer"), "visible"));
      await session.step(139, "And the \"cells\" reading of matrix plot viewer should be 6", () => readingIs(page, "cells", el("matrix plot viewer"), 6));
      await session.step(140, "And the cells of matrix plot viewer should be 2 wide and 3 tall", () => cellsWideTall(page, el("matrix plot viewer"), 2, 3));
      await session.step(141, "And the \"cell viewer type\" reading of matrix plot viewer should be \"Scatter plot\"", () => readingReads(page, "cell viewer type", el("matrix plot viewer"), "Scatter plot"));
      await session.step(142, "And the \"rows shown\" reading of matrix plot viewer should be 1000", () => readingIs(page, "rows shown", el("matrix plot viewer"), 1000));
      await session.step(143, "And the \"cell rows shown of HEIGHT x AGE\" reading of matrix plot viewer should be 872", () => readingIs(page, "cell rows shown of HEIGHT x AGE", el("matrix plot viewer"), 872));
      await session.step(144, "When user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["cellPlotType","Density plot"],["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]));
      await session.step(148, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
      await session.step(149, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
