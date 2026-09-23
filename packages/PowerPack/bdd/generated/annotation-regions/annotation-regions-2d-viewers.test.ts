/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/annotation-regions/annotation-regions-2d-viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart, viewers.density-plot]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/formula-lines.js';
import '../../bindings/home.js';
import '../../bindings/io.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, dragAcrossArea, hasArea, hasNoArea, noErrors, pickFromContextMenu, propertyShouldContain, readingAtLeast, readingIs, readingReads, resizeTo, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Annotation regions on the other two-dimensional viewers", () => {
  const session = feature(test, "features/annotation-regions/annotation-regions-2d-viewers.feature", import.meta.url);
  test("A region drawn on a line chart is keyed to the aggregated Y column", {tag: ["@viewers", "@realizes:viewers.line-chart", "@realizes:viewers.density-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "Given user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"]]);
    await session.step(27, "And user resizes line chart viewer to 800 by 500", () => resizeTo(page, el("line chart viewer"), 800, 500));
    await session.step(28, "Then the \"aggregated\" reading of line chart viewer should be \"true\"", () => readingReads(page, "aggregated", el("line chart viewer"), "true"));
    await session.step(29, "When user picks \"Tools > Draw Annotation Region\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Tools > Draw Annotation Region", el("line chart viewer")));
    await session.step(30, "And user drags across the \"plot\" area of line chart viewer", () => dragAcrossArea(page, "plot", el("line chart viewer")));
    await session.step(31, "Then the \"viewer regions\" reading of line chart viewer should be 1", () => readingIs(page, "viewer regions", el("line chart viewer"), 1));
    await session.step(32, "And \"annotationRegions\" property of line chart viewer should contain \"avg(WEIGHT)\"", () => propertyShouldContain(page, "annotationRegions", el("line chart viewer"), "avg(WEIGHT)"));
    await session.step(33, "And line chart viewer should have a \"region 1\" area", () => hasArea(page, el("line chart viewer"), "region 1"));
    await session.step(34, "When user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(35, "Then no errors should have been logged", () => noErrors(page));
  });
  test("A line chart draws a titled area region on the chart of its Y column", {tag: ["@viewers", "@realizes:viewers.line-chart", "@realizes:viewers.density-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(38, "Given user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"avg(WEIGHT)\",\"header\":\"Adults\",\"area\":[[30.5,60],[60.5,60],[60.5,140],[30.5,140]]}]"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"avg(WEIGHT)\",\"header\":\"Adults\",\"area\":[[30.5,60],[60.5,60],[60.5,140],[30.5,140]]}]"]]);
    await session.step(42, "And user resizes line chart viewer to 800 by 500", () => resizeTo(page, el("line chart viewer"), 800, 500));
    await session.step(43, "Then the \"viewer regions\" reading of line chart viewer should be 1", () => readingIs(page, "viewer regions", el("line chart viewer"), 1));
    await session.step(44, "And the \"regions shown\" reading of line chart viewer should be 1", () => readingIs(page, "regions shown", el("line chart viewer"), 1));
    await session.step(45, "And line chart viewer should have a \"region Adults\" area", () => hasArea(page, el("line chart viewer"), "region Adults"));
    await session.step(46, "And line chart viewer should have a \"region Adults title\" area", () => hasArea(page, el("line chart viewer"), "region Adults title"));
    await session.step(47, "And the \"region Adults title\" area of line chart viewer should lie inside the \"region Adults\" area", () => areaLies(page, "region Adults title", el("line chart viewer"), "inside", "region Adults"));
    await session.step(48, "When user sets \"yColumnNames\" property of line chart viewer to \"WEIGHT, HEIGHT\"", () => setProperty(page, "yColumnNames", el("line chart viewer"), "WEIGHT, HEIGHT"));
    await session.step(49, "Then the \"charts\" reading of line chart viewer should be 2", () => readingIs(page, "charts", el("line chart viewer"), 2));
    await session.step(50, "And line chart viewer should have a \"region Adults\" area", () => hasArea(page, el("line chart viewer"), "region Adults"));
    await session.step(51, "And the \"region Adults title\" area of line chart viewer should lie inside the \"chart WEIGHT\" area", () => areaLies(page, "region Adults title", el("line chart viewer"), "inside", "chart WEIGHT"));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
  test("A line chart band constant on X reserves the strip above its charts", {tag: ["@viewers", "@realizes:viewers.line-chart", "@realizes:viewers.density-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(55, "Given user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT, HEIGHT"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Adults\",\"formula1\":\"${AGE} = 30\",\"formula2\":\"${AGE} = 60\"}]"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT, HEIGHT"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Adults\",\"formula1\":\"${AGE} = 30\",\"formula2\":\"${AGE} = 60\"}]"]]);
    await session.step(59, "Then the \"regions shown\" reading of line chart viewer should be 1", () => readingIs(page, "regions shown", el("line chart viewer"), 1));
    await session.step(60, "And line chart viewer should have a \"region Adults\" area", () => hasArea(page, el("line chart viewer"), "region Adults"));
    await session.step(61, "And line chart viewer should have a \"region Adults title\" area", () => hasArea(page, el("line chart viewer"), "region Adults title"));
    await session.step(62, "And the \"title strip top\" reading of line chart viewer should be at least 1", () => readingAtLeast(page, "title strip top", el("line chart viewer"), 1));
    await session.step(63, "And the \"region Adults title\" area of line chart viewer should lie above the \"region Adults\" area", () => areaLies(page, "region Adults title", el("line chart viewer"), "above", "region Adults"));
    await session.step(64, "And the \"region Adults title\" area of line chart viewer should lie inside the \"chart 1\" area", () => areaLies(page, "region Adults title", el("line chart viewer"), "inside", "chart 1"));
    await session.step(65, "And no errors should have been logged", () => noErrors(page));
  });
  test("A resized line chart keeps its band title", {tag: ["@viewers", "@realizes:viewers.line-chart", "@realizes:viewers.density-plot", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(69, "Given user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT, HEIGHT"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Adults\",\"formula1\":\"${AGE} = 30\",\"formula2\":\"${AGE} = 60\"}]"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT, HEIGHT"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Adults\",\"formula1\":\"${AGE} = 30\",\"formula2\":\"${AGE} = 60\"}]"]]);
      await session.step(73, "Then the \"regions shown\" reading of line chart viewer should be 1", () => readingIs(page, "regions shown", el("line chart viewer"), 1));
      await session.step(74, "And line chart viewer should have a \"region Adults\" area", () => hasArea(page, el("line chart viewer"), "region Adults"));
      await session.step(75, "When user resizes line chart viewer to 800 by 500", () => resizeTo(page, el("line chart viewer"), 800, 500));
      await session.step(76, "Then line chart viewer should have a \"region Adults\" area", () => hasArea(page, el("line chart viewer"), "region Adults"));
      await session.step(77, "And line chart viewer should have a \"region Adults title\" area", () => hasArea(page, el("line chart viewer"), "region Adults title"));
      await session.step(78, "And the \"title strip top\" reading of line chart viewer should be at least 1", () => readingAtLeast(page, "title strip top", el("line chart viewer"), 1));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
  });
  test("A density plot draws a titled area region inside its bins", {tag: ["@viewers", "@realizes:viewers.line-chart", "@realizes:viewers.density-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(82, "Given user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Adults\",\"area\":[[30.5,60],[60.5,60],[60.5,140],[30.5,140]]}]"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Adults\",\"area\":[[30.5,60],[60.5,60],[60.5,140],[30.5,140]]}]"]]);
    await session.step(86, "And user resizes density plot viewer to 800 by 500", () => resizeTo(page, el("density plot viewer"), 800, 500));
    await session.step(87, "Then the \"viewer regions\" reading of density plot viewer should be 1", () => readingIs(page, "viewer regions", el("density plot viewer"), 1));
    await session.step(88, "And the \"regions shown\" reading of density plot viewer should be 1", () => readingIs(page, "regions shown", el("density plot viewer"), 1));
    await session.step(89, "And density plot viewer should have a \"region Adults\" area", () => hasArea(page, el("density plot viewer"), "region Adults"));
    await session.step(90, "And density plot viewer should have a \"region Adults title\" area", () => hasArea(page, el("density plot viewer"), "region Adults title"));
    await session.step(91, "And the \"region Adults title\" area of density plot viewer should lie inside the \"region Adults\" area", () => areaLies(page, "region Adults title", el("density plot viewer"), "inside", "region Adults"));
    await session.step(92, "And the \"region Adults\" area of density plot viewer should lie inside the \"view\" area", () => areaLies(page, "region Adults", el("density plot viewer"), "inside", "view"));
    await session.step(93, "When user sets \"showViewerAnnotationRegions\" property of density plot viewer to \"false\"", () => setProperty(page, "showViewerAnnotationRegions", el("density plot viewer"), "false"));
    await session.step(94, "Then the \"regions shown\" reading of density plot viewer should be 0", () => readingIs(page, "regions shown", el("density plot viewer"), 0));
    await session.step(95, "And density plot viewer should not have a \"region Adults\" area", () => hasNoArea(page, el("density plot viewer"), "region Adults"));
    await session.step(96, "And density plot viewer should not have a \"region Adults title\" area", () => hasNoArea(page, el("density plot viewer"), "region Adults title"));
    await session.step(97, "When user sets \"showViewerAnnotationRegions\" property of density plot viewer to \"true\"", () => setProperty(page, "showViewerAnnotationRegions", el("density plot viewer"), "true"));
    await session.step(98, "Then the \"regions shown\" reading of density plot viewer should be 1", () => readingIs(page, "regions shown", el("density plot viewer"), 1));
    await session.step(99, "And density plot viewer should have a \"region Adults title\" area", () => hasArea(page, el("density plot viewer"), "region Adults title"));
    await session.step(100, "And no errors should have been logged", () => noErrors(page));
  });
});
