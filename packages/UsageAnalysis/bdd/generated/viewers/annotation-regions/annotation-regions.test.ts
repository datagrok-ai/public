/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/annotation-regions/annotation-regions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, viewers.histogram, viewers.line-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, closeContextMenu, dragAcrossArea, hasArea, hasNoArea, menuDoesNotList, menuLists, noErrors, openContextMenu, pickFromAreaContextMenu, pickFromContextMenu, propertyShouldBe, propertyShouldContain, readingIs, readingReads, repainted, rightClickArea, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Annotation regions", () => {
  const session = feature(test, "features/viewers/annotation-regions/annotation-regions.feature", import.meta.url);
  test("Annotation regions", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.histogram", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(25, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["lassoTool","false"]]));
    await session.step(29, "Then the \"viewer regions\" reading of scatter plot viewer should be 0", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 0));
    await session.step(30, "And the \"regions shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "regions shown", el("scatter plot viewer"), 0));
    await session.step(31, "And scatter plot viewer should not have a \"region 1\" area", () => hasNoArea(page, el("scatter plot viewer"), "region 1"));
    await run.scenario("Drawing a rectangle adds a region the viewer paints and reports", async () => {
      await session.step(34, "When user picks \"Tools > Draw Annotation Region\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tools > Draw Annotation Region", el("scatter plot viewer")));
      await session.step(35, "Then the \"region drawing mode\" reading of scatter plot viewer should be \"true\"", () => readingReads(page, "region drawing mode", el("scatter plot viewer"), "true"));
      await session.step(36, "When user drags across the \"view\" area of scatter plot viewer", () => dragAcrossArea(page, "view", el("scatter plot viewer")));
      await session.step(37, "Then the \"viewer regions\" reading of scatter plot viewer should be 1", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 1));
      await session.step(38, "And the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
      await session.step(39, "And the \"dataframe regions\" reading of scatter plot viewer should be 0", () => readingIs(page, "dataframe regions", el("scatter plot viewer"), 0));
      await session.step(40, "And scatter plot viewer should have a \"region 1\" area", () => hasArea(page, el("scatter plot viewer"), "region 1"));
      await session.step(41, "And the \"region 1\" area of scatter plot viewer should be painted", () => areaPainted(page, "region 1", el("scatter plot viewer")));
      await session.step(42, "And scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(43, "And \"annotationRegions\" property of scatter plot viewer should contain \"area\"", () => propertyShouldContain(page, "annotationRegions", el("scatter plot viewer"), "area"));
      await session.step(44, "When user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
      await session.step(45, "Then the \"viewer regions\" reading of scatter plot viewer should be 1", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 1));
      await session.step(46, "And the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hiding the viewer's regions takes the region off the plot but not out of the look", async () => {
      await session.step(50, "When user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")));
      await session.step(51, "And user sets \"showViewerAnnotationRegions\" property of scatter plot viewer to \"false\"", () => setProperty(page, "showViewerAnnotationRegions", el("scatter plot viewer"), "false"));
      await session.step(52, "Then the \"viewer regions\" reading of scatter plot viewer should be 1", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 1));
      await session.step(53, "And the \"regions shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "regions shown", el("scatter plot viewer"), 0));
      await session.step(54, "And scatter plot viewer should not have a \"region 1\" area", () => hasNoArea(page, el("scatter plot viewer"), "region 1"));
      await session.step(55, "And scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(56, "When user sets \"showViewerAnnotationRegions\" property of scatter plot viewer to \"true\"", () => setProperty(page, "showViewerAnnotationRegions", el("scatter plot viewer"), "true"));
      await session.step(57, "Then the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
      await session.step(58, "And scatter plot viewer should have a \"region 1\" area", () => hasArea(page, el("scatter plot viewer"), "region 1"));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Annotation Regions in the Tools menu switches both kinds at once", async () => {
      await session.step(62, "When user picks \"Tools > Show Annotation Regions\" from the context menu of the \"empty space\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tools > Show Annotation Regions", "empty space", el("scatter plot viewer")));
      await session.step(63, "Then \"showViewerAnnotationRegions\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "showViewerAnnotationRegions", el("scatter plot viewer"), "false"));
      await session.step(64, "And \"showDataframeAnnotationRegions\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "showDataframeAnnotationRegions", el("scatter plot viewer"), "false"));
      await session.step(65, "And the \"regions shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "regions shown", el("scatter plot viewer"), 0));
      await session.step(66, "When user picks \"Tools > Show Annotation Regions\" from the context menu of the \"empty space\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Tools > Show Annotation Regions", "empty space", el("scatter plot viewer")));
      await session.step(67, "Then \"showViewerAnnotationRegions\" property of scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "showViewerAnnotationRegions", el("scatter plot viewer"), "true"));
      await session.step(68, "And the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A histogram locks the categorical axis and stores the band as a formula region", async () => {
      await session.step(72, "Given user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"]]));
      await session.step(74, "Then the \"viewer regions\" reading of histogram viewer should be 0", () => readingIs(page, "viewer regions", el("histogram viewer"), 0));
      await session.step(75, "When user picks \"Tools > Draw Annotation Region\" from the context menu of histogram viewer", () => pickFromContextMenu(page, "Tools > Draw Annotation Region", el("histogram viewer")));
      await session.step(76, "And user drags across the \"view\" area of histogram viewer", () => dragAcrossArea(page, "view", el("histogram viewer")));
      await session.step(77, "Then the \"viewer regions\" reading of histogram viewer should be 1", () => readingIs(page, "viewer regions", el("histogram viewer"), 1));
      await session.step(78, "And the \"regions shown\" reading of histogram viewer should be 1", () => readingIs(page, "regions shown", el("histogram viewer"), 1));
      await session.step(79, "And histogram viewer should have a \"region 1\" area", () => hasArea(page, el("histogram viewer"), "region 1"));
      await session.step(80, "And \"annotationRegions\" property of histogram viewer should contain \"formula\"", () => propertyShouldContain(page, "annotationRegions", el("histogram viewer"), "formula"));
      await session.step(81, "And \"annotationRegions\" property of histogram viewer should contain \"${AGE}\"", () => propertyShouldContain(page, "annotationRegions", el("histogram viewer"), "${AGE}"));
      await session.step(82, "When user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
      await session.step(83, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The value axis offers Add Line, Add Band and Add Region", async () => {
      await session.step(86, "When user right-clicks on the \"x axis\" area of histogram viewer", () => rightClickArea(page, "x axis", el("histogram viewer")));
      await session.step(87, "Then the open menu should list \"Annotations > Add Line\"", () => menuLists(page, "Annotations > Add Line"));
      await session.step(88, "And the open menu should list \"Annotations > Add Band\"", () => menuLists(page, "Annotations > Add Band"));
      await session.step(89, "And the open menu should list \"Annotations > Add Region\"", () => menuLists(page, "Annotations > Add Region"));
      await session.step(90, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(91, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A multi-axis line chart offers no Draw Annotation Region", async () => {
      await session.step(94, "Given user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT, HEIGHT"],["multiAxis","false"]]));
      await session.step(98, "When user opens the context menu of line chart viewer", () => openContextMenu(page, el("line chart viewer")));
      await session.step(99, "Then the open menu should list \"Tools > Draw Annotation Region\"", () => menuLists(page, "Tools > Draw Annotation Region"));
      await session.step(100, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(101, "And user sets \"multiAxis\" property of line chart viewer to \"true\"", () => setProperty(page, "multiAxis", el("line chart viewer"), "true"));
      await session.step(102, "And user opens the context menu of line chart viewer", () => openContextMenu(page, el("line chart viewer")));
      await session.step(103, "Then the open menu should not list \"Tools > Draw Annotation Region\"", () => menuDoesNotList(page, "Tools > Draw Annotation Region"));
      await session.step(104, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(105, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
