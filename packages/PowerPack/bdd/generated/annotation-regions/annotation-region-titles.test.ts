/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/annotation-regions/annotation-region-titles.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, viewers.box-plot]
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
import {clearSelection, onlyBetweenSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaShorter, areaTaller, clickArea, hasArea, hoverArea, narrowerRange, noErrors, pointerAway, readingAtLeast, readingIs, resizeTo, setProperties, setProperty, wheelOverArea, widerRange} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies, areaLiesBeside, areaPlacedAsRemembered, rememberAreaPlace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Annotation region titles", () => {
  const session = feature(test, "features/annotation-regions/annotation-region-titles.feature", import.meta.url);
  test("Annotation region titles", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["lassoTool","false"],["markerDefaultSize","2"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Adults\",\"area\":[[30.5,30],[60.5,30],[60.5,180],[30.5,180]]}]"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["lassoTool","false"],["markerDefaultSize","2"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Adults\",\"area\":[[30.5,30],[60.5,30],[60.5,180],[30.5,180]]}]"]]);
    await session.step(29, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(30, "Then the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
    await session.step(31, "And the \"region titles shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "region titles shown", el("scatter plot viewer"), 1));
    await run.scenario("An in-data title sits inside its region and reserves no strip", async () => {
      await session.step(34, "Then scatter plot viewer should have a \"region Adults\" area", () => hasArea(page, el("scatter plot viewer"), "region Adults"));
      await session.step(35, "And scatter plot viewer should have a \"region Adults title\" area", () => hasArea(page, el("scatter plot viewer"), "region Adults title"));
      await session.step(36, "And the \"region Adults title\" area of scatter plot viewer should lie inside the \"region Adults\" area", () => areaLies(page, "region Adults title", el("scatter plot viewer"), "inside", "region Adults"));
      await session.step(37, "And the \"title strip top\" reading of scatter plot viewer should be 0", () => readingIs(page, "title strip top", el("scatter plot viewer"), 0));
      await session.step(38, "And the \"title strip right\" reading of scatter plot viewer should be 0", () => readingIs(page, "title strip right", el("scatter plot viewer"), 0));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking the title selects the region's rows and leaves the title where it was", async () => {
      await session.step(42, "When user remembers the place of the \"region Adults title\" area of scatter plot viewer", () => rememberAreaPlace(page, "region Adults title", el("scatter plot viewer")));
      await session.step(43, "And user hovers over the \"region Adults title\" area of scatter plot viewer", () => hoverArea(page, "region Adults title", el("scatter plot viewer")));
      await session.step(44, "Then the \"regions hovered\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions hovered", el("scatter plot viewer"), 1));
      await session.step(45, "When user clicks on the \"region Adults title\" area of scatter plot viewer", () => clickArea(page, "region Adults title", el("scatter plot viewer")));
      await session.step(46, "Then only rows where \"AGE\" is between 31 and 60 should be selected", () => onlyBetweenSelected(page, "AGE", 31, 60));
      await session.step(47, "And the \"region Adults title\" area of scatter plot viewer should be placed as remembered", () => areaPlacedAsRemembered(page, "region Adults title", el("scatter plot viewer")));
      await session.step(48, "When user clears the row selection", () => clearSelection(page));
      await session.step(49, "And user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(50, "Then the \"region Adults title\" area of scatter plot viewer should be placed as remembered", () => areaPlacedAsRemembered(page, "region Adults title", el("scatter plot viewer")));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Growing the viewer and restoring it puts the title back where it was", async () => {
      await session.step(54, "When user remembers the place of the \"region Adults title\" area of scatter plot viewer", () => rememberAreaPlace(page, "region Adults title", el("scatter plot viewer")));
      await session.step(55, "And user resizes scatter plot viewer to 1000 by 600", () => resizeTo(page, el("scatter plot viewer"), 1000, 600));
      await session.step(56, "Then the \"region Adults title\" area of scatter plot viewer should lie inside the \"region Adults\" area", () => areaLies(page, "region Adults title", el("scatter plot viewer"), "inside", "region Adults"));
      await session.step(57, "When user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
      await session.step(58, "Then the \"region Adults title\" area of scatter plot viewer should be placed as remembered", () => areaPlacedAsRemembered(page, "region Adults title", el("scatter plot viewer")));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Zooming keeps the title inside the region", async () => {
      await session.step(62, "When user scrolls the mouse wheel up over the \"region Adults\" area of scatter plot viewer", () => wheelOverArea(page, "up", "region Adults", el("scatter plot viewer")));
      await session.step(63, "Then scatter plot viewer should show a narrower value range than before", () => narrowerRange(page, el("scatter plot viewer")));
      await session.step(64, "And the \"region Adults title\" area of scatter plot viewer should lie inside the \"region Adults\" area", () => areaLies(page, "region Adults title", el("scatter plot viewer"), "inside", "region Adults"));
      await session.step(65, "When user scrolls the mouse wheel down over the \"region Adults\" area of scatter plot viewer", () => wheelOverArea(page, "down", "region Adults", el("scatter plot viewer")));
      await session.step(66, "Then scatter plot viewer should show a wider value range than before", () => widerRange(page, el("scatter plot viewer")));
      await session.step(67, "And the \"region Adults title\" area of scatter plot viewer should lie inside the \"region Adults\" area", () => areaLies(page, "region Adults title", el("scatter plot viewer"), "inside", "region Adults"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A larger annotation font draws a taller title", async () => {
      await session.step(71, "When user sets \"annotationFont\" property of scatter plot viewer to \"normal normal 20px \\\"Roboto\\\"\"", () => setProperty(page, "annotationFont", el("scatter plot viewer"), "normal normal 20px \"Roboto\""));
      await session.step(72, "Then the \"region Adults title\" area of scatter plot viewer should be taller than before", () => areaTaller(page, "region Adults title", el("scatter plot viewer")));
      await session.step(73, "And the \"region Adults title\" area of scatter plot viewer should lie inside the \"region Adults\" area", () => areaLies(page, "region Adults title", el("scatter plot viewer"), "inside", "region Adults"));
      await session.step(74, "When user sets \"annotationFont\" property of scatter plot viewer to \"normal normal 10px \\\"Roboto\\\"\"", () => setProperty(page, "annotationFont", el("scatter plot viewer"), "normal normal 10px \"Roboto\""));
      await session.step(75, "Then the \"region Adults title\" area of scatter plot viewer should be shorter than before", () => areaShorter(page, "region Adults title", el("scatter plot viewer")));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A band on the X column whose title fits takes the strip above the plot", async () => {
      await session.step(79, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["annotationRegions","[{\"type\":\"formula\",\"header\":\"Adults\",\"formula1\":\"${AGE} = 30\",\"formula2\":\"${AGE} = 60\"}]"]]), [["annotationRegions","[{\"type\":\"formula\",\"header\":\"Adults\",\"formula1\":\"${AGE} = 30\",\"formula2\":\"${AGE} = 60\"}]"]]);
      await session.step(81, "Then the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
      await session.step(82, "And the \"region titles shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "region titles shown", el("scatter plot viewer"), 1));
      await session.step(83, "And the \"title strip top\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "title strip top", el("scatter plot viewer"), 1));
      await session.step(84, "And the \"title strip right\" reading of scatter plot viewer should be 0", () => readingIs(page, "title strip right", el("scatter plot viewer"), 0));
      await session.step(85, "And the \"region Adults title\" area of scatter plot viewer should lie above the \"view\" area", () => areaLies(page, "region Adults title", el("scatter plot viewer"), "above", "view"));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A band on the Y column takes the strip to the right of the plot", async () => {
      await session.step(89, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["annotationRegions","[{\"type\":\"formula\",\"header\":\"Middle\",\"formula1\":\"${WEIGHT} = 80\",\"formula2\":\"${WEIGHT} = 120\"}]"]]), [["annotationRegions","[{\"type\":\"formula\",\"header\":\"Middle\",\"formula1\":\"${WEIGHT} = 80\",\"formula2\":\"${WEIGHT} = 120\"}]"]]);
      await session.step(91, "Then the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
      await session.step(92, "And the \"title strip right\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "title strip right", el("scatter plot viewer"), 1));
      await session.step(93, "And the \"title strip top\" reading of scatter plot viewer should be 0", () => readingIs(page, "title strip top", el("scatter plot viewer"), 0));
      await session.step(94, "And the \"region Middle title\" area of scatter plot viewer should lie to the right of the \"view\" area", () => areaLiesBeside(page, "region Middle title", el("scatter plot viewer"), "right", "view"));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Region titles stay drawn after the formula lines go empty", async () => {
      await session.step(98, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Adults\",\"area\":[[30.5,30],[60.5,30],[60.5,180],[30.5,180]]}]"],["formulaLines","[{\"type\":\"band\",\"title\":\"Young\",\"formula\":\"${AGE} in (18, 28)\",\"orientation\":\"Vertical\",\"column2\":\"WEIGHT\"}]"]]), [["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Adults\",\"area\":[[30.5,30],[60.5,30],[60.5,180],[30.5,180]]}]"],["formulaLines","[{\"type\":\"band\",\"title\":\"Young\",\"formula\":\"${AGE} in (18, 28)\",\"orientation\":\"Vertical\",\"column2\":\"WEIGHT\"}]"]]);
      await session.step(101, "Then the \"formula lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "formula lines", el("scatter plot viewer"), 1));
      await session.step(102, "And the \"title strip top\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "title strip top", el("scatter plot viewer"), 1));
      await session.step(103, "And the \"region titles shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "region titles shown", el("scatter plot viewer"), 1));
      await session.step(104, "When user sets \"formulaLines\" property of scatter plot viewer to \"[]\"", () => setProperty(page, "formulaLines", el("scatter plot viewer"), "[]"));
      await session.step(105, "Then the \"formula lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "formula lines", el("scatter plot viewer"), 0));
      await session.step(106, "And the \"title strip top\" reading of scatter plot viewer should be 0", () => readingIs(page, "title strip top", el("scatter plot viewer"), 0));
      await session.step(107, "And the \"region titles shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "region titles shown", el("scatter plot viewer"), 1));
      await session.step(108, "And scatter plot viewer should have a \"region Adults title\" area", () => hasArea(page, el("scatter plot viewer"), "region Adults title"));
      await session.step(109, "And the \"region Adults title\" area of scatter plot viewer should lie inside the \"region Adults\" area", () => areaLies(page, "region Adults title", el("scatter plot viewer"), "inside", "region Adults"));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A box plot keeps a band title in place through a plot style and a bin count change", async () => {
      await session.step(113, "Given user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["category1ColumnName","RACE"],["valueColumnName","AGE"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Prime\",\"formula1\":\"${AGE} = 30\",\"formula2\":\"${AGE} = 50\"}]"]]), [["category1ColumnName","RACE"],["valueColumnName","AGE"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Prime\",\"formula1\":\"${AGE} = 30\",\"formula2\":\"${AGE} = 50\"}]"]]);
      await session.step(117, "Then the \"regions shown\" reading of box plot viewer should be 1", () => readingIs(page, "regions shown", el("box plot viewer"), 1));
      await session.step(118, "And box plot viewer should have a \"region Prime\" area", () => hasArea(page, el("box plot viewer"), "region Prime"));
      await session.step(119, "And box plot viewer should have a \"region Prime title\" area", () => hasArea(page, el("box plot viewer"), "region Prime title"));
      await session.step(120, "And the \"title strip right\" reading of box plot viewer should be at least 1", () => readingAtLeast(page, "title strip right", el("box plot viewer"), 1));
      await session.step(121, "When user remembers the place of the \"region Prime title\" area of box plot viewer", () => rememberAreaPlace(page, "region Prime title", el("box plot viewer")));
      await session.step(122, "And user sets \"plotStyle\" property of box plot viewer to \"violin\"", () => setProperty(page, "plotStyle", el("box plot viewer"), "violin"));
      await session.step(123, "Then the \"region Prime title\" area of box plot viewer should be placed as remembered", () => areaPlacedAsRemembered(page, "region Prime title", el("box plot viewer")));
      await session.step(124, "When user sets \"bins\" property of box plot viewer to \"20\"", () => setProperty(page, "bins", el("box plot viewer"), "20"));
      await session.step(125, "Then the \"region Prime title\" area of box plot viewer should be placed as remembered", () => areaPlacedAsRemembered(page, "region Prime title", el("box plot viewer")));
      await session.step(126, "When user sets \"plotStyle\" property of box plot viewer to \"box\"", () => setProperty(page, "plotStyle", el("box plot viewer"), "box"));
      await session.step(127, "Then the \"region Prime title\" area of box plot viewer should be placed as remembered", () => areaPlacedAsRemembered(page, "region Prime title", el("box plot viewer")));
      await session.step(128, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A band title that does not fit renders in the data and reserves nothing", async () => {
      await session.step(131, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["annotationRegions","[{\"type\":\"formula\",\"header\":\"Adults of the study population\",\"formula1\":\"${AGE} = 40\",\"formula2\":\"${AGE} = 45\"}]"]]), [["annotationRegions","[{\"type\":\"formula\",\"header\":\"Adults of the study population\",\"formula1\":\"${AGE} = 40\",\"formula2\":\"${AGE} = 45\"}]"]]);
      await session.step(133, "Then the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
      await session.step(134, "And the \"region titles shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "region titles shown", el("scatter plot viewer"), 1));
      await session.step(135, "And the \"title strip top\" reading of scatter plot viewer should be 0", () => readingIs(page, "title strip top", el("scatter plot viewer"), 0));
      await session.step(136, "And the \"region Adults of the study population title\" area of scatter plot viewer should lie inside the \"view\" area", () => areaLies(page, "region Adults of the study population title", el("scatter plot viewer"), "inside", "view"));
      await session.step(137, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Auto Layout off drops the strip and the title moves into the data", async () => {
      await session.step(140, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["annotationRegions","[{\"type\":\"formula\",\"header\":\"Adults\",\"formula1\":\"${AGE} = 30\",\"formula2\":\"${AGE} = 60\"}]"],["autoLayout","true"]]), [["annotationRegions","[{\"type\":\"formula\",\"header\":\"Adults\",\"formula1\":\"${AGE} = 30\",\"formula2\":\"${AGE} = 60\"}]"],["autoLayout","true"]]);
      await session.step(143, "Then the \"title strip top\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "title strip top", el("scatter plot viewer"), 1));
      await session.step(144, "When user sets \"autoLayout\" property of scatter plot viewer to \"false\"", () => setProperty(page, "autoLayout", el("scatter plot viewer"), "false"));
      await session.step(145, "Then the \"title strip top\" reading of scatter plot viewer should be 0", () => readingIs(page, "title strip top", el("scatter plot viewer"), 0));
      await session.step(146, "And the \"region Adults title\" area of scatter plot viewer should lie inside the \"view\" area", () => areaLies(page, "region Adults title", el("scatter plot viewer"), "inside", "view"));
      await session.step(147, "When user sets \"autoLayout\" property of scatter plot viewer to \"true\"", () => setProperty(page, "autoLayout", el("scatter plot viewer"), "true"));
      await session.step(148, "Then the \"title strip top\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "title strip top", el("scatter plot viewer"), 1));
      await session.step(149, "And the \"region Adults title\" area of scatter plot viewer should lie above the \"view\" area", () => areaLies(page, "region Adults title", el("scatter plot viewer"), "above", "view"));
      await session.step(150, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
