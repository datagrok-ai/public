/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/annotation-regions/annotation-region-interaction.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
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
import {pressKey} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, filterTo, noneOfSelected, noneSelected, onlyBetweenSelected, renameColumn, resetFilter, selectedPassFilter, setTableTag, setTableTagText, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, clickAreaHolding, hasArea, hasNoArea, hoverArea, noErrors, pickFromContextMenu, pointerAway, propertyShouldContain, readingIs, readingReads, resizeTo, setProperties, setProperty, showsFewerRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Annotation region interaction", () => {
  const session = feature(test, "features/annotation-regions/annotation-region-interaction.feature", import.meta.url);
  test("Annotation region interaction", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["lassoTool","false"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Outer\",\"area\":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Inner\",\"area\":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}]"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["lassoTool","false"],["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Outer\",\"area\":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Inner\",\"area\":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}]"]]);
    await session.step(23, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(24, "Then the \"viewer regions\" reading of scatter plot viewer should be 2", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 2));
    await session.step(25, "And the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
    await session.step(26, "And scatter plot viewer should have a \"region Outer\" area", () => hasArea(page, el("scatter plot viewer"), "region Outer"));
    await session.step(27, "And scatter plot viewer should have a \"region Inner\" area", () => hasArea(page, el("scatter plot viewer"), "region Inner"));
    await run.scenario("Hovering counts the regions under the pointer", async () => {
      await session.step(30, "When user hovers over the \"region Outer\" area of scatter plot viewer", () => hoverArea(page, "region Outer", el("scatter plot viewer")));
      await session.step(31, "Then the \"regions hovered\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions hovered", el("scatter plot viewer"), 1));
      await session.step(32, "When user hovers over the \"region Inner\" area of scatter plot viewer", () => hoverArea(page, "region Inner", el("scatter plot viewer")));
      await session.step(33, "Then the \"regions hovered\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions hovered", el("scatter plot viewer"), 2));
      await session.step(34, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(35, "Then the \"regions hovered\" reading of scatter plot viewer should be 0", () => readingIs(page, "regions hovered", el("scatter plot viewer"), 0));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click selects the rows the region holds", async () => {
      await session.step(39, "When user clicks on the \"region Outer\" area of scatter plot viewer", () => clickArea(page, "region Outer", el("scatter plot viewer")));
      await session.step(40, "Then only rows where \"AGE\" is between 21 and 60 should be selected", () => onlyBetweenSelected(page, "AGE", 21, 60));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click where regions overlap selects the rows they share", async () => {
      await session.step(44, "When user clicks on the \"region Inner\" area of scatter plot viewer", () => clickArea(page, "region Inner", el("scatter plot viewer")));
      await session.step(45, "Then only rows where \"AGE\" is between 31 and 38 should be selected", () => onlyBetweenSelected(page, "AGE", 31, 38));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control drops a selected region's rows and Shift adds a region's rows", async () => {
      await session.step(49, "When user clicks on the \"region Inner\" area of scatter plot viewer", () => clickArea(page, "region Inner", el("scatter plot viewer")));
      await session.step(50, "Then only rows where \"AGE\" is between 31 and 38 should be selected", () => onlyBetweenSelected(page, "AGE", 31, 38));
      await session.step(51, "When user clicks on the \"region Inner\" area of scatter plot viewer holding Control", () => clickAreaHolding(page, "region Inner", el("scatter plot viewer"), "Control"));
      await session.step(52, "Then no rows should be selected", () => noneSelected(page));
      await session.step(53, "When user clicks on the \"region Inner\" area of scatter plot viewer", () => clickArea(page, "region Inner", el("scatter plot viewer")));
      await session.step(54, "And user clicks on the \"region Outer\" area of scatter plot viewer holding Shift", () => clickAreaHolding(page, "region Outer", el("scatter plot viewer"), "Shift"));
      await session.step(55, "Then only rows where \"AGE\" is between 21 and 60 should be selected", () => onlyBetweenSelected(page, "AGE", 21, 60));
      await session.step(56, "When user clears the row selection", () => clearSelection(page));
      await session.step(57, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click selects nothing the filter has taken out", async () => {
      await session.step(60, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(61, "And user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(62, "Then scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(63, "When user hovers over the \"region Outer\" area of scatter plot viewer", () => hoverArea(page, "region Outer", el("scatter plot viewer")));
      await session.step(64, "Then the \"regions hovered\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions hovered", el("scatter plot viewer"), 1));
      await session.step(65, "When user clicks on the \"region Outer\" area of scatter plot viewer", () => clickArea(page, "region Outer", el("scatter plot viewer")));
      await session.step(66, "Then some rows should be selected", () => someSelected(page));
      await session.step(67, "And every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(68, "And no rows where \"SEX\" is \"M\" should be selected", () => noneOfSelected(page, "SEX", "M"));
      await session.step(69, "When user resets the filter", () => resetFilter(page));
      await session.step(70, "And user clears the row selection", () => clearSelection(page));
      await session.step(71, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A dataframe region is drawn next to the viewer's and hidden on its own switch", async () => {
      await session.step(74, "When user sets the \".annotation-regions\" tag of the table to:", () => setTableTagText(page, ".annotation-regions", "[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Shared\",\"area\":[[65.5,30],[85.5,30],[85.5,180],[65.5,180]]}]"));
      await session.step(78, "Then the \"dataframe regions\" reading of scatter plot viewer should be 1", () => readingIs(page, "dataframe regions", el("scatter plot viewer"), 1));
      await session.step(79, "And the \"regions shown\" reading of scatter plot viewer should be 3", () => readingIs(page, "regions shown", el("scatter plot viewer"), 3));
      await session.step(80, "And scatter plot viewer should have a \"region Shared\" area", () => hasArea(page, el("scatter plot viewer"), "region Shared"));
      await session.step(81, "When user sets \"showDataframeAnnotationRegions\" property of scatter plot viewer to \"false\"", () => setProperty(page, "showDataframeAnnotationRegions", el("scatter plot viewer"), "false"));
      await session.step(82, "Then the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
      await session.step(83, "And scatter plot viewer should not have a \"region Shared\" area", () => hasNoArea(page, el("scatter plot viewer"), "region Shared"));
      await session.step(84, "And scatter plot viewer should have a \"region Outer\" area", () => hasArea(page, el("scatter plot viewer"), "region Outer"));
      await session.step(85, "When user sets \"showDataframeAnnotationRegions\" property of scatter plot viewer to \"true\"", () => setProperty(page, "showDataframeAnnotationRegions", el("scatter plot viewer"), "true"));
      await session.step(86, "Then scatter plot viewer should have a \"region Shared\" area", () => hasArea(page, el("scatter plot viewer"), "region Shared"));
      await session.step(87, "When user sets the \".annotation-regions\" tag of the table to \"\"", () => setTableTag(page, ".annotation-regions", ""));
      await session.step(88, "Then the \"dataframe regions\" reading of scatter plot viewer should be 0", () => readingIs(page, "dataframe regions", el("scatter plot viewer"), 0));
      await session.step(89, "And the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
      await session.step(90, "And scatter plot viewer should not have a \"region Shared\" area", () => hasNoArea(page, el("scatter plot viewer"), "region Shared"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A hidden region is kept in the look but not drawn", async () => {
      await session.step(94, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Outer\",\"area\":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Inner\",\"hidden\":true,\"area\":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}]"]]), [["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Outer\",\"area\":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Inner\",\"hidden\":true,\"area\":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}]"]]);
      await session.step(96, "Then the \"viewer regions\" reading of scatter plot viewer should be 2", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 2));
      await session.step(97, "And the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
      await session.step(98, "And scatter plot viewer should have a \"region Outer\" area", () => hasArea(page, el("scatter plot viewer"), "region Outer"));
      await session.step(99, "And scatter plot viewer should not have a \"region Inner\" area", () => hasNoArea(page, el("scatter plot viewer"), "region Inner"));
      await session.step(100, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Outer\",\"area\":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Inner\",\"area\":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}]"]]), [["annotationRegions","[{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Outer\",\"area\":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{\"type\":\"area\",\"x\":\"AGE\",\"y\":\"WEIGHT\",\"header\":\"Inner\",\"area\":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}]"]]);
      await session.step(102, "Then the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Swapped axes keep the regions and another column drops them", async () => {
      await session.step(106, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xColumnName","WEIGHT"],["yColumnName","AGE"]]), [["xColumnName","WEIGHT"],["yColumnName","AGE"]]);
      await session.step(109, "Then the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
      await session.step(110, "And scatter plot viewer should have a \"region Outer\" area", () => hasArea(page, el("scatter plot viewer"), "region Outer"));
      await session.step(111, "When user sets \"xColumnName\" property of scatter plot viewer to \"HEIGHT\"", () => setProperty(page, "xColumnName", el("scatter plot viewer"), "HEIGHT"));
      await session.step(112, "Then scatter plot viewer should not have a \"region Outer\" area", () => hasNoArea(page, el("scatter plot viewer"), "region Outer"));
      await session.step(113, "And scatter plot viewer should not have a \"region Inner\" area", () => hasNoArea(page, el("scatter plot viewer"), "region Inner"));
      await session.step(114, "And the \"region titles shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "region titles shown", el("scatter plot viewer"), 0));
      await session.step(115, "And the \"viewer regions\" reading of scatter plot viewer should be 2", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 2));
      await session.step(116, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]);
      await session.step(119, "Then the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
      await session.step(120, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column rename follows into the regions", async () => {
      await session.step(123, "When user renames \"AGE\" column to \"AGE (years)\"", () => renameColumn(page, "AGE", "AGE (years)"));
      await session.step(124, "Then \"annotationRegions\" property of scatter plot viewer should contain \"\\\"x\\\":\\\"AGE (years)\\\"\"", () => propertyShouldContain(page, "annotationRegions", el("scatter plot viewer"), "\"x\":\"AGE (years)\""));
      await session.step(125, "And the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
      await session.step(126, "When user renames \"AGE (years)\" column to \"AGE\"", () => renameColumn(page, "AGE (years)", "AGE"));
      await session.step(127, "Then \"annotationRegions\" property of scatter plot viewer should contain \"\\\"x\\\":\\\"AGE\\\"\"", () => propertyShouldContain(page, "annotationRegions", el("scatter plot viewer"), "\"x\":\"AGE\""));
      await session.step(128, "And the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
      await session.step(129, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Escape leaves the drawing mode with the regions untouched", async () => {
      await session.step(132, "When user sets \"showViewerAnnotationRegions\" property of scatter plot viewer to \"false\"", () => setProperty(page, "showViewerAnnotationRegions", el("scatter plot viewer"), "false"));
      await session.step(133, "And user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(134, "And user picks \"Tools > Draw Annotation Region\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tools > Draw Annotation Region", el("scatter plot viewer")));
      await session.step(135, "Then the \"region drawing mode\" reading of scatter plot viewer should be \"true\"", () => readingReads(page, "region drawing mode", el("scatter plot viewer"), "true"));
      await session.step(136, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(137, "Then the \"region drawing mode\" reading of scatter plot viewer should be \"false\"", () => readingReads(page, "region drawing mode", el("scatter plot viewer"), "false"));
      await session.step(138, "And the \"viewer regions\" reading of scatter plot viewer should be 2", () => readingIs(page, "viewer regions", el("scatter plot viewer"), 2));
      await session.step(139, "When user sets \"showViewerAnnotationRegions\" property of scatter plot viewer to \"true\"", () => setProperty(page, "showViewerAnnotationRegions", el("scatter plot viewer"), "true"));
      await session.step(140, "Then the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
      await session.step(141, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
