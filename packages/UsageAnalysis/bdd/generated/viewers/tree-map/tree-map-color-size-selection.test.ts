/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tree-map/tree-map-color-size-selection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tree-map]
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
import {addCategoricalFilter, clearSelection, filterPasses, filterPassesAll, noneSelected, onlyOfAnySelected, onlyOfSelected, selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaNarrower, areaShorter, areaWider, clickArea, clickAreaHolding, hasArea, hasNoArea, loadLayout, noErrors, painted, readingBetween, readingFinite, readingIs, readingReads, repainted, saveLayout, setProperties, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tree map colour, size, selection and filtering", () => {
  const session = feature(test, "features/viewers/tree-map/tree-map-color-size-selection.feature", import.meta.url);
  test("Tree map colour, size, selection and filtering", {tag: ["@journey", "@viewers", "@realizes:viewers.tree-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "And user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE"]]), [["splitByColumnNames","RACE"]]);
    await session.step(24, "Then the \"split columns\" reading of tree map viewer should be \"RACE\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE"));
    await session.step(25, "And the \"leaves\" reading of tree map viewer should be 4", () => readingIs(page, "leaves", el("tree map viewer"), 4));
    await session.step(26, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await session.step(27, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
    await session.step(28, "And tree map viewer should be painted", () => painted(page, el("tree map viewer")));
    await run.scenario("Colouring by AGE scores every leaf and the aggregation moves the score", async () => {
      await session.step(31, "Then the \"color column\" reading of tree map viewer should be \"\"", () => readingReads(page, "color column", el("tree map viewer"), ""));
      await session.step(32, "And the \"color of Caucasian\" reading of tree map viewer should be \"#2ca02c\"", () => readingReads(page, "color of Caucasian", el("tree map viewer"), "#2ca02c"));
      await session.step(33, "When user picks \"AGE\" in the \"color\" column selector of tree map viewer", () => pickInColumnSelector(page, "AGE", "color", el("tree map viewer")));
      await session.step(34, "Then the \"color column\" reading of tree map viewer should be \"AGE\"", () => readingReads(page, "color column", el("tree map viewer"), "AGE"));
      await session.step(35, "And the \"color aggregation\" reading of tree map viewer should be \"avg\"", () => readingReads(page, "color aggregation", el("tree map viewer"), "avg"));
      await session.step(36, "And the \"color score of Caucasian\" reading of tree map viewer should be between 45.6 and 45.7", () => readingBetween(page, "color score of Caucasian", el("tree map viewer"), 45.6, 45.7));
      await session.step(37, "And the \"color score of Asian\" reading of tree map viewer should be between 38.0 and 38.1", () => readingBetween(page, "color score of Asian", el("tree map viewer"), 38, 38.1));
      await session.step(38, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(39, "When user sets \"colorAggrType\" property of tree map viewer to \"max\"", () => setProperty(page, "colorAggrType", el("tree map viewer"), "max"));
      await session.step(40, "Then the \"color score of Caucasian\" reading of tree map viewer should be 89", () => readingIs(page, "color score of Caucasian", el("tree map viewer"), 89));
      await session.step(41, "And the \"color score of Asian\" reading of tree map viewer should be 64", () => readingIs(page, "color score of Asian", el("tree map viewer"), 64));
      await session.step(42, "And the \"color of Caucasian\" reading of tree map viewer should be \"#ff0000\"", () => readingReads(page, "color of Caucasian", el("tree map viewer"), "#ff0000"));
      await session.step(43, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(44, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(45, "When user sets properties of tree map viewer:", () => setProperties(page, el("tree map viewer"), [["colorColumnName",""],["colorAggrType","avg"]]), [["colorColumnName",""],["colorAggrType","avg"]]);
      await session.step(48, "Then the \"color column\" reading of tree map viewer should be \"\"", () => readingReads(page, "color column", el("tree map viewer"), ""));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sizing by WEIGHT replaces the row count with the aggregation", async () => {
      await session.step(52, "Then the \"size column\" reading of tree map viewer should be \"\"", () => readingReads(page, "size column", el("tree map viewer"), ""));
      await session.step(53, "And the \"area of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "area of Caucasian", el("tree map viewer"), 896));
      await session.step(54, "When user sets \"sizeColumnName\" property of tree map viewer to \"WEIGHT\"", () => setProperty(page, "sizeColumnName", el("tree map viewer"), "WEIGHT"));
      await session.step(55, "Then the \"size aggregation\" reading of tree map viewer should be \"sum\"", () => readingReads(page, "size aggregation", el("tree map viewer"), "sum"));
      await session.step(56, "And the \"area of Caucasian\" reading of tree map viewer should be between 71221 and 71222", () => readingBetween(page, "area of Caucasian", el("tree map viewer"), 71221, 71222));
      await session.step(57, "And the \"area of Asian\" reading of tree map viewer should be between 1057 and 1058", () => readingBetween(page, "area of Asian", el("tree map viewer"), 1057, 1058));
      await session.step(58, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(59, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(60, "When user sets \"sizeAggrType\" property of tree map viewer to \"max\"", () => setProperty(page, "sizeAggrType", el("tree map viewer"), "max"));
      await session.step(61, "Then the \"area of Caucasian\" reading of tree map viewer should be 165", () => readingIs(page, "area of Caucasian", el("tree map viewer"), 165));
      await session.step(62, "And the \"area of Asian\" reading of tree map viewer should be between 91.7 and 91.9", () => readingBetween(page, "area of Asian", el("tree map viewer"), 91.7, 91.9));
      await session.step(63, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(64, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(65, "When user sets properties of tree map viewer:", () => setProperties(page, el("tree map viewer"), [["sizeColumnName",""],["sizeAggrType","sum"]]), [["sizeColumnName",""],["sizeAggrType","sum"]]);
      await session.step(68, "Then the \"area of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "area of Caucasian", el("tree map viewer"), 896));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking a rectangle selects exactly the rows it holds, and the band is drawn over it", async () => {
      await session.step(72, "Given user clears the row selection", () => clearSelection(page));
      await session.step(73, "Then the \"selected rows of Caucasian\" reading of tree map viewer should be 0", () => readingIs(page, "selected rows of Caucasian", el("tree map viewer"), 0));
      await session.step(74, "And tree map viewer should not have a \"selection of Caucasian\" area", () => hasNoArea(page, el("tree map viewer"), "selection of Caucasian"));
      await session.step(75, "When user clicks on the \"leaf Caucasian\" area of tree map viewer", () => clickArea(page, "leaf Caucasian", el("tree map viewer")));
      await session.step(76, "Then 896 rows should be selected", () => selectedRowCount(page, 896));
      await session.step(77, "And only rows where \"RACE\" is \"Caucasian\" should be selected", () => onlyOfSelected(page, "RACE", "Caucasian"));
      await session.step(78, "And the \"selected rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "selected rows of Caucasian", el("tree map viewer"), 896));
      await session.step(79, "And the \"selected rows of Asian\" reading of tree map viewer should be 0", () => readingIs(page, "selected rows of Asian", el("tree map viewer"), 0));
      await session.step(80, "And tree map viewer should have a \"selection of Caucasian\" area", () => hasArea(page, el("tree map viewer"), "selection of Caucasian"));
      await session.step(81, "And tree map viewer should not have a \"selection of Asian\" area", () => hasNoArea(page, el("tree map viewer"), "selection of Asian"));
      await session.step(82, "When user clears the row selection", () => clearSelection(page));
      await session.step(83, "Then no rows should be selected", () => noneSelected(page));
      await session.step(84, "And tree map viewer should not have a \"selection of Caucasian\" area", () => hasNoArea(page, el("tree map viewer"), "selection of Caucasian"));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter on the table reshapes the map without re-splitting it", async () => {
      await session.step(88, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(89, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(90, "And the \"rows shown\" reading of tree map viewer should be 447", () => readingIs(page, "rows shown", el("tree map viewer"), 447));
      await session.step(91, "And the \"rows of Caucasian\" reading of tree map viewer should be 416", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 416));
      await session.step(92, "And the \"rows of Other\" reading of tree map viewer should be 14", () => readingIs(page, "rows of Other", el("tree map viewer"), 14));
      await session.step(93, "And the \"rows of Black\" reading of tree map viewer should be 9", () => readingIs(page, "rows of Black", el("tree map viewer"), 9));
      await session.step(94, "And the \"rows of Asian\" reading of tree map viewer should be 8", () => readingIs(page, "rows of Asian", el("tree map viewer"), 8));
      await session.step(95, "And the \"leaves\" reading of tree map viewer should be 4", () => readingIs(page, "leaves", el("tree map viewer"), 4));
      await session.step(96, "And the \"split columns\" reading of tree map viewer should be \"RACE\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE"));
      await session.step(97, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(98, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(99, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(100, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(101, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
      await session.step(102, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shift-click adds a rectangle's rows to the selection and Control-click takes them off", async () => {
      await session.step(106, "Given user clears the row selection", () => clearSelection(page));
      await session.step(107, "When user clicks on the \"leaf Asian\" area of tree map viewer", () => clickArea(page, "leaf Asian", el("tree map viewer")));
      await session.step(108, "Then 15 rows should be selected", () => selectedRowCount(page, 15));
      await session.step(109, "When user clicks on the \"leaf Black\" area of tree map viewer holding Shift", () => clickAreaHolding(page, "leaf Black", el("tree map viewer"), "Shift"));
      await session.step(110, "Then 42 rows should be selected", () => selectedRowCount(page, 42));
      await session.step(111, "And only rows where \"RACE\" is one of \"Asian, Black\" should be selected", () => onlyOfAnySelected(page, "RACE", "Asian, Black"));
      await session.step(112, "And the \"selected rows of Black\" reading of tree map viewer should be 27", () => readingIs(page, "selected rows of Black", el("tree map viewer"), 27));
      await session.step(113, "When user clicks on the \"leaf Asian\" area of tree map viewer holding Control", () => clickAreaHolding(page, "leaf Asian", el("tree map viewer"), "Control"));
      await session.step(114, "Then 27 rows should be selected", () => selectedRowCount(page, 27));
      await session.step(115, "And only rows where \"RACE\" is \"Black\" should be selected", () => onlyOfSelected(page, "RACE", "Black"));
      await session.step(116, "And the \"selected rows of Asian\" reading of tree map viewer should be 0", () => readingIs(page, "selected rows of Asian", el("tree map viewer"), 0));
      await session.step(117, "When user clears the row selection", () => clearSelection(page));
      await session.step(118, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source decides which rows the map is built from", async () => {
      await session.step(121, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(122, "Then the \"rows shown\" reading of tree map viewer should be 447", () => readingIs(page, "rows shown", el("tree map viewer"), 447));
      await session.step(123, "When user sets \"rowSource\" property of tree map viewer to \"All\"", () => setProperty(page, "rowSource", el("tree map viewer"), "All"));
      await session.step(124, "Then the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
      await session.step(125, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(126, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(127, "And user sets \"rowSource\" property of tree map viewer to \"Selected\"", () => setProperty(page, "rowSource", el("tree map viewer"), "Selected"));
      await session.step(128, "Then the \"rows shown\" reading of tree map viewer should be 15", () => readingIs(page, "rows shown", el("tree map viewer"), 15));
      await session.step(129, "And the \"leaves\" reading of tree map viewer should be 1", () => readingIs(page, "leaves", el("tree map viewer"), 1));
      await session.step(130, "When user sets \"rowSource\" property of tree map viewer to \"FilteredSelected\"", () => setProperty(page, "rowSource", el("tree map viewer"), "FilteredSelected"));
      await session.step(131, "Then the \"rows shown\" reading of tree map viewer should be 8", () => readingIs(page, "rows shown", el("tree map viewer"), 8));
      await session.step(132, "When user sets \"rowSource\" property of tree map viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("tree map viewer"), "Filtered"));
      await session.step(133, "Then the \"rows shown\" reading of tree map viewer should be 447", () => readingIs(page, "rows shown", el("tree map viewer"), 447));
      await session.step(134, "When user clears the row selection", () => clearSelection(page));
      await session.step(135, "And user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(136, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(137, "Then the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
      await session.step(138, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer's Filter formula narrows the map and leaves the table alone", async () => {
      await session.step(141, "When user sets \"filter\" property of tree map viewer to \"${AGE} > 40\"", () => setProperty(page, "filter", el("tree map viewer"), "${AGE} > 40"));
      await session.step(142, "Then the \"rows shown\" reading of tree map viewer should be 635", () => readingIs(page, "rows shown", el("tree map viewer"), 635));
      await session.step(143, "And the \"rows of Caucasian\" reading of tree map viewer should be 562", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 562));
      await session.step(144, "And the \"rows of Asian\" reading of tree map viewer should be 6", () => readingIs(page, "rows of Asian", el("tree map viewer"), 6));
      await session.step(145, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(146, "When user sets \"filter\" property of tree map viewer to \"\"", () => setProperty(page, "filter", el("tree map viewer"), ""));
      await session.step(147, "Then the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
      await session.step(148, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The outer margins inset the box the rectangles are laid out in", async () => {
      await session.step(151, "When user sets properties of tree map viewer:", () => setProperties(page, el("tree map viewer"), [["outerMarginLeft","30"],["outerMarginTop","30"],["outerMarginRight","30"],["outerMarginBottom","30"]]), [["outerMarginLeft","30"],["outerMarginTop","30"],["outerMarginRight","30"],["outerMarginBottom","30"]]);
      await session.step(156, "Then the \"view\" area of tree map viewer should be narrower than before", () => areaNarrower(page, "view", el("tree map viewer")));
      await session.step(157, "And the \"view\" area of tree map viewer should be shorter than before", () => areaShorter(page, "view", el("tree map viewer")));
      await session.step(158, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(159, "When user sets properties of tree map viewer:", () => setProperties(page, el("tree map viewer"), [["outerMarginLeft","0"],["outerMarginTop","0"],["outerMarginRight","0"],["outerMarginBottom","0"]]), [["outerMarginLeft","0"],["outerMarginTop","0"],["outerMarginRight","0"],["outerMarginBottom","0"]]);
      await session.step(164, "Then the \"view\" area of tree map viewer should be wider than before", () => areaWider(page, "view", el("tree map viewer")));
      await session.step(165, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the two levels and the colour column back", async () => {
      await session.step(168, "When user sets properties of tree map viewer:", () => setProperties(page, el("tree map viewer"), [["splitByColumnNames","RACE, SEX"],["colorColumnName","AGE"]]), [["splitByColumnNames","RACE, SEX"],["colorColumnName","AGE"]]);
      await session.step(171, "Then the \"levels\" reading of tree map viewer should be 2", () => readingIs(page, "levels", el("tree map viewer"), 2));
      await session.step(172, "And the \"color column\" reading of tree map viewer should be \"AGE\"", () => readingReads(page, "color column", el("tree map viewer"), "AGE"));
      await session.step(173, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(174, "And user clicks on close icon of tree map viewer", () => clickOn(page, el("close icon of tree map viewer")));
      await session.step(175, "Then the open tableview should have 0 tree map viewers", () => viewerCount(page, 0, "tree map"));
      await session.step(176, "When user loads the saved layout", () => loadLayout(page));
      await session.step(177, "Then the open tableview should have 1 tree map viewer", () => viewerCount(page, 1, "tree map"));
      await session.step(178, "And the \"split columns\" reading of tree map viewer should be \"RACE, SEX\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE, SEX"));
      await session.step(179, "And the \"color column\" reading of tree map viewer should be \"AGE\"", () => readingReads(page, "color column", el("tree map viewer"), "AGE"));
      await session.step(180, "And the \"rows of Caucasian | F\" reading of tree map viewer should be 480", () => readingIs(page, "rows of Caucasian | F", el("tree map viewer"), 480));
      await session.step(181, "And the \"color score of Caucasian | F\" reading of tree map viewer should be a finite number", () => readingFinite(page, "color score of Caucasian | F", el("tree map viewer")));
      await session.step(182, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
