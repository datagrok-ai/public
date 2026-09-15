/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tree-map/tree-map-color-size-selection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tree-map]
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
import {clickOn, hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, clearSelection, filterPasses, noneSelected, onlyOfSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hasArea, hasNoArea, noErrors, painted, readingBetween, readingIs, readingReads, repainted, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tree map colour, size, selection and filtering", () => {
  const session = feature(test, "features/viewers/tree-map/tree-map-color-size-selection.feature", import.meta.url);
  test("Tree map colour, size, selection and filtering", {tag: ["@journey", "@viewers", "@realizes:viewers.tree-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE"]]));
    await session.step(19, "Then the \"split columns\" reading of tree map viewer should be \"RACE\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE"));
    await session.step(20, "And the \"leaves\" reading of tree map viewer should be 4", () => readingIs(page, "leaves", el("tree map viewer"), 4));
    await session.step(21, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await session.step(22, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
    await session.step(23, "And tree map viewer should be painted", () => painted(page, el("tree map viewer")));
    await run.scenario("Colouring by AGE scores every leaf and the aggregation moves the score", async () => {
      await session.step(26, "Then the \"color column\" reading of tree map viewer should be \"\"", () => readingReads(page, "color column", el("tree map viewer"), ""));
      await session.step(27, "And the \"color of Caucasian\" reading of tree map viewer should be \"#2ca02c\"", () => readingReads(page, "color of Caucasian", el("tree map viewer"), "#2ca02c"));
      await session.step(28, "When user picks \"AGE\" in the \"color\" column selector of tree map viewer", () => pickInColumnSelector(page, "AGE", "color", el("tree map viewer")));
      await session.step(29, "Then the \"color column\" reading of tree map viewer should be \"AGE\"", () => readingReads(page, "color column", el("tree map viewer"), "AGE"));
      await session.step(30, "And the \"color aggregation\" reading of tree map viewer should be \"avg\"", () => readingReads(page, "color aggregation", el("tree map viewer"), "avg"));
      await session.step(31, "And the \"color score of Caucasian\" reading of tree map viewer should be between 45.6 and 45.7", () => readingBetween(page, "color score of Caucasian", el("tree map viewer"), 45.6, 45.7));
      await session.step(32, "And the \"color score of Asian\" reading of tree map viewer should be between 38.0 and 38.1", () => readingBetween(page, "color score of Asian", el("tree map viewer"), 38, 38.1));
      await session.step(33, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(34, "When user sets \"colorAggrType\" property of tree map viewer to \"max\"", () => setProperty(page, "colorAggrType", el("tree map viewer"), "max"));
      await session.step(35, "Then the \"color score of Caucasian\" reading of tree map viewer should be 89", () => readingIs(page, "color score of Caucasian", el("tree map viewer"), 89));
      await session.step(36, "And the \"color score of Asian\" reading of tree map viewer should be 64", () => readingIs(page, "color score of Asian", el("tree map viewer"), 64));
      await session.step(37, "And the \"color of Caucasian\" reading of tree map viewer should be \"#ff0000\"", () => readingReads(page, "color of Caucasian", el("tree map viewer"), "#ff0000"));
      await session.step(38, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(39, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(40, "When user sets properties of tree map viewer:", () => setProperties(page, el("tree map viewer"), [["colorColumnName",""],["colorAggrType","avg"]]));
      await session.step(43, "Then the \"color column\" reading of tree map viewer should be \"\"", () => readingReads(page, "color column", el("tree map viewer"), ""));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sizing by WEIGHT replaces the row count with the aggregation", async () => {
      await session.step(47, "Then the \"size column\" reading of tree map viewer should be \"\"", () => readingReads(page, "size column", el("tree map viewer"), ""));
      await session.step(48, "And the \"area of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "area of Caucasian", el("tree map viewer"), 896));
      await session.step(49, "When user sets \"sizeColumnName\" property of tree map viewer to \"WEIGHT\"", () => setProperty(page, "sizeColumnName", el("tree map viewer"), "WEIGHT"));
      await session.step(50, "Then the \"size aggregation\" reading of tree map viewer should be \"sum\"", () => readingReads(page, "size aggregation", el("tree map viewer"), "sum"));
      await session.step(51, "And the \"area of Caucasian\" reading of tree map viewer should be between 71221 and 71222", () => readingBetween(page, "area of Caucasian", el("tree map viewer"), 71221, 71222));
      await session.step(52, "And the \"area of Asian\" reading of tree map viewer should be between 1057 and 1058", () => readingBetween(page, "area of Asian", el("tree map viewer"), 1057, 1058));
      await session.step(53, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(54, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(55, "When user sets \"sizeAggrType\" property of tree map viewer to \"max\"", () => setProperty(page, "sizeAggrType", el("tree map viewer"), "max"));
      await session.step(56, "Then the \"area of Caucasian\" reading of tree map viewer should be 165", () => readingIs(page, "area of Caucasian", el("tree map viewer"), 165));
      await session.step(57, "And the \"area of Asian\" reading of tree map viewer should be between 91.7 and 91.9", () => readingBetween(page, "area of Asian", el("tree map viewer"), 91.7, 91.9));
      await session.step(58, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(59, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(60, "When user sets properties of tree map viewer:", () => setProperties(page, el("tree map viewer"), [["sizeColumnName",""],["sizeAggrType","sum"]]));
      await session.step(63, "Then the \"area of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "area of Caucasian", el("tree map viewer"), 896));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clicking a rectangle selects exactly the rows it holds, and the band is drawn over it", async () => {
      await session.step(67, "Given user clears the row selection", () => clearSelection(page));
      await session.step(68, "Then the \"selected rows of Caucasian\" reading of tree map viewer should be 0", () => readingIs(page, "selected rows of Caucasian", el("tree map viewer"), 0));
      await session.step(69, "And tree map viewer should not have a \"selection of Caucasian\" area", () => hasNoArea(page, el("tree map viewer"), "selection of Caucasian"));
      await session.step(70, "When user clicks on the \"leaf Caucasian\" area of tree map viewer", () => clickArea(page, "leaf Caucasian", el("tree map viewer")));
      await session.step(71, "Then 896 rows should be selected", () => selectedRowCount(page, 896));
      await session.step(72, "And only rows where \"RACE\" is \"Caucasian\" should be selected", () => onlyOfSelected(page, "RACE", "Caucasian"));
      await session.step(73, "And the \"selected rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "selected rows of Caucasian", el("tree map viewer"), 896));
      await session.step(74, "And the \"selected rows of Asian\" reading of tree map viewer should be 0", () => readingIs(page, "selected rows of Asian", el("tree map viewer"), 0));
      await session.step(75, "And tree map viewer should have a \"selection of Caucasian\" area", () => hasArea(page, el("tree map viewer"), "selection of Caucasian"));
      await session.step(76, "And tree map viewer should not have a \"selection of Asian\" area", () => hasNoArea(page, el("tree map viewer"), "selection of Asian"));
      await session.step(77, "When user clears the row selection", () => clearSelection(page));
      await session.step(78, "Then no rows should be selected", () => noneSelected(page));
      await session.step(79, "And tree map viewer should not have a \"selection of Caucasian\" area", () => hasNoArea(page, el("tree map viewer"), "selection of Caucasian"));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter on the table reshapes the map without re-splitting it", async () => {
      await session.step(83, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(84, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(85, "And the \"rows shown\" reading of tree map viewer should be 447", () => readingIs(page, "rows shown", el("tree map viewer"), 447));
      await session.step(86, "And the \"rows of Caucasian\" reading of tree map viewer should be 416", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 416));
      await session.step(87, "And the \"rows of Other\" reading of tree map viewer should be 14", () => readingIs(page, "rows of Other", el("tree map viewer"), 14));
      await session.step(88, "And the \"rows of Black\" reading of tree map viewer should be 9", () => readingIs(page, "rows of Black", el("tree map viewer"), 9));
      await session.step(89, "And the \"rows of Asian\" reading of tree map viewer should be 8", () => readingIs(page, "rows of Asian", el("tree map viewer"), 8));
      await session.step(90, "And the \"leaves\" reading of tree map viewer should be 4", () => readingIs(page, "leaves", el("tree map viewer"), 4));
      await session.step(91, "And the \"split columns\" reading of tree map viewer should be \"RACE\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE"));
      await session.step(92, "And tree map viewer should have repainted", () => repainted(page, el("tree map viewer")));
      await session.step(93, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(94, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(95, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(96, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
      await session.step(97, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(98, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
