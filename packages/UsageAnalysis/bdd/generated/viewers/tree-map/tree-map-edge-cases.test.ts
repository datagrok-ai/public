/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tree-map/tree-map-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tree-map, core#1135]
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
import {filterBetween, filterPasses, removeColumn, renameColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasNoArea, loadLayout, noErrors, readingIs, readingReads, saveLayout, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Tree map edge cases: no rows, split, color and size columns removed or renamed, a layout with a missing column", () => {
  const session = feature(test, "features/viewers/tree-map/tree-map-edge-cases.feature", import.meta.url);
  test("A filter that removes every row leaves the map empty without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tree-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE, SEX"]]), [["splitByColumnNames","RACE, SEX"]]);
    await session.step(18, "Then the \"split columns\" reading of tree map viewer should be \"RACE, SEX\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE, SEX"));
    await session.step(19, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
    await session.step(20, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await session.step(23, "When user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(24, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(25, "And the \"rows shown\" reading of tree map viewer should be 0", () => readingIs(page, "rows shown", el("tree map viewer"), 0));
    await session.step(26, "And tree map viewer should not have a \"leaf Caucasian | M\" area", () => hasNoArea(page, el("tree map viewer"), "leaf Caucasian | M"));
    await session.step(27, "When user resets the filter", () => resetFilter(page));
    await session.step(28, "Then the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await session.step(29, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
  });
  test("A split column removed from the table is dropped from the split", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tree-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE, SEX"]]), [["splitByColumnNames","RACE, SEX"]]);
    await session.step(18, "Then the \"split columns\" reading of tree map viewer should be \"RACE, SEX\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE, SEX"));
    await session.step(19, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
    await session.step(20, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await session.step(33, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(34, "Then the \"split columns\" reading of tree map viewer should be \"RACE\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE"));
    await session.step(35, "And the \"leaves\" reading of tree map viewer should be 4", () => readingIs(page, "leaves", el("tree map viewer"), 4));
    await session.step(36, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("A split column renamed is followed by the split", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tree-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE, SEX"]]), [["splitByColumnNames","RACE, SEX"]]);
    await session.step(18, "Then the \"split columns\" reading of tree map viewer should be \"RACE, SEX\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE, SEX"));
    await session.step(19, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
    await session.step(20, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await session.step(40, "When user renames \"SEX\" column to \"GENDER\"", () => renameColumn(page, "SEX", "GENDER"));
    await session.step(41, "Then the \"split columns\" reading of tree map viewer should be \"RACE, GENDER\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE, GENDER"));
    await session.step(42, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout naming a split column the table no longer has is split by the columns that remain", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tree-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE, SEX"]]), [["splitByColumnNames","RACE, SEX"]]);
    await session.step(18, "Then the \"split columns\" reading of tree map viewer should be \"RACE, SEX\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE, SEX"));
    await session.step(19, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
    await session.step(20, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await session.step(46, "When user saves the layout of the current table view", () => saveLayout(page));
    await session.step(47, "And user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(48, "And user loads the saved layout", () => loadLayout(page));
    await session.step(49, "Then the \"split columns\" reading of tree map viewer should be \"RACE\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE"));
    await session.step(50, "And the \"leaves\" reading of tree map viewer should be 4", () => readingIs(page, "leaves", el("tree map viewer"), 4));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the Color column keeps the map working without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tree-map", "@known-failure:60", "@realizes:core#1135"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE, SEX"]]), [["splitByColumnNames","RACE, SEX"]]);
    await session.step(18, "Then the \"split columns\" reading of tree map viewer should be \"RACE, SEX\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE, SEX"));
    await session.step(19, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
    await session.step(20, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await knownFailure(async () => {
      await session.step(57, "When user sets \"colorColumnName\" property of tree map viewer to \"AGE\"", () => setProperty(page, "colorColumnName", el("tree map viewer"), "AGE"));
      await session.step(58, "Then the \"color column\" reading of tree map viewer should be \"AGE\"", () => readingReads(page, "color column", el("tree map viewer"), "AGE"));
      await session.step(59, "When user removes \"AGE\" column", () => removeColumn(page, "AGE"));
      await session.step(60, "Then no errors should have been logged", () => noErrors(page));
      await session.step(61, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
    }, 60);
  });
  test("Removing the Size column keeps the map working without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tree-map", "@known-failure:70", "@realizes:core#1135"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE, SEX"]]), [["splitByColumnNames","RACE, SEX"]]);
    await session.step(18, "Then the \"split columns\" reading of tree map viewer should be \"RACE, SEX\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE, SEX"));
    await session.step(19, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
    await session.step(20, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await knownFailure(async () => {
      await session.step(67, "When user sets \"sizeColumnName\" property of tree map viewer to \"WEIGHT\"", () => setProperty(page, "sizeColumnName", el("tree map viewer"), "WEIGHT"));
      await session.step(68, "Then the \"size column\" reading of tree map viewer should be \"WEIGHT\"", () => readingReads(page, "size column", el("tree map viewer"), "WEIGHT"));
      await session.step(69, "When user removes \"WEIGHT\" column", () => removeColumn(page, "WEIGHT"));
      await session.step(70, "Then no errors should have been logged", () => noErrors(page));
      await session.step(71, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
    }, 70);
  });
});
