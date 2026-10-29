/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tree-map/tree-map-edge-cases.feature
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
import {removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Tree map edge cases: a split column removed", () => {
  const session = feature(test, "features/viewers/tree-map/tree-map-edge-cases.feature", import.meta.url);
  test("A split column removed from the table is dropped from the split", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tree-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE, SEX"]]), [["splitByColumnNames","RACE, SEX"]]);
    await session.step(14, "Then the \"split columns\" reading of tree map viewer should be \"RACE, SEX\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE, SEX"));
    await session.step(15, "And the \"leaves\" reading of tree map viewer should be 8", () => readingIs(page, "leaves", el("tree map viewer"), 8));
    await session.step(16, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await session.step(19, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(20, "Then the \"split columns\" reading of tree map viewer should be \"RACE\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE"));
    await session.step(21, "And the \"leaves\" reading of tree map viewer should be 4", () => readingIs(page, "leaves", el("tree map viewer"), 4));
    await session.step(22, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
    await session.step(23, "And no errors should have been logged", () => noErrors(page));
  });
});
