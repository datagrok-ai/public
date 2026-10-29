/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/matrix-plot/matrix-plot-edge-cases.feature
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
import {addCalculated, deleteSelected, renameColumn, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, noErrors, propertyShouldContain, readingIs, readingReads, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Matrix plot edge cases — no X columns, every row deleted, a constant column, a renamed column", () => {
  const session = feature(test, "features/viewers/matrix-plot/matrix-plot-edge-cases.feature", import.meta.url);
  test("With X Column Names cleared the plot draws no cell and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.matrix-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a matrix plot viewer", () => addViewer(page, "matrix plot"));
    await session.step(17, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
    await session.step(18, "And the \"cells drawn\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 16));
    await session.step(21, "When user sets \"xColumnNames\" property of matrix plot viewer to \"\"", () => setProperty(page, "xColumnNames", el("matrix plot viewer"), ""));
    await session.step(22, "Then the \"x columns\" reading of matrix plot viewer should be 0", () => readingIs(page, "x columns", el("matrix plot viewer"), 0));
    await session.step(23, "And matrix plot viewer should not have a \"cell 1,1\" area", () => hasNoArea(page, el("matrix plot viewer"), "cell 1,1"));
    await session.step(24, "And the \"error\" reading of matrix plot viewer should be \"\"", () => readingReads(page, "error", el("matrix plot viewer"), ""));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
    await session.step(26, "When user sets \"xColumnNames\" property of matrix plot viewer to \"AGE, HEIGHT\"", () => setProperty(page, "xColumnNames", el("matrix plot viewer"), "AGE, HEIGHT"));
    await session.step(27, "Then the \"cells drawn\" reading of matrix plot viewer should be 8", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 8));
  });
  test("Once every row of the table is deleted the plot logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.matrix-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a matrix plot viewer", () => addViewer(page, "matrix plot"));
    await session.step(17, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
    await session.step(18, "And the \"cells drawn\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 16));
    await session.step(30, "When user selects all rows", () => selectAllRows(page));
    await session.step(31, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(32, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(33, "And the \"rows shown\" reading of matrix plot viewer should be 0", () => readingIs(page, "rows shown", el("matrix plot viewer"), 0));
    await session.step(34, "And the \"error\" reading of matrix plot viewer should be \"\"", () => readingReads(page, "error", el("matrix plot viewer"), ""));
    await session.step(35, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column holding one value in every row draws its cells without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.matrix-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a matrix plot viewer", () => addViewer(page, "matrix plot"));
    await session.step(17, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
    await session.step(18, "And the \"cells drawn\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 16));
    await session.step(38, "When user adds a calculated column \"ONE_VALUE\" with formula \"5\"", () => addCalculated(page, "ONE_VALUE", "5"));
    await session.step(39, "And user sets properties of matrix plot viewer:", () => setProperties(page, el("matrix plot viewer"), [["xColumnNames","AGE, ONE_VALUE"],["yColumnNames","AGE, ONE_VALUE"]]), [["xColumnNames","AGE, ONE_VALUE"],["yColumnNames","AGE, ONE_VALUE"]]);
    await session.step(42, "Then the \"cells drawn\" reading of matrix plot viewer should be 4", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 4));
    await session.step(43, "And the \"cell viewer type of ONE_VALUE x ONE_VALUE\" reading of matrix plot viewer should be \"Histogram\"", () => readingReads(page, "cell viewer type of ONE_VALUE x ONE_VALUE", el("matrix plot viewer"), "Histogram"));
    await session.step(44, "And the \"cell rows shown of ONE_VALUE x AGE\" reading of matrix plot viewer should be 1000", () => readingIs(page, "cell rows shown of ONE_VALUE x AGE", el("matrix plot viewer"), 1000));
    await session.step(45, "And the \"error\" reading of matrix plot viewer should be \"\"", () => readingReads(page, "error", el("matrix plot viewer"), ""));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column renamed while the plot shows it is renamed in the plot's columns", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.matrix-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a matrix plot viewer", () => addViewer(page, "matrix plot"));
    await session.step(17, "Then the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
    await session.step(18, "And the \"cells drawn\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 16));
    await session.step(49, "Then matrix plot viewer should have a \"cell HEIGHT x AGE\" area", () => hasArea(page, el("matrix plot viewer"), "cell HEIGHT x AGE"));
    await session.step(50, "When user renames \"HEIGHT\" column to \"HEIGHT_CM\"", () => renameColumn(page, "HEIGHT", "HEIGHT_CM"));
    await session.step(51, "Then \"xColumnNames\" property of matrix plot viewer should contain \"HEIGHT_CM\"", () => propertyShouldContain(page, "xColumnNames", el("matrix plot viewer"), "HEIGHT_CM"));
    await session.step(52, "And matrix plot viewer should not have a \"cell HEIGHT x AGE\" area", () => hasNoArea(page, el("matrix plot viewer"), "cell HEIGHT x AGE"));
    await session.step(53, "And the \"cells drawn\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells drawn", el("matrix plot viewer"), 16));
    await session.step(54, "And matrix plot viewer should have a \"cell HEIGHT_CM x AGE\" area", () => hasArea(page, el("matrix plot viewer"), "cell HEIGHT_CM x AGE"));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
  });
});
