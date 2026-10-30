/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/correlation-plot/correlation-plot-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.correlation-plot, core#1172]
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
import {removeColumn, renameColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, noErrors, readingReads, setProperties} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Correlation plot edge cases", () => {
  const session = feature(test, "features/viewers/correlation-plot/correlation-plot-edge-cases.feature", import.meta.url);
  test("A column renamed while the plot is open keeps its cells and coefficients under the new name", {tag: ["@viewers", "@realizes:viewers.correlation-plot", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a correlation plot viewer", () => addViewer(page, "correlation plot"));
    await session.step(14, "Then the \"x columns\" reading of correlation plot viewer should be \"AGE, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "x columns", el("correlation plot viewer"), "AGE, HEIGHT, WEIGHT, STARTED"));
    await session.step(15, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
    await session.step(19, "When user renames \"AGE\" column to \"AGE2\"", () => renameColumn(page, "AGE", "AGE2"));
    await session.step(20, "Then the \"x columns\" reading of correlation plot viewer should be \"AGE2, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "x columns", el("correlation plot viewer"), "AGE2, HEIGHT, WEIGHT, STARTED"));
    await session.step(21, "And the \"y columns\" reading of correlation plot viewer should be \"AGE2, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "y columns", el("correlation plot viewer"), "AGE2, HEIGHT, WEIGHT, STARTED"));
    await session.step(22, "And correlation plot viewer should have a \"cell HEIGHT x AGE2\" area", () => hasArea(page, el("correlation plot viewer"), "cell HEIGHT x AGE2"));
    await session.step(23, "And correlation plot viewer should not have a \"cell HEIGHT x AGE\" area", () => hasNoArea(page, el("correlation plot viewer"), "cell HEIGHT x AGE"));
    await session.step(24, "And the \"text of cell HEIGHT x AGE2\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE2", el("correlation plot viewer"), "-0.23"));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column removed from narrowed axes keeps the other narrowed columns", {tag: ["@viewers", "@realizes:viewers.correlation-plot", "@edge-cases", "@realizes:core#1172"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a correlation plot viewer", () => addViewer(page, "correlation plot"));
    await session.step(14, "Then the \"x columns\" reading of correlation plot viewer should be \"AGE, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "x columns", el("correlation plot viewer"), "AGE, HEIGHT, WEIGHT, STARTED"));
    await session.step(15, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
    await session.step(29, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT"],["yColumnNames","AGE, HEIGHT, WEIGHT"]]), [["xColumnNames","AGE, HEIGHT, WEIGHT"],["yColumnNames","AGE, HEIGHT, WEIGHT"]]);
    await session.step(32, "Then the \"x columns\" reading of correlation plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => readingReads(page, "x columns", el("correlation plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(33, "When user removes \"WEIGHT\" column", () => removeColumn(page, "WEIGHT"));
    await session.step(34, "Then the \"x columns\" reading of correlation plot viewer should be \"AGE, HEIGHT\"", () => readingReads(page, "x columns", el("correlation plot viewer"), "AGE, HEIGHT"));
    await session.step(35, "And the \"y columns\" reading of correlation plot viewer should be \"AGE, HEIGHT\"", () => readingReads(page, "y columns", el("correlation plot viewer"), "AGE, HEIGHT"));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
  });
});
