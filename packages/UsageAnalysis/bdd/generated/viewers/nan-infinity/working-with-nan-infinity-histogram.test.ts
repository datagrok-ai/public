/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/nan-infinity/working-with-nan-infinity-histogram.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.histogram]
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
import {setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, painted, readingAsRemembered, readingFinite, readingIs, rememberReading, reportsNoError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A histogram over a column that holds a NaN and an Infinity", () => {
  const session = feature(test, "features/viewers/nan-infinity/working-with-nan-infinity-histogram.feature", import.meta.url);
  test("A histogram over a column that holds a NaN and an Infinity", {tag: ["@journey", "@viewers", "@realizes:viewers.histogram"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","WEIGHT"]]), [["valueColumnName","WEIGHT"]]);
    await session.step(15, "Then histogram viewer should report no error", () => reportsNoError(page, el("histogram viewer")));
    await session.step(16, "And the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
    await run.scenario("An Infinity in WEIGHT keeps the value axis finite and the bins as they were", async () => {
      await session.step(19, "When user remembers the \"axis max\" reading of histogram viewer", () => rememberReading(page, "axis max", el("histogram viewer")));
      await session.step(20, "And user sets \"WEIGHT\" column in row 2 to \"Infinity\"", () => setCell(page, "WEIGHT", 2, "Infinity"));
      await session.step(21, "Then the \"axis max\" reading of histogram viewer should be a finite number", () => readingFinite(page, "axis max", el("histogram viewer")));
      await session.step(22, "And the \"axis max\" reading of histogram viewer should be as remembered", () => readingAsRemembered(page, "axis max", el("histogram viewer")));
      await session.step(23, "And the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
      await session.step(24, "And histogram viewer should report no error", () => reportsNoError(page, el("histogram viewer")));
      await session.step(25, "And histogram viewer should be painted", () => painted(page, el("histogram viewer")));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A NaN in HEIGHT leaves the axis finite", async () => {
      await session.step(29, "When user sets \"valueColumnName\" property of histogram viewer to \"HEIGHT\"", () => setProperty(page, "valueColumnName", el("histogram viewer"), "HEIGHT"));
      await session.step(30, "Then the \"rows shown\" reading of histogram viewer should be 1000", () => readingIs(page, "rows shown", el("histogram viewer"), 1000));
      await session.step(31, "When user remembers the \"axis max\" reading of histogram viewer", () => rememberReading(page, "axis max", el("histogram viewer")));
      await session.step(32, "And user sets \"HEIGHT\" column in row 1 to \"NaN\"", () => setCell(page, "HEIGHT", 1, "NaN"));
      await session.step(33, "Then the \"axis min\" reading of histogram viewer should be a finite number", () => readingFinite(page, "axis min", el("histogram viewer")));
      await session.step(34, "And the \"axis max\" reading of histogram viewer should be as remembered", () => readingAsRemembered(page, "axis max", el("histogram viewer")));
      await session.step(35, "And histogram viewer should report no error", () => reportsNoError(page, el("histogram viewer")));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The two values put back leave the axis as it was", async () => {
      await session.step(39, "When user sets \"HEIGHT\" column in row 1 to \"174.705\"", () => setCell(page, "HEIGHT", 1, "174.705"));
      await session.step(40, "And user sets \"WEIGHT\" column in row 2 to \"64\"", () => setCell(page, "WEIGHT", 2, "64"));
      await session.step(41, "And user sets \"valueColumnName\" property of histogram viewer to \"WEIGHT\"", () => setProperty(page, "valueColumnName", el("histogram viewer"), "WEIGHT"));
      await session.step(42, "Then the \"axis max\" reading of histogram viewer should be a finite number", () => readingFinite(page, "axis max", el("histogram viewer")));
      await session.step(43, "And histogram viewer should report no error", () => reportsNoError(page, el("histogram viewer")));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
