/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/form/form-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.form, core#1138]
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
import {uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, makeLastRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Form viewer edge cases", () => {
  const session = feature(test, "features/viewers/form/form-edge-cases.feature", import.meta.url);
  test("The row arrows stop at the first and the last row", {tag: ["@viewers", "@realizes:viewers.form", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(15, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await session.step(19, "When user clicks on the \"next row\" area of form viewer", () => clickArea(page, "next row", el("form viewer")));
    await session.step(20, "Then row 2 should be current", () => currentRowIs(page, 2));
    await session.step(21, "When user clicks on the \"previous row\" area of form viewer", () => clickArea(page, "previous row", el("form viewer")));
    await session.step(22, "Then row 1 should be current", () => currentRowIs(page, 1));
    await session.step(23, "When user clicks on the \"previous row\" area of form viewer", () => clickArea(page, "previous row", el("form viewer")));
    await session.step(24, "Then row 1 should be current", () => currentRowIs(page, 1));
    await session.step(25, "And the \"row\" reading of form viewer should be 1", () => readingIs(page, "row", el("form viewer"), 1));
    await session.step(26, "When user makes the last row current", () => makeLastRowCurrent(page));
    await session.step(27, "Then the \"row\" reading of form viewer should be 1000", () => readingIs(page, "row", el("form viewer"), 1000));
    await session.step(28, "When user clicks on the \"next row\" area of form viewer", () => clickArea(page, "next row", el("form viewer")));
    await session.step(29, "Then row 1000 should be current", () => currentRowIs(page, 1000));
    await session.step(30, "And the \"row\" reading of form viewer should be 1000", () => readingIs(page, "row", el("form viewer"), 1000));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
  test("A file that is not JSON loaded as a form leaves the form as it was and logs no error", {tag: ["@viewers", "@realizes:viewers.form", "@edge-cases", "@known-failure:39", "@realizes:core#1138"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(15, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await knownFailure(async () => {
      await session.step(37, "When user uploads \"fixtures/cars-small.csv\" through folder open icon in form viewer", () => uploadThrough(page, "fixtures/cars-small.csv", el("folder open icon in form viewer")));
      await session.step(38, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    }, 39);
  });
  test("A JSON file that is not a form loaded as a form leaves the form as it was", {tag: ["@viewers", "@realizes:viewers.form", "@edge-cases", "@known-failure:47", "@realizes:core#1138"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(15, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await knownFailure(async () => {
      await session.step(46, "When user uploads \"fixtures/nx/spgi-legend-position-old.layout\" through folder open icon in form viewer", () => uploadThrough(page, "fixtures/nx/spgi-legend-position-old.layout", el("folder open icon in form viewer")));
      await session.step(47, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    }, 47);
  });
});
