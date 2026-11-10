/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-int-formats.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#249]
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
import {setGridColumnFormat} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("An int column shows the named number formats", () => {
  const session = feature(test, "features/viewers/grid/grid-int-formats.feature", import.meta.url);
  test("An int column shows the named number formats", {tag: ["@journey", "@viewers", "@realizes:viewers.grid", "@realizes:core#249"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await run.scenario("An int column takes the float format names", async () => {
      await session.step(12, "When user sets the format of \"AGE\" column in the grid to \"scientific\"", () => setGridColumnFormat(page, "AGE", "scientific"));
      await session.step(13, "Then the \"text of cell 1 of AGE\" reading of grid should be \"2.60E1\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "2.60E1"));
      await session.step(14, "When user sets the format of \"AGE\" column in the grid to \"two digits after comma\"", () => setGridColumnFormat(page, "AGE", "two digits after comma"));
      await session.step(15, "Then the \"text of cell 1 of AGE\" reading of grid should be \"26.00\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26.00"));
      await session.step(16, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
