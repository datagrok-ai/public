/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-inner-grid-drag.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pivot-table]
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragAreaOntoWidget} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — a key column header dragged inside the inner grid", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-inner-grid-drag.feature", import.meta.url);
  test("Dragging the DIS_POP header onto None avg(AGE) logs no error", {tag: ["@viewers", "@realizes:viewers.pivot-table", "@known-failure:20"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(13, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(14, "And the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
    await knownFailure(async () => {
      await session.step(19, "When user drags the \"grid header DIS_POP\" area of pivot table viewer onto the \"grid header None avg(AGE)\" area of pivot table viewer", () => dragAreaOntoWidget(page, "grid header DIS_POP", el("pivot table viewer"), "grid header None avg(AGE)", el("pivot table viewer")));
      await session.step(20, "Then no errors should have been logged", () => noErrors(page));
    }, 20);
  });
});
