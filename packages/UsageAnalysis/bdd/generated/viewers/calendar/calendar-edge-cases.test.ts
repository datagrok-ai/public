/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/calendar/calendar-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.calendar, core#233]
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
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {deleteSelected, filterPasses, filterPassesAll, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, hasNoArea, noErrors, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Calendar edge cases", () => {
  const session = feature(test, "features/viewers/calendar/calendar-edge-cases.feature", import.meta.url);
  test("A calendar over a table with every row deleted counts no rows and logs no error", {tag: ["@viewers", "@realizes:viewers.calendar", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a calendar viewer", () => addViewer(page, "calendar"));
    await session.step(12, "Then the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
    await session.step(16, "When user selects all rows", () => selectAllRows(page));
    await session.step(17, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(18, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(19, "And the \"rows shown\" reading of calendar viewer should be 0", () => readingIs(page, "rows shown", el("calendar viewer"), 0));
    await session.step(20, "And calendar viewer should not have a \"busiest day\" area", () => hasNoArea(page, el("calendar viewer"), "busiest day"));
    await session.step(21, "And no errors should have been logged", () => noErrors(page));
  });
  test("Closing a calendar releases the filter its click applied", {tag: ["@viewers", "@realizes:viewers.calendar", "@edge-cases", "@known-failure:33", "@realizes:core#233"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a calendar viewer", () => addViewer(page, "calendar"));
    await session.step(12, "Then the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
    await knownFailure(async () => {
      await session.step(28, "When user sets \"onClick\" property of calendar viewer to \"Filter\"", () => setProperty(page, "onClick", el("calendar viewer"), "Filter"));
      await session.step(29, "And user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(30, "Then 32 rows should pass the filter", () => filterPasses(page, 32));
      await session.step(31, "When user clicks on close icon of calendar viewer", () => clickOn(page, el("close icon of calendar viewer")));
      await session.step(32, "Then calendar viewer should be absent", () => shouldBe(page, el("calendar viewer"), "absent"));
      await session.step(33, "And all rows should pass the filter", () => filterPassesAll(page));
    }, 33);
  });
});
