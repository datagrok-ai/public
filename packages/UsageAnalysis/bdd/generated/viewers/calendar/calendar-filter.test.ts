/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/calendar/calendar-filter.feature
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
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {addCategoricalFilter, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, noErrors, painted, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The calendar filters the way the other viewers do", () => {
  const session = feature(test, "features/viewers/calendar/calendar-filter.feature", import.meta.url);
  test("The calendar filters the way the other viewers do", {tag: ["@journey", "@viewers", "@realizes:viewers.calendar", "@realizes:core#233"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a calendar viewer", () => addViewer(page, "calendar"));
    await session.step(13, "Then the \"date column\" reading of calendar viewer should be \"STARTED\"", () => readingReads(page, "date column", el("calendar viewer"), "STARTED"));
    await session.step(14, "And the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
    await session.step(15, "And calendar viewer should be painted", () => painted(page, el("calendar viewer")));
    await run.scenario("Closing the calendar lets go of its filter", async () => {
      await session.step(19, "When user sets \"onClick\" property of calendar viewer to \"Filter\"", () => setProperty(page, "onClick", el("calendar viewer"), "Filter"));
      await session.step(20, "And user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(21, "Then 32 rows should pass the filter", () => filterPasses(page, 32));
      await session.step(22, "When user clicks on close icon of calendar viewer", () => clickOn(page, el("close icon of calendar viewer")));
      await session.step(23, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(24, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second click takes the filter off, another part replaces it", async () => {
      await session.step(28, "When user adds a calendar viewer", () => addViewer(page, "calendar"));
      await session.step(29, "And user sets \"onClick\" property of calendar viewer to \"Filter\"", () => setProperty(page, "onClick", el("calendar viewer"), "Filter"));
      await session.step(30, "And user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(31, "Then 32 rows should pass the filter", () => filterPasses(page, 32));
      await session.step(32, "When user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(33, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(34, "When user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(35, "Then 32 rows should pass the filter", () => filterPasses(page, 32));
      await session.step(36, "When user clicks on the \"weekday Sunday\" area of calendar viewer", () => clickArea(page, "weekday Sunday", el("calendar viewer")));
      await session.step(37, "Then 153 rows should pass the filter", () => filterPasses(page, 153));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Leaving On Click = Filter and resetting the filter let go of it", async () => {
      await session.step(42, "When user sets \"onClick\" property of calendar viewer to \"Select\"", () => setProperty(page, "onClick", el("calendar viewer"), "Select"));
      await session.step(43, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(44, "When user sets \"onClick\" property of calendar viewer to \"Filter\"", () => setProperty(page, "onClick", el("calendar viewer"), "Filter"));
      await session.step(45, "And user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(46, "Then 32 rows should pass the filter", () => filterPasses(page, 32));
      await session.step(47, "When user picks \"View > Reset Filter\" from the top menu", () => pickFromTopMenu(page, "View > Reset Filter"));
      await session.step(48, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(49, "When user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(50, "Then 32 rows should pass the filter", () => filterPasses(page, 32));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The calendar's filter and the filter panel apply together, whichever comes first", async () => {
      await session.step(55, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(56, "Then 13 rows should pass the filter", () => filterPasses(page, 13));
      await session.step(57, "When user picks \"View > Reset Filter\" from the top menu", () => pickFromTopMenu(page, "View > Reset Filter"));
      await session.step(58, "And user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(59, "And user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(60, "Then 13 rows should pass the filter", () => filterPasses(page, 13));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
