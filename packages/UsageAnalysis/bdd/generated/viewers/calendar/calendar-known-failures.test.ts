/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/calendar/calendar-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.calendar, core#1221]
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
import {filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {addViewer, clickArea, propertyShouldBe, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Calendar known failures", () => {
  const session = feature(test, "features/viewers/calendar/calendar-known-failures.feature", import.meta.url);
  test("The calendar's click filter stays when a category is toggled in the filter panel (core#1221)", {tag: ["@viewers", "@realizes:viewers.calendar", "@edge-cases", "@known-failure:25", "@realizes:core#1221"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a calendar viewer", () => addViewer(page, "calendar"));
    await session.step(13, "Then the \"date column\" reading of calendar viewer should be \"STARTED\"", () => readingReads(page, "date column", el("calendar viewer"), "STARTED"));
    await session.step(14, "And \"On Click\" property of calendar viewer should be \"Select\"", () => propertyShouldBe(page, "On Click", el("calendar viewer"), "Select"));
    await session.step(15, "And all rows should pass the filter", () => filterPassesAll(page));
    await knownFailure(async () => {
      await session.step(19, "When user sets \"onClick\" property of calendar viewer to \"Filter\"", () => setProperty(page, "onClick", el("calendar viewer"), "Filter"));
      await session.step(20, "And user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(21, "Then 32 rows should pass the filter", () => filterPasses(page, 32));
      await session.step(22, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(23, "And user adds a card for \"SEX\" to the filter panel", () => addCardFor(page, "SEX"));
      await session.step(24, "And user clicks on the \"category F of SEX\" area of filter panel", () => clickArea(page, "category F of SEX", el("filter panel")));
      await session.step(25, "Then 19 rows should pass the filter", () => filterPasses(page, 19));
    }, 25);
  });
});
