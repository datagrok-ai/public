/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/calendar/calendar-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.calendar, core#1221, core#1222]
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
import {columnType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {deleteSelected, filterPasses, openEmptyFilterPanel, removeColumn, rowCount, selectAllRows, setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {addViewer, clickArea, hasNoArea, noErrors, readingIs, readingReads, reportsError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
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
  test("The calendar's click filter stays when a category is toggled in the filter panel, and goes when On Click is switched to Select", {tag: ["@viewers", "@realizes:viewers.calendar", "@edge-cases", "@known-failure:33", "@realizes:core#1221"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a calendar viewer", () => addViewer(page, "calendar"));
    await session.step(12, "Then the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
    await knownFailure(async () => {
      await session.step(27, "When user sets \"onClick\" property of calendar viewer to \"Filter\"", () => setProperty(page, "onClick", el("calendar viewer"), "Filter"));
      await session.step(28, "And user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")));
      await session.step(29, "Then 32 rows should pass the filter", () => filterPasses(page, 32));
      await session.step(30, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(31, "And user adds a card for \"SEX\" to the filter panel", () => addCardFor(page, "SEX"));
      await session.step(32, "And user clicks on the \"category F of SEX\" area of filter panel", () => clickArea(page, "category F of SEX", el("filter panel")));
      await session.step(33, "Then 19 rows should pass the filter", () => filterPasses(page, 19));
      await session.step(34, "When user sets \"onClick\" property of calendar viewer to \"Select\"", () => setProperty(page, "onClick", el("calendar viewer"), "Select"));
      await session.step(35, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    }, 33);
  });
  test("A calendar over a date column whose values are all empty shows a \"No data\" message instead of a grid", {tag: ["@viewers", "@realizes:viewers.calendar", "@edge-cases", "@known-failure:51", "@realizes:core#1222"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a calendar viewer", () => addViewer(page, "calendar"));
    await session.step(12, "Then the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
    await knownFailure(async () => {
      await session.step(41, "Given user opens a table \"empty-dates\" with:", () => openTableOf(page, "empty-dates", [["When","Value"],["2024-01-05","1"],["2024-01-06","2"]]), [["When","Value"],["2024-01-05","1"],["2024-01-06","2"]]);
      await session.step(45, "When user sets \"When\" column in row 1 to \"\"", () => setCell(page, "When", 1, ""));
      await session.step(46, "And user sets \"When\" column in row 2 to \"\"", () => setCell(page, "When", 2, ""));
      await session.step(47, "Then \"When\" column should have type \"datetime\"", () => columnType(page, "When", "datetime"));
      await session.step(48, "When user adds a calendar viewer", () => addViewer(page, "calendar"));
      await session.step(49, "Then the \"date column\" reading of calendar viewer should be \"When\"", () => readingReads(page, "date column", el("calendar viewer"), "When"));
      await session.step(50, "And the \"rows shown\" reading of calendar viewer should be 0", () => readingIs(page, "rows shown", el("calendar viewer"), 0));
      await session.step(51, "And calendar viewer should report the error \"No data\"", () => reportsError(page, el("calendar viewer"), "No data"));
      await session.step(52, "And calendar viewer should not have a \"days\" area", () => hasNoArea(page, el("calendar viewer"), "days"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    }, 51);
  });
  test("Removing the date column of an open calendar shows a \"No data\" message instead of a grid", {tag: ["@viewers", "@realizes:viewers.calendar", "@edge-cases", "@known-failure:59", "@realizes:core#1222"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a calendar viewer", () => addViewer(page, "calendar"));
    await session.step(12, "Then the \"rows shown\" reading of calendar viewer should be 1000", () => readingIs(page, "rows shown", el("calendar viewer"), 1000));
    await knownFailure(async () => {
      await session.step(58, "When user removes \"STARTED\" column", () => removeColumn(page, "STARTED"));
      await session.step(59, "Then calendar viewer should report the error \"No data\"", () => reportsError(page, el("calendar viewer"), "No data"));
      await session.step(60, "And calendar viewer should not have a \"days\" area", () => hasNoArea(page, el("calendar viewer"), "days"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    }, 59);
  });
});
