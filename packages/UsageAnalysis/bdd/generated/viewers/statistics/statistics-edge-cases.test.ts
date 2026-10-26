/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/statistics/statistics-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.stats-viewer]
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
import {addCalculated, filterBetween, filterPasses, removeColumn, renameColumn, resetFilter, selectNoRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, noErrors, readingIs, readingReads, reportsNoError, resizeTo, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Statistics edge cases: no rows counted, one row, columns added, removed and renamed", () => {
  const session = feature(test, "features/viewers/statistics/statistics-edge-cases.feature", import.meta.url);
  test("A filter that removes every row leaves no values counted and no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(18, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(19, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(20, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(21, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(24, "When user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(25, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(26, "And the \"rows shown\" reading of statistics viewer should be 0", () => readingIs(page, "rows shown", el("statistics viewer"), 0));
    await session.step(27, "And the \"values of AGE\" reading of statistics viewer should be \"0\"", () => readingReads(page, "values of AGE", el("statistics viewer"), "0"));
    await session.step(28, "And the \"values of RACE\" reading of statistics viewer should be \"0\"", () => readingReads(page, "values of RACE", el("statistics viewer"), "0"));
    await session.step(29, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
    await session.step(30, "When user resets the filter", () => resetFilter(page));
    await session.step(31, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(32, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
  });
  test("Row Source Selected with nothing selected counts no rows without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(18, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(19, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(20, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(21, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(36, "Given user selects no rows", () => selectNoRows(page));
    await session.step(37, "When user sets \"rowSource\" property of statistics viewer to \"Selected\"", () => setProperty(page, "rowSource", el("statistics viewer"), "Selected"));
    await session.step(38, "Then the \"rows shown\" reading of statistics viewer should be 0", () => readingIs(page, "rows shown", el("statistics viewer"), 0));
    await session.step(39, "And the \"values of AGE\" reading of statistics viewer should be \"0\"", () => readingReads(page, "values of AGE", el("statistics viewer"), "0"));
    await session.step(40, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
    await session.step(41, "When user sets \"rowSource\" property of statistics viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("statistics viewer"), "Filtered"));
    await session.step(42, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
  test("A one-row table counts its one value", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(18, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(19, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(20, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(21, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(46, "Given user opens a table \"one-row\" with:", () => openTableOf(page, "one-row", [["name","x"],["a","42"]]), [["name","x"],["a","42"]]);
    await session.step(49, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(50, "Then the \"rows shown\" reading of statistics viewer should be 1", () => readingIs(page, "rows shown", el("statistics viewer"), 1));
    await session.step(51, "And the \"values of x\" reading of statistics viewer should be \"1\"", () => readingReads(page, "values of x", el("statistics viewer"), "1"));
    await session.step(52, "And the \"nulls of x\" reading of statistics viewer should be \"0\"", () => readingReads(page, "nulls of x", el("statistics viewer"), "0"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
  test("A calculated column added while the viewer is open gets a row", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(18, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(19, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(20, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(21, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(56, "When user adds a calculated column \"AGE2\" with formula \"${AGE} * 2\"", () => addCalculated(page, "AGE2", "${AGE} * 2"));
    await session.step(57, "Then the \"columns shown\" reading of statistics viewer should be 12", () => readingIs(page, "columns shown", el("statistics viewer"), 12));
    await session.step(58, "And the \"values of AGE2\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "values of AGE2", el("statistics viewer"), "1000"));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
  });
  test("A removed column loses its row without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(18, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(19, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(20, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(21, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(62, "Then statistics viewer should have a \"row RACE\" area", () => hasArea(page, el("statistics viewer"), "row RACE"));
    await session.step(63, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(64, "Then the \"columns shown\" reading of statistics viewer should be 10", () => readingIs(page, "columns shown", el("statistics viewer"), 10));
    await session.step(65, "And statistics viewer should not have a \"row RACE\" area", () => hasNoArea(page, el("statistics viewer"), "row RACE"));
    await session.step(66, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(67, "And no errors should have been logged", () => noErrors(page));
  });
  test("A renamed column keeps its row under the new name", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(18, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(19, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(20, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(21, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(70, "When user renames \"AGE\" column to \"AGE_YEARS\"", () => renameColumn(page, "AGE", "AGE_YEARS"));
    await session.step(71, "Then the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(72, "And the \"avg of AGE_YEARS\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE_YEARS", el("statistics viewer"), "45.68"));
    await session.step(73, "And statistics viewer should not have a \"row AGE\" area", () => hasNoArea(page, el("statistics viewer"), "row AGE"));
    await session.step(74, "And no errors should have been logged", () => noErrors(page));
  });
});
