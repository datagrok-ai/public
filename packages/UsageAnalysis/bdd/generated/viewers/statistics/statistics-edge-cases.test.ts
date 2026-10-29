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
import {columnType, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, filterBetween, filterPasses, removeColumn, renameColumn, resetFilter, selectNoRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, noErrors, pickFromAreaContextMenu, readingIs, readingReads, reportsNoError, resizeTo, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

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
    await session.step(29, "And the \"avg of AGE\" reading of statistics viewer should be \"\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), ""));
    await session.step(30, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
    await session.step(31, "When user resets the filter", () => resetFilter(page));
    await session.step(32, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(33, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
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
    await session.step(37, "Given user selects no rows", () => selectNoRows(page));
    await session.step(38, "When user sets \"rowSource\" property of statistics viewer to \"Selected\"", () => setProperty(page, "rowSource", el("statistics viewer"), "Selected"));
    await session.step(39, "Then the \"rows shown\" reading of statistics viewer should be 0", () => readingIs(page, "rows shown", el("statistics viewer"), 0));
    await session.step(40, "And the \"values of AGE\" reading of statistics viewer should be \"0\"", () => readingReads(page, "values of AGE", el("statistics viewer"), "0"));
    await session.step(41, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
    await session.step(42, "When user sets \"rowSource\" property of statistics viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("statistics viewer"), "Filtered"));
    await session.step(43, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
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
    await session.step(47, "Given user opens a table \"one-row\" with:", () => openTableOf(page, "one-row", [["name","x"],["a","42"]]), [["name","x"],["a","42"]]);
    await session.step(50, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(51, "Then the \"rows shown\" reading of statistics viewer should be 1", () => readingIs(page, "rows shown", el("statistics viewer"), 1));
    await session.step(52, "And the \"values of x\" reading of statistics viewer should be \"1\"", () => readingReads(page, "values of x", el("statistics viewer"), "1"));
    await session.step(53, "And the \"nulls of x\" reading of statistics viewer should be \"0\"", () => readingReads(page, "nulls of x", el("statistics viewer"), "0"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
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
    await session.step(57, "When user adds a calculated column \"AGE2\" with formula \"${AGE} * 2\"", () => addCalculated(page, "AGE2", "${AGE} * 2"));
    await session.step(58, "Then the \"columns shown\" reading of statistics viewer should be 12", () => readingIs(page, "columns shown", el("statistics viewer"), 12));
    await session.step(59, "And the \"values of AGE2\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "values of AGE2", el("statistics viewer"), "1000"));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
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
    await session.step(63, "Then statistics viewer should have a \"row RACE\" area", () => hasArea(page, el("statistics viewer"), "row RACE"));
    await session.step(64, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(65, "Then the \"columns shown\" reading of statistics viewer should be 10", () => readingIs(page, "columns shown", el("statistics viewer"), 10));
    await session.step(66, "And statistics viewer should not have a \"row RACE\" area", () => hasNoArea(page, el("statistics viewer"), "row RACE"));
    await session.step(67, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
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
    await session.step(71, "When user renames \"AGE\" column to \"AGE_YEARS\"", () => renameColumn(page, "AGE", "AGE_YEARS"));
    await session.step(72, "Then the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(73, "And the \"avg of AGE_YEARS\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE_YEARS", el("statistics viewer"), "45.68"));
    await session.step(74, "And statistics viewer should not have a \"row AGE\" area", () => hasNoArea(page, el("statistics viewer"), "row AGE"));
    await session.step(75, "And no errors should have been logged", () => noErrors(page));
  });
  test("A numeric column with no values has empty numeric statistics, not 0", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.stats-viewer", "@known-failure:85", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(18, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(19, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(20, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(21, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await knownFailure(async () => {
      await session.step(80, "When user adds a calculated column \"EMPTY\" with formula \"if(${AGE} > 1000, ${AGE}, null)\"", () => addCalculated(page, "EMPTY", "if(${AGE} > 1000, ${AGE}, null)"));
      await session.step(81, "Then \"EMPTY\" column should have type \"int\"", () => columnType(page, "EMPTY", "int"));
      await session.step(82, "And the \"values of EMPTY\" reading of statistics viewer should be \"0\"", () => readingReads(page, "values of EMPTY", el("statistics viewer"), "0"));
      await session.step(83, "And the \"nulls of EMPTY\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "nulls of EMPTY", el("statistics viewer"), "1000"));
      await session.step(84, "And the \"avg of EMPTY\" reading of statistics viewer should be \"\"", () => readingReads(page, "avg of EMPTY", el("statistics viewer"), ""));
      await session.step(85, "And the \"stdev of EMPTY\" reading of statistics viewer should be \"\"", () => readingReads(page, "stdev of EMPTY", el("statistics viewer"), ""));
      await session.step(86, "And the \"min of EMPTY\" reading of statistics viewer should be \"\"", () => readingReads(page, "min of EMPTY", el("statistics viewer"), ""));
      await session.step(87, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    }, 85);
  });
  test("A constant column has its value as min, max and avg and a stdev of 0", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(18, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(19, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(20, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(21, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(91, "When user adds a calculated column \"CONST\" with formula \"5\"", () => addCalculated(page, "CONST", "5"));
    await session.step(92, "Then the \"values of CONST\" reading of statistics viewer should be \"1000\"", () => readingReads(page, "values of CONST", el("statistics viewer"), "1000"));
    await session.step(93, "And the \"min of CONST\" reading of statistics viewer should be \"5.00\"", () => readingReads(page, "min of CONST", el("statistics viewer"), "5.00"));
    await session.step(94, "And the \"max of CONST\" reading of statistics viewer should be \"5.00\"", () => readingReads(page, "max of CONST", el("statistics viewer"), "5.00"));
    await session.step(95, "And the \"avg of CONST\" reading of statistics viewer should be \"5.00\"", () => readingReads(page, "avg of CONST", el("statistics viewer"), "5.00"));
    await session.step(96, "And the \"stdev of CONST\" reading of statistics viewer should be \"0.00\"", () => readingReads(page, "stdev of CONST", el("statistics viewer"), "0.00"));
    await session.step(97, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
    await session.step(98, "And no errors should have been logged", () => noErrors(page));
  });
  test("A Histograms column is dropped when its column is removed", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(18, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(19, "Then the \"rows shown\" reading of statistics viewer should be 1000", () => readingIs(page, "rows shown", el("statistics viewer"), 1000));
    await session.step(20, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(21, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(101, "When user picks \"Histograms > SEX\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Histograms > SEX", "row AGE", el("statistics viewer")));
    await session.step(102, "Then the \"histogram columns\" reading of statistics viewer should be 1", () => readingIs(page, "histogram columns", el("statistics viewer"), 1));
    await session.step(103, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
    await session.step(104, "Then the table should not have a column \"SEX\"", () => hasNoColumn(page, "SEX"));
    await session.step(105, "And the \"histogram columns\" reading of statistics viewer should be 0", () => readingIs(page, "histogram columns", el("statistics viewer"), 0));
    await session.step(106, "And the \"avg of AGE\" reading of statistics viewer should be \"45.68\"", () => readingReads(page, "avg of AGE", el("statistics viewer"), "45.68"));
    await session.step(107, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
    await session.step(108, "And no errors should have been logged", () => noErrors(page));
  });
});
