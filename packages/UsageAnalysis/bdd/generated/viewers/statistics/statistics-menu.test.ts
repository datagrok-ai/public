/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/statistics/statistics-menu.feature
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
import {shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {tableColumns, tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, closeContextMenu, hasArea, hasNoArea, menuDoesNotList, menuLists, noErrors, pickFromAreaContextMenu, readingIs, readingReads, reportsNoError, resizeTo, rightClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains, readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Statistics and Histograms submenus, and what they turn on", () => {
  const session = feature(test, "features/viewers/statistics/statistics-menu.feature", import.meta.url);
  test("The Statistics and Histograms submenus, and what they turn on", {tag: ["@journey", "@viewers", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(25, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(26, "Then the \"stats\" reading of statistics viewer should be \"values, nulls, unique, min, max, avg, med, stdev\"", () => readingReads(page, "stats", el("statistics viewer"), "values, nulls, unique, min, max, avg, med, stdev"));
    await session.step(27, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(28, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
    await run.scenario("The Statistics submenu states which aggregations are on", async () => {
      await session.step(31, "When user right-clicks on the \"row AGE\" area of statistics viewer", () => rightClickArea(page, "row AGE", el("statistics viewer")));
      await session.step(32, "Then the open menu should list \"Statistics > sum\"", () => menuLists(page, "Statistics > sum"));
      await session.step(33, "And \"min\" menu item should be selected", () => shouldBe(page, el("\"min\" menu item"), "selected"));
      await session.step(34, "And \"max\" menu item should be selected", () => shouldBe(page, el("\"max\" menu item"), "selected"));
      await session.step(35, "And \"avg\" menu item should be selected", () => shouldBe(page, el("\"avg\" menu item"), "selected"));
      await session.step(36, "And \"stdev\" menu item should be selected", () => shouldBe(page, el("\"stdev\" menu item"), "selected"));
      await session.step(37, "And \"sum\" menu item should not be selected", () => shouldNotBe(page, el("\"sum\" menu item"), "selected"));
      await session.step(38, "And \"geomean\" menu item should not be selected", () => shouldNotBe(page, el("\"geomean\" menu item"), "selected"));
      await session.step(39, "And \"variance\" menu item should not be selected", () => shouldNotBe(page, el("\"variance\" menu item"), "selected"));
      await session.step(40, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(41, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Statistics > sum adds the column, and picking it again takes it away", async () => {
      await session.step(44, "Then the \"stats\" reading of statistics viewer should not contain \"sum\"", () => readingNotContains(page, "stats", el("statistics viewer"), "sum"));
      await session.step(45, "And statistics viewer should not have a \"header sum\" area", () => hasNoArea(page, el("statistics viewer"), "header sum"));
      await session.step(46, "When user picks \"Statistics > sum\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Statistics > sum", "row AGE", el("statistics viewer")));
      await session.step(47, "Then the \"stats\" reading of statistics viewer should contain \"sum\"", () => readingContains(page, "stats", el("statistics viewer"), "sum"));
      await session.step(48, "And statistics viewer should have a \"header sum\" area", () => hasArea(page, el("statistics viewer"), "header sum"));
      await session.step(49, "And the \"sum of AGE\" reading of statistics viewer should be \"45677.00\"", () => readingReads(page, "sum of AGE", el("statistics viewer"), "45677.00"));
      await session.step(50, "And the \"sum of SEX\" reading of statistics viewer should be \"\"", () => readingReads(page, "sum of SEX", el("statistics viewer"), ""));
      await session.step(51, "When user right-clicks on the \"row AGE\" area of statistics viewer", () => rightClickArea(page, "row AGE", el("statistics viewer")));
      await session.step(52, "Then the open menu should list \"Statistics > sum\"", () => menuLists(page, "Statistics > sum"));
      await session.step(53, "And \"sum\" menu item should be selected", () => shouldBe(page, el("\"sum\" menu item"), "selected"));
      await session.step(54, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(55, "And user picks \"Statistics > sum\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Statistics > sum", "row AGE", el("statistics viewer")));
      await session.step(56, "Then the \"stats\" reading of statistics viewer should not contain \"sum\"", () => readingNotContains(page, "stats", el("statistics viewer"), "sum"));
      await session.step(57, "And statistics viewer should not have a \"header sum\" area", () => hasNoArea(page, el("statistics viewer"), "header sum"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An aggregation dropped from the list takes its column with it", async () => {
      await session.step(61, "Then statistics viewer should have a \"header med\" area", () => hasArea(page, el("statistics viewer"), "header med"));
      await session.step(62, "When user picks \"Statistics > med\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Statistics > med", "row AGE", el("statistics viewer")));
      await session.step(63, "Then the \"stats\" reading of statistics viewer should not contain \"med\"", () => readingNotContains(page, "stats", el("statistics viewer"), "med"));
      await session.step(64, "And statistics viewer should not have a \"header med\" area", () => hasNoArea(page, el("statistics viewer"), "header med"));
      await session.step(65, "And statistics viewer should have a \"header avg\" area", () => hasArea(page, el("statistics viewer"), "header avg"));
      await session.step(66, "When user picks \"Statistics > med\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Statistics > med", "row AGE", el("statistics viewer")));
      await session.step(67, "Then the \"stats\" reading of statistics viewer should contain \"med\"", () => readingContains(page, "stats", el("statistics viewer"), "med"));
      await session.step(68, "And the \"med of AGE\" reading of statistics viewer should be \"45.00\"", () => readingReads(page, "med of AGE", el("statistics viewer"), "45.00"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Histograms submenu offers the categorical columns with fewer than ten categories", async () => {
      await session.step(72, "When user right-clicks on the \"row AGE\" area of statistics viewer", () => rightClickArea(page, "row AGE", el("statistics viewer")));
      await session.step(73, "Then the open menu should list \"Histograms > SEX\"", () => menuLists(page, "Histograms > SEX"));
      await session.step(74, "And the open menu should list \"Histograms > RACE\"", () => menuLists(page, "Histograms > RACE"));
      await session.step(75, "And the open menu should list \"Histograms > DIS_POP\"", () => menuLists(page, "Histograms > DIS_POP"));
      await session.step(76, "And the open menu should list \"Histograms > SEVERITY\"", () => menuLists(page, "Histograms > SEVERITY"));
      await session.step(77, "And the open menu should not list \"Histograms > USUBJID\"", () => menuDoesNotList(page, "Histograms > USUBJID"));
      await session.step(78, "And the open menu should not list \"Histograms > AGE\"", () => menuDoesNotList(page, "Histograms > AGE"));
      await session.step(79, "And the open menu should not list \"Histograms > STARTED\"", () => menuDoesNotList(page, "Histograms > STARTED"));
      await session.step(80, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(81, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Histograms > SEX adds a histogram column, and picking it again removes it", async () => {
      await session.step(84, "Then the \"histogram columns\" reading of statistics viewer should be 0", () => readingIs(page, "histogram columns", el("statistics viewer"), 0));
      await session.step(85, "When user picks \"Histograms > SEX\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Histograms > SEX", "row AGE", el("statistics viewer")));
      await session.step(86, "Then the \"histogram columns\" reading of statistics viewer should be 1", () => readingIs(page, "histogram columns", el("statistics viewer"), 1));
      await session.step(87, "And the \"stats\" reading of statistics viewer should contain \"avg\"", () => readingContains(page, "stats", el("statistics viewer"), "avg"));
      await session.step(88, "And the \"stats\" reading of statistics viewer should contain \"med\"", () => readingContains(page, "stats", el("statistics viewer"), "med"));
      await session.step(89, "And the \"stats\" reading of statistics viewer should contain \"stdev\"", () => readingContains(page, "stats", el("statistics viewer"), "stdev"));
      await session.step(90, "And the \"stats\" reading of statistics viewer should not contain \"sum\"", () => readingNotContains(page, "stats", el("statistics viewer"), "sum"));
      await session.step(91, "When user picks \"Histograms > SEX\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Histograms > SEX", "row AGE", el("statistics viewer")));
      await session.step(92, "Then the \"histogram columns\" reading of statistics viewer should be 0", () => readingIs(page, "histogram columns", el("statistics viewer"), 0));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Open as table hands the statistics over as a table of their own", async () => {
      await session.step(96, "When user picks \"Open as table\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Open as table", "row AGE", el("statistics viewer")));
      await session.step(97, "Then table \"demog-1000 stats\" should be open", () => tableOpen(page, "demog-1000 stats"));
      await session.step(98, "And table \"demog-1000 stats\" should have 11 rows", () => tableRows(page, "demog-1000 stats", 11));
      await session.step(99, "And table \"demog-1000 stats\" should have columns \"__t, __name, values, nulls, unique, min, max, avg, stdev, med\"", () => tableColumns(page, "demog-1000 stats", "__t, __name, values, nulls, unique, min, max, avg, stdev, med"));
      await session.step(100, "When user switches to the \"demog-1000 stats\" table view", () => switchTableView(page, "demog-1000 stats"));
      await session.step(101, "Then the value of \"__name\" column in row 2 should be \"AGE\"", () => valueInRow(page, "__name", 2, "AGE"));
      await session.step(102, "And the value of \"max\" column in row 2 should be \"89\"", () => valueInRow(page, "max", 2, "89"));
      await session.step(103, "And the value of \"__name\" column in row 6 should be \"HEIGHT\"", () => valueInRow(page, "__name", 6, "HEIGHT"));
      await session.step(104, "And the value of \"nulls\" column in row 6 should be \"128\"", () => valueInRow(page, "nulls", 6, "128"));
      await session.step(105, "When user closes the current view", () => closeCurrentView(page));
      await session.step(106, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(107, "Then the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
