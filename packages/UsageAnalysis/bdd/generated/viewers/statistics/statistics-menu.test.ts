/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/statistics/statistics-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.stats-viewer]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, closeContextMenu, hasArea, hasNoArea, menuDoesNotList, menuLists, noErrors, pickFromAreaContextMenu, readingIs, readingReads, reportsNoError, resizeTo, rightClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains, readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Statistics and Histograms submenus, and what they turn on", () => {
  const session = feature(test, "features/viewers/statistics/statistics-menu.feature", import.meta.url);
  test("The Statistics and Histograms submenus, and what they turn on", {tag: ["@journey", "@viewers", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(19, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(20, "Then the \"stats\" reading of statistics viewer should be \"values, nulls, unique, min, max, avg, med, stdev\"", () => readingReads(page, "stats", el("statistics viewer"), "values, nulls, unique, min, max, avg, med, stdev"));
    await session.step(21, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(22, "And statistics viewer should report no error", () => reportsNoError(page, el("statistics viewer")));
    await run.scenario("The Statistics submenu states which aggregations are on", async () => {
      await session.step(25, "When user right-clicks on the \"row AGE\" area of statistics viewer", () => rightClickArea(page, "row AGE", el("statistics viewer")));
      await session.step(26, "Then the open menu should list \"Statistics > sum\"", () => menuLists(page, "Statistics > sum"));
      await session.step(27, "And \"min\" menu item should be selected", () => shouldBe(page, el("\"min\" menu item"), "selected"));
      await session.step(28, "And \"max\" menu item should be selected", () => shouldBe(page, el("\"max\" menu item"), "selected"));
      await session.step(29, "And \"avg\" menu item should be selected", () => shouldBe(page, el("\"avg\" menu item"), "selected"));
      await session.step(30, "And \"stdev\" menu item should be selected", () => shouldBe(page, el("\"stdev\" menu item"), "selected"));
      await session.step(31, "And \"sum\" menu item should not be selected", () => shouldNotBe(page, el("\"sum\" menu item"), "selected"));
      await session.step(32, "And \"geomean\" menu item should not be selected", () => shouldNotBe(page, el("\"geomean\" menu item"), "selected"));
      await session.step(33, "And \"variance\" menu item should not be selected", () => shouldNotBe(page, el("\"variance\" menu item"), "selected"));
      await session.step(34, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(35, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Statistics > sum adds the column, and picking it again takes it away", async () => {
      await session.step(38, "Then the \"stats\" reading of statistics viewer should not contain \"sum\"", () => readingNotContains(page, "stats", el("statistics viewer"), "sum"));
      await session.step(39, "And statistics viewer should not have a \"header sum\" area", () => hasNoArea(page, el("statistics viewer"), "header sum"));
      await session.step(40, "When user picks \"Statistics > sum\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Statistics > sum", "row AGE", el("statistics viewer")));
      await session.step(41, "Then the \"stats\" reading of statistics viewer should contain \"sum\"", () => readingContains(page, "stats", el("statistics viewer"), "sum"));
      await session.step(42, "And statistics viewer should have a \"header sum\" area", () => hasArea(page, el("statistics viewer"), "header sum"));
      await session.step(43, "And the \"sum of AGE\" reading of statistics viewer should be \"45677.00\"", () => readingReads(page, "sum of AGE", el("statistics viewer"), "45677.00"));
      await session.step(44, "And the \"sum of SEX\" reading of statistics viewer should be \"\"", () => readingReads(page, "sum of SEX", el("statistics viewer"), ""));
      await session.step(45, "When user right-clicks on the \"row AGE\" area of statistics viewer", () => rightClickArea(page, "row AGE", el("statistics viewer")));
      await session.step(46, "Then the open menu should list \"Statistics > sum\"", () => menuLists(page, "Statistics > sum"));
      await session.step(47, "And \"sum\" menu item should be selected", () => shouldBe(page, el("\"sum\" menu item"), "selected"));
      await session.step(48, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(49, "And user picks \"Statistics > sum\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Statistics > sum", "row AGE", el("statistics viewer")));
      await session.step(50, "Then the \"stats\" reading of statistics viewer should not contain \"sum\"", () => readingNotContains(page, "stats", el("statistics viewer"), "sum"));
      await session.step(51, "And statistics viewer should not have a \"header sum\" area", () => hasNoArea(page, el("statistics viewer"), "header sum"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An aggregation dropped from the list takes its column with it", async () => {
      await session.step(55, "Then statistics viewer should have a \"header med\" area", () => hasArea(page, el("statistics viewer"), "header med"));
      await session.step(56, "When user picks \"Statistics > med\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Statistics > med", "row AGE", el("statistics viewer")));
      await session.step(57, "Then the \"stats\" reading of statistics viewer should not contain \"med\"", () => readingNotContains(page, "stats", el("statistics viewer"), "med"));
      await session.step(58, "And statistics viewer should not have a \"header med\" area", () => hasNoArea(page, el("statistics viewer"), "header med"));
      await session.step(59, "And statistics viewer should have a \"header avg\" area", () => hasArea(page, el("statistics viewer"), "header avg"));
      await session.step(60, "When user picks \"Statistics > med\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Statistics > med", "row AGE", el("statistics viewer")));
      await session.step(61, "Then the \"stats\" reading of statistics viewer should contain \"med\"", () => readingContains(page, "stats", el("statistics viewer"), "med"));
      await session.step(62, "And the \"med of AGE\" reading of statistics viewer should be \"45.00\"", () => readingReads(page, "med of AGE", el("statistics viewer"), "45.00"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Histograms submenu offers the categorical columns with fewer than ten categories", async () => {
      await session.step(66, "When user right-clicks on the \"row AGE\" area of statistics viewer", () => rightClickArea(page, "row AGE", el("statistics viewer")));
      await session.step(67, "Then the open menu should list \"Histograms > SEX\"", () => menuLists(page, "Histograms > SEX"));
      await session.step(68, "And the open menu should list \"Histograms > RACE\"", () => menuLists(page, "Histograms > RACE"));
      await session.step(69, "And the open menu should list \"Histograms > DIS_POP\"", () => menuLists(page, "Histograms > DIS_POP"));
      await session.step(70, "And the open menu should list \"Histograms > SEVERITY\"", () => menuLists(page, "Histograms > SEVERITY"));
      await session.step(71, "And the open menu should not list \"Histograms > USUBJID\"", () => menuDoesNotList(page, "Histograms > USUBJID"));
      await session.step(72, "And the open menu should not list \"Histograms > AGE\"", () => menuDoesNotList(page, "Histograms > AGE"));
      await session.step(73, "And the open menu should not list \"Histograms > STARTED\"", () => menuDoesNotList(page, "Histograms > STARTED"));
      await session.step(74, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(75, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Histograms > SEX adds a histogram column, and picking it again removes it", async () => {
      await session.step(78, "Then the \"histogram columns\" reading of statistics viewer should be 0", () => readingIs(page, "histogram columns", el("statistics viewer"), 0));
      await session.step(79, "When user picks \"Histograms > SEX\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Histograms > SEX", "row AGE", el("statistics viewer")));
      await session.step(80, "Then the \"histogram columns\" reading of statistics viewer should be 1", () => readingIs(page, "histogram columns", el("statistics viewer"), 1));
      await session.step(81, "And the \"stats\" reading of statistics viewer should contain \"avg\"", () => readingContains(page, "stats", el("statistics viewer"), "avg"));
      await session.step(82, "And the \"stats\" reading of statistics viewer should contain \"med\"", () => readingContains(page, "stats", el("statistics viewer"), "med"));
      await session.step(83, "And the \"stats\" reading of statistics viewer should contain \"stdev\"", () => readingContains(page, "stats", el("statistics viewer"), "stdev"));
      await session.step(84, "And the \"stats\" reading of statistics viewer should not contain \"sum\"", () => readingNotContains(page, "stats", el("statistics viewer"), "sum"));
      await session.step(85, "When user picks \"Histograms > SEX\" from the context menu of the \"row AGE\" area of statistics viewer", () => pickFromAreaContextMenu(page, "Histograms > SEX", "row AGE", el("statistics viewer")));
      await session.step(86, "Then the \"histogram columns\" reading of statistics viewer should be 0", () => readingIs(page, "histogram columns", el("statistics viewer"), 0));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
