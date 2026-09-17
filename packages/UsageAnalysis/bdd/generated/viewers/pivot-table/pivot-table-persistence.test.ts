/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pivot-table]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {pivotedAggregationMatches} from '../../../bindings/pivot-table.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, closeContextMenu, loadLayout, noErrors, pickFromAreaContextMenu, readingDoesNotRead, readingIs, readingReads, saveLayoutToServer, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — a configured cross tab that survives a round trip", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table-persistence.feature", import.meta.url);
  test("Pivot table — a configured cross tab that survives a round trip", {tag: ["@journey", "@viewers", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(16, "And user picks \"Aggregation > med\" from the context menu of the \"aggregate chip avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Aggregation > med", "aggregate chip avg(AGE)", el("pivot table viewer")));
    await session.step(17, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(18, "Then the \"aggregations\" reading of pivot table viewer should be \"med(AGE)\"", () => readingReads(page, "aggregations", el("pivot table viewer"), "med(AGE)"));
    await session.step(19, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
    await session.step(20, "And the aggregated values of pivot table viewer should match \"med(AGE)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "med(AGE)", "DIS_POP", "SEVERITY"));
    await run.scenario("A layout saved on the server brings the non-default aggregation back", async () => {
      await session.step(23, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(24, "And user clicks on close icon of pivot table viewer", () => clickOn(page, el("close icon of pivot table viewer")));
      await session.step(25, "Then pivot table viewer should be absent", () => shouldBe(page, el("pivot table viewer"), "absent"));
      await session.step(26, "When user loads the saved layout", () => loadLayout(page));
      await session.step(27, "Then pivot table viewer should be visible", () => shouldBe(page, el("pivot table viewer"), "visible"));
      await session.step(28, "And the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
      await session.step(29, "And the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
      await session.step(30, "And the \"aggregate\" reading of pivot table viewer should be \"med(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "med(AGE)"));
      await session.step(31, "And the \"aggregations\" reading of pivot table viewer should be \"med(AGE)\"", () => readingReads(page, "aggregations", el("pivot table viewer"), "med(AGE)"));
      await session.step(32, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(33, "And the aggregated values of pivot table viewer should match \"med(AGE)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "med(AGE)", "DIS_POP", "SEVERITY"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title and the inner grid's colour coding travel with the layout", async () => {
      await session.step(37, "When user sets \"Title\" property of pivot table viewer to \"Cross tab\"", () => setProperty(page, "Title", el("pivot table viewer"), "Cross tab"));
      await session.step(38, "And user picks \"Grid > Color Coding > Linear\" from the context menu of the \"grid header None med(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Grid > Color Coding > Linear", "grid header None med(AGE)", el("pivot table viewer")));
      await session.step(39, "Then the \"color coding of None med(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None med(AGE)", el("pivot table viewer"), "Linear"));
      await session.step(40, "And title of pivot table viewer should have text \"Cross tab\"", () => shouldHaveText(page, el("title of pivot table viewer"), "Cross tab"));
      await session.step(41, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(42, "And user clicks on close icon of pivot table viewer", () => clickOn(page, el("close icon of pivot table viewer")));
      await session.step(43, "Then pivot table viewer should be absent", () => shouldBe(page, el("pivot table viewer"), "absent"));
      await session.step(44, "When user loads the saved layout", () => loadLayout(page));
      await session.step(45, "Then pivot table viewer should be visible", () => shouldBe(page, el("pivot table viewer"), "visible"));
      await session.step(46, "And title of pivot table viewer should have text \"Cross tab\"", () => shouldHaveText(page, el("title of pivot table viewer"), "Cross tab"));
      await session.step(47, "And the \"color coding of None med(AGE)\" reading of pivot table viewer should be \"Linear\"", () => readingReads(page, "color coding of None med(AGE)", el("pivot table viewer"), "Linear"));
      await session.step(48, "And the \"color of grid cell 2 of None med(AGE)\" reading of pivot table viewer should not be \"#ffffff\"", () => readingDoesNotRead(page, "color of grid cell 2 of None med(AGE)", el("pivot table viewer"), "#ffffff"));
      await session.step(49, "And the \"aggregations\" reading of pivot table viewer should be \"med(AGE)\"", () => readingReads(page, "aggregations", el("pivot table viewer"), "med(AGE)"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project saved, closed and reopened brings the whole configuration back", async () => {
      await session.step(53, "When user saves the current view as project \"bdd pivot table round trip\"", () => saveAsProject(page, "bdd pivot table round trip"));
      await session.step(54, "And user closes all views", () => closeAllViews(page));
      await session.step(55, "And user opens the \"bdd pivot table round trip\" project", () => openProject(page, "bdd pivot table round trip"));
      await session.step(56, "Then pivot table viewer should be visible", () => shouldBe(page, el("pivot table viewer"), "visible"));
      await session.step(57, "And the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
      await session.step(58, "And the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
      await session.step(59, "And the \"aggregations\" reading of pivot table viewer should be \"med(AGE)\"", () => readingReads(page, "aggregations", el("pivot table viewer"), "med(AGE)"));
      await session.step(60, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(61, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(62, "And title of pivot table viewer should have text \"Cross tab\"", () => shouldHaveText(page, el("title of pivot table viewer"), "Cross tab"));
      await session.step(63, "And the aggregated values of pivot table viewer should match \"med(AGE)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "med(AGE)", "DIS_POP", "SEVERITY"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
