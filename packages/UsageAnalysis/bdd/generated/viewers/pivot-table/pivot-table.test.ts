/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pivot-table/pivot-table.feature
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, hasArea, hasNoArea, noErrors, propertyShouldBe, readingIs, readingReads, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pivot table — the frame around the aggregation", () => {
  const session = feature(test, "features/viewers/pivot-table/pivot-table.feature", import.meta.url);
  test("Pivot table — the frame around the aggregation", {tag: ["@journey", "@viewers", "@realizes:viewers.pivot-table"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(19, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
    await session.step(20, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
    await session.step(21, "And the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
    await run.scenario("A fresh viewer configures itself from the column types", async () => {
      await session.step(24, "Then the \"key columns\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "key columns", el("pivot table viewer"), "DIS_POP"));
      await session.step(25, "And the \"aggregations\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregations", el("pivot table viewer"), "avg(AGE)"));
      await session.step(26, "And the \"data\" reading of pivot table viewer should be \"demog-1000\"", () => readingReads(page, "data", el("pivot table viewer"), "demog-1000"));
      await session.step(27, "And the \"default aggregation\" reading of pivot table viewer should be \"avg\"", () => readingReads(page, "default aggregation", el("pivot table viewer"), "avg"));
      await session.step(28, "And the \"error\" reading of pivot table viewer should be \"\"", () => readingReads(page, "error", el("pivot table viewer"), ""));
      await session.step(29, "And pivot table viewer should have a \"group by chip DIS_POP\" area", () => hasArea(page, el("pivot table viewer"), "group by chip DIS_POP"));
      await session.step(30, "And pivot table viewer should have a \"aggregate chip avg(AGE)\" area", () => hasArea(page, el("pivot table viewer"), "aggregate chip avg(AGE)"));
      await session.step(31, "And pivot table viewer should have a \"pivot chip SEVERITY\" area", () => hasArea(page, el("pivot table viewer"), "pivot chip SEVERITY"));
      await session.step(32, "And \"Group By Column Names\" property of pivot table viewer should be \"DIS_POP\"", () => propertyShouldBe(page, "Group By Column Names", el("pivot table viewer"), "DIS_POP"));
      await session.step(33, "And \"Aggregate Column Names\" property of pivot table viewer should be \"AGE\"", () => propertyShouldBe(page, "Aggregate Column Names", el("pivot table viewer"), "AGE"));
      await session.step(34, "And \"Aggregate Agg Types\" property of pivot table viewer should be \"avg\"", () => propertyShouldBe(page, "Aggregate Agg Types", el("pivot table viewer"), "avg"));
      await session.step(35, "And \"Pivot Column Names\" property of pivot table viewer should be \"SEVERITY\"", () => propertyShouldBe(page, "Pivot Column Names", el("pivot table viewer"), "SEVERITY"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The counts area reports the shape of the aggregation, and the cells hold it", async () => {
      await session.step(39, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(40, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(41, "And the \"limited\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "limited", el("pivot table viewer"), "false"));
      await session.step(42, "And the \"rows shown\" reading of pivot table viewer should be 1000", () => readingIs(page, "rows shown", el("pivot table viewer"), 1000));
      await session.step(43, "And pivot table viewer should have a \"counts\" area", () => hasArea(page, el("pivot table viewer"), "counts"));
      await session.step(44, "And pivot table viewer should have a \"add to workspace\" area", () => hasArea(page, el("pivot table viewer"), "add to workspace"));
      await session.step(45, "And the \"text of grid cell 5 of DIS_POP\" reading of pivot table viewer should be \"RA\"", () => readingReads(page, "text of grid cell 5 of DIS_POP", el("pivot table viewer"), "RA"));
      await session.step(46, "And the \"text of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"52.30\"", () => readingReads(page, "text of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "52.30"));
      await session.step(47, "And the \"text of grid cell 2 of Critical avg(AGE)\" reading of pivot table viewer should be \"29.00\"", () => readingReads(page, "text of grid cell 2 of Critical avg(AGE)", el("pivot table viewer"), "29.00"));
      await session.step(48, "And the \"text of grid cell 1 of Critical avg(AGE)\" reading of pivot table viewer should be \"\"", () => readingReads(page, "text of grid cell 1 of Critical avg(AGE)", el("pivot table viewer"), ""));
      await session.step(49, "And the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"DIS_POP\" pivoted on \"SEVERITY\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "DIS_POP", "SEVERITY"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Header takes the tag rows and the counts away and gives them back", async () => {
      await session.step(53, "When user sets \"Show Header\" property of pivot table viewer to \"false\"", () => setProperty(page, "Show Header", el("pivot table viewer"), "false"));
      await session.step(54, "Then the \"header shown\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "header shown", el("pivot table viewer"), "false"));
      await session.step(55, "And pivot table viewer should not have a \"group by row\" area", () => hasNoArea(page, el("pivot table viewer"), "group by row"));
      await session.step(56, "And pivot table viewer should not have a \"aggregate row\" area", () => hasNoArea(page, el("pivot table viewer"), "aggregate row"));
      await session.step(57, "And pivot table viewer should not have a \"pivot row\" area", () => hasNoArea(page, el("pivot table viewer"), "pivot row"));
      await session.step(58, "And pivot table viewer should not have a \"counts\" area", () => hasNoArea(page, el("pivot table viewer"), "counts"));
      await session.step(59, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(60, "And the \"text of grid cell 5 of DIS_POP\" reading of pivot table viewer should be \"RA\"", () => readingReads(page, "text of grid cell 5 of DIS_POP", el("pivot table viewer"), "RA"));
      await session.step(61, "When user sets \"Show Header\" property of pivot table viewer to \"true\"", () => setProperty(page, "Show Header", el("pivot table viewer"), "true"));
      await session.step(62, "Then the \"header shown\" reading of pivot table viewer should be \"true\"", () => readingReads(page, "header shown", el("pivot table viewer"), "true"));
      await session.step(63, "And pivot table viewer should have a \"group by row\" area", () => hasArea(page, el("pivot table viewer"), "group by row"));
      await session.step(64, "And pivot table viewer should have a \"counts\" area", () => hasArea(page, el("pivot table viewer"), "counts"));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Command Bar takes the history and refresh icons away and gives them back", async () => {
      await session.step(68, "Then pivot table viewer should have a \"history\" area", () => hasArea(page, el("pivot table viewer"), "history"));
      await session.step(69, "And pivot table viewer should have a \"refresh\" area", () => hasArea(page, el("pivot table viewer"), "refresh"));
      await session.step(70, "When user sets \"Show Command Bar\" property of pivot table viewer to \"false\"", () => setProperty(page, "Show Command Bar", el("pivot table viewer"), "false"));
      await session.step(71, "Then the \"command bar shown\" reading of pivot table viewer should be \"false\"", () => readingReads(page, "command bar shown", el("pivot table viewer"), "false"));
      await session.step(72, "And pivot table viewer should not have a \"command bar\" area", () => hasNoArea(page, el("pivot table viewer"), "command bar"));
      await session.step(73, "And pivot table viewer should not have a \"history\" area", () => hasNoArea(page, el("pivot table viewer"), "history"));
      await session.step(74, "And pivot table viewer should not have a \"refresh\" area", () => hasNoArea(page, el("pivot table viewer"), "refresh"));
      await session.step(75, "And pivot table viewer should have a \"group by row\" area", () => hasArea(page, el("pivot table viewer"), "group by row"));
      await session.step(76, "When user sets \"Show Command Bar\" property of pivot table viewer to \"true\"", () => setProperty(page, "Show Command Bar", el("pivot table viewer"), "true"));
      await session.step(77, "Then the \"command bar shown\" reading of pivot table viewer should be \"true\"", () => readingReads(page, "command bar shown", el("pivot table viewer"), "true"));
      await session.step(78, "And pivot table viewer should have a \"history\" area", () => hasArea(page, el("pivot table viewer"), "history"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title bar shows the title and the description obeys its visibility mode", async () => {
      await session.step(82, "When user sets properties of pivot table viewer:", () => setProperties(page, el("pivot table viewer"), [["Show Title","true"],["Title","Cross tab"]]));
      await session.step(85, "Then title of pivot table viewer should have text \"Cross tab\"", () => shouldHaveText(page, el("title of pivot table viewer"), "Cross tab"));
      await session.step(86, "When user sets properties of pivot table viewer:", () => setProperties(page, el("pivot table viewer"), [["Description","Rows by disease"],["Description Visibility Mode","Always"]]));
      await session.step(89, "Then description of pivot table viewer should be visible", () => shouldBe(page, el("description of pivot table viewer"), "visible"));
      await session.step(90, "And description of pivot table viewer should have text \"Rows by disease\"", () => shouldHaveText(page, el("description of pivot table viewer"), "Rows by disease"));
      await session.step(91, "When user sets \"Description Visibility Mode\" property of pivot table viewer to \"Never\"", () => setProperty(page, "Description Visibility Mode", el("pivot table viewer"), "Never"));
      await session.step(92, "Then description of pivot table viewer should be hidden", () => shouldBe(page, el("description of pivot table viewer"), "hidden"));
      await session.step(93, "When user sets properties of pivot table viewer:", () => setProperties(page, el("pivot table viewer"), [["Title",""],["Description",""],["Description Visibility Mode","Auto"]]));
      await session.step(97, "Then the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(98, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the viewer and adding it again brings the same cross tab back", async () => {
      await session.step(101, "When user clicks on close icon of pivot table viewer", () => clickOn(page, el("close icon of pivot table viewer")));
      await session.step(102, "Then pivot table viewer should be absent", () => shouldBe(page, el("pivot table viewer"), "absent"));
      await session.step(103, "When user adds a pivot table viewer", () => addViewer(page, "pivot table"));
      await session.step(104, "Then pivot table viewer should be visible", () => shouldBe(page, el("pivot table viewer"), "visible"));
      await session.step(105, "And the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
      await session.step(106, "And the \"aggregate\" reading of pivot table viewer should be \"avg(AGE)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE)"));
      await session.step(107, "And the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
      await session.step(108, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(109, "And the \"text of grid cell 5 of None avg(AGE)\" reading of pivot table viewer should be \"52.30\"", () => readingReads(page, "text of grid cell 5 of None avg(AGE)", el("pivot table viewer"), "52.30"));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Select columns dialog of the Data row leaves the row as it was", async () => {
      await session.step(113, "When user clicks on the \"data chip demog-1000\" area of pivot table viewer", () => clickArea(page, "data chip demog-1000", el("pivot table viewer")));
      await session.step(114, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(115, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(116, "Then \"Select columns...\" dialog should be absent", () => shouldBe(page, el("\"Select columns...\" dialog"), "absent"));
      await session.step(117, "And the \"data\" reading of pivot table viewer should be \"demog-1000\"", () => readingReads(page, "data", el("pivot table viewer"), "demog-1000"));
      await session.step(118, "And the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
      await session.step(119, "And the \"aggregated rows\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 6));
      await session.step(120, "And the \"aggregated columns\" reading of pivot table viewer should be 6", () => readingIs(page, "aggregated columns", el("pivot table viewer"), 6));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
