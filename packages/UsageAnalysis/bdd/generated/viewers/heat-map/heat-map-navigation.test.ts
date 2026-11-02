/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map-navigation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map]
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
import {makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {hasCurrentRow, noneSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, clickAreaHolding, doubleClickArea, dragZoomOverArea, hasArea, hasNoArea, noErrors, painted, readingBetween, readingDiffers, readingHigher, readingIs, readingLower, readingReads, repainted, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map navigation and grid mode", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map-navigation.feature", import.meta.url);
  test("Heat map navigation and grid mode", {tag: ["@journey", "@viewers", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(19, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(20, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
    await session.step(21, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await run.scenario("A click on a column band makes a row of that column current", async () => {
      await session.step(24, "Given user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(25, "When user clicks on the \"column AGE\" area of heat map viewer", () => clickArea(page, "column AGE", el("heat map viewer")));
      await session.step(26, "Then the \"current column\" reading of heat map viewer should be \"AGE\"", () => readingReads(page, "current column", el("heat map viewer"), "AGE"));
      await session.step(27, "And the \"current row\" reading of heat map viewer should differ from before", () => readingDiffers(page, "current row", el("heat map viewer")));
      await session.step(28, "And the table should have a current row", () => hasCurrentRow(page));
      await session.step(29, "When user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(30, "Then the \"current row\" reading of heat map viewer should be 1", () => readingIs(page, "current row", el("heat map viewer"), 1));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An Alt-drag zooms both sliders in and a double-click on each resets its own axis", async () => {
      await session.step(34, "Then the \"x scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "x scroll span", el("heat map viewer"), 1));
      await session.step(35, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
      await session.step(36, "When user drags a zoom box over the \"column AGE\" area of heat map viewer", () => dragZoomOverArea(page, "column AGE", el("heat map viewer")));
      await session.step(37, "Then the \"x scroll span\" reading of heat map viewer should be lower than before", () => readingLower(page, "x scroll span", el("heat map viewer")));
      await session.step(38, "And the \"y scroll span\" reading of heat map viewer should be lower than before", () => readingLower(page, "y scroll span", el("heat map viewer")));
      await session.step(39, "And the \"row height\" reading of heat map viewer should be higher than before", () => readingHigher(page, "row height", el("heat map viewer")));
      await session.step(40, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(41, "And no rows should be selected", () => noneSelected(page));
      await session.step(42, "When user double-clicks on the \"y scroll slider\" area of heat map viewer", () => doubleClickArea(page, "y scroll slider", el("heat map viewer")));
      await session.step(43, "Then the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
      await session.step(44, "And the \"x scroll span\" reading of heat map viewer should be between 0 and 0.99", () => readingBetween(page, "x scroll span", el("heat map viewer"), 0, 0.99));
      await session.step(45, "When user double-clicks on the \"x scroll slider\" area of heat map viewer", () => doubleClickArea(page, "x scroll slider", el("heat map viewer")));
      await session.step(46, "Then the \"x scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "x scroll span", el("heat map viewer"), 1));
      await session.step(47, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
      await session.step(49, "When user clicks on the \"column AGE\" area of heat map viewer holding Alt", () => clickAreaHolding(page, "column AGE", el("heat map viewer"), "Alt"));
      await session.step(50, "Then the \"y scroll span\" reading of heat map viewer should be 0", () => readingIs(page, "y scroll span", el("heat map viewer"), 0));
      await session.step(51, "And the \"row height\" reading of heat map viewer should be higher than before", () => readingHigher(page, "row height", el("heat map viewer")));
      await session.step(52, "And no rows should be selected", () => noneSelected(page));
      await session.step(53, "When user double-clicks on the \"y scroll slider\" area of heat map viewer", () => doubleClickArea(page, "y scroll slider", el("heat map viewer")));
      await session.step(54, "And user double-clicks on the \"x scroll slider\" area of heat map viewer", () => doubleClickArea(page, "x scroll slider", el("heat map viewer")));
      await session.step(55, "Then the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
      await session.step(56, "And the \"x scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "x scroll span", el("heat map viewer"), 1));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Is Heatmap off draws the same table as a spreadsheet", async () => {
      await session.step(60, "Then heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
      await session.step(61, "And heat map viewer should not have a \"cell 1 of AGE\" area", () => hasNoArea(page, el("heat map viewer"), "cell 1 of AGE"));
      await session.step(62, "When user sets \"isHeatmap\" property of heat map viewer to \"false\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "false"));
      await session.step(63, "Then the \"is heatmap\" reading of heat map viewer should be \"false\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "false"));
      await session.step(64, "And the \"row height\" reading of heat map viewer should be between 20 and 40", () => readingBetween(page, "row height", el("heat map viewer"), 20, 40));
      await session.step(65, "And the \"y scroll span\" reading of heat map viewer should be lower than before", () => readingLower(page, "y scroll span", el("heat map viewer")));
      await session.step(66, "And heat map viewer should have a \"cell 1 of AGE\" area", () => hasArea(page, el("heat map viewer"), "cell 1 of AGE"));
      await session.step(67, "And the \"text of cell 1 of AGE\" reading of heat map viewer should be \"26\"", () => readingReads(page, "text of cell 1 of AGE", el("heat map viewer"), "26"));
      await session.step(68, "And heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
      await session.step(69, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(70, "When user sets \"isHeatmap\" property of heat map viewer to \"true\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "true"));
      await session.step(71, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Is Heatmap on again brings the whole table back on screen", async () => {
      await session.step(80, "Given user sets \"isHeatmap\" property of heat map viewer to \"false\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "false"));
      await session.step(81, "Then the \"row height\" reading of heat map viewer should be between 20 and 40", () => readingBetween(page, "row height", el("heat map viewer"), 20, 40));
      await session.step(82, "When user sets \"isHeatmap\" property of heat map viewer to \"true\"", () => setProperty(page, "isHeatmap", el("heat map viewer"), "true"));
      await session.step(83, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
      await session.step(84, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
      await session.step(85, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
      await session.step(86, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    });
    run.finish();
  });
});
