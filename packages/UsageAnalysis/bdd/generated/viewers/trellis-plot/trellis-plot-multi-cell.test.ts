/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-multi-cell.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
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
import {clickOn, hoverOver, pressKeyIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, clearSelection, filterPasses, filterPassesAll, noneOfFiltered, noneOfSelected, noneSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, clickAreaHolding, noErrors, propertyShouldBe, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains, readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot multi-cell selection and filter", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-multi-cell.feature", import.meta.url);
  test("Trellis plot multi-cell selection and filter", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(19, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await session.step(20, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
    await run.scenario("The inner viewers stay live under None and go inert under Select", async () => {
      await session.step(23, "Then the \"cells inert\" reading of trellis plot viewer should be \"false\"", () => readingReads(page, "cells inert", el("trellis plot viewer"), "false"));
      await session.step(24, "When user sets \"On Click\" property of trellis plot viewer to \"Select\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Select"));
      await session.step(25, "Then the \"cells inert\" reading of trellis plot viewer should be \"true\"", () => readingReads(page, "cells inert", el("trellis plot viewer"), "true"));
      await session.step(26, "And the \"selected cells\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "selected cells", el("trellis plot viewer"), ""));
      await session.step(27, "And no rows should be selected", () => noneSelected(page));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a cell's body selects the whole cell", async () => {
      await session.step(31, "When user clicks on the \"cell body F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(32, "Then 480 rows should be selected", () => selectedRowCount(page, 480));
      await session.step(33, "And no rows where \"SEX\" is \"M\" should be selected", () => noneOfSelected(page, "SEX", "M"));
      await session.step(34, "And the \"selected cells\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "selected cells", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on the only selected cell clears the selection", async () => {
      await session.step(38, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(39, "Then no rows should be selected", () => noneSelected(page));
      await session.step(40, "And the \"selected cells\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "selected cells", el("trellis plot viewer"), ""));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control adds a cell, and takes a selected one back out", async () => {
      await session.step(44, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(45, "And user clicks on the \"cell M | Caucasian\" area of trellis plot viewer holding Control", () => clickAreaHolding(page, "cell M | Caucasian", el("trellis plot viewer"), "Control"));
      await session.step(46, "Then 896 rows should be selected", () => selectedRowCount(page, 896));
      await session.step(47, "And the \"selected cells\" reading of trellis plot viewer should be \"F | Caucasian, M | Caucasian\"", () => readingReads(page, "selected cells", el("trellis plot viewer"), "F | Caucasian, M | Caucasian"));
      await session.step(48, "When user clicks on the \"cell M | Asian\" area of trellis plot viewer holding Control", () => clickAreaHolding(page, "cell M | Asian", el("trellis plot viewer"), "Control"));
      await session.step(49, "Then 904 rows should be selected", () => selectedRowCount(page, 904));
      await session.step(50, "And the \"selected cells\" reading of trellis plot viewer should contain \"M | Asian\"", () => readingContains(page, "selected cells", el("trellis plot viewer"), "M | Asian"));
      await session.step(51, "When user clicks on the \"cell M | Caucasian\" area of trellis plot viewer holding Control", () => clickAreaHolding(page, "cell M | Caucasian", el("trellis plot viewer"), "Control"));
      await session.step(52, "Then 488 rows should be selected", () => selectedRowCount(page, 488));
      await session.step(53, "And the \"selected cells\" reading of trellis plot viewer should not contain \"M | Caucasian\"", () => readingNotContains(page, "selected cells", el("trellis plot viewer"), "M | Caucasian"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shift only adds a cell, Control+Shift only takes one away", async () => {
      await session.step(57, "When user clicks on the \"cell M | Asian\" area of trellis plot viewer holding Shift", () => clickAreaHolding(page, "cell M | Asian", el("trellis plot viewer"), "Shift"));
      await session.step(58, "Then 488 rows should be selected", () => selectedRowCount(page, 488));
      await session.step(59, "When user clicks on the \"cell F | Black\" area of trellis plot viewer holding Shift", () => clickAreaHolding(page, "cell F | Black", el("trellis plot viewer"), "Shift"));
      await session.step(60, "Then 506 rows should be selected", () => selectedRowCount(page, 506));
      await session.step(61, "And the \"selected cells\" reading of trellis plot viewer should contain \"F | Black\"", () => readingContains(page, "selected cells", el("trellis plot viewer"), "F | Black"));
      await session.step(62, "When user clicks on the \"cell M | Black\" area of trellis plot viewer holding Control+Shift", () => clickAreaHolding(page, "cell M | Black", el("trellis plot viewer"), "Control+Shift"));
      await session.step(63, "Then 506 rows should be selected", () => selectedRowCount(page, 506));
      await session.step(64, "And the \"selected cells\" reading of trellis plot viewer should not contain \"M | Black\"", () => readingNotContains(page, "selected cells", el("trellis plot viewer"), "M | Black"));
      await session.step(65, "When user clicks on the \"cell F | Black\" area of trellis plot viewer holding Control+Shift", () => clickAreaHolding(page, "cell F | Black", el("trellis plot viewer"), "Control+Shift"));
      await session.step(66, "Then 488 rows should be selected", () => selectedRowCount(page, 488));
      await session.step(67, "And the \"selected cells\" reading of trellis plot viewer should not contain \"F | Black\"", () => readingNotContains(page, "selected cells", el("trellis plot viewer"), "F | Black"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A selection changed elsewhere drops the outlines", async () => {
      await session.step(71, "When user clears the row selection", () => clearSelection(page));
      await session.step(72, "Then no rows should be selected", () => noneSelected(page));
      await session.step(73, "And the \"selected cells\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "selected cells", el("trellis plot viewer"), ""));
      await session.step(74, "When user clicks on the \"cell M | Asian\" area of trellis plot viewer holding Control", () => clickAreaHolding(page, "cell M | Asian", el("trellis plot viewer"), "Control"));
      await session.step(75, "Then 8 rows should be selected", () => selectedRowCount(page, 8));
      await session.step(76, "And the \"selected cells\" reading of trellis plot viewer should be \"M | Asian\"", () => readingReads(page, "selected cells", el("trellis plot viewer"), "M | Asian"));
      await session.step(77, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(78, "Then no rows should be selected", () => noneSelected(page));
      await session.step(79, "And the \"selected cells\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "selected cells", el("trellis plot viewer"), ""));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Under Filter the cells filter to their union, composed with a filter card", async () => {
      await session.step(83, "When user sets \"On Click\" property of trellis plot viewer to \"Filter\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Filter"));
      await session.step(84, "Then the \"cells inert\" reading of trellis plot viewer should be \"true\"", () => readingReads(page, "cells inert", el("trellis plot viewer"), "true"));
      await session.step(85, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(86, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(87, "And user clicks on the \"cell M | Caucasian\" area of trellis plot viewer holding Control", () => clickAreaHolding(page, "cell M | Caucasian", el("trellis plot viewer"), "Control"));
      await session.step(88, "And user clicks on the \"cell M | Asian\" area of trellis plot viewer holding Control", () => clickAreaHolding(page, "cell M | Asian", el("trellis plot viewer"), "Control"));
      await session.step(89, "Then 904 rows should pass the filter", () => filterPasses(page, 904));
      await session.step(90, "And no rows where \"RACE\" is \"Black\" should pass the filter", () => noneOfFiltered(page, "RACE", "Black"));
      await session.step(91, "When user adds a categorical filter on \"DIS_POP\" keeping \"RA\"", () => addCategoricalFilter(page, "DIS_POP", "RA"));
      await session.step(92, "Then 365 rows should pass the filter", () => filterPasses(page, 365));
      await session.step(93, "When user clicks on the \"cell M | Caucasian\" area of trellis plot viewer holding Control", () => clickAreaHolding(page, "cell M | Caucasian", el("trellis plot viewer"), "Control"));
      await session.step(94, "Then 272 rows should pass the filter", () => filterPasses(page, 272));
      await session.step(95, "And the \"selected cells\" reading of trellis plot viewer should not contain \"M | Caucasian\"", () => readingNotContains(page, "selected cells", el("trellis plot viewer"), "M | Caucasian"));
      await session.step(96, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(97, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(98, "And the \"selected cells\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "selected cells", el("trellis plot viewer"), ""));
      await session.step(99, "When user hovers over \"DIS_POP\" filter card", () => hoverOver(page, el("\"DIS_POP\" filter card")));
      await session.step(100, "And user clicks on close of \"DIS_POP\" filter card", () => clickOn(page, el("close of \"DIS_POP\" filter card")));
      await session.step(101, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on the only filtered cell gives every row back", async () => {
      await session.step(105, "When user clicks on the \"cell M | Asian\" area of trellis plot viewer", () => clickArea(page, "cell M | Asian", el("trellis plot viewer")));
      await session.step(106, "Then 8 rows should pass the filter", () => filterPasses(page, 8));
      await session.step(107, "When user clicks on the \"cell M | Asian\" area of trellis plot viewer", () => clickArea(page, "cell M | Asian", el("trellis plot viewer")));
      await session.step(108, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(109, "And the \"selected cells\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "selected cells", el("trellis plot viewer"), ""));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Switching On Click drops the trellis filter, and None makes the inner viewers live again", async () => {
      await session.step(113, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(114, "Then 480 rows should pass the filter", () => filterPasses(page, 480));
      await session.step(115, "When user sets \"On Click\" property of trellis plot viewer to \"Select\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Select"));
      await session.step(116, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(117, "And no rows should be selected", () => noneSelected(page));
      await session.step(118, "And the \"selected cells\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "selected cells", el("trellis plot viewer"), ""));
      await session.step(119, "When user sets \"On Click\" property of trellis plot viewer to \"None\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(120, "Then the \"cells inert\" reading of trellis plot viewer should be \"false\"", () => readingReads(page, "cells inert", el("trellis plot viewer"), "false"));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
