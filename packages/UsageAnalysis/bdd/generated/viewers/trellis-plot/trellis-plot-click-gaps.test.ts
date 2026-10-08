/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-click-gaps.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, isExpanded, pressKeyIn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, allOfFiltered, clearSelection, filterPasses, filterPassesAll, noneOfFiltered, noneOfSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hasArea, hasNoArea, noErrors, propertyShouldBe, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInnerViewer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot On Click set in the context panel, and what a change keeps", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-click-gaps.feature", import.meta.url);
  test("Trellis plot On Click set in the context panel, and what a change keeps", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(23, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await session.step(24, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
    await run.scenario("Under Select, a change of split column keeps the selection", async () => {
      await session.step(27, "When user sets \"On Click\" property of trellis plot viewer to \"Select\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Select"));
      await session.step(28, "And user clicks on the \"cell M | Asian\" area of trellis plot viewer", () => clickArea(page, "cell M | Asian", el("trellis plot viewer")));
      await session.step(29, "Then 8 rows should be selected", () => selectedRowCount(page, 8));
      await session.step(30, "And no rows where \"SEX\" is \"F\" should be selected", () => noneOfSelected(page, "SEX", "F"));
      await session.step(31, "And no rows where \"RACE\" is \"Caucasian\" should be selected", () => noneOfSelected(page, "RACE", "Caucasian"));
      await session.step(32, "When user sets \"X Column Names\" property of trellis plot viewer to \"CONTROL\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "CONTROL"));
      await session.step(33, "Then trellis plot viewer should not have a \"cell M | Asian\" area", () => hasNoArea(page, el("trellis plot viewer"), "cell M | Asian"));
      await session.step(34, "And trellis plot viewer should have a \"cell false | Asian\" area", () => hasArea(page, el("trellis plot viewer"), "cell false | Asian"));
      await session.step(35, "And 8 rows should be selected", () => selectedRowCount(page, 8));
      await session.step(36, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX"));
      await session.step(37, "And user clears the row selection", () => clearSelection(page));
      await session.step(38, "And user sets \"On Click\" property of trellis plot viewer to \"None\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(39, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Under Filter, a change of the inner viewer keeps the filter, and Escape then drops it", async () => {
      await session.step(42, "When user sets \"On Click\" property of trellis plot viewer to \"Filter\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Filter"));
      await session.step(43, "And user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(44, "Then 480 rows should pass the filter", () => filterPasses(page, 480));
      await session.step(45, "When user picks \"Bar chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Bar chart", el("trellis plot viewer")));
      await session.step(46, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Bar chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Bar chart"));
      await session.step(47, "And the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(48, "And 480 rows should pass the filter", () => filterPasses(page, 480));
      await session.step(49, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(50, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(51, "When user picks \"Scatter plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Scatter plot", el("trellis plot viewer")));
      await session.step(52, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With a filter card on, Escape takes back only the trellis's part", async () => {
      await session.step(56, "When user adds a categorical filter on \"DIS_POP\" keeping \"RA\"", () => addCategoricalFilter(page, "DIS_POP", "RA"));
      await session.step(57, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(58, "When user clicks on the \"cell M | Asian\" area of trellis plot viewer", () => clickArea(page, "cell M | Asian", el("trellis plot viewer")));
      await session.step(59, "Then 1 row should pass the filter", () => filterPasses(page, 1));
      await session.step(60, "And no rows where \"SEX\" is \"F\" should pass the filter", () => noneOfFiltered(page, "SEX", "F"));
      await session.step(61, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(62, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(63, "And all rows where \"DIS_POP\" is \"RA\" should pass the filter", () => allOfFiltered(page, "DIS_POP", "RA"));
      await session.step(64, "When user hovers over \"DIS_POP\" filter card", () => hoverOver(page, el("\"DIS_POP\" filter card")));
      await session.step(65, "And user clicks on close of \"DIS_POP\" filter card", () => clickOn(page, el("close of \"DIS_POP\" filter card")));
      await session.step(66, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(67, "When user sets \"On Click\" property of trellis plot viewer to \"None\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(68, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Set in the context panel, On Click and Row Source correct each other", async () => {
      await session.step(71, "When user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
      await session.step(72, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(73, "When user selects \"Filtered\" in \"Row Source\" property in context panel", () => selectIn(page, "Filtered", el("\"Row Source\" property in context panel")));
      await session.step(74, "And \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(75, "When user selects \"Filter\" in \"On Click\" property in context panel", () => selectIn(page, "Filter", el("\"On Click\" property in context panel")));
      await session.step(76, "Then \"On Click\" property of trellis plot viewer should be \"Filter\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "Filter"));
      await session.step(77, "And \"Row Source\" property of trellis plot viewer should be \"All\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "All"));
      await session.step(78, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(79, "Then 480 rows should pass the filter", () => filterPasses(page, 480));
      await session.step(80, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(81, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(82, "When user selects \"Filtered\" in \"Row Source\" property in context panel", () => selectIn(page, "Filtered", el("\"Row Source\" property in context panel")));
      await session.step(83, "Then \"Row Source\" property of trellis plot viewer should be \"Filtered\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "Filtered"));
      await session.step(84, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(85, "When user sets \"Row Source\" property of trellis plot viewer to \"All\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "All"));
      await session.step(86, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
