/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pc-plot, core#967]
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
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addRangeFilter, filterPasses, filterPassesAll, filterPassesFewer, openFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, eventFired, hasArea, lessInk, listenFor, moreInk, noErrors, pickFromContextMenu, propertyShouldBe, readingAsRemembered, readingBetween, readingHigher, readingIs, readingLower, readingNotAsRemembered, readingReads, readingsEqual, rememberReading, setProperty, showsFewerRows, showsRows, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragRangeHandle} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot in-chart range filter", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-filter.feature", import.meta.url);
  test("PC plot in-chart range filter", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot", "@realizes:core#967"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]), [["Column Names","AGE, HEIGHT, WEIGHT"]]);
    await session.step(22, "Then pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(23, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(24, "And the \"filtering\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("pc plot viewer"), "false"));
    await session.step(25, "And the \"range max of \\\"AGE\\\"\" reading of pc plot viewer should be 89", () => readingIs(page, "range max of \"AGE\"", el("pc plot viewer"), 89));
    await run.scenario("Dragging the AGE max handle filters the table, Reset View gives it back", async () => {
      await session.step(28, "Given user listens for \"d4-pc-plot-reset-view\" event on pc plot viewer", () => listenFor(page, "d4-pc-plot-reset-view", el("pc plot viewer")));
      await session.step(29, "When user hovers over pc plot viewer", () => hoverOver(page, el("pc plot viewer")));
      await session.step(30, "Then pc plot viewer should have a \"range max handle \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "range max handle \"AGE\""));
      await session.step(31, "When user drags the max handle of the \"AGE\" range slider of pc plot viewer by 120 pixels", () => dragRangeHandle(page, "max", "AGE", el("pc plot viewer"), 120));
      await session.step(32, "Then the \"filtering\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "filtering", el("pc plot viewer"), "true"));
      await session.step(33, "And the \"range max of \\\"AGE\\\"\" reading of pc plot viewer should be lower than before", () => readingLower(page, "range max of \"AGE\"", el("pc plot viewer")));
      await session.step(34, "And fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(35, "And pc plot viewer should show fewer rows than before", () => showsFewerRows(page, el("pc plot viewer")));
      await session.step(36, "And the \"lines drawn\" reading of pc plot viewer should be lower than before", () => readingLower(page, "lines drawn", el("pc plot viewer")));
      await session.step(37, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(38, "Then \"d4-pc-plot-reset-view\" event should have fired on pc plot viewer", () => eventFired(page, "d4-pc-plot-reset-view", el("pc plot viewer")));
      await session.step(39, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(40, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(41, "And the \"filtering\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("pc plot viewer"), "false"));
      await session.step(42, "And the \"range max of \\\"AGE\\\"\" reading of pc plot viewer should be 89", () => readingIs(page, "range max of \"AGE\"", el("pc plot viewer"), 89));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Filtered Out Lines draws the rows the slider dropped", async () => {
      await session.step(46, "When user drags the max handle of the \"AGE\" range slider of pc plot viewer by 120 pixels", () => dragRangeHandle(page, "max", "AGE", el("pc plot viewer"), 120));
      await session.step(47, "Then pc plot viewer should show fewer rows than before", () => showsFewerRows(page, el("pc plot viewer")));
      await session.step(48, "And the \"filtered out lines drawn\" reading of pc plot viewer should be 0", () => readingIs(page, "filtered out lines drawn", el("pc plot viewer"), 0));
      await session.step(49, "When user picks \"Filter > Show Filtered Out Lines\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Filter > Show Filtered Out Lines", el("pc plot viewer")));
      await session.step(50, "Then \"Show Filtered Out Lines\" property of pc plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Filtered Out Lines", el("pc plot viewer"), "true"));
      await session.step(51, "And the \"filtered out lines drawn\" reading of pc plot viewer should be higher than before", () => readingHigher(page, "filtered out lines drawn", el("pc plot viewer")));
      await session.step(52, "And pc plot viewer should have more ink than before", () => moreInk(page, el("pc plot viewer")));
      await session.step(53, "When user picks \"Filter > Show Filtered Out Lines\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Filter > Show Filtered Out Lines", el("pc plot viewer")));
      await session.step(54, "Then \"Show Filtered Out Lines\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Filtered Out Lines", el("pc plot viewer"), "false"));
      await session.step(55, "And the \"filtered out lines drawn\" reading of pc plot viewer should be 0", () => readingIs(page, "filtered out lines drawn", el("pc plot viewer"), 0));
      await session.step(56, "And pc plot viewer should have less ink than before", () => lessInk(page, el("pc plot viewer")));
      await session.step(57, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(58, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter card and the in-chart slider compose with AND", async () => {
      await session.step(62, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(63, "And user adds a range filter on \"AGE\" from 30 to 50", () => addRangeFilter(page, "AGE", 30, 50));
      await session.step(64, "Then 494 rows should pass the filter", () => filterPasses(page, 494));
      await session.step(65, "And pc plot viewer should show 494 rows", () => showsRows(page, el("pc plot viewer"), 494));
      await session.step(66, "When user drags the max handle of the \"HEIGHT\" range slider of pc plot viewer by 120 pixels", () => dragRangeHandle(page, "max", "HEIGHT", el("pc plot viewer"), 120));
      await session.step(67, "Then fewer than 494 rows should pass the filter", () => filterPassesFewer(page, 494));
      await session.step(68, "And pc plot viewer should show fewer rows than before", () => showsFewerRows(page, el("pc plot viewer")));
      await session.step(69, "And the \"filtering\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "filtering", el("pc plot viewer"), "true"));
      await session.step(70, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(71, "Then 494 rows should pass the filter", () => filterPasses(page, 494));
      await session.step(72, "And pc plot viewer should show 494 rows", () => showsRows(page, el("pc plot viewer"), 494));
      await session.step(73, "And the \"filtering\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("pc plot viewer"), "false"));
      await session.step(74, "When user drags the max handle of the \"HEIGHT\" range slider of pc plot viewer by 120 pixels", () => dragRangeHandle(page, "max", "HEIGHT", el("pc plot viewer"), 120));
      await session.step(75, "Then fewer than 494 rows should pass the filter", () => filterPassesFewer(page, 494));
      await session.step(76, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(77, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(78, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(79, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(80, "And the \"filtering\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("pc plot viewer"), "false"));
      await session.step(81, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(82, "Then filter panel should be absent", () => shouldBe(page, el("filter panel"), "absent"));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second range filter after a DateTime colour split still filters (GROK-18489)", async () => {
      await session.step(86, "When user sets \"Color\" property of pc plot viewer to \"STARTED\"", () => setProperty(page, "Color", el("pc plot viewer"), "STARTED"));
      await session.step(87, "And user drags the max handle of the \"AGE\" range slider of pc plot viewer by 120 pixels", () => dragRangeHandle(page, "max", "AGE", el("pc plot viewer"), 120));
      await session.step(88, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(89, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(90, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(91, "When user drags the min handle of the \"AGE\" range slider of pc plot viewer by -120 pixels", () => dragRangeHandle(page, "min", "AGE", el("pc plot viewer"), -120));
      await session.step(92, "Then the \"range min of \\\"AGE\\\"\" reading of pc plot viewer should be higher than before", () => readingHigher(page, "range min of \"AGE\"", el("pc plot viewer")));
      await session.step(93, "And fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(94, "When user remembers the \"rows shown\" reading of pc plot viewer", () => rememberReading(page, "rows shown", el("pc plot viewer")));
      await session.step(95, "And user drags the max handle of the \"AGE\" range slider of pc plot viewer by 100 pixels", () => dragRangeHandle(page, "max", "AGE", el("pc plot viewer"), 100));
      await session.step(96, "Then the \"rows shown\" reading of pc plot viewer should not be as remembered", () => readingNotAsRemembered(page, "rows shown", el("pc plot viewer")));
      await session.step(97, "And the \"range max of \\\"AGE\\\"\" reading of pc plot viewer should be lower than before", () => readingLower(page, "range max of \"AGE\"", el("pc plot viewer")));
      await session.step(98, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(99, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(100, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(101, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Changing a histogram's column leaves the plot's filter alone (github-972)", async () => {
      await session.step(104, "When user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"]]), [["Value","AGE"]]);
      await session.step(106, "And user drags the max handle of the \"AGE\" range slider of pc plot viewer by 120 pixels", () => dragRangeHandle(page, "max", "AGE", el("pc plot viewer"), 120));
      await session.step(107, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(108, "When user remembers the \"rows shown\" reading of pc plot viewer", () => rememberReading(page, "rows shown", el("pc plot viewer")));
      await session.step(109, "And user sets \"Value\" property of histogram viewer to \"HEIGHT\"", () => setProperty(page, "Value", el("histogram viewer"), "HEIGHT"));
      await session.step(110, "Then \"Value\" property of histogram viewer should be \"HEIGHT\"", () => propertyShouldBe(page, "Value", el("histogram viewer"), "HEIGHT"));
      await session.step(111, "And the \"rows shown\" reading of pc plot viewer should be as remembered", () => readingAsRemembered(page, "rows shown", el("pc plot viewer")));
      await session.step(112, "And the \"filtering\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "filtering", el("pc plot viewer"), "true"));
      await session.step(113, "When user clicks on close icon of histogram viewer", () => clickOn(page, el("close icon of histogram viewer")));
      await session.step(114, "Then histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
      await session.step(115, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(116, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(117, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The slider of a second plot filters the table both plots draw", async () => {
      await session.step(120, "Given user adds a pc plot viewer", () => addViewer(page, "pc plot"));
      await session.step(121, "When user sets \"Column Names\" property of second pc plot viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "Column Names", el("second pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(122, "Then the open tableview should have 2 pc plot viewers", () => viewerCount(page, 2, "pc plot"));
      await session.step(123, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(124, "When user hovers over second pc plot viewer", () => hoverOver(page, el("second pc plot viewer")));
      await session.step(125, "And user drags the max handle of the \"AGE\" range slider of second pc plot viewer by 120 pixels", () => dragRangeHandle(page, "max", "AGE", el("second pc plot viewer"), 120));
      await session.step(126, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(127, "And the \"filtering\" reading of second pc plot viewer should be \"true\"", () => readingReads(page, "filtering", el("second pc plot viewer"), "true"));
      await session.step(128, "And the \"filtering\" reading of first pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("first pc plot viewer"), "false"));
      await session.step(129, "And the \"rows shown\" and \"lines drawn\" readings of first pc plot viewer should be the same", () => readingsEqual(page, "rows shown", "lines drawn", el("first pc plot viewer")));
      await session.step(130, "And the \"rows shown\" reading of first pc plot viewer should be between 1 and 999", () => readingBetween(page, "rows shown", el("first pc plot viewer"), 1, 999));
      await session.step(131, "When user clicks on close icon of second pc plot viewer", () => clickOn(page, el("close icon of second pc plot viewer")));
      await session.step(132, "Then the open tableview should have 1 pc plot viewer", () => viewerCount(page, 1, "pc plot"));
      await session.step(133, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(134, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the plot releases the filter it contributed", async () => {
      await session.step(139, "When user drags the max handle of the \"AGE\" range slider of pc plot viewer by 120 pixels", () => dragRangeHandle(page, "max", "AGE", el("pc plot viewer"), 120));
      await session.step(140, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(141, "When user clicks on close icon of pc plot viewer", () => clickOn(page, el("close icon of pc plot viewer")));
      await session.step(142, "Then pc plot viewer should be absent", () => shouldBe(page, el("pc plot viewer"), "absent"));
      await session.step(143, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(144, "When user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]), [["Column Names","AGE, HEIGHT, WEIGHT"]]);
      await session.step(146, "Then pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(147, "And the \"filtering\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("pc plot viewer"), "false"));
      await session.step(148, "And no errors should have been logged", () => noErrors(page));
      await session.step(149, "When user hovers over pc plot viewer", () => hoverOver(page, el("pc plot viewer")));
      await session.step(150, "And user drags the max handle of the \"AGE\" range slider of pc plot viewer by 120 pixels", () => dragRangeHandle(page, "max", "AGE", el("pc plot viewer"), 120));
      await session.step(151, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(152, "When user sets \"Transformation\" property of pc plot viewer to '[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"SEX\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"AGE\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"WEIGHT\"}]'", () => setProperty(page, "Transformation", el("pc plot viewer"), "[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"SEX\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"AGE\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"WEIGHT\"}]"));
      await session.step(153, "Then the \"transformed\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "transformed", el("pc plot viewer"), "true"));
      await session.step(154, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(155, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
