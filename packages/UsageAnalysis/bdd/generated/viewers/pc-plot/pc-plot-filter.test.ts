/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pc-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {dragAxisHandle} from '../../../bindings/pc-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addRangeFilter, filterPasses, filterPassesAll, filterPassesFewer, openFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, eventFired, hasArea, lessInk, listenFor, moreInk, noErrors, pickFromContextMenu, propertyShouldBe, readingAsRemembered, readingHigher, readingIs, readingLower, readingNotAsRemembered, readingReads, rememberReading, setProperty, showsFewerRows, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot in-chart range filter", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-filter.feature", import.meta.url);
  test("PC plot in-chart range filter", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]));
    await session.step(20, "Then pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(21, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(22, "And the \"filtering\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("pc plot viewer"), "false"));
    await session.step(23, "And the \"range max of \\\"AGE\\\"\" reading of pc plot viewer should be 89", () => readingIs(page, "range max of \"AGE\"", el("pc plot viewer"), 89));
    await run.scenario("Dragging the AGE max handle filters the table, Reset View gives it back", async () => {
      await session.step(26, "Given user listens for \"d4-pc-plot-reset-view\" event on pc plot viewer", () => listenFor(page, "d4-pc-plot-reset-view", el("pc plot viewer")));
      await session.step(27, "When user hovers over pc plot viewer", () => hoverOver(page, el("pc plot viewer")));
      await session.step(28, "Then pc plot viewer should have a \"range max handle \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "range max handle \"AGE\""));
      await session.step(29, "When user drags the max handle of the \"AGE\" axis range slider of pc plot viewer by 120 pixels", () => dragAxisHandle(page, "max", "AGE", 120));
      await session.step(30, "Then the \"filtering\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "filtering", el("pc plot viewer"), "true"));
      await session.step(31, "And the \"range max of \\\"AGE\\\"\" reading of pc plot viewer should be lower than before", () => readingLower(page, "range max of \"AGE\"", el("pc plot viewer")));
      await session.step(32, "And fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(33, "And pc plot viewer should show fewer rows than before", () => showsFewerRows(page, el("pc plot viewer")));
      await session.step(34, "And the \"lines drawn\" reading of pc plot viewer should be lower than before", () => readingLower(page, "lines drawn", el("pc plot viewer")));
      await session.step(35, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(36, "Then \"d4-pc-plot-reset-view\" event should have fired on pc plot viewer", () => eventFired(page, "d4-pc-plot-reset-view", el("pc plot viewer")));
      await session.step(37, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(38, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(39, "And the \"filtering\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("pc plot viewer"), "false"));
      await session.step(40, "And the \"range max of \\\"AGE\\\"\" reading of pc plot viewer should be 89", () => readingIs(page, "range max of \"AGE\"", el("pc plot viewer"), 89));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Filtered Out Lines draws the rows the slider dropped", async () => {
      await session.step(44, "When user drags the max handle of the \"AGE\" axis range slider of pc plot viewer by 120 pixels", () => dragAxisHandle(page, "max", "AGE", 120));
      await session.step(45, "Then pc plot viewer should show fewer rows than before", () => showsFewerRows(page, el("pc plot viewer")));
      await session.step(46, "And the \"filtered out lines drawn\" reading of pc plot viewer should be 0", () => readingIs(page, "filtered out lines drawn", el("pc plot viewer"), 0));
      await session.step(47, "When user picks \"Filter > Show Filtered Out Lines\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Filter > Show Filtered Out Lines", el("pc plot viewer")));
      await session.step(48, "Then \"Show Filtered Out Lines\" property of pc plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Filtered Out Lines", el("pc plot viewer"), "true"));
      await session.step(49, "And the \"filtered out lines drawn\" reading of pc plot viewer should be higher than before", () => readingHigher(page, "filtered out lines drawn", el("pc plot viewer")));
      await session.step(50, "And pc plot viewer should have more ink than before", () => moreInk(page, el("pc plot viewer")));
      await session.step(51, "When user picks \"Filter > Show Filtered Out Lines\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Filter > Show Filtered Out Lines", el("pc plot viewer")));
      await session.step(52, "Then \"Show Filtered Out Lines\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Filtered Out Lines", el("pc plot viewer"), "false"));
      await session.step(53, "And the \"filtered out lines drawn\" reading of pc plot viewer should be 0", () => readingIs(page, "filtered out lines drawn", el("pc plot viewer"), 0));
      await session.step(54, "And pc plot viewer should have less ink than before", () => lessInk(page, el("pc plot viewer")));
      await session.step(55, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(56, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter card and the in-chart slider compose with AND", async () => {
      await session.step(60, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(61, "And user adds a range filter on \"AGE\" from 30 to 50", () => addRangeFilter(page, "AGE", 30, 50));
      await session.step(62, "Then 494 rows should pass the filter", () => filterPasses(page, 494));
      await session.step(63, "And pc plot viewer should show 494 rows", () => showsRows(page, el("pc plot viewer"), 494));
      await session.step(64, "When user drags the max handle of the \"HEIGHT\" axis range slider of pc plot viewer by 120 pixels", () => dragAxisHandle(page, "max", "HEIGHT", 120));
      await session.step(65, "Then fewer than 494 rows should pass the filter", () => filterPassesFewer(page, 494));
      await session.step(66, "And pc plot viewer should show fewer rows than before", () => showsFewerRows(page, el("pc plot viewer")));
      await session.step(67, "And the \"filtering\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "filtering", el("pc plot viewer"), "true"));
      await session.step(68, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(69, "Then 494 rows should pass the filter", () => filterPasses(page, 494));
      await session.step(70, "And pc plot viewer should show 494 rows", () => showsRows(page, el("pc plot viewer"), 494));
      await session.step(71, "And the \"filtering\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("pc plot viewer"), "false"));
      await session.step(72, "When user drags the max handle of the \"HEIGHT\" axis range slider of pc plot viewer by 120 pixels", () => dragAxisHandle(page, "max", "HEIGHT", 120));
      await session.step(73, "Then fewer than 494 rows should pass the filter", () => filterPassesFewer(page, 494));
      await session.step(74, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(75, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(76, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(77, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(78, "And the \"filtering\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("pc plot viewer"), "false"));
      await session.step(79, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(80, "Then filter panel should be absent", () => shouldBe(page, el("filter panel"), "absent"));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second range filter after a DateTime colour split still filters (GROK-18489)", async () => {
      await session.step(84, "When user sets \"Color\" property of pc plot viewer to \"STARTED\"", () => setProperty(page, "Color", el("pc plot viewer"), "STARTED"));
      await session.step(85, "And user drags the max handle of the \"AGE\" axis range slider of pc plot viewer by 120 pixels", () => dragAxisHandle(page, "max", "AGE", 120));
      await session.step(86, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(87, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(88, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(89, "When user drags the min handle of the \"AGE\" axis range slider of pc plot viewer by -120 pixels", () => dragAxisHandle(page, "min", "AGE", -120));
      await session.step(90, "Then the \"range min of \\\"AGE\\\"\" reading of pc plot viewer should be higher than before", () => readingHigher(page, "range min of \"AGE\"", el("pc plot viewer")));
      await session.step(91, "And fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(92, "When user remembers the \"rows shown\" reading of pc plot viewer", () => rememberReading(page, "rows shown", el("pc plot viewer")));
      await session.step(93, "And user drags the max handle of the \"AGE\" axis range slider of pc plot viewer by 100 pixels", () => dragAxisHandle(page, "max", "AGE", 100));
      await session.step(94, "Then the \"rows shown\" reading of pc plot viewer should not be as remembered", () => readingNotAsRemembered(page, "rows shown", el("pc plot viewer")));
      await session.step(95, "And the \"range max of \\\"AGE\\\"\" reading of pc plot viewer should be lower than before", () => readingLower(page, "range max of \"AGE\"", el("pc plot viewer")));
      await session.step(96, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(97, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(98, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(99, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Changing a histogram's column leaves the plot's filter alone (github-972)", async () => {
      await session.step(102, "When user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"]]));
      await session.step(104, "And user drags the max handle of the \"AGE\" axis range slider of pc plot viewer by 120 pixels", () => dragAxisHandle(page, "max", "AGE", 120));
      await session.step(105, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(106, "When user remembers the \"rows shown\" reading of pc plot viewer", () => rememberReading(page, "rows shown", el("pc plot viewer")));
      await session.step(107, "And user sets \"Value\" property of histogram viewer to \"HEIGHT\"", () => setProperty(page, "Value", el("histogram viewer"), "HEIGHT"));
      await session.step(108, "Then \"Value\" property of histogram viewer should be \"HEIGHT\"", () => propertyShouldBe(page, "Value", el("histogram viewer"), "HEIGHT"));
      await session.step(109, "And the \"rows shown\" reading of pc plot viewer should be as remembered", () => readingAsRemembered(page, "rows shown", el("pc plot viewer")));
      await session.step(110, "And the \"filtering\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "filtering", el("pc plot viewer"), "true"));
      await session.step(111, "When user clicks on close icon of histogram viewer", () => clickOn(page, el("close icon of histogram viewer")));
      await session.step(112, "Then histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
      await session.step(113, "When user picks \"Reset View\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Reset View", el("pc plot viewer")));
      await session.step(114, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the plot releases the filter it contributed", async () => {
      await session.step(118, "When user drags the max handle of the \"AGE\" axis range slider of pc plot viewer by 120 pixels", () => dragAxisHandle(page, "max", "AGE", 120));
      await session.step(119, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(120, "When user clicks on close icon of pc plot viewer", () => clickOn(page, el("close icon of pc plot viewer")));
      await session.step(121, "Then pc plot viewer should be absent", () => shouldBe(page, el("pc plot viewer"), "absent"));
      await session.step(122, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(123, "When user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]));
      await session.step(125, "Then pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(126, "And the \"filtering\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("pc plot viewer"), "false"));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
