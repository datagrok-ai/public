/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot-selection.feature
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
import {lineSelected} from '../../../bindings/pc-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {currentRowIs, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {allOfSelected, clearSelection, filterNotNull, filterPasses, filterPassesAll, noneSelected, resetFilter, selectWhereIs, selectedPassFilter, selectedRowCount, tableColumnIncomplete} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, dragSelectionOverArea, eventFired, hasArea, hoverArea, lessHighlight, lessInk, listenFor, moreHighlight, moreInk, noErrors, pickFromContextMenu, propertyShouldBe, readingHigher, readingIs, repaintedBy, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot selection, current row and mouse-over", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-selection.feature", import.meta.url);
  test("PC plot selection, current row and mouse-over", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(27, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]));
    await session.step(29, "Then pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(30, "And the \"rows selected\" reading of pc plot viewer should be 0", () => readingIs(page, "rows selected", el("pc plot viewer"), 0));
    await session.step(31, "And the \"current row\" reading of pc plot viewer should be 1", () => readingIs(page, "current row", el("pc plot viewer"), 1));
    await session.step(32, "And table \"demog-1000\" should have missing values in \"HEIGHT\" column", () => tableColumnIncomplete(page, "demog-1000", "HEIGHT"));
    await session.step(33, "And pc plot viewer should have a \"band \\\"AGE\\\" - \\\"HEIGHT\\\"\" area", () => hasArea(page, el("pc plot viewer"), "band \"AGE\" - \"HEIGHT\""));
    await session.step(34, "And pc plot viewer should have a \"line of row 1000\" area", () => hasArea(page, el("pc plot viewer"), "line of row 1000"));
    await run.scenario("A Shift-drag over a band selects the polylines it crosses", async () => {
      await session.step(37, "When user filters rows where \"HEIGHT\" is not null", () => filterNotNull(page, "HEIGHT"));
      await session.step(38, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
      await session.step(39, "And pc plot viewer should show 872 rows", () => showsRows(page, el("pc plot viewer"), 872));
      await session.step(40, "When user drags a selection box over the \"band \\\"AGE\\\" - \\\"HEIGHT\\\"\" area of pc plot viewer", () => dragSelectionOverArea(page, "band \"AGE\" - \"HEIGHT\"", el("pc plot viewer")));
      await session.step(41, "Then 872 rows should be selected", () => selectedRowCount(page, 872));
      await session.step(42, "And the \"rows selected\" reading of pc plot viewer should be 872", () => readingIs(page, "rows selected", el("pc plot viewer"), 872));
      await session.step(43, "And the line of row 1000 of pc plot viewer should be selected", () => lineSelected(page, 1000, el("pc plot viewer")));
      await session.step(44, "And every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(45, "And pc plot viewer should show more selection highlight than before", () => moreHighlight(page, el("pc plot viewer")));
      await session.step(46, "When user clears the row selection", () => clearSelection(page));
      await session.step(47, "Then no rows should be selected", () => noneSelected(page));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second Shift-drag adds to the selection instead of replacing it", async () => {
      await session.step(51, "Given 872 rows should pass the filter", () => filterPasses(page, 872));
      await session.step(52, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(53, "Then 15 rows should be selected", () => selectedRowCount(page, 15));
      await session.step(54, "When user drags a selection box over the \"band \\\"AGE\\\" - \\\"HEIGHT\\\"\" area of pc plot viewer", () => dragSelectionOverArea(page, "band \"AGE\" - \"HEIGHT\"", el("pc plot viewer")));
      await session.step(55, "Then 874 rows should be selected", () => selectedRowCount(page, 874));
      await session.step(56, "And all rows where \"RACE\" is \"Asian\" should be selected", () => allOfSelected(page, "RACE", "Asian"));
      await session.step(57, "And the line of row 899 of pc plot viewer should be selected", () => lineSelected(page, 899, el("pc plot viewer")));
      await session.step(58, "And the \"rows selected\" reading of pc plot viewer should be higher than before", () => readingHigher(page, "rows selected", el("pc plot viewer")));
      await session.step(59, "When user clears the row selection", () => clearSelection(page));
      await session.step(60, "And user resets the filter", () => resetFilter(page));
      await session.step(61, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(62, "And no rows should be selected", () => noneSelected(page));
      await session.step(63, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on chart space no polyline passes through clears the selection", async () => {
      await session.step(67, "When user selects rows where \"RACE\" is \"Caucasian\"", () => selectWhereIs(page, "RACE", "Caucasian"));
      await session.step(68, "Then 896 rows should be selected", () => selectedRowCount(page, 896));
      await session.step(69, "And pc plot viewer should show more selection highlight than before", () => moreHighlight(page, el("pc plot viewer")));
      await session.step(70, "When user clicks on the \"empty space\" area of pc plot viewer", () => clickArea(page, "empty space", el("pc plot viewer")));
      await session.step(71, "Then no rows should be selected", () => noneSelected(page));
      await session.step(72, "And the \"rows selected\" reading of pc plot viewer should be 0", () => readingIs(page, "rows selected", el("pc plot viewer"), 0));
      await session.step(73, "And pc plot viewer should show less selection highlight than before", () => lessHighlight(page, el("pc plot viewer")));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a polyline makes that row current", async () => {
      await session.step(77, "Given user listens for \"d4-pc-plot-on-line-clicked\" event on pc plot viewer", () => listenFor(page, "d4-pc-plot-on-line-clicked", el("pc plot viewer")));
      await session.step(78, "And no rows should be selected", () => noneSelected(page));
      await session.step(79, "Then the \"current row\" reading of pc plot viewer should be 1", () => readingIs(page, "current row", el("pc plot viewer"), 1));
      await session.step(80, "When user clicks on the \"line of row 1000\" area of pc plot viewer", () => clickArea(page, "line of row 1000", el("pc plot viewer")));
      await session.step(81, "Then \"d4-pc-plot-on-line-clicked\" event should have fired on pc plot viewer", () => eventFired(page, "d4-pc-plot-on-line-clicked", el("pc plot viewer")));
      await session.step(82, "And the \"current row\" reading of pc plot viewer should be 1000", () => readingIs(page, "current row", el("pc plot viewer"), 1000));
      await session.step(83, "And row 1000 should be current", () => currentRowIs(page, 1000));
      await session.step(84, "And no rows should be selected", () => noneSelected(page));
      await session.step(85, "When user clicks on the \"line of row 999\" area of pc plot viewer", () => clickArea(page, "line of row 999", el("pc plot viewer")));
      await session.step(86, "Then the \"current row\" reading of pc plot viewer should be 999", () => readingIs(page, "current row", el("pc plot viewer"), 999));
      await session.step(87, "When user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(88, "Then the \"current row\" reading of pc plot viewer should be 1", () => readingIs(page, "current row", el("pc plot viewer"), 1));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hovering a polyline makes it the mouse-over row", async () => {
      await session.step(92, "Given user listens for \"d4-pc-plot-on-line-hovered\" event on pc plot viewer", () => listenFor(page, "d4-pc-plot-on-line-hovered", el("pc plot viewer")));
      await session.step(93, "When user hovers over the \"empty space\" area of pc plot viewer", () => hoverArea(page, "empty space", el("pc plot viewer")));
      await session.step(94, "Then the \"hovered row\" reading of pc plot viewer should be 0", () => readingIs(page, "hovered row", el("pc plot viewer"), 0));
      await session.step(95, "When user hovers over the \"line of row 1000\" area of pc plot viewer", () => hoverArea(page, "line of row 1000", el("pc plot viewer")));
      await session.step(96, "Then \"d4-pc-plot-on-line-hovered\" event should have fired on pc plot viewer", () => eventFired(page, "d4-pc-plot-on-line-hovered", el("pc plot viewer")));
      await session.step(97, "And the \"hovered row\" reading of pc plot viewer should be 1000", () => readingIs(page, "hovered row", el("pc plot viewer"), 1000));
      await session.step(98, "When user hovers over the \"empty space\" area of pc plot viewer", () => hoverArea(page, "empty space", el("pc plot viewer")));
      await session.step(99, "Then the \"hovered row\" reading of pc plot viewer should be 0", () => readingIs(page, "hovered row", el("pc plot viewer"), 0));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show All Lines off leaves the current, selected and hovered polylines", async () => {
      await session.step(103, "When user clears the row selection", () => clearSelection(page));
      await session.step(104, "And user hovers over the \"empty space\" area of pc plot viewer", () => hoverArea(page, "empty space", el("pc plot viewer")));
      await session.step(105, "Then the \"hovered row\" reading of pc plot viewer should be 0", () => readingIs(page, "hovered row", el("pc plot viewer"), 0));
      await session.step(106, "And the \"lines drawn\" reading of pc plot viewer should be 1000", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1000));
      await session.step(107, "When user sets \"Show All Lines\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show All Lines", el("pc plot viewer"), "false"));
      await session.step(108, "Then the \"lines drawn\" reading of pc plot viewer should be 1", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1));
      await session.step(109, "And pc plot viewer should have less ink than before", () => lessInk(page, el("pc plot viewer")));
      await session.step(110, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(111, "Then the \"lines drawn\" reading of pc plot viewer should be 16", () => readingIs(page, "lines drawn", el("pc plot viewer"), 16));
      await session.step(112, "And pc plot viewer should have more ink than before", () => moreInk(page, el("pc plot viewer")));
      await session.step(113, "When user hovers over the \"line of row 1000\" area of pc plot viewer", () => hoverArea(page, "line of row 1000", el("pc plot viewer")));
      await session.step(114, "Then the \"hovered row\" reading of pc plot viewer should be 1000", () => readingIs(page, "hovered row", el("pc plot viewer"), 1000));
      await session.step(115, "And the \"lines drawn\" reading of pc plot viewer should be 17", () => readingIs(page, "lines drawn", el("pc plot viewer"), 17));
      await session.step(116, "When user clears the row selection", () => clearSelection(page));
      await session.step(117, "Then the \"lines drawn\" reading of pc plot viewer should be 2", () => readingIs(page, "lines drawn", el("pc plot viewer"), 2));
      await session.step(118, "When user hovers over the \"empty space\" area of pc plot viewer", () => hoverArea(page, "empty space", el("pc plot viewer")));
      await session.step(119, "And user sets \"Show All Lines\" property of pc plot viewer to \"true\"", () => setProperty(page, "Show All Lines", el("pc plot viewer"), "true"));
      await session.step(120, "Then the \"lines drawn\" reading of pc plot viewer should be 1000", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1000));
      await session.step(121, "And pc plot viewer should have more ink than before", () => moreInk(page, el("pc plot viewer")));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Selection menu drives the lines and the overlay", async () => {
      await session.step(125, "When user clears the row selection", () => clearSelection(page));
      await session.step(126, "And user hovers over the \"empty space\" area of pc plot viewer", () => hoverArea(page, "empty space", el("pc plot viewer")));
      await session.step(127, "Then the \"lines drawn\" reading of pc plot viewer should be 1000", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1000));
      await session.step(128, "When user picks \"Selection > Show All Lines\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Selection > Show All Lines", el("pc plot viewer")));
      await session.step(129, "Then \"Show All Lines\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Show All Lines", el("pc plot viewer"), "false"));
      await session.step(130, "And pc plot viewer should have less ink than before", () => lessInk(page, el("pc plot viewer")));
      await session.step(131, "When user picks \"Selection > Show Current Line\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Selection > Show Current Line", el("pc plot viewer")));
      await session.step(132, "Then \"Show Current Line\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Current Line", el("pc plot viewer"), "false"));
      await session.step(133, "And pc plot viewer should have repainted by at least 200 pixels", () => repaintedBy(page, el("pc plot viewer"), 200));
      await session.step(134, "When user picks \"Selection > Show Current Line\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Selection > Show Current Line", el("pc plot viewer")));
      await session.step(135, "Then \"Show Current Line\" property of pc plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Current Line", el("pc plot viewer"), "true"));
      await session.step(136, "And pc plot viewer should have repainted by at least 200 pixels", () => repaintedBy(page, el("pc plot viewer"), 200));
      await session.step(137, "When user picks \"Selection > Show Mouse Over Line\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Selection > Show Mouse Over Line", el("pc plot viewer")));
      await session.step(138, "Then \"Show Mouse Over Line\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Mouse Over Line", el("pc plot viewer"), "false"));
      await session.step(139, "When user picks \"Selection > Show Mouse Over Line\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Selection > Show Mouse Over Line", el("pc plot viewer")));
      await session.step(140, "Then \"Show Mouse Over Line\" property of pc plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Mouse Over Line", el("pc plot viewer"), "true"));
      await session.step(141, "When user picks \"Selection > Show All Lines\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Selection > Show All Lines", el("pc plot viewer")));
      await session.step(142, "Then \"Show All Lines\" property of pc plot viewer should be \"true\"", () => propertyShouldBe(page, "Show All Lines", el("pc plot viewer"), "true"));
      await session.step(143, "And the \"lines drawn\" reading of pc plot viewer should be 1000", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1000));
      await session.step(144, "And pc plot viewer should have more ink than before", () => moreInk(page, el("pc plot viewer")));
      await session.step(145, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
