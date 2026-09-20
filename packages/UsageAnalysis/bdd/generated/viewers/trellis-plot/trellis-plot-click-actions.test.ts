/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-click-actions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
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
import {clickOn, hoverOver, pressKeyIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, allOfFiltered, allOfSelected, filterPasses, filterPassesAll, noneOfFiltered, noneSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, clickAreaHolding, eventFired, hasArea, listenFor, noErrors, propertyShouldBe, readingAsRemembered, readingIs, readingReads, readingsDiffer, rememberReading, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot click actions", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-click-actions.feature", import.meta.url);
  test("Trellis plot click actions", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 14, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]));
    await session.step(20, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await session.step(21, "And \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
    await run.scenario("Switching On Click alone draws nothing new", async () => {
      await session.step(24, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(25, "And user sets \"On Click\" property of trellis plot viewer to \"Select\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Select"));
      await session.step(26, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(27, "And the \"cell signature F | Caucasian\" and \"cell signature M | Asian\" readings of trellis plot viewer should differ", () => readingsDiffer(page, "cell signature F | Caucasian", "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(28, "And the \"columns\" reading of trellis plot viewer should be 2", () => readingIs(page, "columns", el("trellis plot viewer"), 2));
      await session.step(29, "And the \"rows\" reading of trellis plot viewer should be 4", () => readingIs(page, "rows", el("trellis plot viewer"), 4));
      await session.step(30, "And the \"current cell\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "current cell", el("trellis plot viewer"), ""));
      await session.step(31, "And no rows should be selected", () => noneSelected(page));
      await session.step(32, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell click selects exactly its rows", async () => {
      await session.step(36, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(37, "Then 480 rows should be selected", () => selectedRowCount(page, 480));
      await session.step(38, "And the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(39, "And the \"columns\" reading of trellis plot viewer should be 2", () => readingIs(page, "columns", el("trellis plot viewer"), 2));
      await session.step(40, "And the \"rows\" reading of trellis plot viewer should be 4", () => readingIs(page, "rows", el("trellis plot viewer"), 4));
      await session.step(41, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(42, "And trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An inner type change does not disturb the selection", async () => {
      await session.step(46, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Bar chart\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Bar chart"));
      await session.step(47, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Bar chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Bar chart"));
      await session.step(48, "And 480 rows should be selected", () => selectedRowCount(page, 480));
      await session.step(49, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Scatter plot\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(50, "Then 480 rows should be selected", () => selectedRowCount(page, 480));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Another cell replaces the selection", async () => {
      await session.step(54, "When user clicks on the \"cell M | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell M | Caucasian", el("trellis plot viewer")));
      await session.step(55, "Then 416 rows should be selected", () => selectedRowCount(page, 416));
      await session.step(56, "And the \"current cell\" reading of trellis plot viewer should be \"M | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Caucasian"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control adds another cell's rows to the selection", async () => {
      await session.step(60, "When user clicks on the \"cell M | Asian\" area of trellis plot viewer holding Control", () => clickAreaHolding(page, "cell M | Asian", el("trellis plot viewer"), "Control"));
      await session.step(61, "Then 424 rows should be selected", () => selectedRowCount(page, 424));
      await session.step(62, "And the \"current cell\" reading of trellis plot viewer should be \"M | Asian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Asian"));
      await session.step(63, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(64, "Then no rows should be selected", () => noneSelected(page));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control on a cell the filter panel emptied still adds the cell's rows", async () => {
      await session.step(68, "Given \"Row Source\" property of trellis plot viewer should be \"Filtered\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "Filtered"));
      await session.step(69, "When user sets \"Row Source\" property of trellis plot viewer to \"All\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "All"));
      await session.step(70, "And user clicks on the \"cell M | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell M | Caucasian", el("trellis plot viewer")));
      await session.step(71, "Then 416 rows should be selected", () => selectedRowCount(page, 416));
      await session.step(72, "When user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(73, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(74, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(75, "And trellis plot viewer should have a \"cell F | Caucasian\" area", () => hasArea(page, el("trellis plot viewer"), "cell F | Caucasian"));
      await session.step(76, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer holding Control", () => clickAreaHolding(page, "cell F | Caucasian", el("trellis plot viewer"), "Control"));
      await session.step(77, "Then 896 rows should be selected", () => selectedRowCount(page, 896));
      await session.step(78, "And all rows where \"RACE\" is \"Caucasian\" should be selected", () => allOfSelected(page, "RACE", "Caucasian"));
      await session.step(79, "And the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(80, "When user adds a categorical filter on \"SEX\" keeping \"F, M\"", () => addCategoricalFilter(page, "SEX", "F, M"));
      await session.step(81, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(82, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(83, "And user sets \"Row Source\" property of trellis plot viewer to \"Filtered\"", () => setProperty(page, "Row Source", el("trellis plot viewer"), "Filtered"));
      await session.step(84, "Then no rows should be selected", () => noneSelected(page));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Control on an empty cell adds nothing but takes the current cell", async () => {
      await session.step(88, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","RACE"],["Y Column Names","SEVERITY"],["Pack Categories","false"]]));
      await session.step(92, "Then the \"cells\" reading of trellis plot viewer should be 20", () => readingIs(page, "cells", el("trellis plot viewer"), 20));
      await session.step(93, "And the \"cells drawn\" reading of trellis plot viewer should be 17", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 17));
      await session.step(94, "And trellis plot viewer should have a \"cell Asian | Critical\" area", () => hasArea(page, el("trellis plot viewer"), "cell Asian | Critical"));
      await session.step(95, "When user clicks on the \"cell Caucasian | Critical\" area of trellis plot viewer", () => clickArea(page, "cell Caucasian | Critical", el("trellis plot viewer")));
      await session.step(96, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
      await session.step(97, "When user clicks on the \"cell Asian | Critical\" area of trellis plot viewer holding Control", () => clickAreaHolding(page, "cell Asian | Critical", el("trellis plot viewer"), "Control"));
      await session.step(98, "Then 5 rows should be selected", () => selectedRowCount(page, 5));
      await session.step(99, "And the \"current cell\" reading of trellis plot viewer should be \"Asian | Critical\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "Asian | Critical"));
      await session.step(100, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(101, "Then no rows should be selected", () => noneSelected(page));
      await session.step(102, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Pack Categories","true"]]));
      await session.step(106, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell click filters to exactly its rows", async () => {
      await session.step(110, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(111, "And user remembers the \"cell signature M | Asian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(112, "And user sets \"On Click\" property of trellis plot viewer to \"Filter\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Filter"));
      await session.step(113, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(114, "And the \"cell signature M | Asian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(115, "And \"Row Source\" property of trellis plot viewer should be \"All\"", () => propertyShouldBe(page, "Row Source", el("trellis plot viewer"), "All"));
      await session.step(116, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(117, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(118, "Then 480 rows should pass the filter", () => filterPasses(page, 480));
      await session.step(119, "And no rows where \"SEX\" is \"M\" should pass the filter", () => noneOfFiltered(page, "SEX", "M"));
      await session.step(120, "And no rows where \"RACE\" is \"Asian\" should pass the filter", () => noneOfFiltered(page, "RACE", "Asian"));
      await session.step(121, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Escape drops the trellis contribution", async () => {
      await session.step(125, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(126, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A filter card and a cell click compose", async () => {
      await session.step(130, "When user adds a categorical filter on \"DIS_POP\" keeping \"RA\"", () => addCategoricalFilter(page, "DIS_POP", "RA"));
      await session.step(131, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(132, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(133, "Then 271 rows should pass the filter", () => filterPasses(page, 271));
      await session.step(134, "And no rows where \"SEX\" is \"M\" should pass the filter", () => noneOfFiltered(page, "SEX", "M"));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Changing a split column drops the trellis contribution and keeps the card", async () => {
      await session.step(138, "When user sets \"X Column Names\" property of trellis plot viewer to \"CONTROL\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "CONTROL"));
      await session.step(139, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(140, "And all rows where \"DIS_POP\" is \"RA\" should pass the filter", () => allOfFiltered(page, "DIS_POP", "RA"));
      await session.step(141, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX"));
      await session.step(142, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(143, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the card restores every row", async () => {
      await session.step(146, "When user hovers over \"DIS_POP\" filter card", () => hoverOver(page, el("\"DIS_POP\" filter card")));
      await session.step(147, "And user clicks on close of \"DIS_POP\" filter card", () => clickOn(page, el("close of \"DIS_POP\" filter card")));
      await session.step(148, "Then \"DIS_POP\" filter card should be absent", () => shouldBe(page, el("\"DIS_POP\" filter card"), "absent"));
      await session.step(149, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(150, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On Click None leaves both channels silent", async () => {
      await session.step(153, "When user sets \"On Click\" property of trellis plot viewer to \"None\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(154, "And user clicks on the \"cell M | Black\" area of trellis plot viewer", () => clickArea(page, "cell M | Black", el("trellis plot viewer")));
      await session.step(155, "Then no rows should be selected", () => noneSelected(page));
      await session.step(156, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(157, "And the \"current cell\" reading of trellis plot viewer should be \"M | Black\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Black"));
      await session.step(158, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Both trellis events fire off one cell", async () => {
      await session.step(161, "Given user listens for \"d4-trellis-plot-current-cell-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-current-cell-changed", el("trellis plot viewer")));
      await session.step(162, "And user listens for \"d4-trellis-plot-inner-viewer-clicked\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-inner-viewer-clicked", el("trellis plot viewer")));
      await session.step(163, "When user clicks on the \"cell body F | Black\" area of trellis plot viewer", () => clickArea(page, "cell body F | Black", el("trellis plot viewer")));
      await session.step(164, "And user clicks on the \"cell F | Black\" area of trellis plot viewer", () => clickArea(page, "cell F | Black", el("trellis plot viewer")));
      await session.step(165, "Then \"d4-trellis-plot-current-cell-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-current-cell-changed", el("trellis plot viewer")));
      await session.step(166, "And \"d4-trellis-plot-inner-viewer-clicked\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-inner-viewer-clicked", el("trellis plot viewer")));
      await session.step(167, "And the \"current cell\" reading of trellis plot viewer should be \"F | Black\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Black"));
      await session.step(168, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
