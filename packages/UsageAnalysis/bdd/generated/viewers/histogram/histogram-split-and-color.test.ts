/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/histogram/histogram-split-and-color.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.histogram, core#998, public#130]
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
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, colorCategorical, colorOff, filterIsExactly, filterPasses, filterPassesAll, noneSelected, onlyOfSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaPainted, areasDiffer, clickArea, enterIntoArea, eventFired, hasArea, hasNoArea, hoverArea, legendLists, lessInk, listenFor, moreInk, noErrors, pointerAway, readingIs, readingReads, repaintedBy, resizeTo, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Histogram split, stacking and color coding", () => {
  const session = feature(test, "features/viewers/histogram/histogram-split-and-color.feature", import.meta.url);
  test("Histogram split, stacking and color coding", {tag: ["@journey", "@viewers", "@realizes:viewers.histogram", "@realizes:core#998", "@realizes:public#130"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"]]), [["Value","AGE"]]);
    await session.step(17, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(18, "And user clicks on settings icon of histogram viewer", () => clickOn(page, el("settings icon of histogram viewer")));
    await session.step(19, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(20, "And histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
    await run.scenario("Color coding by a numerical column", async () => {
      await session.step(23, "Then \"Color\" property in context panel should be enabled", () => shouldBe(page, el("\"Color\" property in context panel"), "enabled"));
      await session.step(24, "And \"Color Aggr Type\" property in context panel should be disabled", () => shouldBe(page, el("\"Color Aggr Type\" property in context panel"), "disabled"));
      await session.step(25, "When user sets \"Color\" property of histogram viewer to \"HEIGHT\"", () => setProperty(page, "Color", el("histogram viewer"), "HEIGHT"));
      await session.step(26, "Then histogram viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("histogram viewer"), 2000));
      await session.step(27, "And the \"bin 1\" area of histogram viewer should contain the color \"#FF0000\"", () => areaColor(page, "bin 1", el("histogram viewer"), "#FF0000"));
      await session.step(28, "And the \"bin 17\" area of histogram viewer should contain the color \"#0000FF\"", () => areaColor(page, "bin 17", el("histogram viewer"), "#0000FF"));
      await session.step(29, "And the \"bin 1\" and \"bin 17\" areas of histogram viewer should be painted in different colors", () => areasDiffer(page, "bin 1", "bin 17", el("histogram viewer")));
      await session.step(30, "And \"Color Aggr Type\" property in context panel should be enabled", () => shouldBe(page, el("\"Color Aggr Type\" property in context panel"), "enabled"));
      await session.step(31, "And \"Invert Color Scheme\" property in context panel should be enabled", () => shouldBe(page, el("\"Invert Color Scheme\" property in context panel"), "enabled"));
      await session.step(32, "When user sets \"Color Aggr Type\" property of histogram viewer to \"stdev\"", () => setProperty(page, "Color Aggr Type", el("histogram viewer"), "stdev"));
      await session.step(33, "Then histogram viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("histogram viewer"), 500));
      await session.step(34, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Color Aggr Type","avg"],["Invert Color Scheme","true"]]), [["Color Aggr Type","avg"],["Invert Color Scheme","true"]]);
      await session.step(37, "Then histogram viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("histogram viewer"), 2000));
      await session.step(38, "And the \"bin 1\" area of histogram viewer should contain the color \"#0000FF\"", () => areaColor(page, "bin 1", el("histogram viewer"), "#0000FF"));
      await session.step(39, "And the \"bin 17\" area of histogram viewer should contain the color \"#FF0000\"", () => areaColor(page, "bin 17", el("histogram viewer"), "#FF0000"));
      await session.step(40, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Invert Color Scheme","false"],["Color",""]]), [["Invert Color Scheme","false"],["Color",""]]);
      await session.step(43, "Then histogram viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("histogram viewer"), 2000));
      await session.step(44, "And \"Color Aggr Type\" property in context panel should be disabled", () => shouldBe(page, el("\"Color Aggr Type\" property in context panel"), "disabled"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A split disables color coding and draws lines instead of bars", async () => {
      await session.step(48, "When user sets \"Color\" property of histogram viewer to \"HEIGHT\"", () => setProperty(page, "Color", el("histogram viewer"), "HEIGHT"));
      await session.step(49, "And user sets \"Split\" property of histogram viewer to \"RACE\"", () => setProperty(page, "Split", el("histogram viewer"), "RACE"));
      await session.step(50, "Then \"Color\" property in context panel should be disabled", () => shouldBe(page, el("\"Color\" property in context panel"), "disabled"));
      await session.step(51, "And \"Color Aggr Type\" property in context panel should be disabled", () => shouldBe(page, el("\"Color Aggr Type\" property in context panel"), "disabled"));
      await session.step(52, "And \"Invert Color Scheme\" property in context panel should be disabled", () => shouldBe(page, el("\"Invert Color Scheme\" property in context panel"), "disabled"));
      await session.step(53, "And histogram viewer should not have a \"bin 8\" area", () => hasNoArea(page, el("histogram viewer"), "bin 8"));
      await session.step(54, "And histogram viewer should have a \"line Caucasian\" area", () => hasArea(page, el("histogram viewer"), "line Caucasian"));
      await session.step(55, "And histogram viewer should have a \"line Asian\" area", () => hasArea(page, el("histogram viewer"), "line Asian"));
      await session.step(56, "When user sets \"Split\" property of histogram viewer to \"\"", () => setProperty(page, "Split", el("histogram viewer"), ""));
      await session.step(57, "Then \"Color\" property in context panel should be enabled", () => shouldBe(page, el("\"Color\" property in context panel"), "enabled"));
      await session.step(58, "And histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
      await session.step(59, "When user sets \"Color\" property of histogram viewer to \"\"", () => setProperty(page, "Color", el("histogram viewer"), ""));
      await session.step(60, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A split line answers to the pointer as a whole category", async () => {
      await session.step(63, "Given user listens for \"d4-histogram-mouse-over-line\" event on histogram viewer", () => listenFor(page, "d4-histogram-mouse-over-line", el("histogram viewer")));
      await session.step(64, "And user listens for \"d4-histogram-select-line\" event on histogram viewer", () => listenFor(page, "d4-histogram-select-line", el("histogram viewer")));
      await session.step(65, "When user sets \"Split\" property of histogram viewer to \"SEX\"", () => setProperty(page, "Split", el("histogram viewer"), "SEX"));
      await session.step(66, "Then histogram viewer should have a \"line F point\" area", () => hasArea(page, el("histogram viewer"), "line F point"));
      await session.step(67, "And histogram viewer should have a \"line M point\" area", () => hasArea(page, el("histogram viewer"), "line M point"));
      await session.step(68, "When user hovers over the \"line F point\" area of histogram viewer", () => hoverArea(page, "line F point", el("histogram viewer")));
      await session.step(69, "Then \"d4-histogram-mouse-over-line\" event should have fired on histogram viewer", () => eventFired(page, "d4-histogram-mouse-over-line", el("histogram viewer")));
      await session.step(70, "When user clicks on the \"line F point\" area of histogram viewer", () => clickArea(page, "line F point", el("histogram viewer")));
      await session.step(71, "Then \"d4-histogram-select-line\" event should have fired on histogram viewer", () => eventFired(page, "d4-histogram-select-line", el("histogram viewer")));
      await session.step(72, "And only rows where \"SEX\" is \"F\" should be selected", () => onlyOfSelected(page, "SEX", "F"));
      await session.step(73, "And 553 rows should be selected", () => selectedRowCount(page, 553));
      await session.step(74, "When user clears the row selection", () => clearSelection(page));
      await session.step(75, "And user moves the pointer away from histogram viewer", () => pointerAway(page, el("histogram viewer")));
      await session.step(76, "And user sets \"Split\" property of histogram viewer to \"\"", () => setProperty(page, "Split", el("histogram viewer"), ""));
      await session.step(77, "Then no rows should be selected", () => noneSelected(page));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Split Stack brings the bars back as one segment per category", async () => {
      await session.step(81, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Split","SEX"],["Split Stack","true"]]), [["Split","SEX"],["Split Stack","true"]]);
      await session.step(84, "Then histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
      await session.step(85, "And histogram viewer should have a \"bin 8 | F\" area", () => hasArea(page, el("histogram viewer"), "bin 8 | F"));
      await session.step(86, "And histogram viewer should have a \"bin 8 | M\" area", () => hasArea(page, el("histogram viewer"), "bin 8 | M"));
      await session.step(87, "And the \"bin 8 | F\" and \"bin 8 | M\" areas of histogram viewer should be painted in different colors", () => areasDiffer(page, "bin 8 | F", "bin 8 | M", el("histogram viewer")));
      await session.step(88, "And the legend of histogram viewer should list 2 items", () => legendLists(page, el("histogram viewer"), 2));
      await session.step(89, "When user sets \"Show Values\" property of histogram viewer to \"true\"", () => setProperty(page, "Show Values", el("histogram viewer"), "true"));
      await session.step(90, "Then histogram viewer should have a \"bin labels\" area", () => hasArea(page, el("histogram viewer"), "bin labels"));
      await session.step(91, "And the \"bin labels\" area of histogram viewer should be painted", () => areaPainted(page, "bin labels", el("histogram viewer")));
      await session.step(92, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Show Values",""],["Split Stack","false"],["Split",""]]), [["Show Values",""],["Split Stack","false"],["Split",""]]);
      await session.step(96, "Then histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Normalize Values decides what the vertical axis is scaled to", async () => {
      await session.step(100, "Given histogram viewer should have a \"y axis\" area", () => hasArea(page, el("histogram viewer"), "y axis"));
      await session.step(101, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Split","SEX"],["Normalize Values","true"]]), [["Split","SEX"],["Normalize Values","true"]]);
      await session.step(104, "Then histogram viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("histogram viewer"), 500));
      await session.step(105, "And histogram viewer should not have a \"y axis\" area", () => hasNoArea(page, el("histogram viewer"), "y axis"));
      await session.step(106, "When user sets \"Normalize Values\" property of histogram viewer to \"false\"", () => setProperty(page, "Normalize Values", el("histogram viewer"), "false"));
      await session.step(107, "Then histogram viewer should have a \"y axis\" area", () => hasArea(page, el("histogram viewer"), "y axis"));
      await session.step(108, "And the \"y axis max\" reading of histogram viewer should be 53", () => readingIs(page, "y axis max", el("histogram viewer"), 53));
      await session.step(109, "And histogram viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("histogram viewer"), 500));
      await session.step(110, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Normalize Values","true"],["Split",""]]), [["Normalize Values","true"],["Split",""]]);
      await session.step(113, "Then histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Distribution lines and markers under a split", async () => {
      await session.step(117, "When user sets \"Split\" property of histogram viewer to \"SEX\"", () => setProperty(page, "Split", el("histogram viewer"), "SEX"));
      await session.step(118, "And user sets \"Show Distribution Lines\" property of histogram viewer to \"true\"", () => setProperty(page, "Show Distribution Lines", el("histogram viewer"), "true"));
      await session.step(119, "Then histogram viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("histogram viewer"), 500));
      await session.step(120, "And histogram viewer should have more ink than before", () => moreInk(page, el("histogram viewer")));
      await session.step(121, "When user sets \"Show Markers\" property of histogram viewer to \"false\"", () => setProperty(page, "Show Markers", el("histogram viewer"), "false"));
      await session.step(122, "Then histogram viewer should have less ink than before", () => lessInk(page, el("histogram viewer")));
      await session.step(123, "When user sets \"Spline Tension\" property of histogram viewer to \"5\"", () => setProperty(page, "Spline Tension", el("histogram viewer"), "5"));
      await session.step(124, "Then histogram viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("histogram viewer"), 500));
      await session.step(125, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Spline Tension","0"],["Show Markers","true"],["Show Distribution Lines","false"],["Split",""]]), [["Spline Tension","0"],["Show Markers","true"],["Show Distribution Lines","false"],["Split",""]]);
      await session.step(130, "Then histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
      await session.step(131, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A range narrowed under a split keeps the table filter valid", async () => {
      await session.step(134, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Split","RACE"],["Show Range Inputs","true"]]), [["Split","RACE"],["Show Range Inputs","true"]]);
      await session.step(137, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(138, "When user enters \"30\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "30", "range min input", el("histogram viewer")));
      await session.step(139, "Then 861 rows should pass the filter", () => filterPasses(page, 861));
      await session.step(140, "And the filter should pass exactly the rows where \"AGE\" is between 30 and 89", () => filterIsExactly(page, "AGE", 30, 89));
      await session.step(141, "And histogram viewer should show 861 rows", () => showsRows(page, el("histogram viewer"), 861));
      await session.step(142, "And the legend of histogram viewer should list 4 items", () => legendLists(page, el("histogram viewer"), 4));
      await session.step(143, "And histogram viewer should have a \"line Caucasian\" area", () => hasArea(page, el("histogram viewer"), "line Caucasian"));
      await session.step(144, "When user enters \"18\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "18", "range min input", el("histogram viewer")));
      await session.step(145, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(146, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Show Range Inputs","false"],["Split",""]]), [["Show Range Inputs","false"],["Split",""]]);
      await session.step(149, "Then histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
      await session.step(150, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A split with a single bin can be hovered without errors (core#998)", async () => {
      await session.step(154, "When user sets \"Split\" property of histogram viewer to \"SEX\"", () => setProperty(page, "Split", el("histogram viewer"), "SEX"));
      await session.step(155, "Then histogram viewer should have a \"line F\" area", () => hasArea(page, el("histogram viewer"), "line F"));
      await session.step(156, "When user sets \"Bins\" property of histogram viewer to \"1\"", () => setProperty(page, "Bins", el("histogram viewer"), "1"));
      await session.step(157, "Then the \"bins shown\" reading of histogram viewer should be 1", () => readingIs(page, "bins shown", el("histogram viewer"), 1));
      await session.step(158, "And histogram viewer should not have a \"line F\" area", () => hasNoArea(page, el("histogram viewer"), "line F"));
      await session.step(159, "When user hovers over the \"view\" area of histogram viewer", () => hoverArea(page, "view", el("histogram viewer")));
      await session.step(160, "And user moves the pointer away from histogram viewer", () => pointerAway(page, el("histogram viewer")));
      await session.step(161, "Then the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
      await session.step(162, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Bins","20"],["Split",""]]), [["Bins","20"],["Split",""]]);
      await session.step(165, "Then the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
      await session.step(166, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Split lines take the split column's category colours (public#130)", async () => {
      await session.step(170, "When user colors \"SEX\" column categorically:", () => colorCategorical(page, "SEX", [["F","#00FF00"],["M","#FF00FF"]]), [["F","#00FF00"],["M","#FF00FF"]]);
      await session.step(173, "And user sets \"Split\" property of histogram viewer to \"SEX\"", () => setProperty(page, "Split", el("histogram viewer"), "SEX"));
      await session.step(174, "Then histogram viewer should have a \"line F\" area", () => hasArea(page, el("histogram viewer"), "line F"));
      await session.step(175, "And the \"line F\" area of histogram viewer should contain the color \"#00FF00\"", () => areaColor(page, "line F", el("histogram viewer"), "#00FF00"));
      await session.step(176, "And the \"line M\" area of histogram viewer should contain the color \"#FF00FF\"", () => areaColor(page, "line M", el("histogram viewer"), "#FF00FF"));
      await session.step(177, "When user sets \"Split\" property of histogram viewer to \"\"", () => setProperty(page, "Split", el("histogram viewer"), ""));
      await session.step(178, "And user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(179, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
