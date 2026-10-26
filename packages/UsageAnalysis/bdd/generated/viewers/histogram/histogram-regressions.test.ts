/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/histogram/histogram-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.histogram, core#998, public#2301, public#130]
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
import {addCalculated, clearSelection, colorCategorical, colorOff, filterPasses, filterPassesAll, noneOfSelected, removeColumn, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, clickArea, doubleClickArea, dragAreaBy, hasArea, hasNoArea, hoverArea, noErrors, pointerAway, readingIs, readingReads, resizeTo, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Histogram regressions from closed issues", () => {
  const session = feature(test, "features/viewers/histogram/histogram-regressions.feature", import.meta.url);
  test("Histogram regressions from closed issues", {tag: ["@journey", "@viewers", "@realizes:viewers.histogram", "@realizes:core#998", "@realizes:public#2301", "@realizes:public#130"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"]]), [["Value","AGE"]]);
    await session.step(19, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(20, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(21, "And the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
    await session.step(22, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
    await run.scenario("Under the viewer's own Filter a bin click selects only the rows the bar counts (core#998)", async () => {
      await session.step(27, "When user sets \"Filter\" property of histogram viewer to \"${SEX} == \\\"F\\\"\"", () => setProperty(page, "Filter", el("histogram viewer"), "${SEX} == \"F\""));
      await session.step(28, "Then histogram viewer should show 553 rows", () => showsRows(page, el("histogram viewer"), 553));
      await session.step(29, "When user clicks on the \"bin 8\" area of histogram viewer", () => clickArea(page, "bin 8", el("histogram viewer")));
      await session.step(30, "Then 53 rows should be selected", () => selectedRowCount(page, 53));
      await session.step(31, "And no rows where \"SEX\" is \"M\" should be selected", () => noneOfSelected(page, "SEX", "M"));
      await session.step(32, "When user clears the row selection", () => clearSelection(page));
      await session.step(33, "And user sets \"Filter\" property of histogram viewer to \"\"", () => setProperty(page, "Filter", el("histogram viewer"), ""));
      await session.step(34, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A split with a single bin can be hovered without errors (core#998)", async () => {
      await session.step(40, "When user sets \"Split\" property of histogram viewer to \"SEX\"", () => setProperty(page, "Split", el("histogram viewer"), "SEX"));
      await session.step(41, "Then histogram viewer should have a \"line F\" area", () => hasArea(page, el("histogram viewer"), "line F"));
      await session.step(42, "When user sets \"Bins\" property of histogram viewer to \"1\"", () => setProperty(page, "Bins", el("histogram viewer"), "1"));
      await session.step(43, "Then the \"bins shown\" reading of histogram viewer should be 1", () => readingIs(page, "bins shown", el("histogram viewer"), 1));
      await session.step(44, "And histogram viewer should not have a \"line F\" area", () => hasNoArea(page, el("histogram viewer"), "line F"));
      await session.step(45, "When user hovers over the \"view\" area of histogram viewer", () => hoverArea(page, "view", el("histogram viewer")));
      await session.step(46, "And user moves the pointer away from histogram viewer", () => pointerAway(page, el("histogram viewer")));
      await session.step(47, "Then the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
      await session.step(49, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Bins","20"],["Split",""]]), [["Bins","20"],["Split",""]]);
      await session.step(52, "Then the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
    });
    await run.scenario("A logarithmic X axis over a column with no positive value is reported, not thrown (core#998)", async () => {
      await session.step(57, "When user adds a calculated column \"NEG_AGE\" with formula \"0 - ${AGE}\"", () => addCalculated(page, "NEG_AGE", "0 - ${AGE}"));
      await session.step(58, "And user sets \"Value\" property of histogram viewer to \"NEG_AGE\"", () => setProperty(page, "Value", el("histogram viewer"), "NEG_AGE"));
      await session.step(59, "Then histogram viewer should have a \"bin 1\" area", () => hasArea(page, el("histogram viewer"), "bin 1"));
      await session.step(60, "When user sets \"X Axis Type\" property of histogram viewer to \"logarithmic\"", () => setProperty(page, "X Axis Type", el("histogram viewer"), "logarithmic"));
      await session.step(61, "Then the \"error\" reading of histogram viewer should be \"no positive values to show on the logarithmic X axis\"", () => readingReads(page, "error", el("histogram viewer"), "no positive values to show on the logarithmic X axis"));
      await session.step(62, "And histogram viewer should not have a \"bin 1\" area", () => hasNoArea(page, el("histogram viewer"), "bin 1"));
      await session.step(63, "When user hovers over the \"view\" area of histogram viewer", () => hoverArea(page, "view", el("histogram viewer")));
      await session.step(64, "Then no errors should have been logged", () => noErrors(page));
      await session.step(65, "When user sets \"X Axis Type\" property of histogram viewer to \"linear\"", () => setProperty(page, "X Axis Type", el("histogram viewer"), "linear"));
      await session.step(66, "Then the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
      await session.step(67, "And histogram viewer should have a \"bin 1\" area", () => hasArea(page, el("histogram viewer"), "bin 1"));
      await session.step(68, "When user sets \"Value\" property of histogram viewer to \"AGE\"", () => setProperty(page, "Value", el("histogram viewer"), "AGE"));
      await session.step(69, "And user removes \"NEG_AGE\" column", () => removeColumn(page, "NEG_AGE"));
      await session.step(70, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Y Min above the Y Max is reported, not thrown (core#998)", async () => {
      await session.step(75, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Y Min","80"],["Y Max","20"]]), [["Y Min","80"],["Y Max","20"]]);
      await session.step(78, "Then the \"error\" reading of histogram viewer should be \"Y Max should be at least 80\"", () => readingReads(page, "error", el("histogram viewer"), "Y Max should be at least 80"));
      await session.step(79, "And histogram viewer should not have a \"bin 8\" area", () => hasNoArea(page, el("histogram viewer"), "bin 8"));
      await session.step(80, "When user hovers over the \"view\" area of histogram viewer", () => hoverArea(page, "view", el("histogram viewer")));
      await session.step(81, "Then no errors should have been logged", () => noErrors(page));
      await session.step(82, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Y Min",""],["Y Max",""]]), [["Y Min",""],["Y Max",""]]);
      await session.step(85, "Then the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
      await session.step(86, "And histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The min handle dragged past the max handle leaves a valid range (public#2301)", async () => {
      await session.step(92, "When user hovers over the \"view\" area of histogram viewer", () => hoverArea(page, "view", el("histogram viewer")));
      await session.step(93, "Then histogram viewer should have a \"range min handle\" area", () => hasArea(page, el("histogram viewer"), "range min handle"));
      await session.step(94, "When user drags the \"range min handle\" area of histogram viewer by 450 pixels to the right", () => dragAreaBy(page, "range min handle", el("histogram viewer"), 450, "right"));
      await session.step(95, "Then the \"range min\" reading of histogram viewer should be 89", () => readingIs(page, "range min", el("histogram viewer"), 89));
      await session.step(96, "And the \"range max\" reading of histogram viewer should be 89", () => readingIs(page, "range max", el("histogram viewer"), 89));
      await session.step(97, "And 1 row should pass the filter", () => filterPasses(page, 1));
      await session.step(98, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
      await session.step(99, "And no errors should have been logged", () => noErrors(page));
      await session.step(100, "When user double-clicks on the \"range slider\" area of histogram viewer", () => doubleClickArea(page, "range slider", el("histogram viewer")));
      await session.step(101, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(102, "And the \"range min\" reading of histogram viewer should be 18", () => readingIs(page, "range min", el("histogram viewer"), 18));
    });
    await run.scenario("Split lines take the split column's category colours (public#130)", async () => {
      await session.step(107, "When user colors \"SEX\" column categorically:", () => colorCategorical(page, "SEX", [["F","#00FF00"],["M","#FF00FF"]]), [["F","#00FF00"],["M","#FF00FF"]]);
      await session.step(110, "And user sets \"Split\" property of histogram viewer to \"SEX\"", () => setProperty(page, "Split", el("histogram viewer"), "SEX"));
      await session.step(111, "Then histogram viewer should have a \"line F\" area", () => hasArea(page, el("histogram viewer"), "line F"));
      await session.step(112, "And the \"line F\" area of histogram viewer should contain the color \"#00FF00\"", () => areaColor(page, "line F", el("histogram viewer"), "#00FF00"));
      await session.step(113, "And the \"line M\" area of histogram viewer should contain the color \"#FF00FF\"", () => areaColor(page, "line M", el("histogram viewer"), "#FF00FF"));
      await session.step(114, "When user sets \"Split\" property of histogram viewer to \"\"", () => setProperty(page, "Split", el("histogram viewer"), ""));
      await session.step(115, "And user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(116, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
