/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/histogram/histogram-range-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.histogram]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {rangeInputReads} from '../../../bindings/histogram.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactly, filterPasses, filterPassesAll, filterPassesFewer} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areasDiffer, doubleClickArea, dragAreaToArea, enterIntoArea, hasArea, hasNoArea, hoverArea, loadLayout, noErrors, propertyShouldBe, readingHigher, readingIs, readingLower, readingReads, repaintedBy, resizeTo, saveLayoutToServer, setProperties, setProperty, showsFewerRows, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Histogram range filter", () => {
  const session = feature(test, "features/viewers/histogram/histogram-range-filter.feature", import.meta.url);
  test("Histogram range filter", {tag: ["@journey", "@viewers", "@realizes:viewers.histogram"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]));
    await session.step(17, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(18, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(19, "And histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(20, "And the \"range min\" reading of histogram viewer should be 18", () => readingIs(page, "range min", el("histogram viewer"), 18));
    await session.step(21, "And the \"range max\" reading of histogram viewer should be 89", () => readingIs(page, "range max", el("histogram viewer"), 89));
    await session.step(22, "And the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
    await session.step(23, "And histogram viewer should have a \"range min input\" area", () => hasArea(page, el("histogram viewer"), "range min input"));
    await session.step(24, "And histogram viewer should have a \"range max input\" area", () => hasArea(page, el("histogram viewer"), "range max input"));
    await session.step(25, "And histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
    await run.scenario("Typed bounds filter the table and zoom the axis to them", async () => {
      await session.step(28, "When user enters \"30\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "30", "range min input", el("histogram viewer")));
      await session.step(29, "Then 861 rows should pass the filter", () => filterPasses(page, 861));
      await session.step(30, "And the \"range min\" reading of histogram viewer should be 30", () => readingIs(page, "range min", el("histogram viewer"), 30));
      await session.step(31, "And the \"axis min\" reading of histogram viewer should be higher than before", () => readingHigher(page, "axis min", el("histogram viewer")));
      await session.step(32, "And the \"bins shown\" reading of histogram viewer should be lower than before", () => readingLower(page, "bins shown", el("histogram viewer")));
      await session.step(33, "When user enters \"60\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "60", "range max input", el("histogram viewer")));
      await session.step(34, "Then 708 rows should pass the filter", () => filterPasses(page, 708));
      await session.step(35, "And the filter should pass exactly the rows where \"AGE\" is between 30 and 60", () => filterIsExactly(page, "AGE", 30, 60));
      await session.step(36, "And histogram viewer should show 708 rows", () => showsRows(page, el("histogram viewer"), 708));
      await session.step(37, "And the \"range max\" reading of histogram viewer should be 60", () => readingIs(page, "range max", el("histogram viewer"), 60));
      await session.step(38, "And the \"bins shown\" reading of histogram viewer should be lower than before", () => readingLower(page, "bins shown", el("histogram viewer")));
      await session.step(39, "And histogram viewer should have a \"range bar\" area", () => hasArea(page, el("histogram viewer"), "range bar"));
      await session.step(40, "When user enters \"18\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "18", "range min input", el("histogram viewer")));
      await session.step(41, "And user enters \"89\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "89", "range max input", el("histogram viewer")));
      await session.step(42, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(43, "And the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
      await session.step(44, "And the \"range min\" reading of histogram viewer should be 18", () => readingIs(page, "range min", el("histogram viewer"), 18));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A stacked split keeps the range filter", async () => {
      await session.step(48, "When user enters \"30\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "30", "range min input", el("histogram viewer")));
      await session.step(49, "And user enters \"60\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "60", "range max input", el("histogram viewer")));
      await session.step(50, "Then 708 rows should pass the filter", () => filterPasses(page, 708));
      await session.step(51, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Split","SEX"],["Split Stack","true"]]));
      await session.step(54, "Then 708 rows should pass the filter", () => filterPasses(page, 708));
      await session.step(55, "And histogram viewer should show 708 rows", () => showsRows(page, el("histogram viewer"), 708));
      await session.step(56, "And histogram viewer should have a \"bin 8 | F\" area", () => hasArea(page, el("histogram viewer"), "bin 8 | F"));
      await session.step(57, "And histogram viewer should have a \"bin 8 | M\" area", () => hasArea(page, el("histogram viewer"), "bin 8 | M"));
      await session.step(58, "And the \"bin 8 | F\" and \"bin 8 | M\" areas of histogram viewer should be painted in different colors", () => areasDiffer(page, "bin 8 | F", "bin 8 | M", el("histogram viewer")));
      await session.step(59, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Split",""],["Split Stack","false"]]));
      await session.step(62, "And user enters \"18\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "18", "range min input", el("histogram viewer")));
      await session.step(63, "And user enters \"89\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "89", "range max input", el("histogram viewer")));
      await session.step(64, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(65, "And histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Zoom To Range decides whether the axis follows the range", async () => {
      await session.step(69, "When user enters \"60\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "60", "range min input", el("histogram viewer")));
      await session.step(70, "Then 171 rows should pass the filter", () => filterPasses(page, 171));
      await session.step(71, "And the \"axis min\" reading of histogram viewer should be higher than before", () => readingHigher(page, "axis min", el("histogram viewer")));
      await session.step(72, "And the \"bins shown\" reading of histogram viewer should be lower than before", () => readingLower(page, "bins shown", el("histogram viewer")));
      await session.step(73, "When user sets \"Zoom To Range\" property of histogram viewer to \"false\"", () => setProperty(page, "Zoom To Range", el("histogram viewer"), "false"));
      await session.step(74, "Then the \"axis min\" reading of histogram viewer should be 18", () => readingIs(page, "axis min", el("histogram viewer"), 18));
      await session.step(75, "And the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
      await session.step(76, "And histogram viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("histogram viewer"), 500));
      await session.step(77, "And 171 rows should pass the filter", () => filterPasses(page, 171));
      await session.step(78, "When user sets \"Zoom To Range\" property of histogram viewer to \"true\"", () => setProperty(page, "Zoom To Range", el("histogram viewer"), "true"));
      await session.step(79, "Then the \"axis min\" reading of histogram viewer should be higher than before", () => readingHigher(page, "axis min", el("histogram viewer")));
      await session.step(80, "And the \"bins shown\" reading of histogram viewer should be lower than before", () => readingLower(page, "bins shown", el("histogram viewer")));
      await session.step(81, "When user enters \"18\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "18", "range min input", el("histogram viewer")));
      await session.step(82, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(83, "And the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Normalize To Filter scales the bars to the bins the range leaves", async () => {
      await session.step(87, "When user enters \"60\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "60", "range min input", el("histogram viewer")));
      await session.step(88, "Then 171 rows should pass the filter", () => filterPasses(page, 171));
      await session.step(89, "And the \"y axis max\" reading of histogram viewer should be 68", () => readingIs(page, "y axis max", el("histogram viewer"), 68));
      await session.step(90, "When user sets \"Normalize To Filter\" property of histogram viewer to \"false\"", () => setProperty(page, "Normalize To Filter", el("histogram viewer"), "false"));
      await session.step(91, "Then the \"y axis max\" reading of histogram viewer should be 99", () => readingIs(page, "y axis max", el("histogram viewer"), 99));
      await session.step(92, "And histogram viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("histogram viewer"), 500));
      await session.step(93, "When user sets \"Normalize To Filter\" property of histogram viewer to \"true\"", () => setProperty(page, "Normalize To Filter", el("histogram viewer"), "true"));
      await session.step(94, "Then the \"y axis max\" reading of histogram viewer should be 68", () => readingIs(page, "y axis max", el("histogram viewer"), 68));
      await session.step(95, "When user enters \"18\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "18", "range min input", el("histogram viewer")));
      await session.step(96, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(97, "And the \"y axis max\" reading of histogram viewer should be 99", () => readingIs(page, "y axis max", el("histogram viewer"), 99));
      await session.step(98, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A maximum below the minimum collapses the filter and reports it", async () => {
      await session.step(101, "When user enters \"40\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "40", "range min input", el("histogram viewer")));
      await session.step(102, "Then 665 rows should pass the filter", () => filterPasses(page, 665));
      await session.step(103, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
      await session.step(104, "When user enters \"20\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "20", "range max input", el("histogram viewer")));
      await session.step(105, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(106, "And the \"error\" reading of histogram viewer should be \"max should be greater than min\"", () => readingReads(page, "error", el("histogram viewer"), "max should be greater than min"));
      await session.step(107, "And histogram viewer should not have a \"bin 8\" area", () => hasNoArea(page, el("histogram viewer"), "bin 8"));
      await session.step(108, "When user enters \"60\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "60", "range max input", el("histogram viewer")));
      await session.step(109, "Then 512 rows should pass the filter", () => filterPasses(page, 512));
      await session.step(110, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
      await session.step(111, "And histogram viewer should have a \"bin 8\" area", () => hasArea(page, el("histogram viewer"), "bin 8"));
      await session.step(112, "When user enters \"18\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "18", "range min input", el("histogram viewer")));
      await session.step(113, "And user enters \"89\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "89", "range max input", el("histogram viewer")));
      await session.step(114, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bounds outside the column extent are clamped, the typed text is kept", async () => {
      await session.step(118, "When user enters \"40\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "40", "range min input", el("histogram viewer")));
      await session.step(119, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(120, "And the \"range min\" reading of histogram viewer should be 40", () => readingIs(page, "range min", el("histogram viewer"), 40));
      await session.step(121, "When user enters \"-999\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "-999", "range min input", el("histogram viewer")));
      await session.step(122, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(123, "And the \"range min\" reading of histogram viewer should be 18", () => readingIs(page, "range min", el("histogram viewer"), 18));
      await session.step(124, "And the range min input of histogram viewer should read \"-999\"", () => rangeInputReads(page, "min", "-999"));
      await session.step(125, "When user enters \"60\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "60", "range max input", el("histogram viewer")));
      await session.step(126, "Then 847 rows should pass the filter", () => filterPasses(page, 847));
      await session.step(127, "And the \"range min\" reading of histogram viewer should be 18", () => readingIs(page, "range min", el("histogram viewer"), 18));
      await session.step(128, "When user enters \"999\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "999", "range max input", el("histogram viewer")));
      await session.step(129, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(130, "And the \"range max\" reading of histogram viewer should be 89", () => readingIs(page, "range max", el("histogram viewer"), 89));
      await session.step(131, "And the range max input of histogram viewer should read \"999\"", () => rangeInputReads(page, "max", "999"));
      await session.step(132, "When user enters \"18\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "18", "range min input", el("histogram viewer")));
      await session.step(133, "And user enters \"89\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "89", "range max input", el("histogram viewer")));
      await session.step(134, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(135, "And the range min input of histogram viewer should read \"18\"", () => rangeInputReads(page, "min", "18"));
      await session.step(136, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A handle dragged with the pointer moves the range, a double click resets it", async () => {
      await session.step(139, "When user hovers over the \"view\" area of histogram viewer", () => hoverArea(page, "view", el("histogram viewer")));
      await session.step(140, "Then histogram viewer should have a \"range min handle\" area", () => hasArea(page, el("histogram viewer"), "range min handle"));
      await session.step(141, "And histogram viewer should have a \"range max handle\" area", () => hasArea(page, el("histogram viewer"), "range max handle"));
      await session.step(142, "When user drags the \"range max handle\" area of histogram viewer to the \"bin 12\" area", () => dragAreaToArea(page, "range max handle", el("histogram viewer"), "bin 12"));
      await session.step(143, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(144, "And the \"range max\" reading of histogram viewer should be lower than before", () => readingLower(page, "range max", el("histogram viewer")));
      await session.step(145, "And histogram viewer should show fewer rows than before", () => showsFewerRows(page, el("histogram viewer")));
      await session.step(146, "And histogram viewer should have a \"range bar\" area", () => hasArea(page, el("histogram viewer"), "range bar"));
      await session.step(147, "When user double-clicks on the \"range slider\" area of histogram viewer", () => doubleClickArea(page, "range slider", el("histogram viewer")));
      await session.step(148, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(149, "And the \"range max\" reading of histogram viewer should be 89", () => readingIs(page, "range max", el("histogram viewer"), 89));
      await session.step(150, "And the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
      await session.step(151, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Filtering Enabled off releases the table and keeps the zoom", async () => {
      await session.step(154, "When user enters \"40\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "40", "range min input", el("histogram viewer")));
      await session.step(155, "Then 665 rows should pass the filter", () => filterPasses(page, 665));
      await session.step(156, "And the \"bins shown\" reading of histogram viewer should be 14", () => readingIs(page, "bins shown", el("histogram viewer"), 14));
      await session.step(157, "When user sets \"Filtering Enabled\" property of histogram viewer to \"false\"", () => setProperty(page, "Filtering Enabled", el("histogram viewer"), "false"));
      await session.step(158, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(159, "And the \"range min\" reading of histogram viewer should be 40", () => readingIs(page, "range min", el("histogram viewer"), 40));
      await session.step(160, "And the \"bins shown\" reading of histogram viewer should be 14", () => readingIs(page, "bins shown", el("histogram viewer"), 14));
      await session.step(161, "When user sets \"Filtering Enabled\" property of histogram viewer to \"true\"", () => setProperty(page, "Filtering Enabled", el("histogram viewer"), "true"));
      await session.step(162, "Then 665 rows should pass the filter", () => filterPasses(page, 665));
      await session.step(163, "When user enters \"18\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "18", "range min input", el("histogram viewer")));
      await session.step(164, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(165, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the histogram back at the full range", async () => {
      await session.step(168, "When user enters \"40\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "40", "range min input", el("histogram viewer")));
      await session.step(169, "And user enters \"60\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "60", "range max input", el("histogram viewer")));
      await session.step(170, "Then 512 rows should pass the filter", () => filterPasses(page, 512));
      await session.step(171, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(172, "And user clicks on close icon of histogram viewer", () => clickOn(page, el("close icon of histogram viewer")));
      await session.step(173, "Then histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
      await session.step(174, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(175, "When user loads the saved layout", () => loadLayout(page));
      await session.step(176, "Then histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
      await session.step(177, "And \"Value\" property of histogram viewer should be \"AGE\"", () => propertyShouldBe(page, "Value", el("histogram viewer"), "AGE"));
      await session.step(178, "And histogram viewer should have a \"range min input\" area", () => hasArea(page, el("histogram viewer"), "range min input"));
      await session.step(179, "And the \"range min\" reading of histogram viewer should be 18", () => readingIs(page, "range min", el("histogram viewer"), 18));
      await session.step(180, "And the \"range max\" reading of histogram viewer should be 89", () => readingIs(page, "range max", el("histogram viewer"), 89));
      await session.step(181, "And the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
      await session.step(182, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(183, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
