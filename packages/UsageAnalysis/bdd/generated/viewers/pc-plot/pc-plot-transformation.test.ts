/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot-transformation.feature
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
import {axesShouldBe} from '../../../bindings/pc-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addRangeFilter, filterPasses, filterPassesAll, openFilterPanel, rowCount, selectFirstRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, errorBalloonText, hasArea, hasNoArea, noBalloons, noErrors, painted, readingAsRemembered, readingIs, readingNotAsRemembered, readingReads, rememberReading, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot transformations", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-transformation.feature", import.meta.url);
  test("PC plot transformations", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(29, "And user adds a pc plot viewer", () => addViewer(page, "pc plot"));
    await session.step(30, "Then the \"transformed\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "transformed", el("pc plot viewer"), "false"));
    await session.step(31, "And the \"axes\" reading of pc plot viewer should be 10", () => readingIs(page, "axes", el("pc plot viewer"), 10));
    await session.step(32, "And pc plot viewer should show 100 rows", () => showsRows(page, el("pc plot viewer"), 100));
    await session.step(33, "And the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
    await session.step(34, "And user remembers the \"axis order\" reading of pc plot viewer", () => rememberReading(page, "axis order", el("pc plot viewer")));
    await run.scenario("A pivot replaces the axes with the pivot column's categories", async () => {
      await session.step(37, "When user sets \"Transformation\" property of pc plot viewer to '[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"Chemist 521\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"pivot\",\"colName\":\"Series\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"count\",\"colName\":\"Id\"}]'", () => setProperty(page, "Transformation", el("pc plot viewer"), "[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"Chemist 521\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"pivot\",\"colName\":\"Series\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"count\",\"colName\":\"Id\"}]"));
      await session.step(38, "Then the \"transformed\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "transformed", el("pc plot viewer"), "true"));
      await session.step(39, "And the \"axes\" reading of pc plot viewer should be 5", () => readingIs(page, "axes", el("pc plot viewer"), 5));
      await session.step(40, "And the axes of pc plot viewer should be \"Triazoles, Diazabicyclooctane, Pyrrolidines, , Aminopiperidines\"", () => axesShouldBe(page, el("pc plot viewer"), "Triazoles, Diazabicyclooctane, Pyrrolidines, , Aminopiperidines"));
      await session.step(41, "And pc plot viewer should have an \"axis \\\"Triazoles\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis \"Triazoles\""));
      await session.step(42, "And the \"rows\" reading of pc plot viewer should be 13", () => readingIs(page, "rows", el("pc plot viewer"), 13));
      await session.step(43, "And pc plot viewer should show 13 rows", () => showsRows(page, el("pc plot viewer"), 13));
      await session.step(44, "And the \"axis order\" reading of pc plot viewer should not be as remembered", () => readingNotAsRemembered(page, "axis order", el("pc plot viewer")));
      await session.step(45, "And pc plot viewer should be painted", () => painted(page, el("pc plot viewer")));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clearing the transformation puts the original axes back", async () => {
      await session.step(49, "When user sets \"Transformation\" property of pc plot viewer to \"\"", () => setProperty(page, "Transformation", el("pc plot viewer"), ""));
      await session.step(50, "Then the \"transformed\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "transformed", el("pc plot viewer"), "false"));
      await session.step(51, "And the \"axes\" reading of pc plot viewer should be 10", () => readingIs(page, "axes", el("pc plot viewer"), 10));
      await session.step(52, "And the \"axis order\" reading of pc plot viewer should be as remembered", () => readingAsRemembered(page, "axis order", el("pc plot viewer")));
      await session.step(53, "And pc plot viewer should not have an \"axis \\\"Triazoles\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis \"Triazoles\""));
      await session.step(54, "And pc plot viewer should show 100 rows", () => showsRows(page, el("pc plot viewer"), 100));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An aggregation draws the aggregated columns", async () => {
      await session.step(58, "When user sets \"Transformation\" property of pc plot viewer to '[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"Series\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"Average Mass\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"TPSA\"}]'", () => setProperty(page, "Transformation", el("pc plot viewer"), "[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"Series\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"Average Mass\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"TPSA\"}]"));
      await session.step(59, "Then the \"transformed\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "transformed", el("pc plot viewer"), "true"));
      await session.step(60, "And the axes of pc plot viewer should be \"avg(Average Mass), avg(TPSA)\"", () => axesShouldBe(page, el("pc plot viewer"), "avg(Average Mass), avg(TPSA)"));
      await session.step(61, "And the \"rows\" reading of pc plot viewer should be 5", () => readingIs(page, "rows", el("pc plot viewer"), 5));
      await session.step(62, "And pc plot viewer should show 5 rows", () => showsRows(page, el("pc plot viewer"), 5));
      await session.step(63, "And pc plot viewer should be painted", () => painted(page, el("pc plot viewer")));
      await session.step(64, "When user sets \"Transformation\" property of pc plot viewer to \"\"", () => setProperty(page, "Transformation", el("pc plot viewer"), ""));
      await session.step(65, "Then the \"transformed\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "transformed", el("pc plot viewer"), "false"));
      await session.step(66, "And the \"axis order\" reading of pc plot viewer should be as remembered", () => readingAsRemembered(page, "axis order", el("pc plot viewer")));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A transformation that does not parse raises a balloon and leaves the axes alone", async () => {
      await session.step(70, "When user sets \"Transformation\" property of pc plot viewer to \"not a transformation\"", () => setProperty(page, "Transformation", el("pc plot viewer"), "not a transformation"));
      await session.step(71, "Then an error balloon containing \"Invalid transformation\" should have been shown", () => errorBalloonText(page, "Invalid transformation"));
      await session.step(72, "And the \"error\" reading of pc plot viewer should be \"Invalid transformation\"", () => readingReads(page, "error", el("pc plot viewer"), "Invalid transformation"));
      await session.step(73, "And the \"transformed\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "transformed", el("pc plot viewer"), "false"));
      await session.step(74, "And the \"axis order\" reading of pc plot viewer should be as remembered", () => readingAsRemembered(page, "axis order", el("pc plot viewer")));
      await session.step(75, "And pc plot viewer should show 100 rows", () => showsRows(page, el("pc plot viewer"), 100));
      await session.step(76, "And pc plot viewer should be painted", () => painted(page, el("pc plot viewer")));
      await session.step(77, "When user sets \"Transformation\" property of pc plot viewer to '[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"Series\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"Average Mass\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"TPSA\"}]'", () => setProperty(page, "Transformation", el("pc plot viewer"), "[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"Series\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"Average Mass\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"TPSA\"}]"));
      await session.step(78, "Then the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
      await session.step(79, "And the \"transformed\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "transformed", el("pc plot viewer"), "true"));
      await session.step(80, "When user sets \"Transformation\" property of pc plot viewer to \"\"", () => setProperty(page, "Transformation", el("pc plot viewer"), ""));
      await session.step(81, "Then the \"axis order\" reading of pc plot viewer should be as remembered", () => readingAsRemembered(page, "axis order", el("pc plot viewer")));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing the filter panel with a transformation applied breaks nothing (GROK-18091)", async () => {
      await session.step(85, "Given the package autostarts have completed", () => autostartsCompleted(page));
      await session.step(86, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(87, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(88, "When user sets \"Transformation\" property of pc plot viewer to '[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"Series\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"Average Mass\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"TPSA\"}]'", () => setProperty(page, "Transformation", el("pc plot viewer"), "[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"Series\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"Average Mass\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"TPSA\"}]"));
      await session.step(89, "Then the \"transformed\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "transformed", el("pc plot viewer"), "true"));
      await session.step(90, "And pc plot viewer should show 5 rows", () => showsRows(page, el("pc plot viewer"), 5));
      await session.step(91, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(92, "Then filter panel should be absent", () => shouldBe(page, el("filter panel"), "absent"));
      await session.step(93, "And pc plot viewer should be visible", () => shouldBe(page, el("pc plot viewer"), "visible"));
      await session.step(94, "And pc plot viewer should be painted", () => painted(page, el("pc plot viewer")));
      await session.step(95, "And the \"transformed\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "transformed", el("pc plot viewer"), "true"));
      await session.step(96, "And the \"rows\" reading of pc plot viewer should be 5", () => readingIs(page, "rows", el("pc plot viewer"), 5));
      await session.step(97, "And the table should have 100 rows", () => rowCount(page, 100));
      await session.step(98, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(99, "When user sets \"Transformation\" property of pc plot viewer to \"\"", () => setProperty(page, "Transformation", el("pc plot viewer"), ""));
      await session.step(100, "Then the \"axis order\" reading of pc plot viewer should be as remembered", () => readingAsRemembered(page, "axis order", el("pc plot viewer")));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With a transformation, Reset filters restores the rows but drops the selection (GROK-17306)", async () => {
      await session.step(105, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(106, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]));
      await session.step(108, "Then pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(109, "When user sets \"Transformation\" property of pc plot viewer to '[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"SEX\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"pivot\",\"colName\":\"DIS_POP\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"WEIGHT\"}]'", () => setProperty(page, "Transformation", el("pc plot viewer"), "[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"SEX\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"pivot\",\"colName\":\"DIS_POP\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"WEIGHT\"}]"));
      await session.step(110, "Then the \"transformed\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "transformed", el("pc plot viewer"), "true"));
      await session.step(111, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(112, "And user adds a range filter on \"AGE\" from 30 to 50", () => addRangeFilter(page, "AGE", 30, 50));
      await session.step(113, "Then 494 rows should pass the filter", () => filterPasses(page, 494));
      await session.step(114, "When user selects the first 10 rows", () => selectFirstRows(page, 10));
      await session.step(115, "Then 10 rows should be selected", () => selectedRowCount(page, 10));
      await session.step(116, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(117, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(118, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(119, "And 10 rows should be selected", () => selectedRowCount(page, 10));
    }, {knownFailure: true});
    run.finish();
  });
});
