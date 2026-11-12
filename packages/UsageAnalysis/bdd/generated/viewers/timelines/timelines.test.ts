/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/timelines/timelines.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.timelines]
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
import {clickOn, doubleClickOn, hoverOver, isExpanded, selectIn, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, clickArea, clickLegendItem, clickLegendItemHolding, legendLists, noErrors, painted, pickFromContextMenu, propertyShouldBe, readingAsRemembered, readingIs, readingLowerThanRemembered, rememberReading, repainted, takeSnapshot, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {doesNotOfferColumn, offersColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Timelines legend, legend visibility, legend clicks, the table filter and Reset View", () => {
  const session = feature(test, "features/viewers/timelines/timelines.feature", import.meta.url);
  test("Timelines legend, legend visibility, legend clicks, the table filter and Reset View", {tag: ["@journey", "@viewers", "@realizes:charts.viewer.timelines"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(15, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(16, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(17, "And Files---App-Data---Charts tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Charts tree node inside browse tree")));
    await session.step(18, "When user double-clicks on Files---App-Data---Charts---ae.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Charts---ae.csv tree node inside browse tree")));
    await session.step(19, "Then the \"ae\" view should be current", () => viewIsCurrent(page, "ae"));
    await session.step(20, "And the table should have 143 rows", () => rowCount(page, 143));
    await session.step(21, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(22, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(23, "When user clicks on first \"Timelines\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Timelines\" card in \"Add Viewer\" dialog")));
    await session.step(24, "Then timelines viewer should be visible", () => shouldBe(page, el("timelines viewer"), "visible"));
    await session.step(25, "And timelines viewer should be bound to table \"ae\"", () => boundTable(page, el("timelines viewer"), "ae"));
    await session.step(26, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(27, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(28, "And user clicks on settings icon of timelines viewer", () => clickOn(page, el("settings icon of timelines viewer")));
    await session.step(29, "Given \"Color\" category in context panel is expanded", () => isExpanded(page, el("\"Color\" category in context panel")));
    await session.step(30, "Then \"Color\" property in context panel should be visible", () => shouldBe(page, el("\"Color\" property in context panel"), "visible"));
    await run.scenario("The viewer draws, and the legend follows Legend Visibility (GROK-20800)", async () => {
      await session.step(33, "Then no errors should have been logged", () => noErrors(page));
      await session.step(34, "And the \"intervals\" reading of timelines viewer should be 143", () => readingIs(page, "intervals", el("timelines viewer"), 143));
      await session.step(35, "And \"Color\" property in context panel should offer the columns \"STUDYID, DOMAIN, USUBJID, AESPID, AETERM, AELLT, AELLTCD, AEDECOD, AEPTCD, AEHLT, AEHLTCD, AEHLGT, AEHLGTCD, AEBODSYS, AEBDSYCD, AESOC, AESOCCD, AESEV, AEACN, AEREL, AEOUT, AESTDTC\"", () => offersColumns(page, el("\"Color\" property in context panel"), "STUDYID, DOMAIN, USUBJID, AESPID, AETERM, AELLT, AELLTCD, AEDECOD, AEPTCD, AEHLT, AEHLTCD, AEHLGT, AEHLGTCD, AEBODSYS, AEBDSYCD, AESOC, AESOCCD, AESEV, AEACN, AEREL, AEOUT, AESTDTC"));
      await session.step(36, "And \"Color\" property in context panel should not offer the column \"AESTDY\"", () => doesNotOfferColumn(page, el("\"Color\" property in context panel"), "AESTDY"));
      await session.step(37, "And \"Color\" property in context panel should not offer the column \"AESEQ\"", () => doesNotOfferColumn(page, el("\"Color\" property in context panel"), "AESEQ"));
      await session.step(38, "When user selects \"AESOC\" in \"Color\" property in context panel", () => selectIn(page, "AESOC", el("\"Color\" property in context panel")));
      await session.step(39, "Then \"Color\" property of timelines viewer should be \"AESOC\"", () => propertyShouldBe(page, "Color", el("timelines viewer"), "AESOC"));
      await session.step(40, "And the legend of timelines viewer should list 15 items", () => legendLists(page, el("timelines viewer"), 15));
      await session.step(41, "And legend of timelines viewer should be visible", () => shouldBe(page, el("legend of timelines viewer"), "visible"));
      await session.step(42, "Given \"Legend\" category in context panel is expanded", () => isExpanded(page, el("\"Legend\" category in context panel")));
      await session.step(43, "When user selects \"Always\" in \"Legend Visibility\" property in context panel", () => selectIn(page, "Always", el("\"Legend Visibility\" property in context panel")));
      await session.step(44, "Then legend of timelines viewer should be visible", () => shouldBe(page, el("legend of timelines viewer"), "visible"));
      await session.step(45, "And the legend of timelines viewer should list 15 items", () => legendLists(page, el("timelines viewer"), 15));
      await session.step(46, "When user selects \"Never\" in \"Legend Visibility\" property in context panel", () => selectIn(page, "Never", el("\"Legend Visibility\" property in context panel")));
      await session.step(47, "Then legend of timelines viewer should be hidden", () => shouldBe(page, el("legend of timelines viewer"), "hidden"));
      await session.step(48, "When user selects \"Auto\" in \"Legend Visibility\" property in context panel", () => selectIn(page, "Auto", el("\"Legend Visibility\" property in context panel")));
      await session.step(49, "Then legend of timelines viewer should be visible", () => shouldBe(page, el("legend of timelines viewer"), "visible"));
      await session.step(50, "And the legend of timelines viewer should list 15 items", () => legendLists(page, el("timelines viewer"), 15));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend clicks narrow the viewer without blanking it, and so does a filter card (GROK-19033, GROK-19535, GROK-18608)", async () => {
      await session.step(54, "When user selects \"AESOC\" in \"Color\" property in context panel", () => selectIn(page, "AESOC", el("\"Color\" property in context panel")));
      await session.step(55, "Then the legend of timelines viewer should list 15 items", () => legendLists(page, el("timelines viewer"), 15));
      await session.step(56, "And the \"intervals\" reading of timelines viewer should be 143", () => readingIs(page, "intervals", el("timelines viewer"), 143));
      await session.step(57, "When user clicks on \"SKIN AND SUBCUTANEOUS TISSUE DISORDERS\" item in the legend of timelines viewer", () => clickLegendItem(page, "SKIN AND SUBCUTANEOUS TISSUE DISORDERS", el("timelines viewer")));
      await session.step(58, "Then the \"intervals\" reading of timelines viewer should be 35", () => readingIs(page, "intervals", el("timelines viewer"), 35));
      await session.step(59, "And 143 rows should pass the filter", () => filterPasses(page, 143));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
      await session.step(61, "When user clicks on \"NERVOUS SYSTEM DISORDERS\" item in the legend of timelines viewer holding Control", () => clickLegendItemHolding(page, "NERVOUS SYSTEM DISORDERS", el("timelines viewer"), "Control"));
      await session.step(62, "Then the \"intervals\" reading of timelines viewer should be 68", () => readingIs(page, "intervals", el("timelines viewer"), 68));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
      await session.step(64, "When user clicks on \"NERVOUS SYSTEM DISORDERS\" item in the legend of timelines viewer holding Control", () => clickLegendItemHolding(page, "NERVOUS SYSTEM DISORDERS", el("timelines viewer"), "Control"));
      await session.step(65, "Then the \"intervals\" reading of timelines viewer should be 35", () => readingIs(page, "intervals", el("timelines viewer"), 35));
      await session.step(66, "When user clicks on \"SKIN AND SUBCUTANEOUS TISSUE DISORDERS\" item in the legend of timelines viewer", () => clickLegendItem(page, "SKIN AND SUBCUTANEOUS TISSUE DISORDERS", el("timelines viewer")));
      await session.step(67, "Then the \"intervals\" reading of timelines viewer should be 143", () => readingIs(page, "intervals", el("timelines viewer"), 143));
      await session.step(68, "And \"SKIN AND SUBCUTANEOUS TISSUE DISORDERS\" legend item in legend of timelines viewer should not be selected", () => shouldNotBe(page, el("\"SKIN AND SUBCUTANEOUS TISSUE DISORDERS\" legend item in legend of timelines viewer"), "selected"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
      await session.step(70, "Given \"Data\" category in context panel is expanded", () => isExpanded(page, el("\"Data\" category in context panel")));
      await session.step(71, "When user takes a snapshot of timelines viewer", () => takeSnapshot(page, el("timelines viewer")));
      await session.step(72, "And user selects \"AESEV\" in \"Split By\" property in context panel", () => selectIn(page, "AESEV", el("\"Split By\" property in context panel")));
      await session.step(73, "Then \"Split By\" property of timelines viewer should be \"AESEV\"", () => propertyShouldBe(page, "Split By", el("timelines viewer"), "AESEV"));
      await session.step(74, "And timelines viewer should have repainted", () => repainted(page, el("timelines viewer")));
      await session.step(75, "And the legend of timelines viewer should list 15 items", () => legendLists(page, el("timelines viewer"), 15));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
      await session.step(77, "When user clicks on \"CARDIAC DISORDERS\" item in the legend of timelines viewer", () => clickLegendItem(page, "CARDIAC DISORDERS", el("timelines viewer")));
      await session.step(78, "Then \"CARDIAC DISORDERS\" legend item in legend of timelines viewer should be selected", () => shouldBe(page, el("\"CARDIAC DISORDERS\" legend item in legend of timelines viewer"), "selected"));
      await session.step(79, "And the \"intervals\" reading of timelines viewer should be 11", () => readingIs(page, "intervals", el("timelines viewer"), 11));
      await session.step(80, "And 143 rows should pass the filter", () => filterPasses(page, 143));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
      await session.step(82, "When user takes a snapshot of timelines viewer", () => takeSnapshot(page, el("timelines viewer")));
      await session.step(83, "And user selects \"USUBJID\" in \"Split By\" property in context panel", () => selectIn(page, "USUBJID", el("\"Split By\" property in context panel")));
      await session.step(84, "Then \"Split By\" property of timelines viewer should be \"USUBJID\"", () => propertyShouldBe(page, "Split By", el("timelines viewer"), "USUBJID"));
      await session.step(85, "And timelines viewer should have repainted", () => repainted(page, el("timelines viewer")));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
      await session.step(87, "When user clicks on \"CARDIAC DISORDERS\" item in the legend of timelines viewer", () => clickLegendItem(page, "CARDIAC DISORDERS", el("timelines viewer")));
      await session.step(88, "Then the \"intervals\" reading of timelines viewer should be 143", () => readingIs(page, "intervals", el("timelines viewer"), 143));
      await session.step(89, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(90, "And user clicks on the \"category SEVERE of AESEV\" area of filter panel", () => clickArea(page, "category SEVERE of AESEV", el("filter panel")));
      await session.step(91, "Then 7 rows should pass the filter", () => filterPasses(page, 7));
      await session.step(92, "And the \"intervals\" reading of timelines viewer should be 7", () => readingIs(page, "intervals", el("timelines viewer"), 7));
      await session.step(93, "And the \"lanes\" reading of timelines viewer should be 4", () => readingIs(page, "lanes", el("timelines viewer"), 4));
      await session.step(94, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(95, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(96, "Then the \"intervals\" reading of timelines viewer should be 143", () => readingIs(page, "intervals", el("timelines viewer"), 143));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Reset View brings back every interval after a zoom", async () => {
      await session.step(100, "When user remembers the \"intervals\" reading of timelines viewer", () => rememberReading(page, "intervals", el("timelines viewer")));
      await session.step(101, "And user scrolls the mouse wheel up 5 times over the \"view\" area of timelines viewer", () => wheelOverAreaTimes(page, "up", 5, "view", el("timelines viewer")));
      await session.step(102, "Then the \"intervals\" reading of timelines viewer should be lower than remembered", () => readingLowerThanRemembered(page, "intervals", el("timelines viewer")));
      await session.step(103, "When user picks \"Reset View\" from the context menu of timelines viewer", () => pickFromContextMenu(page, "Reset View", el("timelines viewer")));
      await session.step(104, "Then the \"intervals\" reading of timelines viewer should be as remembered", () => readingAsRemembered(page, "intervals", el("timelines viewer")));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
