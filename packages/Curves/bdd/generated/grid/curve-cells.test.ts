/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/grid/curve-cells.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/curves.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe, shouldNotOfferChoice, shouldOfferChoice} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {currentColumnIs, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelOpen, openDataset, packageInstalled, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridPins, scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Curve cells in the grid and their Context Panel", () => {
  const session = feature(test, "features/grid/curve-cells.feature", import.meta.url);
  test("Curve cells in the grid and their Context Panel", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(17, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(18, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("The curve columns are detected as fitted curves", async () => {
      await session.step(21, "Then the table should have 30 rows", () => rowCount(page, 30));
      await session.step(22, "And \"multiple prefit\" column should have semantic type \"fit\"", () => columnSemType(page, "multiple prefit", "fit"));
      await session.step(23, "And \"multiple styled series\" column should have semantic type \"fit\"", () => columnSemType(page, "multiple styled series", "fit"));
      await session.step(24, "And \"styled proprtional with IC50\" column should have semantic type \"fit\"", () => columnSemType(page, "styled proprtional with IC50", "fit"));
      await session.step(25, "And \"log min/max candlesticks\" column should have semantic type \"fit\"", () => columnSemType(page, "log min/max candlesticks", "fit"));
      await session.step(26, "And \"custom fit function\" column should have semantic type \"fit\"", () => columnSemType(page, "custom fit function", "fit"));
      await session.step(27, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell with two curves shows its panes, the fit functions and the summary inputs", async () => {
      await session.step(30, "When user scrolls the grid to the \"multiple styled series\" column", () => scrollGridTo(page, "multiple styled series"));
      await session.step(31, "And user clicks on the \"cell 2 of multiple styled series\" area of grid", () => clickArea(page, "cell 2 of multiple styled series", el("grid")));
      await session.step(32, "Then the current column should be \"multiple styled series\"", () => currentColumnIs(page, "multiple styled series"));
      await session.step(33, "And \"Options\" pane in context panel should be visible", () => shouldBe(page, el("\"Options\" pane in context panel"), "visible"));
      await session.step(34, "And \"Chart\" pane in context panel should be visible", () => shouldBe(page, el("\"Chart\" pane in context panel"), "visible"));
      await session.step(35, "And \"Fit\" pane in context panel should be visible", () => shouldBe(page, el("\"Fit\" pane in context panel"), "visible"));
      await session.step(36, "Given \"Options\" pane in context panel is expanded", () => isExpanded(page, el("\"Options\" pane in context panel")));
      await session.step(37, "Then \"Level\" input in context panel should have value \"Column\"", () => shouldHaveValue(page, el("\"Level\" input in context panel"), "Column"));
      await session.step(38, "And \"Level\" input in context panel should offer the choice \"Dataframe\"", () => shouldOfferChoice(page, el("\"Level\" input in context panel"), "Dataframe"));
      await session.step(39, "And \"Level\" input in context panel should offer the choice \"Cell\"", () => shouldOfferChoice(page, el("\"Level\" input in context panel"), "Cell"));
      await session.step(40, "And \"Fit Function\" input in context panel should offer the choice \"linear\"", () => shouldOfferChoice(page, el("\"Fit Function\" input in context panel"), "linear"));
      await session.step(41, "And \"Fit Function\" input in context panel should offer the choice \"sigmoid\"", () => shouldOfferChoice(page, el("\"Fit Function\" input in context panel"), "sigmoid"));
      await session.step(42, "And \"Fit Function\" input in context panel should offer the choice \"log-linear\"", () => shouldOfferChoice(page, el("\"Fit Function\" input in context panel"), "log-linear"));
      await session.step(43, "And \"Fit Function\" input in context panel should offer the choice \"exponential\"", () => shouldOfferChoice(page, el("\"Fit Function\" input in context panel"), "exponential"));
      await session.step(44, "And \"Fit Function\" input in context panel should offer the choice \"4pl-regression\"", () => shouldOfferChoice(page, el("\"Fit Function\" input in context panel"), "4pl-regression"));
      await session.step(45, "And \"Fit Function\" input in context panel should offer the choice \"4pl-dose-response\"", () => shouldOfferChoice(page, el("\"Fit Function\" input in context panel"), "4pl-dose-response"));
      await session.step(46, "And \"Fit Function\" input in context panel should not offer the choice \"\"", () => shouldNotOfferChoice(page, el("\"Fit Function\" input in context panel"), ""));
      await session.step(47, "And \"Droplines\" input in context panel should contain text \"IC10\"", () => shouldContainText(page, el("\"Droplines\" input in context panel"), "IC10"));
      await session.step(48, "And \"Droplines\" input in context panel should contain text \"IC90\"", () => shouldContainText(page, el("\"Droplines\" input in context panel"), "IC90"));
      await session.step(49, "And \"Statistics mode\" input in context panel should be visible", () => shouldBe(page, el("\"Statistics mode\" input in context panel"), "visible"));
      await session.step(50, "And \"Aggregation\" input in context panel should be visible", () => shouldBe(page, el("\"Aggregation\" input in context panel"), "visible"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell with one curve has nothing to summarise", async () => {
      await session.step(54, "When user scrolls the grid to the \"styled proprtional with IC50\" column", () => scrollGridTo(page, "styled proprtional with IC50"));
      await session.step(55, "And user clicks on the \"cell 1 of styled proprtional with IC50\" area of grid", () => clickArea(page, "cell 1 of styled proprtional with IC50", el("grid")));
      await session.step(56, "Then the current column should be \"styled proprtional with IC50\"", () => currentColumnIs(page, "styled proprtional with IC50"));
      await session.step(57, "Given \"Options\" pane in context panel is expanded", () => isExpanded(page, el("\"Options\" pane in context panel")));
      await session.step(58, "Then \"Fit Function\" input in context panel should be visible", () => shouldBe(page, el("\"Fit Function\" input in context panel"), "visible"));
      await session.step(59, "And \"Statistics mode\" input in context panel should not be visible", () => shouldNotBe(page, el("\"Statistics mode\" input in context panel"), "visible"));
      await session.step(60, "And \"Aggregation\" input in context panel should not be visible", () => shouldNotBe(page, el("\"Aggregation\" input in context panel"), "visible"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A double-click opens the chart in the Edit chart dialog, and Close closes it", async () => {
      await session.step(64, "When user scrolls the grid to the \"multiple prefit\" column", () => scrollGridTo(page, "multiple prefit"));
      await session.step(65, "And user double-clicks on the \"cell 3 of multiple prefit\" area of grid", () => doubleClickArea(page, "cell 3 of multiple prefit", el("grid")));
      await session.step(66, "Then \"Edit chart\" dialog should be visible", () => shouldBe(page, el("\"Edit chart\" dialog"), "visible"));
      await session.step(67, "When user clicks on \"Close\" button in \"Edit chart\" dialog", () => clickOn(page, el("\"Close\" button in \"Edit chart\" dialog")));
      await session.step(68, "Then \"Edit chart\" dialog should be absent", () => shouldBe(page, el("\"Edit chart\" dialog"), "absent"));
      await session.step(69, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The filter panel and a pinned curve column take the curves without errors", async () => {
      await session.step(73, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(74, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
      await session.step(76, "When user scrolls the grid to the \"multiple prefit\" column", () => scrollGridTo(page, "multiple prefit"));
      await session.step(77, "And user picks \"Pin > Pin Column\" from the context menu of the \"header multiple prefit\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Column", "header multiple prefit", el("grid")));
      await session.step(78, "Then the grid should pin the columns \"multiple prefit\"", () => gridPins(page, "multiple prefit"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
