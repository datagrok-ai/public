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
import {clickOn, isExpanded, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentColumnIs, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelOpen, openDataset, packageInstalled, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridPins, scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Curve cells in the grid and their Context Panel", () => {
  const session = feature(test, "features/grid/curve-cells.feature", import.meta.url);
  test("Curve cells in the grid and their Context Panel", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(21, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("A cell with two curves shows its panes and the level its options are set at", async () => {
      await session.step(25, "When user scrolls the grid to the \"multiple styled series\" column", () => scrollGridTo(page, "multiple styled series"));
      await session.step(26, "And user clicks on the \"cell 2 of multiple styled series\" area of grid", () => clickArea(page, "cell 2 of multiple styled series", el("grid")));
      await session.step(27, "Then the current column should be \"multiple styled series\"", () => currentColumnIs(page, "multiple styled series"));
      await session.step(28, "And \"Options\" pane in context panel should be visible", () => shouldBe(page, el("\"Options\" pane in context panel"), "visible"));
      await session.step(29, "And \"Chart\" pane in context panel should be visible", () => shouldBe(page, el("\"Chart\" pane in context panel"), "visible"));
      await session.step(30, "And \"Fit\" pane in context panel should be visible", () => shouldBe(page, el("\"Fit\" pane in context panel"), "visible"));
      await session.step(31, "Given \"Options\" pane in context panel is expanded", () => isExpanded(page, el("\"Options\" pane in context panel")));
      await session.step(32, "Then \"Level\" input in context panel should have value \"Column\"", () => shouldHaveValue(page, el("\"Level\" input in context panel"), "Column"));
      await session.step(33, "And \"Statistics mode\" input in context panel should be visible", () => shouldBe(page, el("\"Statistics mode\" input in context panel"), "visible"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A double-click opens the chart in the Edit chart dialog, and Close closes it", async () => {
      await session.step(37, "When user scrolls the grid to the \"multiple prefit\" column", () => scrollGridTo(page, "multiple prefit"));
      await session.step(38, "And user double-clicks on the \"cell 3 of multiple prefit\" area of grid", () => doubleClickArea(page, "cell 3 of multiple prefit", el("grid")));
      await session.step(39, "Then \"Edit chart\" dialog should be visible", () => shouldBe(page, el("\"Edit chart\" dialog"), "visible"));
      await session.step(40, "When user clicks on \"Close\" button in \"Edit chart\" dialog", () => clickOn(page, el("\"Close\" button in \"Edit chart\" dialog")));
      await session.step(41, "Then \"Edit chart\" dialog should be absent", () => shouldBe(page, el("\"Edit chart\" dialog"), "absent"));
      await session.step(42, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The filter panel and a pinned curve column take the curves without errors", async () => {
      await session.step(46, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(47, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
      await session.step(49, "When user scrolls the grid to the \"multiple prefit\" column", () => scrollGridTo(page, "multiple prefit"));
      await session.step(50, "And user picks \"Pin > Pin Column\" from the context menu of the \"header multiple prefit\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Column", "header multiple prefit", el("grid")));
      await session.step(51, "Then the grid should pin the columns \"multiple prefit\"", () => gridPins(page, "multiple prefit"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
