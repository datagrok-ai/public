/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/grid/curve-options.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/curves.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {enterInto, isExpanded, selectIn, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentColumnIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelOpen, openDataset, packageInstalled, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Curve options and the level they are set at", () => {
  const session = feature(test, "features/grid/curve-options.feature", import.meta.url);
  test("Curve options and the level they are set at", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(20, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("A Column level fit function reaches the other cells of its column, not another column", async () => {
      await session.step(24, "When user scrolls the grid to the \"styled proprtional with IC50\" column", () => scrollGridTo(page, "styled proprtional with IC50"));
      await session.step(25, "And user clicks on the \"cell 1 of styled proprtional with IC50\" area of grid", () => clickArea(page, "cell 1 of styled proprtional with IC50", el("grid")));
      await session.step(26, "Then the current column should be \"styled proprtional with IC50\"", () => currentColumnIs(page, "styled proprtional with IC50"));
      await session.step(27, "Given \"Options\" pane in context panel is expanded", () => isExpanded(page, el("\"Options\" pane in context panel")));
      await session.step(28, "Then \"Level\" input in context panel should have value \"Column\"", () => shouldHaveValue(page, el("\"Level\" input in context panel"), "Column"));
      await session.step(29, "And \"Fit Function\" input in context panel should have value \"sigmoid\"", () => shouldHaveValue(page, el("\"Fit Function\" input in context panel"), "sigmoid"));
      await session.step(30, "When user selects \"linear\" in \"Fit Function\" input in context panel", () => selectIn(page, "linear", el("\"Fit Function\" input in context panel")));
      await session.step(31, "And user clicks on the \"cell 5 of styled proprtional with IC50\" area of grid", () => clickArea(page, "cell 5 of styled proprtional with IC50", el("grid")));
      await session.step(32, "Then the current column should be \"styled proprtional with IC50\"", () => currentColumnIs(page, "styled proprtional with IC50"));
      await session.step(33, "And \"Fit Function\" input in context panel should have value \"linear\"", () => shouldHaveValue(page, el("\"Fit Function\" input in context panel"), "linear"));
      await session.step(34, "When user scrolls the grid to the \"custom fit function\" column", () => scrollGridTo(page, "custom fit function"));
      await session.step(35, "And user clicks on the \"cell 1 of custom fit function\" area of grid", () => clickArea(page, "cell 1 of custom fit function", el("grid")));
      await session.step(36, "Then the current column should be \"custom fit function\"", () => currentColumnIs(page, "custom fit function"));
      await session.step(37, "And \"Fit Function\" input in context panel should have value \"Polynomial\"", () => shouldHaveValue(page, el("\"Fit Function\" input in context panel"), "Polynomial"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Cell level title stays in its cell", async () => {
      await session.step(41, "When user scrolls the grid to the \"styled proprtional with IC50\" column", () => scrollGridTo(page, "styled proprtional with IC50"));
      await session.step(42, "And user clicks on the \"cell 1 of styled proprtional with IC50\" area of grid", () => clickArea(page, "cell 1 of styled proprtional with IC50", el("grid")));
      await session.step(43, "Then the current column should be \"styled proprtional with IC50\"", () => currentColumnIs(page, "styled proprtional with IC50"));
      await session.step(44, "And \"Title\" input in context panel should have value \"Dose-Response curves\"", () => shouldHaveValue(page, el("\"Title\" input in context panel"), "Dose-Response curves"));
      await session.step(45, "When user selects \"Cell\" in \"Level\" input in context panel", () => selectIn(page, "Cell", el("\"Level\" input in context panel")));
      await session.step(46, "And user enters \"Only this one\" into \"Title\" input in context panel", () => enterInto(page, "Only this one", el("\"Title\" input in context panel")));
      await session.step(47, "And user clicks on the \"cell 3 of styled proprtional with IC50\" area of grid", () => clickArea(page, "cell 3 of styled proprtional with IC50", el("grid")));
      await session.step(48, "Then the current column should be \"styled proprtional with IC50\"", () => currentColumnIs(page, "styled proprtional with IC50"));
      await session.step(49, "And \"Title\" input in context panel should have value \"Dose-Response curves\"", () => shouldHaveValue(page, el("\"Title\" input in context panel"), "Dose-Response curves"));
      await session.step(50, "When user clicks on the \"cell 1 of styled proprtional with IC50\" area of grid", () => clickArea(page, "cell 1 of styled proprtional with IC50", el("grid")));
      await session.step(51, "Then \"Title\" input in context panel should have value \"Only this one\"", () => shouldHaveValue(page, el("\"Title\" input in context panel"), "Only this one"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
