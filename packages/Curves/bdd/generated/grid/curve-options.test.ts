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
import {isExpanded, selectIn, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentColumnIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, closeAllViews, contextPanelOpen, layoutsDeleted, noProjectOnServer, openDataset, openProjectWithTable, packageInstalled, saveAsProject, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {applyNamedLayout, clickArea, noErrors, saveNamedLayout} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Curve options and the level they are set at", () => {
  const session = feature(test, "features/grid/curve-options.feature", import.meta.url);
  test("Curve options and the level they are set at", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(20, "And no project named \"CurvesOptions-{run}\" is on the server", () => noProjectOnServer(page, session.text("CurvesOptions-{run}")));
    await session.step(21, "And the layouts named \"curves-options-{run}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("curves-options-{run}")));
    await session.step(22, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(23, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("A Column level fit function reaches the other cells of its column", async () => {
      await session.step(26, "When user scrolls the grid to the \"styled proprtional with IC50\" column", () => scrollGridTo(page, "styled proprtional with IC50"));
      await session.step(27, "And user clicks on the \"cell 1 of styled proprtional with IC50\" area of grid", () => clickArea(page, "cell 1 of styled proprtional with IC50", el("grid")));
      await session.step(28, "Then the current column should be \"styled proprtional with IC50\"", () => currentColumnIs(page, "styled proprtional with IC50"));
      await session.step(29, "Given \"Options\" pane in context panel is expanded", () => isExpanded(page, el("\"Options\" pane in context panel")));
      await session.step(30, "Then \"Level\" input in context panel should have value \"Column\"", () => shouldHaveValue(page, el("\"Level\" input in context panel"), "Column"));
      await session.step(31, "And \"Fit Function\" input in context panel should have value \"sigmoid\"", () => shouldHaveValue(page, el("\"Fit Function\" input in context panel"), "sigmoid"));
      await session.step(32, "When user selects \"linear\" in \"Fit Function\" input in context panel", () => selectIn(page, "linear", el("\"Fit Function\" input in context panel")));
      await session.step(33, "And user clicks on the \"cell 5 of styled proprtional with IC50\" area of grid", () => clickArea(page, "cell 5 of styled proprtional with IC50", el("grid")));
      await session.step(34, "Then the current column should be \"styled proprtional with IC50\"", () => currentColumnIs(page, "styled proprtional with IC50"));
      await session.step(35, "And \"Fit Function\" input in context panel should have value \"linear\"", () => shouldHaveValue(page, el("\"Fit Function\" input in context panel"), "linear"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Column level fit function is kept by a layout applied to the table opened again", async () => {
      await session.step(39, "When user saves the layout of the current table view as \"curves-options-{run}\"", () => saveNamedLayout(page, session.text("curves-options-{run}")));
      await session.step(40, "And user closes all views", () => closeAllViews(page));
      await session.step(41, "Given user opens curves dataset", () => openDataset(page, ds("curves")));
      await session.step(42, "And the context panel is open", () => contextPanelOpen(page));
      await session.step(43, "When user applies the layout \"curves-options-{run}\" to the current table view", () => applyNamedLayout(page, session.text("curves-options-{run}")));
      await session.step(44, "And user scrolls the grid to the \"styled proprtional with IC50\" column", () => scrollGridTo(page, "styled proprtional with IC50"));
      await session.step(45, "And user clicks on the \"cell 4 of styled proprtional with IC50\" area of grid", () => clickArea(page, "cell 4 of styled proprtional with IC50", el("grid")));
      await session.step(46, "Then the current column should be \"styled proprtional with IC50\"", () => currentColumnIs(page, "styled proprtional with IC50"));
      await session.step(47, "Given \"Options\" pane in context panel is expanded", () => isExpanded(page, el("\"Options\" pane in context panel")));
      await session.step(48, "Then \"Fit Function\" input in context panel should have value \"linear\"", () => shouldHaveValue(page, el("\"Fit Function\" input in context panel"), "linear"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Column level fit function is kept by a project", async () => {
      await session.step(52, "When user saves the current view as project \"CurvesOptions-{run}\"", () => saveAsProject(page, session.text("CurvesOptions-{run}")));
      await session.step(53, "And user closes all views", () => closeAllViews(page));
      await session.step(54, "And user opens the \"CurvesOptions-{run}\" project and waits for its table", () => openProjectWithTable(page, session.text("CurvesOptions-{run}")));
      await session.step(55, "And user scrolls the grid to the \"styled proprtional with IC50\" column", () => scrollGridTo(page, "styled proprtional with IC50"));
      await session.step(56, "And user clicks on the \"cell 4 of styled proprtional with IC50\" area of grid", () => clickArea(page, "cell 4 of styled proprtional with IC50", el("grid")));
      await session.step(57, "Then the current column should be \"styled proprtional with IC50\"", () => currentColumnIs(page, "styled proprtional with IC50"));
      await session.step(58, "Given \"Options\" pane in context panel is expanded", () => isExpanded(page, el("\"Options\" pane in context panel")));
      await session.step(59, "Then \"Fit Function\" input in context panel should have value \"linear\"", () => shouldHaveValue(page, el("\"Fit Function\" input in context panel"), "linear"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
