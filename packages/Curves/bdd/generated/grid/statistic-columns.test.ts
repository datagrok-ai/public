/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/grid/statistic-columns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/curves.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, selectIn, shouldBe, shouldOfferChoice} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType, columnsExactly, distinctValues} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {currentColumnIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, closeAllViews, contextPanelOpen, openDataset, packageInstalled, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Fit statistics of a curve cell become columns", () => {
  const session = feature(test, "features/grid/statistic-columns.feature", import.meta.url);
  test("Fit statistics of a curve cell become columns", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(20, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("The IC50 of a one-curve column is added right after it, filled with numbers", async () => {
      await session.step(24, "When user scrolls the grid to the \"styled proprtional with IC50\" column", () => scrollGridTo(page, "styled proprtional with IC50"));
      await session.step(25, "And user clicks on the \"cell 1 of styled proprtional with IC50\" area of grid", () => clickArea(page, "cell 1 of styled proprtional with IC50", el("grid")));
      await session.step(26, "Then the current column should be \"styled proprtional with IC50\"", () => currentColumnIs(page, "styled proprtional with IC50"));
      await session.step(27, "Given \"Fit\" pane in context panel is expanded", () => isExpanded(page, el("\"Fit\" pane in context panel")));
      await session.step(28, "Then \"IC50\" input in context panel should be visible", () => shouldBe(page, el("\"IC50\" input in context panel"), "visible"));
      await session.step(29, "And \"R²\" input in context panel should be visible", () => shouldBe(page, el("\"R²\" input in context panel"), "visible"));
      await session.step(30, "When user clicks on first \"Calculate ic50 for the whole column\" icon in context panel", () => clickOn(page, el("first \"Calculate ic50 for the whole column\" icon in context panel")));
      await session.step(31, "Then the table should have the columns \"smiles, multiple prefit, multiple styled series, styled proprtional with IC50, styled proprtional with IC50 1 ic50, log min/max candlesticks, custom fit function\"", () => columnsExactly(page, "smiles, multiple prefit, multiple styled series, styled proprtional with IC50, styled proprtional with IC50 1 ic50, log min/max candlesticks, custom fit function"));
      await session.step(32, "And \"styled proprtional with IC50 1 ic50\" column should have type \"double\"", () => columnType(page, "styled proprtional with IC50 1 ic50", "double"));
      await session.step(33, "And \"styled proprtional with IC50 1 ic50\" column should have at least 10 distinct values", () => distinctValues(page, "styled proprtional with IC50 1 ic50", 10));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The average IC50 of a two-curve column is added right after it", async () => {
      await session.step(37, "Given user closes all views", () => closeAllViews(page));
      await session.step(38, "And user opens curves dataset", () => openDataset(page, ds("curves")));
      await session.step(39, "And the context panel is open", () => contextPanelOpen(page));
      await session.step(40, "When user scrolls the grid to the \"multiple styled series\" column", () => scrollGridTo(page, "multiple styled series"));
      await session.step(41, "And user clicks on the \"cell 2 of multiple styled series\" area of grid", () => clickArea(page, "cell 2 of multiple styled series", el("grid")));
      await session.step(42, "Then the current column should be \"multiple styled series\"", () => currentColumnIs(page, "multiple styled series"));
      await session.step(43, "Given \"Fit\" pane in context panel is expanded", () => isExpanded(page, el("\"Fit\" pane in context panel")));
      await session.step(44, "Then \"Series\" input in context panel should offer the choice \"all\"", () => shouldOfferChoice(page, el("\"Series\" input in context panel"), "all"));
      await session.step(45, "And \"Series\" input in context panel should offer the choice \"aggregated\"", () => shouldOfferChoice(page, el("\"Series\" input in context panel"), "aggregated"));
      await session.step(46, "When user selects \"aggregated\" in \"Series\" input in context panel", () => selectIn(page, "aggregated", el("\"Series\" input in context panel")));
      await session.step(47, "And user selects \"avg\" in \"Aggregation type\" input in context panel", () => selectIn(page, "avg", el("\"Aggregation type\" input in context panel")));
      await session.step(48, "And user clicks on first \"Calculate ic50 avg for the whole column\" icon in context panel", () => clickOn(page, el("first \"Calculate ic50 avg for the whole column\" icon in context panel")));
      await session.step(49, "Then the table should have the columns \"smiles, multiple prefit, multiple styled series, multiple styled series avg ic50, styled proprtional with IC50, log min/max candlesticks, custom fit function\"", () => columnsExactly(page, "smiles, multiple prefit, multiple styled series, multiple styled series avg ic50, styled proprtional with IC50, log min/max candlesticks, custom fit function"));
      await session.step(50, "And \"multiple styled series avg ic50\" column should have type \"double\"", () => columnType(page, "multiple styled series avg ic50", "double"));
      await session.step(51, "And \"multiple styled series avg ic50\" column should have at least 10 distinct values", () => distinctValues(page, "multiple styled series avg ic50", 10));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
