/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/curves-legend-cells.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {everyCellLegend, hideGridColumn, noCellLegend, onlyCurrentCellLegend, sizeCurveCells, wheelIntoView} from '../../bindings/curves.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, expand, selectIn, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelOpen, openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Show Legend of curves in grid cells is set for the table, a column or one cell", () => {
  const session = feature(test, "features/guides/curves-legend-cells.feature", import.meta.url);
  test("Hide the legend for the column, then show it in one cell", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And user opens assay-curves dataset", () => openDataset(page, ds("assay-curves")));
    await session.step(14, "And the grid of the current view does not show the \"Molecule\" column", () => hideGridColumn(page, "Molecule"));
    await session.step(15, "And the grid of the current view has rows 290 px high and the \"Fitted Curve\" column 340 px wide", () => sizeCurveCells(page, 290, "Fitted Curve", 340));
    await session.step(16, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(17, "When user clicks on the \"cell 1 of Fitted Curve\" area of grid", () => clickArea(page, "cell 1 of Fitted Curve", el("grid")));
    await session.step(18, "And user expands \"Options\" pane in context panel", () => expand(page, el("\"Options\" pane in context panel")));
    await session.step(19, "Then \"Level\" input in context panel should have value \"Column\"", () => shouldHaveValue(page, el("\"Level\" input in context panel"), "Column"));
    await session.step(20, "And every curve cell in view of \"Fitted Curve\" column should draw a legend", () => everyCellLegend(page, "Fitted Curve"));
    await session.step(21, "When user scrolls \"Show Legend\" input in context panel into view", () => wheelIntoView(page, el("\"Show Legend\" input in context panel")));
    await session.step(22, "And user unchecks \"Show Legend\" input in context panel", () => uncheck(page, el("\"Show Legend\" input in context panel")));
    await session.step(23, "Then no curve cell in view of \"Fitted Curve\" column should draw a legend", () => noCellLegend(page, "Fitted Curve"));
    await session.step(24, "When user scrolls \"Level\" input in context panel into view", () => wheelIntoView(page, el("\"Level\" input in context panel")));
    await session.step(25, "And user selects \"Cell\" in \"Level\" input in context panel", () => selectIn(page, "Cell", el("\"Level\" input in context panel")));
    await session.step(26, "And user scrolls \"Show Legend\" input in context panel into view", () => wheelIntoView(page, el("\"Show Legend\" input in context panel")));
    await session.step(27, "And user checks \"Show Legend\" input in context panel", () => check(page, el("\"Show Legend\" input in context panel")));
    await session.step(28, "Then only the current curve cell of \"Fitted Curve\" column should draw a legend", () => onlyCurrentCellLegend(page, "Fitted Curve"));
  });
});
