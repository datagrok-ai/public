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
import {columnWidth} from '../../bindings/curves.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, expand, scrollIntoView, selectIn, shouldBe, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, setProperties} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Show Legend of curves in grid cells is set for the table, a column or one cell", () => {
  const session = feature(test, "features/guides/curves-legend-cells.feature", import.meta.url);
  test("Hide the legend for the column, then show it in one cell", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens assay-curves dataset", () => openDataset(page, ds("assay-curves")));
    await session.step(18, "And user removes \"Molecule\" column", () => removeColumn(page, "Molecule"));
    await session.step(19, "And user sets properties of grid:", () => setProperties(page, el("grid"), [["Row Height","290"],["Col Header Height","30"]]), [["Row Height","290"],["Col Header Height","30"]]);
    await session.step(22, "And the \"Fitted Curve\" column of the grid is 340 px wide", () => columnWidth(page, "Fitted Curve", 340));
    await session.step(23, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(24, "When user clicks on the \"cell 1 of Fitted Curve\" area of grid", () => clickArea(page, "cell 1 of Fitted Curve", el("grid")));
    await session.step(25, "And user expands \"Options\" pane in context panel", () => expand(page, el("\"Options\" pane in context panel")));
    await session.step(26, "Then \"Level\" input in context panel should have value \"Column\"", () => shouldHaveValue(page, el("\"Level\" input in context panel"), "Column"));
    await session.step(27, "When user scrolls \"Show Legend\" input in context panel into view", () => scrollIntoView(page, el("\"Show Legend\" input in context panel")));
    await session.step(28, "And user unchecks \"Show Legend\" input in context panel", () => uncheck(page, el("\"Show Legend\" input in context panel")));
    await session.step(29, "And user scrolls \"Level\" input in context panel into view", () => scrollIntoView(page, el("\"Level\" input in context panel")));
    await session.step(30, "And user selects \"Cell\" in \"Level\" input in context panel", () => selectIn(page, "Cell", el("\"Level\" input in context panel")));
    await session.step(31, "And user scrolls \"Show Legend\" input in context panel into view", () => scrollIntoView(page, el("\"Show Legend\" input in context panel")));
    await session.step(32, "And user checks \"Show Legend\" input in context panel", () => check(page, el("\"Show Legend\" input in context panel")));
    await session.step(33, "And user clicks on the \"cell 2 of Fitted Curve\" area of grid", () => clickArea(page, "cell 2 of Fitted Curve", el("grid")));
    await session.step(34, "Then \"Show Legend\" input in context panel should be unchecked", () => shouldBe(page, el("\"Show Legend\" input in context panel"), "unchecked"));
    await session.step(35, "When user clicks on the \"cell 1 of Fitted Curve\" area of grid", () => clickArea(page, "cell 1 of Fitted Curve", el("grid")));
    await session.step(36, "Then \"Show Legend\" input in context panel should be checked", () => shouldBe(page, el("\"Show Legend\" input in context panel"), "checked"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
});
