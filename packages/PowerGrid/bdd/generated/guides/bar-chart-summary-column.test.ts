/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/bar-chart-summary-column.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Chart several properties in one grid column", () => {
  const session = feature(test, "features/guides/bar-chart-summary-column.feature", import.meta.url);
  test("Add a bar chart summary column and choose the columns it draws", {tag: ["@guide", "@help:visualize/viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(17, "When user picks \"Add > Summary Columns > Bar Chart\" from the context menu of the \"cell 2 of Id\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Bar Chart", "cell 2 of Id", el("grid")));
    await session.step(18, "And user clicks on the \"header Bar Chart\" area of grid", () => clickArea(page, "header Bar Chart", el("grid")));
    await session.step(19, "And user clicks on Columns input in context panel", () => clickOn(page, el("Columns input in context panel")));
    await session.step(20, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(21, "When user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
    await session.step(22, "And user toggles the \"Average Mass\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Average Mass", el("\"Select columns...\" dialog")));
    await session.step(23, "And user toggles the \"TPSA\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "TPSA", el("\"Select columns...\" dialog")));
    await session.step(24, "And user toggles the \"Num Rotatable Bonds\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "Num Rotatable Bonds", el("\"Select columns...\" dialog")));
    await session.step(25, "And user toggles the \"NIBR logP\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "NIBR logP", el("\"Select columns...\" dialog")));
    await session.step(26, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(27, "Then Columns input in context panel should contain text \"(4) Average Mass, TPSA, Num Rotatable Bonds, NIBR logP\"", () => shouldContainText(page, el("Columns input in context panel"), "(4) Average Mass, TPSA, Num Rotatable Bonds, NIBR logP"));
    await session.step(28, "When user clicks on the \"cell 3 of Bar Chart\" area of grid", () => clickArea(page, "cell 3 of Bar Chart", el("grid")));
    await session.step(29, "Then context panel should contain text \"TPSA: 92.5\"", () => shouldContainText(page, el("context panel"), "TPSA: 92.5"));
  });
});
