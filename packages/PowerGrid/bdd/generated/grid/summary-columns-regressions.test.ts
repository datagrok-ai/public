/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/grid/summary-columns-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hasArea, hasNoArea, noErrors, pickFromAreaContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {checkedInColumnList, toggleInColumnList, uncheckedInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid summary columns — regressions of closed issues", () => {
  const session = feature(test, "features/grid/summary-columns-regressions.feature", import.meta.url);
  test("A summary column stays when another column is hidden and shown again (public#1093)", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "Then the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
    await session.step(19, "When user picks \"Add > Summary Columns > Sparklines\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Sparklines", "cell 2 of USUBJID", el("grid")));
    await session.step(20, "Then the \"cell type of Sparklines\" reading of grid should be \"sparkline\"", () => readingReads(page, "cell type of Sparklines", el("grid"), "sparkline"));
    await session.step(21, "When user picks \"Hide\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Hide", "header AGE", el("grid")));
    await session.step(22, "Then grid should not have a \"header AGE\" area", () => hasNoArea(page, el("grid"), "header AGE"));
    await session.step(23, "And the \"column order\" reading of grid should be \"USUBJID, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines\"", () => readingReads(page, "column order", el("grid"), "USUBJID, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines"));
    await session.step(24, "When user picks \"Order or Hide Columns...\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Order or Hide Columns...", "header SEX", el("grid")));
    await session.step(25, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
    await session.step(26, "When user toggles the \"AGE\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "AGE", el("Order or Hide Columns dialog")));
    await session.step(27, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
    await session.step(28, "Then grid should have a \"header AGE\" area", () => hasArea(page, el("grid"), "header AGE"));
    await session.step(29, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines"));
    await session.step(30, "And the \"cell type of Sparklines\" reading of grid should be \"sparkline\"", () => readingReads(page, "cell type of Sparklines", el("grid"), "sparkline"));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
  test("A summary column itself is hidden and shown again in Order or Hide Columns (public#1093)", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "Then the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
    await session.step(34, "When user picks \"Add > Summary Columns > Sparklines\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Sparklines", "cell 2 of USUBJID", el("grid")));
    await session.step(35, "Then grid should have a \"header Sparklines\" area", () => hasArea(page, el("grid"), "header Sparklines"));
    await session.step(36, "When user picks \"Order or Hide Columns...\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Order or Hide Columns...", "header SEX", el("grid")));
    await session.step(37, "Then Order or Hide Columns dialog should be visible", () => shouldBe(page, el("Order or Hide Columns dialog"), "visible"));
    await session.step(38, "And the \"Sparklines\" column should be checked in the column list of Order or Hide Columns dialog", () => checkedInColumnList(page, "Sparklines", el("Order or Hide Columns dialog")));
    await session.step(39, "When user toggles the \"Sparklines\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "Sparklines", el("Order or Hide Columns dialog")));
    await session.step(40, "Then the \"Sparklines\" column should not be checked in the column list of Order or Hide Columns dialog", () => uncheckedInColumnList(page, "Sparklines", el("Order or Hide Columns dialog")));
    await session.step(41, "And grid should not have a \"header Sparklines\" area", () => hasNoArea(page, el("grid"), "header Sparklines"));
    await session.step(42, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
    await session.step(43, "When user toggles the \"Sparklines\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "Sparklines", el("Order or Hide Columns dialog")));
    await session.step(44, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
    await session.step(45, "Then grid should have a \"header Sparklines\" area", () => hasArea(page, el("grid"), "header Sparklines"));
    await session.step(46, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, Sparklines"));
    await session.step(47, "And the \"cell type of Sparklines\" reading of grid should be \"sparkline\"", () => readingReads(page, "cell type of Sparklines", el("grid"), "sparkline"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
});
