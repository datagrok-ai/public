/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/control-comparisons.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.analyze.group-comparison.control-comparisons]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {tableColumnComplete, tableColumns, tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, noBalloons, noErrors, painted, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Control comparisons", () => {
  const session = feature(test, "features/analyze/control-comparisons.feature", import.meta.url);
  test("Running the default comparison docks a box plot and the comparison statistics", {tag: ["@eda", "@realizes:ml.menu.analyze.group-comparison.control-comparisons"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(15, "When user picks \"ML > Analyze > Group Comparison > Control Comparisons...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Group Comparison > Control Comparisons..."));
    await session.step(16, "Then \"Control comparisons\" dialog should be visible", () => shouldBe(page, el("\"Control comparisons\" dialog"), "visible"));
    await session.step(17, "And editor of Category input in \"Control comparisons\" dialog should have text \"RACE\"", () => shouldHaveText(page, el("editor of Category input in \"Control comparisons\" dialog"), "RACE"));
    await session.step(18, "And Control input in \"Control comparisons\" dialog should have value \"Asian\"", () => shouldHaveValue(page, el("Control input in \"Control comparisons\" dialog"), "Asian"));
    await session.step(19, "And editor of Feature input in \"Control comparisons\" dialog should have text \"AGE\"", () => shouldHaveText(page, el("editor of Feature input in \"Control comparisons\" dialog"), "AGE"));
    await session.step(20, "And Alpha input in \"Control comparisons\" dialog should have value \"0.05\"", () => shouldHaveValue(page, el("Alpha input in \"Control comparisons\" dialog"), "0.05"));
    await session.step(21, "And Run button in \"Control comparisons\" dialog should be enabled", () => shouldBe(page, el("Run button in \"Control comparisons\" dialog"), "enabled"));
    await session.step(22, "When user clicks on Run button in \"Control comparisons\" dialog", () => clickOn(page, el("Run button in \"Control comparisons\" dialog")));
    await session.step(23, "Then \"Control comparisons\" dialog should be hidden", () => shouldBe(page, el("\"Control comparisons\" dialog"), "hidden"));
    await session.step(24, "And box plot viewer should be visible", () => shouldBe(page, el("box plot viewer"), "visible"));
    await session.step(25, "And description of box plot viewer should be visible", () => shouldBe(page, el("description of box plot viewer"), "visible"));
    await session.step(26, "And description of box plot viewer should contain text \"Asian\"", () => shouldContainText(page, el("description of box plot viewer"), "Asian"));
    await session.step(27, "And box plot viewer should be painted", () => painted(page, el("box plot viewer")));
    await session.step(28, "And table \"Control comparisons result\" should be open", () => tableOpen(page, "Control comparisons result"));
    await session.step(29, "And table \"Control comparisons result\" should have 3 rows", () => tableRows(page, "Control comparisons result", 3));
    await session.step(30, "And table \"Control comparisons result\" should have columns \"Conclusion, Group, n, Mean diff, 95% CI low, 95% CI high, t, df, p (raw), p (adj), Hedges' g\"", () => tableColumns(page, "Control comparisons result", "Conclusion, Group, n, Mean diff, 95% CI low, 95% CI high, t, df, p (raw), p (adj), Hedges' g"));
    await session.step(31, "And table \"Control comparisons result\" should have no missing values in \"Mean diff\" column", () => tableColumnComplete(page, "Control comparisons result", "Mean diff"));
    await session.step(32, "And table \"Control comparisons result\" should have no missing values in \"95% CI low\" column", () => tableColumnComplete(page, "Control comparisons result", "95% CI low"));
    await session.step(33, "And table \"Control comparisons result\" should have no missing values in \"95% CI high\" column", () => tableColumnComplete(page, "Control comparisons result", "95% CI high"));
    await session.step(34, "And table \"Control comparisons result\" should have no missing values in \"t\" column", () => tableColumnComplete(page, "Control comparisons result", "t"));
    await session.step(35, "And table \"Control comparisons result\" should have no missing values in \"df\" column", () => tableColumnComplete(page, "Control comparisons result", "df"));
    await session.step(36, "And table \"Control comparisons result\" should have no missing values in \"p (raw)\" column", () => tableColumnComplete(page, "Control comparisons result", "p (raw)"));
    await session.step(37, "And table \"Control comparisons result\" should have no missing values in \"p (adj)\" column", () => tableColumnComplete(page, "Control comparisons result", "p (adj)"));
    await session.step(38, "And table \"Control comparisons result\" should have no missing values in \"Hedges' g\" column", () => tableColumnComplete(page, "Control comparisons result", "Hedges' g"));
    await session.step(39, "And second grid viewer should be bound to table \"Control comparisons result\"", () => boundTable(page, el("second grid viewer"), "Control comparisons result"));
    await session.step(40, "And the \"text of cell 1 of Group\" reading of second grid viewer should be \"Black\"", () => readingReads(page, "text of cell 1 of Group", el("second grid viewer"), "Black"));
    await session.step(41, "And the \"text of cell 1 of n\" reading of second grid viewer should be \"157\"", () => readingReads(page, "text of cell 1 of n", el("second grid viewer"), "157"));
    await session.step(42, "And the \"text of cell 2 of Group\" reading of second grid viewer should be \"Caucasian\"", () => readingReads(page, "text of cell 2 of Group", el("second grid viewer"), "Caucasian"));
    await session.step(43, "And the \"text of cell 2 of n\" reading of second grid viewer should be \"5266\"", () => readingReads(page, "text of cell 2 of n", el("second grid viewer"), "5266"));
    await session.step(44, "And the \"text of cell 3 of Group\" reading of second grid viewer should be \"Other\"", () => readingReads(page, "text of cell 3 of Group", el("second grid viewer"), "Other"));
    await session.step(45, "And the \"text of cell 3 of n\" reading of second grid viewer should be \"354\"", () => readingReads(page, "text of cell 3 of n", el("second grid viewer"), "354"));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
});
