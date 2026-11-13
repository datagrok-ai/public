/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/hypothesis-testing.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {tableColumns, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Two-sample t-test compares a feature between the two groups of a category", () => {
  const session = feature(test, "features/analyze/hypothesis-testing.feature", import.meta.url);
  test("The dialog opens on the two-group column, a three-group one disables Run, and Run docks the box plot and the result table", {tag: ["@eda"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens a table \"groups\" with:", () => openTableOf(page, "groups", [["grp","site","v"],["a","x","1"],["a","y","2"],["a","z","3"],["b","x","11"],["b","y","12"],["b","z","13"]]), [["grp","site","v"],["a","x","1"],["a","y","2"],["a","z","3"],["b","x","11"],["b","y","12"],["b","z","13"]]);
    await session.step(24, "When user picks \"ML > Analyze > Group Comparison > T-test...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Group Comparison > T-test..."));
    await session.step(25, "Then \"Two-sample t-test\" dialog should be visible", () => shouldBe(page, el("\"Two-sample t-test\" dialog"), "visible"));
    await session.step(26, "And editor of Category input in \"Two-sample t-test\" dialog should have text \"grp\"", () => shouldHaveText(page, el("editor of Category input in \"Two-sample t-test\" dialog"), "grp"));
    await session.step(27, "And editor of Feature input in \"Two-sample t-test\" dialog should have text \"v\"", () => shouldHaveText(page, el("editor of Feature input in \"Two-sample t-test\" dialog"), "v"));
    await session.step(28, "When user selects \"site\" in Category input in \"Two-sample t-test\" dialog", () => selectIn(page, "site", el("Category input in \"Two-sample t-test\" dialog")));
    await session.step(29, "Then Run button in \"Two-sample t-test\" dialog should be disabled", () => shouldBe(page, el("Run button in \"Two-sample t-test\" dialog"), "disabled"));
    await session.step(30, "When user selects \"grp\" in Category input in \"Two-sample t-test\" dialog", () => selectIn(page, "grp", el("Category input in \"Two-sample t-test\" dialog")));
    await session.step(31, "Then Run button in \"Two-sample t-test\" dialog should be enabled", () => shouldBe(page, el("Run button in \"Two-sample t-test\" dialog"), "enabled"));
    await session.step(32, "When user clicks on Run button in \"Two-sample t-test\" dialog", () => clickOn(page, el("Run button in \"Two-sample t-test\" dialog")));
    await session.step(33, "Then description of box plot viewer should contain text \"\\\"v\\\" differs between \\\"a\\\" and \\\"b\\\"\"", () => shouldContainText(page, el("description of box plot viewer"), "\"v\" differs between \"a\" and \"b\""));
    await session.step(34, "And table \"T-test result\" should have 1 row", () => tableRows(page, "T-test result", 1));
    await session.step(35, "And table \"T-test result\" should have columns \"Conclusion, t, df, p-value, Mean difference, 95% CI low, 95% CI high, Cohen's d, Hedges' g\"", () => tableColumns(page, "T-test result", "Conclusion, t, df, p-value, Mean difference, 95% CI low, 95% CI high, Cohen's d, Hedges' g"));
    await session.step(36, "And second grid viewer should be bound to table \"T-test result\"", () => boundTable(page, el("second grid viewer"), "T-test result"));
    await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
});
