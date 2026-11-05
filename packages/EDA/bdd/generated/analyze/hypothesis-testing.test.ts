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
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {tableColumns, tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, noBalloons, noErrors, painted, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Two-sample t-test compares a feature between the two groups of a category", () => {
  const session = feature(test, "features/analyze/hypothesis-testing.feature", import.meta.url);
  test("Two-sample t-test compares a feature between the two groups of a category", {tag: ["@eda", "@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens a table \"groups\" with:", () => openTableOf(page, "groups", [["grp","site","v","w"],["a","x","1","5"],["a","y","2","7"],["a","z","3","6"],["b","x","11","6"],["b","y","12","5"],["b","z","13","7"]]), [["grp","site","v","w"],["a","x","1","5"],["a","y","2","7"],["a","z","3","6"],["b","x","11","6"],["b","y","12","5"],["b","z","13","7"]]);
    await run.scenario("The dialog opens on the two-group column, a three-group one disables Run, and Run docks the box plot and the result table", async () => {
      await session.step(32, "When user picks \"ML > Analyze > Group Comparison > T-test...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Group Comparison > T-test..."));
      await session.step(33, "Then \"Two-sample t-test\" dialog should be visible", () => shouldBe(page, el("\"Two-sample t-test\" dialog"), "visible"));
      await session.step(34, "And editor of Category input in \"Two-sample t-test\" dialog should have text \"grp\"", () => shouldHaveText(page, el("editor of Category input in \"Two-sample t-test\" dialog"), "grp"));
      await session.step(35, "And editor of Feature input in \"Two-sample t-test\" dialog should have text \"v\"", () => shouldHaveText(page, el("editor of Feature input in \"Two-sample t-test\" dialog"), "v"));
      await session.step(36, "And Alpha input in \"Two-sample t-test\" dialog should have value \"0.05\"", () => shouldHaveValue(page, el("Alpha input in \"Two-sample t-test\" dialog"), "0.05"));
      await session.step(37, "And \"Full report\" input in \"Two-sample t-test\" dialog should be checked", () => shouldBe(page, el("\"Full report\" input in \"Two-sample t-test\" dialog"), "checked"));
      await session.step(38, "And Run button in \"Two-sample t-test\" dialog should be enabled", () => shouldBe(page, el("Run button in \"Two-sample t-test\" dialog"), "enabled"));
      await session.step(39, "When user selects \"site\" in Category input in \"Two-sample t-test\" dialog", () => selectIn(page, "site", el("Category input in \"Two-sample t-test\" dialog")));
      await session.step(40, "Then Run button in \"Two-sample t-test\" dialog should be disabled", () => shouldBe(page, el("Run button in \"Two-sample t-test\" dialog"), "disabled"));
      await session.step(41, "When user selects \"grp\" in Category input in \"Two-sample t-test\" dialog", () => selectIn(page, "grp", el("Category input in \"Two-sample t-test\" dialog")));
      await session.step(42, "Then Run button in \"Two-sample t-test\" dialog should be enabled", () => shouldBe(page, el("Run button in \"Two-sample t-test\" dialog"), "enabled"));
      await session.step(43, "When user clicks on Run button in \"Two-sample t-test\" dialog", () => clickOn(page, el("Run button in \"Two-sample t-test\" dialog")));
      await session.step(44, "Then \"Two-sample t-test\" dialog should be hidden", () => shouldBe(page, el("\"Two-sample t-test\" dialog"), "hidden"));
      await session.step(45, "And box plot viewer should be visible", () => shouldBe(page, el("box plot viewer"), "visible"));
      await session.step(46, "And description of box plot viewer should contain text \"\\\"v\\\" differs between \\\"a\\\" and \\\"b\\\"\"", () => shouldContainText(page, el("description of box plot viewer"), "\"v\" differs between \"a\" and \"b\""));
      await session.step(47, "And box plot viewer should be painted", () => painted(page, el("box plot viewer")));
      await session.step(48, "And table \"T-test result\" should be open", () => tableOpen(page, "T-test result"));
      await session.step(49, "And table \"T-test result\" should have 1 row", () => tableRows(page, "T-test result", 1));
      await session.step(50, "And table \"T-test result\" should have columns \"Conclusion, t, df, p-value, Mean difference, 95% CI low, 95% CI high, Cohen's d, Hedges' g\"", () => tableColumns(page, "T-test result", "Conclusion, t, df, p-value, Mean difference, 95% CI low, 95% CI high, Cohen's d, Hedges' g"));
      await session.step(51, "And second grid viewer should be bound to table \"T-test result\"", () => boundTable(page, el("second grid viewer"), "T-test result"));
      await session.step(52, "And the \"text of cell 1 of Conclusion\" reading of second grid viewer should be \"Significant\"", () => readingReads(page, "text of cell 1 of Conclusion", el("second grid viewer"), "Significant"));
      await session.step(53, "And table \"groups\" should have 6 rows", () => tableRows(page, "groups", 6));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Overlapping groups are reported as not differing", async () => {
      await session.step(58, "When user closes all views", () => closeAllViews(page));
      await session.step(59, "And user opens a table \"groups\" with:", () => openTableOf(page, "groups", [["grp","site","v","w"],["a","x","1","5"],["a","y","2","7"],["a","z","3","6"],["b","x","11","6"],["b","y","12","5"],["b","z","13","7"]]), [["grp","site","v","w"],["a","x","1","5"],["a","y","2","7"],["a","z","3","6"],["b","x","11","6"],["b","y","12","5"],["b","z","13","7"]]);
      await session.step(67, "When user picks \"ML > Analyze > Group Comparison > T-test...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Group Comparison > T-test..."));
      await session.step(68, "And user selects \"w\" in Feature input in \"Two-sample t-test\" dialog", () => selectIn(page, "w", el("Feature input in \"Two-sample t-test\" dialog")));
      await session.step(69, "And user clicks on Run button in \"Two-sample t-test\" dialog", () => clickOn(page, el("Run button in \"Two-sample t-test\" dialog")));
      await session.step(70, "Then box plot viewer should be visible", () => shouldBe(page, el("box plot viewer"), "visible"));
      await session.step(71, "And description of box plot viewer should contain text \"\\\"w\\\" doesn't differ between \\\"a\\\" and \\\"b\\\"\"", () => shouldContainText(page, el("description of box plot viewer"), "\"w\" doesn't differ between \"a\" and \"b\""));
      await session.step(72, "And the \"text of cell 1 of Conclusion\" reading of second grid viewer should be \"Not significant\"", () => readingReads(page, "text of cell 1 of Conclusion", el("second grid viewer"), "Not significant"));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Student with Full report off docks only the box plot", async () => {
      await session.step(76, "When user closes all views", () => closeAllViews(page));
      await session.step(77, "And user opens a table \"groups\" with:", () => openTableOf(page, "groups", [["grp","site","v","w"],["a","x","1","5"],["a","y","2","7"],["a","z","3","6"],["b","x","11","6"],["b","y","12","5"],["b","z","13","7"]]), [["grp","site","v","w"],["a","x","1","5"],["a","y","2","7"],["a","z","3","6"],["b","x","11","6"],["b","y","12","5"],["b","z","13","7"]]);
      await session.step(85, "When user picks \"ML > Analyze > Group Comparison > T-test...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Group Comparison > T-test..."));
      await session.step(86, "And user selects \"Student\" in Method input in \"Two-sample t-test\" dialog", () => selectIn(page, "Student", el("Method input in \"Two-sample t-test\" dialog")));
      await session.step(87, "And user unchecks \"Full report\" input in \"Two-sample t-test\" dialog", () => uncheck(page, el("\"Full report\" input in \"Two-sample t-test\" dialog")));
      await session.step(88, "And user clicks on Run button in \"Two-sample t-test\" dialog", () => clickOn(page, el("Run button in \"Two-sample t-test\" dialog")));
      await session.step(89, "Then box plot viewer should be visible", () => shouldBe(page, el("box plot viewer"), "visible"));
      await session.step(90, "And description of box plot viewer should contain text \"\\\"v\\\" differs between \\\"a\\\" and \\\"b\\\"\"", () => shouldContainText(page, el("description of box plot viewer"), "\"v\" differs between \"a\" and \"b\""));
      await session.step(91, "And the open tableview should have 1 grid viewer", () => viewerCount(page, 1, "grid"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
