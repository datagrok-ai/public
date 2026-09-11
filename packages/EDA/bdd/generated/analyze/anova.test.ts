/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/anova.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.analyze.group-comparison.anova]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {tableColumns, tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, noBalloons, noErrors, painted, propertyShouldContain} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("One-way ANOVA", () => {
  const session = feature(test, "features/analyze/anova.feature", import.meta.url);
  test("The dialog opens on a category, a feature and a significance level", {tag: ["@eda", "@realizes:ml.menu.analyze.group-comparison.anova"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(16, "When user picks \"ML > Analyze > Group Comparison > ANOVA...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Group Comparison > ANOVA..."));
    await session.step(17, "Then \"ANOVA\" dialog should be visible", () => shouldBe(page, el("\"ANOVA\" dialog"), "visible"));
    await session.step(18, "And editor of Category input in \"ANOVA\" dialog should have text \"RACE\"", () => shouldHaveText(page, el("editor of Category input in \"ANOVA\" dialog"), "RACE"));
    await session.step(19, "And editor of Feature input in \"ANOVA\" dialog should have text \"AGE\"", () => shouldHaveText(page, el("editor of Feature input in \"ANOVA\" dialog"), "AGE"));
    await session.step(20, "And Alpha input in \"ANOVA\" dialog should have value \"0.05\"", () => shouldHaveValue(page, el("Alpha input in \"ANOVA\" dialog"), "0.05"));
    await session.step(21, "And Run button in \"ANOVA\" dialog should be enabled", () => shouldBe(page, el("Run button in \"ANOVA\" dialog"), "enabled"));
  });
  test("Running it docks a box plot with the conclusion and the table of the test", {tag: ["@eda", "@realizes:ml.menu.analyze.group-comparison.anova"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(24, "When user picks \"ML > Analyze > Group Comparison > ANOVA...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Group Comparison > ANOVA..."));
    await session.step(25, "And user clicks on Run button in \"ANOVA\" dialog", () => clickOn(page, el("Run button in \"ANOVA\" dialog")));
    await session.step(26, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(27, "And \"ANOVA\" dialog should be hidden", () => shouldBe(page, el("\"ANOVA\" dialog"), "hidden"));
    await session.step(28, "And box plot viewer should be visible", () => shouldBe(page, el("box plot viewer"), "visible"));
    await session.step(29, "And \"Description\" property of box plot viewer should contain \"doesn't affect\"", () => propertyShouldContain(page, "Description", el("box plot viewer"), "doesn't affect"));
    await session.step(30, "And \"Description\" property of box plot viewer should contain \"p = 0.176\"", () => propertyShouldContain(page, "Description", el("box plot viewer"), "p = 0.176"));
    await session.step(31, "And box plot viewer should be painted", () => painted(page, el("box plot viewer")));
    await session.step(32, "And table \"ANOVA result\" should be open", () => tableOpen(page, "ANOVA result"));
    await session.step(33, "And table \"ANOVA result\" should have 1 row", () => tableRows(page, "ANOVA result", 1));
    await session.step(34, "And table \"ANOVA result\" should have columns \"Conclusion, Source of variance, F, df₁, df₂, F-critical, p-value\"", () => tableColumns(page, "ANOVA result", "Conclusion, Source of variance, F, df₁, df₂, F-critical, p-value"));
    await session.step(35, "And second grid viewer should be bound to table \"ANOVA result\"", () => boundTable(page, el("second grid viewer"), "ANOVA result"));
    await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
});
