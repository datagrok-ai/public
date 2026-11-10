/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/data/impute-missing-values.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnsExactly} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("ML > Impute Missing Values fills the missing values of a column from its nearest neighbours", () => {
  const session = feature(test, "features/data/impute-missing-values.feature", import.meta.url);
  test("ML > Impute Missing Values fills the missing values of a column from its nearest neighbours", {tag: ["@eda", "@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens a table \"gaps\" with:", () => openTableOf(page, "gaps", [["a","b"],["1","10"],["2","20"],["5",""],["9","90"],["10","100"]]), [["a","b"],["1","10"],["2","20"],["5",""],["9","90"],["10","100"]]);
    await run.scenario("The dialog opens on the column with missing values, and RUN fills it in place", async () => {
      await session.step(23, "When user picks \"ML > Impute Missing Values...\" from the top menu", () => pickFromTopMenu(page, "ML > Impute Missing Values..."));
      await session.step(24, "Then \"k-NN Imputation\" dialog should be visible", () => shouldBe(page, el("\"k-NN Imputation\" dialog"), "visible"));
      await session.step(25, "And Impute input in \"k-NN Imputation\" dialog should contain text \"(1)\"", () => shouldContainText(page, el("Impute input in \"k-NN Imputation\" dialog"), "(1)"));
      await session.step(26, "When user clicks on RUN button in \"k-NN Imputation\" dialog", () => clickOn(page, el("RUN button in \"k-NN Imputation\" dialog")));
      await session.step(27, "Then \"b\" column should have no missing values", () => columnComplete(page, "b"));
      await session.step(28, "And the table should have the columns \"a, b\"", () => columnsExactly(page, "a, b"));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A table with no missing values gets an info balloon and no dialog", async () => {
      await session.step(32, "Given user opens a table \"full\" with:", () => openTableOf(page, "full", [["a","b"],["1","10"],["2","20"]]), [["a","b"],["1","10"],["2","20"]]);
      await session.step(36, "When user picks \"ML > Impute Missing Values...\" from the top menu", () => pickFromTopMenu(page, "ML > Impute Missing Values..."));
      await session.step(37, "Then an info balloon containing \"No missing values\" should have been shown", () => infoBalloonText(page, "No missing values"));
      await session.step(38, "And \"k-NN Imputation\" dialog should be absent", () => shouldBe(page, el("\"k-NN Imputation\" dialog"), "absent"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
