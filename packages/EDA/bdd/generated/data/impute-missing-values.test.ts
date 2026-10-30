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
import {clickOn, enterInto, shouldBe, shouldContainText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnIncomplete, columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("ML > Impute Missing Values fills the missing values of a column from its nearest neighbours", () => {
  const session = feature(test, "features/data/impute-missing-values.feature", import.meta.url);
  test("The dialog opens on the column with missing values, and RUN fills it in place", {tag: ["@eda"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens a table \"gaps\" with:", () => openTableOf(page, "gaps", [["a","b"],["1","10"],["2","20"],["5",""],["9","90"],["10","100"]]), [["a","b"],["1","10"],["2","20"],["5",""],["9","90"],["10","100"]]);
    await session.step(28, "When user picks \"ML > Impute Missing Values...\" from the top menu", () => pickFromTopMenu(page, "ML > Impute Missing Values..."));
    await session.step(29, "Then \"k-NN Imputation\" dialog should be visible", () => shouldBe(page, el("\"k-NN Imputation\" dialog"), "visible"));
    await session.step(30, "And Impute input in \"k-NN Imputation\" dialog should contain text \"(1)\"", () => shouldContainText(page, el("Impute input in \"k-NN Imputation\" dialog"), "(1)"));
    await session.step(31, "And Using input in \"k-NN Imputation\" dialog should contain text \"(2)\"", () => shouldContainText(page, el("Using input in \"k-NN Imputation\" dialog"), "(2)"));
    await session.step(32, "And Neighbors input in \"k-NN Imputation\" dialog should have value \"4\"", () => shouldHaveValue(page, el("Neighbors input in \"k-NN Imputation\" dialog"), "4"));
    await session.step(33, "And \"In-place\" input in \"k-NN Imputation\" dialog should be checked", () => shouldBe(page, el("\"In-place\" input in \"k-NN Imputation\" dialog"), "checked"));
    await session.step(34, "And RUN button in \"k-NN Imputation\" dialog should be enabled", () => shouldBe(page, el("RUN button in \"k-NN Imputation\" dialog"), "enabled"));
    await session.step(35, "When user clicks on RUN button in \"k-NN Imputation\" dialog", () => clickOn(page, el("RUN button in \"k-NN Imputation\" dialog")));
    await session.step(36, "Then \"k-NN Imputation\" dialog should be hidden", () => shouldBe(page, el("\"k-NN Imputation\" dialog"), "hidden"));
    await session.step(37, "And \"b\" column should have no missing values", () => columnComplete(page, "b"));
    await session.step(38, "And the value of \"b\" column in row 3 should be \"55\"", () => valueInRow(page, "b", 3, "55"));
    await session.step(39, "And the value of \"b\" column in row 2 should be \"20\"", () => valueInRow(page, "b", 2, "20"));
    await session.step(40, "And the value of \"a\" column in row 3 should be \"5\"", () => valueInRow(page, "a", 3, "5"));
    await session.step(41, "And the table should have the columns \"a, b\"", () => columnsExactly(page, "a, b"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
  test("With In-place off and one neighbour, RUN adds the filled copy and leaves the column as it was", {tag: ["@eda"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens a table \"gaps\" with:", () => openTableOf(page, "gaps", [["a","b"],["1","10"],["2","20"],["5",""],["9","90"],["10","100"]]), [["a","b"],["1","10"],["2","20"],["5",""],["9","90"],["10","100"]]);
    await session.step(45, "When user picks \"ML > Impute Missing Values...\" from the top menu", () => pickFromTopMenu(page, "ML > Impute Missing Values..."));
    await session.step(46, "And user enters \"1\" into Neighbors input in \"k-NN Imputation\" dialog", () => enterInto(page, "1", el("Neighbors input in \"k-NN Imputation\" dialog")));
    await session.step(47, "And user unchecks \"In-place\" input in \"k-NN Imputation\" dialog", () => uncheck(page, el("\"In-place\" input in \"k-NN Imputation\" dialog")));
    await session.step(48, "And user clicks on RUN button in \"k-NN Imputation\" dialog", () => clickOn(page, el("RUN button in \"k-NN Imputation\" dialog")));
    await session.step(49, "Then a new column \"b(copy)\" should have been added", () => newColumnNamed(page, "b(copy)"));
    await session.step(50, "And the value of \"b(copy)\" column in row 3 should be \"20\"", () => valueInRow(page, "b(copy)", 3, "20"));
    await session.step(51, "And the value of \"b(copy)\" column in row 4 should be \"90\"", () => valueInRow(page, "b(copy)", 4, "90"));
    await session.step(52, "And \"b\" column should have missing values", () => columnIncomplete(page, "b"));
    await session.step(53, "And the table should have the columns \"a, b, b(copy)\"", () => columnsExactly(page, "a, b, b(copy)"));
    await session.step(54, "And no errors should have been logged", () => noErrors(page));
  });
  test("Nothing to impute disables RUN", {tag: ["@eda"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens a table \"gaps\" with:", () => openTableOf(page, "gaps", [["a","b"],["1","10"],["2","20"],["5",""],["9","90"],["10","100"]]), [["a","b"],["1","10"],["2","20"],["5",""],["9","90"],["10","100"]]);
    await session.step(57, "When user picks \"ML > Impute Missing Values...\" from the top menu", () => pickFromTopMenu(page, "ML > Impute Missing Values..."));
    await session.step(58, "And user clicks on Impute input in \"k-NN Imputation\" dialog", () => clickOn(page, el("Impute input in \"k-NN Imputation\" dialog")));
    await session.step(59, "And user toggles the \"b\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "b", el("\"Select columns...\" dialog")));
    await session.step(60, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(61, "Then RUN button in \"k-NN Imputation\" dialog should be disabled", () => shouldBe(page, el("RUN button in \"k-NN Imputation\" dialog"), "disabled"));
    await session.step(62, "When user clicks on CANCEL button in \"k-NN Imputation\" dialog", () => clickOn(page, el("CANCEL button in \"k-NN Imputation\" dialog")));
    await session.step(63, "Then \"k-NN Imputation\" dialog should be hidden", () => shouldBe(page, el("\"k-NN Imputation\" dialog"), "hidden"));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
  test("A table with no missing values gets an info balloon and no dialog", {tag: ["@eda"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens a table \"gaps\" with:", () => openTableOf(page, "gaps", [["a","b"],["1","10"],["2","20"],["5",""],["9","90"],["10","100"]]), [["a","b"],["1","10"],["2","20"],["5",""],["9","90"],["10","100"]]);
    await session.step(67, "Given user opens a table \"full\" with:", () => openTableOf(page, "full", [["a","b"],["1","10"],["2","20"]]), [["a","b"],["1","10"],["2","20"]]);
    await session.step(71, "When user picks \"ML > Impute Missing Values...\" from the top menu", () => pickFromTopMenu(page, "ML > Impute Missing Values..."));
    await session.step(72, "Then an info balloon containing \"No missing values\" should have been shown", () => infoBalloonText(page, "No missing values"));
    await session.step(73, "And \"k-NN Imputation\" dialog should be absent", () => shouldBe(page, el("\"k-NN Imputation\" dialog"), "absent"));
    await session.step(74, "And no errors should have been logged", () => noErrors(page));
  });
});
