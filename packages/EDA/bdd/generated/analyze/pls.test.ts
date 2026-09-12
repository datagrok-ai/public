/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/pls.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.analyze.pls]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Partial least squares regression", () => {
  const session = feature(test, "features/analyze/pls.feature", import.meta.url);
  test("Partial least squares regression", {tag: ["@journey", "@eda", "@realizes:ml.menu.analyze.pls"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await run.scenario("The dialog opens on what to predict, the predictors and the components", async () => {
      await session.step(16, "When user picks \"ML > Analyze > PLS...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > PLS..."));
      await session.step(17, "Then \"PLS\" dialog should be visible", () => shouldBe(page, el("\"PLS\" dialog"), "visible"));
      await session.step(18, "And editor of Predict input in \"PLS\" dialog should have text \"price\"", () => shouldHaveText(page, el("editor of Predict input in \"PLS\" dialog"), "price"));
      await session.step(19, "And editor of Using input in \"PLS\" dialog should contain text \"(15)\"", () => shouldContainText(page, el("editor of Using input in \"PLS\" dialog"), "(15)"));
      await session.step(20, "And Components input in \"PLS\" dialog should have value \"3\"", () => shouldHaveValue(page, el("Components input in \"PLS\" dialog"), "3"));
      await session.step(21, "And Quadratic input in \"PLS\" dialog should be unchecked", () => shouldBe(page, el("Quadratic input in \"PLS\" dialog"), "unchecked"));
      await session.step(22, "And RUN button in \"PLS\" dialog should be enabled", () => shouldBe(page, el("RUN button in \"PLS\" dialog"), "enabled"));
    });
    await run.scenario("Every column as a predictor includes the predicted one, which RUN does not accept", async () => {
      await session.step(25, "When user clicks on editor of Using input in \"PLS\" dialog", () => clickOn(page, el("editor of Using input in \"PLS\" dialog")));
      await session.step(26, "And user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(27, "Then \"Select columns...\" dialog should contain text \"16 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "16 checked"));
      await session.step(28, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(29, "Then editor of Using input in \"PLS\" dialog should contain text \"(16)\"", () => shouldContainText(page, el("editor of Using input in \"PLS\" dialog"), "(16)"));
      await session.step(30, "And RUN button in \"PLS\" dialog should be disabled", () => shouldBe(page, el("RUN button in \"PLS\" dialog"), "disabled"));
    });
    await run.scenario("With price taken back out, RUN adds three PLS components", async () => {
      await session.step(33, "When user clicks on editor of Using input in \"PLS\" dialog", () => clickOn(page, el("editor of Using input in \"PLS\" dialog")));
      await session.step(34, "And user types \"price\" into Search input in \"Select columns...\" dialog", () => typeInto(page, "price", el("Search input in \"Select columns...\" dialog")));
      await session.step(35, "Then the \"text of cell 16 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"price\"", () => readingReads(page, "text of cell 16 of __name", el("grid viewer in \"Select columns...\" dialog"), "price"));
      await session.step(36, "When user clicks on the \"cell 16 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 16 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(37, "Then the \"text of cell 16 of x\" reading of grid viewer in \"Select columns...\" dialog should be \"false\"", () => readingReads(page, "text of cell 16 of x", el("grid viewer in \"Select columns...\" dialog"), "false"));
      await session.step(38, "And \"Select columns...\" dialog should contain text \"15 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "15 checked"));
      await session.step(39, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(40, "Then editor of Using input in \"PLS\" dialog should contain text \"(15)\"", () => shouldContainText(page, el("editor of Using input in \"PLS\" dialog"), "(15)"));
      await session.step(41, "And RUN button in \"PLS\" dialog should be enabled", () => shouldBe(page, el("RUN button in \"PLS\" dialog"), "enabled"));
      await session.step(42, "When user clicks on RUN button in \"PLS\" dialog", () => clickOn(page, el("RUN button in \"PLS\" dialog")));
      await session.step(43, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(44, "And \"PLS\" dialog should be hidden", () => shouldBe(page, el("\"PLS\" dialog"), "hidden"));
      await session.step(45, "And 3 new columns should have been added", () => newColumnsCount(page, 3));
      await session.step(46, "And the table should have a column \"PLS1\"", () => hasColumn(page, "PLS1"));
      await session.step(47, "And the table should have a column \"PLS2\"", () => hasColumn(page, "PLS2"));
      await session.step(48, "And the table should have a column \"PLS3\"", () => hasColumn(page, "PLS3"));
      await session.step(49, "And \"PLS1\" column should have no missing values", () => columnComplete(page, "PLS1"));
      await session.step(50, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
