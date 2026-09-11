/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/pareto-front.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [eda.viewer.pareto-front]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Pareto front viewer", () => {
  const session = feature(test, "features/pareto-front.feature", import.meta.url);
  test("The viewer picks the column of unique values as its label", {tag: ["@eda", "@realizes:eda.viewer.pareto-front"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(15, "When user picks \"ML > Pareto Front...\" from the top menu", () => pickFromTopMenu(page, "ML > Pareto Front..."));
    await session.step(16, "Then pareto front viewer should be visible", () => shouldBe(page, el("pareto front viewer"), "visible"));
    await session.step(17, "And \"Label Columns\" property of pareto front viewer should be \"model\"", () => propertyShouldBe(page, "Label Columns", el("pareto front viewer"), "model"));
    await session.step(18, "And \"Minimize\" property of pareto front viewer should be \"highway.mpg, price\"", () => propertyShouldBe(page, "Minimize", el("pareto front viewer"), "highway.mpg, price"));
    await session.step(19, "And no errors should have been logged", () => noErrors(page));
  });
  test("On demog the unique subject id is the label", {tag: ["@eda", "@realizes:eda.viewer.pareto-front"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(22, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(23, "When user picks \"ML > Pareto Front...\" from the top menu", () => pickFromTopMenu(page, "ML > Pareto Front..."));
    await session.step(24, "Then pareto front viewer should be visible", () => shouldBe(page, el("pareto front viewer"), "visible"));
    await session.step(25, "And \"Label Columns\" property of pareto front viewer should be \"USUBJID\"", () => propertyShouldBe(page, "Label Columns", el("pareto front viewer"), "USUBJID"));
  });
  test("Without a column of unique values the label stays empty", {tag: ["@eda", "@realizes:eda.viewer.pareto-front"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(28, "Given user opens iris dataset", () => openDataset(page, ds("iris")));
    await session.step(29, "When user picks \"ML > Pareto Front...\" from the top menu", () => pickFromTopMenu(page, "ML > Pareto Front..."));
    await session.step(30, "Then pareto front viewer should be visible", () => shouldBe(page, el("pareto front viewer"), "visible"));
    await session.step(31, "And \"Minimize\" property of pareto front viewer should be \"Petal.Length, Petal.Width\"", () => propertyShouldBe(page, "Minimize", el("pareto front viewer"), "Petal.Length, Petal.Width"));
    await session.step(32, "And \"Label Columns\" property of pareto front viewer should be \"\"", () => propertyShouldBe(page, "Label Columns", el("pareto front viewer"), ""));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
  });
});
