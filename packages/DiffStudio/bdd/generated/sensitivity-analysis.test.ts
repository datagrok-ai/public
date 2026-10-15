/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sensitivity-analysis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [diffstudio.model.bioreactor]
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openLibraryModel} from '../bindings/diff-studio.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldBeSwitchedOn, switchOff, switchOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {viewHoldsViewers, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Sensitivity analysis over a model", () => {
  const session = feature(test, "features/sensitivity-analysis.feature", import.meta.url);
  test("Sensitivity analysis over a model", {tag: ["@journey", "@diffstudio", "@realizes:diffstudio.model.bioreactor"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens the \"Bioreactor\" model of the Diff Studio library", () => openLibraryModel(page, "Bioreactor"));
    await run.scenario("The model is the one the ribbon acts on", async () => {
      await session.step(18, "Then the \"Bioreactor\" view should be current", () => viewIsCurrent(page, "Bioreactor"));
      await session.step(19, "And Sensitivity ribbon item should be visible", () => shouldBe(page, el("Sensitivity ribbon item"), "visible"));
      await session.step(20, "And Fit ribbon item should be visible", () => shouldBe(page, el("Fit ribbon item"), "visible"));
    });
    await run.scenario("Sensitivity opens a view of its own", async () => {
      await session.step(23, "When user clicks on Sensitivity ribbon item", () => clickOn(page, el("Sensitivity ribbon item")));
      await session.step(24, "Then the \"Bioreactor - comparison\" view should be current", () => viewIsCurrent(page, "Bioreactor - comparison"));
    });
    await run.scenario("A parameter is chosen by the switch beside its input", async () => {
      await session.step(27, "When user switches on \"FFox\" input", () => switchOn(page, el("\"FFox\" input")));
      await session.step(28, "Then \"FFox min\" input should be switched on", () => shouldBeSwitchedOn(page, el("\"FFox min\" input")));
      await session.step(29, "And \"FFox max\" input should be visible", () => shouldBe(page, el("\"FFox max\" input"), "visible"));
      await session.step(30, "When user switches off \"FFox min\" input", () => switchOff(page, el("\"FFox min\" input")));
      await session.step(31, "Then \"FFox\" input should be visible", () => shouldBe(page, el("\"FFox\" input"), "visible"));
    });
    await run.scenario("Running the analysis over three parameters puts its viewers on screen", async () => {
      await session.step(34, "When user switches on \"FFox\" input", () => switchOn(page, el("\"FFox\" input")));
      await session.step(35, "And user switches on \"FKox\" input", () => switchOn(page, el("\"FKox\" input")));
      await session.step(36, "And user switches on \"FFred\" input", () => switchOn(page, el("\"FFred\" input")));
      await session.step(37, "And user clicks on \"Run\" icon", () => clickOn(page, el("\"Run\" icon")));
      await session.step(38, "Then the current view should hold at least 4 viewers", () => viewHoldsViewers(page, 4));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
