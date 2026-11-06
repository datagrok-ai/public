/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/flows/plan-to-analysis.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {finishedRunning, noSavedPlan, pickSavedPlan, planIsSaved} from '../../bindings/pmx.js';
import {openApp} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("From the plan to the analysis", () => {
  const session = feature(test, "features/flows/plan-to-analysis.feature", import.meta.url);
  test("From the plan to the analysis", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await run.scenario("The saved plan is what the analysis compares with", async () => {
      await session.step(6, "Given user is logged in", () => loggedIn(page));
      await session.step(7, "And no plan named \"BDD flow\" is saved", () => noSavedPlan(page, "BDD flow"));
      await session.step(8, "And user opens PMX Studio", () => openApp(page));
      await session.step(9, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(10, "When user clicks on \"Save plan\" button", () => clickOn(page, el("\"Save plan\" button")));
      await session.step(11, "And user enters \"BDD flow\" into \"Name\" input in \"Save plan\" dialog", () => enterInto(page, "BDD flow", el("\"Name\" input in \"Save plan\" dialog")));
      await session.step(12, "And user clicks on \"Save\" button in \"Save plan\" dialog", () => clickOn(page, el("\"Save\" button in \"Save plan\" dialog")));
      await session.step(13, "Then a plan named \"BDD flow\" should be saved", () => planIsSaved(page, "BDD flow"));
      await session.step(14, "When user clicks on \"3 Analyze\" button", () => clickOn(page, el("\"3 Analyze\" button")));
      await session.step(15, "Then analyze start should be visible", () => shouldBe(page, el("analyze start"), "visible"));
      await session.step(16, "When user clicks on \"Analyze an example study\" button", () => clickOn(page, el("\"Analyze an example study\" button")));
      await session.step(17, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(18, "When user picks the saved plan \"BDD flow\"", () => pickSavedPlan(page, "BDD flow"));
      await session.step(19, "Then stale note should be visible", () => shouldBe(page, el("stale note"), "visible"));
      await session.step(20, "When user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(21, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(22, "And study error should be hidden", () => shouldBe(page, el("study error"), "hidden"));
      await session.step(23, "And stale note should be hidden", () => shouldBe(page, el("stale note"), "hidden"));
      await session.step(24, "When user clicks on \"Plan vs observed\" tab", () => clickOn(page, el("\"Plan vs observed\" tab")));
      await session.step(25, "Then results pane should contain text \"plan: BDD flow\"", () => shouldContainText(page, el("results pane"), "plan: BDD flow"));
      await session.step(26, "And \"Effect of each treatment\" section should be visible", () => shouldBe(page, el("\"Effect of each treatment\" section"), "visible"));
    });
    run.finish();
  });
});
